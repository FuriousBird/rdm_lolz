import numpy as np
import utils.constants as cst
import datetime as dt

# https://youtu.be/MfT0YCWPev8 - listen to this to cope with the terrible code

class Noeud:
    """
    A class representing a point in 3D space.
    It does not handle any file saving/formatting
    """
    def __init__(self, loc):
        self.loc = np.array(loc, dtype=float)
    def __add__(self, other):
        return Noeud(self.loc + other.loc)
    def __repr__(self):
        return f"Point({self.loc})"
    def __sub__(self, other):
        return Noeud(self.loc - other.loc)
    def __mul__(self, scalar):
        return Noeud(self.loc * scalar)
    def __truediv__(self, scalar):
        return Noeud(self.loc / scalar)
    def close(self, other, threshold=1e-6):
        return self.distance(other) < threshold
    def distance(self, other):
        return np.linalg.norm(self.loc - other.loc)
    def __xor__(self, other):
        return Noeud(np.linalg.cross(self.loc, other.loc))
    def norm(self):
        return np.linalg.norm(self.loc)
    def normalize(self):
        norm = self.norm()
        if norm == 0:
            raise ValueError("Cannot normalize a zero vector")
        return Noeud(self.loc / norm)
    def frmt(self):
        return f"  {self.idx + 1} {self.loc[0]:.11E} {self.loc[1]:.11E} {self.loc[2]:.11E}\n"

class Poutre:
    """A class representing a beam defined by two points in 3D space.
    It does not handle any file saving/formatting
    """
    def __init__(self, point1, point2):
        self.point1 = point1
        self.point2 = point2
    def __iter__(self, deg=0, keep=False): #can you even iterate over a beam?
        if deg == 0:
            if keep:
                yield self.point1 
                yield self.point2
            else:
                yield Noeud(self.point1.loc)
                yield Noeud(self.point2.loc)
        else:
            for i in range(deg+2):
                if not keep:
                    yield Noeud(self.point1.loc + (self.point2.loc - self.point1.loc) * (i / (deg + 1)))
                    continue
                if i == 0:
                    yield self.point1
                    continue
                if i == deg + 1:
                    yield self.point2

    def __repr__(self):
        return f"Poutre({self.point1}, {self.point2})"
    def length(self):
        return self.point1.distance(self.point2)
    def vecR(self, other=Noeud((0, 0, 1))):
        return ((self.point2-self.point1)^other).normalize()
    def frmt(self):
        n1 = self.point1.idx + 1
        n2 = self.point2.idx + 1
        vec = self.vecR().loc
        return f"   {n1:2d} {cst.POUTRE_TYPE}     {n1:2d}    {n2:2d} {vec[0]:.11E} {vec[1]:.11E} {vec[2]:.11E} 11 11\n"

class Block:
    def __init__(self, name):
        self.name = name
    #iteration
    def __iter__(self):
        yield self.name


class Fichier:
    def __init__(self, *args):
        self.sections = [i for i in args if isinstance(i, Block)]
        self.meta = {}
        self.init_meta() #THE META OF THE HOLY CROSS
        for sect in self.sections:
            self.bind(sect)
    def init_meta(self):
        self.meta[cst.POINTS_COUNT_KEY] = 0
        self.meta[cst.POUTRES_COUNT_KEY] = 0
        self.meta[cst.LIAISONS_COUNT_KEY] = 0
    def bind(self, section):
        if isinstance(section, Block):
            section.parentfile = self
            self.sections.append(section)
            if hasattr(section, '_setup'):
                section._setup(self)
        else:
            raise TypeError("Only Block instances can be bound to Fichier")
    def __iter__(self):
        return iter(self.sections)
    def write(self, file, info=cst.DEFAULT_INFO):
        file.write(info)
        file.write(cst.FILE_START)
        file.write(cst.DEFAULT_FILE_HEADER.format(""))
        nodesblocks = [s for s in self.sections if isinstance(s, Noeuds)]
        KEY = cst.POINTS_COUNT_KEY
        file.write(f"noeuds {self.meta.get(KEY, 0)}\n")
        for node in nodesblocks:
            node.write(file)
        file.write(f"  0\n////\n")

        poutresblocks = [s for s in self.sections if isinstance(s, Poutres)]
        KEY = cst.POUTRES_COUNT_KEY
        file.write(f"poutres {self.meta.get(KEY, 0)}\n")
        for poutre in poutresblocks:
            poutre.write(file)
        file.write(f"  0\n////\n")
        liaisonsblocks = [s for s in self.sections if isinstance(s, Liaisons)]

        KEY = cst.LIAISONS_COUNT_KEY
        file.write(f"liaisons {self.meta.get(KEY, 0)}\n")
        for liaison in liaisonsblocks:
            liaison.write(file)
        file.write(f"  0\n////\n")
        file.write(f"fin du fichier\n")
        file.write(f"{cst.FILE_END}")

class Noeuds(Block):
    def __init__(self, *args, autojoin=False, closed=False):
        self.name="nodes"
        self.nodes = [i for i in args if isinstance(i, Noeud)]
        self.autojoin = autojoin
        self.closed = closed
    def _setup(self, parentfile:Fichier):
        self.parentfile = parentfile
        for point in self.nodes:
            point.idx = parentfile.meta[cst.POINTS_COUNT_KEY]
            parentfile.meta[cst.POINTS_COUNT_KEY] += 1
        if self.autojoin:
            poutres = [Poutre(self.nodes[i], self.nodes[i+1]) for i in range(len(self.nodes)-1)]
            if self.closed:
                poutres.append(Poutre(self.nodes[-1], self.nodes[0]))
            parentfile.bind(Poutres(*poutres))
    def build_quadtree(self):
        from scipy.spatial import cKDTree
        coords = np.array([point.loc for point in self.nodes])
        self.tree = cKDTree(coords)
    def query_pnt(self, point, k=1):
        if not hasattr(self, 'tree'):
            self.build_quadtree()
        dist, idx = self.tree.query(point.loc, k=k)
        return [self.nodes[i] for i in idx], dist
    def __iter__(self):
        return iter(self.nodes)
    def write(self, file):
        file.write(f"{self.name} ( {len(self.nodes):d} )\n")
        for i in range(len(self.nodes)):
            node = self.nodes[i]
            file.write(node.frmt()) #  4  1.00000000000E+000  1.00000000000E+000  0.00000000000E+000
        
class Poutres(Block):
    def __init__(self, *args):
        self.name = "poutres"
        self.poutres: list[Poutre] = [i for i in args if isinstance(i, Poutre)]
    def _setup(self, parentfile):
        self.parentfile = parentfile
        for poutre in self.poutres:
            poutre.idx = parentfile.meta[cst.POUTRES_COUNT_KEY] + 1
            parentfile.meta[cst.POUTRES_COUNT_KEY] += 1
    def __iter__(self):
        return iter(self.poutres)
    def write(self, file):
        for poutre in self.poutres:
            # Example: 1 RIRI     7    1  1.00000000000E+000  0.00000000000E+000  0.00000000000E+000 11 11
            file.write(poutre.frmt())

class Liaisons(Block):
    class Elastique:
        def __init__(self, param1:list, stiffness, point ):
            self.name = "liaisons"
            self.point = Noeud(param1)
            self.type = param1
            self.stiffness = stiffness
            self.point = point
        def __repr__(self):
            return f"Elastique({self.point}, {self.type}, {self.stiffness:.11E})"
        def format(self):
            return f"elastique {self.type} {self.stiffness:.11E} {self.point.idx}\n"
    class Imposer:
        def __init__(self, param1:list, valeur, point):
            self.point = Noeud(param1)
            self.type = param1
            self.valeur = valeur
            self.point = point
        def __repr__(self):
            return f"Imposer({self.point}, {self.type}, {self.valeur:.11E})"
        def format(self):
            return f"imposer {self.type} {self.valeur:.11E} {self.point.idx}\n"
    class Encastremement:
        def __init__(self, point):
            self.point = Noeud(point)
        def __repr__(self):
            return f"Encastremement({self.point})"
        def format(self):
            return f"encastrement {self.point.idx}\n"
    class Rotule:
        def __init__(self, point):
            self.point = Noeud(point)
        def __repr__(self):
            return f"Rotule({self.point})"
        def format(self):
            return f"rotule {self.point.idx}\n"
    
    def __init__(self, *args):
        self.name = "liaisons"
        self.liaisons = [i for i in args if isinstance(i, Noeud)]
    def _setup(self, parentfile):
        self.parentfile = parentfile
    def write(self, file):
        for liaison in self.liaisons:
            if isinstance(liaison, Liaisons.Elastique):
                file.write(liaison.format())
            elif isinstance(liaison, Liaisons.Imposer):
                file.write(liaison.format())
            elif isinstance(liaison, Liaisons.Encastremement):
                file.write(liaison.format())
            elif isinstance(liaison, Liaisons.Rotule):
                file.write(liaison.format())

class Composant:
    def __init__(self):
        self.meta = {}
    def __repr__(self):
        return f"Composant({self.meta})"
    def _setup(self, parentfile):
        self.parentfile = parentfile
    def bind(self, section):
        return self.parentfile.bind(section)
    def resolve(self):
        """This method should be overridden by subclasses.
        It computes the necessary blocks to be added to the parent file.
        """
        raise NotImplementedError("Subclasses should implement this method.")