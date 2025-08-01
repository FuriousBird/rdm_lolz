import numpy as np
import utils.constants as cst

class Point:
    """
    A class representing a point in 3D space.
    It does not handle any file saving/formatting
    """
    def __init__(self, loc):
        self.loc = np.array(loc, dtype=float)
    def __add__(self, other):
        return Point(self.loc + other.loc)
    def __repr__(self):
        return f"Point({self.loc})"
    def __sub__(self, other):
        return Point(self.loc - other.loc)
    def __mul__(self, scalar):
        return Point(self.loc * scalar)
    def __truediv__(self, scalar):
        return Point(self.loc / scalar)
    def close(self, other, threshold=1e-6):
        return self.distance(other) < threshold
    def distance(self, other):
        return np.linalg.norm(self.loc - other.loc)
    def __xor__(self, other):
        return Point(np.linalg.cross(self.loc, other.loc))
    def norm(self):
        return np.linalg.norm(self.loc)
    def normalize(self):
        norm = self.norm()
        if norm == 0:
            raise ValueError("Cannot normalize a zero vector")
        return Point(self.loc / norm)

class Poutre:
    """A class representing a beam defined by two points in 3D space.
    It does not handle any file saving/formatting
    """
    def __init__(self, point1, point2):
        self.point1 = point1
        self.point2 = point2
    def __iter__(self):
        yield self.point1
        yield self.point2
    def __repr__(self):
        return f"Poutre({self.point1}, {self.point2})"
    def length(self):
        return self.point1.distance(self.point2)
    def vecR(self, other=Point((0, 0, 1))):
        return ((self.point2-self.point1)^other).normalize()

class Block:
    def __init__(self, name):
        self.name = name
    #iteration
    def __iter__(self):
        yield self.name

class Literal(Block):
    def __init__(self, text):
        super().__init__(text)
        self.text = text
    def 

class Fichier:
    def __init__(self, *args):
        self.sections = [i for i in args if isinstance(i, Block)]
        self.meta = {}
        for sect in self.sections:
            sect.parentfile = self
            if hasattr(sect, '_setup'):
                sect._setup(self)
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
        file.write(cst.DEFAULT_FILE_HEADER)
        for section in self.sections:
            for line in section:
                file.write(line + "\n")
        file.write(cst.FILE_END)

class Points(Block):
    def __init__(self, *args):
        self.name="nodes"
        self.points = [i for i in args if isinstance(i, Point)]
    def _setup(self, parentfile):
        parentfile.meta[cst.POINTS_COUNT_KEY] = parentfile.meta.get(cst.POINTS_COUNT_KEY, 0)
        self.parentfile = parentfile
        for point in self.points:
            point.idx = parentfile.meta[cst.POINTS_COUNT_KEY]
            parentfile.meta[cst.POINTS_COUNT_KEY] += 1
    def build_quadtree(self):
        from scipy.spatial import cKDTree
        coords = np.array([point.loc for point in self.points])
        self.tree = cKDTree(coords)
    def query_pnt(self, point, k=1):
        if not hasattr(self, 'tree'):
            self.build_quadtree()
        dist, idx = self.tree.query(point.loc, k=k)
        return [self.points[i] for i in idx], dist
    def __iter__(self):
        return iter(self.points)
    def write(self, file):
        file.write(f"{self.name} ( {len(self.points):d} )\n")
        for i in range(len(self.points)):
            point = self.points[i]
            file.write(f"  {i+1} {point.loc[0]:.11E} {point.loc[1]:.11E} {point.loc[2]:.11E}\n") #  4  1.00000000000E+000  1.00000000000E+000  0.00000000000E+000
        file.write(f"  0")

class Poutres(Block):
    def __init__(self, *args):
        self.name = "poutres"
        self.poutres: list[Poutre] = [i for i in args if isinstance(i, Poutre)]
    def _setup(self, parentfile):
        parentfile.meta[cst.POUTRES_COUNT_KEY] = parentfile.meta.get(cst.POUTRES_COUNT_KEY, 0)
        self.parentfile = parentfile
        for poutre in self.poutres:
            poutre.idx = parentfile.meta[cst.POUTRES_COUNT_KEY] + 1
            parentfile.meta[cst.POUTRES_COUNT_KEY] += 1
    def __iter__(self):
        return iter(self.poutres)
    def write(self, file):
        file.write(f"${self.name} ( {len(self.poutres)} )\n")
        for i, poutre in enumerate(self.poutres):
            # Example: 1 RIRI     7    1  1.00000000000E+000  0.00000000000E+000  0.00000000000E+000 11 11
            n1 = poutre.point1.idx + 1
            n2 = poutre.point2.idx + 1
            vec = poutre.vecR().loc
            file.write(f"   {i+1} {cst.POUTRE_TYPE}     {n1:2d}    {n2:2d} {vec[0]:.11E} {vec[1]:.11E} {vec[2]:.11E} 11 11\n")
        file.write("   0\n")

class Liaisons(Block):
    class Elastique:
        def __init__(self, param1:list, stiffness, point ):
            self.name = "liaisons"
            self.point = Point(param1)
            self.type = param1
            self.stiffness = stiffness
            self.point = point
        def __repr__(self):
            return f"Elastique({self.point}, {self.type}, {self.stiffness})"
        def format(self):
            return f"elastique {self.type} {self.stiffness:.11E} {self.point.idx}\n"
    class Imposer:
        def __init__(self, param1:list, valeur, point):
            self.point = Point(param1)
            self.type = param1
            self.valeur = valeur
            self.point = point
        def __repr__(self):
            return f"Imposer({self.point}, {self.type})"
        def format(self):
            return f"imposer {self.type} {self.valeur:.11E} {self.point.idx}\n"
    class Encastremement:
        def __init__(self, point):
            self.point = Point(point)
        def __repr__(self):
            return f"Encastremement({self.point})"
        def format(self):
            return f"encastrement {self.point.idx}\n"
    class Rotule:
        def __init__(self, point):
            self.point = Point(point)
        def __repr__(self):
            return f"Rotule({self.point})"
        def format(self):
            return f"rotule {self.point.idx}\n"
    
    def __init__(self, *args):
        self.name = "liaisons"
        self.liaisons = [i for i in args if isinstance(i, Point)]
    def _setup(self, parentfile):
        pass
    def write(self, file):
        file.write(f"${self.name}\n")
        for liaison in self.liaisons:
            if isinstance(liaison, self.Elastique):
                file.write(liaison.format())
            elif isinstance(liaison, self.Imposer):
                file.write(liaison.format())
            elif isinstance(liaison, self.Encastremement):
                file.write(liaison.format())
            elif isinstance(liaison, self.Rotule):
                file.write(liaison.format())
        file.write("  0\n")