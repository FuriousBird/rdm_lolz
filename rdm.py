import numpy as np
from utils.porfile import Composant, Noeuds, Liaisons, Noeud, Poutre, Poutres, Empty, Fichier

edgeL = 1e-3
N = 5
A,B = N, N

# lines are x, columns are y
def matches(i,j, A,B):
    return not((i==0 and j==0) or (i==A-1 and j==B-1) or (i==0 and j==B-1) or (i==A-1 and j==0))

points = []
pointidx = {}
for i in range(A):
    for j in range(B):
        location = [i * edgeL, j * edgeL]
        if matches(i,j,A,B):
            points.append(Noeud(location))
            pointidx[(i, j)] = len(points) - 1

class Ellipsoid(Composant):
    def __init__(self, center, Yvec, Zvec):
        self.meta = {
            "center": center,
            "Yvec": Yvec,
            "Zvec": Zvec
        }
    def resolve(self, parentfile):
        #this computes the necessary blocks to be added to the parent file


# edges = []
# for i in range(A):
#     for j in range(B):
#         if (i, j) in pointidx:
#             idx = pointidx[(i, j)]
#             if i < A - 1 and (i + 1, j) in pointidx:
#                 edges.append([pointidx[(i, j)], pointidx[(i + 1, j)]])
#             if j < B - 1 and (i, j + 1) in pointidx:
#                 edges.append([pointidx[(i, j)], pointidx[(i, j + 1)]])


fichier = Fichier(
    Noeuds(*points, autojoin=True, closed=True),
    # Poutres(*[Poutre(*edge) for edge in edges]),
)

