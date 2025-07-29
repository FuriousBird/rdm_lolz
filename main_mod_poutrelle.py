import numpy as np

HEAD = """
   RDM - Ossatures
   Calcul des Structures par la Méthode des Éléments Finis

   Version  - 7.04 - 30 janvier 2019

   Utilisateur : GC
  
$debut du fichier
$version
7.04
$SI unites
$nom du fichier
test.por
$date
28/7/2025
$heure
12/25/45
$ossature
plancher
"""

nodeinf = []
nodeidx = {}
N = 5
dLi,dLj = 80e-3,5*1e-2
Lj = 1200*1e-3 #m
A,B = 1,int(Lj//dLj)+2 #+1 for edge to node cnt and +1 for the rest
rest = Lj%dLj
tol_rest = 5e-3
coeff = 1e2
if rest<=tol_rest:
    B=B-1
def iscorner(i,j, A,B):
    return False
    return (i==0) * (j==0) or (i==A-1)*(j==B-1) or (i==A-1)*(j==0) or (i==0)*(j==B-1)

for i in range(A):
    for j in range(B):
        if iscorner(i,j,A,B):
            continue
        ix = len(nodeinf)
        if j==B-1 and rest>tol_rest:
            nodeinf.append([(i*dLi*coeff, ((j-1)*dLj+rest)*coeff, 0), i,j, ix])
        else:
            nodeinf.append([(i*dLi*coeff, j*dLj*coeff, 0), i,j, ix])
        nodeidx[(i,j)] = ix 

nodeout = f"$Noeuds ({len(nodeinf)})\n"
for inode, node in enumerate(nodeinf):
    nodeout += "   {:d}  {: .11E}  {: .11E}  {: .11E}\n".format(node[-1]+1, *node[0])
    
nodeout+="   0\n////\n"

def isedgeok(node1, node2,A,B):
    i1,j1 = node1
    i2,j2 = node2
    if i1>=A or i2>=A or j1>=B or j2>=B:
        return False
    return not iscorner(i1,j1,A,B) and not iscorner(i2,j2,A,B)
    
#build the grid, offet all indices by the 

def frmtedge(edge):
    return "   {:d} {}     {:d}    {:d}  {: .11E}  {: .11E}  {: .11E} {:d} {:d}".format(*edge)
    

edgeinf = []
offset = 0
tot = 1
for i in range(A):
    for j in range(B):
        
        do_edge = isedgeok((i,j), (i+1,j), A,B)

        if do_edge and False:
            ix1 = nodeidx[(i,j)]
            ix2 = nodeidx[(i+1,j)]
            vecCol = np.subtract(nodeinf[ix2][0], nodeinf[ix1][0])
            vecR = np.cross(vecCol, [0,0,1])
            vecR = vecR/np.linalg.norm(vecR)
            edge = [tot, "RIRI", ix1+1, ix2+1, *vecR, 11, 11]
            edgeinf.append(edge)
            tot += 1
        do_edge = isedgeok((i,j), (i,j+1), A,B)

        if do_edge:
            ix1 = nodeidx[(i,j)]
            ix2 = nodeidx[(i,j+1)]
            vecCol = np.subtract(nodeinf[ix2][0], nodeinf[ix1][0])
            vecR = np.cross(vecCol, [0,0,1])
            vecR = vecR/np.linalg.norm(vecR)
            edge = [tot, "RIRI", ix1+1, ix2+1, *vecR, 11, 11]
            edgeinf.append(edge)
            tot += 1

edgeout = f"$poutres ({len(edgeinf)})\n"
for iedge, edge in enumerate(edgeinf):
    edgeout += "   {:d} {} {:d} {:d} {: .11E}  {: .11E}  {: .11E} {:d} {:d}\n".format(*edge)
    
edgeout+="   0\n////\n"

def matches_link(node, A,B):
    i,j = node[1:3]
    return j==0 or j==B-1



linksinf = [node for node in nodeinf if matches_link(node,A,B)]

linkout = "$liaisons ({:d})\n".format(len(linksinf))
for ilink, link in enumerate(linksinf):
    linkout+="encastrement {:d}\n".format(link[-1]+1)
linkout+="   0\n////\n"

MID = """$sections
11
TYPE PARAMETREE
NOM *Rectangle creux
DESIGNATION *LY = 30.0 LZ = 30.0 e = 2.0 mm
LOGO 7
DIMENSIONS 3
 3.000000E-002
 3.000000E-002
 2.000000E-003
AIRE  2.24000000000E-004
IYY  2.94186666667E-008
IZZ  2.94186666667E-008
WPY  2.35598032304E-006
WPZ  2.35598033881E-006
TORSION  4.54939167005E-008
KYY  0.4328590
KZZ  0.4328585
IWW  6.54107820417E-015
///
0
$materiaux
11
NOM Aluminium
MOD  6.750E+010
POI 0.3400
MAS 2700.00
DIL  2.4000E-005
LIM  3.000E+007
///
0
"""


END = """///
$gpesanteur
10.000
$cas de charges
0
////
$modes propres
nombre 1
methode sous_espace
precision 1.00000E-002
decalage_spectral 0.00000E+000
////
$maillage
20
$fin du fichier
"""
        
out = HEAD+nodeout+edgeout+MID+linkout+END

print(out)

print(f"L = {Lj}")
print(f"{len(nodeinf):d} nodes")
print(f"{len(edgeinf):d} edges")
print(f"{len(linksinf):d} node constraints")
