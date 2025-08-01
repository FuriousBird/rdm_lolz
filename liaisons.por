
   RDM - Ossatures
   Calcul des Structures par la M�thode des �l�ments Finis

   Version  - 7.04 - 5 mars 2020

   Utilisateur : elian

$debut du fichier
$version
7.04
$SI unites
$nom du fichier
contraintes.por
$date
29/7/2025
$heure
20/31/13
$ossature
spatiale
$symetrie
xys  0.000000E+000 5
xzs  4.000000E+000 4
yzs  0.000000E+000 4
////
$noeuds ( 7 )
   1  0.00000000000E+000  1.00000000000E+000  0.00000000000E+000
   2  0.00000000000E+000  2.00000000000E+000  0.00000000000E+000
   3  0.00000000000E+000  3.00000000000E+000  0.00000000000E+000
   4  0.00000000000E+000  4.00000000000E+000  0.00000000000E+000
   5  0.00000000000E+000  5.00000000000E+000  0.00000000000E+000
   6  0.00000000000E+000  6.00000000000E+000  0.00000000000E+000
   7  0.00000000000E+000  0.00000000000E+000  0.00000000000E+000
   0
$poutres ( 6 )
   1 RIRI     7    1  1.00000000000E+000  0.00000000000E+000  0.00000000000E+000 11 11
   2 RIRI     2    1 -1.00000000000E+000  0.00000000000E+000  0.00000000000E+000 11 11
   3 RIRI     2    3  1.00000000000E+000  0.00000000000E+000  0.00000000000E+000 11 11
   4 RIRI     3    4  1.00000000000E+000  0.00000000000E+000  0.00000000000E+000 11 11
   5 RIRI     4    5  1.00000000000E+000  0.00000000000E+000  0.00000000000E+000 11 11
   6 RIRI     5    6  1.00000000000E+000  0.00000000000E+000  0.00000000000E+000 11 11
   0
$sections
11
TYPE QUELCONQUE
NOM *
DESIGNATION *
LOGO 0
AIRE  0.00000000000E+000
IYY  0.00000000000E+000
IZZ  0.00000000000E+000
alpha  0.00000000000E+000
WPY  0.00000000000E+000
WPZ  0.00000000000E+000
TORSION  0.00000000000E+000
KYY  1.0000000
KZZ  1.0000000
IWW  0.00000000000E+000
YCISAILLEMENT  0.00000000000E+000
ZCISAILLEMENT  0.00000000000E+000
BTY  0.00000000000E+000
BTZ  0.00000000000E+000
BTW  0.00000000000E+000
///
0
$materiaux
11
NOM Acier
MOD  2.100E+011
POI 0.3000
MAS 7800.00
DIL  1.3000E-005
LIM  2.500E+008
///
0
$liaisons ( 14 )
elastique 1  1.000000E+000 3
elastique 2  2.000000E+000 3
elastique 3  3.000000E+000 3
elastique 4  4.000000E+000 3
elastique 5  5.000000E+000 3
elastique 6  6.000000E+000 3
imposer 1  1.000000E-001 2
imposer 2  2.000000E-001 2
imposer 3  3.000000E-001 2
imposer 4  1.500000E-001 2
imposer 5  2.500000E-001 2
imposer 6  3.500000E-001 2
encastrement 7
encastrement 1
///
$gpesanteur
10.000
$modes propres
nombre 1
methode sous_espace
precision 1.00000E-002
decalage_spectral 0.00000E+000
////
$maillage
20
$fin du fichier
