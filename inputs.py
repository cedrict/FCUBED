###############################################################################
#
#  FFFF  CCCC  U   U  BBB   EEEE  DDD      C.Thieulot
#  F     C     U   U  B  B  E     D  D     F.Gueydan
#  FFF   C     U   U  BBB   EEEE  D  D     A.Lemaitre
#  F     C     U   U  B  B  E     D  D
#  F     CCCC  UUUUU  BBB   EEEE  DDD
#
###############################################################################

from  constants_and_tools import *

avrg=3

# -1: pure shear
# -2: simple shear
# -3: solvi
# -4: poiseuille nl

# 1: clast

experiment=-3

###################################################################################################

if experiment==-1:
   from inputs_m1 import *

if experiment==-2:
   from inputs_m2 import *

if experiment==-3:
   from inputs_m3 import *

if experiment==1:
   from inputs_p1 import *

###################################################################################################
