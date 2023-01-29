###############################################################################
#
#  FFFF  CCCC  U   U  BBB   EEEE  DDD      C.Thieulot
#  F     C     U   U  B  B  E     D  D     F.Gueydan
#  FFF   C     U   U  BBB   EEEE  D  D     A.Lemaitre
#  F     C     U   U  B  B  E     D  D
#  F     CCCC  UUUUU  BBB   EEEE  DDD
#
###############################################################################

from constants_and_tools import *
from inputs import *

#@jit(nopython=True)
def gravity_model(x,y):

    #---------------------------
    if experiment<0: #benchmarks

       if experiment==-1 or experiment==-2: # shear
          gx=0
          gy=0

       if experiment==-3: #solvi
          gx=0
          gy=0

       if experiment==-4: #solkz
          gx=0
          gy=1

       if experiment==-5: #viscoplastic block
          gx=0
          gy=0

    #--------------------------
    elif experiment==1: # clast

       gx=0
       gy=0

    else:
       exit("experiment unknown in material model")

    return gx,gy 

###############################################################################
