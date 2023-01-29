###############################################################################
#
#  FFFF  CCCC  U   U  BBB   EEEE  DDD      C.Thieulot
#  F     C     U   U  B  B  E     D  D     F.Gueydan
#  FFF   C     U   U  BBB   EEEE  D  D     A.Lemaitre
#  F     C     U   U  B  B  E     D  D
#  F     CCCC  UUUUU  BBB   EEEE  DDD
#
###############################################################################

import numpy as np

def analytical_solution(x,y,experiment):

    #-----------------------------
    if experiment==-1: #pure shear
       vx=1-x
       vy=-0.5+y
       p=0

    #-----------------------------
    if experiment==-2: #simple shear
       vx=2*(y-1/2)
       vy=0
       p=0

    #-----------------------------
    if experiment==-3: #SolVi
       min_eta = 1.
       max_eta = 1.e3
       epsilon = 1.
       A=min_eta*(max_eta-min_eta)/(max_eta+min_eta)
       r_inclusion=0.2
       r2_inclusion=r_inclusion*r_inclusion
       r2=x*x+y*y
       # phi, psi, dphi are complex
       z=x+y*1j
       if r2<r2_inclusion:
          phi=0+0.*1j
          dphi=0+0.*1j
          psi=-4*epsilon*(max_eta*min_eta/(min_eta+max_eta))*z
          visc=1e3
       else:
          phi=-2*epsilon*A*r2_inclusion/z
          dphi=-phi/z
          psi=-2*epsilon*(min_eta*z+A*r2_inclusion*r2_inclusion/(z*z*z))
          visc=1.

       v = (phi-z*np.conjugate(dphi)-np.conjugate(psi))/(2.*visc)
       vx=v.real
       vy=v.imag
       p=-2*epsilon*dphi.real

    #-----------------------------
    if experiment==-4: #SolKz

       import solkz
       vx,vy,p=solkz.SolKzSolution(x,y)

    #-----------------------------
    if experiment==-5: #viscoplastic block
       vx=0
       vy=0
       p=0

    else:
       exit("experiment unknown in analytical_solution")

    return vx,vy,p

###############################################################################
