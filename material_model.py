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
def material_model(x,y,ee,T,imat,iter,plastic_strain_marker,pressure,pf):

    #---------------------------
    if experiment<0: #benchmarks

       if experiment==-1 or experiment==-2: # shear
          val=1
          is_plastic=False
          yield_DP=0
          strain_level=0
          rho=0

       if experiment==-3: #solvi
          if (np.sqrt(x*x+y*y) < 0.2):
             val=1e3
          else:
             val=1.
          is_plastic=False
          yield_DP=0
          strain_level=0
          rho=0

       if experiment==-4: #solkz
          val=np.exp(13.8155*y)
          is_plastic=False
          yield_DP=0
          strain_level=0
          rho=np.sin(2*y)*np.cos(3*np.pi*x)

       if experiment==-5: #visco-plastic block

          val=1e17 #both air and inclusion
          is_plastic=False
          yield_DP=0
          strain_level=0
          rho=2700

          if imat==1:
             phi=37/180*np.pi
             c=1e8
             #pressure=0
             sr=max(1e-19,ee) # minimum strain rate
             yield_DP=max(0,pressure*np.sin(phi)+c*np.cos(phi))
             val=yield_DP/2/sr
             is_plastic=True
             rho=2700
          
          val=min(val,1e23)
          val=max(val,1e17)

    #--------------------------
    elif experiment==1: # clast

       sr=max(1e-19,ee) # minimum strain rate

       eta_dsl=A_values[imat-1]**(-1/n_values[imat-1])*sr**(1/n_values[imat-1]-1)*\
               np.exp(Q_values[imat-1]/Rgas/T/n_values[imat-1])*MPa

       #strain weakening
       if plastic_strain_marker <eps1:
          phi_sw=phi_values[imat-1]
          c_sw=cohesion_values[imat-1]
          strain_level=0
       elif plastic_strain_marker<eps2:
          phi_sw=phi_values[imat-1]*(weakening_factor_phi-1)/(eps2-eps1)*(plastic_strain_marker-eps1)+phi_values[imat-1] 
          c_sw=cohesion_values[imat-1]*(weakening_factor_cohesion-1)/(eps2-eps1)*(plastic_strain_marker-eps1)+cohesion_values[imat-1] 
          strain_level=(plastic_strain_marker-eps1)/(eps2-eps1)
       else:
          phi_sw=phi_values[imat-1]*weakening_factor_phi
          c_sw=cohesion_values[imat-1]*weakening_factor_cohesion
          strain_level=1

       yield_DP=background_pressure*(1-pf_coefficient)*np.sin(phi_sw)+c_sw*np.cos(phi_sw)

       eta_pl=yield_DP/2/sr
       val=min(eta_pl,eta_dsl)
       if eta_pl<eta_dsl: 
          is_plastic=1 
       else:
          is_plastic=0 

       #viscosity cutoffs
       val=min(val,1e26)
       val=max(val,1e18)

       rho=0

    else:
       exit("experiment unknown in material model")

    return val,is_plastic,yield_DP,strain_level,rho

###############################################################################
