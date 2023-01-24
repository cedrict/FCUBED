from constants_and_tools import *
from inputs import *

#------------------------------------------------------------------------------

#@jit(nopython=True)
def viscosity(x,y,ee,T,imat,iter,plastic_strain_marker):

    if experiment<0:

       if experiment==-1: # shear
          val=1
          is_plastic=False
          yield_vM=0
          strain_level=0

    else:
       #val=1e20
       #return val,0,0,0

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

       yield_vM=background_pressure*(1-pf_coefficient)*np.sin(phi_sw)+c_sw*np.cos(phi_sw)

       eta_pl=yield_vM/2/ee
       val=min(eta_pl,eta_dsl)
       if eta_pl<eta_dsl: 
          is_plastic=1 
       else:
          is_plastic=0 

       #viscosity cutoffs
       val=min(val,1e26)
       val=max(val,1e18)

    return val,is_plastic,yield_vM,strain_level

