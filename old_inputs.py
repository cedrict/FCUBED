from constants import *
import numpy as np
import os

"""
Test = choix du répertoire
!!! Préciser le nombre d'itération (niter) par istep 
!!! Préciser le nombre de nelx pour la résolution
Inclusion = présence ou pas de l'inclusion
Size_inclusion & Shape_inclusion = géométrie de l'inclusion
Seed = Nombre et position des défauts dans l'inclusion
Plastic_strain_seed = Valeur de la déformation plastique dans le(s) défaut(s)
!!! Controlé par V_PS_seed : random à comme maximum V_PS_seed
Roughness = présence d'une rugosité
use_plasticity = présence ou non d'une déformation plastique dans l'inclusion
Weakening = présence ou pas d'un adoucissement 
"""
                        
Size = 'M'                          # S, M or L                                # AM
Shape = 'Circle'                    # Circle or Ellipsis                       # AM
Roughness = 'No'                    # Rgh : Yes or No                          # AM

Amount = 'Many'                     # Without, Single, Many or Only            # AM
Location = 'CT'                       # Random, Definite                       # AM
Values = 'Definite'                 # No, Random, Included , Definite  
                     
use_plasticity=True                 # Plastic : True or False                                         # AM
Weakening = 'With'                  # Weak : With or Without                                          # AM

"""
Size        
Shape       
Roughness

Seed_proportion         
Seed_localisation       
Values_Seed             

use_plasticity          
Weakening               

Other
"""   
#################################################################################
# Repertory choice
#################################################################################


#################################################################################



#################################################################################
# Weakening - without weakening = 1
#################################################################################
if Weakening == 'With' :
    weakening_factor_phi=0.1
    weakening_factor_cohesion=0.1
    
if Weakening == 'Without' :
    weakening_factor_phi=1
    weakening_factor_cohesion=1
#################################################################################

#################################################################################
# medium value of plastic strain in seed
#################################################################################
Denominator_seed = 2

PS_seed = 0.75
PS_seed_calc = PS_seed/np.sqrt(2)
PS_seed_variation = 0.1*PS_seed
PS_seed_variation_calc = PS_seed_variation/np.sqrt(2)

#################################################################################

#background strainrate
strainrate=1e-15
tfinal=15*Myr

#mat 1: A_ab=398.1, n_ab=3, Q_ab=356e3
#mat 2: A_cpx=6.80220477, Q_cpx=534e3, n_cpx=5.52

# INPUT PARAMETERS
strainrate=1e-15
depthkm=40  #AM
tempdegC=510 # AM
pf=0.9 # AM
grainsize=100 #FG
velofact=1  # for increasing velocity at timestep 10 #FG

niter=250
Lx=2e-2 # horizontal extent of the domain in m 
Ly=1e-2 # vertical extent of the domain in m
nelx = 80             #number of elements in horizontal direction
nely = int(nelx*Ly/Lx) #number of elements in vertical direction
v0=strainrate*Ly # bc velocity so that shear strain rate is 10^-15

#inclusion parameters
x_inclusion=Lx/2
y_inclusion=Ly/2

################################      Geometry      ###############################
if Size == 'L' and Shape == 'Circle' :
    a_inclusion=Ly/2.5
    b_inclusion=Ly/2.5

if Size == 'L' and Shape == 'Ellipsis' :
    a_inclusion=Ly/2.5
    b_inclusion=Ly/1.3

if Size == 'M' and Shape == 'Circle' :
    a_inclusion=Ly/3.8
    b_inclusion=Ly/3.8

if Size == 'M' and Shape == 'Ellipsis' :
    a_inclusion=Ly/3.8
    b_inclusion=Ly/1.9

if Size == 'S' and Shape == 'Circle' :                                                     #
    a_inclusion=Ly/5
    b_inclusion=Ly/5

if Size == 'S' and Shape == 'Ellipsis' :
    a_inclusion=Ly/5
    b_inclusion=Ly/2.5
################################################################################

angle_inclusion=-np.pi/4

#################################################################################
#                                     Roughness                                 #
#################################################################################
#################################################################################
#################################################################################

#time stepping
nstep=100000
CFL_nb=0.5

#markers parameters
nmarker_per_dim=6
avrg=3

#nonlinear iterations parameters
tol=1e-4

##########################################################################################
#                                   Name and location                                    #
##########################################################################################

#ddir="/Users/fred/1-RECHERCHE/NUMERIQUE/SimpleShear2D-CThieulot/INCLUSION/Run-2mai2022/"
##########################################################################################

##########################################################################################    


Variation = 'Niter='+str(niter)


ddir="./"
  
if not os.path.isdir(ddir):os.mkdir(ddir)
if tempdegC == 510 and depthkm == 40 and strainrate == 1e-15 and pf == 0.9 and matrix == 'Ab_GSI' and inclusion == 'DiW_GSI' : #AM
    if Roughness == 'Yes':
        if Values == 'Included' :
            folder='RM'+'__'+str(Variation)+'__Inc='+str(Shape[0:3])+'+'+str(Size)+'+'+'Rgh'+'__'+str(Amount)+'Seed_in'+str(Location)+'='+'-'+str(PS_seed)+'-'+str(Values)+'='+str(PS_seed_variation)+'/'   #AM
        else : 
            folder='RM'+'__'+str(Variation)+'__Inc='+str(Shape[0:3])+'+'+str(Size)+'+'+'Rgh'+'__'+str(Amount)+'Seed_in'+str(Location)+'='+'-'+str(PS_seed)+'-'+str(Values)+'/'   #AM
    
    if Roughness == 'No' :
        if Values == 'Included' :
            folder='RM'+'__'+str(Variation)+'__Inc='+str(Shape[0:3])+'+'+str(Size)+'__'+str(Amount)+'Seed_in'+str(Location)+'='+'-'+str(PS_seed)+'-'+str(Values)+'='+str(PS_seed_variation)+'/'   #AM
        else : 
            folder='RM'+'__'+str(Variation)+'__Inc='+str(Shape[0:3])+'+'+str(Size)+'__'+str(Amount)+'Seed_in'+str(Location)+'='+'-'+str(PS_seed)+'-'+str(Values)+'/'   #AM
else :  #AM
    folder='M'+'__'+str(Variation)+'__Inc='+str(Shape[0:3])+'+'+str(Size)+'+'+'Rgh'+'__'+str(Amount)+'Seed_in'+str(Location)+'='+'-'+str(PS_seed)+'-'+str(Values)+'='+str(PS_seed_variation)+'/'   #AM

output_folder=ddir+folder
if not os.path.isdir(output_folder):os.mkdir(output_folder)

print(output_folder)  

###########################################################################################
###########################################################################################







count_matrice = 0
count_inclusion = 0
count_seed = 0


print ('background_pressure=',background_pressure)


eta_ref=1e20 # purely numerical parameter - do not change
