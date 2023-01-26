from constants_and_tools import *

#geometry
Lx=2e-2 # horizontal extent of the domain in m 
Ly=1e-2 # vertical extent of the domain in m
nelx = 30             #number of elements in horizontal direction
nely = int(nelx*Ly/Lx) #number of elements in vertical direction

#clast
x_inclusion=Lx/2
y_inclusion=Ly/2
a_inclusion=Ly/3.8
b_inclusion=a_inclusion
angle_inclusion=0
roughness=False

#step velocity for bc
t1=100*year
t2=200*year
velofact=1

niter=2
nstep=1
CFL_nb=0.25
nmarker_per_dim=7
eta_ref=1e20
tol=1e-8
tfinal=15*Myr
every_png=1
every_vtu=1

background_strainrate=1e-15
v_ref=background_strainrate*Ly # bc velocity so that shear strain rate is 10^-15

#ductile rheology
matrix='Ab_GSI' 
if matrix=='An_GSI':
      Am=398.1 ; Qm=356e3; nm=3
if matrix=='An_GSS':
      Am=xx ; Qm=xx; nm=xx
if matrix=='Ab_GSI':
      Am=2511.886 ; Qm=332e3; nm=3
if matrix=='Ab_GSS':
      grainsize=100 #microns
      Am=7943.2823*grainsize**-3 ; Qm=193.1e3; nm=1

inclusion='DiW_GSI' 
if inclusion == 'DiW_GSI':
      Ai=6.80220477 ; Qi = 534e3 ; ni= 3

A_values=[Am,Ai] 
Q_values=[Qm,Qi] 
n_values=[nm,ni] 

#plastic rheology
cohesion_values=[20e6,20e6] 
phi_values=[30/180*np.pi,30/180*np.pi] 
eps1=0
eps2=1
weakening_factor_phi=0.1
weakening_factor_cohesion=0.1

depth_km=50
depth=depth_km*km
background_pressure=3000*9.81*depth 
background_temperature_C=700
background_temperature=background_temperature_C+273
pf_coefficient=0.7

nseed=20
wseed=a_inclusion/8 #width
aseed=0.75          #amplitude

#fluids
use_fluid=True
beta=1e-10
eta_fluid=1.33e-4
p_ref=background_pressure*pf_coefficient
phi0=0.01
phi_max=0.1
K0=1e-26
Hcoeff=5e-8
H0=0
Hmax=1e-5 #grossse louche 

#--------------------

path='.'

if not use_fluid:
   if velofact==1: #no velocity jump
         output_folder=path+'/S-T='+str(background_temperature_C)+\
                              '-P='+str(depth_km)+\
                              '-pf='+str(pf_coefficient)+'/'
   else:
         output_folder=path+'/SV-T='+str(background_temperature_C)+\
                               '-P='+str(depth_km)+\
                              '-pf='+str(pf_coefficient)+'/'
else:
   if velofact==1: #no velocity jump
         output_folder=path+'/F-T='+str(background_temperature_C)+\
                              '-P='+str(depth_km)+\
                             '-pf='+str(pf_coefficient)+'/'
   else:
         output_folder=path+'/FV-T='+str(background_temperature_C)+\
                               '-P='+str(depth_km)+\
                              '-pf='+str(pf_coefficient)+'/'

