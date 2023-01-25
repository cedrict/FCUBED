from inputs import *
from constants_and_tools import *
from analytical_solutions import *

def define_bc_V(Lx,Ly,NV,bc_fix,bc_val,xV,yV,experiment,total_time):

    if experiment==-1: #pure shear

       for i in range(0, NV):
           #left boundary 
           if xV[i]/Lx<eps:
              bc_fix[i*ndofV] = True ; bc_val[i*ndofV] = v_ref/Ly
           #right boundary 
           if xV[i]/Lx>1-eps:
              bc_fix[i*ndofV] = True ; bc_val[i*ndofV] = -v_ref/Ly
           #bottom boundary 
           if yV[i]/Ly<eps:
              bc_fix[i*ndofV+1] = True ; bc_val[i*ndofV+1] = -v_ref/Lx
           #top boundary 
           if yV[i]/Ly>1-eps:
              bc_fix[i*ndofV+1] = True ; bc_val[i*ndofV+1] = v_ref/Lx
       #end for

    if experiment==-2: #simple shear

       for i in range(0,NV):
           #left boundary 
           if xV[i]/Lx<eps:
              bc_fix[i*ndofV+1] = True ; bc_val[i*ndofV+1] = 0 # vy
           #right boundary 
           if xV[i]/Lx>1-eps:
              bc_fix[i*ndofV+1] = True ; bc_val[i*ndofV+1] = 0 # vy
           #bottom boundary 
           if yV[i]/Ly<eps:
              bc_fix[i*ndofV+0] = True ; bc_val[i*ndofV  ] = -v_ref # vx
              bc_fix[i*ndofV+1] = True ; bc_val[i*ndofV+1] = 0   # vy
           #top boundary 
           if yV[i]/Ly>1-eps:
              bc_fix[i*ndofV+0] = True ; bc_val[i*ndofV  ] = +v_ref # vx
              bc_fix[i*ndofV+1] = True ; bc_val[i*ndofV+1] = 0   # vy
       #end for

    if experiment==-3: #solvi

       for i in range(0, NV):
           ui,vi,pi=analytical_solution(xV[i],yV[i],experiment)
           if xV[i]<eps:
              bc_fix[i*ndofV+0]   = True ; bc_val[i*ndofV+0] = ui
              bc_fix[i*ndofV+1]   = True ; bc_val[i*ndofV+1] = vi
           if xV[i]>(Lx-eps):
              bc_fix[i*ndofV+0]   = True ; bc_val[i*ndofV+0] = ui
              bc_fix[i*ndofV+1]   = True ; bc_val[i*ndofV+1] = vi
           if yV[i]<eps:
              bc_fix[i*ndofV+0]   = True ; bc_val[i*ndofV+0] = ui
              bc_fix[i*ndofV+1]   = True ; bc_val[i*ndofV+1] = vi
           if yV[i]>(Ly-eps):
              bc_fix[i*ndofV+0]   = True ; bc_val[i*ndofV+0] = ui
              bc_fix[i*ndofV+1]   = True ; bc_val[i*ndofV+1] = vi


    elif experiment==1: # clast under simple shear

       for i in range(0, NV):
           #left boundary 
           if xV[i]/Lx<eps:
              bc_fix[i*ndofV+1] = True ; bc_val[i*ndofV+1] = 0 # vy
           #right boundary 
           if xV[i]/Lx>1-eps:
              bc_fix[i*ndofV+1] = True ; bc_val[i*ndofV+1] = 0 # vy
           #bottom boundary 
           if yV[i]/Ly<eps:
              bc_fix[i*ndofV+0] = True ; bc_val[i*ndofV  ] = -v_ref*bc_factor(t1,t2,velofact,total_time) # vx
              bc_fix[i*ndofV+1] = True ; bc_val[i*ndofV+1] = 0   # vy
           #top boundary 
           if yV[i]/Ly>1-eps:
              bc_fix[i*ndofV+0] = True ; bc_val[i*ndofV  ] = +v_ref*bc_factor(t1,t2,velofact,total_time) # vx
              bc_fix[i*ndofV+1] = True ; bc_val[i*ndofV+1] = 0   # vy
       #end for

    else:

       exit("experiment unknown in define_bc")

###############################################################################

def define_bc_Pf(Lx,Ly,NPf,bc_fix_Pf,bc_val_Pf,x,y,experiment):

    for i in range(0,NPf):
        if y[i]/Ly<eps:
           bc_fix_Pf[i]=True ; bc_val_Pf[i]=p_ref
        if y[i]/Ly>(1-eps):
           bc_fix_Pf[i]=True ; bc_val_Pf[i]=p_ref
    #end for







#------------------------------------------------------------------------------
# function which returns a scalar which then mutiplies v_ref applied on the boundaries
# arguments are x,y coordinates, the time step integer, and the real time (in seconds).
#------------------------------------------------------------------------------

def bc_factor(t1,t2,velofact,total_time):
    if total_time<t1:
       factor=1
    elif total_time<t2:
       factor=velofact
    else:
       factor=1
    return factor


