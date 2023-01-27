###############################################################################
#
#  FFFF  CCCC  U   U  BBB   EEEE  DDD      C.Thieulot
#  F     C     U   U  B  B  E     D  D     F.Gueydan
#  FFF   C     U   U  BBB   EEEE  D  D     A.Lemaitre
#  F     C     U   U  B  B  E     D  D
#  F     CCCC  UUUUU  BBB   EEEE  DDD
#
###############################################################################
# if istep==0 create file, otherwise append data to it

import numpy as np

def write_history(output_folder,total_time,istep,u,v,ee,\
                  swarm_mat,\
                  swarm_total_strain_eff,\
                  swarm_plastic_strain_eff,\
                  swarm_sw_level,\
                  swarm_tau_eff,\
                  swarm_plastic_strain_eff0):

    if istep==0:
       histfile=open(output_folder+'hist.ascii',"w")
       histfile.write("# col  1: istep       \n")
       histfile.write("# col  2: total_time  \n")
       histfile.write("# col  3: min(u)      \n")
       histfile.write("# col  4: max(u)      \n")
       histfile.write("# col  5: min(v)      \n")
       histfile.write("# col  6: max(v)      \n")
       histfile.write("# col  7: min(ee)     \n")
       histfile.write("# col  8: max(ee)     \n")
       histfile.write("# col  9: ratio_tau   \n")
       histfile.write("# col 10: faulting_TS \n")
       histfile.write("# col 11: faulting_PS \n")
    else:
       histfile=open(output_folder+'hist.ascii',"a")

    #------------------------------------------------------
    # faulting with total_strain_eff  

    count_inclusion = 0
    count_matrice = 0
    
    for i in range (len(swarm_mat)):
        if swarm_mat[i] == 2:
            count_inclusion+=1
            inclusion = np.where(swarm_mat == 2)
           
        else :
            count_matrice+=1
            matrice = np.where(swarm_mat == 1)
            
    tau_incl = np.mean(swarm_tau_eff[inclusion])
    tau_matr = np.mean(swarm_tau_eff[matrice])
    tau_bulk = np.mean(swarm_tau_eff)

    ratio_tau = tau_incl/tau_matr

    #------------------------------------------------------
    # faulting with total_strain_eff 

    fault_TS = 0
    nofault_TS = 0
    
    for i in range (len(swarm_total_strain_eff)):
        if swarm_mat[i] == 2 and swarm_total_strain_eff[i] > 0.9 :
            fault_TS+=1
           
        if swarm_mat[i] == 2 and swarm_total_strain_eff[i] <= 0.9 :
            nofault_TS+=1
            
    faulting_TS = (fault_TS/(fault_TS+nofault_TS))*100

    #------------------------------------------------------
    # faulting with plastic_strain_eff  

    fault_PS = 0
    nofault_PS = 0
    
    for i in range (len(swarm_plastic_strain_eff)):
        if swarm_mat[i] == 2 and swarm_plastic_strain_eff[i] >= 1 :
            fault_PS+=1
           
        if swarm_mat[i] == 2 and swarm_plastic_strain_eff[i] < 1 :
            nofault_PS+=1
            
    faulting_PS = (fault_PS/(fault_PS+nofault_PS))*100

    #------------------------------------------------------

    histfile.write("%d %e %e %e %e %e %e %e %e %e %e\n" %(\
                    istep,      \
                    total_time, \
                    np.min(u),  \
                    np.max(u),  \
                    np.min(v),  \
                    np.max(v),  \
                    np.min(ee), \
                    np.max(ee), \
                    ratio_tau,  \
                    faulting_TS,\
                    faulting_PS))
        

    histfile.close()

###############################################################################
