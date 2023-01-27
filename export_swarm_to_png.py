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
import matplotlib.pyplot as plt

def export_swarm_to_png(Lx,Ly,swarm_x,swarm_y,\
                              swarm_eta,\
                              swarm_mat,\
                              swarm_ee,\
                              swarm_total_strain_eff,\
                              swarm_tau_eff,\
                              swarm_plastic_strain_eff,\
                              swarm_plastic_strain_eff0,\
                              output_folder,istep):

    plt.subplots(nrows = 3, ncols = 3,figsize = (18,10))
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.8, wspace=0.45, hspace= 0.3)

    #materials
    plt.subplot(3,3,1, adjustable='datalim')
    sc=plt.scatter(swarm_x,swarm_y,c=swarm_mat,s=10,cmap=plt.cm.get_cmap('coolwarm'))
    plt.colorbar(sc,anchor = (2,0.5), pad = -0.18, shrink=0.75)
    plt.ylim(0,0.010)
    plt.xlim(0,0.020)
    plt.xticks([], [])
    plt.yticks([], [])
    plt.title('materials')
    
    #effective strain rate
    plt.subplot(3,3,2, adjustable='datalim')
    sc=plt.scatter(swarm_x,swarm_y,c=np.log10(swarm_ee),s=10,cmap=plt.cm.get_cmap('jet'))
    plt.colorbar(sc,anchor = (2,0.5), pad = -0.18, shrink=0.75)
    plt.ylim(0,0.010)
    plt.xlim(0,0.020)
    plt.xticks([], [])
    plt.yticks([], [])
    plt.title('$\dot{\epsilon}_{eff}$')

    #total strain
    plt.subplot(3,3,3, adjustable='datalim')
    sc=plt.scatter(swarm_x,swarm_y,c=swarm_total_strain_eff,s=10,vmin = 0, vmax = 1,cmap=plt.cm.get_cmap('Greens'))
    plt.colorbar(sc,anchor = (2,0.5), pad = -0.18, shrink=0.75)
    plt.ylim(0,0.010)
    plt.xlim(0,0.020)
    plt.xticks([], [])
    plt.yticks([], [])
    plt.title('$\epsilon_{eff}$')

    #tau eff
    plt.subplot(3,3,4, adjustable='datalim')
    sc=plt.scatter(swarm_x,swarm_y,c=swarm_tau_eff,s=10,cmap=plt.cm.get_cmap('YlOrRd'))
    plt.colorbar(sc,anchor = (2,0.5), pad = -0.18, shrink=0.75)
    plt.ylim(0,0.010)
    plt.xlim(0,0.020)
    plt.xticks([], [])
    plt.yticks([], [])
    plt.title(chr(945)+'$_{eff}$')

    #plastic strain
    plt.subplot(3,3,5, adjustable='datalim')
    sc=plt.scatter(swarm_x,swarm_y,c=swarm_plastic_strain_eff,s=10,vmin = 0, vmax = 1.5,cmap=plt.cm.get_cmap('RdYlBu_r'))
    plt.colorbar(sc,anchor = (2,0.5), pad = -0.18, shrink=0.75)
    plt.ylim(0,0.010)
    plt.xlim(0,0.020)
    plt.xticks([], [])
    plt.yticks([], [])
    plt.title('plastic strain eff')
    
    #plastic strain 0
    plt.subplot(3,3,6, adjustable='datalim')
    sc=plt.scatter(swarm_x,swarm_y,c=swarm_plastic_strain_eff0,s=10,vmin = 0, vmax = 1.5,cmap=plt.cm.get_cmap('RdYlBu_r'))
    plt.colorbar(sc,anchor = (2,0.5), pad = -0.18, shrink=0.75)
    plt.ylim(0,0.010)
    plt.xlim(0,0.020)
    plt.xticks([], [])
    plt.yticks([], [])
    plt.title('plastic strain eff sans seed')

    #viscosity 
    plt.subplot(3,3,7, adjustable='datalim')
    sc=plt.scatter(swarm_x,swarm_y,c=np.log10(swarm_eta),s=10,cmap=plt.cm.get_cmap('gist_rainbow'))
    plt.colorbar(sc,anchor = (2,0.5), pad = -0.18, shrink=0.75)
    plt.ylim(0,0.010)
    plt.xlim(0,0.020)
    plt.xticks([], [])
    plt.yticks([], [])
    plt.title('$\eta$')

    #viscosity 
    #plt.subplot(3,3,8, adjustable='datalim')
    #sc=plt.scatter(swarm_x,swarm_y,c=np.log10(swarm_eta),s=10,cmap=plt.cm.get_cmap('gist_rainbow'))
    #plt.colorbar(sc,anchor = (2,0.5), pad = -0.18, shrink=0.75)
    #plt.ylim(0,0.010)
    #plt.xlim(0,0.020)
    #plt.xticks([], [])
    #plt.yticks([], [])
    #plt.title('$\eta$')
    
    #viscosity 
    #plt.subplot(3,3,9, adjustable='datalim')
    #sc=plt.scatter(swarm_x,swarm_y,c=np.log10(swarm_eta),s=10,cmap=plt.cm.get_cmap('gist_rainbow'))
    #plt.colorbar(sc,anchor = (2,0.5), pad = -0.18, shrink=0.75)
    #plt.ylim(0,0.010)
    #plt.xlim(0,0.020)
    #plt.xticks([], [])
    #plt.yticks([], [])
    #plt.title('$\eta$')

    filename = output_folder+'img_swarm_{:04d}.png'.format(istep)
    plt.savefig(filename,bbox_inches='tight')

###############################################################################
