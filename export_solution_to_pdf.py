import numpy as np
import matplotlib.pyplot as plt

def onePlot(axes,fig,variable, plotX, plotY, title, labelX, labelY, extVal, limitX, limitY, colorMap):
    im = axes[plotX][plotY].imshow(np.flipud(variable),extent=extVal, cmap=colorMap, interpolation="nearest")
    axes[plotX][plotY].set_title(title,fontsize=10, y=1.01)

    if (limitX != 0.0):
       axes[plotX][plotY].set_xlim(0,limitX)

    if (limitY != 0.0):
       axes[plotX][plotY].set_ylim(0,limitY)

    axes[plotX][plotY].set_xlabel(labelX)
    axes[plotX][plotY].set_ylabel(labelY)
    fig.colorbar(im,ax=axes[plotX][plotY])
    return

#############################################################################################

def export_solution_to_pdf(Lx,Ly,nnx,nny,nelx,nely,x,y,u,v,q,exx,eyy,exy,ee,eta,pf,output_folder,istep):

   u_temp=np.reshape(u,(nny,nnx))
   v_temp=np.reshape(v,(nny,nnx))
   p_temp=np.reshape(q,(nny,nnx))
   exx_temp=np.reshape(exx,(nny,nnx))
   eyy_temp=np.reshape(eyy,(nny,nnx))
   exy_temp=np.reshape(exy,(nny,nnx))
   ee_temp=np.reshape((ee),(nny,nnx))
   eta_temp=np.reshape(np.log10(eta),(nely,nelx))
   pf_temp=np.reshape(pf,(nny,nnx))

   fig,axes = plt.subplots(nrows=3,ncols=3,figsize=(18,18))

   extent=(np.amin(x),np.amax(x),np.amin(y),np.amax(y))

   onePlot(axes,fig,u_temp,   0, 0, "$v_x$",                  "x", "y", extent,  Lx,  Ly, 'Spectral_r')
   onePlot(axes,fig,v_temp,   0, 1, "$v_y$",                  "x", "y", extent,  Lx,  Ly, 'Spectral_r')
   onePlot(axes,fig,p_temp,   0, 2, "$p$",                    "x", "y", extent,  Lx,  Ly, 'RdGy_r')
   onePlot(axes,fig,exx_temp, 1, 0, "$\dot{\epsilon}_{xx}$",  "x", "y", extent,  Lx,  Ly, 'viridis')
   onePlot(axes,fig,eyy_temp, 1, 1, "$\dot{\epsilon}_{yy}$",  "x", "y", extent,  Lx,  Ly, 'Spectral_r')
   onePlot(axes,fig,exy_temp, 1, 2, "$\dot{\epsilon}_{xy}$",  "x", "y", extent,  Lx,  Ly, 'RdGy_r')
   onePlot(axes,fig,ee_temp,  2, 0, "$\dot{\epsilon}_{eff}$", "x", "y", extent,  Lx,  Ly, 'Spectral_r')
   onePlot(axes,fig,eta_temp, 2, 1, "$\eta$",                 "x", "y", extent,  Lx,  Ly, 'Spectral_r')
   onePlot(axes,fig,pf_temp,  2, 2, "$p_f$",                  "x", "y", extent,  Lx,  Ly, 'RdGy_r')

   plt.subplots_adjust(hspace=0.5)

   filename = output_folder+'img_solution_{:04d}.png'.format(istep)
   plt.savefig(filename, bbox_inches='tight')
   #plt.show()
   return 

#############################################################################################

def export_swarm_to_pdf(Lx,Ly,swarm_x,\
                              swarm_y,\
                              swarm_eta,\
                              swarm_mat,\
                              swarm_ee,\
                              swarm_total_strain_eff,\
                              swarm_tau_eff,\
                              swarm_plastic_strain_eff,\
                              output_folder,istep):


    cm = plt.cm.get_cmap('RdYlBu')

    fig = plt.figure(figsize=(18,18))

    #effective strain rate
    ax1 = fig.add_subplot(3,3,1, adjustable='datalim')
    sc=ax1.scatter(swarm_x,swarm_y,c=np.log10(swarm_ee),s=10,cmap=cm)
    plt.colorbar(sc)
    plt.title('$\dot{\epsilon}_{eff}$')

    #total strain
    ax2 = fig.add_subplot(3,3,2, adjustable='datalim')
    sc=ax2.scatter(swarm_x,swarm_y,c=swarm_total_strain_eff,s=10,cmap=cm)
    plt.colorbar(sc)
    plt.title('total strain eff')

    #tau eff
    ax3 = fig.add_subplot(3,3,3, adjustable='datalim')
    sc=ax3.scatter(swarm_x,swarm_y,c=swarm_tau_eff,s=10,cmap=cm)
    plt.colorbar(sc)
    plt.title('$\tau_{eff}$')

    #plastic strain
    ax4 = fig.add_subplot(3,3,4, adjustable='datalim')
    sc=ax4.scatter(swarm_x,swarm_y,c=swarm_plastic_strain_eff,s=10,cmap=cm)
    plt.colorbar(sc)
    plt.title('plastic strain eff')

    #viscosity 
    ax5 = fig.add_subplot(3,3,5, adjustable='datalim')
    sc=ax5.scatter(swarm_x,swarm_y,c=np.log10(swarm_eta),s=10,cmap=cm)
    plt.colorbar(sc)
    plt.title('$\eta$')

    #materials
    ax6 = fig.add_subplot(3,3,6)
    sc=ax6.scatter(swarm_x,swarm_y,c=swarm_mat,s=10,cmap=cm)
    plt.colorbar(sc)
    plt.title('material')

    #plt.subplot(3,3,9)
    #sc=plt.scatter(swarm_x,swarm_y,c=swarm_mat,s=10,cmap=cm)
    #plt.colorbar(sc)
    #plt.title('material')

    filename = output_folder+'img_swarm_{:04d}.png'.format(istep)
    plt.savefig(filename,bbox_inches='tight')


