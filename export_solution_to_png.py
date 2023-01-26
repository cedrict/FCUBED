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

def export_solution_to_png(Lx,Ly,nnx,nny,nelx,nely,x,y,u,v,q,exx,eyy,exy,ee,eta,\
                           Pf,K,phi,H,output_folder,istep):

   u_temp=np.reshape(u,(nny,nnx))
   v_temp=np.reshape(v,(nny,nnx))
   p_temp=np.reshape(q,(nny,nnx))
   exx_temp=np.reshape(exx,(nny,nnx))
   eyy_temp=np.reshape(eyy,(nny,nnx))
   exy_temp=np.reshape(exy,(nny,nnx))
   ee_temp=np.reshape((ee),(nny,nnx))
   eta_temp=np.reshape(np.log10(eta),(nely,nelx))
   Pf_temp=np.reshape(Pf,(nny,nnx))
   phi_temp=np.reshape(phi,(nely,nelx))
   K_temp=np.reshape(K,(nely,nelx))
   H_temp=np.reshape(H,(nely,nelx))

   fig,axes = plt.subplots(nrows=4,ncols=3,figsize=(18,18))

   extent=(np.amin(x),np.amax(x),np.amin(y),np.amax(y))

   onePlot(axes,fig,u_temp,   0, 0, "$v_x$",                  "x", "y", extent,  Lx,  Ly, 'Spectral_r')
   onePlot(axes,fig,v_temp,   0, 1, "$v_y$",                  "x", "y", extent,  Lx,  Ly, 'Spectral_r')
   onePlot(axes,fig,p_temp,   0, 2, "$p$",                    "x", "y", extent,  Lx,  Ly, 'RdGy_r')
   onePlot(axes,fig,exx_temp, 1, 0, "$\dot{\epsilon}_{xx}$",  "x", "y", extent,  Lx,  Ly, 'viridis')
   onePlot(axes,fig,eyy_temp, 1, 1, "$\dot{\epsilon}_{yy}$",  "x", "y", extent,  Lx,  Ly, 'Spectral_r')
   onePlot(axes,fig,exy_temp, 1, 2, "$\dot{\epsilon}_{xy}$",  "x", "y", extent,  Lx,  Ly, 'RdGy_r')
   onePlot(axes,fig,ee_temp,  2, 0, "$\dot{\epsilon}_{eff}$", "x", "y", extent,  Lx,  Ly, 'Spectral_r')
   onePlot(axes,fig,eta_temp, 2, 1, "$\eta$",                 "x", "y", extent,  Lx,  Ly, 'Spectral_r')
   onePlot(axes,fig,Pf_temp,  2, 2, "$P_f$",                  "x", "y", extent,  Lx,  Ly, 'RdGy_r')
   onePlot(axes,fig,K_temp,   3, 0, "$K$",                    "x", "y", extent,  Lx,  Ly, 'Spectral_r')
   onePlot(axes,fig,phi_temp, 3, 1, "$\phi$",                 "x", "y", extent,  Lx,  Ly, 'Spectral_r')
   onePlot(axes,fig,H_temp,   3, 2, "$H",                     "x", "y", extent,  Lx,  Ly, 'RdGy_r')

   plt.subplots_adjust(hspace=0.5)

   filename = output_folder+'img_solution_{:04d}.png'.format(istep)
   plt.savefig(filename, bbox_inches='tight')
   return 

###############################################################################
