###############################################################################
#
#  FFFF  CCCC  U   U  BBB   EEEE  DDD      C.Thieulot
#  F     C     U   U  B  B  E     D  D     F.Gueydan
#  FFF   C     U   U  BBB   EEEE  D  D     A.Lemaitre
#  F     C     U   U  B  B  E     D  D
#  F     CCCC  UUUUU  BBB   EEEE  DDD
#
###############################################################################

from inputs import *
import numpy as np
import random
import jatten

def make_clast(nmarker,swarm_x,swarm_y,swarm_mat,swarm_plastic_strainxx,\
               swarm_plastic_strainyy,swarm_plastic_strainxy):
   
   if experiment==-1 or experiment==-2 or experiment==-3 or experiment==-4:
      swarm_mat[:]=1

   if experiment==-5: #viscoplastic-block
      #1: block
      #2: air
      #3: inclusion
      swarm_mat[:]=2

      for im in range(0,nmarker):
          if swarm_y[im]>25e3 and swarm_y[im]<75e3: 
             swarm_mat[im]=1

      for im in range(0,nmarker):
          if abs(swarm_x[im]-Lx/2)<12.5e3/2 and abs(swarm_y[im]-Ly/2)<12.5e3/2:
             swarm_mat[im]=3


   if experiment==1: # clast

      for im in range (0,nmarker):
          swarm_mat[im]=1
          xxi=swarm_x[im]-x_inclusion
          yyi=swarm_y[im]-y_inclusion
          xxxi=xxi*np.cos(-angle_inclusion)-yyi*np.sin(-angle_inclusion)
          yyyi=xxi*np.sin(-angle_inclusion)+yyi*np.cos(-angle_inclusion)

          if not roughness:
             if xxxi**2/a_inclusion**2+yyyi**2/b_inclusion**2<1: # ellipse
                swarm_mat[im]=2
          else:
             A1=0.05
             A2=0.08
             A3=0.01
             k1=3.33
             k2=7.77
             k3=11.11
             thetai=np.arctan(yyyi/xxxi)
             xi1=random.uniform(0.5,1)
             xi2=random.uniform(0.5,1)
             xi3=random.uniform(0.5,1)
             if np.sqrt(xxxi**2/a_inclusion**2+yyyi**2/b_inclusion**2)\
                <1+xi1*A1*np.sin(k1*thetai)+xi2*A2*np.cos(k2*thetai)+xi3*A3*np.cos(k3*thetai):
                swarm_mat[im]=2

      # seeds

      #for iseed in range(0,nseed):
      #    rs=random.uniform(0,1)*a_inclusion
      #    ts=random.uniform(0,1)*2*np.pi
      #    xs=x_inclusion+rs*np.cos(ts)  
      #    ys=y_inclusion+rs*np.sin(ts)
      #for iseed in range(0,nseed):
      #    rs=random.uniform(-1,1)*a_inclusion
      #    ts=random.uniform(-1,1)*a_inclusion
      #    xs=x_inclusion+rs
      #    ys=y_inclusion+ts

      kpoisson=30
      nseed_wish=2*nseed
      avrgdist=np.sqrt((a_inclusion*2)**2/nseed_wish)/1.125
      nseedP,xseed,yseed = jatten.PoissonDisc(kpoisson,avrgdist,a_inclusion*2,a_inclusion*2)
      xseed=np.array(xseed).astype(np.float64)
      xseed[:]+=x_inclusion-a_inclusion
      yseed=np.array(yseed).astype(np.float64)
      yseed[:]+=y_inclusion-a_inclusion

      counter=0
      for iseed in range(0,nseedP):
          if np.sqrt((xseed[iseed]-x_inclusion)**2+(yseed[iseed]-y_inclusion)**2)<a_inclusion-wseed:
             counter+=1 
             for im in range(0,nmarker):
                 if swarm_mat[im]==2:
                    dist=np.sqrt((swarm_x[im]-xseed[iseed])**2+(swarm_y[im]-yseed[iseed])**2)
                    if dist<wseed:
                       swarm_plastic_strainxx[im]+=aseed*(np.cos(dist*np.pi/wseed)+1)/2 /np.sqrt(2)
                       swarm_plastic_strainxy[im]+=aseed*(np.cos(dist*np.pi/wseed)+1)/2 /np.sqrt(2)
                       swarm_plastic_strainyy[im]+=aseed*(np.cos(dist*np.pi/wseed)+1)/2 /np.sqrt(2)
                    #end if
                 #end if
             #end for
          #end if
      #end for
       
      print('real nseed in clast:',counter)

      counter=0
      counter2=0
      for im in range(0,nmarker):
          if swarm_mat[im]==2:
             counter+=1
             if swarm_plastic_strainxx[im]>0:
                counter2+=1
          #end if
      #end for
      print('nmarker in clast:',counter)
      print('nmarker in seeds:',counter2)

###############################################################################
