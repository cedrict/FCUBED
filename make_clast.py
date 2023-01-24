from inputs import *
import numpy as np
import random
import jatten

def make_clast(nmarker,swarm_x,swarm_y,swarm_mat,swarm_plastic_strainxx,swarm_plastic_strainyy,swarm_plastic_strainxy):

   if experiment==-1:
      swarm_mat[:]=1


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


######################################################################################################################################
#######################################################################################################################################



#######################################################################################################################################
# Tests FG multi_inclusion

    # xxi=swarm_x[im]-x_inclusion-Lx/4
    # yyi=swarm_y[im]-y_inclusion+Ly/4
    # xxxi=xxi*np.cos(-angle_inclusion)-yyi*np.sin(-angle_inclusion)
    # yyyi=xxi*np.sin(-angle_inclusion)+yyi*np.cos(-angle_inclusion)
    # thetai=np.arctan(yyyi/xxxi)
    # xi1=random.uniform(0.5,1)
    # xi2=random.uniform(0.5,1)
    # xi3=random.uniform(0.5,1)
    # if np.sqrt(xxxi**2/a_inclusion**2+yyyi**2/b_inclusion**2)\
    #    <1+xi1*A1*np.sin(k1*thetai)+xi2*A2*np.cos(k2*thetai)+xi3*A3*np.cos(k3*thetai):
    #    swarm_mat[im]=2
    # xxi=swarm_x[im]-x_inclusion+Lx/4
    # yyi=swarm_y[im]-y_inclusion-Ly/4
    # xxxi=xxi*np.cos(-angle_inclusion)-yyi*np.sin(-angle_inclusion)
    # yyyi=xxi*np.sin(-angle_inclusion)+yyi*np.cos(-angle_inclusion)
    # thetai=np.arctan(yyyi/xxxi)
    # xi1=random.uniform(0.5,1)
    # xi2=random.uniform(0.5,1)
    # xi3=random.uniform(0.5,1)
    # if np.sqrt(xxxi**2/a_inclusion**2+yyyi**2/b_inclusion**2)\
    #    <1+xi1*A1*np.sin(k1*thetai)+xi2*A2*np.cos(k2*thetai)+xi3*A3*np.cos(k3*thetai):
    #    swarm_mat[im]=2

#######################################################################################################################################   






#################################################################################################
#             random seed                                                                       #
#################################################################################################
#start = time.time()

#for im in range (0,nmarker):
    
# Seed in all inclusion           
                              #
#    if Amount == 'Only' and swarm_mat[im]==2:
#        if Values_seed == 'Definite' :
#            swarm_plastic_strainxx[im]=PS_seed_calc
#            swarm_plastic_strainxy[im]=PS_seed_calc
#            swarm_plastic_strainyy[im]=PS_seed_calc
#            count_seed += 1
            
#        if Values_seed == 'Included' :
#            swarm_plastic_strainxx[im]=random.uniform(PS_seed_calc-PS_seed_variation_calc,PS_seed_calc+PS_seed_variation_calc)
#            swarm_plastic_strainxy[im]=random.uniform(PS_seed_calc-PS_seed_variation_calc,PS_seed_calc+PS_seed_variation_calc)
#            swarm_plastic_strainyy[im]=random.uniform(PS_seed_calc-PS_seed_variation_calc,PS_seed_calc+PS_seed_variation_calc)
#            count_seed += 1  
            
#        if Values_seed == 'Random' :
#            swarm_plastic_strainxx[im]=random.uniform(0,PS_seed_calc)
#            swarm_plastic_strainxy[im]=random.uniform(0,PS_seed_calc)
#            swarm_plastic_strainyy[im]=random.uniform(0,PS_seed_calc)
#            count_seed += 1
            
# Punctual seed in inclusion                                        #

#    if Amount == 'Many' and swarm_mat[im]==2 and (swarm_iel[im]%11==0 or swarm_iel[im]%17==0 or swarm_iel[im]%29==0):
#    #if Amount == 'Many' and swarm_mat[im]==2 and (swarm_iel[im]%Location==0):
       
#        if Values == 'Definite' :
#            swarm_plastic_strainxx[im]=PS_seed_calc
#            swarm_plastic_strainxy[im]=PS_seed_calc
#            swarm_plastic_strainyy[im]=PS_seed_calc
#            count_seed += 1
        
#        if Values == 'Included' :
#            swarm_plastic_strainxx[im]=random.uniform(PS_seed_calc-PS_seed_variation_calc,PS_seed_calc+PS_seed_variation_calc)
#            swarm_plastic_strainxy[im]=random.uniform(PS_seed_calc-PS_seed_variation_calc,PS_seed_calc+PS_seed_variation_calc)
#            swarm_plastic_strainyy[im]=random.uniform(PS_seed_calc-PS_seed_variation_calc,PS_seed_calc+PS_seed_variation_calc)
#            count_seed += 1
        
#        if Values == 'Random' :
#            swarm_plastic_strainxx[im]=random.uniform(0,PS_seed_calc)
#            swarm_plastic_strainxy[im]=random.uniform(0,PS_seed_calc)
#            swarm_plastic_strainyy[im]=random.uniform(0,PS_seed_calc)
#            count_seed += 1
           
          
        
        # A single seed at the core of inclusion                         #
        # element in core of inclusion : ((nel/2) -/+ (nelx/2))
        # Ajouter un +/- 1 si vous voulez que ce soit l'élément à gauche de ceux qui plot

#    if Amount == 'Single' and swarm_mat[im]==2 and (swarm_iel[im]%((nel/2)+(nelx/2))==0):  
        
#        if Values == 'Definite' :
#            swarm_plastic_strainxx[im]=PS_seed_calc
#            swarm_plastic_strainxy[im]=PS_seed_calc
#            swarm_plastic_strainyy[im]=PS_seed_calc
#            count_seed += 1
            
#        if Values == 'Included' :
#            swarm_plastic_strainxx[im]=random.uniform(PS_seed_calc-PS_seed_variation_calc,PS_seed_calc+PS_seed_variation_calc)
#            swarm_plastic_strainxy[im]=random.uniform(PS_seed_calc-PS_seed_variation_calc,PS_seed_calc+PS_seed_variation_calc)
#            swarm_plastic_strainyy[im]=random.uniform(PS_seed_calc-PS_seed_variation_calc,PS_seed_calc+PS_seed_variation_calc)
#            count_seed += 1               
            
#        if Values == 'Random' :
#            swarm_plastic_strainxx[im]=random.uniform(0,PS_seed_calc)
#            swarm_plastic_strainxy[im]=random.uniform(0,PS_seed_calc)
#            swarm_plastic_strainyy[im]=random.uniform(0,PS_seed_calc)
#            count_seed += 1
       
                  
#print('seed = ',count_seed)   
#prop_seed = (count_seed/count_inclusion)*100
#print ('CONCLUSION = ', prop_seed,'% de seed dans inclusion')
#print("paint swarm strain in inclusion: %.3f s" % (time.time() - start))

#######################################################################################################################################
#######################################################################################################################################





