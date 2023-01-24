import numpy as np
import math as math
import sys as sys
import scipy
import scipy.sparse as sps
from scipy.sparse.linalg.dsolve import linsolve
from scipy.sparse import csr_matrix, lil_matrix 
import time as time
import os
from numpy import linalg as LA
import matplotlib.pyplot as plt

###############################################################################

from basis_functions import * 
from constants_and_tools import *
from inputs import *
from material_model import *
from define_bc import *
from export_solution_to_vtu import *
from export_swarm_to_vtu import *
from make_clast import *
from export_solution_to_pdf import *

###############################################################################

#sys.stdout = open(output_folder+'log.txt', 'w')

print("-----------------------------")
print("----------- fcubed ----------")
print("-----------------------------")

    
nnx=2*nelx+1  # number of elements, x direction
nny=2*nely+1  # number of elements, y direction

NV=nnx*nny           # number of velocity nodes
nel=nelx*nely        # number of elements, total
NP=(nelx+1)*(nely+1) # number of pressure nodes

NfemV=NV*ndofV   # number of velocity dofs
NfemP=NP*ndofP   # number of pressure dofs
Nfem=NfemV+NfemP # total number of dofs

hx=Lx/nelx # size of element in x direction
hy=Ly/nely # size of element in y direction

nmarker_per_element=nmarker_per_dim**2
nmarker=nel*nmarker_per_element

#quadrature parameters
nqperdim=3
qcoords=[-np.sqrt(3./5.),0.,np.sqrt(3./5.)]
qweights=[5./9.,8./9.,5./9.]

###############################################################################
#opening/creating files/folders 

if not os.path.isdir(output_folder):
   print('The results folder '+output_folder+' does not exist. Creating a new one..')
   print("------------------------------")
   os.mkdir(output_folder)
else:
   print('The results folder '+output_folder+' already exists!')
   print("------------------------------")

convfile=open(output_folder+'conv.ascii',"w")

###############################################################################

print("nelx",nelx)
print("nely",nely)
print("nel",nel)
print("nnx=",nnx)
print("nny=",nny)
print("NV=",NV)
print("NP=",NP)
print("NfemV=",NfemV)
print("NfemP=",NfemP)
print("nmarker_per_dim=",nmarker_per_dim)
print("nmarker=",nmarker)
print("------------------------------")

rVnodes=[-1,+1,+1,-1, 0,+1, 0,-1,0]
sVnodes=[-1,-1,+1,+1,-1,0,+1,0,0]

###############################################################################
# grid point and icon setup
###############################################################################
start = time.time()

xV=np.zeros(NV,dtype=np.float64)  # x coordinates
yV=np.zeros(NV,dtype=np.float64)  # y coordinates
iconV=np.zeros((mV,nel),dtype=np.int32)

counter = 0
for j in range(0,nny):
    for i in range(0,nnx):
        xV[counter]=i*hx/2.
        yV[counter]=j*hy/2.
        counter += 1
    #end for
#end for

counter = 0
for j in range(0,nely):
    for i in range(0,nelx):
        iconV[0,counter]=(i)*2+1+(j)*2*nnx -1
        iconV[1,counter]=(i)*2+3+(j)*2*nnx -1
        iconV[2,counter]=(i)*2+3+(j)*2*nnx+nnx*2 -1
        iconV[3,counter]=(i)*2+1+(j)*2*nnx+nnx*2 -1
        iconV[4,counter]=(i)*2+2+(j)*2*nnx -1
        iconV[5,counter]=(i)*2+3+(j)*2*nnx+nnx -1
        iconV[6,counter]=(i)*2+2+(j)*2*nnx+nnx*2 -1
        iconV[7,counter]=(i)*2+1+(j)*2*nnx+nnx -1
        iconV[8,counter]=(i)*2+2+(j)*2*nnx+nnx -1
        counter += 1
    #end for
#end for

print("velocity connectivity & nodes: %.3f s" % (time.time() - start))

#################################################################
# pressure connectivity array
#################################################################
start = time.time()

xP=np.zeros(NP,dtype=np.float64)     # x coordinates
yP=np.zeros(NP,dtype=np.float64)     # y coordinates
iconP=np.zeros((mP,nel),dtype=np.int32)

counter = 0
for j in range(0,nely):
    for i in range(0,nelx):
        iconP[0,counter]=i+j*(nelx+1)
        iconP[1,counter]=i+1+j*(nelx+1)
        iconP[2,counter]=i+1+(j+1)*(nelx+1)
        iconP[3,counter]=i+(j+1)*(nelx+1)
        counter += 1
    #end for
#end for

counter = 0
for j in range(0, nely+1):
    for i in range(0, nelx+1):
        xP[counter]=i*Lx/float(nelx)
        yP[counter]=j*Ly/float(nely)
        counter += 1

print("pressure connectivity & nodes: %.3f s" % (time.time() - start))

#################################################################
# compute area of elements
# This is a good test because it uses the quadrature points and 
# weights as well as the shape functions. If any area comes out
# negative or zero, or if the sum does not equal to the area of the 
# whole domain then there is a major problem which needs to 
# be addressed before FE are set into motion.
#################################################################
start = time.time()

area=np.zeros(nel,dtype=np.float64) 
NNNV=np.zeros(mV,dtype=np.float64)       
dNNNVdr=np.zeros(mV,dtype=np.float64)   
dNNNVds= np.zeros(mV,dtype=np.float64) 

for iel in range(0,nel):
    for iq in range(0,nqperdim):
        for jq in range(0,nqperdim):
            rq=qcoords[iq]
            sq=qcoords[jq]
            weightq=qweights[iq]*qweights[jq]
            NNNV[0:mV]=NNV(rq,sq)
            dNNNVdr[0:mV]=dNNVdr(rq,sq)
            dNNNVds[0:mV]=dNNVds(rq,sq)
            jcb=np.zeros((2,2),dtype=np.float64)
            for k in range(0,mV):
                jcb[0,0] += dNNNVdr[k]*xV[iconV[k,iel]]
                jcb[0,1] += dNNNVdr[k]*yV[iconV[k,iel]]
                jcb[1,0] += dNNNVds[k]*xV[iconV[k,iel]]
                jcb[1,1] += dNNNVds[k]*yV[iconV[k,iel]]
            jcob = np.linalg.det(jcb)
            area[iel]+=jcob*weightq
        if area[iel]<0: 
           for k in range(0,mV):
               print (xV[iconV[k,iel]],yV[iconV[k,iel]])
        #end for
    #end for
#end for

print("     -> area (m,M) %.6e %.6e " %(np.min(area),np.max(area)))
print("     -> total area meas %.8e " %(area.sum()))
print("     -> total area anal %.8e " %(Lx*Ly))

print("compute elements areas: %.3f s" % (time.time() - start))

###############################################################################
# swarm (=all the particles) setup
###############################################################################
start = time.time()

swarm_x=np.zeros(nmarker,dtype=np.float64)                  # x coordinates   
swarm_y=np.zeros(nmarker,dtype=np.float64)                  # y coordinates 
swarm_mat=np.zeros(nmarker,dtype=np.int8)                   # type of material 
swarm_paint=np.zeros(nmarker,dtype=np.float64)              # paint to show the virtual deformation of the grid
swarm_r=np.zeros(nmarker,dtype=np.float64)                  # reduced coordinates r
swarm_s=np.zeros(nmarker,dtype=np.float64)                  # reduced coordinates s
swarm_u=np.zeros(nmarker,dtype=np.float64)                  # velocity x
swarm_v=np.zeros(nmarker,dtype=np.float64)                  # velocity y
swarm_exx=np.zeros(nmarker,dtype=np.float64)                # strain rate xx
swarm_exy=np.zeros(nmarker,dtype=np.float64)                # strain rate yy
swarm_eyy=np.zeros(nmarker,dtype=np.float64)                # strain rate xy
swarm_ee=np.zeros(nmarker,dtype=np.float64)                 # effective strain rate
swarm_total_strainxx=np.zeros(nmarker,dtype=np.float64)     # total strain xx
swarm_total_strainxy=np.zeros(nmarker,dtype=np.float64)     # total strain yy
swarm_total_strainyy=np.zeros(nmarker,dtype=np.float64)     # total strain xy
swarm_total_strain_eff=np.zeros(nmarker,dtype=np.float64)   # effective strain 
swarm_plastic_strainxx=np.zeros(nmarker,dtype=np.float64)   # plastic strain xx
swarm_plastic_strainxy=np.zeros(nmarker,dtype=np.float64)   # plastic strain yy
swarm_plastic_strainyy=np.zeros(nmarker,dtype=np.float64)   # plastic strain xy
swarm_plastic_strain_eff=np.zeros(nmarker,dtype=np.float64) # effective plastic strain 
swarm_plastic_strain_eff0=np.zeros(nmarker,dtype=np.float64) # effective plastic strain initial
swarm_tauxx=np.zeros(nmarker,dtype=np.float64)              # dev stress xx
swarm_tauxy=np.zeros(nmarker,dtype=np.float64)              # dev stress yy
swarm_tauyy=np.zeros(nmarker,dtype=np.float64)              # dev stress xy
swarm_tau_eff=np.zeros(nmarker,dtype=np.float64)            # effective dev stress
swarm_iel=np.zeros(nmarker,dtype=np.int32)                  # element identity
swarm_eta=np.zeros(nmarker,dtype=np.float64)                # viscosity
swarm_p_dyn=np.zeros(nmarker,dtype=np.float64)              # pressure 
swarm_yield=np.zeros(nmarker,dtype=np.float64)              # yield value
swarm_is_plastic=np.zeros(nmarker,dtype=np.int32)           # plastic deformation active 
swarm_tau_angle=np.zeros(nmarker,dtype=np.float64)          # principal angle dev stress 
swarm_sigmaxx=np.zeros(nmarker,dtype=np.float64)            # full stress xx
swarm_sigmaxy=np.zeros(nmarker,dtype=np.float64)            # full stress yy
swarm_sigmayy=np.zeros(nmarker,dtype=np.float64)            # full stress xy
swarm_sigma_angle=np.zeros(nmarker,dtype=np.float64)        # principal angle full stress
swarm_sigma1=np.zeros(nmarker,dtype=np.float64)             # principal stress
swarm_sigma2=np.zeros(nmarker,dtype=np.float64)             # principal stress 
swarm_sw_level=np.zeros(nmarker,dtype=np.float64)           # level of strain weaking betw. 0 and 1

counter=0
for iel in range(0,nel):
    x1=xV[iconV[0,iel]] ; y1=yV[iconV[0,iel]]
    x2=xV[iconV[1,iel]] ; y2=yV[iconV[1,iel]]
    x3=xV[iconV[2,iel]] ; y3=yV[iconV[2,iel]]
    x4=xV[iconV[3,iel]] ; y4=yV[iconV[3,iel]]
    for j in range(0,nmarker_per_dim):
        for i in range(0,nmarker_per_dim):
            r=-1.+i*2./nmarker_per_dim + 1./nmarker_per_dim
            s=-1.+j*2./nmarker_per_dim + 1./nmarker_per_dim
            swarm_r[counter]=r
            swarm_s[counter]=s
            N1=0.25*(1-r)*(1-s)
            N2=0.25*(1+r)*(1-s)
            N3=0.25*(1+r)*(1+s)
            N4=0.25*(1-r)*(1+s)
            swarm_x[counter]=N1*x1+N2*x2+N3*x3+N4*x4
            swarm_y[counter]=N1*y1+N2*y2+N3*y3+N4*y4
            swarm_iel[counter]=iel
            counter+=1
        #end for 
    #end for 
#end for 

print("swarm setup: %.3f s" % (time.time() - start))

###############################################################################
# assign material id to markers, roughness, seeds, ... 
###############################################################################
start = time.time()

make_clast(nmarker,swarm_x,swarm_y,swarm_mat,swarm_plastic_strainxx,swarm_plastic_strainyy,swarm_plastic_strainxy)


print("create clast on swarm: %.3f s" % (time.time() - start))

###############################################################################
# paint markers 
###############################################################################
start = time.time()

for im in range (0,nmarker):
    for i in range(0,11):
        if abs(swarm_x[im]-i*Lx/10)<hx/4:
            swarm_paint[im]=1
    #end for 
    for i in range(0,5):
        if abs(swarm_y[im]-i*Ly/4)<hy/4:
            swarm_paint[im]=1
    #end for 
#end for 

print("paint swarm: %.3f s" % (time.time() - start))

###################################################################################################
###################################################################################################
# time stepping loop 
###################################################################################################
###################################################################################################

total_time=0
    
umem=np.zeros(NV,dtype=np.float64)      
vmem=np.zeros(NV,dtype=np.float64)    
u   =np.zeros(NV,dtype=np.float64)   
v   =np.zeros(NV,dtype=np.float64)    
exx = np.zeros(NV,dtype=np.float64)  
eyy = np.zeros(NV,dtype=np.float64)  
exy = np.zeros(NV,dtype=np.float64)  
        
# since all elements are identical, the jacobian and derived 
# products are all equal and can also be computed analytically
jcb   = np.array([[hx/2,0],[0,hy/2]],dtype=np.float64) 
jcob=hx*hy/4
jcbi  = np.array([[2/hx,0],[0,2/hy]],dtype=np.float64) 

for istep in range(0,nstep):

    print ('||******************************************************||')
    print ('||    istep= %i ' %istep)
    print ('||******************************************************||')

    ###########################################################################
    # define boundary conditions
    ###########################################################################
    start = time.time()

    bc_fix=np.zeros(NfemV,dtype=bool)        # boundary condition, yes/no
    bc_val=np.zeros(NfemV,dtype=np.float64)  # boundary condition, value

    define_bc(Lx,Ly,NV,bc_fix,bc_val,xV,yV,experiment,total_time)

    print("      boundary conditions: %.3f s" % (time.time() - start))

    ###########################################################################
    # localise markers 
    ###########################################################################
    start = time.time()

    nmarker_in_element=np.zeros(nel,dtype=np.int16)
    list_of_markers_in_element=np.zeros((2*nmarker_per_element,nel),dtype=np.int32)
    for im in range(0,nmarker):
        ielx=int(swarm_x[im]/Lx*nelx)
        iely=int(swarm_y[im]/Ly*nely)
        iel=nelx*(iely)+ielx
        list_of_markers_in_element[nmarker_in_element[iel],iel]=im
        nmarker_in_element[iel]+=1
        swarm_iel[im]=iel
    #end for

    print("     localise markers: %.3f s" % (time.time() - start))

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # non linear iterations
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    for iter in range(0,niter):

        print('     ----------------------------------')
        print('     ----- iteration------',iter,'----------')
        print('     ----------------------------------')

        #################################################################
        # compute elemental averagings 
        #################################################################
        start = time.time()
        eta_elemental=np.zeros(nel,dtype=np.float64)
        for im in range(0,nmarker):
            iel=swarm_iel[im]
            NNNV[0:mV]=NNV(swarm_r[im],swarm_s[im])
            swarm_exx[im]=sum(NNNV[0:mV]*exx[iconV[0:mV,iel]])
            swarm_eyy[im]=sum(NNNV[0:mV]*eyy[iconV[0:mV,iel]])
            swarm_exy[im]=sum(NNNV[0:mV]*exy[iconV[0:mV,iel]])
            swarm_ee[im]=np.sqrt(0.5*(swarm_exx[im]**2+swarm_eyy[im]**2+2*swarm_exy[im]**2) ) 
            swarm_eta[im],swarm_is_plastic[im],swarm_yield[im],dummy=\
            viscosity(swarm_x[im],swarm_y[im],swarm_ee[im],\
                      background_temperature,swarm_mat[im],iter,swarm_plastic_strain_eff[im])

            if abs(avrg)==1 : # arithmetic
               eta_elemental[iel]     +=swarm_eta[im]
            if abs(avrg)==2: # geometric
               eta_elemental[iel]     +=np.log10(swarm_eta[im])
            if abs(avrg)==3: # harmonic
               eta_elemental[iel]     +=1/swarm_eta[im]
        #end for
        if abs(avrg)==1:
           eta_elemental[:]/=nmarker_in_element[:]
        if abs(avrg)==2:
           eta_elemental[:]=10.**(eta_elemental[:]/nmarker_in_element[:])
        if abs(avrg)==3:
           eta_elemental[:]=nmarker_in_element[:]/eta_elemental[:]

        print("          -> nmarker_in_elt(m,M) %.5e %.5e " %(np.min(nmarker_in_element),np.max(nmarker_in_element)))
        print("          -> eta_elemental (m,M) %.5e %.5e " %(np.min(eta_elemental),np.max(eta_elemental)))

        print("     markers onto grid: %.3f s" % (time.time() - start))

        #################################################################
        # build FE matrix
        # [ K G ][u]=[f]
        # [GT 0 ][p] [h]
        #################################################################
        start = time.time()

        A_sparse = lil_matrix((Nfem,Nfem),dtype=np.float64)

        f_rhs   = np.zeros(NfemV,dtype=np.float64)         # right hand side f 
        h_rhs   = np.zeros(NfemP,dtype=np.float64)         # right hand side h 
        b_mat   = np.zeros((3,ndofV*mV),dtype=np.float64)  # gradient matrix B 
        N_mat   = np.zeros((3,ndofP*mP),dtype=np.float64)  # matrix N
        NNNV    = np.zeros(mV,dtype=np.float64)            # shape functions V
        NNNP    = np.zeros(mP,dtype=np.float64)            # shape functions P
        dNNNVdx = np.zeros(mV,dtype=np.float64)            # shape functions derivatives
        dNNNVdy = np.zeros(mV,dtype=np.float64)            # shape functions derivatives
        dNNNVdr = np.zeros(mV,dtype=np.float64)            # shape functions derivatives
        dNNNVds = np.zeros(mV,dtype=np.float64)            # shape functions derivatives
        RV      = np.zeros(NfemV,dtype=np.float64)         # nonlinear residual V
        RP      = np.zeros(NfemP,dtype=np.float64)         # nonlinear residual P
        p       = np.zeros(NfemP,dtype=np.float64)         # pressure solution vector

        c_mat   = np.array([[2,0,0],[0,2,0],[0,0,1]],dtype=np.float64) 
        

        for iel in range(0,nel):

            f_el =np.zeros((mV*ndofV),dtype=np.float64)
            K_el =np.zeros((mV*ndofV,mV*ndofV),dtype=np.float64)
            G_el=np.zeros((mV*ndofV,mP*ndofP),dtype=np.float64)
            h_el=np.zeros((mP*ndofP),dtype=np.float64)
            NNNNP= np.zeros(mP*ndofP,dtype=np.float64)   

            V_el=np.zeros((mV*ndofV),dtype=np.float64)
            P_el=np.zeros((mP*ndofP),dtype=np.float64)
            RV_el=np.zeros((mV*ndofV),dtype=np.float64)
            RP_el=np.zeros((mP*ndofP),dtype=np.float64)

            # element vector containing V and P solution
            for k in range(0,mV):
                V_el[2*k+0]=u[iconV[k,iel]]
                V_el[2*k+1]=v[iconV[k,iel]]
            for k in range(0,mP):
                P_el[k]=p[iconP[k,iel]]
            

            for iq in range(0,nqperdim):
                for jq in range(0,nqperdim):

                    # position & weight of quad. point
                    rq=qcoords[iq]
                    sq=qcoords[jq]
                    weightq=qweights[iq]*qweights[jq]

                    NNNV[0:mV]=NNV(rq,sq)
                    dNNNVdr[0:mV]=dNNVdr(rq,sq)
                    dNNNVds[0:mV]=dNNVds(rq,sq)
                    NNNP[0:mP]=NNP(rq,sq)

                    # calculate jacobian matrix
                    # replaced by analytical values above

                    # compute dNdx & dNdy
                    for k in range(0,mV):
                        dNNNVdx[k]=jcbi[0,0]*dNNNVdr[k]
                        dNNNVdy[k]=jcbi[1,1]*dNNNVds[k]

                    # construct b_mat matrix
                    for i in range(0,mV):
                        b_mat[0:3, 2*i:2*i+2] = [[dNNNVdx[i],0.       ],
                                                 [0.        ,dNNNVdy[i]],
                                                 [dNNNVdy[i],dNNNVdx[i]]]

                    # compute elemental a_mat matrix
                    K_el+=b_mat.T.dot(c_mat.dot(b_mat))*eta_elemental[iel]*weightq*jcob

                    # compute elemental rhs vector
                    #for i in range(0,mV):
                    #    f_el[ndofV*i  ]+=NNNV[i]*jcob*weightq*gx*density(xq,yq)
                    #    f_el[ndofV*i+1]+=NNNV[i]*jcob*weightq*gy*density(xq,yq)

                    for i in range(0,mP):
                        N_mat[0,i]=NNNP[i]
                        N_mat[1,i]=NNNP[i]
                        N_mat[2,i]=0.

                    G_el-=b_mat.T.dot(N_mat)*weightq*jcob

                    #NNNNP[:]+=NNNP[:]*jcob*weightq

                # end for jq
            # end for iq


            # impose b.c. 
            for k1 in range(0,mV):
                for i1 in range(0,ndofV):
                    ikk=ndofV*k1          +i1
                    m1 =ndofV*iconV[k1,iel]+i1
                    if bc_fix[m1]:
                       K_ref=K_el[ikk,ikk] 
                       for jkk in range(0,mV*ndofV):
                           f_el[jkk]-=K_el[jkk,ikk]*bc_val[m1]
                           K_el[ikk,jkk]=0
                           K_el[jkk,ikk]=0
                       K_el[ikk,ikk]=K_ref
                       f_el[ikk]=K_ref*bc_val[m1]
                       h_el[:]-=G_el[ikk,:]*bc_val[m1]
                       G_el[ikk,:]=0
                #end for
            #end for

            #----------------------------------
            # compute nonlinear residual
            RV_el=K_el.dot(V_el)+G_el.dot(P_el)-f_el
            RP_el=G_el.T.dot(V_el)-h_el

            # scaling of blocks 
            G_el*=eta_ref/Ly
            h_el*=eta_ref/Ly

            # assemble FEM matrix blocks and right hand side rhs
            for k1 in range(0,mV):
                for i1 in range(0,ndofV):
                    ikk=ndofV*k1          +i1
                    m1 =ndofV*iconV[k1,iel]+i1
                    for k2 in range(0,mV):
                        for i2 in range(0,ndofV):
                            jkk=ndofV*k2          +i2
                            m2 =ndofV*iconV[k2,iel]+i2
                            A_sparse[m1,m2] += K_el[ikk,jkk]
                        #end for
                    #end for
                    for k2 in range(0,mP):
                        jkk=k2
                        m2 =iconP[k2,iel]
                        A_sparse[m1,NfemV+m2]+=G_el[ikk,jkk]
                        A_sparse[NfemV+m2,m1]+=G_el[ikk,jkk]
                    #end for
                    f_rhs[m1]+=f_el[ikk]
                    RV[m1]+=RV_el[ikk]
                #end for
            #end for
            for k2 in range(0,mP):
                m2=iconP[k2,iel]
                h_rhs[m2]+=h_el[k2]
                RP[m2]+=RP_el[k2]
            #end for
        #end for iel

        print("     build FE matrix: %.3f s" % (time.time() - start))

        ######################################################################
        # assemble K, G, GT, f, h into A and rhs
        ######################################################################
        start = time.time()

        rhs=np.zeros(Nfem,dtype=np.float64)   # right hand side of Ax=b
        rhs[0:NfemV]=f_rhs
        rhs[NfemV:Nfem]=h_rhs

        sparse_matrix=A_sparse.tocsr()

        print("     assemble rhs, convert to csr: %.3f s" % (time.time() - start))

        ######################################################################
        # solve system
        ######################################################################
        start = time.time()

        sol=sps.linalg.spsolve(sparse_matrix,rhs)

        print("     solve time: %.3f s" % (time.time() - start))

        ######################################################################
        # put solution into separate x,y velocity arrays
        ######################################################################
        start = time.time()

        u,v=np.reshape(sol[0:NfemV],(NV,2)).T
        p=sol[NfemV:Nfem]*(eta_ref/Ly)

        print("          -> u (m,M) %.4e %.4e (cm/year)" %(np.min(u)/cm*year,np.max(u)/cm*year))
        print("          -> v (m,M) %.4e %.4e (cm/year)" %(np.min(v)/cm*year,np.max(v)/cm*year))
        print("          -> p (m,M) %.4e %.4e (Pa)     " %(np.min(p),np.max(p)))

        print("     split vel into u,v: %.3f s" % (time.time() - start))

        ######################################################################
        #normalise pressure
        ######################################################################
        start = time.time()

        avrg_p=0
        for iel in range(0,nel):
            for iq in range(0,nqperdim):
                for jq in range(0,nqperdim):
                    rq=qcoords[iq]
                    sq=qcoords[jq]
                    weightq=qweights[iq]*qweights[jq]
                    NNNP[0:mP]=NNP(rq,sq)
                    pq=NNNP.dot(p[iconP[0:mP,iel]])
                    jcob=hx*hy/4
                    avrg_p+=pq*jcob*weightq

        print('          -> avrg_p=',avrg_p)

        p-=avrg_p/Lx/Ly

        print("          -> p (m,M) %.4e %.4e (Pa)     " %(np.min(p),np.max(p)))
            
        print("     pressure normalisation: %.3f s" % (time.time() - start))

        ######################################################################
        # compute nodal strain rate 
        ######################################################################
        start = time.time()
 
        exx = np.zeros(NV,dtype=np.float64)  
        eyy = np.zeros(NV,dtype=np.float64)  
        exy = np.zeros(NV,dtype=np.float64)  
        ee  = np.zeros(NV,dtype=np.float64)  
        ccc = np.zeros(NV,dtype=np.float64)  
 
        for iel in range(0,nel):
            for k in range(0,mV):
                rq = rVnodes[k]
                sq = sVnodes[k]
                inode=iconV[k,iel]
                NNNV[0:mV]=NNV(rq,sq)
                dNNNVdr[0:mV]=dNNVdr(rq,sq)
                dNNNVds[0:mV]=dNNVds(rq,sq)
                #jcb=np.zeros((2,2),dtype=np.float64)
                #for k in range(0,mV):
                #    jcb[0,0]+=dNNNVdr[k]*xV[iconV[k,iel]]
                #    jcb[0,1]+=dNNNVdr[k]*yV[iconV[k,iel]]
                #    jcb[1,0]+=dNNNVds[k]*xV[iconV[k,iel]]
                #    jcb[1,1]+=dNNNVds[k]*yV[iconV[k,iel]]
                #jcbi=np.linalg.inv(jcb)
                for k in range(0,mV):
                    dNNNVdx[k]=jcbi[0,0]*dNNNVdr[k]+jcbi[0,1]*dNNNVds[k]
                    dNNNVdy[k]=jcbi[1,0]*dNNNVdr[k]+jcbi[1,1]*dNNNVds[k]
                ccc[inode]+=1
                exx[inode]+=dNNNVdx.dot(u[iconV[:,iel]])
                eyy[inode]+=dNNNVdy.dot(v[iconV[:,iel]])
                exy[inode]+=dNNNVdx.dot(v[iconV[:,iel]])*0.5+dNNNVdy.dot(u[iconV[:,iel]])*0.5
            #end for
        #end for
        exx[:]/=ccc[:]
        eyy[:]/=ccc[:]
        exy[:]/=ccc[:]
        ee[:]=np.sqrt(0.5*(exx[:]*exx[:]+eyy[:]*eyy[:])+exy[:]*exy[:])
 
        print("          -> exx (m,M) %.4e %.4e " %(np.min(exx),np.max(exx)))
        print("          -> eyy (m,M) %.4e %.4e " %(np.min(eyy),np.max(eyy)))
        print("          -> exy (m,M) %.4e %.4e " %(np.min(exy),np.max(exy)))
 
        print("     compute nodal strain rate: %.3f s" % (time.time() - start))

        ######################################################################
        # convergence criterion is based on difference between two consecutively
        # obtained velocity fields, normalised by the boundary condition velocity

        chi_u=LA.norm(u-umem,2)/v_ref # vx convergence indicator
        chi_v=LA.norm(v-vmem,2)/v_ref # vy convergence indicator

        if iter==0: RV0=1 #LA.norm(RV,2)

        print('          -> Nonlinear residual V,P:',LA.norm(RV,2)/RV0,LA.norm(RP,2))

        print('          -> convergence u,v: %.3e %.3e | tol= %.2e' %(chi_u,chi_v,tol))

        convfile.write("%f %10e %10e %10e %e %e \n" %(istep+iter/200,chi_u,chi_v,LA.norm(RV,2)/RV0,LA.norm(RP,2),tol))
        convfile.flush()

        umem[:]=u[:]
        vmem[:]=v[:]

        if chi_u<tol and chi_v<tol:
           print('     ***************')
           print('     ***converged***')
           print('     ***************')
           break

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #end for iter
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    print('     ----------------------------------')
    print('     end of nl iteration ')
    print('     ----------------------------------')

    ######################################################################
    # compute timestep using CFL condition 
    ######################################################################

    dt=CFL_nb*min(hx,hy)/max(max(abs(u)),max(abs(v)))

    total_time+=dt

    print("     -> dt= %.3e yr" %(dt/year))
    print("     -> time= %.3e yr" %(total_time/year))

    ######################################################################
    # advect markers 
    # a simple one-step Euler method is used so timesteps should
    # be kept to rather low values, i.e. CFL_nb<=0.25
    # Periodic boundary conditions are implemented on the markers.
    ######################################################################
    start = time.time()

    for im in range(0,nmarker):
           ielx=int(swarm_x[im]/Lx*nelx)
           iely=int(swarm_y[im]/Ly*nely)
           iel=nelx*iely+ielx
           swarm_r[im]=-1.+2*(swarm_x[im]-xV[iconV[0,iel]])/hx
           swarm_s[im]=-1.+2*(swarm_y[im]-yV[iconV[0,iel]])/hy
           NNNV[0:mV]=NNV(swarm_r[im],swarm_s[im])
           um=sum(NNNV[0:mV]*u[iconV[0:mV,iel]])
           vm=sum(NNNV[0:mV]*v[iconV[0:mV,iel]])
           exxm=sum(NNNV[0:mV]*exx[iconV[0:mV,iel]])
           eyym=sum(NNNV[0:mV]*eyy[iconV[0:mV,iel]])
           exym=sum(NNNV[0:mV]*exy[iconV[0:mV,iel]])
           #assign velocity to marker
           swarm_u[im]=um 
           swarm_v[im]=vm
           #update its position
           swarm_x[im]+=um*dt
           swarm_y[im]+=vm*dt
           #update its strain tensor components
           swarm_total_strainxx[im]+=exxm*dt
           swarm_total_strainyy[im]+=eyym*dt
           swarm_total_strainxy[im]+=exym*dt
           if swarm_is_plastic[im]==1: 
              swarm_plastic_strainxx[im]+=exxm*dt
              swarm_plastic_strainyy[im]+=eyym*dt
              swarm_plastic_strainxy[im]+=exym*dt

           swarm_total_strain_eff[im]  =effective(swarm_total_strainxx[im],swarm_total_strainyy[im],swarm_total_strainxy[im])
           swarm_plastic_strain_eff[im]=effective(swarm_plastic_strainxx[im],swarm_plastic_strainyy[im],swarm_plastic_strainxy[im])

           #assign strain rate tensor components
           swarm_exx[im]=exxm
           swarm_eyy[im]=eyym
           swarm_exy[im]=exym
           #assign effective strain rate and viscosity
           swarm_ee[im]=effective(swarm_exx[im],swarm_eyy[im],swarm_exy[im]) 
           swarm_eta[im],swarm_is_plastic[im],swarm_yield[im],swarm_sw_level[im]=\
           viscosity(swarm_x[im],swarm_y[im],swarm_ee[im],
                     background_temperature,swarm_mat[im],iter,swarm_plastic_strain_eff[im])
           #assign dev stress values
           swarm_tauxx[im]=2*exxm*swarm_eta[im]
           swarm_tauyy[im]=2*eyym*swarm_eta[im]
           swarm_tauxy[im]=2*exym*swarm_eta[im]
           #assign pressure
           NNNP[0:mP]=NNP(swarm_r[im],swarm_s[im])
           swarm_p_dyn[im]=NNNP.dot(p[iconP[0:mP,iel]])

           if swarm_x[im]>Lx: swarm_x[im]-=Lx #periodic b.c. on right side
           if swarm_x[im]<0:  swarm_x[im]+=Lx #periodic b.c. on left side
       #end for
    #end for

    swarm_sigmaxx[:]=-swarm_p_dyn[:]+swarm_tauxx[:]
    swarm_sigmayy[:]=-swarm_p_dyn[:]+swarm_tauyy[:]
    swarm_sigmaxy[:]=               +swarm_tauxy[:]
    swarm_sigma_angle[:]=0.5*np.arctan(2*swarm_sigmaxy[:]/(swarm_sigmayy[:]-swarm_sigmaxx[:])) 

    swarm_tau_eff[:]=effective(swarm_tauxx[:],swarm_tauyy[:],swarm_tauxy[:])

    swarm_tau_angle[:]=0.5*np.arctan(2*swarm_tauxy[:]/(swarm_tauyy[:]-swarm_tauxx[:])) 

    swarm_sigma1[:]=(swarm_sigmaxx[:]+swarm_sigmayy[:])/2. \
                   + np.sqrt( (swarm_sigmaxx[:]-swarm_sigmayy[:])**2/4 +swarm_sigmaxy[:]**2 ) 
    swarm_sigma2[:]=(swarm_sigmaxx[:]+swarm_sigmayy[:])/2. \
                   - np.sqrt( (swarm_sigmaxx[:]-swarm_sigmayy[:])**2/4 +swarm_sigmaxy[:]**2 ) 

    print("     advect markers: %.3f s" % (time.time() - start))

    if istep==0: swarm_plastic_strain_eff0[:]=swarm_plastic_strain_eff[:]

    #####################################################################
    # interpolate pressure onto velocity grid points for paraview output
    #####################################################################
    start = time.time()

    q=np.zeros(NV,dtype=np.float64)
    counter=np.zeros(NV,dtype=np.float64)

    for iel in range(0,nel):
        q[iconV[0,iel]]=p[iconP[0,iel]]
        q[iconV[1,iel]]=p[iconP[1,iel]]
        q[iconV[2,iel]]=p[iconP[2,iel]]
        q[iconV[3,iel]]=p[iconP[3,iel]]
        q[iconV[4,iel]]=(p[iconP[0,iel]]+p[iconP[1,iel]])*0.5
        q[iconV[5,iel]]=(p[iconP[1,iel]]+p[iconP[2,iel]])*0.5
        q[iconV[6,iel]]=(p[iconP[2,iel]]+p[iconP[3,iel]])*0.5
        q[iconV[7,iel]]=(p[iconP[3,iel]]+p[iconP[0,iel]])*0.5
        q[iconV[8,iel]]=(p[iconP[0,iel]]+p[iconP[1,iel]]+p[iconP[2,iel]]+p[iconP[3,iel]])*0.25

    print("     project p on Q2: %.3f s" % (time.time() - start))

    ###########################################################################
    # export solution to vtu
    ###########################################################################
    start = time.time()

    if istep%every_vtu==0:
       export_solution_to_vtu(NV,nel,xV,yV,iconV,u,v,q,eta_elemental,exx,eyy,exy,ee,output_folder,istep)

    print("     export solution to vtu: %.3f s" % (time.time() - start))

    ###########################################################################
    # export solution to vtu
    ###########################################################################
    start = time.time()

    if istep%every_vtu==0:
       export_swarm_to_vtu(nmarker,swarm_mat,\
                                swarm_paint,\
                                swarm_total_strainxx,\
                                swarm_total_strainyy,\
                                swarm_total_strainxy,\
                                swarm_total_strain_eff,\
                                swarm_plastic_strainxx,\
                                swarm_plastic_strainyy,\
                                swarm_plastic_strainxy,\
                                swarm_plastic_strain_eff,\
                                swarm_plastic_strain_eff0,\
                                swarm_exx,\
                                swarm_eyy,\
                                swarm_exy,\
                                swarm_ee,\
                                swarm_sw_level,\
                                swarm_is_plastic,\
                                swarm_tauxx,\
                                swarm_tauyy,\
                                swarm_tauxy,\
                                swarm_tau_eff,\
                                swarm_tau_angle,\
                                swarm_sigma_angle,\
                                swarm_sigma1,\
                                swarm_sigma2,\
                                swarm_eta,\
                                swarm_p_dyn,\
                                swarm_u,\
                                swarm_v,\
                                swarm_x,\
                                swarm_y,\
                                output_folder,istep)

    print("     export swarm to vtu: %.3f s" % (time.time() - start))

    ###########################################################################
    # export solution to ascii
    ###########################################################################

    #np.savetxt('velocity.ascii',np.array([xV,yV,u,v]).T,header='# x,y,u,v')
    #np.savetxt('pressure_aft.ascii',np.array([xP,yP,p]).T,header='# x,y,p')
    np.savetxt('markers.ascii',np.array([swarm_x,swarm_y,swarm_mat,swarm_eta,swarm_plastic_strain_eff]).T)

    ###########################################################################
    # export solution to pdf via matplotlib
    ###########################################################################
    start = time.time()

    if istep%every_pdf==0:
       export_solution_to_pdf(Lx,Ly,nnx,nny,nelx,nely,xV,yV,u,v,q,exx,eyy,exy,ee,\
                              eta_elemental,q,output_folder,istep)

    print("     export solution to pdf: %.3f s" % (time.time() - start))

    ###########################################################################
    # export swarm to pdf via matplotlib
    ###########################################################################
    start = time.time()

    if istep%every_pdf==11110:
       export_swarm_to_pdf(Lx,Ly,swarm_x,\
                                 swarm_y,\
                                 swarm_eta,\
                                 swarm_mat,\
                                 swarm_ee,\
                                 swarm_total_strain_eff,\
                                 swarm_tau_eff,\
                                 swarm_plastic_strain_eff,\
                                 output_folder,istep)

    print("     export swarm to pdf: %.3f s" % (time.time() - start))

    ###########################################################################
        
    if total_time>tfinal:
       print('     ***tfinal reached**')
       break

###################################################################################################
# end for istep
###################################################################################################

print("-----------------------------")
print("------------the end----------")
print("-----------------------------")