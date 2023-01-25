import numpy as np
from constants_and_tools import *
from basis_functions import * 
from analytical_solutions import *

def compute_errors(Lx,Ly,nel,hx,xV,yV,u,v,p,iconV,iconP,experiment):

   errv=0
   errp=0
   vrms=0

   for iel in range(0,nel):
       for iq in range(0,nqperdim):
           for jq in range(0,nqperdim):
               rq=qcoords[iq]
               sq=qcoords[jq]
               weightq=qweights[iq]*qweights[jq]
               NNNV=NNV(rq,sq)
               dNNNVdr=dNNVdr(rq,sq)
               dNNNVds=dNNVds(rq,sq)
               NNNP=NNP(rq,sq)
               jcb=np.zeros((2,2),dtype=np.float64)
               for k in range(0,mV):
                   jcb[0,0] += dNNNVdr[k]*xV[iconV[k,iel]]
                   jcb[0,1] += dNNNVdr[k]*yV[iconV[k,iel]]
                   jcb[1,0] += dNNNVds[k]*xV[iconV[k,iel]]
                   jcb[1,1] += dNNNVds[k]*yV[iconV[k,iel]]
               xq=0
               yq=0
               for k in range(0,mV):
                   xq+=NNNV[k]*xV[iconV[k,iel]]
                   yq+=NNNV[k]*yV[iconV[k,iel]]
               jcob = np.linalg.det(jcb)

               uq=NNNV.dot(u[iconV[0:mV,iel]])
               vq=NNNV.dot(v[iconV[0:mV,iel]])
               pq=NNNP.dot(p[iconP[0:mP,iel]])
               vrms+=(uq**2+vq**2)*weightq*jcob

               a,b,c=analytical_solution(xq,yq,experiment)
               errv+=(uq-a)**2*weightq*jcob+(vq-b)**2*weightq*jcob
               errp+=(pq-c)**2*weightq*jcob

           #end for
       #end for
   #end for

   vrms=np.sqrt(vrms/(Lx*Ly))
   errv=np.sqrt(errv/(Lx*Ly))
   errp=np.sqrt(errp/(Lx*Ly))

   print("     errv,errp,vrms,hx-> ",errv,errp,vrms,hx)


