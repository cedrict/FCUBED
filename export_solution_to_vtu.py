###############################################################################
#
#  FFFF  CCCC  U   U  BBB   EEEE  DDD      C.Thieulot
#  F     C     U   U  B  B  E     D  D     F.Gueydan
#  FFF   C     U   U  BBB   EEEE  D  D     A.Lemaitre
#  F     C     U   U  B  B  E     D  D
#  F     CCCC  UUUUU  BBB   EEEE  DDD
#
###############################################################################
# the 9-node Q2 elt does not exist in vtk, but the 8-node one does, ie type=23. 

from constants_and_tools import *
from inputs import *
from analytical_solutions import *

def export_solution_to_vtu(NV,nel,xV,yV,iconV,u,v,q,eta_elemental,rho_elemental,exx,eyy,exy,ee,\
                           Pf,phi,K,plastic_strain_eff,H,u_darcy,v_darcy,output_folder,istep):

    filename=output_folder+'solution_{:04d}.vtu'.format(istep)
    vtufile=open(filename,"w")
    vtufile.write("<VTKFile type='UnstructuredGrid' version='0.1' byte_order='BigEndian'> \n")
    vtufile.write("<UnstructuredGrid> \n")
    vtufile.write("<Piece NumberOfPoints=' %5d ' NumberOfCells=' %5d '> \n" %(NV,nel))
    #####
    vtufile.write("<Points> \n")
    vtufile.write("<DataArray type='Float32' NumberOfComponents='3' Format='ascii'> \n")
    for i in range(0,NV):
        vtufile.write("%10e %10e %10e \n" %(xV[i],yV[i],0.))
    vtufile.write("</DataArray>\n")
    vtufile.write("</Points> \n")
    #####
    vtufile.write("<CellData Scalars='scalars'>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='viscosity' Format='ascii'> \n")
    for iel in range (0,nel):
        vtufile.write("%10e\n" % eta_elemental[iel]) 
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='density' Format='ascii'> \n")
    for iel in range (0,nel):
        vtufile.write("%10e\n" % rho_elemental[iel]) 
    vtufile.write("</DataArray>\n")

    #--
    if use_fluid:
       vtufile.write("<DataArray type='Float32' NumberOfComponents='3' Name='velocity Darcy' Format='ascii'> \n")
       for iel in range (0,nel):
           vtufile.write("%e %e %e\n" % (u_darcy[iel],v_darcy[iel],0)) 
       vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='plastic strain (eff)' Format='ascii'> \n")
    for iel in range (0,nel):
        vtufile.write("%10e\n" % plastic_strain_eff[iel]) 
    vtufile.write("</DataArray>\n")
    #--
    if use_fluid:
       vtufile.write("<DataArray type='Float32' Name='H' Format='ascii'> \n")
       for iel in range (0,nel):
           vtufile.write("%10e\n" % H[iel]) 
       vtufile.write("</DataArray>\n")
       #--
       vtufile.write("<DataArray type='Float32' Name='porosity (phi)' Format='ascii'> \n")
       for iel in range (0,nel):
           vtufile.write("%10e\n" % phi[iel]) 
       vtufile.write("</DataArray>\n")
       #--
       vtufile.write("<DataArray type='Float32' Name='permeability (K)' Format='ascii'> \n")
       for iel in range (0,nel):
           vtufile.write("%10e\n" % K[iel]) 
       vtufile.write("</DataArray>\n")
    #--
    vtufile.write("</CellData>\n")
    #####
    vtufile.write("<PointData Scalars='scalars'>\n")
    #--
    vtufile.write("<DataArray type='Float32' NumberOfComponents='3' Name='velocity (m/s)' Format='ascii'> \n")
    for i in range(0,NV):
        vtufile.write("%10e %10e %10e \n" %(u[i],v[i],0.))
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' NumberOfComponents='3' Name='velocity (cm/year)' Format='ascii'> \n")
    for i in range(0,NV):
        vtufile.write("%10e %10e %10e \n" %(u[i]/cm*year,v[i]/cm*year,0.))
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='pressure' Format='ascii'> \n")
    for i in range(0,NV):
        vtufile.write("%10e \n" %q[i])
    vtufile.write("</DataArray>\n")
    #--
    if use_fluid:
       vtufile.write("<DataArray type='Float32' Name='pore fluid pressure (Pf)' Format='ascii'> \n")
       for i in range(0,NV):
           vtufile.write("%15e \n" %Pf[i])
       vtufile.write("</DataArray>\n")
       #--
       vtufile.write("<DataArray type='Float32' Name='pore fluid pressure (Pf/p_ref)' Format='ascii'> \n")
       for i in range(0,NV):
           vtufile.write("%15e \n" %(Pf[i]/p_ref))
       vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='exx' Format='ascii'> \n")
    for i in range (0,NV):
        vtufile.write("%10e\n" % exx[i])
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='eyy' Format='ascii'> \n")
    for i in range (0,NV):
        vtufile.write("%10e\n" % eyy[i])
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='exy' Format='ascii'> \n")
    for i in range (0,NV):
        vtufile.write("%10e\n" % exy[i])
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='ee' Format='ascii'> \n")
    for i in range (0,NV):
        vtufile.write("%10e\n" % ee[i])
    vtufile.write("</DataArray>\n")

    if experiment<0: 
       #--
       vtufile.write("<DataArray type='Float32' NumberOfComponents='3' Name='velocity (analytical)' Format='ascii'> \n")
       for i in range (0,NV):
           ui,vi,pi=analytical_solution(xV[i],yV[i],experiment)
           vtufile.write("%e %e %e\n" % (ui,vi,0)) 
       vtufile.write("</DataArray>\n")
       #--
       vtufile.write("<DataArray type='Float32' Name='pressure (analytical)' Format='ascii'> \n")
       for i in range (0,NV):
           ui,vi,pi=analytical_solution(xV[i],yV[i],experiment)
           vtufile.write("%e\n" % pi) 
       vtufile.write("</DataArray>\n")
    #--
    vtufile.write("</PointData>\n")
    #####
    vtufile.write("<Cells>\n")
    #--
    vtufile.write("<DataArray type='Int32' Name='connectivity' Format='ascii'> \n")
    for iel in range (0,nel):
        vtufile.write("%d %d %d %d %d %d %d %d\n" %(iconV[0,iel],iconV[1,iel],iconV[2,iel],\
                                                    iconV[3,iel],iconV[4,iel],iconV[5,iel],\
                                                    iconV[6,iel],iconV[7,iel]))
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Int32' Name='offsets' Format='ascii'> \n")
    for iel in range (0,nel):
        vtufile.write("%d \n" %((iel+1)*8))
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Int32' Name='types' Format='ascii'>\n")
    for iel in range (0,nel):
        vtufile.write("%d \n" %23)
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("</Cells>\n")
    #####
    vtufile.write("</Piece>\n")
    vtufile.write("</UnstructuredGrid>\n")
    vtufile.write("</VTKFile>\n")
    vtufile.close()

###############################################################################
