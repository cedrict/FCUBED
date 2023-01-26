###############################################################################
#
#  FFFF  CCCC  U   U  BBB   EEEE  DDD      C.Thieulot
#  F     C     U   U  B  B  E     D  D     F.Gueydan
#  FFF   C     U   U  BBB   EEEE  D  D     A.Lemaitre
#  F     C     U   U  B  B  E     D  D
#  F     CCCC  UUUUU  BBB   EEEE  DDD
#
###############################################################################

from constants_and_tools import *

def export_swarm_to_vtu(nmarker,swarm_mat,\
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
                                output_folder,istep):

    filename = output_folder+'swarm_{:04d}.vtu'.format(istep)
    vtufile=open(filename,"w")
    vtufile.write("<VTKFile type='UnstructuredGrid' version='0.1' byte_order='BigEndian'> \n")
    vtufile.write("<UnstructuredGrid> \n")
    vtufile.write("<Piece NumberOfPoints=' %5d ' NumberOfCells=' %5d '> \n" %(nmarker,nmarker))
    vtufile.write("<PointData Scalars='scalars'>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='mat' Format='ascii'>\n")
    for i in range(0,nmarker):
        vtufile.write("%3e \n" %swarm_mat[i])
    vtufile.write("</DataArray>\n")
    #--
    #vtufile.write("<DataArray type='Float32' Name='iel' Format='ascii'>\n")
    #for i in range(0,nmarker):
    #    vtufile.write("%3e \n" %swarm_iel[i])
    #vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='paint' Format='ascii'>\n")
    for i in range(0,nmarker):
        vtufile.write("%3e \n" %swarm_paint[i])
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='total strain (xx)' Format='ascii'>\n")
    for i in range(0,nmarker):
        vtufile.write("%3e \n" %swarm_total_strainxx[i])
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='total strain (yy)' Format='ascii'>\n")
    for i in range(0,nmarker):
        vtufile.write("%3e \n" %swarm_total_strainyy[i])
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='total strain (xy)' Format='ascii'>\n")
    for i in range(0,nmarker):
        vtufile.write("%3e \n" %swarm_total_strainxy[i])
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='total strain (eff.)' Format='ascii'>\n")
    for i in range(0,nmarker):
        vtufile.write("%3e \n" % swarm_total_strain_eff[i] )
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='plastic strain (xx)' Format='ascii'>\n")
    for i in range(0,nmarker):
        vtufile.write("%3e \n" %swarm_plastic_strainxx[i])
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='plastic strain (yy)' Format='ascii'>\n")
    for i in range(0,nmarker):
        vtufile.write("%3e \n" %swarm_plastic_strainyy[i])
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='plastic strain (xy)' Format='ascii'>\n")
    for i in range(0,nmarker):
        vtufile.write("%3e \n" %swarm_plastic_strainxy[i])
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='plastic strain (eff.)' Format='ascii'>\n")
    for i in range(0,nmarker):
        vtufile.write("%3e \n" % swarm_plastic_strain_eff[i] )
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='plastic strain (eff.) noninitial' Format='ascii'>\n")
    for i in range(0,nmarker):
        vtufile.write("%3e \n" % (swarm_plastic_strain_eff[i]-swarm_plastic_strain_eff0[i]))
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='strain rate (xx)' Format='ascii'>\n")
    for i in range(0,nmarker):
        vtufile.write("%3e \n" %swarm_exx[i])
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='strain rate (yy)' Format='ascii'>\n")
    for i in range(0,nmarker):
        vtufile.write("%3e \n" %swarm_eyy[i])
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='strain rate (xy)' Format='ascii'>\n")
    for i in range(0,nmarker):
        vtufile.write("%3e \n" %swarm_exy[i])
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='strain rate (eff.)' Format='ascii'>\n")
    for i in range(0,nmarker):
        vtufile.write("%3e \n" % swarm_ee[i] )
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='strain weakening level' Format='ascii'>\n")
    for i in range(0,nmarker):
        vtufile.write("%3e \n" % swarm_sw_level[i] )
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='is_plastic' Format='ascii'>\n")
    for i in range(0,nmarker):
        vtufile.write("%3e \n" % swarm_is_plastic[i] )
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='tauxx (MPa)' Format='ascii'>\n")
    for i in range(0,nmarker):
        vtufile.write("%3e \n" %(swarm_tauxx[i]/MPa))
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='tauyy (MPa)' Format='ascii'>\n")
    for i in range(0,nmarker):
        vtufile.write("%3e \n" %(swarm_tauyy[i]/MPa))
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='tauxy (MPa)' Format='ascii'>\n")
    for i in range(0,nmarker):
        vtufile.write("%3e \n" %(swarm_tauxy[i]/MPa))
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='tau (eff.) (MPa)' Format='ascii'>\n")
    for i in range(0,nmarker):
        vtufile.write("%3e \n" % (swarm_tau_eff[i]/MPa))
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='tau angle' Format='ascii'>\n")
    for i in range(0,nmarker):
        vtufile.write("%3e \n" % (swarm_tau_angle[i]/np.pi*180))
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='sigma angle (abs)' Format='ascii'>\n")
    for i in range(0,nmarker):
        vtufile.write("%3e \n" % (np.abs(swarm_sigma_angle[i])/np.pi*180))
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='sigma 1 (MPa)' Format='ascii'>\n")
    for i in range(0,nmarker):
        vtufile.write("%e \n" % (swarm_sigma1[i]/MPa))
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='sigma 2 (MPa)' Format='ascii'>\n")
    for i in range(0,nmarker):
        vtufile.write("%e \n" % (swarm_sigma2[i]/MPa))
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='sigma_1 (dir)' NumberOfComponents='3' Format='ascii'> \n")
    for i in range(0,nmarker):
        vtufile.write("%10e %10e %10e \n" %( np.cos(swarm_sigma_angle[i]),np.sin(swarm_sigma_angle[i]),0) )
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='sigma_2 (dir)' NumberOfComponents='3' Format='ascii'> \n")
    for i in range(0,nmarker):
        vtufile.write("%10e %10e %10e \n" %( np.cos(swarm_sigma_angle[i]+np.pi/2),np.sin(swarm_sigma_angle[i]+np.pi/2),0) )
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='viscosity (Pa.s)' Format='ascii'>\n")
    for i in range(0,nmarker):
        vtufile.write("%10e \n" %swarm_eta[i])
    vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='pressure dyn. (MPa)' Format='ascii'>\n")
    for i in range(0,nmarker):
        vtufile.write("%10e \n" %(swarm_p_dyn[i]/MPa))
    vtufile.write("</DataArray>\n")
    #--
    #vtufile.write("<DataArray type='Float32' Name='yield value(MPa)' Format='ascii'>\n")
    #for i in range(0,nmarker):
    #    vtufile.write("%10e \n" %(swarm_tau_eff[i]/MPa-swarm_yield[i]/MPa))
    #vtufile.write("</DataArray>\n")
    #--
    #vtufile.write("<DataArray type='Float32' Name='r,s,t' NumberOfComponents='3' Format='ascii'>\n")
    #for i in range(0,nmarker):
    #    vtufile.write("%5e %5e %5e \n" %(swarm_r[i],swarm_s[i],0.))
    #vtufile.write("</DataArray>\n")
    #--
    vtufile.write("<DataArray type='Float32' Name='velocity (cm/year)' NumberOfComponents='3' Format='ascii'>\n")
    for i in range(0,nmarker):
        vtufile.write("%5e %5e %5e \n" %(swarm_u[i]/cm*year,swarm_v[i]/cm*year,0.))
    vtufile.write("</DataArray>\n")
    #--
    #vtufile.write("<DataArray type='Float32' NumberOfComponents='3' Name='velocity-v_simpleshear (cm/year)' Format='ascii'> \n")
    #for i in range(0,nmarker):
    #    vtufile.write("%10e %10e %10e \n" %((swarm_u[i]-(2*v0/Ly*swarm_y[i]-v0))/cm*year ,swarm_v[i]/cm*year,0.))
    #vtufile.write("</DataArray>\n")

    #--
    vtufile.write("</PointData>\n")
    vtufile.write("<Points> \n")
    vtufile.write("<DataArray type='Float32' NumberOfComponents='3' Format='ascii'>\n")
    for i in range(0,nmarker):
        vtufile.write("%10e %10e %10e \n" %(swarm_x[i],swarm_y[i],0.))
    vtufile.write("</DataArray>\n")
    vtufile.write("</Points> \n")
    vtufile.write("<Cells>\n")
    vtufile.write("<DataArray type='Int32' Name='connectivity' Format='ascii'> \n")
    for i in range(0,nmarker):
        vtufile.write("%d " % i)
    vtufile.write("</DataArray>\n")
    vtufile.write("<DataArray type='Int32' Name='offsets' Format='ascii'> \n")
    for i in range(0,nmarker):
        vtufile.write("%d " % (i+1))
    vtufile.write("</DataArray>\n")
    vtufile.write("<DataArray type='Int32' Name='types' Format='ascii'>\n")
    for i in range(0,nmarker):
        vtufile.write("%d " % 1)
    vtufile.write("</DataArray>\n")
    vtufile.write("</Cells>\n")
    vtufile.write("</Piece>\n")
    vtufile.write("</UnstructuredGrid>\n")
    vtufile.write("</VTKFile>\n")
    vtufile.close()

###############################################################################
