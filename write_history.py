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

def write_history(output_folder,total_time,istep,u,v,ee):

    if istep==0:
       histfile=open(output_folder+'hist.ascii',"w")
       histfile.write("# col 1: istep      \n")
       histfile.write("# col 2: total_time \n")
       histfile.write("# col 3: min(u)     \n")
       histfile.write("# col 4: max(u)     \n")
       histfile.write("# col 5: min(v)     \n")
       histfile.write("# col 6: max(v)     \n")
       histfile.write("# col 7: min(ee)    \n")
       histfile.write("# col 8: max(ee)    \n")
    else:
       histfile=open(output_folder+'hist.ascii',"a")

    histfile.write("%d %e %e %e %e %e %e %e\n" %(\
                    istep,     \
                    total_time,\
                    np.min(u), \
                    np.max(u), \
                    np.min(v), \
                    np.max(v), \
                    np.min(ee),\
                    np.max(ee)))

    histfile.close()

###############################################################################
