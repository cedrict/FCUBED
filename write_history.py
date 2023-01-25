import numpy as np

def write_history(output_folder,total_time,istep,u,v):

    if istep==0:
       histfile=open(output_folder+'hist.ascii',"w")
       histfile.write("# col 1: istep \n")
       histfile.write("# col 2: total_time \n")
       histfile.write("# col 3: min(u) \n")
       histfile.write("# col 4: max(u) \n")
       histfile.write("# col 5: min(v) \n")
       histfile.write("# col 6: max(v) \n")
    else:
       histfile=open(output_folder+'hist.ascii',"a")
    histfile.write("%d %e %e %e %e %e\n" %(\
                    istep,     \
                    total_time,\
                    np.min(u), \
                    np.max(u), \
                    np.min(v), \
                    np.max(v)))
    histfile.close()
