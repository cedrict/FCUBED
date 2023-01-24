import numpy as np

#------------------------------------------------------------------------------
# basic physical constants
#------------------------------------------------------------------------------

Rgas=8.3145
Tkelvin=273
cm=0.01
km=1e3
year=3600*24*365.25
MPa=1e6
Myr=1e6*year
eps=1.e-10

mV=9     # number of velocity nodes making up an element
mP=4     # number of pressure nodes making up an element
ndofV=2  # number of velocity degrees of freedom per node
ndofP=1  # number of pressure degrees of freedom 

###############################################################################

def effective(xx,yy,xy):
    return np.sqrt(0.5*(xx**2+yy**2)+xy**2)

###############################################################################

