#!/opt/local/bin/python

####################################################################
# Purpose: Create xyz files of icosahedral cluster consisting of 1 #
#          or two different atom types                             #
####################################################################
#                Elke Fasshauer                                    #
####################################################################

import numpy as np
import math

##################Input Variables ##################################
atcore = 'Xe' # atomtype of the core atoms
atouter = 'Ar' # atomtype of the outer shells

rcore =  2.16 # radius of core atoms 
router = 2.07 # radius of outer shell atoms

n_core = 2 #number of atoms for the longest edge
n_outer = 0

phi = (1 + math.sqrt(5))/2

####################### Functions used #############################

####################### Programme ##################################

testvec = np.array([2,5,3])
avec    = np.array([1,4,2.5])

sp      = np.vdot(testvec,avec)
norm    = np.linalg.norm(avec)

##################################

for i in (1,n_core):
    kante  = rcore * (i-1)

# Die Ecken des Ikosaeders der entsprechenden Groesse
    ecke1  = np.array([           0,    -kante/2, kante/2*phi])
    ecke2  = np.array([           0,     kante/2, kante/2*phi])
    ecke3  = np.array([           0,    -kante/2,-kante/2*phi])
    ecke4  = np.array([           0,     kante/2,-kante/2*phi])
    ecke5  = np.array([     kante/2,-kante/2*phi,           0])
    ecke6  = np.array([     kante/2, kante/2*phi,           0])
    ecke7  = np.array([    -kante/2,-kante/2*phi,           0])
    ecke8  = np.array([    -kante/2, kante/2*phi,           0])
    ecke9  = np.array([ kante/2*phi,           0,     kante/2])
    ecke10 = np.array([ kante/2*phi,           0,    -kante/2])
    ecke11 = np.array([-kante/2*phi,           0,     kante/2])
    ecke12 = np.array([-kante/2*phi,           0,    -kante/2])

    ecken = np.vstack((ecke1,ecke2,ecke3,ecke4,ecke5,ecke6,ecke7,\
                       ecke8,ecke9,ecke10,ecke11,ecke12))
    print ecken
