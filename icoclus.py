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

# Definition of the surfaces
surf1  = np.array([ 1, 2,11])
surf2  = np.array([ 1, 2, 9])
surf3  = np.array([ 1, 7,11])
surf4  = np.array([ 1, 5, 7])
surf5  = np.array([ 1, 5, 9])
surf6  = np.array([ 2, 6, 8])
surf7  = np.array([ 2, 8,11])
surf8  = np.array([ 2, 6, 9])
surf9  = np.array([ 3, 4,10])
surf10 = np.array([ 3, 4,12])
surf11 = np.array([ 3, 5,10])
surf12 = np.array([ 3, 5, 7])
surf13 = np.array([ 3, 7,12])
surf14 = np.array([ 4, 6, 8])
surf15 = np.array([ 4, 6,10])
surf16 = np.array([ 4, 8,12])
surf17 = np.array([ 5, 9,10])
surf18 = np.array([ 6, 9,10])
surf19 = np.array([ 7,11,12])
surf20 = np.array([ 8,11,12])

surfaces = np.vstack((surf1,surf2,surf3,surf4,surf5,surf6,surf7,\
                      surf8,surf9,surf10,surf11,surf12,surf13,\
                      surf14,surf15,surf16,surf17,surf18,surf19,surf20))
print surfaces

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
    print ecken[11]
