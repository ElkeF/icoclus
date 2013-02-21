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

n_core = 3 #number of atoms for the longest edge
n_outer = 0

################## Definitions #####################################

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

central = np.array([0,0,0])
coords  = central

####################### Functions used #############################

def vec2str(vec):
    return "   ".join([str(v) for v in vec])

####################### Programme ##################################


##################################
#####  Build the core icosahedra #
##################################
for i in range (2,n_core+1):
    kante  = 2* rcore * (i-1)

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

    latest = ecken
#    print latest

    for j in range (0,19):
        vec1 = ecken[surfaces[j,0] -1]
        vec2 = ecken[surfaces[j,1] -1]
        vec3 = ecken[surfaces[j,2] -1]

#        print ' '.join(map(str, vec1))
#        print ' '.join(map(str, vec2))
#        print ' '.join(map(str, vec3))
        normkante = (vec2-vec1) / np.linalg.norm(vec2-vec1)
        normlauf  = (vec3-vec2) / np.linalg.norm(vec3-vec2)

        if (i > 2):

            for k in range (1,i):
                kantatom = vec1 + (k * normkante * 2 * rcore)
                latest = np.vstack((latest,kantatom))
#                print vec1
#                print kantatom
                
                for l in range (1,k+1):
                    flatom = kantatom + l * normlauf * 2 * rcore
                    latest = np.vstack((latest,flatom))
                    print kantatom
                    print normlauf
                    print flatom

#    print latest
    coords = np.vstack((coords,latest))

#########################################
# Write Output
#########################################




xyz_1st  = [vec2str(coord) for coord in coords]

lines_1st = []
for coord in xyz_1st:
    line = '%s    %s' %(atcore,coord)
    lines_1st.append(line)
print_1st = '\n'.join(lines_1st)


outlist = [print_1st]

#print outlist
outlines = '\n'.join(outlist)
#print outlines

outfile = open("%s%s_ico.xyz" %(atouter,atcore), mode="w")
outfile.writelines(outlines)

