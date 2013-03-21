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
atcore = 'Ar' # atomtype of the core atoms
atouter = 'Ne' # atomtype of the outer shells

rcore =  1.88 # radius of core atoms 
router = 1.54 # radius of outer shell atoms

n_core = 8 #number of atoms for the longest edge
n_outer = raw_input('How many layers of atoms do you want to have? ')
n_outer = int(n_outer)
#n_outer = 1

################## Definitions #####################################

phi = (1 + math.sqrt(5))/2
scale = 2.0 / math.sqrt(1+phi**2) # passt fuer einheitliche Atome

thres = 1e-10

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

#def unique_rows(a):
#    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
#    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))

def unique_rows(a):
#    order = np.lexsort(a.T)
#    a = a[order]
    a = np.around(a,decimals=10)
    a = a[np.lexsort(a.T)]
#    a = a[a[:,2].argsort()]
#    a = a[a[:,1].argsort()]
#    a = a[a[:,0].argsort()]
#    print a
    diff = np.diff(a, axis=0)
    ui = np.ones(len(a), 'bool')
    ui[1:] = (diff > thres).any(axis=1) 
    return a[ui]

####################### Programme ##################################


##################################
#####  Build the core icosahedra #
##################################
for i in range (2,n_core+1):
    kante  = 2* rcore * (i-1) * scale

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

    for j in range (0,20):
        vec1 = ecken[surfaces[j,0] -1]
        vec2 = ecken[surfaces[j,1] -1]
        vec3 = ecken[surfaces[j,2] -1]

        #print j+1
        #print surfaces[j,0], surfaces[j,1], surfaces[j,2]

#        print ' '.join(map(str, vec1))
#        print ' '.join(map(str, vec2))
#        print ' '.join(map(str, vec3))
        normkante = (vec2-vec1) / np.linalg.norm(vec2-vec1)
        normlauf  = (vec3-vec2) / np.linalg.norm(vec3-vec2)

        if (i > 2):

            for k in range (1,i):
                kantatom = vec1 + (k * normkante * 2 * rcore * scale)
                latest = np.vstack((latest,kantatom))
#                print vec1
#                print kantatom
                
                for l in range (1,k+1):
                    flatom = kantatom + l * normlauf * 2 * rcore * scale
                    latest = np.vstack((latest,flatom))
                    #print kantatom


# Entferne Duplikate innerhalb der Liste
    unique = unique_rows(latest)

# vereine die Koordianten der letzten Schicht mit allen anderen
    coords = np.vstack((coords,unique))



##############################################
##### Outer Atom Layers ######################
##############################################
for i in range (1,n_outer+1):
    kante  = (rcore * (2*n_core-1) + (2*i-1) *router) * scale

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

    for j in range (0,20):
        vec1 = ecken[surfaces[j,0] -1]
        vec2 = ecken[surfaces[j,1] -1]
        vec3 = ecken[surfaces[j,2] -1]

        #print j+1
        #print surfaces[j,0], surfaces[j,1], surfaces[j,2]

#        print ' '.join(map(str, vec1))
#        print ' '.join(map(str, vec2))
#        print ' '.join(map(str, vec3))
        normkante = (vec2-vec1) / np.linalg.norm(vec2-vec1)
        normlauf  = (vec3-vec2) / np.linalg.norm(vec3-vec2)
        atdist    = kante / (n_core+i-1)

        for k in range (1,n_core + i):
            kantatom = vec1 + (k * normkante * atdist)
            latest = np.vstack((latest,kantatom))
#            print vec1
#            print kantatom

            for l in range (1,k+1):
                flatom = kantatom + l * normlauf * atdist
                latest = np.vstack((latest,flatom))
                #print kantatom


# Entferne Duplikate innerhalb der Liste
    unique = unique_rows(latest)

    if i == 1:
        coords2nd = unique
    else:
# vereine die Koordianten der letzten Schicht mit allen anderen
        coords2nd = np.vstack((coords2nd,unique))





#########################################
# Write Output
#########################################




xyz_1st  = [vec2str(coord) for coord in coords]
xyz_2nd  = [vec2str(coord) for coord in coords2nd]

lines_1st = []
for coord in xyz_1st:
    line = '%s    %s' %(atcore,coord)
    lines_1st.append(line)
print_1st = '\n'.join(lines_1st)

lines_2nd = []
for coord in xyz_2nd:
    line = '%s    %s' %(atouter,coord)
    lines_2nd.append(line)
print_2nd = '\n'.join(lines_2nd)


no_core_atoms = len(lines_1st)
no_outer_atoms = len(lines_2nd)
no_atoms = no_core_atoms + no_outer_atoms


outlist = [print_1st,print_2nd]

#print outlist
outlines = '\n'.join(outlist)
#print outlines

outfile = open("%s%s_ico_%d.xyz" %(atouter,atcore,no_atoms), mode="w")
outfile.writelines(outlines)

