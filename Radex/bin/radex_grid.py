#! /usr/bin/python
#
import math
import os
import time
#
# Run a series of Radex models to estimate temperature & density
# from observed ratios of HCO+ 1-0/3-2, 3-2/4-3, and 4-3/7-6 lines.
#
# Floris van der Tak, 27oct05; revised 28jun06
#
# Grid boundaries
#
tmin = 10.0  # minimum kinetic temperature (K)
tmax = 200.0 # maximum kinetic temperature (K)
nmin = 1e3   # minimum H2 density (cm^-3)
nmax = 1e8   # maximum H2 density (cm^-3)
#
# Parameters to keep constant
#
tbg   = 2.73 # background radiation temperature
cdmol = 1e12 # low enough to be optically thin
dv    = 1.0  # line width (km/s)
#
# Numerical parameters
#
ntemp = 30  # number of temperature points
ndens = 30  # number of density points
bw    = 0.001 # "bandwidth": free spectral range around line
#
# No user changes needed below this point.
#
def write_input(infile,tkin,nh2):
    infile.write(mole+'.dat\n')
    infile.write('radex.out\n')
    infile.write(str(flow*(1-bw))+' '+str(fupp/(1-bw))+'\n')
    infile.write(str(tkin)+'\n')
    infile.write('1\n')
    infile.write('H2\n')
    infile.write(str(nh2)+'\n')
    infile.write(str(tbg)+'\n')
    infile.write(str(cdmol)+'\n')
    infile.write(str(dv)+'\n')

def read_radex(outfile):
    line  = outfile.readline()
    words = line.split()
    while (words[1] != "T(kin)"):
        line  = outfile.readline()
        words = line.split()
    temp  = float(words[-1])
    line  = outfile.readline()
    words = line.split()
    dens  = float(words[-1])
    while (words[-1] != "FLUX"):
        line  = outfile.readline()
        words = line.split()
    line  = outfile.readline()
    line  = outfile.readline()
    words = line.split()
    ftmp  = float(words[4])
    while ((ftmp < flow*(1-bw)) or (ftmp > flow/(1-bw))):
        line  = outfile.readline()
        words = line.split()
        ftmp  = float(words[4])
    low   = float(words[-2])
    line  = outfile.readline()
    words = line.split()
    ftmp  = float(words[4])
    while ((ftmp < fupp*(1-bw)) or (ftmp > fupp/(1-bw))):
        line  = outfile.readline()
        words = line.split()
        ftmp  = float(words[4])
    upp   = float(words[-2])
    ratio = low / upp
    return temp,dens,ratio
 
# Begin of main program

start = time.time()

mole = 'hco+'
acts = ([89.2,267.6,'1-0_3-2.dat'],[267.6,356.7,'3-2_4-3.dat'],[356.7,624.2,'4-3_7-6.dat'])

for act in acts:
    flow = act[0]
    fupp = act[1]
    gfil = act[2]
    infile = open('radex.inp','w')
    print "Starting ",gfil

    for itemp in range(ntemp+1):
        for idens in range(ndens+1):

            temp = tmin*((tmax/tmin)**(float(itemp)/ntemp))
            dens = nmin*((nmax/nmin)**(float(idens)/ndens))
            write_input(infile,temp,dens)
            if (itemp == ntemp and idens == ndens):
                infile.write('0\n')
                infile.close()
            else:
                infile.write('1\n')

    os.system('radex < radex.inp > /dev/null')
    grid = open(gfil,'w')
    fmt  = '%10.3e %10.3e %10.3e \n'

    outfile  = open('radex.out')

    rmin = 100
    rmax = 0.1

    for itemp in range(ntemp+1):
        for idens in range(ndens+1):

            temp = tmin*((tmax/tmin)**(float(itemp)/ntemp))
            dens = nmin*((nmax/nmin)**(float(idens)/ndens))
            temp,dens,ratio = read_radex(outfile)

            if (ratio > 0.0):
                if (ratio < rmin):
                    rmin = ratio
                if (ratio > rmax):
                    rmax = ratio
            grid.write(fmt %(temp, math.log10(dens), ratio))

    grid.close()
    outfile.close()

    print "Min, max:", rmin, rmax

stop = time.time()
dure = stop - start
print "Run time = ",dure, "seconds"
