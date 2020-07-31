

## ------------- Dont edit anything between here and next line! ----------


import numpy as np
import matplotlib.pyplot as plt
	
def autocorr_3d(xyz):
	speedup = 1 ## keep this VERY low! 100 gives worthless Far-IR. I strongly suggest keeping it at 1.
	datalen = len(xyz)
	outlen = int(datalen/2)
	
	## Create arrays to store the calculated datapoints.
	cfun = np.zeros(outlen)
	cfunx = np.zeros(outlen)
	cfuny = np.zeros(outlen)
	cfunz = np.zeros(outlen)
	cfunxpy = np.zeros(outlen)
	cfunxmy = np.zeros(outlen)
	cfunxpz = np.zeros(outlen)
	cfunypz = np.zeros(outlen)
	cfunxpypz = np.zeros(outlen)
	cfunxmypz = np.zeros(outlen)
	
	## this loops over the different time differences.
	for i in range(datapoints_to_calculate):
		if i%(progress_update*10)==0 or i in [0,progress_update]: print(i)
		
		## define these 6 so the remaining code is more legible. 
		## each of the variables is a 1D array. The 0's are the unshifted
		## datapoints, the 1's are the shifted datapoints
		x0, y0, z0 = xyz[:-(i+1),0], xyz[:-(i+1),1], xyz[:-(i+1),2]
		x1, y1, z1 = xyz[i:-1,0], xyz[i:-1,1], xyz[i:-1,2]
		devis = datalen-i
		
		## Actually calculating the ACF datapoints. Using np.dot() is faster
		## than manually looping over the datapoints.
		cfun[i] = (np.dot(x0,x1) + np.dot(y0,y1) + np.dot(z0,z1))/devis
		cfunx[i] = np.dot(x0,x1)/devis
		cfuny[i] = np.dot(y0,y1)/devis
		cfunz[i] = np.dot(z0,z1)/devis
		cfunxmy[i] = 0.5 * np.dot(x0-y0,x1-y1)/devis ## *0.5 for normalization!
		cfunxpy[i] = 0.5 * np.dot(x0+y0,x1+y1)/devis
		cfunxpz[i] = 0.5 * np.dot(x0+z0,x1+z1)/devis
		cfunypz[i] = 0.5 * np.dot(y0+z0,y1+z1)/devis
		cfunypz[i] = 0.5 * np.dot(y0+z0,y1+z1)/devis
		cfunxpypz[i] = 0.33 * np.dot(x0+y0+z0,x1+y1+z1)/devis
		cfunxmypz[i] = 0.33 * np.dot(x0-y0+z0,x1-y1+z1)/devis
			
	return cfun, cfunx, cfuny, cfunz, cfunxmy, cfunxpy, cfunxpz, cfunypz, cfunxpypz, cfunxmypz
	

## -------------- Dont edit anything above this line! --------------

	
## the script only calculates the y-values of the ACF. The X values are
## taken from the gromacs generated ACF function.
testdat = np.genfromtxt("../XVGforIR/md_21_6_SnI.xvg",skip_header=17, skip_footer = 1)

## what data to calculate the ACF over (assumed file format is gromacs output
## dipole file)
input = np.genfromtxt("../XVGforIR/Mtot_raw_21_6_totcorr_All.xvg",skip_header=27)

## where to write the ACF to (and specify file name).
outname = "../md_21_6_"
outnamefin = "All"


## The gromacs files from this project are very lengthy. The full files could
## be read/used, but the resulting spectra are very noisy. 
datapoints_to_calculate = 3000


## -------------- Dont edit anything below this line! --------------


## create all the output files
outnamefin = outnamefin+".xvg"
fouttot = open(outname+"xyz_"+outnamefin,"w")
foutx = open(outname+"x_"+outnamefin,"w")
fouty = open(outname+"y_"+outnamefin,"w")
foutz = open(outname+"z_"+outnamefin,"w")
foutxmy = open(outname+"x-y_"+outnamefin,"w")
foutxpy = open(outname+"x+y_"+outnamefin,"w")
foutxpz = open(outname+"x+z_"+outnamefin,"w")
foutypz = open(outname+"y+z_"+outnamefin,"w")
foutxmypz = open(outname+"x-y+z_"+outnamefin,"w")
foutxpypz = open(outname+"x+y+z_"+outnamefin,"w")

## calculate this to have the command prompt indicate the progress through
## the calculation
progress_update = round(datapoints_to_calculate/100)

## split the input files into the different components in a suitable format
inpt = input[:,0]
inptot = input[:,4]
inpxyz = input[:,1:4]

## calculate the ACFs
outtot, outx, outy, outz, outxmy, outxpy, outxpz, outypz, outxpypz, outxmypz = autocorr_3d(inpxyz)

## This code generates files to be read by the FarIR.py script. That script
## can also read Gromacs-calculated ACFs, so this script has to generate
## files with the same format. Gromacs has 17 lines with info, which are
## not used anyways.
for i in range(17):
	fouttot.write("this is filler line " + str(i+1) + "\n")
	foutx.write("this is filler line " + str(i+1) + "\n")
	fouty.write("this is filler line " + str(i+1) + "\n")
	foutz.write("this is filler line " + str(i+1) + "\n")
	foutxmy.write("this is filler line " + str(i+1) + "\n")
	foutxpy.write("this is filler line " + str(i+1) + "\n")
	foutxpz.write("this is filler line " + str(i+1) + "\n")
	foutypz.write("this is filler line " + str(i+1) + "\n")
	foutxmypz.write("this is filler line " + str(i+1) + "\n")
	foutxpypz.write("this is filler line " + str(i+1) + "\n")
	
## write the actual data. Testdat provides the x-coordinates, the second part
## provides the y-coordinates.
for i in range(outtot.size):
	fouttot.write(str(testdat[i,0])+ "  " + str(outtot[i]) + "\n")
	foutx.write(str(testdat[i,0])+ "  " + str(outx[i]) + "\n")
	fouty.write(str(testdat[i,0])+ "  " + str(outy[i]) + "\n")
	foutz.write(str(testdat[i,0])+ "  " + str(outz[i]) + "\n")
	foutxmy.write(str(testdat[i,0])+ "  " + str(outxmy[i]) + "\n")
	foutxpy.write(str(testdat[i,0])+ "  " + str(outxpy[i]) + "\n")
	foutxpz.write(str(testdat[i,0])+ "  " + str(outxpz[i]) + "\n")
	foutypz.write(str(testdat[i,0])+ "  " + str(outypz[i]) + "\n")
	foutxmypz.write(str(testdat[i,0])+ "  " + str(outxmypz[i]) + "\n")
	foutxpypz.write(str(testdat[i,0])+ "  " + str(outxpypz[i]) + "\n")
	

