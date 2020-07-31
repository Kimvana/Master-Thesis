

## ------------- Dont edit anything between here and next line! ----------


import numpy as np
import matplotlib.pyplot as plt
import math

## Calculating the actual IR spectrum. 
def fourtrans(data_in):
	## do the Fourier transform
	func_tempx = np.fft.rfftfreq(desired_ndatapoints,delta_t)
	func_tempy = np.real(np.fft.rfft(data_in[:,1]))
	## subtracting the last term from all values makes sure that the function
	## approaches zero at larger values
	func_tempy -= func_tempy[-1]
	func_IR = np.transpose(np.stack((func_tempx,func_tempy)))
	
	## each y-coordinate has to be multiplied with the square of the frequency
	func_IR[:,1]*= func_IR[:,0]**2
	
	## scaling the x-axis so the frequencies are correct
	func_IR[:,0]*= IR_x_scale_factor
	
	## terminating the series to only contain the datapoints needed for
	## drawing the graph (only interested in low-wavenumber part of spectrum)
	func_delta_cm=func_IR[1,0]-func_IR[0,0]
	needed_datapoints = 1+ math.ceil(desired_cm/func_delta_cm)
	func_IR = func_IR[:needed_datapoints,:]

	return func_IR

## Plot all the graphs with custom colours
def make_IRc(varname,data_in,make_IR_label,color,boxsize):
	globals()[varname] = np.genfromtxt(folder_sample+data_in+".xvg",skip_header=17, max_rows=desired_ndatapoints)
	globals()["delta_t"] = globals()[varname][1,0]-globals()[varname][0,0]
	globals()[varname+"_IR"] = fourtrans(globals()[varname])
	globals()[varname+"_IR"][:,1] /= boxsize
	plt.plot(globals()[varname+"_IR"][:,0], globals()[varname+"_IR"][:,1], label = make_IR_label, c=color)

## Plot all the graphs with default colours
def make_IR(varname,data_in,make_IR_label):
	globals()[varname] = np.genfromtxt(folder_sample+data_in+".xvg",skip_header=17, max_rows=desired_ndatapoints)
	globals()["delta_t"] = globals()[varname][1,0]-globals()[varname][0,0]
	globals()[varname+"_IR"] = fourtrans(globals()[varname])
	globals()[varname+"_IR"][:,1] /= boxsize
	plt.plot(globals()[varname+"_IR"][:,0], globals()[varname+"_IR"][:,1], label = make_IR_label)


## -------------- Dont edit anything above this line! --------------


## amount of data points to read from the ACF files. Those from gromacs
## contain many datapoints. The ones from 'ACF calculator.py' only contain
## a few, with appended zeros. The value used here should not exceed the one
## used in 'ACF calculator.py'. 3000 seems to give good results.
desired_ndatapoints = 3000

## highest wavenumber to plot in the graph.
desired_cm = 500

## 33.35641 if timestep in gromacs is 0.5 fs (saved every 5 fs)
IR_x_scale_factor = 33.35641

## Export the data of the plot (for use in e.g. origin)
folder_sample = "../XVGforIR/md_"
output = open("../farIR.txt","w")

## The file names of data to plot
input_data=["21_6_All","25_6_All","26_6_All","7_6_nonorm_All","27_6_All","29_6_All"]

## The names of the datasets to be used in the legend
graph_label = ["0 FA","1 FA","2 FA","3 FA","4 FA","Bulk"]

## The colours the graphs should have.
custom_colors = True ## set to false if default are desired
colors=["green","blue","purple","red","orange","black"]

## Axis names
Xaxis_label = "frequency (cm-1)"

## Plot title 
plot_title = "FASnI3 (nBA) Far-IR spectrum for different directions (nBA only)"

## The desired scaling factors. Define as a list of 1s when none are desired
boxsizes = [29.859,27.899,18.4766,22.184,26.599,16.0094]


## -------------- Dont edit anything below this line! --------------


## behind-the-scenes variable names (in case special changes must be made
## to a single dataset). The first one, "zero", is being used!
variable_names=["zero","one","two","three","four","hh","hi","hg","hl"]
ndata = len(input_data)
variable_names = variable_names[:ndata]

## plt doesn't clean up after itself... just housekeeping.
plt.close('all')

## instruct plt to plot the data
for i, e in enumerate(variable_names):
	if custom_colors: make_IRc(e, input_data[i], graph_label[i], colors[i], boxsizes[i])
	else: make_IR(e, input_data[i], graph_label[i] , boxsizes[i])

## write the data for use in other software
for i in range(zero_IR.shape[0]):
	output.write(str(zero_IR[i,0]))
	for varname in variable_names:
		output.write("\t"+str(globals()[varname+"_IR"][i,1]))
	output.write("\n")

## show the plot	
plt.xlim(0, desired_cm)
plt.xlabel(Xaxis_label)
plt.title(plot_title)
plt.legend()
plt.show()

## plt doesn't clean up after itself... just housekeeping.
plt.close()
plt.close('all')
print("quitting!")
quit()

