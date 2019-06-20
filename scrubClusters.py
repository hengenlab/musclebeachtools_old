import numpy as np
import os
import glob
import musclebeachtools as mb

def scrubClusters(datafile, rawdatadir, multi_probe, probeNumber):
	"""
	runs through all the cells in a dataset and calls the checkqual function so you can set the quality 

	start_cell is the first cell you scrub, it runs from the start_cell to the end 
	"""
	os.chdir(datafile)
	file_list=mb.makeFileList(datafile, rawdatadir, multi_probe, probeNumber=probeNumber)
	clusters=np.unique(file_list[0])
	num_cells = len(clusters)

	sq = np.sort(glob.glob("*scrubbed*.npy"))
	start_cell = 0
	if len(sq) > 0:
		print ("at least one of these cells has been scrubbed before")
		sq_array = file_list[6]
		first_unscrubbed = np.where(np.isnan(sq_array))[0]
		print("the first unscrubbed cell is index {}".format(first_unscrubbed[0]))
		start_cell= int(input("what cell index would you like to start scrubbing at? From 0 - {}".format(num_cells)))

	cells = []
	for c in np.arange(start_cell, num_cells):
		cell=mb.neuron(datadir =datafile, rawdatadir = rawdatadir, clust_idx = clusters[c], file_list = file_list)
		cells.append(cell)

	for c in cells:
		c.checkqual(scrub_cell = True)

	#prints the quality stats for the data set once you get to the end
	#basically how many of each quality cell there was
	quals=np.zeros(4)
	#need to restate this line in case there wasn't a scrubbed quality array to begin with
	sq = np.sort(glob.glob("*scrubbed*.npy"))
	scrubbed_qual_array = np.load(sq[0])
	for q in range(len(scrubbed_qual_array)):
		quals[int(scrubbed_qual_array[q])] = quals[int(scrubbed_qual_array[q])] +1

	print("Quality Statistics:\n1: {} \n2: {} \n3: {}".format(quals[1], quals[2], quals[3]))




	
	










