import numpy as np
import os
import glob
import musclebeachtools as mb

def scrubClusters(datafile, start_cell=0):
	os.chdir(datafile)
	clustfiles = np.sort(glob.glob("*spike_clusters*.npy"))
	curr_clusts=[np.load(clustfiles[i]) for i in range(1)]
	clusters=np.unique(curr_clusts)
	num_cells = len(clusters)
	sq = np.sort(glob.glob("*scrubbed*.npy"))

	if len(sq) > 0:
		print ("at least one of these cells has been scrubbed before")

	for i in range(start_cell, num_cells):
		cell=mb.neuron(datafile=datafile, cell_idx=i)
		cell.checkqual(save_update=True)

	#add stats here
	quals=np.zeros(4)
	sq = np.sort(glob.glob("*scrubbed*.npy"))
	scrubbed_qual_array = np.load(sq[0])
	for q in range(len(scrubbed_qual_array)):
		quals[int(scrubbed_qual_array[q])] = quals[int(scrubbed_qual_array[q])] +1

	print("\n1: {} \n2: {} \n 3: {}".format(quals[1], quals[2], quals[3]))