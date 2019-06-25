import musclebeachtools as mbt
import numpy as np 
import os
def bin_FR_mat(datadir, binsz = 40, qual = 1, rawdatadir = [], multi_probe = False, probeNumber = 1, start_block = 0, end_block = 1):
	''' INPUTS:
		datadir: t folder for desired clustering job
		binsz: desired binsize (ms)
		qual: what quality to include (if more than one enter as list)
		rawdatadir: data with sleep info if applicable
		multi_probe: True or False
		probeNumber: if multiprobe, which one
		start_block/end_block: if choosing a particular block from clustering 
		OUTPUT:
		matrix with columns as different timepoints and rows as different neurons
		'''
	os.chdir(datadir)
	qual_list = np.load('scrubbed_quality.npy')
	if np.size(qual) > 1:
		cell_list = []
		for q in qual:
			cell_list.append(np.where(qual_list == q)[0])
		cell_list = np.concatenate(cell_list)
	else:
		cell_list = np.where(qual_list == qual)[0]
	cell_list = np.sort(cell_list)

	filelist = mbt.makeFileList(datadir)
	cells = []
	for a,i in enumerate(cell_list):
	    print(a, ' ', i)
	    cells.append(mbt.neuron(datadir, datatype='npy', cell_idx=i, file_list = filelist, rawdatadir = rawdatadir, start_block = start_block, end_block = end_block))
	    
	#make a function that spikes up waking spike times and sleeping spike times
	spks = mbt.getspikes(cells, 0, 3600*24) # gets spiketimes
	data_T = mbt.spiketimes_to_spikewords(spks,0,3600*24,binsz,1) #binarizes spikes
	data = data_T.T #transposes matrix so each column is a time bin and each row is a cell
	np.save('cr_mat.npy',data)

	return data, cells
