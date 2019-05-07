import numpy as np
import os
import glob
import musclebeachtools as mb

def scrubClusters(datafile, start_cell=0):
	"""
	runs through all the cells in a dataset and calls the checkqual function so you can set the quality 

	start_cell is the first cell you scrub, it runs from the start_cell to the end 
	"""
	os.chdir(datafile)
	file_list=makeFileList(datafile)
	clusters=np.unique(file_list[0])
	num_cells = len(clusters)

	sq = np.sort(glob.glob("*scrubbed*.npy"))

	if len(sq) > 0:
		print ("at least one of these cells has been scrubbed before")
		sq_array = file_list[6]
		first_unscrubbed = sq_array.index(np.nan)
		print("the first unscrubbed cell is index {}".format(first_unscrubbed))
		start_cell= input("what cell index would you like to start scrubbing at? From 0 - {}".format(num_cells))


	for i in range(start_cell, num_cells):
		cell=mb.neuron(datafile=datafile, cell_idx=i, file_list=file_list)
		cell.checkqual()

	#prints the quality stats for the data set once you get to the end
	#basically how many of each quality cell there was
	quals=np.zeros(4)
	#need to restate this line in case there wasn't a scrubbed quality array to begin with
	sq = np.sort(glob.glob("*scrubbed*.npy"))
	scrubbed_qual_array = np.load(sq[0])
	for q in range(len(scrubbed_qual_array)):
		quals[int(scrubbed_qual_array[q])] = quals[int(scrubbed_qual_array[q])] +1

	print("Quality Statistics:\n1: {} \n2: {} \n3: {}".format(quals[1], quals[2], quals[3]))



"""
sorts and loads all the files in a directory in the same way that neuon_class does
then puts each loaded array into a list and returns that list

this way if we're running through mulitple cells at a time it'll only need to load once
#current format is:
#index:      data:
#0				curr_clusts
#1				curr_spikes
#2				peak_chs
#3				waveform
#4				amplitude
#5				auto_qual_array		
#6				scrubbed_qual_array
#7				scrubbed quality file name
	"""

def makeFileList(datafile, silicon=False, start_day=0, end_day=1):

	try:
		os.chdir(datafile)
	except FileNotFoundError:
		print("*** Data File does not exist *** check the path")
		return

	file_list=[]
	length = np.zeros(end_day)

		
   #SORTS DATA FILES
	if(silicon):
		ch = input("What probe would you like to look at?")
		f="*chg_"+str(ch)+"*"
		channelFiles=np.sort(glob.glob(f))
		#sorts spikes and clusters
		spikefiles = [channelFiles[i] for i in range(len(channelFiles)) if channelFiles[i] in np.sort(glob.glob("*spike_times*.npy"))]
		#print("spike_files: ", spikefiles)
		clustfiles = [channelFiles[i] for i in range(len(channelFiles)) if channelFiles[i] in np.sort(glob.glob("*spike_clusters*.npy"))]
		#sorts any peak channel files found in folder
		peakfiles = [channelFiles[i] for i in range(len(channelFiles)) if (channelFiles[i] in np.sort(glob.glob("*peakchannel*.npy")) or channelFiles[i] in np.sort(glob.glob("*max_channel*.npy")))]
		#sorts wavefiles in two forms, named "waveform" or "templates"
		wavefiles=[channelFiles[i] for i in range(len(channelFiles)) if channelFiles[i] in np.sort(glob.glob("*waveform*.npy"))]
		templates_all=[channelFiles[i] for i in range(len(channelFiles)) if channelFiles[i] in np.sort(glob.glob("*templates*.npy"))]
		#since there are multiple files that have "templates" at the end this pulls out only the ones we want
		templates_wf=[fn for fn in templates_all if fn not in glob.glob("*spike*.npy") and fn not in glob.glob("*similar*.npy") and fn not in glob.glob("*number_of_*.npy") and fn not in glob.glob("*templates_in_clust.npy")]
		#checks for amplitude files
		amplitude_files = [channelFiles[i] for i in range(len(channelFiles)) if channelFiles[i] in np.sort(glob.glob("*amplitudes*.npy"))]
		#this checks for an automated quality array from the clustering algorithm
		aq = [channelFiles[i] for i in range(len(channelFiles)) if channelFiles[i] in np.sort(glob.glob("*qual*.npy"))]
		#looks for scrubbed quality and loads if possible
		sq = [channelFiles[i] for i in range(len(channelFiles)) if channelFiles[i] in np.sort(glob.glob("*scrubbed*.npy"))]
	#pulls data if not silicon
	else:
		#sorts spikes and clusters
		spikefiles = np.sort(glob.glob("*spike_times*.npy"))
		clustfiles = np.sort(glob.glob("*spike_clusters*.npy"))
		#sorts any peak channel files found in folder
		peakfiles = np.sort(glob.glob("*max_channel*.npy"))
		if len(peakfiles)==0:
			peakfiles=np.sort(glob.glob("*peakchannel*.npy"))
		#peakfiles.extend(np.sort(glob.glob("*max_channel*.npy")))
		#sorts wavefiles in two forms, named "waveform" or "templates"
		wavefiles = np.sort(glob.glob("*waveform*.npy"))
		templates_all=np.sort(glob.glob("*template*.npy"))
		#since there are multiple files that have "templates" at the end this pulls out only the ones we want
		templates_wf=[fn for fn in templates_all if fn not in glob.glob("*spike*.npy") and fn not in glob.glob("*similar*.npy") and fn not in glob.glob("*number_of_*.npy") and fn not in glob.glob("*templates_in_clust.npy")]
		#checks for amplitude files
		amplitude_files = np.sort(glob.glob("*amplitudes*.npy"))
		#this checks for an automated quality array from the clustering algorithm
		aq = np.sort(glob.glob("*qual*.npy"))
		#looks for scrubbed quality and loads if possible
		sq = np.sort(glob.glob("*scrubbed*.npy"))

	has_peak_files = not len(peakfiles)==0
	has_twf = not (len(wavefiles)==0 and len(templates_wf)==0)
	has_aqual = not len(aq)==0
	has_squal = not len(sq)==0

	try:
		print("Loading files...")
		curr_clust = [np.load(clustfiles[i]) for i in range(start_day, end_day)]
		curr_spikes = [np.load(spikefiles[i])+length[i] for i in range(start_day, end_day)]
		file_list.append(curr_clust)
		file_list.append(curr_spikes)
		if has_peak_files:
			peak_ch = np.concatenate([np.load(peakfiles[i])for i in range(start_day, end_day)])
			file_list.append(peak_ch)
		else:
			file_list.append([])
		if has_twf and len(wavefiles)==0:
			w=np.load(templates_wf[0])
			file_list.append(w)
		elif has_twf:
			w=np.load(wavefiles[0])
			file_list.append(w)
		else:
			file_list.append([])
		if len(amplitude_files)>0:
			amps=np.load(amplitude_files[0])
			file_list.append(amps)
			#self.amplitudes=np.load(amplitude_files[0])
		else:
			file_list.append([])
		if has_aqual:
			aqa=np.load(aq[0])[peak_ch != 0]
			file_list.append(aqa)
		else:
			file_list.append([])
			#self.auto_qual_array = np.load(aq[0])[peak_ch != 0]
		if has_squal:
			sqa=np.load(sq[0])
			file_list.append(sqa)
			#self.scrubbed_qual_array = np.load(sq[0])
			file_list.append(sq[0])
		else:
			file_list.append([])
	except IndexError:
		print("files do not exist for that day range")

	return file_list

	
	










