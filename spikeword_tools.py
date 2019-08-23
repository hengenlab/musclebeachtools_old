#SPIKEWORDS_9418

import musclebeachtools as mbt
import numpy as np
import matplotlib.pyplot as plt
import os
import datetime
import seaborn as sns
import autoreload
import pandas as pd
import h5py
import re
import csv
import time
import pickle

def bintoint(row):
    sw_s = "".join(row.astype(str))
    sw_s = sw_s.replace('.','')
    #sw = int(sw_s,2)
    return(sw_s)

def readhdf5datafile(datdir, animalnames):
    print('LOADING CELLS')
    fn_end = '.hdf5'
    datfile = [datdir+i+fn_end for i in animalnames]

    # negpos time
    RSUxlims = (0.4, 1)
    FSxlims = (0, 0.4)

    # tailslope
    RSUylims = (0, 400)
    FSylims = (-300, 0.0)

    # Index the neurons of interest
    allcells = mbt.indexdata(datfile, qual=(1, 2), contstat=1, dep=1)  # negpos=FSxlims, slope=FSylims

    global all_nrn
    all_nrn = []
    for i in allcells[0][0]:
    	var     = [mbt.neuron((datfile[0]), i)]
    	all_nrn.append(var)
    
    all_nrn     = sum(all_nrn[:],[])
    all_net 	= mbt.nrn_network( all_nrn ) 
    return all_nrn

def getspikes(neuron_list,startime,stoptime):
	n_cells = np.shape(neuron_list)[0]
	spiketimes_allcells = list()
	for i in np.arange(n_cells):
		print('Getting spiketimes for cell ' + str(i))
		#get all spiketimes for cell
		spiketimes = neuron_list[i].time
		spiketimes = spiketimes[(spiketimes>startime)&(spiketimes<stoptime)]
		spiketimes_allcells.append(spiketimes)
	return (spiketimes_allcells)

def shuffle_spikes(spiketimes,starttime,stoptime):
	shuffled_spikes = []
	for i in np.arange(len(spiketimes)):
		spikestoshuffle = spiketimes[i]
		n_spikes = len(spiketimes[i])
		randspikes = np.random.uniform(size = n_spikes, low = starttime, high = stoptime)
		shuffled_spikes.append(np.sort(randspikes))
	return shuffled_spikes

def spiketimes_to_spikewords(spiketimes,startime,stoptime,binsize,binarize): 
    ### ARGUMENTS
    #spiketimes - list of spiketime arrays where each element of list is for a different neuron (e.g. the output of getspikes())
    #startime,stoptime - bounds of time to convert to spikewords (in seconds)
    #binsize - size of bin in milliseconds
    #binarize - '1' to convert to 0's and 1's, '0' to keep as histogram counts
    ### RETURNS
    #array of spikewords with each column as a cell and each rows as time in bins 
    
    #get sec_time to bin conversion factor
    #startime in bins
    startime_ms = startime * 1000
    stoptime_ms = stoptime * 1000
    binrange = np.arange(start = startime_ms,stop = stoptime_ms+1, step = binsize)
    n_cells = len(spiketimes)

    spikewords_array = np.zeros([n_cells,binrange.shape[0]-1])
    for i in range(n_cells):
        spiketimes_cell = np.asarray(spiketimes)[i] * 1000. #spiketimes in seconds * 1000 msec/sec
        counts, bins = np.histogram(spiketimes_cell,bins = binrange)
        if binarize == 1:
	        #binarize the counts
	        counts[counts>0] = 1
        # print(counts.astype(np.int))
        spikewords_array[i,:] = counts
    return(spikewords_array.astype(np.int8).T)

def spikeword_distributions(spikeword_array,save_output,spikewordfn):
	print('CONVERTING SPIKE WORDS, COUNTING SPIKE WORDS')
	print("shape(sw) ",np.shape(spikeword_array))
	wordslist = np.apply_along_axis(bintoint, axis=1, arr = spikeword_array)
	unique, counts = np.unique(wordslist, return_counts=True)
	dict(zip(unique,counts))
	pdf = {'spikeword': unique,'counts': counts}
	swdf = pd.DataFrame(data=pdf)
	#swdf = swdf.sort_values('counts',axis=0,ascending=False)
	#swdf = swdf.iloc[:,1:3]
	if save_output==1:
		swdf.to_csv(spikewordfn)
	# 72618: added line to generate tseries of spikewords
	tsdf = pd.DataFrame(data=wordslist)
	if save_output==1:
		tsdf.to_csv(spikewordfn[:-4] + '_tseries.csv')
	return swdf, tsdf


###Not sure this one works yet.  Need to test it
def spikeword_count_FRcorrection(neuron_list,spiketimes,startime,stoptime,binsize,swdf):
	n_cells = np.shape(neuron_list)[0]

	#firing rates of each cell
	fr = np.zeros(n_cells)
	for i in np.arange(n_cells):
		print('Getting spiketimes for cell ' + str(i))
		#fr[i] = np.shape(spiketimes[i])[0]/(np.max(spiketimes[i])-np.min(spiketimes[i])) #in seconds
		fr[i] = np.shape(spiketimes[i])[0]/(stoptime-startime) #in seconds

	#convert FR from Hz to spikes/bin
	# spikes/s * (1s/1000 ms) * (10 ms/1 tbin)
	netfrhits = (fr/1000)*binsize
	netfrhits[netfrhits>1]=1 #if probability > 1 then it is inevitable
	netfrmisses = 1-netfrhits
	
	n_spkwrds = len(swdf['spikeword'])
	prob_sw=np.ones(n_spkwrds)
	for i in np.arange(n_spkwrds):
		spkwrdex = swdf['spikeword'][i]
		sw_current = np.array(list(spkwrdex),dtype=int)
		for j in np.arange(n_cells):
		    if sw_current[j] == 1:
		        prob_sw[i] = prob_sw[i] * netfrhits[j]
		    elif sw_current[j] == 0:
		        prob_sw[i] = prob_sw[i] * netfrmisses[j]

	n_bins = (stoptime-startime)*1000*(1/binsize)
	sw_correction = prob_sw * n_bins
	sw_correction = np.round(sw_correction)
	sw_correction = sw_correction.astype(int)
	#make: counts observed minus counts expected
	# swdf['counts_adj'] = swdf['counts']-sw_correction
	return sw_correction
