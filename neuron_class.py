import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import h5py
import time
import seaborn as sns
import matplotlib
import glob
import math
import neuraltoolkit.ntk_ecube as ntk
import datetime as dt
#matplotlib.use('TkAgg')
import pdb

class neuron(object):
    '''Create instances (objects) of the class neuron '''
    def __init__(self, datafile, datatype='npy', cell_idx = 0, start_day=0, end_day=1, silicon=False,probenumber=1,fs=25000,file_list=[]):
        if datatype == 'npy':

            print('You are using WashU data')
            # fs = 25000
            print('Sampling rate {:10.2f}'.format(fs))
            print('working on neuron '+str(cell_idx)+'\n')

            # going to folder that contains processed data if it exists
            try:
                os.chdir(datafile)
            except FileNotFoundError:
                print("*** Data File does not exist *** check the path")
                #return

            if len(file_list)>0:
                print("using file list")
                has_peak_files  = not len(file_list[2]) == 0
                has_twf         = not len(file_list[3]) == 0
                has_aqual       = not len(file_list[5]) == 0
                has_squal       = not len(file_list[6]) == 0
                #has_amp= not len(file_list[4])==0

                curr_clust  = file_list[0]
                curr_spikes = file_list[1]
                if(has_peak_files):
                    peak_ch = file_list[2]
                    if np.size(np.where(peak_ch == -1)[0]) == 0:
                        print('You are using older data, changing indexing')
                        peak_ch     = peak_ch-1
                        curr_clust  = curr_clust[0]-1
                        curr_spikes = curr_spikes[0]-1
                if(has_twf):
                    w = file_list[3]
                #if(has_amp):
                #    self.amplitudes= file_list[4]
                if(has_aqual):
                    self.auto_qual_array=file_list[5]
                if(has_squal):
                    self.scrubbed_qual_array=file_list[6]
                    self._scqu_file=file_list[7]
                length = np.zeros(end_day)
                spikefiles = np.sort(glob.glob("*spike_times*.npy"))

            else:

               #SORTS DATA FILES
                if(silicon):
                    #ch  = input("What probe would you like to look at?")
                    ch = probenumber
                    f   = "*chg_"+str(ch)+"*"
                    channelFiles = np.sort(glob.glob(f))
                    #sorts spikes and clusters
                    spikefiles  = [channelFiles[i] for i in range(len(channelFiles)) if channelFiles[i] in np.sort(glob.glob("*spike_times*.npy"))]
                    #print("spike_files: ", spikefiles)
                    clustfiles  = [channelFiles[i] for i in range(len(channelFiles)) if channelFiles[i] in np.sort(glob.glob("*spike_clusters*.npy"))]
                    #sorts any peak channel files found in folder
                    peakfiles   = [channelFiles[i] for i in range(len(channelFiles)) if (channelFiles[i] in np.sort(glob.glob("*peakchannel*.npy")) or channelFiles[i] in np.sort(glob.glob("*max_channel*.npy")))]
                    #sorts wavefiles in two forms, named "waveform" or "templates"
                    wavefiles = [channelFiles[i] for i in range(len(channelFiles)) if channelFiles[i] in np.sort(glob.glob("*waveform*.npy"))]
                    templates_all = [channelFiles[i] for i in range(len(channelFiles)) if channelFiles[i] in np.sort(glob.glob("*templates*.npy"))]
                    #since there are multiple files that have "templates" at the end this pulls out only the ones we want
                    templates_wf = [fn for fn in templates_all if fn not in glob.glob("*spike*.npy") and fn not in glob.glob("*similar*.npy") and fn not in glob.glob("*number_of_*.npy") and fn not in glob.glob("*templates_in_clust.npy")]
                    #checks for amplitude files
                    #amplitude_files = [channelFiles[i] for i in range(len(channelFiles)) if channelFiles[i] in np.sort(glob.glob("*amplitudes*.npy"))]
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
                    #amplitude_files = np.sort(glob.glob("*amplitudes*.npy"))
                    #this checks for an automated quality array from the clustering algorithm
                    aq = np.sort(glob.glob("*qual*.npy"))
                    #looks for scrubbed quality and loads if possible
                    sq = np.sort(glob.glob("*scrubbed*.npy"))


                has_peak_files = not len(peakfiles)==0
                has_twf = not (len(wavefiles)==0 and len(templates_wf)==0)
                has_aqual = not len(aq)==0
                has_squal = not len(sq)==0

                #eventually length will be calculated differently but for now length is just 0
                #see original neuron_class for length calculation with keys
                length = np.zeros(end_day)

                #LOADS DATA
                try:
                    print("Loading files...")
                    curr_clust = [np.load(clustfiles[i]) for i in range(start_day, end_day)]
                    curr_spikes = [np.load(spikefiles[i])+length[i] for i in range(start_day, end_day)]
                    if has_peak_files:
                        peak_ch = np.concatenate([np.load(peakfiles[i])for i in range(start_day, end_day)])
                        if np.size(np.where(peak_ch == -1)[0]) == 0:
                            print('You are using older data, changing indexing')
                            peak_ch = peak_ch-1
                            curr_clust = curr_clust[0]-1
                            curr_spikes = curr_spikes[0]-1
                    if has_twf and len(wavefiles)==0:
                        w=np.load(templates_wf[0])
                    elif has_twf:
                        w=np.load(wavefiles[0])
                    # if len(amplitude_files)>0:
                    #     # pdb.set_trace()
                    #     tempamps = np.load(amplitude_files[0])

                    if has_aqual:
                        # self.auto_qual_array = np.load(aq[0])[peak_ch >= 0]
                        self.auto_qual_array = np.load(aq[0])
                    if has_squal:
                        self.scrubbed_qual_array = np.load(sq[0])
                        self._scqu_file = sq[0]
                except IndexError:
                    print("files do not exist for that day range")

            if end_day-start_day > 1:
                print("this cannot be done yet")
                #can't do this yet, see past code for an idea of how it will be done
            else:
                if np.shape(curr_clust)[0] == 1:
                    curr_clust = curr_clust[0]
                if np.shape(curr_spikes)[0] == 1:
                    curr_spikes = curr_spikes[0]
                #pulls out the unique cluster numbers
                self.unique_clusters = np.unique(curr_clust)
                clust_idx = self.unique_clusters[int(cell_idx)]

                #pulls out all indices of that cluster spiking based on cluster index and spike clusters
                spk_idx = np.where(curr_clust == clust_idx)[0]
                #loads all times at those indicies
                spiketimes = curr_spikes[spk_idx]/fs
                #amplitudes = [int(tempamps[i]) for i in spk_idx]
                #self.amplitudes =  np.asarray(amplitudes)
                #if there are peak channels this loads them into the instance variables
                #eventually this will be determined by day so that's why it's loaded here
                if has_peak_files:
                    peak_ch_no0 = np.delete(peak_ch,np.where(peak_ch == -1)[0],0)
                    self.peak_chans = peak_ch_no0[int(cell_idx)]

            if has_twf:
                w_real = np.delete(w,np.where(peak_ch == -1)[0],0)
                self.wf_temp = w_real[cell_idx]
                bottom      = np.argmin(self.wf_temp)
                top         = np.argmax(self.wf_temp[bottom:]) + bottom
                np_samples  = top - bottom

                seconds_per_sample  = 1/fs
                ms_per_sample       = seconds_per_sample*1e3
                self.neg_pos_time   = np_samples * ms_per_sample
                self.cell_idx = cell_idx
                if self.neg_pos_time >=0.4:
                    self.cell_type  = 'RSU'
                if self.neg_pos_time <0.4:
                    self.cell_type = 'FS'

            self._has_twf = has_twf
            self._has_squal = has_squal

            #sets the time instance variables to all the spike times
            self.time = np.concatenate(spiketimes)
            #sets on an off time based on first and last numbers in the spike time array
            self.onTime = np.array([self.time[0]])
            self.offTime = np.array([self.time[-1]])
            startTimeIndex = spikefiles[0].find("times_") + 6
            if startTimeIndex==-1:
                startTimeIndex = spikefiles[0].find("fs") + 2
            startTimeIndexEnd = spikefiles[0].find("-timee_")
            if startTimeIndexEnd ==-1:
                startTimeIndexEnd = spikefiles[0].find("-fs")
            clocktimeSIdx = spikefiles[0].find("_2019-") # this was hardcoded to 2018. changed to 2019. need to make flexible. KBH 5/11/19
            self.ecube_start_time = int(spikefiles[0][startTimeIndex : startTimeIndexEnd])
            self.clocktime_start_time = spikefiles[0][clocktimeSIdx: clocktimeSIdx+20]

            print("\nData set information:\nPeak Channels: {}\nWaveform Template: {}\nNumber of clusters: {}".format(has_peak_files, has_twf, len(self.unique_clusters)))


            #QUALITY
            if has_squal:
                last_idx=None
                for n in range(len(self.scrubbed_qual_array)):
                    if np.isnan(self.scrubbed_qual_array[n]):
                        last_idx=n
                        break
                if last_idx != None:
                    print("First unscrubbed cell index: ", last_idx)
                if np.isnan(self.scrubbed_qual_array[cell_idx]):
                    if has_aqual:
                        self.quality = self.auto_qual_array[int(clust_idx)]
                        print("\nScrubbed: NO")
                        print("Quality rating (automated): ", self.quality)
                else:
                    self.quality = self.scrubbed_qual_array[int(cell_idx)]
                    print('\nScrubbed: YES')
                    print("Quality set at: ", self.quality)
            else:
                if has_aqual:
                    self.quality = self.auto_qual_array[int(clust_idx)]
                    print("\nQuality rating (automated): ", self.quality)
                else:
                    print("\nNo automated or scrubbed quality rating, check the quality to set the quality manually")
                    self.quality = 0

            # #SLEEP-WAKE
            # fname = '{}*SleepStates*.npy'.format(rawdatadir)
            # sleepFiles = glob.glob(fname)
            # sleepFiles = np.sort(sleepFiles)
            # sleep_states = np.zeros(np.shape(self.time))
            # for i, t in enumerate(self.time):
            #     hr = math.ceil(t/3600)
            #     if hr > len(sleepFiles):
            #         print("Sleep states recorded through hour {}".format(hr-1))
            #         break
            #     sleeps = np.load(sleepFiles[hr-1])
            #     for idx, state in enumerate(sleeps):
            #         if t > (idx*4)+(900*(hr-1)) and t<(idx+1)*4 + (900*(hr-1)):
            #             sleep_states[i] = state
            #             break   
            # self.sleep_states = sleep_states            

        #Brandeis data
        else:
            print('You are using Brandeis data')
            print('Building neuron number {}'.format(cell_idx))
            f               = h5py.File(datafile,'r')
            try:
                self.animal     = np.array(f['neurons/neuron_'+str(cell_idx)+'/animal'][:], dtype=np.int8).tostring().decode("ascii")
            except:
                print('Stuck at line 19 in neuron_class.py')
                Tracer()()

            self.HDF5_tag   = (datafile,cell_idx)
            self.deprived   = f['neurons/neuron_'+str(cell_idx)+'/deprived'][0]
            self.channel    = np.int(f['neurons/neuron_'+str(cell_idx)+'/channel'][0])
            self.cont_stat  = f['neurons/neuron_'+str(cell_idx)+'/cont_stat'].value#bool(f['neurons/neuron_'+str(cell_idx)+'/cont_stat'][0])
            self.halfwidth  = f['neurons/neuron_'+str(cell_idx)+'/halfwidth'][0]
            self.idx        = f['neurons/neuron_'+str(cell_idx)+'/idx'][:]
            self.meantrace  = f['neurons/neuron_'+str(cell_idx)+'/meantrace'][:]
            self.scaledWF   = f['neurons/neuron_'+str(cell_idx)+'/scaledWF'][:]

            if np.size(f['neurons/neuron_'+str(cell_idx)+'/offTime'].shape) > 0:

                self.offTime    = f['neurons/neuron_'+str(cell_idx)+'/offTime'][:]
                self.onTime     = f['neurons/neuron_'+str(cell_idx)+'/onTime'][:]

            elif np.size(f['neurons/neuron_'+str(cell_idx)+'/offTime'].shape) == 0:
                self.offTime    = []
                self.onTime     = []

            try:
                self.qflag      = f['neurons/neuron_'+str(cell_idx)+'/qflag'][0]
                self.score      = f['neurons/neuron_'+str(cell_idx)+'/score'][:]
            except:
                pass

            self.tailSlope  = f['neurons/neuron_'+str(cell_idx)+'/tailSlope'][0]
            self.trem       = f['neurons/neuron_'+str(cell_idx)+'/trem'][0]

            self.time       = f['neurons/neuron_'+str(cell_idx)+'/time'][:]
            self.quality    = np.int(f['neurons/neuron_'+str(cell_idx)+'/quality'][0])
            self.cell_idx   = cell_idx

            # calculate the half-width correctly:
            sampling_rate   = 24414.06
            og_samples      = np.shape(self.meantrace)[0]

            if og_samples == 33 or og_samples == 21:
                interp_factor = 1
            elif og_samples == 59 or og_samples == 97 or og_samples == 100:
                interp_factor = 3
            elif og_samples == 91 or og_samples == 161:
                interp_factor = 5

            bottom      = np.argmin(self.meantrace)
            top         = np.argmax(self.meantrace[bottom:]) + bottom
            np_samples  = top - bottom

            seconds_per_sample  = 1/(interp_factor*sampling_rate)
            ms_per_sample       = seconds_per_sample*1e3
            self.neg_pos_time   = np_samples * ms_per_sample
        self.datatype = datatype
        self.clust_idx = clust_idx


    def spikeTimeToClockTime(self):
        '''creates a new instance variable that consists of the spike times converted to clock times'''
        yr = int(self.clocktime_start_time[1:5])
        month = int(self.clocktime_start_time[6:8])
        day = int(self.clocktime_start_time[9:11])
        hr = int(self.clocktime_start_time[12:14])
        minu = int(self.clocktime_start_time[15:17])
        sec = int(self.clocktime_start_time[18:20])

        t = dt.datetime(yr,month,day,hr,minu,sec,00)
        self.clock_spk_times = []
        for time in self.time:
            tdelt = dt.timedelta(seconds = time)
            tnew = t+tdelt
            strTime = "{}-{}-{}___{}:{}:{}:{}".format(tnew.year, tnew.month, tnew.day, tnew.hour, tnew.minute, tnew.second, tnew.microsecond)
            self.clock_spk_times.append(strTime)



    def plotFR(self, axes = None, savehere=None,counter=None, binsz = 3600):
        # Plot the firing rate of the neuron in 1h bins, and add grey bars to indicate dark times. This is all based on the standard approach for our recordings and will have to be updated to accomodate any majorly different datasets (e.g. very short recordings or other L/D arrangements)
        edges   = np.arange(0,max(self.time)+binsz,binsz)
        bins    = np.histogram(self.time,edges)
        hzcount = bins[0]
        hzcount = hzcount/binsz
        hzcount[hzcount==0] = 'NaN'
        xbins   = bins[1]
        xbins   = xbins/binsz

        plt.ion()
        if axes:
            currentAxis = axes
        else:
            with sns.axes_style("white"):
                fig1        = plt.figure()
                currentAxis = plt.gca()
                plt.plot(xbins[:-1],hzcount)

        plt.gca().set_xlim([0,edges[-1]/binsz])
        ylim    = plt.gca().get_ylim()[1]

        # make an array of the lights-off times
        if self.datatype == 'hdf5':
            lt_off = np.arange(12*3600,max(self.time)+12*3600,24*3600)
            # cycle through the lights-off times and plot a transparent grey bar to indicate the dark hours
            for p in lt_off/binsz:
                currentAxis.add_patch(patches.Rectangle((p, 0), 12, ylim, facecolor="grey",alpha=0.1, edgecolor="none"))

        # deal with on/off times
        if np.size(self.onTime) == 1:
            # if there is only one on/off time, just plot as dashed lines
            plt.plot((self.onTime[0]/binsz,self.onTime[0]/binsz),(0,ylim),'g--')
            plt.plot((self.offTime[0]/binsz,self.offTime[0]/binsz),(0,ylim),'r--')

        elif np.size(self.onTime) > 1:
            # in this case, plot the start and end of the cell as a dashed line, and add a red shaded box around any remaining periods of "off" times
            count = 0
            for ee in self.offTime[:-1]/binsz:
                count += 1
                wdth = self.onTime[count]/binsz - ee

            plt.plot((self.onTime[0]/binsz,self.onTime[0]/binsz),(0,ylim),'g--')
            plt.plot((self.offTime[-1]/binsz,self.offTime[-1]/binsz),(0,ylim),'r--')

        elif np.size(self.onTime) == 0:
            print('I did not find any on/off times for this cell.')


        # if self.deprived:
        #     plt.text(12,0.7*ylim,'Deprived')
        # else:
        #     plt.text(12,0.7*ylim,'Control')


        plt.ion()
        sns.set()
        sns.despine()
        plt.gca().set_xlabel('Time (hours)')
        plt.gca().set_ylabel('Firing rate (Hz)')
        plt.show()

        if savehere:
            fig1.savefig(savehere+'NRN_'+str(counter)+'_FRplot.pdf')

    def isi_hist(self, start = 0, end = False, isi_thresh = 0.1, nbins = 101):
        '''Return a histogram of the interspike interval (ISI) distribution. This is a method built into the class "neuron", so input is self (). This is typically used to evaluate whether a spike train exhibits a refractory period and is thus consistent with a single unit or a multi-unit recording. '''
        # For a view of how much this cell is like a "real" neuron, calculate the ISI distribution between 0 and 100 msec. This function will plot the bar histogram of that distribution and calculate the percentage of ISIs that fall under 2 msec.
        if end == False:
            end = max(self.time)
        idx = np.where(np.logical_and(self.time>=start, self.time<=end))[0]
        ISI = np.diff(self.time[idx])
        edges = np.linspace(0,isi_thresh,nbins)
        hist_isi        = np.histogram(ISI,edges)
        contamination   = 100*(sum(hist_isi[0][0:int((0.1/isi_thresh)*(nbins-1)/50)])/sum(hist_isi[0]))
        contamination   = round(contamination,2)
        cont_text       = 'Contamination is {} percent.' .format(contamination)

        plt.ion()
        with sns.axes_style("white"):
            fig1            = plt.figure()
            ax              = fig1.add_subplot(111)
            ax.bar(edges[1:]*1000-0.5, hist_isi[0],color='#6a79f7')
            ax.set_ylim(bottom = 0)
            ax.set_xlim(left = 0)
            ax.set_xlabel('ISI (ms)')
            ax.set_ylabel('Number of intervals')
            ax.text(30,0.7*ax.get_ylim()[1],cont_text)
        sns.despine()
        return ISI

    def crosscorr(self,friend=None,dt=1e-3,tspan=1.0,nsegs=None):
        # if "friend" argument is not give, compute autocorrelation
        if friend is None:
            print('Computing autocorrelation for cell {:d}.'.format(self.cell_idx))
            print('  Parameters:\n\tTime span: {} ms\n\tBin step: {:.2f} ms'.format(int(tspan*1000),dt*1000))
            # select spike timestamps within on/off times for this cell
            stimes1     = self.time[ (self.time>self.onTime[0]) & (self.time<self.offTime[-1]) ]
            # remove spikes at the edges (where edges are tspan/2)
            t_start     = time.time()
            subset      = [ (stimes1 > stimes1[0]+tspan/2) & (stimes1 < stimes1[-1] - tspan/2) ]
            # line above returns an array of booleans. Convert to indices
            subset      = np.where(subset)[1]
            # Take a subset of indices. We want "nsegs" elements, evenly spaced. "segindx" therefore contains nsegs spike indices, evenly spaced between the first and last one in subset
            if nsegs is None:
                nsegs   = np.int(np.ceil(np.max(self.time/120)))
            print('\tUsing {:d} segments.'.format(nsegs))
            segindx     = np.ceil( np.linspace( subset[0], subset[-1], nsegs) )
            # The spikes pointed at by the indices in "segindx" are our reference spikes for autocorrelation calculation

            # initialize timebins
            timebins    = np.arange(0,tspan+dt,dt)
            # number of bins
            nbins       = timebins.shape[0] - 1

            # initialize autocorrelation array
            ACorrs      = np.zeros((nsegs,2*nbins-1),float)

            # ref is the index of the reference spike
            for i,ref in enumerate(segindx):
                ref = int(ref)
                # "t" is the timestamp of reference spike
                t = stimes1[ref]
                # find indices of spikes between t and t+tspan
                spikeindx = np.where((stimes1>t) & (stimes1 <= t+tspan))[0]
                # get timestamps of those and subtract "t" to get "tau" for each spike
                # "tau" is the lag between each spike
                # divide by "dt" step to get indices of bins in which those spikes fall
                spikebins = np.ceil((stimes1[spikeindx] - t)/dt)
                if spikebins.any():
                    # if any spikes were found using this method, create a binary array of spike presence in each bin
                    bintrain    = np.zeros(nbins,int)
                    bintrain[spikebins.astype(int)-1] = 1
                    # the auto-correlation is the convolution of that binary sequence with itself
                    # mode="full" ensures np.correlate uses convolution to compute correlation
                    ACorrs[i,:] = np.correlate(bintrain,bintrain,mode="full")

            # sum across all segments to get auto-correlation across dataset
            Y = np.sum(ACorrs,0)
            # remove artifactual central peak
            Y[nbins-1] = 0

            # measure time elapsed for benchmarking
            elapsed = time.time() - t_start
            print('Elapsed time: {:.2f} seconds'.format(elapsed))

            # plotting ----------------------
            # -------------------------------
            plt.ion()
            fig = plt.figure(facecolor='white')
            fig.suptitle('Auto-correlation, cell {:d}'.format(self.cell_idx))
            ax2 = fig.add_subplot(211, frame_on=False)

            ax2.bar( 1000*np.arange(-tspan+dt,tspan,dt), Y, width =  0.5, color = 'k' )
            xlim2       = 100 # in milliseconds
            tickstep2   = 20 # in milliseconds
            ax2.set_xlim(0,xlim2)
            ax2_ticks   = [i for i in range(0,xlim2+1,tickstep2)]
            ax2_labels  = [str(i) for i in ax2_ticks]
            ax2.set_xticks(ax2_ticks)
            ax2.set_xticklabels(ax2_labels)
            ax2.set_ylabel('Counts')

            ax3         = fig.add_subplot(212, frame_on=False)
            ax3.bar( 1000*np.arange(-tspan+dt,tspan,dt), Y, width =  0.5, color = 'k' )
            xlim3       = int(tspan*1000) # milliseconds - set to tspan
            tickstep3   = int(xlim3/5) # milliseconds
            ax3_ticks   = [i for i in range(-xlim3,xlim3+1,tickstep3)]
            ax3_labels  = [str(i) for i in ax3_ticks]
            ax3.set_xticks(ax3_ticks)
            ax3.set_xticklabels(ax3_labels)
            ax3.set_ylabel('Counts')
            ax3.set_xlabel('Lag (ms)')

            fig.show()
            plt.show()

            # -------------------------------

        # if friend is an instance of class "neuron", compute cross-correlation
        elif isinstance(friend, neuron):
            print('Computing cross correlation between cells {:d} and {:d}.'.format(self.cell_idx,friend.cell_idx))
            print('  Parameters:\n\tTime span: {} ms\n\tBin step: {:.2f} ms'.format(int(tspan*1000),dt*1000))
            # select spike timestamps within on/off times for self cell
            stimes1 = self.time[ (self.time>self.onTime[0]) & (self.time<self.offTime[-1]) ]
            # select spike timestamps within on/off times for self cell
            stimes2 = friend.time[ (friend.time>friend.onTime[0]) & (friend.time<friend.offTime[-1]) ]
            # start timer for benchmarking
            t_start = time.time()
            # remove spikes at the edges (where edges are tspan/2)
            subset1      = [ (stimes1 > stimes1[0]+tspan/2) & (stimes1 < stimes1[-1] - tspan/2) ]
            subset2      = [ (stimes2 > stimes2[0]+tspan/2) & (stimes2 < stimes2[-1] - tspan/2) ]
            # line above returns an array of booleans. Convert to indices
            subset1      = np.where(subset1)[1]
            subset2      = np.where(subset2)[1]
            # check to see if nsegs is user provided or default
            if nsegs is None:
                nsegs1  = np.int(np.ceil(np.max(self.time/120)))
                nsegs2  = np.int(np.ceil(np.max(friend.time/120)))
                nsegs   = max(nsegs1,nsegs2)
            print('\tUsing {:d} segments.'.format(nsegs))
            # Take a subset of indices. We want "nsegs" elements, evenly spaced. "segindx" therefore contains nsegs spike indices, evenly spaced between the first and last one in subset
            segindx1     = np.ceil(np.linspace(subset1[0], subset1[-1], nsegs))
            segindx2     = np.ceil(np.linspace(subset2[0], subset2[-1], nsegs))

            # The spikes pointed at by the indices in "segindx" are our reference spikes for autocorrelation calculation

            # initialize timebins
            timebins = np.arange(0, tspan+dt,   dt)
            # number of bins
            nbins    = timebins.shape[0] - 1

            # initialize autocorrelation array
            XCorrs = np.zeros((nsegs,2*nbins-1),float)

            # ref is the index of the reference spike
            for i, ref in enumerate(segindx1):
                ref = int(ref)
                # "t" is the timestamp of reference spike
                t = stimes1[ref]
                # find indices of spikes between t and t+tspan, for cell SELF
                spikeindx1 = np.where((stimes1>t) & (stimes1 <= t+tspan))[0]
                # same thing but for cell FRIEND
                spikeindx2 = np.where((stimes2>t) & (stimes2 <= t+tspan))[0]
                # get timestamps of those and subtract "t" to get "tau" for each spike
                # "tau" is the lag between each spike
                spikebins1 = np.ceil((stimes1[spikeindx1] - t)/dt)
                spikebins2 = np.ceil((stimes2[spikeindx2] - t)/dt)
                if spikebins1.any() & spikebins2.any():
                    # binary sequence for cell SELF
                    bintrain1  = np.zeros(nbins, int)
                    bintrain1[spikebins1.astype(int)-1] = 1
                    # binary sequence for cell FRIEND
                    bintrain2   = np.zeros(nbins, int)
                    bintrain2[spikebins2.astype(int)-1] = 1
                    XCorrs[i, :] = np.correlate(bintrain1, bintrain2, mode="full")

            Y = np.sum(XCorrs,0)

            elapsed     = time.time() - t_start
            print('Elapsed time: {:.2f} seconds'.format(elapsed))
            plt.ion()
            figx        = plt.figure(facecolor='white')
            figx.suptitle('Cross-correlation, cells {:d} and {:d}'.format(self.cell_idx,friend.cell_idx))

            ax2         = figx.add_subplot(211, frame_on=False)

            ax2.bar( 1000*np.arange(-tspan+dt,tspan,dt), Y, width =  0.5, color = 'k' )
            xlim2       = int(200) # in milliseconds
            tickstep2   = int(xlim2/5) # in milliseconds
            ax2.set_xlim(-xlim2,xlim2)
            ax2_ticks   = [i for i in range(-xlim2,xlim2+1,tickstep2)]
            ax2_labels  = [str(i) for i in ax2_ticks]
            #ax2_labels = str(ax2_labels)
            ax2.set_xticks(ax2_ticks)
            ax2.set_xticklabels(ax2_labels)
            ax2.set_ylabel('Counts')

            ax3         = figx.add_subplot(212, frame_on=False)
            ax3.bar( 1000*np.arange(-tspan+dt,tspan,dt), Y, width =  0.5, color = 'k' )
            xlim3       = int(tspan*1000) # milliseconds - set to tspan
            tickstep3   = int(xlim3/5) # milliseconds
            ax3_ticks   = [i for i in range(-xlim3,xlim3+1,tickstep3)]
            ax3_labels  = [str(i) for i in ax3_ticks]
            ax3.set_xticks(ax3_ticks)
            ax3.set_xticklabels(ax3_labels)
            ax3.set_ylabel('Counts')
            ax3.set_xlabel('Lag (ms)')

            figx.show()
            plt.show()

        elif not isinstance(friend, neuron):
            # error out
            raise TypeError('ERROR: Second argument to crosscorr should be an instance of ''neuron'' class.')


    def checkqual(self, save_update=True):

        #binsz set as elapsed time/100
        elapsed = self.time[-1] - self.time[0]
        binsz = elapsed/100

        if np.size(np.shape(self.offTime)) == 2:
            offts = np.squeeze(self.offTime)
        else:
            offts = self.offTime

        if np.size(np.shape(self.onTime)) == 2:
            onts = np.squeeze(self.onTime)
        else:
            onts = self.onTime

        try:
            ontimes         = self.time[ (self.time>onts[0]) & (self.time<offts[-1]) ]
        except:
            Tracer()()

        #determine contamination and set up ISI histogram
        isis            = np.diff(ontimes)
        edges           = np.linspace(0,1.0,1001)
        oldedges        = np.linspace(0,0.1,101)
        hist_isi        = np.histogram(isis,edges)
        oldhist         = np.histogram(isis,oldedges)
        contamination   = 100*(sum(hist_isi[0][0:2])/sum(hist_isi[0]))
        contamination   = round(contamination,2)
        oldcont         = 100*(sum(oldhist[0][0:2])/sum(oldhist[0]))
        oldcont         = round(oldcont,2)
        cont_text       = 'Contamination: {} percent.' .format(contamination)
        oldconttext     = 'Prior measure: {} percent.'.format(oldcont)



        # Check contamination by bin size and then return a descriptive statistic of the contamination in bins of bin size:
        hrcont  = np.repeat(999.0, 260)
        hct = 0
        tcount  = 0
        for ee in onts:

            iedge   = np.arange( ee, offts[tcount], binsz)

            for i in iedge[:-1]:
                tmpisi          = np.diff(self.time[ (self.time>i) & (self.time<(i+binsz)) ]) # ISIs/binsz
                hist_isi        = np.histogram(tmpisi, edges)
                hrcont[hct]     = 100*(sum(hist_isi[0][0:2])/sum(hist_isi[0]))
                hct += 1

            tcount += 1

        hrcont          = np.delete( hrcont,[np.where(hrcont==999.0)] )

        if self.time[-1] < binsz:
            print ("Not enough data for firing rate, check bin size")
            newcont_text    = 'Mean Contamination by hour: --not enough data--'
            newcont_text2   = 'Median Contamination by hour: --not enough data--'
        else:
            newcont_text    = 'Mean Contamination by bin size: {0:.{1}} percent.' .format(np.nanmean(hrcont),4)
            newcont_text2   = 'Median Contamination by bin size: {0:.{1}} percent.' .format(np.nanmedian(hrcont),4)

        plt.ion()
        plt.rcParams['font.family'] = 'serif'
        fig8            = plt.figure(8,figsize=(12, 6), dpi=100,facecolor='white')
        sns.set(style="ticks")
        sns.despine()


        numGraphs=2

        if self._has_twf:
            numGraphs = 3


        # PLOT THE ISI HISTOGRAM:
        ax1             = plt.subplot(1,numGraphs,1,frame_on=False)
        # create ISI historgram over entire period of cell on, and plot that here
        tmpisi          = np.diff(self.time[ (self.time>self.onTime) & (self.time<self.offTime) ]) # ISIs/binsz
        hist_isi        = np.histogram(tmpisi, edges)
        ax1.bar(edges[1:]*1000-0.5, hist_isi[0],color=(0.7, 0.7, 0.7))
        ax1.set_xlim(0,100)
        ax1.set_ylim(bottom = 0)
        ax1.set_xlabel('ISIs (msec)')
        ax1.set_ylabel('$Number of intervals$')
        ax1.text(5,0.90*ax1.get_ylim()[1], oldconttext, fontsize="small")
        ax1.text(5,0.85*ax1.get_ylim()[1], cont_text, fontsize="small")
        ax1.text(5,0.8*ax1.get_ylim()[1], newcont_text, fontsize="small")
        ax1.text(5,0.75*ax1.get_ylim()[1], newcont_text2, fontsize="small")


        # PLOT THE MEAN TRACE (template waveform):
        if self._has_twf:
            ax1.text(0,.96, "Cell type: " + self.cell_type, transform=ax1.transAxes, color='black', fontsize='medium')
            sns.set(style="ticks")
            sns.despine()
            ax2             = plt.subplot(1,numGraphs,2,frame_on=False)
            if self.datatype == 'hdf5':
                ax2.plot(self.meantrace)
            else:
                 ax2.plot(self.wf_temp)
            ax2.set_ylabel('$Millivolts$')
        else:
            ax1.text(0,.96,"No waveform template found", transform=ax1.transAxes, color='red')

        quality_rating = "Quality set at: %s" % self.quality if self.quality!=None else "Quality NOT RATED yet"
        cell_idx_info = "cell index: %s" % self.cell_idx
        ax1.text(0,1,quality_rating, transform=ax1.transAxes, color='green')
        ax1.text(.7, 1, cell_idx_info, transform=ax1.transAxes)
        # if self.wf:
        #     ax1.text(0,.96, "Cell type: " + self.cell_type, transform=ax1.transAxes, color='black', fontsize='medium')

        # PLOT THE FIRING RATE TRACE
        edges   = np.arange(0,max(self.time)+2*binsz,binsz)
        bins    = np.histogram(self.time,edges)
        hzcount = bins[0]
        hzcount = hzcount/binsz
        hzcount[hzcount==0] = 'NaN'
        xbins   = bins[1]
        xbins   = xbins/binsz
        ax3     = plt.subplot(1,numGraphs,numGraphs,frame_on=False)
        with sns.axes_style("white"):
            ax3.plot(xbins[:-1],hzcount)

        ax3.set_xlim([0,edges[-1]/3600])
        ylim    = ax3.get_ylim()[1]

        if self.datatype == 'hdf5':
            # make an array of the lights-off times
            lt_off  = np.arange(12*3600,max(self.time)+12*3600,24*3600)

            # cycle through the lights-off times and plot a transparent grey bar to indicate the dark hours
            for p in lt_off/3600:
                ax3.add_patch(patches.Rectangle((p, 0), 12, ylim, facecolor="grey",alpha=0.1, edgecolor="none"))


        # deal with on/off times
        if np.size(onts) == 1:
            # if there is only one on/off time, just plot as dashed lines
            ax3.plot((onts[0]/3600,onts[0]/3600),(0,ylim),'g--')
            ax3.plot((offts[0]/3600,offts[0]/3600),(0,ylim),'r--')
        else:
            # in this case, plot the start and end of the cell as a dashed line, and add a red shaded box around any remaining periods of "off" times
            count = 0
            for ee in offts[:-1]/3600:
                count += 1
                wdth = onts[count]/3600 - ee
                ax3.add_patch(patches.Rectangle((ee, 0), wdth, ylim, facecolor="red",alpha=0.5, edgecolor="none"))

            ax3.plot((onts[0]/3600,onts[0]/3600),(0,ylim),'g--')
            ax3.plot((offts[-1]/3600,offts[-1]/3600),(0,ylim),'r--')


        # if self.deprived == 1:
        #     ax3.text(12,0.7*ylim,'Deprived')
        # else:
        #     ax3.text(12,0.7*ylim,'Control')

        sns.set()
        sns.despine()
        ax3.set_xlabel('Time (hours)')
        ax3.set_ylabel('Firing rate (Hz)')

        plt.subplots_adjust(wspace=.3)


        fig8.show()

        g = None
        g = input('What is the quality? ')

        #QUALITY UPDATES
        flag = None
        while flag is None:
            time.sleep(0.1)
            if g in (['1', '2', '3']):
                #currently assuming there might not be a quality file
                if save_update:
                    #always update quality if the save_update flag is true
                    self.quality=g

                    if self._has_squal:
                        #if there is already a scrubbed quality array - update the index of the cell
                        self.scrubbed_qual_array[self.cell_idx] = self.quality
                        np.save(self._scqu_file, self.scrubbed_qual_array)
                    else:
                        #make array if there isn't one already

                        self.scrubbed_qual_array = np.full(np.shape(self.unique_clusters), np.NaN)
                        self.scrubbed_qual_array[self.cell_idx] = self.quality
                        np.save("scrubbed_quality.npy", self.scrubbed_qual_array)

                flag = 1
            else:
                g = input('What is the quality? ')


        plt.close(fig8)

        # Convert g to number and compare this to the current quality. if it's different, overwrite the HDF5 file and update the quality field there.
        if self.datatype == 'hdf5':
            if self.quality != float(g):
                f1          = h5py.File(self.HDF5_tag[0], 'r+')
                data        = f1['neurons/neuron_'+str(self.HDF5_tag[1])+'/quality']
                data[...]   = float(g)
                f1.close()
                print('Updated cell quality from {} to {}. Thanks, boss!'.format(self.quality, g))

        return(g)

    def updatecontstat(self, t0 = 36*3600, t1 = 144*3600):

        try:
            cstat = self.cont_stat[0]
        except:
            cstat = self.cont_stat

        if np.size(np.shape(self.onTime))>1:
            onTs    = np.min(np.squeeze(self.onTime))
            offTs   = np.max(np.squeeze(self.offTime))
        else:
            onTs    = np.min(self.onTime)
            offTs   = np.max(self.offTime)

        try:
            tcheck = onTs <= t0 and offTs >= t1
        except:
            Tracer()()

        try:
            if not tcheck and cstat == 1:
            # change this cell to NOT continuous
                f1          = h5py.File(self.HDF5_tag[0], 'r+')
                data        = f1['neurons/neuron_'+str(self.HDF5_tag[1])+'/cont_stat']
                data[...]   = float(0)
                f1.close()
                print('Changed this cell from continuous to not.')

            elif tcheck and cstat == 0:
                # change this cell to continuous
                f1          = h5py.File(self.HDF5_tag[0], 'r+')
                data        = f1['neurons/neuron_'+str(self.HDF5_tag[1])+'/cont_stat']
                data[...]   = float(1)
                f1.close()
                print('Changed this cell to continuous.')
            else:
                print('Continuous status stays as {}.'.format(cstat))
        except:
            print('Caught at 535')
            Tracer()()

    def setonoff(self,keepprior=1, binsz = 3600):
        '''Plot the firing rate of the neuron in 1h bins and ask the user to click on the figure to indicate the on time and then again for the off time. If KEEPPRIOR is set as 1, this will preserve any currently existent on/off times that are between the newly indicated on/off time. If KEEPPRIOR is set as 0, any currently saved on/off times will be erased and replaced with whatever the user selects when running this code.'''

        global ix, ix2, cid
        ix  = None
        ix2 = None

        def __onclick(event):
            global ix, ix2, cid

            if ix is None:
                ix, iy = event.xdata, event.ydata
                print ( ' You clicked at x = {} hours. '.format(ix) )
            elif ix >0:
                ix2, iy = event.xdata, event.ydata
                print ( ' You clicked at x = {} hours. '.format(ix2) )



            bbb = plt.gcf()
            bbb.canvas.mpl_disconnect(cid)

            return ix, ix2

        def __pltfr(self):
            edges   = np.arange(0,max(self.time)+binsz,binsz)
            bins    = np.histogram(self.time,edges)
            hzcount = bins[0]
            hzcount = hzcount/3600
            hzcount[hzcount==0] = 'NaN'
            xbins   = bins[1]
            xbins   = xbins/3600

            plt.ion()
            with sns.axes_style("white"):
                fig1    = plt.figure()
                ax1     = plt.subplot(1,1,1)
                ax1.plot(xbins[:-1],hzcount)

            ax1.set_xlim(0,edges[-1]/3600)

            ylim    = ax1.get_ylim()[1]
            ylim    = plt.gca().get_ylim()[1]

            # make an array of the lights-off times
            lt_off = np.arange(12*3600,max(self.time)+12*3600,24*3600)
            #currentAxis = plt.gca()

            # cycle through the lights-off times and plot a transparent grey bar to indicate the dark hours
            for p in lt_off/3600:
                ax1.add_patch(patches.Rectangle((p, 0), 12, ylim, facecolor="grey",alpha=0.1, edgecolor="none"))

            # if self.deprived:
            #     plt.text(12,0.7*ylim,'Deprived')
            # else:
            #     plt.text(12,0.7*ylim,'Control')

            sns.set()
            sns.despine()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('Firing rate (Hz)')
            plt.show()

            # Set the on time for this cell:

            xlim    = plt.xlim()[1]
            temptxt = ax1.text(xlim*0.1, ylim*0.9, 'Click to set the on time:')

            global cid
            cid = fig1.canvas.mpl_connect('button_press_event', __onclick)


            return fig1


        #global ix, ix2
        fign = __pltfr(self)

        # set the on time  - - - - - - - - - - - -
        while ix is None:
            fign.canvas.flush_events()
            time.sleep(0.05)

        newax   = fign.get_axes()[0]
        ylim    = newax.get_ylim()[1]
        xlim    = newax.get_xlim()[1]

        newax.set_ylim(0, ylim )
        newax.plot((ix,ix),(0,ylim),'g--')

        ontime  = ix

        # And now deal with the off time: - - - - -
        txt = newax.texts[1]
        txt.set_text('And now select the off time:')

        fign.show()
        fign.canvas.flush_events()

        cid = fign.canvas.mpl_connect('button_press_event', __onclick)

        while ix2 is None:
            fign.canvas.flush_events()
            time.sleep(0.05)

        offtime = ix2
        txt.set_visible(False)
        newax.plot((ix2,ix2),(0,ylim),'r--')

        fign.show()
        fign.canvas.flush_events()

        # overwrite the old on/off time data with the updated information you just provided:
        f1      = h5py.File(self.HDF5_tag[0], 'r+')

        ontime  *= 3600
        offtime *= 3600

        if keepprior == 1:

            if np.shape(f1['neurons/neuron_'+str(self.HDF5_tag[1])+'/onTime']) != ():
                # Case in which there are already on/off times (there _should_ be).
                onTs    = np.squeeze(f1['neurons/neuron_'+str(self.HDF5_tag[1])+'/onTime'][:])
                offTs   = np.squeeze(f1['neurons/neuron_'+str(self.HDF5_tag[1])+'/offTime'][:])

                if np.size(onTs) > 1:
                    onTs    = onTs[1:]
                    offTs   = offTs[0:-1]
                    # get rid of any remaining on/off times prior to the new on time
                    onTs    = np.delete(onTs, onTs[ np.where( onTs<ontime )])
                    offTs   = np.delete(offTs, offTs[ np.where( offTs<ontime )])

                    # get rid of any remaining on/off times following the new off time
                    onTs    = np.delete(onTs, onTs[ np.where( onTs>offtime )])
                    offTs   = np.delete(offTs, offTs[ np.where( offTs>offtime )])
                elif np.size(onTs) == 1:
                    onTs    = []
                    offTs   = []
                elif np.size(onTs) == 0:
                    offTs   = []
                else:
                    print('line 660 neuron class')
                    Tracer()()

                ontdat  = np.append(ontime, onTs)
                ontdat  = np.squeeze([ontdat])
                offtdat = np.append(offTs, offtime)
                offtdat = np.squeeze([offtdat])

                del f1['neurons/neuron_'+str(self.HDF5_tag[1])+'/onTime']
                del f1['neurons/neuron_'+str(self.HDF5_tag[1])+'/offTime']

                f1.create_dataset('neurons/neuron_'+str(self.HDF5_tag[1])+'/onTime', data = [ontdat] )
                f1.create_dataset('neurons/neuron_'+str(self.HDF5_tag[1])+'/offTime', data = [offtdat] )

            elif np.shape(f1['neurons/neuron_'+str(self.HDF5_tag[1])+'/onTime']) == ():
                # Case in which there is no on/off time
                ontdat  = ontime
                offtdat = offtime

                del f1['neurons/neuron_'+str(self.HDF5_tag[1])+'/onTime']
                del f1['neurons/neuron_'+str(self.HDF5_tag[1])+'/offTime']

                f1.create_dataset('neurons/neuron_'+str(self.HDF5_tag[1])+'/onTime', data = [ontdat] )
                f1.create_dataset('neurons/neuron_'+str(self.HDF5_tag[1])+'/offTime', data = [offtdat] )




        elif keepprior == 0:
            # this is the case in which you want to WIPE the prior on/off times and start with a clean slate.
            ontdat  = ontime
            offtdat = offtime

            del f1['neurons/neuron_'+str(self.HDF5_tag[1])+'/onTime']
            del f1['neurons/neuron_'+str(self.HDF5_tag[1])+'/offTime']

            f1.create_dataset('neurons/neuron_'+str(self.HDF5_tag[1])+'/onTime', data = [ontdat] )
            f1.create_dataset('neurons/neuron_'+str(self.HDF5_tag[1])+'/offTime', data = [offtdat] )

        f1.close()


        f1          = h5py.File(self.HDF5_tag[0], 'r')
        try:
            newontime   = f1['neurons/neuron_'+str(self.HDF5_tag[1])+'/onTime'][:]
            newofftime  = f1['neurons/neuron_'+str(self.HDF5_tag[1])+'/offTime'][:]
            print('Updated to on: {} and off: {} '.format(newontime/3600, newofftime/3600))
            f1.close()

        except:
            print('Problem writing on/off times to disc. Line 696 in neuron_class.py')
            Tracer()()

        plt.close(fign)

        print('Saved new on/off times to disk.')
