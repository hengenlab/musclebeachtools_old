import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import h5py
import time
import seaborn as sns
import glob
import math
import datetime as dt

class neuron(object):
    """
           :param datafile: directory where the clustering output is
           :param rawdatadir: directory where the raw data and the sleep states are stored
           :param datatype: 'npy' for washu 'hp5' for brandeis
           :param cell_idx: cell number to look at in relation to the total number of clusters found (different from the cluster index)
           :param start_block: what cluster block to begin loading
           :param end_block: what cluster block to end at
           :param clust_idx: this is the cluster number found in the unique clusters array,
                            if you're building a neuron from database information then you will know this number and not the cell_idx.
                            This is a better identifier than cell_idx as it correspondes with Phy, use this when possible.
           :param multi_probe: True if the recording was multiple probe
           :param probenumber: what probe to look at for multi_probe recordings
           :param fs: sampling rate
           :param file_list: a list returned by the makeFileList() function. If loading more than one cell use this function to avoid loading the same file multiple times, makes the code much faster.
    """
    def __init__(self, datadir, rawdatadir=False, datatype= 'npy', cell_idx = 0, start_block=0, clust_idx = False, end_block=1, multi_probe=False, probenumber=1, fs=25000, file_list=[]):

        def pull_files(fileName, also = None, butNot = None):
            '''
            :param fileName: wildcard name of file to pull from directory
            :param also: if there are more names to pull and save in the same list
            :param butNot: if there is a type to discard from those names. for example 'scrubbed' could be discarded when pulling 'quality' to only pull the phy qual output
            :return: list of file names
            '''
            if also is not None:
                return [f for f in files_full if f in glob.glob(fileName) or f in glob.glob(also)]
            if butNot is not None:
                return [f for f in files_full if f in glob.glob(fileName) and f not in glob.glob(butNot)]
            return [f for f in files_full if f in glob.glob(fileName)]
        def load_files(files):
            '''
            :param files: list of file names
            :return: the loaded file(s)
            '''
            return [np.load(files[i]) for i in range(start_block, end_block)]

        if datatype == 'npy': # washU data



            print('You are using WashU data')
            if not clust_idx:
                print('working on cell index '+str(cell_idx))
            else:
                print('working on cluster number '+str(clust_idx))


            # going to folder that contains processed data if it exists
            if datadir not in os.getcwd():
                try:
                    os.chdir(os.path.expanduser('~'))
                    os.chdir(os.path.realpath(datadir))
                except FileNotFoundError:
                    print("*** Data File does not exist *** check the path")
                    return
            files_present = []
            if len(file_list)>0: #if there is a file list load the files from that
                dat = file_list[9]
                files_present = file_list[8]
                print('using file list')
                curr_clust = file_list[0]
                curr_spikes = file_list[1]
                if dat['max_channel']:
                    peak_ch = file_list[2]
                templates = file_list[3]
                if dat['qual']:
                    qual_array = file_list[5]
                if dat['amplitudes']:
                    amps = file_list[4]
                block_labels = file_list[10]

            else: # load the files manually
                # SORTS DATA FILES
                if multi_probe:
                    ch = probenumber
                    f   = "*chg_"+str(ch)+"*.npy"
                    files_full = np.sort(glob.glob(f))
                else:
                    files_full = np.sort(glob.glob('*.npy'))
                possible_files = ['amplitudes', 'waveform','qual', 'spline']
                dat = {}
                for f in possible_files:
                    g = glob.glob('*{}*'.format(f))
                    dat[f] = len(g) > 0
                spikefiles = pull_files('*spike_times*')  # for file name purposes
                clustfiles = pull_files('*spike_clusters*')
                peakfiles = pull_files('*peak*', also = '*max_channel*')
                templatefiles = pull_files('*waveform.npy')
                qual = pull_files('*qual*', butNot = '*scrubbed*')
                ampfiles = pull_files('*amplitudes*')
                splinefiles = pull_files('*spline*')

                block_labels = [i[:i.find('spike')] for i in spikefiles] # block labels from all blocks


                dat['max_channel'] = len(peakfiles) > 0

                # LOADS DATA
                try:
                    print("Loading files...")

                    curr_clust = load_files(clustfiles)
                    curr_spikes = load_files(spikefiles)
                    if dat['max_channel']:
                        peak_ch = load_files(peakfiles)
                        files_present.append('max/peak channels')
                    if dat['spline']:
                        templates = load_files(splinefiles)
                        files_present.append('spine waveform *yea fix this theres definitly a better name for this')
                    elif dat['waveform']:
                        if len(glob.glob('*mean_waveform.npy')) > 0:
                            templates = load_files(glob.glob('*mean_waveform.npy'))
                            files_present.append('mean waveform')
                        else:
                            templates = load_files(templatefiles)
                            files_present.append('template waveform')
                    if dat['qual']:
                        qual_array = load_files(qual)
                        files_present.append('automated cluster quality')
                    if dat['amplitudes']:
                        amps = load_files(ampfiles)
                        files_present.append('amplitudes')
                except IndexError:
                    print("files do not exist for that day range")

            startTimeIndex = [b.find("times_") + 6 for b in block_labels]
            if startTimeIndex == -1:
                startTimeIndex = [b.find("fs") + 2 for b in block_labels]
            startTimeIndexEnd = [b.find("-timee_") for b in block_labels]
            endTimeIndexStart = [b.find("-timee_") + 7 for b in block_labels]
            endTimeIndexEnd = [b.find("_length") for b in block_labels]
            if startTimeIndexEnd == -1:
                startTimeIndexEnd = [b.find("-fs") for b in block_labels]
            clocktimeSIdx = [b.find("int16_") + 6 for b in block_labels]
            clocktimeEIdx = [b.find('-P') - 1 for b in block_labels]
            self.ecube_start_time = [int(b[startTimeIndex[i]: startTimeIndexEnd[i]]) for i, b in enumerate(block_labels)]
            self.ecube_end_time = [int(b[endTimeIndexStart[i]: endTimeIndexEnd[i]]) for i, b in enumerate(block_labels)]
            self.clocktime_start_time = [b[clocktimeSIdx[i]: clocktimeEIdx[i]] for i, b in enumerate(block_labels)]

            lengths = [self.ecube_end_time[i] - self.ecube_start_time[i] for i,b in enumerate(block_labels)]
            lengths = np.insert(lengths, 0, 0)


            if end_block-start_block > 1:
                print("this cannot be done yet")
                self.unique_clusters = [np.unique(curr_clust[i]) for i in range(start_block, end_block)]
                if not clust_idx: #establish both cluster and cell index if only one was given
                    clust_idx = self.unique_clusters[0][int(cell_idx)]
                    self.clust_idx = clust_idx
                    self.cell_idx = cell_idx
                else:
                    clust_idx=int(clust_idx)
                    cell_idx = np.where(self.unique_clusters[0] == clust_idx)[0]
                    cell_idx = cell_idx[0]
                    self.clust_idx = clust_idx
                    self.cell_idx = cell_idx
                # Tracking
                total_keys = np.load(block_labels[0] + 'keys.npy')
                KEYS = total_keys[clust_idx]

                # make instance variables from this
                peak_channel = np.zeros(end_block-start_block)
                neg_pos_time = np.zeros(end_block-start_block)
                mean_amplitude = np.zeros(end_block-start_block)
                waveform = list(np.zeros(end_block-start_block))
                spike_times = list(np.zeros(end_block - start_block))
                for d in np.arange(end_block-start_block):
                    temp_idx = int(KEYS[d]-1)

                    if temp_idx == '0':
                        print(f'this cluster was not found on block {d} of the {end_block - start_block} you are loading')
                    else:
                        if dat['max_channel']:
                            peak_channel[d] = peak_ch[d][temp_idx]
                        if dat['waveform']:
                            temp = templates[d][temp_idx]
                            bottom = np.argmin(temp)
                            top = np.argmax(temp[bottom:]) + bottom
                            np_samples = top - bottom
                            seconds_per_sample = 1 / fs
                            ms_per_sample = seconds_per_sample * 1e3
                            neg_pos_time[d] = np_samples * ms_per_sample
                            mean_amplitude[d] = np.abs(np.amin(temp))

                            # in the future it might be good to later take the mean neg_pos and do cell type from that but lets go with this for now
                            if neg_pos_time[d] >= 0.4:
                                self.cell_type = 'RSU'
                            if neg_pos_time[d] < 0.4:
                                self.cell_type = 'FS'
                            waveform[d] = temp
                        if len(glob.glob('*mean_waveform.npy')) > 0:
                            print('taking mean amplitude from mean waveform')
                            mean_waveform = load_files(sorted(glob.glob('*mean_waveform.npy')))[d]
                            mean_waveform = mean_waveform[clust_idx]
                            mean_amplitude[d] = np.abs(np.amin(mean_waveform))

                        # SPIKE STUFF
                        spk_idx = np.where(curr_clust[d] == temp_idx)[0]
                        spike_times[d] = curr_spikes[d][spk_idx] / fs
                        spike_times[d] = spike_times[d] + int(lengths[d])
                self.peak_channel = peak_channel
                self.neg_pos_time = neg_pos_time
                self.mean_amplitude = mean_amplitude
                self.waveform = waveform
                self.spike_times = spike_times


                #offset is going to be the difference in the ecube time in the block_labels
                # use keys to find cluster index for that block

                # find the time of the overlap
                # convert the overlap to seconds, cut that part out of the spike array
                # add the overlap time to the spikes in the rest of the spike times array

                #concatenate at the end


            else:
                # do peak stuff
                if dat['max_channel']:
                    if np.size(np.where(peak_ch == -1)[0]) == 0:
                        print('You are using older data, changing indexing')
                        peak_ch = peak_ch - 1
                        curr_clust = curr_clust - 1
                        curr_spikes = curr_spikes - 1

                # do spike stuff
                self.unique_clusters = np.unique(curr_clust)
                if not clust_idx: #establish both cluster and cell index if only one was given
                    clust_idx = self.unique_clusters[int(cell_idx)]
                    self.clust_idx = clust_idx
                    self.cell_idx = cell_idx
                else:
                    clust_idx=int(clust_idx)
                    cell_idx = np.where(self.unique_clusters == clust_idx)[0]
                    cell_idx = cell_idx[0]
                    self.clust_idx = clust_idx
                    self.cell_idx = cell_idx

                # load the file with the tracking keys
                # index into the start_block and find the new cluster number

                spk_idx = np.where(curr_clust == clust_idx)[0]
                spiketimes = curr_spikes[spk_idx]/fs
                self.time = np.concatenate(spiketimes)
                if dat['max_channel']:
                    self.peak_channel = peak_ch[clust_idx]
                # do template stuff
                if dat['waveform']:
                    temp = templates[clust_idx]
                    bottom = np.argmin(temp)
                    top = np.argmax(temp[bottom:]) + bottom
                    np_samples = top - bottom
                    seconds_per_sample = 1 / fs
                    ms_per_sample = seconds_per_sample * 1e3
                    self.neg_pos_time = np_samples * ms_per_sample
                    self.mean_amplitude = np.abs(np.amin(temp))
                    if self.neg_pos_time >= 0.4:
                        self.cell_type = 'RSU'
                    if self.neg_pos_time < 0.4:
                        self.cell_type = 'FS'
                    self.waveform = temp

                if len(glob.glob('*mean_waveform.npy')) > 0:
                    print('taking mean amplitude from mean waveform')
                    mean_waveform = load_files(sorted(glob.glob('*mean_waveform.npy')))[0]
                    mean_waveform = mean_waveform[clust_idx]
                    self.mean_amplitude = np.abs(np.amin(mean_waveform))

                # do quality stuff
                scrubbed_files = glob.glob(f'scrubbed_quality_{start_block}.npy')
                if len(scrubbed_files) > 0:
                    scrubbed_qual_array = np.load(scrubbed_files[0])
                    dat['scrubbed'] = True
                else:
                    scrubbed_quality = np.full(np.shape(self.unique_clusters), np.NaN)
                    np.save(f'scrubbed_quality_{start_block}.npy', scrubbed_quality)
                    self.scrubbed_qual_array = scrubbed_quality
                    dat['scrubbed'] = False
                if dat['qual']:
                    if dat['scrubbed']:
                        self.scrubbed_qual_array = scrubbed_qual_array
                        last_idx = None
                        for n in range(len(self.scrubbed_qual_array)):
                            if np.isnan(self.scrubbed_qual_array[n]):
                                last_idx = n
                                break
                        if last_idx is not None:
                            print("First unscrubbed cell index: ", last_idx)
                        if np.isnan(self.scrubbed_qual_array[cell_idx]):
                            self.quality_array = qual_array
                            self.quality = qual_array[clust_idx]
                        else:
                            self.quality = self.scrubbed_qual_array[cell_idx]
                    else:
                        self.quality_array = qual_array
                        self.quality = qual_array[clust_idx]
                else:
                    print("There is no quality rating for any of these cells. Run 'checkqual()' and change the save_update flag to True if you'd like to start a quality array")
                    self.quality = 0
                if dat['amplitudes']:
                    # this is not correct bc phy is not correct for the amplitudes
                    self.amplitudes = amps[spk_idx]

            self._dat = dat
            self.onTime = np.array([self.time[0]])
            self.offTime = np.array([self.time[-1]])
            # SLEEP-WAKE
            if rawdatadir and len(file_list)==0:
                fname = '{}*SleepStates*.npy'.format(rawdatadir)
                files = glob.glob(fname)
                numHrs = math.ceil((self.offTime - self.onTime)/3600)
                baseName = files[0][:files[0].find('SleepStates')+11]
                hour_list = [int(files[i][files[i].find('SleepStates')+11:files[i].find('SleepStates')+13]) for i in np.arange(np.size(files))]
                hour_list = np.sort(hour_list)
                sleepFiles=[]
                for i in hour_list:
                    sf = baseName+str(i)+'.npy'
                    sleepFiles.append(sf)
                sleep_states = np.zeros((2,0))
                for idx, f in enumerate(sleepFiles):
                    sleeps = np.load(f)
                    timestamps = (np.nonzero(np.diff(sleeps))[0]+1)*4
                    time_ind = (timestamps/4)-1
                    states = sleeps[time_ind.astype(int)]
                    timestamps = timestamps+(3600*idx)
                    s = np.array(states)
                    t = np.stack((timestamps,s))
                    sleep_states = np.concatenate((sleep_states,t), axis =1)
                    last = idx
                self.behavior = sleep_states
                files_present.append('You have data for the following SLEEP STATES: {}'.format(hour_list))
                if np.size(hour_list) < numHrs:
                    files_present.append('PLEASE NOTE: you do not have sleep states for the entire block!')
            elif rawdatadir:
                self.behavior = file_list[7]
            #TIME STUFF

            self.start_block = start_block
            self.end_block = end_block

            self.directory = datadir
            self.probe_number = probenumber

            if len(file_list) == 0:
                print("Data set information: \nThis clustering output contains:")
                s = '\t'+files_present[0]
                for f in files_present[1:]:
                    s += '\n\t{}'.format(f)
                print(s+'\nRecording started at: {} \nNumber of clusters: {}'.format(self.clocktime_start_time, len(self.unique_clusters)))
            if dat['qual']:
                print(f'Cell quality: {self.quality}')
            print('\n')

        #Brandeis data
        else:
            print('You are using Brandeis data')
            print('Building neuron number {}'.format(cell_idx))
            f               = h5py.File(datadir, 'r')
            try:
                self.animal     = np.array(f['neurons/neuron_'+str(cell_idx)+'/animal'][:], dtype=np.int8).tostring().decode("ascii")
            except:
                print('Stuck at line 19 in neuron_class.py')
                Tracer()()

            self.HDF5_tag   = (datadir, cell_idx)
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
        uScore = [i for i, a in enumerate(self.clocktime_start_time) if a=='_']
        dashes = [i for i, a in enumerate(self.clocktime_start_time) if a=='-']
        yr = int(self.clocktime_start_time[0:4])
        month = int(self.clocktime_start_time[dashes[0]+1:dashes[1]])
        day = int(self.clocktime_start_time[dashes[1]+1:uScore[0]])
        hr = int(self.clocktime_start_time[uScore[0]+1:dashes[2]])
        minu = int(self.clocktime_start_time[dashes[2]+1:dashes[3]])
        sec = int(self.clocktime_start_time[dashes[3]+1:])

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


    def checkqual(self, scrub_cell=False):

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

        if self._dat['waveform']:
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
        if self._dat['waveform']:
            ax1.text(0,.96, "Cell type: " + self.cell_type, transform=ax1.transAxes, color='black', fontsize='medium')
            sns.set(style="ticks")
            sns.despine()
            ax2             = plt.subplot(1,numGraphs,2,frame_on=False)
            if self.datatype == 'hdf5':
                ax2.plot(self.meantrace)
            else:
                ax2.plot(self.waveform)
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
            if g in (['1', '2', '3', '4']):
                self.quality = g
                #currently assuming there might not be a quality file
                if scrub_cell:
                    #always update quality if the scrub_cell flag is true
                    self.scrubbed_qual_array[self.cell_idx] = self.quality
                    np.save(f"scrubbed_quality_{self.start_block}.npy", self.scrubbed_qual_array)

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
