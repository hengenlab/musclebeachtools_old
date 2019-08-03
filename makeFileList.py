import numpy as np
import os
import glob
import math
import neuraltoolkit as ntk 

def makeFileList(datadir, file_startclust = False, rawdatadir=False, multi_probe = False, start_block = 0, end_block = 1, probeNumber = False, fs = 25000, tracked = False):
    """
    sorts and loads all the files in a directory in the same way that neuon_class does
    then puts each loaded array into a list and returns that list

    this way if we're running through multiple cells at a time it'll only need to load once
    if a file doesn't exist it's loaded as 'nan'
    current format is:
    index:      data:
    0				curr_clusts
    1				curr_spikes
    2				peak_chs
    3				waveform(mean waveform, template waveform, or spline depending on what exists in the file)
    4				amplitude
    5				auto_qual_array
    6				scrubbed_qual_array
    7               behavior
    8               files_present
    9               dat dictionary
    10              block label
    """

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
        :return: the loaded file
        '''
        return [np.load(files[i]) for i in range(start_block, end_block)]

    if datadir not in os.getcwd():
        try:
            os.chdir(os.path.expanduser('~'))
            os.chdir(os.path.realpath(datadir))
        except FileNotFoundError:
            print("*** Data File does not exist *** check the path")
            return


    file_list = list(np.zeros(11))
    files_present = []


    # pulls all files in a directory according to probe number
    if multi_probe:
        ch = probeNumber
        f = "*chg_" + str(ch) + "*.npy"
        files_full = np.sort(glob.glob(f))
    else:
        files_full = np.sort(glob.glob('*.npy'))
    idx = files_full[0].find('chg_')
    baseName = files_full[0][:idx + 6]
    possible_files = ['amplitudes', 'waveform', 'qual', 'spline']

    dat = {}
    for f in possible_files: # determines if there are certain files to load
        g = glob.glob('*{}*'.format(f))
        dat[f] = len(g) > 0

    spikefiles = pull_files('*spike_times*')
    clustfiles = pull_files('*spike_clusters*')
    peakfiles = pull_files('*max_channel*')
    templatefiles = pull_files('*waveform.npy')
    qual = pull_files('*qual*', butNot = '*scrubbed*')
    ampfiles = pull_files('*amplitudes*')
    splinefiles = pull_files('*spline*')

    dat['max_channel'] = len(peakfiles) > 0


    try: # LOADS DATA
        print("Loading files...")
        curr_clust = load_files(clustfiles)
        file_list[0] = curr_clust
        curr_spikes = load_files(spikefiles)
        file_list[1] = curr_spikes
        if dat['max_channel']:
            peak_ch = load_files(peakfiles)
            file_list[2] = peak_ch
            files_present.append('max/peak channels')
        if dat['spline']:
            templates = load_files(splinefiles)
            file_list[3] = templates
            files_present.append('smoothed waveform')
        elif dat['waveform']:
            if len(glob.glob('*mean_waveform.npy')) > 0:
                templates = load_files(glob.glob('*mean_waveform.npy'))
                file_list[3] = templates
                files_present.append('mean waveform')
            else:
                templates = load_files(templatefiles)
                file_list[3] = templates
                files_present.append('template waveform')
        if dat['qual']:
            qual_array = load_files(qual)
            file_list[5] = qual_array
            files_present.append('automated cluster quality')
        if dat['amplitudes']:
            amps = load_files(ampfiles)
            file_list[4] = amps
            files_present.append('amplitudes')
    except:
        print('there was an error loading the files')
        return
    if rawdatadir:
        file_startdir = np.sort(glob.glob(rawdatadir + 'Head*.bin'))[0]        
        t1 = ntk.load_raw_binary_gain_chmap(file_startdir, 512, 'PCB_tetrode', nprobes=8, t_only=1)
        t2 = ntk.load_raw_binary_gain_chmap(rawdatadir +file_startclust, 512, 'PCB_tetrode', nprobes=8, t_only=1)
        cluster_time = (t2-t1)/1e9
        fname = '{}*SleepStates*.npy'.format(rawdatadir)
        files = glob.glob(fname)
        numHrs = int(np.round((curr_spikes[-1] - curr_spikes[0]) / (3600*fs)))
        baseName = files[0][:files[0].find('SleepStates')+11]
        hour_list = [int(files[i][files[i].find('SleepStates')+11:files[i].find('.')]) for i in np.arange(np.size(files))]
        hour_list = np.sort(hour_list)
        start_hour = math.floor(cluster_time/3600)
        these_hours = hour_list[np.where(hour_list >= start_hour)[0]]
        contcheck = np.unique(np.diff(these_hours))
        if np.size(contcheck) > 1:
            ('sleep state hours are not continuous, taking only continuous ones')
            stop = np.where(np.diff(these_hours) != 1)[0][0]
            these_hours = these_hours[:stop]
        if np.size(these_hours) > numHrs:
            these_hours = these_hours[:numHrs]
        sleepFiles=[]
        for i in these_hours:
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
        ss_datstart = (these_hours[0]-1)*3600
        ss_offset = ss_datstart-cluster_time
        sleep_states[0] = sleep_states[0]+ss_offset
        null_idx = np.where(sleep_states[0]<0)[0]
        sleep_states = np.delete(sleep_states, null_idx, 1)
        behavior = sleep_states
        file_list[7] = behavior
        files_present.append('You have data for the following SLEEP STATES: {}'.format(these_hours))
        if np.size(hour_list) < numHrs:
            files_present.append('PLEASE NOTE: you do not have sleep states for the entire block!')
    file_list[8] = files_present
    file_list[9]=dat

    clocktimeSIdx = spikefiles[0].find("int16_") + 6
    clocktimeEIdx = spikefiles[0].find('-P') - 1
    clocktime_start_time = spikefiles[0][clocktimeSIdx: clocktimeEIdx]

    if (end_block - start_block) > 1:
        keys = np.load(glob.glob('*keys*')[0])
        unique_clusters = np.unique(np.where(keys!=0)[0])
    else:
        unique_clusters = np.unique(curr_clust[0])

    block_label = [i[:i.find('spike')] for i in spikefiles]
    file_list[10]=block_label

    print("Data set information: \nThis clustering output contains:")
    s = '\t' + files_present[0]
    for f in files_present[1:]:
        s += '\n\t{}'.format(f)
    print(s + '\nRecording started at: {} \nMinimum number of clusters: {}'.format(clocktime_start_time,
                                                                           len(unique_clusters)))
    return file_list
