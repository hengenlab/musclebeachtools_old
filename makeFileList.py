import numpy as np
import os
import glob
import math


def makeFileList(datafile, rawdatadir=False, multi_probe = False, start_block = 0, end_block = 1, probeNumber = False, fs = 25000):
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
    if datafile not in os.getcwd():
        try:
            os.chdir(os.path.expanduser('~'))
            os.chdir(os.path.realpath(datafile))
        except FileNotFoundError:
            print("*** Data File does not exist *** check the path")
            return


    file_list = list(np.zeros(11))
    files_present = []

    if multi_probe:
        ch = probeNumber
        f = "*chg_" + str(ch) + "*.npy"
        files_full = np.sort(glob.glob(f))
    else:
        files_full = np.sort(glob.glob('*.npy'))
    idx = files_full[0].find('chg_')
    baseName = files_full[0][:idx + 6]
    possible_files = ['amplitudes', 'waveform', 'qual', 'spline', 'scrubbed']

    dat = {}
    for f in possible_files:
        g = glob.glob('*{}*'.format(f))
        dat[f] = len(g) > 0

    spikefiles = [f for f in files_full if f in glob.glob('*spike_times*')]
    clustfiles = [f for f in files_full if f in glob.glob('*spike_clusters*')]
    peakfiles = [f for f in files_full if f in glob.glob('*peak*') or f in glob.glob('*max_channel*')]
    templatefiles = [f for f in files_full if f in glob.glob('*waveform.npy')]
    qual = [f for f in files_full if f in glob.glob('*qual*') if f not in glob.glob('*scrubbed*')]
    ampfiles = [f for f in files_full if f in glob.glob('*amplitudes*')]
    splinefiles = [f for f in files_full if f in glob.glob('*spline*')]

    dat['max_channel'] = len(peakfiles) > 0

    # LOADS DATA
    try:
        print("Loading files...")
        curr_clust = [np.load(clustfiles[i]) for i in range(start_block, end_block)][0]
        file_list[0] = curr_clust
        curr_spikes = [np.load(spikefiles[i]) for i in range(start_block, end_block)][0]
        file_list[1] = curr_spikes
        if dat['max_channel']:
            peak_ch = [np.load(peakfiles[i]) for i in range(start_block, end_block)][0]
            file_list[2] = peak_ch
            files_present.append('max/peak channels')
        if dat['spline']:
            templates = [np.load(splinefiles[i]) for i in range(start_block, end_block)][0]
            file_list[3] = templates
            files_present.append('smoothed waveform')
        elif dat['waveform']:
            if len(glob.glob('*mean_waveform.npy')) > 0:
                templates = [np.load(glob.glob('*mean_waveform.npy')[i]) for i in range(start_block, end_block)][0]
                file_list[3] = templates
                files_present.append('mean waveform')
            else:
                templates = [np.load(templatefiles[i]) for i in range(start_block, end_block)][0]
                file_list[3] = templates
                files_present.append('template waveform')
        if dat['qual']:
            qual_array = [np.load(qual[i]) for i in range(start_block, end_block)][0]
            file_list[5] = qual_array
            files_present.append('automated cluster quality')
        # if dat['scrubbed']:
        #     scrubbed_qual_array = np.load('scrubbed_quality.npy')
        #     file_list[6] = scrubbed_qual_array
        #     files_present.append('scrubbed cluster quality')
        if dat['amplitudes']:
            amps = [np.load(ampfiles[i]) for i in range(start_block, end_block)][0]
            file_list[4] = amps
            files_present.append('amplitudes')
    except:
        print('there was an error loading the files')
        return
    if rawdatadir:
        fname = '{}*SleepStates*.npy'.format(rawdatadir)
        files = glob.glob(fname)
        numHrs = np.round((curr_spikes[-1] - curr_spikes[0]) / (3600*fs))
        baseName = files[0][:files[0].find('SleepStates') + 11]
        sleepFiles = []
        for i in range(numHrs):
            file = baseName + str(i + 1) + '.npy'
            sleepFiles.append(file)
        sleep_states = np.zeros((2, 0))
        for idx, f in enumerate(sleepFiles):
            sleeps = np.load(f)
            timestamps = (np.nonzero(np.diff(sleeps))[0] + 1) * 4
            time_ind = (timestamps / 4) - 1
            states = sleeps[time_ind.astype(int)]
            timestamps = timestamps + (900 * idx)
            s = np.array(states)
            t = np.stack((timestamps, s))
            sleep_states = np.concatenate((sleep_states, t), axis = 1)
            last = idx
        behavior = sleep_states
        files_present.append('SLEEP STATES through hour {}'.format(last + 1))
    file_list[8] = files_present
    file_list[9]=dat

    clocktimeSIdx = spikefiles[0].find("int16_") + 6
    clocktimeEIdx = spikefiles[0].find('-P') - 1
    clocktime_start_time = spikefiles[0][clocktimeSIdx: clocktimeEIdx]

    unique_clusters = np.unique(curr_clust)

    block_label = spikefiles[0][:spikefiles[0].find('spike')]
    file_list[10]=block_label

    print("Data set information: \nThis clustering output contains:")
    s = '\t' + files_present[0]
    for f in files_present[1:]:
        s += '\n\t{}'.format(f)
    print(s + '\nRecording started at: {} \nNumber of clusters: {}'.format(clocktime_start_time,
                                                                           len(unique_clusters)))
    return file_list
