import numpy as np
import os
import glob
import musclebeachtools as mb


def scrubClusters(datafile, rawdatadir, multi_probe, probeNumber,
                  start_block=0, end_block=1):

    '''
    Function to set quality of cells in a dataset using checkqual function
    start_cell is the first cell you scrub to the end

    scrubClusters(datafile, rawdatadir, multi_probe, probeNumber,
                  start_block=0, end_block=1):

    Parameters
    ----------
    datafile : directory where the clustering output is
    rawdatadir : directory where the raw data and the sleep states are stored
    multi_probe : True if the recording was multiple probe
    probeNumber : what probe to look at for multi_probe recordings
    start_block : what cluster block to begin loading
    end_block : what cluster block to end at

    Returns
    -------

    Raises
    ------

    See Also
    --------

    Notes
    -----

    Examples
    --------
    scrubClusters(datafile, rawdatadir, multi_probe, probeNumber,
                  start_block=0, end_block=1):

    '''

    os.chdir(datafile)
    file_list = mb.makeFileList(datafile, rawdatadir, multi_probe,
                                probeNumber=probeNumber)
    clusters = np.unique(file_list[0])
    num_cells = len(clusters)

    sq = glob.glob("*scrubbed_quality_{:d}.npy".format(start_block))
    start_cell = 0
    if len(sq) > 0:
        print("At least one of these cells has been scrubbed before")
        sq_array = np.load(sq[0])
        first_unscrubbed = np.where(np.isnan(sq_array))[0]
        print("The first unscrubbed cell is index {}"
              .format(first_unscrubbed[0]))
        start_cell = int(input(
            "what cell index would you like to start scrubbing at? From 0-{}"
            .format(num_cells)))

    for c in np.arange(start_cell, num_cells):
        cell = mb.neuron(datadir=datafile, rawdatadir=rawdatadir,
                         clust_idx=clusters[c], file_list=file_list,
                         start_block=start_block, end_block=end_block)
        cell.checkqual(scrub_cell=True)

    # prints the quality stats for the data set once you get to the end
    # basically how many of each quality cell there was
    quals = np.zeros(5)
    # need to restate this line in case there wasn't a scrubbed quality array
    # to begin with
    # scrubbed_qual_array = np.load(f'scrubbed_quality_{start_block}.npy')
    scrubbed_qual_array = np.load('scrubbed_quality_{:d}.npy'
                                  .format(start_block))
    for q in range(len(scrubbed_qual_array)):
        quals[int(scrubbed_qual_array[q])] = \
            quals[int(scrubbed_qual_array[q])] + 1

    print("Quality Statistics:\n1: {} \n2: {} \n3: {} \n4: {}"
          .format(quals[1], quals[2], quals[3], quals[4]))
