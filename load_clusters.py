import musclebeachtools as mbt
import numpy as np
import os
import glob



def load_clusters(datadir, filt, tracked = False, file_startclust = False,rawdatadir=False, multi_probe=False, probeNumber=1, start_block=0, end_block=1, scrubbed = False):
    """

    :param datadir: directory with KS2 output (t folder)
    :param rawdatadir: directory with sleep scoring data
    :param multi_probe: is it multiprobbe?
    :param probeNumber: what is the probeNumber (1 if it is single probe)
    :param start_block: first block
    :param end_block: last block
    :return: list of cells
    """
    neurons = []
    os.chdir(datadir)
    fileList = mbt.makeFileList(datadir, file_startclust, rawdatadir, multi_probe, start_block, end_block, probeNumber)
    if tracked:
        keys = np.load(glob.glob('*keys*')[0])
        key_clusters = np.unique(np.where(keys != 0)[0])

        for i, clust_idx in enumerate(key_clusters): # this code picks the quality array to go off of for that cell, since every cell isn't on the first day you can't just go based on that
            if len(np.where(keys[clust_idx] == 0)) > 0:
                continuous = False
                first_block_present = np.where(keys[clust_idx] != 0)[0][0]
                relative_block = first_block_present - start_block
                if relative_block >= end_block:
                    break
                else:
                    qual_list = fileList[5][relative_block]
            else:
                qual_list = fileList[5][0]

            if qual_list[clust_idx] in filt:
                neurons.append(mbt.neuron(datadir, fileList,  start_block=start_block, end_block=end_block, rawdatadir=rawdatadir, clust_idx = clust_idx, probenumber=probeNumber))

        # do the other cells that aren't in the key file
        for block in range(end_block - start_block):
            clusts = np.unique(fileList[0][block])
            clusters = [c for c in clusts if c not in key_clusters]

            for i, clust_idx in enumerate(clusters):
                qual_list = fileList[5][block]
                if qual_list[clust_idx] in filt:
                    neurons.append(mbt.neuron(datadir, fileList,  start_block=start_block, end_block=end_block, rawdatadir=rawdatadir, clust_idx = clust_idx, probenumber=probeNumber))

        return neurons

    else:
        unique_clusters = np.unique(fileList[0])
        neurons = []
        if scrubbed:
            qual_list = np.load('scrubbed_quality_' + str(start_block) + '.npy')
            print(f'Here is the quality list (scrubbed): {qual_list}')
            for i, clust_idx in enumerate(unique_clusters):
                if qual_list[i] in filt:
                    neurons.append(mbt.neuron(datadir, fileList,  start_block=start_block, end_block=end_block, rawdatadir=rawdatadir, clust_idx = clust_idx, probenumber=probeNumber))
            return neurons
        else:
            qual_list = fileList[5][0]
            print(f'Here is the quality list: {qual_list}')
            for i, clust_idx in enumerate(unique_clusters):
                if qual_list[clust_idx] in filt:
                    neurons.append(mbt.neuron(datadir, fileList,  start_block=start_block, end_block=end_block, rawdatadir=rawdatadir, clust_idx = clust_idx, probenumber=probeNumber))
            return neurons
