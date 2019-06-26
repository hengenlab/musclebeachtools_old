import musclebeachtools as mbt
import numpy as np



def load_clusters(datadir, rawdatadir=False, multi_probe=False, probeNumber=1, start_block=0, end_block=1, filter):
    """

    :param datadir: directory with KS2 output (t folder)
    :param rawdatadir: directory with sleep scoring data
    :param multi_probe: is it multiprobbe?
    :param probeNumber: what is the probeNumber (1 if it is single probe)
    :param start_block: first block
    :param end_block: last block
    :return: list of cells
    """
    fileList = mbt.makeFileList(datadir, rawdatadir, multi_probe, start_block, end_block, probeNumber)
    unique_clusters = np.unique(fileList[0])
    qual = fileList[5]
    neurons = []
    for clust_idx in unique_clusters:
        if qual[clust_idx] in filter:
            neurons.append(mbt.neuron(datadir,rawdatadir, clust_idx=clust_idx,start_block = start_block, end_block = end_block, file_list=fileList))
    return neurons