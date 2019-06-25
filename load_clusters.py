import musclebeachtools as mbt
import numpy as np



def load_clusters(datadir, rawdatadir=False, multi_probe=False, probeNumber=1, start_block=0, end_block=1):
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

    neurons = []
    for clust_idx in unique_clusters:
        neurons.append(mbt.neuron(datadir,rawdatadir, clust_idx=clust_idx, file_list=fileList))

    return neurons