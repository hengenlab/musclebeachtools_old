import pymysql
# connect to the clustering database
rootpwd = "%6m5kq2FymMXy5t3"

db = pymysql.connect(user   = "root",
                            passwd  = rootpwd,
                            host    = "localhost",
                            database= "clusteringdb")

#import mysql.connector as my
# db = my.connect(
#     host        = "localhost",
#     user        = "root",
#     passwd      = rootpwd,
#     database    = "clusteringdb"
# )
cursor = db.cursor()

# - - - - - - - - - - -  delete table - - - - - - - - - - - - - - - - -
def deltable(tablestring, cursor = cursor):
    '''Delete tables from the clusteringdb.
    Inputs:
    TABLESTRING: can be either 'implant_db' or 'clusters'

    Outputs:
    None. Prints tables before and after execution of function.'''

    cursor.execute("SHOW TABLES")
    tables_0 = cursor.fetchall()

    if tablestring is 'implant_db':

        cursor = db.cursor()

        [ cursor.execute(killcode) for killcode in ("DROP TABLE clusters", "DROP TABLE implant_db") ]

    elif tablestring is 'clusters':
        cursor.execute("DROP TABLE clusters")

    cursor.execute("SHOW TABLES")
    tables_1 = cursor.fetchall()

    print('Tables at start: {}\nTables end at: {}'.format(tables_0, tables_1))

def createimplanttable(db = db):
    '''Create the implant_db table. This should NOT be used except during development.
    May be called after the deltable function. '''
    # ---------------------------- create table for implant/region info ------------
    cursor = db.cursor()
    cursor.execute( "CREATE TABLE implant_db ( animal_id VARCHAR(255), experiment_id VARCHAR(255), species VARCHAR(255), sex VARCHAR(1), region VARCHAR(255), strain VARCHAR(255), genotype VARCHAR(255), daqsys VARCHAR(255), nchan TINYINT, chan_range VARCHAR(255), n_implant_sites TINYINT, implant_date VARCHAR(255), expt_start VARCHAR(255), expt_end VARCHAR(255), age_t0 TINYINT, surgeon VARCHAR(10), video_binary TINYINT, light_binary TINYINT, sound_binary TINYINT, sleep_state_binary TINYINT, implant_coordinates VARCHAR(255), electrode VARCHAR(255), headstage VARCHAR(255) ) "     )

    cursor.execute("ALTER TABLE implant_db ADD COLUMN implant_id INTEGER AUTO_INCREMENT PRIMARY KEY FIRST")
    print('Created table "implant_db" in the {} database'.format(db.db))

def createclusterstable():
    '''Create the clusters table. This should NOT be used except during development.
    May be called after the deltable function. '''

    cursor.execute( "CREATE TABLE clusters ( quality TINYINT, neg_pos_t TINYINT, half_width TINYINT, slope_falling TINYINT, mean_amplitude TINYINT, fr TINYINT, cluster_number TINYINT, duration SMALLINT, clustering_t0 VARCHAR(255), algorithm VARCHAR(255), track_key VARCHAR(255), implant_id INTEGER, block_label VARCHAR(255), folder_location VARCHAR(255) )" )

    # add a column for cluster barcodes and make it the primary key and make it first.
    cursor.execute("ALTER TABLE clusters ADD COLUMN barcode DOUBLE NOT NULL AUTO_INCREMENT PRIMARY KEY FIRST")

    # make foreign keys in the clusters table
    cursor.execute("ALTER TABLE clusters ADD FOREIGN KEY(implant_id) REFERENCES implant_db(implant_id) ")
    print('Created table "clusters" in the {} database'.format(db.db))

def __implantgui():
    ''' This function asks a user to input the implant/region info about each electrode used in chronic electrophys recordings. The function will write a csv file into the folder that contains the relevant dataset. The csv file will be uploaded into the lab's mysql database into the implant_db table in the clusters database.'''
    import os
    from traits.api import HasTraits, Str, Enum, Range, Directory
    from traitsui.api import View, Item, Handler  #Group,  #EnumEditor
    import numpy as np

    # set up dictionaries for GUI field enumeration
    specieslist = {

        'Unknown'   : ['unknown'],
        'mouse'     : ['m'],
        'rat'       : ['r'],
    }

    surgeonlist = {

        'EAB'   : ['lit'],
        'KBH'   : ['kbh'],
        'LIT'   : ['lit'],
        'SCF'   : ['scf'],
        'SJB'   : ['sjb'],
    }

    sexlist = {
        'unknown' : ['unknown'],
        'female': ['f'],
        'male'  : ['m'],

    }

    regionlist = {

        'Unknown'           : ['unknown'],
        'Basal Forebrain'   : ['basalforebrain'],
        'CA1'               : ['ca1'],
        'DG'                : ['dg'],
        'Endo ctx'          : ['endoctx'],
        'Ento ctx'          : ['entoctx'],
        'HC'                : ['hc'],
        'LGN'               : ['lgn'],
        'M1'                : ['m1'],
        'm2'                : ['m2'],
        'perirh ctx'        : ['perirh ctx'],
        'NAc'               : ['nac'],
        'RSC'               : ['rsc'],
        'SCN'               : ['scn'],
        'S1'                : ['s1'],
        'Subiculum'         : ['subiculum'],
        'V1m'               : ['v1m'],
        'V1b'               : ['v1b'],
        'V2'                : ['v2'],

    }

    # can conditionally display these based on species later on...
    strainlist = {

        'Unknown'       : ['unknown'],
        'Long Evans'    : ['le'],
        'C57 Black 6'   : ['c57b6'],

    }

    genotypelist = {

        'Unknown'   : ['unknown'],
        'WT'        : ['wt'],
        'Dnmt3a'    : ['dnmt3a'],
        'nf1'       : ['nf1'],
        'P301S/E4'  : ['p301se4'],
        'PV-cre'    : ['pvcre'],
        'Shank3'    : ['shank3'],
        'VIP-cre'   : ['vipcre'],

    }

    daqsyslist = {

        'Unknown'   : ['unknown'],
        'ecube'     : ['ecube'],
        'intan'     : ['intan']

    }

    tflist = {

        'yes' : [1],
        'no'  : [0]

    }

    electrodelist = {

        'Unknown'       : ['unknown'],
        'tetrode'       : ['tetrode'],
        'stereotrode'   : ['stereotrode'],
        'single wire'   : ['single'],
        'carbon'        : ['carbon'],
        'MIT silicon'   : ['mitsilicon'],
        'UCLA silicon'  : ['uclasilicon'],
        'neuropixel'    : ['neuropixel'],
    }

    headstagelist = {

        'Unknown' : ['unknown'],
        'intan16' : ['intan16'],
        'intan32' : ['intan32'],
        'intan64' : ['intan64'],
        'HS640'   : ['hs640'],

    }


    class implantupload(HasTraits):
        """
        """
        masterpath      = Directory(os.getcwd())
        animalid        = Str ('ex. ABC12345')
        experiment_id   = Str ('ex. 0010101')
        species         = Enum(list(specieslist.keys())[0], list(specieslist.keys()))
        sex             = Enum(list(sexlist.keys())[0],list(sexlist.keys()))
        region          = Enum(list(regionlist.keys())[0],list(regionlist.keys()))
        strain          = Enum(list(strainlist.keys())[0],list(strainlist.keys()))
        genotype        = Enum(list(genotypelist.keys())[0],list(genotypelist.keys()))
        daqsys          = Enum(list(daqsyslist.keys())[0],list(daqsyslist.keys()))
        nchan           = Range(low = 1, high = 640)
        chanrange       = Str ('ex. 65-128')
        nsites          = Range(low = 1, high = 20)
        implant_date    = Str ('MMDDYYYY')
        exptstart       = Str ('MMDDYYYY')
        exptend         = Str ('MMDDYYYY')
        aget0           = Range(low = 0, high = 730)
        surgeon         = Str ('initials, ex. ABC')
        videobin        = Enum(list(tflist.keys())[0],list(tflist.keys()))
        lightbin        = Enum(list(tflist.keys())[0],list(tflist.keys()))
        soundbin        = Enum(list(tflist.keys())[0],list(tflist.keys()))
        swbin           = Enum(list(tflist.keys())[0],list(tflist.keys()))
        implantcoord    = Str('ex. -2.3,0.4')
        electrode_type  = Enum(list(electrodelist.keys())[0],list(electrodelist.keys()))
        hstype          = Enum(list(headstagelist.keys())[0],list(headstagelist.keys()))

        view = View(
            #Item(name = 'topdir'),
            Item(name='masterpath',label='Directory'),
            #Group(Item(name='topdir',show_label=False)),
            Item(name = 'animalid'),
            Item(name = 'experiment_id'),
            Item(name = 'species'),
            Item(name = 'sex'),
            Item(name = 'region' ),
            Item(name = 'strain'),
            Item(name = 'genotype'),
            Item(name = 'daqsys'),
            Item(name = 'nchan'),
            Item(name = 'chanrange'),
            Item(name = 'nsites'),
            Item(name = 'implant_date'),
            Item(name = 'exptstart'),
            Item(name = 'exptend'),
            Item(name = 'aget0'),
            Item(name = 'surgeon'),
            Item(name = 'videobin'),
            Item(name = 'lightbin'),
            Item(name = 'soundbin'),
            Item(name = 'swbin'),
            Item(name = 'electrode_type'),
            Item(name = 'hstype'),

            title = 'Implant Information.',
            buttons = ['OK'],
            resizable = True,

        )

    # Create the demo:
    igui = implantupload()

    # Run the demo (if invoked from the command line):
    if __name__ == '__main__':
        igui.configure_traits()

    return igui
            # NEED TO GET THE INFORMATION OUT OF THE GUI OUTPUT AND PUT INTO DICTIONARY

def submit_implant():
    import os
    import os.path as op
    import shutil
    import csv
    import pandas as pd
    #from __future__ import absolute_import
    import numpy as np

    g = __implantgui()

    target_val_pair = {
            "animal_id": animal_id,
            "experiment_id" : experiment_id,
            "species" : species,
            "sex" : sex,
            "region" : region,
            "strain" : strain,
            "genotype" : genotype,
            "daqsys" : daqsys,
            "nchan" : nchan,
            "chan_range" : chanrange,
            "n_implant_sites" : n_sites,
            "implant_date" : implant_date,
            "expt_start" : expt_start,
            "expt_end" : expt_end,
            "age_t0" : age_t0,
            "surgeon" : surgeon,
            "video_binary" : video_binary,
            "light_binary" : light_binary,
            "sound_binary" : sound_binary,
            "sleep_state_binary" : sw_binary,
            "implant_coordinates" : implant_coordinates,
            "electrode" : electrode_type,
            "headstage" : headstage
            }


    targets = tuple( [*target_val_pair] )
    values  = tuple( [*target_val_pair.values()] )

    def escape_name(s):
        """Escape name to avoid SQL injection and keyword clashes.

        Doubles embedded backticks, surrounds the whole in backticks.

        Note: not security hardened, caveat emptor.

        """
        return '`{}`'.format(s.replace('`', '``'))

    names = list(target_val_pair)
    cols = ', '.join(map(escape_name, names))  # assumes the keys are *valid column names*.
    placeholders = ', '.join(['%({})s'.format(name) for name in names])

    query = 'INSERT INTO implant_db ({}) VALUES ({})'.format(cols, placeholders)
    cursor.execute(query, target_val_pair)
    uniqueid = cursor.execute('SELECT last_insert_id()')

    # add the fodler location and the implant barcode ID (unique, generated on commit)
    target_val_pair.update({'location':dir_phy, 'implant_id':uniqueid})
    print('Added implant information to the implant_db table in the clusteringdb database.')
    # write to a pandas dataframe, use this to write to a .csv file easily.
    df = pd.DataFrame.from_dict(data=target_val_pair, orient='index')
    fn = dir_phy + '/' + animal_id + '_' + region + '_' + str(n_sites) + '_sites.csv'
    (pd.DataFrame.from_dict(data=target_val_pair, orient='index').to_csv(fn, header=False))
    print('Wrote implant information to .csv file  {}'.format(fn))

    return uniqueid, chanrange, dir_phy
    # NOW YOU'LL NEED TO CALL THE CLUSTERCRAWL FUNCTION. WRITE SCRIPT TO first
    # call the submit implant, then go to the cluster crawl. user should have
    # control over this, but it's unclear when you would run one and not the other

def clustercrawl(topdir,ucode,chans):
    '''Crawl through all of the subfolders in the top directory - these should contain
        all of the clusters by time bin (24h?). Info will be scraped from the datafiles
        and added to the clusters database.

        INPUTS:

            TOPDIR:top level directory for an animal. This should contian all of
                the folders that hold the binned cluster data.
            UCODE: unique code for the implant being processed. This is returned
                by SUBMIT_IMPLANT and will serve as the foreign key in the clusters
                table (tying the cluster info to the implant info).
            CHANS: channel range corresponding to the UCODE referenced implant
                that is being processed.

        OUTPUTS:
        (none)

        under development.
        kbh 4/11/19
        '''
    import numpy as np
    import glob
    import os
    #cwd = os.getcwd() # save the current folder so that you can return the
    #user to this as soon as the code finishes.

    # figure out how many blocks (days?) of recording you're dealing with:
    blocks      = glob.glob(topdir +'/'+ '*_unique_clusters.npy')
    # extract the block's file name key
    blocks  = [blocks[i][0 : blocks[i].find('unique_clusters.npy')] for i in np.arange(0,np.shape(blocks)[0])]
    # extract the label component that is unique to this block and common to
    # all components related to this block
    #blocklabel  = blocks[0][0:blocks[0].find('unique_clusters.npy')])

    # set the directory
    folder_location = topdir

    for block in blocks:
        print(block)
        #DO THE thing

def cluststats(blocklabel):
    # extract the clustering algorithm
    clustal = np.load(blocklabel + 'algorithm.npy')
    # figure out t0 for the block in YYMMDD_HHMMSS
    unders = [i for i, letter in enumerate(blocklabel) if letter == '_']
    t0 = blocklabel[unders[1]+1:unders[3]]
    # - - - - - - - - - caluculate block duration  - - - - - - - -
    # get the sample rate.
    samplerate = np.squeeze( np.load( blocklabel + 'sampling_rate.npy') )
    # read out the number of sample points in the block from the file name
    ltext   = blocks[0].find('length_') + 7
    dur     = np.int( blocks[0][ ltext : blocks[0].find( '_', ltext ) ] )
    # convert sample points to time in seconds. might be better to store
    # sample points and sample rate so nothing gets lost in rounding.
    dur     /= samplerate
    # Load the block idx. Extract cluster labels and counts
    block_idx       = np.load(blocklabel + 'spike_clusters.npy')
    unique, counts  = np.unique(block_idx, return_counts=True)
    # get cluster quality. The default from KS2 is frequently junk
    qs = np.load( blocklabel + 'basicclusterqual.npy' )
    qs = qs[unique]

    # calculate firing rates.
    clust_frs = counts/dur
    # load amplitudes of all spikes
    amps = np.load(blocklabel + 'amplitudes.npy')
    # use list comprehension to calculate the mean amplitude by idx
    mean_amps = [np.mean ( amps[block_idx == i] ) for i in unique  ]

    # - - - get waveform information - - - - - - - - -
    wfs = np.load(blocklabel + 'template_waveform.npy')

    negpos  = np.zeros(unique.size)
    halfs   = np.zeros(unique.size)
    falling = np.zeros(unique.size)
    count   = 0
    for i in unique:
        # calculate peak min to max time
        negpos[count]   =  (wfs[i][wfs[i].argmin():-1].argmax()) * (1/samplerate)
        # calculate the spike halfwidth
        halfs[count]    = np.sum( wfs[i] < 0.5 * wfs[i].min() ) * (1/samplerate)
        # calculate falling phase slope (appears rising on the extracellular WF)
        # this is calculated as microvolts per millisecond
        falling[count]  = (wfs[i][28] - wfs[i][24])/ (4*(1/samplerate)*1000)


        count +=1

    #WRITE THE THINGS YOU NEED TO A PANDAS DATAFRAME AND RETURN THAT

    #convert chans from comma separated string to integers
    chans = [int(i) for i in chans.split(',')]
