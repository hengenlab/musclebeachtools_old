import os
import os.path as op
import numpy as np
import pymysql
import shutil
import csv
import pandas as pd
import glob
from traits.api import HasTraits, Str, Enum, Range, Directory
from traitsui.api import View, Item, Handler
import pdb

def connectclusterdb (usr, pwd):
    ''' CONNECTCLUSTERDB. Connect to the clusteringdb database.
    Inputs:
        USER: username
        PWD: password

    Outputs:
        DB.CURSOR: "A database cursor is an identifier associated with a group
            of rows. It is, in a sense, a pointer to the current row in a buffer.
            You must use a cursor in the following cases: Statements that return
            more than one row of data from the database server: A SELECT statement
            requires a select cursor."
        DB: Database connection.
     '''
    # connect to the clustering database
    # pwd = "%6m5kq2FymMXy5t3"
    # usr = "root"
    db = pymysql.connect(user   = usr,
                                passwd  = pwd,
                                host    = "localhost",
                                database= "clusteringdb")

    return db.cursor(), db

# this is how you'll interact with the database
# - - - - - - - - - - -  delete table - - - - - - - - - - - - - - - - -
def deltable(tablestring, cursor, db):
    '''DELTABLE. Delete tables from the clusteringdb. If you select to delete
        the 'implant_db' table, you will automatically trigger deletion of the
        clusters database. This is not the case if you delete the clusters
        database (implant_db will remain intact). This is because the clusters
        table has a foreign key that intersects with the implant_db table.
    Inputs:
        TABLESTRING: can be either 'implant_db' or 'clusters'
        CURSOR: database interacting cursor. Run ' [cursor, db] = connectclusterdb(pwd) ' to generate.

    Outputs:
        None. Prints tables before and after execution of function.'''

    cursor.execute("SHOW TABLES")
    tables_0 = cursor.fetchall()

    if tablestring is 'implant_db':

        [ cursor.execute(killcode) for killcode in ("DROP TABLE clusters", "DROP TABLE implant_db") ]

    elif tablestring is 'clusters':

        cursor.execute("DROP TABLE clusters")

    cursor.execute("SHOW TABLES")
    tables_1 = cursor.fetchall()

    print('Tables at start: {}\nTables end at: {}'.format(tables_0, tables_1))
    print('Changes apply to the {} database'.format(db.db))

def createimplanttable(cursor,db):
    '''Create the implant_db table. This should NOT be used except during development.
    May be called after the deltable function. '''
    # ---------------------------- create table for implant/region info ------------
    #cursor = db.cursor()
    cursor.execute( "CREATE TABLE implant_db ( animal_id VARCHAR(255), experiment_id VARCHAR(255), species VARCHAR(255), sex VARCHAR(255), region VARCHAR(255), strain VARCHAR(255), genotype VARCHAR(255), daqsys VARCHAR(255), nchan TINYINT, changroup TINYINT, chan_range VARCHAR(255), n_implant_sites TINYINT, implant_date VARCHAR(255), expt_start VARCHAR(255), expt_end VARCHAR(255), age_t0 SMALLINT, surgeon VARCHAR(10), video_binary TINYINT, light_binary TINYINT, sound_binary TINYINT, sleep_state_binary TINYINT, implant_coordinates VARCHAR(255), electrode VARCHAR(255), headstage VARCHAR(255) ) "     )

    cursor.execute("ALTER TABLE implant_db ADD COLUMN implant_id INTEGER NOT NULL AUTO_INCREMENT PRIMARY KEY FIRST")
    print('Created table "implant_db" in the {} database'.format(db.db))

def createclusterstable(cursor,db):
    '''Create the clusters table. This should NOT be used except during development.
    May be called after the deltable function. '''

    cursor.execute( "CREATE TABLE clusters ( quality TINYINT, neg_pos_t SMALLINT, half_width SMALLINT, slope_falling MEDIUMINT, mean_amplitude SMALLINT, fr SMALLINT, cluster_number TINYINT, duration SMALLINT, clustering_t0 VARCHAR(255), algorithm VARCHAR(255), implant_id INTEGER, tracklinks VARCHAR(255), block_label VARCHAR(255), folder_location VARCHAR(255) )" )

    # add a column for cluster barcodes and make it the primary key and make it first.
    cursor.execute("ALTER TABLE clusters ADD COLUMN barcode DOUBLE NOT NULL AUTO_INCREMENT PRIMARY KEY FIRST")

    # make foreign keys in the clusters table
    cursor.execute("ALTER TABLE clusters ADD FOREIGN KEY(implant_id) REFERENCES implant_db(implant_id) ")
    print('Created table "clusters" in the {} database'.format(db.db))

def resetfordev():
    '''Quickly clear and recreate the tables during development. '''
    cursor,db = connectclusterdb ('root','%6m5kq2FymMXy5t3')
    cursor.execute("SHOW TABLES")
    cursor.fetchall()
    cursor.execute("DROP TABLE clusters")
    cursor.execute("DROP TABLE implant_db")

    createimplanttable(cursor,db)
    createclusterstable(cursor,db)

def __implantgui():
    ''' This function asks a user to input the implant/region info about each electrode used in chronic electrophys recordings. The function will write a csv file into the folder that contains the relevant dataset. The csv file will be uploaded into the lab's mysql database into the implant_db table in the clusters database.'''

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
        'BLA'               : ['bla'],
        'CA1'               : ['ca1'],
        'CeA'               : ['cea'],
        'DG'                : ['dg'],
        'Endo ctx'          : ['endoctx'],
        'Ento ctx'          : ['entoctx'],
        'HC'                : ['hc'],
        'LGN'               : ['lgn'],
        'M1'                : ['m1'],
        'm2'                : ['m2'],
        'OFC'               : ['OFC'],
        'perirh ctx'        : ['perirh ctx'],
        'NAc'               : ['nac'],
        'RSC'               : ['rsc'],
        'SCN'               : ['scn'],
        'S1'                : ['s1'],
        'S2'                : ['s2'],
        'Subiculum'         : ['subiculum'],
        'V1m'               : ['v1m'],
        'V1b'               : ['v1b'],
        'V2'                : ['v2'],

    }

    strainlist = {

        'Unknown'       : ['unknown'],
        'Long Evans'    : ['le'],
        'Sprague Dawley': ['sd'],
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
        'HS64'    : ['hs64'],
        'HS640'   : ['hs640'],

    }

    class implantupload(HasTraits):
        """ IMPLANTUPLOAD: Class for traitsui GUI creation and subsequent datastorage.
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
        changroup       = Range (low = 1, high = 10)
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
            Item(name='masterpath',label='Directory'),
            Item(name = 'animalid'),
            Item(name = 'experiment_id'),
            Item(name = 'species'),
            Item(name = 'sex'),
            Item(name = 'region' ),
            Item(name = 'strain'),
            Item(name = 'genotype'),
            Item(name = 'daqsys'),
            Item(name = 'nchan'),
            Item(name = 'changroup'),
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
            Item(name = 'implantcoord'),
            Item(name = 'electrode_type'),
            Item(name = 'hstype'),

            title = 'Implant Information.',
            buttons = ['OK'],
            resizable = True,
            scrollable = True,

        )

    # Create the GUI:
    igui = implantupload()

    # Run the GUI (if invoked from the command line):
    if __name__ == '__main__':
        igui.configure_traits()

    return igui


def __escape_name(s):
    """ Code copied from internet source: formats string data from
    target_val_pair for error-free uploading into mySQL database.

    'Escape name to avoid SQL injection and keyword clashes.
    Doubles embedded backticks, surrounds the whole in backticks.
    Note: not security hardened, caveat emptor.''

    """
    return '`{}`'.format(s.replace('`', '``'))

def submit_implant(g,cursor,db):
    '''SUBMIT_IMPLANT Takes the data structure output from the GUI __implantgui
        and writes the data contents into the clusteringdb in the implant_db
        table. This will automatically call the clustercrawl function to crawl
        through the directory selected by the user (via the GUI), extract the
        cluster information, and write those data into the clusters table.

        Inputs:
            G: this is the output structure from __implantgui.
            CURSOR: Database cursor.
            DB: Database connection.

        Outputs:
            uniqueid,
            g.changroup,
            g.masterpath
            '''
    # convert the binary fields to 0 and 1 integers for proper formatting
    d = np.array([g.videobin, g.lightbin, g.soundbin,g.swbin])
    d = [int(i) for i in d == 'yes']

    # create a dictionary of each of the targets (names) and the correspoding data
    target_val_pair = {
            "animal_id": g.animalid,
            "experiment_id" : g.experiment_id,
            "species" : g.species,
            "sex" : g.sex,
            "region" : g.region,
            "strain" : g.strain,
            "genotype" : g.genotype,
            "daqsys" : g.daqsys,
            "nchan" : g.nchan,
            "changroup" : g.changroup,
            "chan_range" : g.chanrange,
            "n_implant_sites" : g.nsites,
            "implant_date" : g.implant_date,
            "expt_start" : g.exptstart,
            "expt_end" : g.exptend,
            "age_t0" : g.aget0,
            "surgeon" : g.surgeon,
            "video_binary" : d[0],
            "light_binary" : d[1],
            "sound_binary" : d[2],
            "sleep_state_binary" : d[3],
            "implant_coordinates" : g.implantcoord,
            "electrode" : g.electrode_type,
            "headstage" :g.hstype
            }

    # convert dictionary target and value information into tuples
    targets = tuple( [*target_val_pair] )
    #values  = tuple( [*target_val_pair.values()] )
    # automatically rewrite the target names into the proper format for mysql
    cols = ', '.join(map(__escape_name, targets))  # assumes the keys are *valid column names*.
    placeholders = ', '.join(['%({})s'.format(name) for name in targets])
    # submit to the implants_db table.
    query = 'INSERT INTO implant_db ({}) VALUES ({})'.format(cols, placeholders)

    cursor.execute(query, target_val_pair)
    uniqueid = cursor.lastrowid

    db.commit()
    print('Added implant information to the implant_db table in the clusteringdb database.')

    # add the folder location and the implant barcode ID (unique, generated on
    # commit) to the dictionary
    target_val_pair.update({'location':g.masterpath, 'implant_id':uniqueid})

    # write to a pandas dataframe, use this to write to a .csv file easily.
    df = pd.DataFrame.from_dict(data=target_val_pair, orient='index')
    fn = g.masterpath + '/' + g.animalid + '_' + g.region + '_' + str(g.nsites) + '_sites.csv'
    (pd.DataFrame.from_dict(data=target_val_pair, orient='index').to_csv(fn, header=False))
    print('Wrote implant information to .csv file  {}'.format(fn))

    return uniqueid, g.changroup, g.masterpath
    # This information should be sent to the clustercrawl function. clustercrawl
    # will automatically calculate/detect cluster metadata by implant (channel
    # group) and block (time) and write to another table in the database.

def clustercrawl(topdir,ucode,channel_group,cursor,db):
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
        '''
    # extract the base filename of the individual time-blocks (each is associated)
    # with many files. There will be files per block and per channel group (same
    # as implant number). i.e., you could have 30 blocks, each containing n channel
    # groups. You'll need to loop across time-blocks and only consider the cluster
    # from the channel group contained within the implant submission that called
    # this function.
    blocks      = glob.glob(topdir +'/'+ '*_chg_' +
        str(channel_group) +'_unique_clusters.npy')

    blocks  = [blocks[i][0 : blocks[i].find('unique_clusters.npy')] for i in np.arange(0,np.shape(blocks)[0])]


    for block in blocks:
        print('Extracting metadata on clusters from {}'.format(block))
        # get the info that needs to go into the clustersdb table
        blockdf = cluststats(block, ucode, topdir)

        # submit cluster metadata from block to the database.
        submitclusters(blockdf,cursor,db)


def cluststats(blocklabel, uniqid, clustdir):
    '''CLUSTSTATS Extract spiking statistics and measurements from clusters. Info
    will be returned as a pandas dataframe.

    Inputs:
        BLOCKLABEL: The string that is used to identify all of the files corres-
            ponding to the time-block and channel group (implant number) being
            analyzed.
        UNIQID: The unique implant barcode generated by uploading implant info
            to the implant_db table. That upload *should* be the event that
            triggers automated extraction and uploading of cluster data. UNIQID
            is used to link clusters in the clusters table to the correct implant
            info in the implant_db table.
        CLUSTDIR: The directory in which the cluster datafiles are stored.

    Outputs:
        DF: A pandas dataframe containing all of the cluster metadata. This
            *should* typically be passed to SUBMITCLUSTERS.
        '''
    # extract the clustering algorithm
    f = open(blocklabel + 'algorithm.txt','r')
    clustal = f.readlines()[0]
    f.close()
    # figure out t0 for the block in YYMMDD_HHMMSS
    unders = [i for i, letter in enumerate(blocklabel) if letter == '_']
    t0 = blocklabel[unders[1]+1:unders[3]]
    # get the block name independent of the directory information
    blockname = blocklabel[ blocklabel.rfind('/') + 1 : -1 ]

    # - - - - - - - - - caluculate block duration  - - - - - - - -
    # get the sample rate.
    samplerate = np.squeeze( np.load( blocklabel + 'sampling_rate.npy') )
    # read out the number of sample points in the block from the file name
    ltext   = blocklabel.find('length_') + 7
    dur     = np.int( blocklabel[ ltext : blocklabel.find( '_', ltext ) ] )
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

    neg_pos_t   = np.zeros(unique.size)
    halfs       = np.zeros(unique.size)
    falling     = np.zeros(unique.size)
    count       = 0
    for i in unique:
        # calculate peak min to max time
        neg_pos_t[count]=  (wfs[i][wfs[i].argmin():-1].argmax()) * (1/samplerate)
        # calculate the spike halfwidth
        halfs[count]    = np.sum( wfs[i] < 0.5 * wfs[i].min() ) * (1/samplerate)
        # calculate falling phase slope (appears rising on the extracellular WF)
        # this is calculated as microvolts per millisecond
        falling[count]  = (wfs[i][28] - wfs[i][24])/ (4*(1/samplerate)*1000)

        count +=1

    trackvar = 'none'
    # Write data into a pandas dataframe.
    df = pd.DataFrame({ 'quality' : np.squeeze(qs), 'neg_pos_t' :
     neg_pos_t, 'half_width' : halfs,  'slope_falling' : falling , 'mean_amplitude': mean_amps,
     'fr' : clust_frs, 'cluster_number' : unique, 'duration' : np.tile(dur, unique.size),
     'clustering_t0' : np.tile(t0, unique.size), 'algorithm' : np.tile(clustal, unique.size),
     'implant_id' : np.tile(uniqid,unique.size),'tracklinks':trackvar, 'block_label' : np.tile(blockname, unique.size),
     'folder_location' : np.tile(clustdir, unique.size) })

    return df

def submitclusters(clstrdf,cursor,db):
    '''SUBMITCLUSTERS Function called by other scripts that submits cluster data
    to the clusters table.

    Inputs:
        CLSTRDF: Pandas dataframe created by cluststats.
    Outputs:
        None.
        '''
    clstrdic        = clstrdf.to_dict(orient='records') # convert rows to dictionaries
    targets         = tuple( [ *clstrdic[0] ] )
    cols            = ', '.join(map(__escape_name, targets))
    placeholders    = ', '.join(['%({})s'.format(name) for name in targets])
      # assumes the keys are *valid column names*.

    for i in np.arange(0,clstrdf.shape[0]):

        tempd = clstrdic[i]

        query = 'INSERT INTO clusters ({}) VALUES ({})'.format(cols, placeholders)
        cursor.execute(query, tempd)
        uniqueid = cursor.lastrowid #execute('SELECT last_insert_id()')
        db.commit()
        print('Submitted cluster {} information to the clusters table in the clusteringdb database as clust no. {}.'.format(i,uniqueid))

    # write pandas dataframe to a .csv file in the loc folder.
    # same as mpath above
    fn  = clstrdf['folder_location'][0] + '/' + 'metadata_' + clstrdf['block_label'][0] + '_.csv'
    clstrdf.to_csv(fn)
    print('Wrote info from {} clusters to .csv file {}'.format(np.shape(clstrdf)[0],fn))


def upload_implant(user,pwd):
    '''SUBMIT_IMPLANT Top level function that calls subscripts required for
    a user to submit and upload information about an implant, and subsequently
    trigger the automated detection and uploading of information about the
    clusters (neurons) recorded on that implant.

    Inputs:
        CURSOR: Database cursor.
        DB: Database connection.
        '''
    #connect to the clusteringdb
    # pwd = "%6m5kq2FymMXy5t3"
    # usr = "root"
    cursor, db = connectclusterdb (user, pwd)
    # call the implant info GUI and collect relevant information.
    g = __implantgui()
    # format and pass information from GUI to the implant_db table
    uniqueid, chgroup, mpath = submit_implant(g,cursor,db)
    # crawl through the clustering output and create entries in the clusters table
    # that correspond to the uniqueid of the implant submitted by the user.
    clustercrawl(mpath,uniqueid,chgroup,cursor,db)

    db.close()
    cursor.close()




    # LINK THIS TO BUILDING NEURON CLASS INSTANCES
    # WRITE TOOLS FOR DELETING ITEMS FROM THE DATABASE (SEARCH, RETURN BARCODES, DELETE THOSE)

def searchclusters():
    # WRITE A SUITE OF TOOLS FOR SEARCHING AND RETURNING ITEMS IN THE DATABASES
    # Create a connection object
    cursor,db = connectclusterdb ('root','%6m5kq2FymMXy5t3')#cursorclass = pymysql.cursors.DictCursor)

    query = (
            " SELECT clusters.barcode FROM clusters JOIN implant_db ON "
            "clusters.implant_id = implant_db.implant_id WHERE "
            "clusters.mean_amplitude > 10 AND clusters.slope_falling > 0 "
            " AND implant_db.species = 'rat' AND implant_db.region IN ('CA1','HC') "
            )

    cursor.execute(query)
    ## fetching all records from the 'cursor' object
    records = cursor.fetchall()
    for record in records:
        print(record)

    # unpack records into a list then a string.
    res = [rec for rec, in records]
    res2 = ", ".join(str(x) for x in res)


    query = "  SELECT implant_db.animal_id FROM clusters JOIN implant_db ON clusters.implant_id = implant_db.implant_id WHERE clusters.barcode IN ({}) ".format(res2)

        query = 'INSERT INTO implant_db ({}) VALUES ({})'.format(cols, placeholders)
# # "where" search:
    # query = "SELECT barcode FROM clusters WHERE mean_amplitude > 10 AND clusters.implant_id = implant_db.implant_id )  "
    # records = cursor.fetchall()
    #
    #
    #
    #
    # (SELECT salesman_id
    #  FROM salesman
    #  WHERE name='Paul Adam');
    #
    #
    #
    # # SHOW  RECORDS IN TABLE - - - - - - - - - - -
    # retreive one column only
    query = "SELECT quality FROM clusters"
    # retreive all columns
    query = "SELECT * FROM clusters"
    #retrieve some columns
    query = "SELECT quality, mean_amplitude FROM clusters"
    # show the column np_samples
    query = "DESCRIBE clusters"
    ## getting records from the table
    cursor.execute(query)
    ## fetching all records from the 'cursor' object
    records = cursor.fetchall()
    ## Showing the data
    for record in records:
        print(record)
    # - - - - - - - - - - - - - - - - - - - - - - - - -