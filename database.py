#import mysql.connector as my
import pymysql
# connect to the clustering database
rootpwd = "%6m5kq2FymMXy5t3"

db = pymysql.connect(user   = "root",
                            passwd  = rootpwd,
                            host    = "localhost",
                            database= "clusteringdb")


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

def createimplanttable():
    '''Create the implant_db table. This should NOT be used except during development.
    May be called after the deltable function. '''
    # ---------------------------- create table for implant/region info ------------
    cursor = db.cursor()
    cursor.execute( "CREATE TABLE implant_db ( animal_id VARCHAR(255), experiment_id VARCHAR(255), species VARCHAR(255), strain VARCHAR(255), sex VARCHAR(1), genotype VARCHAR(255), region VARCHAR(255), daqsys VARCHAR(255), electrode VARCHAR(255), nchan TINYINT, chan_range VARCHAR(255), n_implant_sites TINYINT, implant_date VARCHAR(255), expt_start VARCHAR(255), expt_end VARCHAR(255), age_t0 TINYINT, surgeon VARCHAR(10), video_binary TINYINT, light_binary TINYINT, sound_binary TINYINT, sleep_state_binary TINYINT, implant_coordinates VARCHAR(255), headstage VARCHAR(255) ) "     )

    cursor.execute("ALTER TABLE implant_db ADD COLUMN implant_id INTEGER AUTO_INCREMENT PRIMARY KEY FIRST")

def createclusterstable():
    '''Create the clusters table. This should NOT be used except during development.
    May be called after the deltable function. '''

    cursor.execute( "CREATE TABLE clusters ( quality TINYINT, neg_pos_t TINYINT, half_width TINYINT, slope_falling TINYINT, mean_amplitude TINYINT, fr TINYINT, cluster_number TINYINT, duration SMALLINT, clustering_date VARCHAR(255), algorithm VARCHAR(255), region VARCHAR(255), expt_id TINYINT, track_key VARCHAR(255), implant_id INTEGER, folder_location VARCHAR(255) )" )

    # add a column for cluster barcodes and make it the primary key and make it first.
    cursor.execute("ALTER TABLE clusters ADD COLUMN barcode DOUBLE NOT NULL AUTO_INCREMENT PRIMARY KEY FIRST")

    # make foreign keys in the clusters table
    cursor.execute("ALTER TABLE clusters ADD FOREIGN KEY(implant_id) REFERENCES implant_db(implant_id) ")



def submit_implant():
    ''' This function asks a user to input the implant/region info about each electrode used in chronic electrophys recordings. The function will write a csv file into the folder that contains the relevant dataset. The csv file will be uploaded into the lab's mysql database into the implant_db table in the clusters database.'''
    from tkinter import filedialog
    from tkinter import Tk
    from tkinter import messagebox
    import os
    import os.path as op
    import shutil
    import csv
    # Tkinter
    root = Tk()
    root.withdraw()

    # Change to output directory
    chances_to_give = 3
    i = 0
    while i < chances_to_give + 1:
        i = i + 1
        dir_phy = filedialog.askdirectory( title = "Select output directory")
        yes_no = messagebox.askyesno( title = "Are You Sure?",
                                     message = Output directory : '+str(dir_phy) )
        if yes_no:
            print( "Current working directory is ", dir_phy )
            #os.chdir(dir_phy)
            break
        else:
            if i == chances_to_give:
                raise Exception( 'Error: output directory not selected' )


    # these are the fields that the user needs to manually enter prior to uploading anything into the database.
    # implant (primary key) # THIS IS THE PRIMARY KEY AND WILL AUTOMATICALLY POPULATE WITH A UNIQUE VALUE
    animal_id           = input ('Animal ID (initials plus 5 digits, ex: abc12345)  ')
    experiment_id       = input ('Experiment ID code (look up in lab file, or create)  ')
    species             = input ('Species (rat or mouse, r/m)  ')
    sex                 = input ('Sex (m or f)  ')
    region              = input ('Implant region (common lab abbreviation, e.g., v1, v1m, scn, rsc, s1)  ')
    strain              = input ('Animal Strain (ex. c57b6, le (long evans), sd (sprague dawley))  ')
    genotype            = input ('Animal genotype (ex. wt, shank3+-, shank3 --)  ')
    daqsys              = input ('Data Acq. System (ecube or intan)  ')
    nchan               = int( input ('Number of channels in this implant, e.g., 64  ') )
    chanrange           = input ('What is the channel range of this implant? (e.g. 65-128 as 65,128)  ')
    n_sites             = int( input ('Total number of implant sites in this animal  ') )
    implant_date        = input ('Implant date (MMDDYYYY)  ')
    expt_start          = input ('Experiment start (MMDDYYYY)  ')
    expt_end            = input ('Experiment end (MMDDYYYY)  ')
    age_t0              = int( input ('Age at start of experiment (postnatal day e.g., 60)  ') )
    surgeon             = input ('Surgeon (3 letter initial, e.g., abc)  ')
    video_binary        = int( input ('Video? (yes video = 1, no video = 0)  ') )
    light_binary        = int( input ('L/D timestamps? (1 or 0)  ') )
    sound_binary        = int( input ('Microphone data? (1 or 0)  ') )
    sw_binary           = int( input ('Sleep scoring? (1 or 0)  ') )
    implant_coordinates = input ('Stereotaxic coordinates (lateral, then bregma), e.g., -1.0, 2.2  ')
    electrode_type      = input ('Electrode type e.g., mit silicon, carbon, tetrode  ')
    headstage           = input ('Headstage e.g., hs640, hs64, intan32, intan64  ')

    # write in a check for this. if the user wants to, s/he should be able to go back and edit before committing the info to the database.


    # upload this into the database:
    #  - - - - - - - - - - - - add data to a table - - - - - - - - - - - - -
    ## defining the Query

    # make sure that your target and value pairing match up with the fields that exist in the database.
    # the "target" in quotes is the field name in the database. the value (not in quotes) is the name
    # of the variable you're passing to the database.
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



    # look at everything in the table, like, all of the data
    query = "SELECT * FROM implant_db"

    ## getting records from the table
    cursor.execute(query)

    ## fetching all records from the 'cursor' object
    records = cursor.fetchall()

    #query = "INSERT INTO implant_db (animal_id, experiment_id, species, sex, region, strain, genotype, daqsys, nchan, chan_range, n_implant_sites, implant_date, expt_start, expt_end, age_t0, surgeon, video_binary, light_binary, sound_binary, sleep_state_binary, implant_coordinates, electrode, headstage ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s )"
    ## storing values in a variable

    #values = (animal_id, experiment_id, species, sex, region, strain, genotype, daqsys, nchan, chanrange, n_sites, implant_date, expt_start, expt_end, age_t0, surgeon, video_binary, light_binary, sound_binary, sw_binary, implant_coordinates, electrode_type, headstage)

    ## executing the query with values



csvData = [ ['animal_id', animal_id],['species', species], ['sex', sex], ['region', region], ['strain', strain], ['genotype', genotype], ['daq_system', daq_system], ['n_channels', nchan], ['n_sites', n_sites], ['implant_date', implant_date], ['expt_start', expt_start], ['expt_end', expt_end], ['age_t0', age_t0], ['surgeon', surgeon], ['video_binary', video_binary], ['light_binary', light_binary], ['sound_binary', sound_binary], ['sw_binary', sw_binary], ['implant_coordinates', implant_coordinates], ['electrode', electrode_type], ['headstage', headstage] ]



#- - - - - - - - - - - - - - - - - insert data into table - - - - - - - - - - -
sql_insert_query = " INSERT INTO `clusters`(`id`, `name`, `birth_date`, `age`) VALUES (1,'Scott','2018-04-11', 26) "

with open('person.csv', 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerows(csvData)

csvFile.close()
