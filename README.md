# ~~Please do not use this musclebeachtools, currently we are doing tests~~
###### Actually, please use this anyway, if you have errors slack me - sahara

# musclebeachtools



## Installation

### Download musclebeachtools
git clone https://github.com/hengenlab/musclebeachtools.git  
Enter your username and password

#### Windows
My Computer > Properties > Advanced System Settings > Environment Variables >  
In system variables, create a new variable  
    Variable name  : PYTHONPATH  
    Variable value : location where musclebeachtools is located  
    Click OK


#### Linux
If you are using bash shell  
In terminal open .barshrc or .bash_profile  
add this line  
`export PYTHONPATH=/location_of_musclebeachtools:$PYTHONPATH`


#### Mac
If you are using bash shell  
In terminal `cd ~/` 
then open  `.profile` using your favourite text editor  
add this line  
`export PYTHONPATH=/location_of_musclebeachtools:$PYTHONPATH`


### Dependencies
Run the line below to install dependencies   
`conda install ipython seaborn numpy scipy h5py colour pymysql traits traitsui`   
`conda install -c https://conda.binstar.org/menpo opencv`   

### Testing
To test it open powershell/terminal  
    ipython  
    `import musclebeachtools as mbt`

# Neuron Class

Init 
------
`cell = mbt.neuron( datadir = '/YOUR_PATH/', rawdatadir='/SW_PATH/', datatype='npy', cell_idx = 0, start_block=0, clust_idx = 0, end_block=1, multi_probe=False, probenumber=1, fs=25000, file_list=[])`
### <a id="init-params"></a>Parameters 
- datafile: Path to the clustering output, in string format. IF YOU USED PHY (and split clusters), MAKE SURE THIS IS THE FINAL FOLDER
- rawdatadir: Path to the sleep-wake data if you have it. Keep the default to False if you don't have this information
- datatype: type of files found in the path. Defaults to 'npy', the datatype for WashU data
- cell_idx: Cluster number relative to the total number of clusters found. Different from cluster index.
- start_block: first block from the clustering output you want to load
- clust_idx: The cluster index that is associated with that cluster. This is how everything is indexed, it corresponds to the number in Phy as well. This is the ideal way to identify a cell, use this number if you know it.
- end_block: last block (exclusive) that you want to load. For example if you clustered 2 24 hour blocks with one hour of overlap then you would say start_block is 0 and end_block is 2
- multi_probe: bool value (True or False) indicating if the reading is multi_probe or not. Defaults to False. The only change is how the data is pulled from the file, everything looks the same from the user side.
- probenumber: if your recording is multi_probe you will cluster cells relative to their probe number so make sure this is correct. Defaults to 1 if your recording is single probe.
- fs: sampling rate. Defualts to 25000
- file_list: Defaults to an empty list. Run makeFileList() and put the output here. If a list of files is passed in the code simply sets the variables based on the list rather than loading all the files every time a new cell is made. 

### Output 
When you initialize a neuron there will be a series of outputs that vary depending on how your data was clustered.
1. What data you're using (Washu or Brandeis)
2. What cluster you're working on (cluster index if given)
3. Loading files... indicates that the code is loading the files, this step will take the longest.
4. If you passed in a list of files it will print "using file list"
4. Files present in the datafile
    - since there have been many iterations of clustering this lets you know what information has been loaded into the neuron object.
5. Quality Rating 
    - this breaks down the statistics about the cluster qualities that are present in that clustering output as well as the cell you just loaded
6. Start date of the recording. (usually) in YYYY-MM-DD_mm-ss-ms but its dependent on the file name
6. Sleep-Wake information.
    - if you gave a directory for sleep_wake scoring this will tell you how many hours of information it loaded. Sleep wake states are loaded under .behavior
    - eventually i'll do a comprehensive list of the variables connected to a neuron object. 
    - for now, if you load sleep wake data, the states are stored based on transition points. A 2-D array is set up with the transition second in the first row, and the previous state in the second row. So if the first entry into behavior is [300][3] then for the first 300 seconds the animal was in state 3 (REM sleep)
    
Check Quality 
------
`cell.checkqual(save_update=True)`
### Function
This method's purpose is to give a quick overview of the neuron with the essential information. When you run the method a large figure will show up with 2 or 3 plots. It will show an ISI Histogram on the left and a Firing Rate plot on the right. If there is a template waveform file it will show that in the middle. 

At the top of the figure the cell type will be printed if there is a waveform plot. It will also print the current quality of the neuron. 

On the ISI Histogram it has the mean contamination by bin size (bin size is calculated by total elapsed time of the reading divided by 100) as well as total contamination.

The Firing Rate graph axis are determined by bin size as well. 

### Parameters 
- scrub_cell: This is a bool value (defaulted to False) that indicates whether the inputed quality rating should be saved after the method completes. If this is set to True then the quality inputed by the user will override the current quality of the neuron and save it to the scrubbed quality array and file.

### Setting The Quality 
When the figure is shown, the user will be prompted in the command line to enter a quality rating for the cell. If the save_update flag is True and there is already a scrubbed quality file then this file will be updated and saved to the directory with the other files. If there is no scrubbed quality file then one will be made and saved to the same directory. If the save_update is False then nothing is saved and the quality of the cell is not updated. 

It will only accept qualities of 1, 2, 3, 4 and will keep asking until you enter a valid option.

It will prompt you even if the scrub_cell value is False, you don't have to force the code to quit. 

Clock Time Spike Times
------
### Function
This method creates an instance variable that's equivalent to the spike times in clock time format. 

It doesn't return anything, only creates a new field of the neuron object

Common Warnings and Errors
------
#### `"*** Data File does not exist *** check the path"`
- The path you sent into the init function was incorrect. 
- This error is also common when loading multiple clusters at once because the system doesn't like 'going' to a directory it's already in. slide into my dms if you get this error and you're sure the path is correct and you're connected to the correct servers.
    - before you message me try putting a slash in front of your path, or removing the slash and trying again.

#### `"files do not exist for that day range"`
- Make sure the "block range" you entered is valid for your animal. For example if you only clustered one 24 hour chunk of time and you try to ask for the 3rd block (start_blocky=3, end_blockd=4), it will error in this way.

#### `"this cannot be done yet"`
- You tried to input a day range greater than one
- tracking will be a thing, stay tuned for all the errors this code decides to throw when that happens.

#### `Index out of bounds exception`
- This is most likely due to inputting a probe number that is too large when looking at Multi-Probe Data

# MakeFileList

# load_clusters

This function takes in parameters that will be passed into neuron_class. (See that [documentation](#init-parameters) for explanation of the variables)
1. Filter: A list of qualities that should be loaded. For example if you entered [1,3] only cells of qualities 1 and 3 will be loaded (ignores noise).
#### Return
This function returns a list of neurons that passed the specified filter

# scrubClusters






