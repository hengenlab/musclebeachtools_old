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
export PYTHONPATH=/location_of_musclebeachtools:$PYTHONPATH


#### Mac
If you are using bash shell  
In terminal cd ~/  
then open  .profile using your favourite text editor  
add this line  
export PYTHONPATH=/location_of_musclebeachtools:$PYTHONPATH


### Testing
To test it open powershell/terminal  
    ipython  
    import musclebeachtools as mbt

