# skCSD scripts

### This code is associated with the paper from Cserpán et al., "Revealing the distribution of transmembrane currents along the dendritic tree of a neuron from extracellular recordings". eLife, 2017. http://dx.doi.org/10.7554/eLife.29384

The scripts can be devided into 3 main groups:
1. Simulation of membrane currents and virtual recodigs of extracellular potentials
2. skCSD reconstruction of membrane currents from extracellular potentials using the morhology of the cell
3. Visualization

These various parts have different software dependencies and aims so they will be covered separately.


## Simulation of membrane currents and virtual recodigs of extracellular potentials

### Dependencies:
Software dependencies
LFPy with NEURON module installed (see: http://lfpy.github.io/information.html)
R  and RStudio

R package dependencies: rgl, gWidgets, scatterplot3d, MASS (To install there, start RStudio and type 
> install.packages(c("rgl", "gWidgets", "scatterplot3d", "MASS")))


### How to run:
Start RStudio and start the folowing scriopt in the main folder of the scripts  and type 
> source("lfpy_D14.R")

Script to create dataset for validation. Loads morphology files then runs NEURON mechanism via LFPy to get biologically realistic electrophysiological data and extracellular potential recording on the selected grid  setup.


### Test cases:

### Required input files:

### Additional used files and what they do

File dependencies:

- '/simulation/ElcoordsDomi14.txt'  electrode  coordinates for experimenal data
- '/simulation/morphology/ballstick.swc'  ballstick morp: hology file
- '/simulation/morphology/villa.swc'  Y-shaped neuron morphology file
- '/simulation/morphology/morpho1.swc'  unknown morphology !!!!!!!!!!!!!!!!!!!!!
- '/simulation/morphology/neuron_agasbogas.swc'  file not found , but in the script, should remove it!!!!!!!!!!!!!!
- '/simulation/morphology/Mainen_swcLike.swc'  Mainen neuron morphology,  not in the paper
- '/simulation/morphology/retina_ganglion.swc' small ganglion, not  in paper
- '/simulation/morphology/Badea2011Fig2Du.CNG.swc'  retinal ganglion cell
- '/simulation/morphology/DomiCell.swc' morphology of pyramidal cell from experiment
- '/simulation/morphology/active.hoc' inserting active ion-channels to soma
- '/alprogik/sCSDFun.R' calculating sCSD (for ballstick neurons only)

#### Python scripts to run LFPy running NEURON:

- 'simulation/LFP_calc.py'#_osc_sine.py' #random inputs all over the cell
- 'simulation/LFP_Y_symmetric.py'  #2 synaptic inputs 
- 'simulation/LFPymod_example6.py' # examply from LFPy page
- 'simulation/LFP_calc_sine.py' #oscillatory inputs
- 'simulation/LFP_calc_constInj.py' #constant current injection to soma
   


## skCSD reconstruction of membrane currents from extracellular potentials using the morhology of the cell

#### To run the  skCSD method:

There are the following R package dependencies: 
scatterplot3d, foreach, doMC, MASS
To install run:
> install.packages(c(" foreach", "doMC")))


'kernel.R' is the main script to run the skCD method with selected parameters, to run:
>source('kernel.R') #working drectory should be set to the main folder of the program

##### Required input files for in case of a simulational/experimental setup:

- segcoordinates.txt (coordinates of segments from NEURON simulations, each row in the membcurr file belongs to a segment in the same order)
- elcoord_x_y_z (electrode coordinates x for each electrode, then y for each then z)
- myLFP (potentials at the electrode locations)
- connections.txt contains the connection info of the segments, which is connected to which other  
- membcurr (ground truth membrane currents)
- seglength: length of each segment

For an example for these files please check the 'simulation/cell_1' pr 'simulation/gang_9x9_200' folders.

Additional scripts (not required for running):
- 'utils/kernel_basisfunctions_Regularizal.R': implementation of the skCSD method
if you want to run the skCSD method on several datasets parallelly or want to run it without GUI:
- 'RunMultiple.R': running for multiple dataset without GUI with default parameters,  working for the simulated data structures (elcoord_x_y_z, myLFP,segcoordinates.txt, connections.txt, membcurr, seglength )
- 'utils/kernel_withoutGUI.R': setting ranges of parameter and calculatins skCSD 



#### Test cases



## Visualization

The scripts for plotting the  reults of the skCSD reconstructions are to be found in the Figure folder.
For running the scripts generating the figures certain R libraries need to be installed and file pathes need to be set correctly. For each figure a folder containing the script to make the figure and in most of the cases the dataset used is provided except for Figure 7., in which case it would mean GBs of data. 
The folders contain files, which were not used for the figures as well, but selection of only the used ones would take a lot of effort. For this the reason is, that in most cases the used files containg the best skCSD reconstructions are selected automatically by the sripts based on the values of reconstruction error which is stored in an other file, hence the data folders contain the reconstruction of skCSD distribution for every possible parameter combination.
Link to data: https://www.dropbox.com/sh/aqljnsps9xvtrhx/AABJDn3LzRpwTQF4ukHh8__Ua?dl=0












--------------------------------------------------------
Some extra features:
kernel_withoutGUI_SelectedElect.R: select electrode  distribution in different  ways
'utils/kernel_basisfunctions_Regularizal.R' - there is an option to use  background activity, description of output files?


