# skCSD scripts

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

/simulation/ElcoordsDomi14.txt : electrode  coordinates for experimenal data
'/simulation/morphology/ballstick.swc' : ballstick morp: hology file
'/simulation/morphology/villa.swc' : Y-shaped neuron morphology file
'/simulation/morphology/morpho1.swc' : unknown morphology !!!!!!!!!!!!!!!!!!!!!
'/simulation/morphology/neuron_agasbogas.swc' : file not found , but in the script, should remove it!!!!!!!!!!!!!!
'/simulation/morphology/Mainen_swcLike.swc' : Mainen neuron morphology,  not in the paper
'/simulation/morphology/retina_ganglion.swc': small ganglion, not  in paper
'/simulation/morphology/Badea2011Fig2Du.CNG.swc':  retinal ganglion cell
'/simulation/morphology/DomiCell.swc': morphology of pyramidal cell from experiment
'/simulation/morphology/active.hoc': inserting active ion-channels to soma
'/alprogik/sCSDFun.R': calculating sCSD (for ballstick neurons only)

Python scripts to run LFPy running NEURON:

'simulation/LFP_calc.py'#_osc_sine.py' #random inputs all over the cell
'simulation/LFP_Y_symmetric.py'  #2 synaptic inputs 
'simulation/LFPymod_example6.py' # examply from LFPy page
'simulation/LFP_calc_sine.py' #oscillatory inputs
'simulation/LFP_calc_constInj.py' #constant current injection to soma
   


## skCSD reconstruction of membrane currents from extracellular potentials using the morhology of the cell

running the  skCSD method:
Package dependencies: 
scatterplot3d, foreach, doMC, MASS 
'RunMultiple.R': running for multiple dataset without GUI with default parameters,  working for the simulated data structures (elcoord_x_y_z, myLFP,segcoordinates.txt, connections.txt, membcurr, seglength )
'utils/kernel_withoutGUI.R': setting ranges of parameter and calculatins skCSD 
'utils/kernel_basisfunctions_Regularizal.R': 



Test cases

Required input files


## Visualization








--------------------------------------------------------
Some extra features:
kernel_withoutGUI_SelectedElect.R: select electrode  distribution in different  ways
'utils/kernel_basisfunctions_Regularizal.R' - there is an option to use  background activity, description of output files?



