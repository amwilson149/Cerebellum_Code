This folder contains the MATLAB 2015a code that runs a time evolution simulation from P3 connectivity conditions,
with the goal of reaching a connectivity pattern like the observed P7 connectivity.

It also contains the files necessary to run the code. The files "P3_Observed_PC_Connectivity_Synapse_Numbers.mat" and
"P7_Observed_PC_Connectivity_Synapse_Numbers" contain the connectivity matrices for the climbing fiber axons that innervate the P3
and P7 seed cells, respectively (i.e. these matrices contain the numbers of synapses formed by these climbing fibers onto all of the
Purkinje cells in the image volumes).  The script also makes use of the file "createPatches.m", a script submitted on MATLAB Central,
in order to visualize histograms produced by this code (https://www.mathworks.com/matlabcentral/answers/uploaded_files/30414/createPatches.m). All of these files should be in the same
directory as the MATLAB .m file in order for it to run as-is.

Finally, this folder contains the output when the code is run with the parameters in the script, with the quantile analysis turned on.

The comments within the code guide the user about where to change the simulation parameters (i.e. probability of synapse removal and
the exponential dependence of synapse addition on the number of synapses that already exist between a climbing fiber and Purkinje cell)
within the code. They also indicate which blocks of code can be uncommented in order to run additional analyses on the simulation.


