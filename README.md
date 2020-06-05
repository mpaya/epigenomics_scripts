# epigenomics_scripts

The present set of scripts represent a draft for a complete analytical workflow for ChIP-seq and RNA-seq data. It contains a configuration file, where most specific instructions, such as location of files and program parameters, are manually indicated, and a set of scripts divided into runner and auxiliary scripts. Each runner script contains, as necessary, instructions to install the required software, and instructions to find files and provide them to the auxiliary script. These scripts were designed to work in CESGA. As such, modules were the preferred option for software usage, but when not available, packages are installed with Anaconda3-2019.10 or as stand-alone.

The initial script creates the environment according to the instructions indicated in the configuration file. Next scripts run different processes within that location. Scripts are prepared to run consecutively working with the configuration file and the results from previous scripts. 
* Each runner job includes a call to the required software, the identification of the auxiliary script/s as well as source and desination directories, and specific parameters if required. These jobs can be selected to run either locally or in a HPC environment. Each submision is recorded for future track.
Runner scripts include functions to find source files for the experiment in progress, They also reconstruct the names for the destination directories.

* Auxiliary scripts contain complex instructions to run multiple software and processes, and are written either in bash or R. They are called either from a runner script or another auxiliary script. They contain running instructions for a slurm cluster. Specific running parameters not included in the configuration file may be specified here.
Specificaly, on read mapping, two software options were tested, bowtie2 and stampy, as such two folders are indicated on subsequent steps to continue the analysis of each result. At the step of peak calling, an optimization program is used and the user may select the optimized gap size on epic2 that better suite their data. Both intermediate files and tables/figures are created on these steps.

