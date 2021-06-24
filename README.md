# Epigenomics scripts

Complete analytical workflow for ChIP-seq and RNA-seq data designed to work in Finis Terrae II server  at “Centro de Supercomputación de Galicia” (CESGA). The scripts contain a configuration file, where most specific instructions, such as location of files and program parameters, are manually indicated, and a set of scripts divided into runner and auxiliary scripts. Each runner script contains, as necessary, instructions to install the required software, and instructions to find files and provide them to the auxiliary script. Modules were the preferred option for software usage, but when not available, packages are installed with Anaconda3-2019.10 or as stand-alone.

The initial script creates the environment according to the instructions indicated in the configuration file. Next scripts run different processes within that location. Scripts are prepared to run consecutively working with the configuration file and the results from previous scripts.
* 	Each runner job includes a call to the required software, the identification of the auxiliary script/s as well as source and destination directories, and specific parameters if required. These jobs can be selected to run either locally or in a HPC environment. Each submission is recorded for future track. Runner scripts include functions to find source files for the experiment in progress, They also reconstruct the names for the destination directories.

* Auxiliary scripts contain complex instructions to run multiple software and processes, and are written either in bash or R. They are called either from a runner script or another auxiliary script. They contain running instructions for a slurm cluster. Specific running parameters not included in the configuration file may be specified here. 

## Usage
Scripts don't require extra options. A configuration file is provided with options used in multiple steps. A metadata file with sample names, file names and grouping conditions must be provided by the user, as specified in the configuration file.

By running `00_start.sh`, directories for the analysis are created for each experiment (by default chipseq and rnaseq) and a copy of these scripts placed at the root for modifications to remain on its own copy. A folder named `ARCHIVE` should contain raw data for the experiment. The folder `analysis` will contain the results of the scripts.

Each runner script contains instructions to read the configuration file and metadata as needed to retrieve inputs. Additionally, instructions are provided to install software not available in CESGA. When called, each runner script finds input data and calls specific scripts that contain the analytical code.
