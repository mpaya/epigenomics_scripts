# Epigenomics scripts

Complete analytical workflow for ChIP-seq and RNA-seq data designed to work in Finis Terrae II server  at “Centro de Supercomputación de Galicia” (CESGA). The scripts contain a configuration file, where most specific instructions, such as location of files and program parameters, are manually indicated, and a set of scripts divided into runner and auxiliary scripts. Each runner script contains, as necessary, instructions to install the required software, and instructions to find files and provide them to the auxiliary script. Modules were the preferred option for software usage, but when not available, packages are installed with Anaconda3-2019.10 or as stand-alone.

The initial script creates the environment according to the instructions indicated in the configuration file. Next scripts run different processes within that location. Scripts are prepared to run consecutively working with the configuration file and the results from previous scripts.
* 	Each runner job includes a call to the required software, the identification of the auxiliary script/s as well as source and destination directories, and specific parameters if required. These jobs can be selected to run either locally or in a HPC environment. Each submission is recorded for future track. Runner scripts include functions to find source files for the experiment in progress, They also reconstruct the names for the destination directories. Logs from each running submission are stored.

* Auxiliary scripts contain complex instructions to run multiple software and processes, and are written either in bash or R. They are called either from a runner script or another auxiliary script. They contain running instructions for a slurm cluster. Specific running parameters not included in the configuration file may be specified here. 

## Usage
A set of scripts in numerical order are provided which perform a complete analysis for ChIP-Seq and RNA-Seq data. The starting point is creating the file system where the results will be stored. In order to do this, the file `config.txt` contains two fields, `basedir` and `refdir`, that are used as root directory for the oncoming project, and the location of resources for the reference genome. When this is ready, run:
```bash
bash 00_start.sh
```
Directories for the analysis are created for each experiment (by default chipseq and rnaseq) and a copy of these scripts placed at the root for modifications to remain on its own copy. A folder named `ARCHIVE` should contain raw data for the experiment. The folder `analysis` will contain the results benerated by the scripts.

The first step is quality control of raw reads. File names are selected from the metadata file, which contains columns for sample name, sample, For/Rev file names and grouping variables. Further specifications in `config.txt`. These files are searched in the `ARCHIVE` folder created in the first step. As parameters and options are provided through `config.txt`, it does not require additional options. The following command will submit a job to a slurm queue.
```bash
bash 01_fastqc_raw.sh
```
To run locally, the executing command may be specified (according to the extension, `bash`).
```bash
bash 01_fastqc_raw.sh bash
```

The next step is raw read trimming. It also find files in `ARCHIVE` but stores results in a newly created folder within `analysis`. 
```bash
bash 02_skewer.sh
```

Following steps find files in directories created along the analysis. The last step, `13_multiqc.sh`, gathers results and generate reports and files ready for publication.
