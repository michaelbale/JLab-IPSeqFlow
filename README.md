# JLab-Flow
=================================================
      C U T & R U N  P I P E L I N E v0.5
=================================================
Author: Michael J. Bale (mib4004@med.cornell.edu)
  Usage:
  The typical command for running the pipeline is as follows:
  nextflow run michaelbale/JLab-Flow --input 'project/*_R{1,2}.fastq.gz' --genome mouse -profile singularity
  Mandatory arguments:
    --input                       Path to input data (must be surrounded with quotes)
    --genome                      Name of Genomes reference (current supported: mouse -- mm10, human -- hg38)
 -profile                      Name of package manager system (available: docker, singularity, conda); for WCM default -- singularity is recommended, but conda works while docker does not. For minimal - use conda.

  Options:
    --executorConfig              Path to config file that contains specifics for execution. Will default to WCM SCU-specific parameters. N.B. for local running use --executorConfig conf/minimal.config
    --singleSample                Specifies that the input is a single sample and will not generate a PCA graph
    --PCATitle                    Title to be included in PCA graph; must be surrounded with quotes
    --catLanes                    Tells CnRFlow to take input files and concatenate lanes into single fastq files
    --name                        Project Name; cannot have whitespace characters
    --addBEDFilesProfile          Path to csv file with info on additional BED files for generating rStart-rEnd profile plots; csv format: rName,BEDPath
    --addBEDFilesRefPoint         Path to csv file with info on additional BED files for generating rStart +/- region profile plots; csv format: pName,BEDPath,PlusMinus,pLabel
    --workDir                     Name of folder to output concatenated fastq files to (not used unless --catLanes)
    --outdir                      Name of folder to output all results to
    --genomeAssets                Home directory of where genome-specific files are. Defaults to /athena/josefowiczlab/scratch/szj2001/JLab-Flow_Genomes

