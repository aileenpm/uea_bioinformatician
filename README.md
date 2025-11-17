# READ ME
This repository stores my code for the interview questions for the Senior Research Associate Bioinformatician role with the Norwich Medical School. 

# Repository contents
## Task 1
**interview_test.pdf** : My code for the first pre-interview task analysing the sock problem.


## Task 2
### HPC_Script

| Script name             | Function                                                                                                                                                |
|-------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------|
| master.slurm            | Main script which runs sample pre-processing, quality checks, genome alignment and subsequent miRNA identification                                      |
| run-multiqc.slurm       | Simple script to collate all FASTQC reports for raw and trimmed samples using MultiQC.                                                                  |
| run-master.sh           | Handler script which handles batch processing from a list of samples (sample_paths.txt) and submits jobs for alignment and MultiQC to the SLURM cluster |

To submit parallelised jobs for RNA-seq alignment, I use the command `bash run-master.sh list-of-fastq-files.txt`.

#### Summary of master.slurm
![Summary of miRNA identification pipeline and tools/databases referenced.](pipeline.png)

### timeseries_clustering

| File                      | Description                                                                                    |
|---------------------------|------------------------------------------------------------------------------------------------|
| timeseries_clustering.pdf | A walkthrough of how I performed fuzzy c-means clustering on miRNA time series expression data |
| mfuzz_ro18.jpeg           | Plot of clustered miRNAs produced from above analysis                                          |

#### Clustering output
![Timeseries clustering plot](timeseries_clustering/mfuzz_ro18.jpeg)
