# high-throughput BLAST

## Idea

Instead of performing thousands of blast search sequentially, if you have access to HPC or a computer with big RAM and CPU, you can parallelize the process and reduce the time needed. For example, it will be possible to do blast for 20,000 protein sequences in around 130 minutes. The idea is to split the input sequence file into multiple pieces and run the BLAST search simultaneously on each of the split file.

## Requirements

- [conda](https://docs.conda.io/en/latest/)
- [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
- [seqkit](https://bioinf.shenwei.me/seqkit/)
- [blast](https://anaconda.org/bioconda/blast)
- awk

In this file, I assume you are using HPC and slurm as workload manager.

### Installation

You can set it up whatever you like. I assume you have conda installed and explain everything after that.

Installing snakemake, blast and seqkit using [mamba](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html) (which is very similar to conda).
```
conda create -n high_throughput_blast -c conda-forge mamba
conda activate high_throughput_blast
mamba install -c conda-forge -c bioconda -n high_throughput_blast snakemake blast seqkit
```

## Protocol

1. I assume you have already moved your protein sequence file that you want to use as query to a server. You also have a set of sequences that you want to use as local database. Here, I use *Arabidopsis thaliana* [1](https://phytozome-next.jgi.doe.gov/info/Athaliana_Araport11) protein sequences to make the database and *Mesotaenium endlicherianum* [2](https://mesotaenium.uni-goettingen.de/index.html) protein sequences as my query. If the sequence IDs contain extra information beside IDs (separated by space or other delimeters), please remove them. Sometimes they cause error. You can do that via a simple regex.

Feel free to change sequence type, database, query, blast settings, etc.

2. Create a local blast database. You can do it via a bash file like this.
```
#!/bin/bash
#SBATCH -p medium # this is a partition I want to use. If you do not need it remove this.
#SBATCH -C scratch # this is a filesystem I want to access. 
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -o job-%J.out

# You should give the absolute path to your fasta file to make a local database.
# You should to execute this code in the same folder as fasta file if you prefer to use relative paths
makeblastdb -in Athaliana_447_Araport11.protein_primaryTranscriptOnly.fa -dbtype prot
```

3. Split query fasta files into many pieces using seqkit.

```
#!/bin/bash
#SBATCH -p medium # this is a partition I want to use. If you do not need it remove this.
#SBATCH -C scratch # this is a filesystem I want to access. 
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -o job-%J.out

# With -s you can define how many sequence are allowed per file after splitting. Here, we create one file per sequence.
seqkit split2 Me1_v2.release.gff3.pep.fasta -s 1
```
Let's assume I renamed the output folder that contains lots of inputs to `splited`.

4. Create a file named `Snakefile` and add these lines to it.

```
import os
# Please replace `PATH_TO_THE_FOLDER` with the proper path in the below line.
project_dir = "PATH_TO_THE_FOLDER"
PEPS, = glob_wildcards(os.path.join(project_dir, "splited/{pep}.fasta"))
# If you used a different fasta file for your local database, update the line below.
DB = os.path.join(project_dir, "Athaliana_447_Araport11.protein_primaryTranscriptOnly.fa")
######################################################################################################

rule all:
    input:
        expand(os.path.join(project_dir, "splited/{pep}.outfmt6",pep=PEPS))

rule run_blast:
    input:
        pep = os.path.join(project_dir, splited/{pep}.fasta")
    output:
        os.path.join(project_dir, splited/{pep}.outfmt6")
    threads: 15
    log:
        os.path.join(project_dir, splited/logs/{pep}_blast.log")
    shell:
        """
        blastp -db {DB} -query {input.pep} -outfmt 6 -out {output} -num_threads 15 -evalue 1e-7 -max_target_seqs 50
        """
```

5. Create a file named `cluster.yaml` and add these lines.

```
__default__:
  time: "04:00:00"
  partition: "medium"
  nodes: 1
  ntasks: 1
  ncpus: 1
  memory: "8G"
  node_properties: scratch
  job-name: "{rule}"
  output: "slurm-%j-%x.out"
  error: "slurm-%j-%x.err"
run_blast:
  time: "02:00:00"
  memory: "32G"
  ncpus: 17
```

6. Create a bash file with a name you like and add these lines.

```
#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH -c 2
#SBATCH --mem-per-cpu=3G
#SBATCH -o outfile-%J
#SBATCH -C scratch
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=USERNAME@DOMAIN.com # if you like to get an email when the process starts and ends.

SLURM_ARGS="-p {cluster.partition} -N {cluster.nodes} -n {cluster.ntasks} -c {cluster.ncpus} -t {cluster.time} -J {cluster.job-name} -o {cluster.output} -e {cluster.error} --mem={cluster.memory} -C {cluster.node_properties}"

# you can change the parameter after -j to increase or decrease the maximum number of jobs that are active simultanously.
snakemake -j 100 -pr --use-conda --cluster-config cluster.yaml --cluster "sbatch $SLURM_ARGS"
```

7. Move to the output folder and concatenate all outputs into one single fasta file:

```
cat *.outfmt6 >> blast_results.outfmt6
```

8. You can get the best blast hit using this.
```
awk '!x[$1]++' blast_results.outfmt6 > blast_results_best_blast_hit.outfmt6
```
## References
1. https://phytozome-next.jgi.doe.gov/info/Athaliana_Araport11
2. https://mesotaenium.uni-goettingen.de/index.html
