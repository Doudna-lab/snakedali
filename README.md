![screenshot](figures/logos/snakedali_logo.png)

![GitHub Tag](https://img.shields.io/github/v/tag/Doudna-lab/snakedali)
![GitHub top language](https://img.shields.io/github/languages/top/Doudna-lab/snakedali)
![GitHub commit activity](https://img.shields.io/github/commit-activity/m/Doudna-lab/snakedali)
![GitHub repo size](https://img.shields.io/github/repo-size/Doudna-lab/snakedali)
![GitHub Downloads (all assets, all releases)](https://img.shields.io/github/downloads/Doudna-lab/snakedali/total)


## OVERVIEW

<details>
<summary>1 What is Snakedali?</summary>

Snakedali is the Snakemake implementation of the multithreaded version of DaliLite v5 to align PDB queries to a pre-built Alphafold database. 
It introduces automated input handling and a unified report that aggregates all queries and hits in a single .xlsx.


</details>

<details>
<summary>2 Citation</summary>
<ul>

Yoon, P.H., Zhang, Z., Loi, K.J., Adler, B.A., Lahiri, A., Vohra, K., Shi, H., Rabelo, D.B., Trinidad, M., Boger, R.S. and Al-Shimary, M.J., 2024. Structure-guided discovery of ancestral CRISPR-Cas13 ribonucleases. Science, p.eadq0553.


</ul>
</details>

## GETTING STARTED

<details>
<summary>3 Dependencies</summary>
<ul>

<details>
<summary>3.1 Anaconda </summary>
<ul>

 - Install Miniconda:
 - Download the installer at: https://docs.conda.io/projects/miniconda/en/latest/
 
   ```
   bash Miniconda3-latest-<your-OS>.sh
   ```
  - Set up and update conda: 
    ```
    conda update --all
    conda config --set channel_priority strict
    ```
</ul>
</details>

<details>
<summary>3.2 Snakemake </summary>
<ul>

- Snakemake can be installed directly via Anaconda:

  ```
  conda install -n base -c conda-forge mamba
  ```
</ul>
</details>

<details>
<summary>3.3 DaliLite v5 </summary>
<ul>

- Snakedali requires a valid installation of DaliLite v5
- Snakedali defaults to the parallelized implementation of DaliLite

```
wget http://ekhidna2.biocenter.helsinki.fi/dali/DaliLite.v5.tar.gz
tar -zxvf DaliLite.v5.tar.gz
cd /home/you/DaliLite.v5/bin
make clean
make parallel
```
 - Future updates are planned to include a non-parallel version of Snakedali
 - More details on how to acquire the program can be found on the software's page:
   - http://ekhidna2.biocenter.helsinki.fi/dali/README.v5.html#install

</ul>
</details>

<details>
<summary>3.4 AWS CLI </summary>
<ul>

- Install the AWS CLI to download the snakedali database directly from the command-line interface
- A comprehensive guide on the installation is provided at the AWS page:
 - Prerequisites: https://docs.aws.amazon.com/cli/latest/userguide/getting-started-prereqs.html
 - Installation Guide: https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html

</ul>
</details>

</ul>
</details>

<details>
<summary>4 Installation</summary>
<ul>

<details>
<summary>4.1 Database Download</summary>
<ul>
  - With AWS CLI installed (see Section 3.4)
  - Download the pre-built database:

```
aws s3 cp s3://snakedali.db/pdb_files_DAT.tar.gz <your_local_path>
tar zxf <your_local_path>/pdb_files_DAT.tar.gz
```


</ul>
</details>

<details>
<summary>4.2 Standard Git</summary>
<ul>

  - Clone repository files
```
git clone https://github.com/Doudna-lab/nidali.git
```
</ul>
</details>

<details>
<summary>4.3 Git LFS</summary>
<ul>
  
  - Two singularity/apptainer containers are provided in this repository
  - Although these are support files which are <b><u>not</u></b> integrated to the pipeline, they could be useful for users who may be facing issues when trying to get DaliLite installed in unsupported machines.
  - These large files will be indexed upon cloning and will take a small amount of storage. 
  - The user can then download them with Git LFS in case they need the containerized version.
  
  - 1.1 Install Git LFS to pull apptainer containers

  -1.1.1 Linux Install
```
apt install git-lfs
git lfs install
```

  -1.1.2 macOS Install
```
brew install git-lfs
git lfs install
```

  -1.1.3 Pull apptainer containers
```
git lfs pull
```
</ul>
</details>

</ul>
</details>


<details>
<summary>5 Snakedali Pipeline Setup</summary>
<ul>

<details>
<summary>5.1 Run Configuration</summary>
<ul>

 - Each Snakedali run can be customized based on the `configuration file`: `config/dali_template.yaml`
 - This file can be replicated, and each subsequent modified yaml file is associated with one Snakedali run.
 - From the configuration file users are expected to set up:
   - In-/Output paths for the run
   - pre-built database path
   - query name(s)
   - default DaliLite v5 binary folder path

</ul>
</details>

<details>
<summary>5.2 Snakemake Profile</summary>
<ul>
  
  - Snakedali was designed to work with `(Sun Grid Engine) SGE` job scheduler
  - The Snakemake profile can be modified to accommodate other schedulers: `/profile/config.yaml`
  
  - The default profile includes:
    - cluster job submission: `qsub -l h_rt={cluster.time} -j y -pe smp 4 -cwd`
    - cluster config path: `config/cluster.yaml`
    - rerun triggers: `mtime`
    - n jobs limit: `600`
    - latency-wait: `120` 
    - reason: `True`
    - rerun-incomplete: `True`
    - show-failed-logs: `True`
    - keep-going: `True`
    - printshellcmds: `True`
    - jobname: `{rule}.{jobid}`
    - jobs: `600`       
   
  - Make sure to adjust the parameters above according to the house rules of your HPC.

</ul>
</details>

</ul>
</details>



<details>
<summary>6 Run Snakedali</summary>
<ul>

 - Once the necessary inputs have been set up in the `configuration file`, Snakedali shall be called as in:

```
snakemake --snakefile snakedali.smk --configfile config/dali_template.yaml --profile profile/
```

</ul>
</details>


<details>
<summary>7 DALI + TCOFFEE Integration</summary>
<ul>
 - The DALI + TCOFFEE workflow is broken down into two parts.
  - 1. The first script, dali_out_to_fasta.py is a python script that takes in a DALI.txt output (that is, the search results for a DALI query against a database in DALI alignment format) and converts them into individual pariwise alignment files in FASTA format.
  - 2. The second script is a wrapper that calls on the first script to take in an entire directory of DALI.txt output files to convert them into directories with FASTA format alignments. This script then calls TCOFFEE to merge the FASTA format alignments into multiple sequence alignments. One alignment is generated per DALI.txt output (that is, one DALI query searched against a database) such that there is an alignment generated for every single query. The script invokes TCOFFEE one more time to merge all such alignments into one final multiple-sequence alignment.


</ul>
</details>
