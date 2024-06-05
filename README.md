![screenshot](figures/logos/snakedali_logo.png)

## OVERVIEW

<details>
<summary>1- What is Snakedali?</summary>

Snakedali is the Snakemake implementation of DaliLite v5 to align PDB queries to a pre-built Alphafold database.


</details>

<details>
<summary>2- Citation</summary>
<ul>

'Structure Guided Discovery of Ancestral CRISPR-Cas13 Ribonucleases '

</ul>
</details>

## GETTING STARTED

<details>
<summary>3. Dependencies</summary>
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

</ul>
</details>

<details>
<summary>4 Installation</summary>
<ul>

<details>
<summary>4.1 Database Download</summary>
<ul>

  - Download the pre-built database through __


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
  
  - A singularity/apptainer container is provided in this repository
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


## Snakedali PIPELINE STRUCTURE

<details>
<summary>Snakemake Profile</summary>
<ul>
  
  - Snakedali was designed to work with `(Sun Grid Engine) SGE` job scheduler
  - Set up the Snakemake profile: `/profile/config.yaml`
  
  - The default profile includes:
    - cluster job submission: `qsub -l h_rt={cluster.time} -j y -pe smp 4 -cwd`
    - cluster config path: `config/cluster.yaml`
    - rerun triggers: `mtime`
    - n jobs limit: `600`
   
  - Make sure to adjust the parameters above according to the house rules of your HPC.
  - If using the containers, make sure to uncomment the #singularity-args line on `profile/config.yaml`

</ul>
</details>
