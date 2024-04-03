# New Improved Dali

## CLONING THIS REPOSITORY

<details>
<summary>1- Standard Git</summary>
<ul>

  Clone repository files
```
git clone https://github.com/Doudna-lab/nidali.git
```
</details>

<details>
<summary>2- Git LFS</summary>
<ul>
  -1.1 Install Git LFS to pull apptainer containers
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
</details>

## NIDALI PIPELINE STRUCTURE

<details>
<summary>Snakemake Profile</summary>
<ul>
  
  - NIDALI is currently set up to work on `@dev1.wynton.ucsf.edu`
  - Set up the Snakemake profile: `/profile/config.yaml`
  - The default profile includes:
    - cluster job submission: `qsub -l h_rt={cluster.time} -j y -pe smp 4 -cwd`
    - cluster config path: `config/cluster.yaml`
    - rerun triggers: `mtime`
    - singularity arguments: `--bind /usr/lib64/openmpi/ --bind /scratch --bind /wynton/home/doudna/bellieny-rabelo/nidali_output --bind /wynton/home/doudna/bellieny-rabelo/nidali_db`
    - n jobs limit: `400`

</details>
