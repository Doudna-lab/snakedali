![screenshot](figures/logos/snakedali_logo.png)

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

</details>
