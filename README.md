# HSA
<details>
<summary><b>Help documentation</b> </summary>

    usage: hsa.py [-h] -p PARM_FILE -t TRAJ_FILE -l LIG_FILE -f NUM_FRAMES
              [-c CLUSTER_FILE] [-d DIST] -o OUTPUT_PREFIX
    
    HSA Analysis of waters (TIP3P model) in protein pocket
    
    optional arguments:
      -h, --help            show this help message and exit
      -c CLUSTER_FILE, --cluster_file CLUSTER_FILE
                            Input crystal waters to analyse
      -d DIST, --dist DIST  Input distance from ligand to analyse (if no crystal
                            waters specified)
      -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                            Output prefix
    
    required arguments:
      -p PARM_FILE, --parm_file PARM_FILE
                            Input parameter file
      -t TRAJ_FILE, --traj_file TRAJ_FILE
                            Input trajectory file
      -l LIG_FILE, --lig_file LIG_FILE
                            Input ligand file (PDB)
      -f NUM_FRAMES, --num_frames NUM_FRAMES
                            Input number of frames

</details>

## Install SSTMap using conda
Either run the below lines or set up conda environment using .yml file provided
```
conda config --add channels omnia
conda config --add channels solvationtools
conda install sstmap
```
## Running test system 
Download trajectory from following link into working directory: https://drive.google.com/drive/folders/1dXeWMGwhTl2kC8hOqZBqS2nk2HEoQvNY?usp=drive_link

Run bash script as follows
```
sh ./test_hsa.sh
```
