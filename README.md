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
                            waters specified). Default = 5
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
Reference summary of HSA results provided for comparison. HSA will also write a folder of more detailed results for deeper analysis.
## Known list of issues
- bad_arrray_new_length: No waters found around specified cluster
    - Check distance of waters in cluster file to ligand file: should be less than distance specified in -d flag
    - There may be no waters in specified cluster centre. Check trajectory file accordingly
 
- Segmentation fault: numpy might need to be downgraded. Latest supported version is 1.17

- Output prefix needs to be changed every run to prevent overriding of files

- Index out of range: Might occur if specifying 1 cluster only.

- Crystal waters go crazy in periodic box: Check original PDB file if there are waters found in between protein and ligand. Move to the end of the PDB file before rurnning MD
    - Side note: Restraints should be implemented on heavy atoms of protein during MD.
