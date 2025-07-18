# HSA
<details>
<summary><b>Help documentation</b> </summary>

    usage: hsa.py [-h] -p PARM_FILE -t TRAJ_FILE -l LIG_FILE -f NUM_FRAMES -w
              WATER_MODEL -s START_FRAME [-c CLUSTER_FILE] [-d DIST] -o
              OUTPUT_PREFIX

    HSA Analysis of waters in protein pocket
    
    optional arguments:
      -h, --help            show this help message and exit
      -c CLUSTER_FILE, --cluster_file CLUSTER_FILE
                            Input crystal waters to analyse
      -d DIST, --dist DIST  Input distance from ligand to analyse (if no crystal
                            waters specified). Default = 10
    
    required arguments:
      -p PARM_FILE, --parm_file PARM_FILE
                            Input parameter file
      -t TRAJ_FILE, --traj_file TRAJ_FILE
                            Input trajectory file
      -l LIG_FILE, --lig_file LIG_FILE
                            Input ligand file (PDB)
      -f NUM_FRAMES, --num_frames NUM_FRAMES
                            Total number of frames to process
      -w WATER_MODEL, --water_model WATER_MODEL
                            Water model used. Supported models are: TIP3P,
                            TIP4PEW, TIP4P, TIP5P, TIP3PFW, SPCE, SPCFW
      -s START_FRAME, --start_frame START_FRAME
                            Starting frame
      -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                            Output prefix

</details>

## Install SSTMap using conda
Either run the below lines or set up conda environment using .yml file provided
```
conda config --add channels omnia
conda config --add channels solvationtools
conda install sstmap
```
## Running SSTMap
Required input files:
- Amber parameter file (.prmtop)
- Trajectory file (.dcd)
- Ligand file (.pdb)

Pointers to note for input files:
- Restraints should be implemented on heavy atoms of protein during MD to generate trajectory file
- PDB file used to run MD should have waters at the end of the file/waters added separately into Modeller
- Distance of all waters indicated in cluster file to the ligand should be less than distance indicated by distance flag
- dH component is calculated with respect to TIP3P water bulk energy. dH and dG values of different water models would need to be adjusted accordingly 
- Output prefix needs to be changed every run to prevent potential overwriting of files

## HSA Outputs
In depth explanations of all output columns from HSA can be found at this link: http://sstmap.org/2017/05/09/undestanding-output/

Briefly, dH is calculated from the Etot column relative to bulk energy of the water model used. dG is then calculate from dH - TStot. Note that water bulk entropy is taken as 0.

Other interesting outputs from HSA include occupancy, number of hydrogen bonds as well as residues that form hydrogen bonds to each water. A PDB file of the probable configurations of waters in each cluster is also output.

## Post SSTMap HSA analysis
Probability of water beinng conserved and displaced is calculated based on Bayes' formula, with terms trained on SZMap dataset. Classification of conserved/displaced is done based on comparison of probability. 

Further determination of replaceable waters from conserved waters is done based on a magic cutoff value for dH of -2.5. dH values more negative than -2.5 are stable.

Note: SZMap dataset is relatively small, hence generalisability of this classifier may be limited.

## Running test system 
Download trajectory from following link into test directory: https://drive.google.com/drive/folders/1dXeWMGwhTl2kC8hOqZBqS2nk2HEoQvNY?usp=drive_link

Run bash script as follows
```
sh ./test_hsa.sh
```
Reference summary of HSA results provided for comparison. HSA will also write a folder of more detailed results for deeper analysis.
## Known list of issues
- bad_arrray_new_length/Segmentation fault (core dumped)/Aborted (core dumped)
    - Check distance of waters in cluster file to ligand file: should be less than distance specified in -d flag
    - There may be no waters in specified cluster centre. Check trajectory file accordingly (especially if periodic box was used for MD, may need to autoimage using cpptraj first)
    - Numpy might need to be downgraded. Latest supported version is 1.17

- Index out of range: Might occur if specifying 1 cluster only.

- Crystal waters go crazy in periodic box: Check original PDB file if there are waters found in between protein and ligand. Move to the end of the PDB file before rurnning MD
    - Side note: Restraints should be implemented on heavy atoms of protein during MD.
