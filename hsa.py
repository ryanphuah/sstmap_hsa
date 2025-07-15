
import parmed as pmd
import argparse
import subprocess
import os
import numpy as np
import pandas as pd

def protein_quick_fix(parm_file,traj_file,lig_file,num_frames:str,output_prefix,cluster_file=None,dist=5):
    # Run HSA
    if cluster_file:
        result=subprocess.run("run_hsa -i "+parm_file+" -t "+traj_file+" -l "+lig_file+" -f "+num_frames+" -d "+dist+" -c "+cluster_file+" -o "+output_prefix,shell=True)
    else: 
        # result=subprocess.run(["run_hsa", "-i", {parm_file},"-t",{traj_file},"-l",{lig_file},"-f",{num_frames},"-d",{dist},"-o",{output_prefix}])
        
        result=subprocess.run("run_hsa -i "+parm_file+" -t "+traj_file+" -l "+lig_file+" -f "+num_frames+" -d "+dist+" -o "+output_prefix,shell=True)
    if result.returncode != 0:
        print("Error running HSA:")
        print(result.stderr)
    else:
        subprocess.run("mv SSTMap_HSA SSTMap_HSA_"+output_prefix,shell=True)
        print(f"HSA successful. Raw outputs written to SSTMAP_HSA_{output_prefix}")

    # Post HSA analysis
    # Classifying based on entropy/enthalpy
    with open('SSTMap_HSA_'+output_prefix+'/'+output_prefix+'_hsa_summary.txt', 'r') as f:
        rawhsa=pd.read_csv(f, sep=" ")
        rawhsa['dH']=rawhsa['Etot']+9.53
        rawhsa['dG']=rawhsa['Etot']-rawhsa['TStot']+9.53
    category=[]
    for i,j,k in zip(rawhsa['TStot'],rawhsa['dH'],range(len(rawhsa))):
        if i>-5.15: cat='Displaced'
        elif i<-5.1 and abs(j)<2.5: cat='Potential Replace'
        else: cat='Conserved'
        category.append(cat)
    rawhsa['HSA Category']=category
    
    # Matching to crystal waters
    if cluster_file:
        colspecs = [(0, 6), (6, 11), (12, 16), (16, 17), (17, 20), (21, 22), (22, 26),
            (26, 27), (30, 38), (38, 46), (46, 54), (54, 60), (60, 66), (76, 78),
            (78, 80)]
        names = ['ATOM', 'serial', 'name', 'altloc', 'resname', 'chainid', 'resseq',
         'icode', 'x', 'y', 'z', 'occupancy', 'tempfactor', 'element', 'charge']
        with open(cluster_file, 'r') as f:
            paperwats=pd.read_fwf(f,names=names, colspecs=colspecs)
            paperwats=paperwats[paperwats['ATOM']=='HETATM']
            paperwats=paperwats[['resseq','x','y','z','tempfactor']].dropna()
            paperwats=paperwats.reset_index()
        rawhsa['ResID']=paperwats['resseq']
        hsa_results=rawhsa[['ResID','x','y','z','occupancy','dH','TStot','dG','HSA Category']]
        hsa_results.to_csv(output_prefix+'_HSA_results.txt',index=True, sep=' ')
    else: 
        hsa_results=rawhsa[['x','y','z','occupancy','dH','TStot','dG','HSA Category']]
        hsa_results.to_csv(output_prefix+'_HSA_results.txt', index=True,sep=' ')

def main():
    parser= argparse.ArgumentParser(description='HSA Analysis of waters (TIP3P model) in protein pocket')
    parser.add_argument('-p','--parm_file',help='Input parameter file', required=True)
    parser.add_argument('-t','--traj_file',help='Input trajectory file',required=True)
    parser.add_argument('-l','--lig_file',help='Input ligand file',required=True)
    parser.add_argument('-f','--num_frames',help='Input number of frames',required=True,type=str)
    parser.add_argument('-c','--cluster_file',help='Input crystal waters to analyse',default=None)
    parser.add_argument('-d','--dist',help='Input distance from ligand to analyse (if no crystal waters specified)',default=5)
    parser.add_argument('-o','--output_prefix',help='Output prefix',required=True)

    # parser.add_argument('pdb_out',help='Output pdb',default=None)
    # parser.add_argument('water_model',help='Water Model (cSPCE/cTIP3P)')
    args = parser.parse_args()
    protein_quick_fix(**vars(args))

if __name__=="__main__":
    main()
    