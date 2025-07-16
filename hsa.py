
import parmed as pmd
import argparse
import subprocess
import os
import numpy as np
import pandas as pd
import scipy
import pickle as pl
import csv

def hsa(parm_file,traj_file,lig_file,num_frames:str,output_prefix,cluster_file=None,dist=10):
    # Check for clashing file names in cwd
    if os.path.isdir('SSTMap_HSA_'+output_prefix): 
        print('Error, directory already exists. Please change output prefix')
        return

    # Run HSA
    if cluster_file:
        result=subprocess.run("run_hsa -i "+parm_file+" -t "+traj_file+" -l "+lig_file+" -f "+num_frames+" -d "+dist+" -c "+cluster_file+" -o "+output_prefix,shell=True)
    else: 
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
    probc=pd.DataFrame(columns=['Prob Conserved'])
    probd=pd.DataFrame(columns=['Prob Displaced'])

    # Trained Bayes classifier from SZMap paper
    with open('displaced_kde.pkl','rb') as f: displaced_kde=pl.load(f)
    with open('conserved_kde.pkl','rb') as f: conserved_kde=pl.load(f)
    def prob(wat):
        d=36/54
        c=18/54
        xd=displaced_kde.pdf(wat)
        xc=conserved_kde.pdf(wat)
        x=xd*d+xc*c
        dx=xd*d/x
        dc=xc*c/x
        return dc,dx
    
    for i,j,k in zip(rawhsa['TStot'],rawhsa['dH'],range(len(rawhsa))):
        c,d=prob(i)
        probc.loc[len(probc)]=c
        probd.loc[len(probd)]=d
        if d>c: cat='Displaced'
        elif abs(j)<2.5: cat='Replaceable'
        else: cat='Conserved'
        category.append(cat)
    rawhsa['HSA_Category']=category
    rawhsa['Conserved_Probability']=probc
    rawhsa['Displaced_Probability']=probd

    # Matching to crystal waters if clusters specified
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
        hsa_results=rawhsa[['ResID','x','y','z','occupancy','dH','TStot','dG','Conserved_Probability','Displaced_Probability','HSA_Category']]
        hsa_results.to_csv(output_prefix+'_HSA_results.txt',index=True, sep=' ',quoting=csv.QUOTE_NONE)
    else: 
        hsa_results=rawhsa[['x','y','z','occupancy','dH','TStot','dG','Conserved Probability','Displaced Probability','HSA Category']]
        hsa_results.to_csv(output_prefix+'_HSA_results.txt', index=True,sep=' ',quoting=csv.QUOTE_NONE)



def main():
    parser= argparse.ArgumentParser(description='HSA Analysis of waters (TIP3P model) in protein pocket')
    required = parser.add_argument_group('required arguments')
    required.add_argument('-p','--parm_file',help='Input parameter file', required=True)
    required.add_argument('-t','--traj_file',help='Input trajectory file',required=True)
    required.add_argument('-l','--lig_file',help='Input ligand file (PDB)',required=True)
    required.add_argument('-f','--num_frames',help='Input number of frames',required=True,type=str)
    parser.add_argument('-c','--cluster_file',help='Input crystal waters to analyse',default=None)
    parser.add_argument('-d','--dist',help='Input distance from ligand to analyse (if no crystal waters specified). Default = 10',default=10)
    parser.add_argument('-o','--output_prefix',help='Output prefix',required=True)
    args = parser.parse_args()
    hsa(**vars(args))

if __name__=="__main__":
    main()
    