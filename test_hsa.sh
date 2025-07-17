cd test
python ../hsa.py -p 1h0s_protein_crystalwat_init_5ns.prmtop -t 1h0s_5ns_autoimaged.dcd -l 1h0s_ligand_lastframealigned.pdb -f 5000 -s 1 -d 10 -o 1h0s -c 1h0s_paperwats.pdb -w TIP3P
python ../hsa.py -p 5std_protein_crystalwat_init_5ns.prmtop -t 5std_5ns_autoimaged.dcd -l 5std_ligand_lastframealigned.pdb -f 5000 -s 1 -d 10 -o 5std -c 5std_paperwats.pdb -w TIP3P
python ../hsa.py -p 2qwj_protein_crystalwat_init_5ns.prmtop -t 2qwj_5ns_autoimaged.dcd -l 2qwj_ligand_lastframealigned.pdb -f 5000 -s 1 -d 10 -o 2qwj -c 2qwj_paperwats.pdb -w TIP3P
