# bt5340-assignment1

***Run order***
fetch_pdb.py --> angles_calc.py

***fetch_pdb.py***
Reads the list of pdbs from the server_data_pdbs file and fetches the files from the RCSB server. 
All PDB files are stored in a directory called 'pdb_files' (has been added to .gitignore).
Wall-time for this script can vary based on the internet connection speed. 

***angles_calc.py***
This reads the pdb files and their corresponding chain-ids. 
It prunes the pdbs present in the pdb_files directory to keep only the data for the chain of interest. 
Phi-psi angles are calculated and stored in csv files in the 'angles' directory with the file naming format PDBID_CHAIN.csv. 
These csv files contain the data in the format: (phi, psi, residue) 

Approximate wall-time for angles_calc.py is around 60 min. 

