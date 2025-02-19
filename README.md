# bt5340-assignment1

***Run order***:   
fetch_pdb.py --> angles_calc.py --> probability_distribution.py
                                --> extract_7residues.py --> ramachandran_plot.py

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

Both scripts after running generate a _debug.log_ file which has logs for starting the execution of a process and errors (if any) for each function call. 

***extract_7residues.py***
This is very similar to the calculation of the phi psi angles in angles_calc.py but only the angles of Isoleucins are calculated which have at least 3 residus in front and after it.
Script is very similar to angles_calc.py. For computational efficiency both scripts can be combined. 
Additional information is extracted and stored: 
    charge: whether at least one of the 6 neighbouring residues is Aspartate, Arginine, Glutamate or Lysine
    secondary structure information: checks with the unpruned pdb file whether the Isolycine of the 7 residue sequence is part of a helix or a strand

***ramachandran_plot.py***
This code generates and stores histograms of the phi psi angles of Isoleucin for different datasets (all angles, only Is present in the 7 residue sequence) and their subsets based on charge of the neighbouring residues and the secondary structure.

***probability_distribution.py***
Histograms of phi and psi angles are created according to a random and a uniform distribution. 
The probability density of these reference datasets and the actual data for Isoleucin is compared with the Kullback-Liebler divergence.

to do:
inspect 9IJC.pdb

