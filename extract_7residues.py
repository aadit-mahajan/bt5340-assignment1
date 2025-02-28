from biopandas import pdb as ppdb
import os
import logging
import concurrent.futures 
import Bio.PDB
from Bio.PDB import PDBParser
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import warnings

logging.basicConfig(level=logging.INFO, filename='debug.log', filemode='w', format='%(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def secondary_structure(pdb_path, chain_id):
    '''Extract secondary structure elements for a chain'''
    
    # extract info on secondary structure
    data = ppdb.PandasPdb().read_pdb(pdb_path)
    helix_raw = data.df["OTHERS"][data.df["OTHERS"]["record_name"] == "HELIX"]
    sheet_raw = data.df["OTHERS"][data.df["OTHERS"]["record_name"] == "SHEET"]

    if helix_raw.empty == False:
        helix_raw = helix_raw["entry"]
        # np array storing chain id, start and end index of every helix
        helix = []
        for row in helix_raw:
            content = row.rstrip().split()
            #check chain id
            if content[3] == chain_id and content[6] == chain_id:
                # very unclean solution to delete letter that is sometimes at the end of int in helix section 
                if content[4][-1].isalpha():
                    content[4] = content[4][:-1]
                if content[7][-1].isalpha():
                    content[7] = content[7][:-1]
                helix.append([int(content[4]), int(content[7])])
        helix = np.array(helix, dtype=int)
    else:
        helix = None

    if sheet_raw.empty == False:
        sheet_raw = sheet_raw["entry"]
        # np array storing chain id, start and end index of every sheet
        sheet = []
        for row in sheet_raw:
            content = row.rstrip().split()
            if content[4] == chain_id and content[7] == chain_id:
                # problem eg in file 2H5C, 5ZE8, 6FU9, 7ZGM, 8AH3, 8AYC, 8QQC
                if content[5][-1].isalpha():
                    content[5] = content[5][:-1]
                # needed for one where content[8]='86A-1'
                if content[8].startswith('-'):
                    content[8] = content[8]
                else:
                    content[8] = content[8].split('-')[0]
                if content[8][-1].isalpha():
                    content[8] = content[8][:-1]
                sheet.append([int(content[5]), int(content[8])])
        sheet = np.array(sheet, dtype=int)
    else:
        sheet = None

    return helix, sheet
    
def get_angles_7res(pdb_id, file_path, full_pdb_dir, chain_id, res = 'I'):
    try:
        print(f'starting {pdb_id}')
        full_pdb_path = f'{full_pdb_dir}/{pdb_id}.pdb'
        
        angles = []
        charged = []
        sec_struc = []
        
        charged_res = ['R', 'K', 'D', 'E']
        helix, sheet = secondary_structure(full_pdb_path, chain_id)
                
        for model in Bio.PDB.PDBParser().get_structure(pdb_id, file_path):
            for chain in model:
                if chain.id == chain_id:
                    poly = Bio.PDB.Polypeptide.Polypeptide(chain)
                    phiPsi_list = np.asarray(poly.get_phi_psi_list())
                    seq = np.array(poly.get_sequence(), dtype=str)
                    #ensure that 3 residues before and after X exist in the sequence
                    idx = np.where(seq[3:-3] == res)[0]
                    idx += 3
                    angles.extend(phiPsi_list[idx])
                    
                    
                    for i in idx:
                        # decide whether charged environment or not
                        # return True or False
                        test_charge = np.isin(seq[i-3:i+4], charged_res)
                        if np.any(test_charge):
                            charged.append(True)
                        else:
                            charged.append(False)                       
                    
                        # get secondary structure of sequence
                        if helix is not None:
                            inside_helix = [1 if helix[r, 0] < i < helix[r, 1] else 0 for r in range(np.shape(helix)[0])]
                            sum_helix = np.sum(inside_helix)
                        else:
                            sum_helix = 0
                        
                        if sheet is not None:
                            inside_sheet = [1 if sheet[r, 0] < i < sheet[r, 1] else 0 for r in range(np.shape(sheet)[0])]
                            sum_sheet = np.sum(inside_sheet)
                        else:
                            sum_sheet = 0
                            
                        if sum_helix > 0 and sum_sheet > 0:
                            sec_struc.append('hs')
                        elif sum_helix > 0:
                            sec_struc.append('h')
                        elif sum_sheet > 0:
                            sec_struc.append('s')
                        else:
                            sec_struc.append('c')         
                        
        return (angles, charged, sec_struc)
    except Exception as e:
        logging.error(f'Error in getting angles of 7 res sequence of {pdb_id}: {e}')

def get_all_angles_7res(pdb_files_dir, full_pdb_files_dir, output_file, res = 'I'):
    angles = []
    charged = []
    sec_struct = []
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = executor.map(get_angles_7res, 
                     [file.split('_')[0] for file in os.listdir(pdb_files_dir)], 
                     [f'{pdb_files_dir}/{file}' for file in os.listdir(pdb_files_dir)], 
                     [full_pdb_files_dir for _ in range(len(os.listdir(pdb_files_dir)))],
                     [file.split('_')[1].split('.')[0] for file in os.listdir(pdb_files_dir)], 
                     [res for _ in range(len(os.listdir(pdb_files_dir)))]
                      )
    
    # concatenate the results
    for result in results:
        angles.extend(result[0])
        charged.extend(result[1])
        sec_struct.extend(result[2])
    data_out = pd.DataFrame(angles, columns=['phi', 'psi'])
    data_out["charged"] = charged
    data_out["sec_struct"] = sec_struct
    data_out.to_csv(output_file, index=False)
    print('Saved:', output_file)
  

if __name__ == '__main__':
    
    warnings.filterwarnings("ignore")
    
    ANGLES_DATA_OUT_DIR = "./angles"
    PDB_FILES_DIR = "./pdb_files"
    CLEANED_PDB_DIR = "./pruned_pdb_files"
    ALL_LENGTHS = []
    
    logger.info('Starting extracting 7 residue sequences')
    get_all_angles_7res(CLEANED_PDB_DIR, PDB_FILES_DIR, "I_angles_7res.csv")
    logger.info('All angles of 7 residue sequence calculated')
    
