import os
import Bio.PDB
import numpy as np
import pandas as pd
import warnings
from biopandas import pdb as ppdb
import concurrent.futures 
import time
import logging

logging.basicConfig(level=logging.INFO, filename='debug.log', filemode='w', format='%(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# metadata = np.array([pdb_ids, chain_ids]).T
# print('PDBs:', len(pdb_ids))
# print(pdb_ids[:5])

def clean_pdb(pdb_file_path, chain_id, output_dir):
    '''
    extract only the chain of interest from the pdb file
    '''
    try:
        aa_list = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']

        pdb_id = os.path.basename(pdb_file_path).split('.')[0]
        data = ppdb.PandasPdb().read_pdb(pdb_file_path)
        clean_data = ppdb.PandasPdb()

        clean_data.df['ATOM'] = data.df['ATOM'][data.df['ATOM']['chain_id'] == chain_id]
        # print(clean_data.df['ATOM']['residue_name'])
        clean_data.df['ATOM'] = clean_data.df['ATOM'][clean_data.df['ATOM']['residue_name'].isin(aa_list)]

        output_path = os.path.join(output_dir, f'{pdb_id}_{chain_id}.pdb')
        ppdb.PandasPdb.to_pdb(clean_data, output_path)

        return clean_data
    
    except Exception as e:
        logging.error(f'Error in cleaning pdb at path {pdb_file_path}: {e}')

def get_angles(pdb_id, file_path, chain_id):
    '''
    Calculates the phi and psi angles for each residue in the chain using the BioPython library

    '''
    try:
        angles = []
        for model in Bio.PDB.PDBParser().get_structure(pdb_id, file_path):
            for chain in model:
                if chain.id == chain_id:
                    poly = Bio.PDB.Polypeptide.Polypeptide(chain)
                    # print ("Model %s Chain %s" % (str(model.id), str(chain.id)))
                    phiPsi_list = poly.get_phi_psi_list()
                    seq = poly.get_sequence()

                    data = []
                    for i in range(len(phiPsi_list)):
                        phi, psi = phiPsi_list[i]
                        # ignoring the entries which are HETATM/ANISOU (these turn up as None vals)
                        if phi is None and psi is None:
                            continue
                        elif phi is None:
                            phi = 0
                        elif psi is None:
                            psi = 0
                        data.append([phi, psi, seq[i]])
                    angles.extend(data)
                    break

        return angles
    
    except Exception as e:
        logging.error(f'Error in getting angles: {e}')

def generate_angle_data(pdb_file_path, output_dir):
    '''
    orchestrating angle data calculation
    '''
    try:
        if os.path.exists(output_dir):
            pdb_id = os.path.basename(pdb_file_path).split('.')[0].split('_')[0]
            chain = os.path.basename(pdb_file_path).split('.')[0].split('_')[1]
            angles = get_angles(pdb_id, file_path=pdb_file_path, chain_id=chain)

            if len(angles) > 0:
                phi_angles, psi_angles, sequence = zip(*angles)
                data = pd.DataFrame({'phi': phi_angles, 'psi': psi_angles, 'sequence': sequence})
                data_file_name = os.path.join(output_dir, pdb_id + '.csv')
                data.to_csv(data_file_name, index=False)
                print('Saved:', data_file_name)
            else:
                print('No angles found for:', pdb_id)
        else:
            print('Output directory not found:', output_dir)
    except Exception as e:
        logging.error(f'Error in generating angle data: {e}')

def prune_pdb_files(file_path, output_dir):
    '''
    Orchestrating the cleaning of pdb files
    '''
    try:
        for file in os.listdir(file_path):
            if file.endswith('.pdb'):
                pdb_id = file.split('.')[0]
                chain_id = chain_ids[pdb_ids.index(pdb_id)]
                clean_pdb(f'{PDB_FILES_DIR}/{file}', chain_id, output_dir=output_dir)
                print('Cleaned:', pdb_id)
    except Exception as e:
        logging.error(f'Error in pruning pdb files: {e}')

    # prune_pdb_files('./pdb_files')
    # generate_angle_data('./pdb_files/4N2P.pdb', 'angles')

    # Expected time for cleaning: 26.33480230967204 minutes (sequential)
    # Expected time for angle calculation: 16.18500550587972 minutes (sequential)

def get_angles_res(pdb_id, file_path, chain_id, res = 'I'):
    angles = []
    for model in Bio.PDB.PDBParser().get_structure(pdb_id, file_path):
        for chain in model:
            if chain.id == chain_id:
                poly = Bio.PDB.Polypeptide.Polypeptide(chain)
                phiPsi_list = np.asarray(poly.get_phi_psi_list())
                seq = np.array(poly.get_sequence(), dtype=str)
                idx = np.where(seq == res)[0]
                angles.extend(phiPsi_list[idx])

    return angles

def get_all_angles_res(pdb_files_dir, output_file, res = 'I'):
    data = []
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = executor.map(get_angles_res, 
                     [file.split('_')[0] for file in os.listdir(pdb_files_dir)], 
                     [f'{pdb_files_dir}/{file}' for file in os.listdir(pdb_files_dir)], 
                     [file.split('_')[1].split('.')[0] for file in os.listdir(pdb_files_dir)], 
                     [res for _ in range(len(os.listdir(pdb_files_dir)))]
                      )
    
    # concatenate the results
    for result in results:
        data.extend(result)

    data = pd.DataFrame(data, columns=['phi', 'psi'])
    data.to_csv(output_file, index=False)
    print('Saved:', output_file)

def check_missing_files(pdb_list, file_dir):
    '''
    Check for missing files in the output directory
    '''
    processed_files = os.listdir(file_dir)
    processed_files = [file.split('.')[0] for file in processed_files]
    diff = list(set(pdb_list) - set(processed_files))

    return diff

if __name__ == '__main__':

    warnings.filterwarnings("ignore") # this was done only for the sake of the pandaspdb I/O warnings. 

    # directories
    ANGLES_DATA_OUT_DIR = "./angles"
    PDB_FILES_DIR = "./pdb_files"
    CLEANED_PDB_DIR = "./pruned_pdb_files"

    # create directories if they don't exist
    os.makedirs(CLEANED_PDB_DIR, exist_ok=True)
    os.makedirs(ANGLES_DATA_OUT_DIR, exist_ok=True)

    pdb_ids = []
    chain_ids = []
    with open('server_data_pdbs', 'r') as f:     # server_data_pdbs is the file containing the pdb_ids and chain_ids (extracted from Dunbrack's server)
        for line in f:
            pdb_ids.append(line.strip()[:4])
            chain_ids.append(line.strip()[4])

    t1 = time.time()
    logger.info('Starting the process')

    # cleaning pdb files (parallelly executed for speed up)
    with concurrent.futures.ThreadPoolExecutor() as executor:
        executor.map(prune_pdb_files, [PDB_FILES_DIR], [CLEANED_PDB_DIR])

    logger.info('Pruning done')
    logger.info(f'Files in pruned_pdb_files: {len(os.listdir(CLEANED_PDB_DIR))}')
    
    pdb_file_list = os.listdir(CLEANED_PDB_DIR)
    print('Total files:', len(pdb_file_list))
    
    # angle calculation (parallelly executed for speed up)
    logger.info('Starting angle calculation')
    with concurrent.futures.ThreadPoolExecutor() as executor:

        executor.map(generate_angle_data, 
            [f'{CLEANED_PDB_DIR}/{file}' for file in os.listdir(CLEANED_PDB_DIR)], 
            ['angles' for _ in range(len(os.listdir(CLEANED_PDB_DIR)))]
            )

    # confirm that that all files have been processed
    print('Files in angles:', len(os.listdir(ANGLES_DATA_OUT_DIR)))
    print('Files in pruned_pdb_files:', len(os.listdir(CLEANED_PDB_DIR)))
    if len(os.listdir(ANGLES_DATA_OUT_DIR)) == len(os.listdir(CLEANED_PDB_DIR)):
        print('All files have been processed\n')
        logger.info('All files have been processed')
    else:
        print('Some files are missing.')
        logger.warning('All files not processed. There maybe a mismatch in the number of files processed due to redundancy in the data from the server!')

    # retrying missing files
        diff = check_missing_files(pdb_ids, ANGLES_DATA_OUT_DIR)
        for pdb_id in diff:
            chain_id = chain_ids[pdb_ids.index(pdb_id)]
            clean_pdb(f'{PDB_FILES_DIR}/{pdb_id}.pdb', chain_id, output_dir=CLEANED_PDB_DIR)
            print('Cleaned:', pdb_id)
            generate_angle_data(f'{CLEANED_PDB_DIR}/{pdb_id}_{chain_id}.pdb', ANGLES_DATA_OUT_DIR)
            print('Processed:', pdb_id)

    t2 = time.time()
    print('Time taken:', (t2-t1)/60, 'minutes')
    print('All done!\n-----------------')



