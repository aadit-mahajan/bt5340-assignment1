import os
import requests
import concurrent.futures

def fetch_pdb_files(pdb_id, out_dir):
    '''
    Using the pdb_id, fetch the pdb file from the RCSB PDB database
    '''
    try:
        filename = f'{pdb_id}.pdb'
        url = f'https://files.rcsb.org/download/{filename}'
        response = requests.get(url)
        with open(os.path.join(out_dir, filename), 'wb') as f:
            f.write(response.content)
    except Exception as e:
        print(e)

# checking difference in downloaded and required PDBs
def check_missing_pdbs(): 
    downloaded_pdbs = os.listdir(OUT_DIR)
    downloaded_pdbs = [pdb.split('.')[0] for pdb in downloaded_pdbs]

    missing_pdbs = list(set(pdb_ids) - set(downloaded_pdbs))
    print(f'Number of missing PDBs: {len(missing_pdbs)}')
    return missing_pdbs

if __name__ == '__main__':
    with open('server_data_pdbs', 'r') as f:
        lines = f.readlines()

    pdb_ids = [line.strip()[:4] for line in lines]    # pdb_id is in the first 4 characters of each line. 
    pdb_ids = list(set(pdb_ids))

    print(f'Number of PDBs to fetch: {len(pdb_ids)}')

    OUT_DIR = './pdb_files'     
    os.makedirs(OUT_DIR, exist_ok=True)

    # executing parallelly to speed up I/O
    with concurrent.futures.ThreadPoolExecutor() as executor:
        executor.map(fetch_pdb_files, pdb_ids, [OUT_DIR]*len(pdb_ids))

    # retrying missing PDBs
    missing_pdbs = check_missing_pdbs()
    for pdb_id in missing_pdbs:
        fetch_pdb_files(pdb_id, OUT_DIR)
    
    print('All done!')


