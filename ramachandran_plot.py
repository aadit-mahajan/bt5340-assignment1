import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os



def ramachandran_plot(phi, psi, save_file, title=""):
    # convert from rad to deg
    phi = np.rad2deg(phi)
    psi = np.rad2deg(psi)
    
    # bins of width 3°x3°
    bin_values = np.linspace(-180, 180, 121)
    
    # plot 2D histogram
    plt.hist2d(phi, psi, bins=[bin_values, bin_values], cmap='gist_ncar_r')
    plt.xlabel(r"$\phi$ [deg]")
    plt.ylabel(r"$\psi$ [deg]")
    plt.title(title)
    plt.savefig(save_file)
    plt.show()

if __name__ == '__main__':
    
    FIG_DIR = "./plots"
    os.makedirs(FIG_DIR, exist_ok=True)
    
    # Ramachandran Plot for all ILEs
    all_angles = np.genfromtxt("I_angles.csv", dtype=float, delimiter=",", skip_header=1)
    ramachandran_plot(all_angles[:, 0], all_angles[:, 1], f'{FIG_DIR}/all_ILE.png', "all I")

    # Isoleucins included in the 7 residue sequence abcIdef
    data_7res = pd.read_csv("I_angles_7res.csv")
    # entire data
    ramachandran_plot(data_7res["phi"], data_7res["psi"], f"{FIG_DIR}/all_ILE_7res.png", "I in 7 res")
    # charged res in sequence
    idx_charge = data_7res[data_7res["charged"] == True].index
    ramachandran_plot(data_7res.iloc[idx_charge]["phi"], data_7res.iloc[idx_charge]["psi"], 
                      f"{FIG_DIR}/charged_7res.png", "charged environment")
    # no charged res in sequence
    idx_nocharge = data_7res[data_7res["charged"] == False].index
    ramachandran_plot(data_7res.iloc[idx_nocharge]["phi"], data_7res.iloc[idx_nocharge]["psi"], 
                      f"{FIG_DIR}/noncharged_7res.png", "non charged environment")
    # indices based on secondary structure
    idx_helix = data_7res[data_7res["sec_struct"] == 'h'].index
    idx_sheet = data_7res[data_7res["sec_struct"] == 's'].index
    idx_coil = data_7res[data_7res["sec_struct"] == 'c'].index
    idx_hs = data_7res[data_7res["sec_struct"] == 'hs'].index
    # I in helix
    ramachandran_plot(data_7res.iloc[idx_helix.union(idx_hs)]["phi"], data_7res.iloc[idx_helix.union(idx_hs)]["psi"], 
                    f"{FIG_DIR}/helix_7res.png", "helix")
    # I in sheet
    ramachandran_plot(data_7res.iloc[idx_sheet.union(idx_hs)]["phi"], data_7res.iloc[idx_sheet.union(idx_hs)]["psi"], 
                    f"{FIG_DIR}/sheet_7res.png", "sheet")
    # I in either helix or sheet
    ramachandran_plot(data_7res.iloc[idx_sheet.union(idx_helix).union(idx_hs)]["phi"], 
                    data_7res.iloc[idx_sheet.union(idx_helix).union(idx_hs)]["psi"], f"{FIG_DIR}/noncoil_7res.png", "helix and sheet")
    # I in coil
    ramachandran_plot(data_7res.iloc[idx_coil]["phi"], data_7res.iloc[idx_coil]["psi"], f"{FIG_DIR}/coil_7res.png", "coiled")

