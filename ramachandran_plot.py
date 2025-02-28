import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

def ramachandran_plot(phi, psi, save_file, title="", norm_factor=1):
    # convert from rad to deg
    phi = np.rad2deg(phi)
    psi = np.rad2deg(psi)
    
    # bins of width 3°x3°
    bin_values = np.linspace(-180, 180, 121)
    
    # plot 2D histogram
    # plt.hist2d(phi, psi, bins=[bin_values, bin_values], cmap='gist_ncar_r')
    fig = plt.figure()
    h, xedges, yedges = np.histogram2d(phi, psi, bins=[bin_values, bin_values], density=True)
    h_normalized = h*norm_factor
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    im = plt.imshow(h_normalized.T, extent=extent, origin='lower', cmap='gist_ncar_r', vmin=0, vmax=1/norm_factor)
    fig.colorbar(im)
    plt.xlabel(r"$\phi$ [deg]")
    plt.ylabel(r"$\psi$ [deg]")
    plt.title(title)
    plt.savefig(save_file)
    # plt.show()

if __name__ == '__main__':
    
    FIG_DIR = "./plots"
    os.makedirs(FIG_DIR, exist_ok=True)
    
    all_angles = np.genfromtxt("I_angles.csv", dtype=float, delimiter=",", skip_header=1)
    # calculating common normalization factor for all plots using the all plot
    bin_values = np.linspace(-180, 180, 121)
    h, xedges, yedges = np.histogram2d(all_angles[:, 0], all_angles[:, 1], bins=[bin_values, bin_values], density=True)
    norm_factor = 1/np.max(h)

    # Ramachandran Plot for all ILEs
    ramachandran_plot(all_angles[:, 0], all_angles[:, 1], f'{FIG_DIR}/all_ILE.png', "all I", norm_factor=norm_factor)

    # Isoleucins included in the 7 residue sequence abcIdef
    data_7res = pd.read_csv("I_angles_7res.csv")
    # entire data
    ramachandran_plot(data_7res["phi"], data_7res["psi"], f"{FIG_DIR}/all_ILE_7res.png", "I in 7 res", norm_factor=norm_factor)
    # charged res in sequence
    idx_charge = data_7res[data_7res["charged"] == True].index
    ramachandran_plot(data_7res.iloc[idx_charge]["phi"], data_7res.iloc[idx_charge]["psi"], 
                      f"{FIG_DIR}/charged_7res.png", "charged environment", norm_factor=norm_factor)
    # no charged res in sequence
    idx_nocharge = data_7res[data_7res["charged"] == False].index
    ramachandran_plot(data_7res.iloc[idx_nocharge]["phi"], data_7res.iloc[idx_nocharge]["psi"], 
                      f"{FIG_DIR}/noncharged_7res.png", "non charged environment", norm_factor=norm_factor)
    # indices based on secondary structure
    idx_helix = data_7res[data_7res["sec_struct"] == 'h'].index
    idx_sheet = data_7res[data_7res["sec_struct"] == 's'].index
    idx_coil = data_7res[data_7res["sec_struct"] == 'c'].index
    idx_hs = data_7res[data_7res["sec_struct"] == 'hs'].index
    # I in helix
    ramachandran_plot(data_7res.iloc[idx_helix.union(idx_hs)]["phi"], data_7res.iloc[idx_helix.union(idx_hs)]["psi"], 
                    f"{FIG_DIR}/helix_7res.png", "helix", norm_factor=norm_factor)
    # I in sheet
    ramachandran_plot(data_7res.iloc[idx_sheet.union(idx_hs)]["phi"], data_7res.iloc[idx_sheet.union(idx_hs)]["psi"], 
                    f"{FIG_DIR}/sheet_7res.png", "sheet", norm_factor=norm_factor)
    # I in either helix or sheet
    ramachandran_plot(data_7res.iloc[idx_sheet.union(idx_helix).union(idx_hs)]["phi"], 
                    data_7res.iloc[idx_sheet.union(idx_helix).union(idx_hs)]["psi"], f"{FIG_DIR}/noncoil_7res.png", "helix and sheet", norm_factor=norm_factor)
    
    # I in coil
    ramachandran_plot(data_7res.iloc[idx_coil]["phi"], data_7res.iloc[idx_coil]["psi"], f"{FIG_DIR}/coil_7res.png", "coiled", norm_factor=norm_factor)

