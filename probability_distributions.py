
import numpy as np
from scipy.special import kl_div

'''
This file loads the phi and psi angles of all Isoleucins in the database and compares
the probability distribution to a random and a uniform one using Kullback-Leibler divergence.
'''

# load data
all_angles = np.genfromtxt("I_angles.csv", dtype=float, delimiter=",", skip_header=1)
# convert from radiant to degree
all_angles = np.rad2deg(all_angles)

#generate histogram
bin_values = np.linspace(-180, 180, 121)
hist = np.histogram2d(all_angles[:, 0], all_angles[:, 1], bins=[bin_values, bin_values])

# random distribution
np.random.seed(0)
data_random = np.random.rand(10**6, 2)
data_random = data_random * 360 - 180
hist_random = np.histogram2d(data_random[:, 0], data_random[:, 1], bins=[bin_values, bin_values])

# uniform distribution
hist_uniform = np.ones_like(hist[0])

# normalize distributions
hist_norm = hist[0] / np.sum(hist[0])
hist_random_norm = hist_random[0] / np.sum(hist_random[0])
hist_uniform_norm = hist_uniform / np.sum(hist_uniform)

# Kullback-Leibler divergence 
kl_random = kl_div(hist_norm, hist_random_norm)
kl_random = np.sum(kl_random)
print(f"KL divergence, random distribution: {kl_random}")

kl_uniform = kl_div(hist_norm, hist_uniform_norm)
kl_uniform = np.sum(kl_uniform)
print(f"KL divergence, uniform distribution: {kl_uniform}")