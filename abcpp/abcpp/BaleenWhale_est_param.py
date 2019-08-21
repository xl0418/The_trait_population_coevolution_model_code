import os,sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
sys.path.append('C:/Liang/abcpp_ms/abcpp')
import seaborn as sns
def argsort2D(X):
    return np.arange(len(X))[:, np.newaxis], np.argsort(X, axis=1)


def normalized_norm(x, y):
    diff_norm = np.linalg.norm(x - y, axis=1)
    max_err = np.nanmax(diff_norm)
    return diff_norm / max_err


dir_path = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/'
data_dir = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/result_cluster/'

obs_file = dir_path + 'treedata/'

with open(obs_file+'extantspecieslabels.csv') as csv_file:
    csv1_reader = csv.reader(csv_file, delimiter=',')
    extantlabels = list(csv1_reader)

with open(dir_path+'slater_length_data.csv') as csv_file:
    csv2_reader = csv.reader(csv_file, delimiter=',')
    lengthdata = list(csv2_reader)

extantlabels_array = np.array(extantlabels)[1:,1]
lengthdata_array = np.array(lengthdata)
length_index = []
gamma_list = []
a_list = []
nv_list = []
distance_list = []
ratio_dis = []
valid_length = []
timescaling_list = []
heri_list = []
count = 0
particle_size = 100
K=10e8
nu=1/(100*K)

for label in extantlabels_array:
    length_index.append(np.where(lengthdata_array[:,0]==label)[0][0])

logTL = lengthdata_array[length_index,1].astype(np.float)

timescale_vec = [20000,40000,80000]
heritability_vec = [1,0.5]


for timescaling_index in range(3):
        for heritability_index in range(2):
            print(count)
            count += 1
            timescaling = timescale_vec[timescaling_index]
            heritability = heritability_vec[heritability_index]

            data_name = data_dir + 'BWest_t%i_h%i.npy' % (int(timescaling), int(heritability_index)+1)
            est_data = np.load(data_name).item()
            generation = len(est_data['gamma'])
            population = len(est_data['gamma'][0])

            # last_fitness = est_data['fitness'][generation - 1]
            # q5 = np.argsort(last_fitness)[-population // 20]  # best 5%
            # fit_index = np.where(last_fitness > last_fitness[q5])[0]

            # mean of all samples
            # gamma_mean = np.mean(est_data['gamma'][generation-1])
            # a_mean = np.mean(est_data['a'][generation-1])
            # nv_mean = np.mean(est_data['nu'][generation-1])
            # mean of the top 5% samples

            gamma_vec = est_data['gamma'][generation-1]*1e7
            a_vec = est_data['a'][generation-1]*1e4
            gamma_list.append(gamma_vec.tolist())
            a_list.append(a_vec.tolist())
            timescaling_list.append(np.repeat(timescaling,population))
            heri_list.append(np.repeat(heritability,population))



gamma_list_flat = [item for sublist in gamma_list for item in sublist]
a_list_flat = [item for sublist in a_list for item in sublist]

timescaling_list_flat = [item for sublist in timescaling_list for item in sublist]
heri_list_flat = [item for sublist in heri_list for item in sublist]


ss_list = {'gamma':gamma_list_flat,'alpha':a_list_flat,'s':timescaling_list_flat,
            'heri':heri_list_flat}
ss_df = pd.DataFrame(ss_list)

def hexbin(x, y, color, **kwargs):
    cmap = sns.light_palette(color, as_cmap=True)
    plt.hexbin(x, y, gridsize=15, cmap=cmap, **kwargs)

g = sns.FacetGrid(ss_df, col="heri",  row="s",margin_titles=True)
g.map(hexbin, "gamma", "alpha", extent=[9, 11, 2, 3])
g.set_xlabels('$\gamma $ $(10^{-7})$')
g.set_ylabels('$\\alpha $ $(10^{-4})$')
axes = g.axes.flatten()
axes[0].set_title("$h^{2}=0.5$")
axes[1].set_title("$h^{2}=1$")
