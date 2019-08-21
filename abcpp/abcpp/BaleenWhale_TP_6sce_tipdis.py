import os,sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
sys.path.append('C:/Liang/abcpp_ms/abcpp')
from dvtraitsim_shared import DVTreeData, DVParamLiang
import dvtraitsim_cpp as dvcpp
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
heritability_list =[]
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

            length = 10 ** logTL
            obsZ = sorted(length)

            with open(dir_path + 'result_cluster/Est/predictsim%i_t5.csv' % count) as csv_file:
                csv2_reader = csv.reader(csv_file, delimiter=',')
                simZ = list(csv2_reader)
                del simZ[0]

            Z = np.array([float(i)  for sublist in simZ for i in sublist]).reshape((len(simZ),15))
            i, j = argsort2D(Z)
            Z = Z[i, j]
            tp_distance = np.linalg.norm(Z - obsZ, axis=1)
            total_distance = np.linalg.norm(Z,axis=1)
            print(np.mean(tp_distance))
            distance_list.append(tp_distance)
            ratio_dis.append(tp_distance/total_distance)
            timescaling_list.append(np.repeat(timescaling,len(simZ)))
            heritability_list.append(np.repeat(heritability,len(simZ)))


distance_list_flat = [item for sublist in distance_list for item in sublist]
ratio_dis_flat = [item for sublist in ratio_dis for item in sublist]

timescaling_list_flat = [item for sublist in timescaling_list for item in sublist]
heritability_list_flat = [item for sublist in heritability_list for item in sublist]


ss_list = {'distance':distance_list_flat,'timescale':timescaling_list_flat,
            'heritability':heritability_list_flat,'ratio':ratio_dis_flat}
ss_df = pd.DataFrame(ss_list)


f, axes = plt.subplots(1, 1,figsize = (9,12),sharex=True,sharey=True)
ax1 = sns.boxplot(x="timescale", y="distance",
                 hue="heritability", palette=["m", "g"],
                 data=ss_df,
                  linewidth=1, ax=axes[0], showfliers=False)
ax1.title.set_text('Dividing scalar: 1')


handles_gamma, labels_gamma = ax1.get_legend_handles_labels()

for ax in [ax1]:
    ax.set_ylabel('')
    ax.set_xlabel('')
    ax.legend_.remove()

axes[0].set_ylabel('Distance',fontsize=15)
axes[0].yaxis.set_label_position("left")


handles = handles_gamma #[ item for subhandle in [handles_gamma,handles_a,handles_nu] for item in subhandle]
labels = ['$h^2 = 0.5$', '$h^2 = 1$']
f.text(0.5, 0.04, 'Time scaling parameters', ha='center',fontsize=15)
# f.text(0.04, 0.5, 'Estimates', va='center', rotation='vertical',fontsize=15)
l = plt.legend(handles, labels, bbox_to_anchor=(0.65, 0.95), loc=2, borderaxespad=0.)
l.get_frame().set_linewidth(0.0)
