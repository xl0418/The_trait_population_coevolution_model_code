import sys, os
import csv
sys.path.append('C:/Liang/abcpp_ms5/abcpp')
from Dvtraitsim_TVP import DVSimTVP
from dvtraitsim_shared import DVTreeData, DVParamLiang
from Dvtraitsim_TV import DVSimTV
from Dvtraitsim_TVM import DVSimTVM

import numpy as np
import matplotlib.pyplot as plt



r = 1  # growth rate
Vmax_TVP = 332
Vmax_TV = 92.8
Vmax_TVM = 48.0

gamma_TVP = 9.954e-7
gamma_TV = 4.072e-6
gamma_TVM = 4.311e-9

a_TVP = 2.365e-4
a_TV = 2.437e-5
a_TVM = 8.577e-4

nu_TVP = 3.038e-3
nu_TV = 2.933e-3
nu_TVM = 4.955e-3

scalar = 20000
K=10e5
K_TVM = 1e12

# let's try to find a true simulation:
model_label = ['TVP','TV','TVM']

#full tree and pruned tree directory
dir_path = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/'

files = dir_path + 'treedata/'

meantrait = 1300

# trait evolution plot

count = 0

td = DVTreeData(path=files, scalar=scalar)
for replicate in range(100):
    print('Simulating TVP model...')
    obs_param_TVP = DVParamLiang(gamma=gamma_TVP, a=a_TVP, K=K,h=1, nu=nu_TVP, r=r, theta=meantrait,V00=.1,V01=.1,
                                 Vmax=Vmax_TVP, inittrait=meantrait, initpop=1e5,
                        initpop_sigma=10.0, break_on_mu=False)
    simresult_TVP = DVSimTVP(td,obs_param_TVP)
    if simresult_TVP['sim_time'] == td.sim_evo_time:
        pic_TVP = 0
        break
    else:
        pic_TVP=1

for replicate in range(100):
    print('Simulating TV model...')

    obs_param_TV = DVParamLiang(gamma=gamma_TV, a=a_TV, K=K,h=1, nu=nu_TV, r=r, theta=meantrait,V00=.1,V01=.1,
                                Vmax=Vmax_TV, inittrait=meantrait, initpop=1e5,
                        initpop_sigma=10.0, break_on_mu=False)
    simresult_TV = DVSimTV(td,obs_param_TV)
    if simresult_TV['sim_time'] == td.sim_evo_time:
        pic_TV = 0
        break
    else:
        pic_TV=1

for replicate in range(100):
    print('Simulating TVM model...')

    obs_param_TVM = DVParamLiang(gamma=gamma_TVM, a=a_TVM, K=K_TVM,h=1, nu=nu_TVM, r=r, theta=meantrait,V00=.1,V01=.1,
                                 Vmax=Vmax_TVM, inittrait=meantrait, initpop=1e5,
                        initpop_sigma=10.0, break_on_mu=False)
    simresult_TVM = DVSimTVM(td,obs_param_TVM)
    if simresult_TVM['sim_time'] == td.sim_evo_time:
        pic_TVM = 0
        break
    else:
        pic_TVM=1
# if pic==0:
evo_time, total_species = simresult_TVP['N'].shape
evo_time = evo_time - 1
trait_TVP = simresult_TVP['Z']
trait_TV = simresult_TV['Z']
row_TV,col_TV = np.where(trait_TV == 0)
trait_TV[row_TV,col_TV] = np.nan
trait_TVM = simresult_TVM['Z']
row_TVM,col_TVM = np.where(trait_TVM == 0)
trait_TVM[row_TVM,col_TVM] = np.nan

num_lines = total_species
x = np.arange(evo_time+1)

f1, axes1 = plt.subplots(3, 1, figsize=(9, 9),sharey=True,sharex=True) #

labels = []
for i in range(1, num_lines + 1):
    axes1[0].plot(x, trait_TVP[:, i - 1])
    axes1[1].plot(x, trait_TV[:, i - 1])
    axes1[2].plot(x, trait_TVM[:, i - 1])

for img in range(3):
    axes1[img].xaxis.set_major_locator(plt.NullLocator())
    axes1[img].ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    axes1[img].set_ylabel(model_label[img])
    axes1[img].yaxis.set_label_position("right")


def argsort2D(X):
    return np.arange(len(X))[:, np.newaxis], np.argsort(X, axis=1)

with open(files+'extantspecieslabels.csv') as csv_file:
    csv1_reader = csv.reader(csv_file, delimiter=',')
    extantlabels = list(csv1_reader)

with open(dir_path+'slater_length_data.csv') as csv_file:
    csv2_reader = csv.reader(csv_file, delimiter=',')
    lengthdata = list(csv2_reader)

extantlabels_array = np.array(extantlabels)[1:,1]
lengthdata_array = np.array(lengthdata)
length_index = []
for label in extantlabels_array:
    length_index.append(np.where(lengthdata_array[:,0]==label)[0][0])

logTL = lengthdata_array[length_index,1].astype(np.float)
length = 10**logTL
obsZ = length
s = np.argsort(obsZ)
obsZ = obsZ[s]

trait_TVP_sorted = np.array(sorted(trait_TVP[-1]))
diff_norm_TVP = np.linalg.norm(trait_TVP_sorted - obsZ)

trait_TV_sorted = np.array(sorted(trait_TV[-1]))
diff_norm_TV = np.linalg.norm(trait_TV_sorted - obsZ)

np.linalg.norm(trait_TV_sorted - trait_TVP_sorted)
dir_fig = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/traittree'


gap_tvp = np.array([trait_TVP_sorted[i+1]-trait_TVP_sorted[i] for i in range(len(trait_TVP_sorted)-1)])
gap_tv = np.array([trait_TV_sorted[i+1]-trait_TV_sorted[i] for i in range(len(trait_TV_sorted)-1)])
gap_obs = np.array([obsZ[i+1]-obsZ[i] for i in range(len(obsZ)-1)])
#
# f1.savefig(dir_fig+'TP2w.png')
# plt.close(f1)

