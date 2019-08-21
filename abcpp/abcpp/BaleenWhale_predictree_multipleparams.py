import os,sys
import platform
if platform.system()=='Windows':
    sys.path.append('C:/Liang/abcpp_ms5/abcpp')
elif platform.system()=='Darwin':
    sys.path.append('/Users/dudupig/Documents/GitHub/Code/Pro2/Python_p2')
import numpy as np
from dvtraitsim_shared import DVTreeData, DVParamLiang
import dvtraitsim_cpp as dvcpp
import pandas as pd
import csv

def argsort2D(X):
    return np.arange(len(X))[:, np.newaxis], np.argsort(X, axis=1)



dir_path = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/'
data_dir = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/result_cluster/result_0617/'

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
vm_list = []
theta_list = []
distance_list = []
valid_length = []
timescaling_list = []
heritability_list =[]
count = 8
particle_size = 1000
K=10e5
nu=1/(100*K)

for label in extantlabels_array:
    length_index.append(np.where(lengthdata_array[:,0]==label)[0][0])

logTL = lengthdata_array[length_index,1].astype(np.float)
empirical_data = {'species':extantlabels_array,'TL':logTL}
ed_df = pd.DataFrame(empirical_data)
ed_df.to_csv(data_dir+'emp.csv', sep=',', index=False)

timescale_vec = [20000,40000,60000,80000,100000]
heritability_vec = [0.5,1.0]
for timescaling_index in range(5):
    for heritability_index in range(2):
        print(count)
        count += 1
        timescaling = timescale_vec[timescaling_index]
        heritability = heritability_vec[heritability_index]

        data_name = data_dir + 'BWest_t%i_h%i.npy' % (int(timescaling),int(heritability_index))
        est_data = np.load(data_name).item()
        generation = len(est_data['gamma'])
        population = len(est_data['gamma'][0])

        last_fitness = est_data['fitness'][generation-1]
        q5 = np.argsort(last_fitness)[-population // 20]  # best 5%
        fit_index = np.where(last_fitness > last_fitness[q5])[0]

        # mean of all samples
        # gamma_mean = np.mean(est_data['gamma'][generation-1])
        # a_mean = np.mean(est_data['a'][generation-1])
        # nv_mean = np.mean(est_data['nu'][generation-1])
        # mean of the top 5% samples
        gamma_mean = np.mean(est_data['gamma'][generation - 1][fit_index])
        a_mean = np.mean(est_data['a'][generation - 1][fit_index])
        nv_mean = np.mean(est_data['nu'][generation - 1][fit_index])
        vm_mean = np.mean(est_data['vm'][generation - 1][fit_index])
        theta_mean = np.mean(est_data['theta'][generation - 1][fit_index])
        gamma_list.append(gamma_mean)
        a_list.append(a_mean)
        nv_list.append(nv_mean)
        vm_list.append(vm_mean)
        theta_list.append(theta_mean)

        print('Count: %i; gamma:%.3e; alpha: %.3e; nv: %.3e; Vm: %f theta: %f ...' % (count,gamma_mean,a_mean,nv_mean,
                                                                            vm_mean,theta_mean))
        # random test

        length = 10 ** logTL
        obsZ = sorted(length)
        meantrait = np.mean(obsZ)
        td = DVTreeData(path=obs_file, scalar=timescaling)

        param = DVParamLiang(gamma=gamma_mean, a=a_mean, K=K,h=np.sqrt(heritability), nu=nv_mean, r=1, theta=theta_mean, V00 = .5, V01=.5, Vmax=vm_mean, inittrait=meantrait, initpop=1e5,
         initpop_sigma = 10.0, break_on_mu=False)
        params = np.tile(param, (particle_size, 1))  # duplicate

        predictsim = dvcpp.DVSimTVP(td, params)

        valid = np.where(predictsim['sim_time'] == td.sim_evo_time)[0]

        Z = predictsim['Z'][valid]
        Z = np.nan_to_num(Z)
        Z_df = pd.DataFrame(Z)
        savefilename = data_dir+'predictsim%i_vm.csv' % count
        Z_df.to_csv(savefilename,sep=',',index=False)


gamma_array = np.array(gamma_list)
a_array = np.array(a_list)
ratio = a_array/gamma_array