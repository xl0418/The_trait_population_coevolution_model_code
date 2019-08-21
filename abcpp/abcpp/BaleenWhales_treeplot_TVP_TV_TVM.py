import sys, os
import platform
sys.path.append('C:/Liang/abcpp_ms5/abcpp')
from Dvtraitsim_TVP import DVSimTVP
from dvtraitsim_shared import DVTreeData, DVParamLiang
from Dvtraitsim_TV import DVSimTV
from Dvtraitsim_TVM import DVSimTVM

import numpy as np
import matplotlib.pyplot as plt



theta = 0  # optimum of natural selection
r = 1  # growth rate
Vmax = 100
scalar = 2000
K=10e5
K_TVM = 1e12
nu=1e-4
timegap = 1

# let's try to find a true simulation:


#full tree and pruned tree directory
dir_path = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/'

files = dir_path + 'treedata/'

meantrait = 1300

# trait evolution plot

gamma_vec = np.array([0, 1e-10,1e-9,1e-8,1e-7,1e-6])
a_vec = np.array([1e-5, 5e-5,1e-4,5e-4,1e-3,5e-3])

gamma_vec_TV = np.array([0, 1e-10,1e-9,1e-8,1e-7,1e-6])
a_vec_TV = np.array([1e-3, 5e-3,1e-2,5e-2,1e-1,5e-1])
row_gamma = len(gamma_vec)
count = 0

td = DVTreeData(path=files, scalar=scalar)


f1, axes1 = plt.subplots(row_gamma, row_gamma, figsize=(9, 9),sharey=True,sharex=True) #
f2, axes2 = plt.subplots(row_gamma, row_gamma, figsize=(9, 9),sharey=True,sharex=True) #
f3, axes3 = plt.subplots(row_gamma, row_gamma, figsize=(9, 9),sharey=True,sharex=True) #


label_gamma = (['$\gamma=0$','$\gamma=1e-10$','$\gamma=1e-9$','$\gamma=1e-8$','$\gamma=1e-7$','$\gamma=1e-6$'])
label_a = (['$\\alpha=1e-5$','$\\alpha=5e-5$','$\\alpha=1e-4$','$\\alpha=5e-4$','$\\alpha=1e-3$','$\\alpha=5e-3$'])

for index_g in range(len(gamma_vec)):
    gamma1=gamma_vec[index_g]
    for index_a in range(len(a_vec)):
        a=a_vec[index_a]
        print( count)
        for replicate in range(100):
            print('Simulating TVP model...')
            obs_param_TVP = DVParamLiang(gamma=gamma1, a=a, K=K,h=1, nu=nu, r=r, theta=meantrait,V00=.1,V01=.1, Vmax=Vmax, inittrait=meantrait, initpop=1e5,
                                initpop_sigma=10.0, break_on_mu=False)
            simresult_TVP = DVSimTVP(td,obs_param_TVP)
            if simresult_TVP['sim_time'] == td.sim_evo_time:
                pic_TVP = 0
                break
            else:
                pic_TVP=1

        for replicate in range(100):
            print('Simulating TV model...')

            obs_param_TV = DVParamLiang(gamma=gamma1, a=a, K=K,h=1, nu=nu, r=r, theta=meantrait,V00=.1,V01=.1, Vmax=100, inittrait=meantrait, initpop=1e5,
                                initpop_sigma=10.0, break_on_mu=False)
            simresult_TV = DVSimTV(td,obs_param_TV)
            if simresult_TV['sim_time'] == td.sim_evo_time:
                pic_TV = 0
                break
            else:
                pic_TV=1

        for replicate in range(100):
            print('Simulating TVM model...')

            obs_param_TVM = DVParamLiang(gamma=gamma1, a=a, K=K_TVM,h=1, nu=nu, r=r, theta=meantrait,V00=.1,V01=.1, Vmax=Vmax, inittrait=meantrait, initpop=1e5,
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
        trait_TVM = simresult_TVM['Z']

        num_lines = total_species
        x = np.arange(evo_time/timegap+1)

        labels = []
        for i in range(1, num_lines + 1):
            axes1[index_g,index_a].plot(x, trait_TVP[::timegap, i - 1])
            axes2[index_g,index_a].plot(x, trait_TV[::timegap, i - 1])
            axes3[index_g,index_a].plot(x, trait_TVM[::timegap, i - 1])


        axes1[index_g, index_a].xaxis.set_major_locator(plt.NullLocator())
        axes2[index_g, index_a].xaxis.set_major_locator(plt.NullLocator())
        axes3[index_g, index_a].xaxis.set_major_locator(plt.NullLocator())


        # axes[index_g, index_a].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
        # axes[index_g, index_a].set_yscale('log')
        axes1[index_g, index_a].ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        axes2[index_g, index_a].ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        axes3[index_g, index_a].ticklabel_format(style='sci', axis='y', scilimits=(0, 0))


        if count in range(0, row_gamma):
            axes1[index_g, index_a].title.set_text(label_a[count])
            axes2[index_g, index_a].title.set_text(label_a[count])
            axes3[index_g, index_a].title.set_text(label_a[count])


        if count in ([5, 11, 17, 23, 29, 35]):
            axes1[index_g, index_a].set_ylabel(label_gamma[int(count / row_gamma)])
            axes1[index_g, index_a].yaxis.set_label_position("right")
            axes2[index_g, index_a].set_ylabel(label_gamma[int(count / row_gamma)])
            axes2[index_g, index_a].yaxis.set_label_position("right")
            axes3[index_g, index_a].set_ylabel(label_gamma[int(count / row_gamma)])
            axes3[index_g, index_a].yaxis.set_label_position("right")

        count += 1


dir_fig = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/traittree'

#
# f1.savefig(dir_fig+'TP2w.png')
# plt.close(f1)

