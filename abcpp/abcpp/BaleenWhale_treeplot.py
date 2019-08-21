import sys, os
import platform
if platform.system()=='Windows':
    sys.path.append('C:/Liang/abcpp_ms5/abcpp')
elif platform.system()=='Darwin':
    sys.path.append('/Users/dudupig/Documents/GitHub/Code/Pro2/Python_p2')
from dvtraitsim_py import DVSimLiang
from dvtraitsim_shared import DVTreeData, DVParamLiang
import numpy as np
import matplotlib.pyplot as plt

theta = 0  # optimum of natural selection
r = 1  # growth rate
Vmax = 1
scalar = 2000
K=10e5
nu=1e-4
timegap = 20

# let's try to find a true simulation:


#full tree and pruned tree directory
dir_path = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/'

files = dir_path + 'treedata/'



# trait evolution plot

gamma_vec = np.array([0, 1e-10,1e-9,1e-8,1e-7,1e-6])
a_vec = np.array([0, 1e-4,1e-3,1e-2,1e-1,1])
row_gamma = len(gamma_vec)
count = 0

td = DVTreeData(path=files, scalar=scalar)


f1, axes1 = plt.subplots(row_gamma, row_gamma, figsize=(9, 9),sharey=True,sharex=True) #


label_gamma = (['$\gamma=0$','$\gamma=1e-10$','$\gamma=1e-9$','$\gamma=1e-8$','$\gamma=1e-7$','$\gamma=1e-6$'])
label_a = (['$\\alpha=0$','$\\alpha=1e-4$','$\\alpha=1e-3$','$\\alpha=1e-2$','$\\alpha=1e-1$','$\\alpha=1$'])

for index_g in range(len(gamma_vec)):
    gamma1=gamma_vec[index_g]
    for index_a in range(len(a_vec)):
        a=a_vec[index_a]
        print( count)
        for replicate in range(100):
            obs_param = DVParamLiang(gamma=gamma1, a=a, K=K,h=1, nu=nu, r=r, theta=theta,V00=.1,V01=.1, Vmax=100, inittrait=0, initpop=500,
                                initpop_sigma=10.0, break_on_mu=False)
            simresult = DVSimLiang(td,obs_param)
            if simresult['sim_time'] == td.sim_evo_time:
                pic = 0
                break
            else:
                pic=1
        # if pic==0:
        evo_time, total_species = simresult['N'].shape
        evo_time = evo_time - 1
        trait_RI_dr = simresult['Z']
        population_RI_dr = simresult['N']
        population_RI_dr = population_RI_dr.astype(float)
        population_RI_dr[population_RI_dr==0] = np.nan
        V_dr = simresult['V']
        num_lines = total_species
        x = np.arange(evo_time/timegap)

        labels = []
        for i in range(1, num_lines + 1):
            axes1[index_g,index_a].plot(x, trait_RI_dr[::timegap, i - 1])


        # axes[index_g, index_a].yaxis.set_major_locator(plt.NullLocator())
        axes1[index_g, index_a].xaxis.set_major_locator(plt.NullLocator())


        # axes[index_g, index_a].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
        # axes[index_g, index_a].set_yscale('log')
        axes1[index_g, index_a].ticklabel_format(style='sci', axis='y', scilimits=(0, 0))


        if count in range(0, row_gamma):
            axes1[index_g, index_a].title.set_text(label_a[count])


        if count in ([5, 11, 17, 23, 29, 35]):
            axes1[index_g, index_a].set_ylabel(label_gamma[int(count / row_gamma)])
            axes1[index_g, index_a].yaxis.set_label_position("right")

        count += 1


dir_fig = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/traittree'


f1.savefig(dir_fig+'TP2w.png')
plt.close(f1)

