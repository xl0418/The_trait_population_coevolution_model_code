import os
import sys
import platform
if platform.system()=='Windows':
    sys.path.append('C:/Liang/abcpp_master/abcpp')
elif platform.system()=='Darwin':
    sys.path.append('/Users/dudupig/Documents/GitHub/Code/Pro2/Python_p2')
import numpy as np
from dvtraitsim_shared import DVTreeData, DVParam
import dvtraitsim_cpp as dvcpp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kendalltau
sns.set(style="white")
#
singlesim = False
gamma_vec = np.array([0,0.001,0.01,0.1,0.5,1])
a_vec = gamma_vec

for no_tree in range(1,2):
    print(no_tree)
# no_tree= 5
    dir_path = 'c:/Liang/Googlebox/Research/Project2'
    files = dir_path + '/treesim_newexp/example%d/' % no_tree

    td = DVTreeData(path=files, scalar=20000)
    K = 10e8
    nu=1/(100*K)
    num = 11
    trait_w = []
    trait_v = []
    pop_w = []

    for gamma in gamma_vec:
        for a in a_vec:
            # let's try to find a true simulation:
            obs_param = DVParam(gamma=gamma, a=a, K=K, nu=nu, r=1, theta=0, Vmax=1, inittrait=0, initpop=500,
                                initpop_sigma=10.0, break_on_mu=False)
            trait_data = ()
            population_data = ()
            traitvar_data = ()
            for loop in range(1,num):
                str = 'gamma = %.3f; a = %.3f; loop = %d' % (gamma,a,loop)
                print(str)
                par_obs = np.array([gamma, a])
                simresult = dvcpp.DVSim(td, obs_param)
                if simresult['sim_time'] == td.sim_evo_time:
                    trait_tips = simresult['Z']
                    population_tips = simresult['N']
                    traitvar_tips = simresult['V']
                    # empirical data for trait and population
                    index_extant = np.where(~np.isnan(trait_tips))[0]
                    trait_tips = trait_tips[index_extant]
                    population_tips = population_tips[index_extant]
                    traitvar_tips = traitvar_tips[index_extant]
                    assert len(trait_tips) == len(population_tips)
                    assert len(trait_tips) == len(traitvar_tips)
                    trait_data = np.append(trait_data, trait_tips)
                    population_data = np.append(population_data, population_tips)
                    traitvar_data = np.append(traitvar_data, traitvar_tips)
                    if singlesim: break
                if simresult['sim_time'] < td.sim_evo_time:
                    print('Jump to the next loop')

            trait_w.append(trait_data)
            trait_v.append(traitvar_data)
            pop_w.append(population_data)

    num_tips = len(trait_tips)

    normed_trait = []
    normed_traitvar = []
    normed_pop = []

    for i in range(0,36):
        if len(trait_w[i])==0:
            normed_trait.append([0])
            normed_traitvar.append([0])
            normed_pop.append([0])
        else:
            normed_trait.append((trait_w[i]- np.min(trait_w[i])) / (np.max(trait_w[i]) - np.min(trait_w[i])))
            normed_traitvar.append((trait_v[i]- np.min(trait_v[i])) / (np.max(trait_v[i]) - np.min(trait_v[i])))
            normed_pop.append((pop_w[i] - np.min(pop_w[i])) / (np.max(pop_w[i]) - np.min(pop_w[i])))


    count = 0
    label_a = (['a=0','a=.001','a=.01','a=.1','a=.5','a=1'])
    label_gamma = (['$\gamma$=0','$\gamma$=.001','$\gamma$=.01','$\gamma$=.1','$\gamma$=.5','$\gamma$=1'])

    # Set up the matplotlib figure
    f1, axes1 = plt.subplots(6, 6, figsize=(9, 9), sharex=True, sharey=True) #, sharex=True, sharey=True
    f2, axes2 = plt.subplots(6, 6, figsize=(9, 9), sharex=True, sharey=True)
    f3, axes3 = plt.subplots(6, 6, figsize=(9, 9), sharex=True, sharey=True)
    # Rotate the starting point around the cubehelix hue circle
    # for ax,ax2,ax3, s in zip(axes.flat,axes2.flat,axes3.flat, np.linspace(0, 3,36)):
        # Create a cubehelix colormap to use with kdeplot
    for index_g in range(len(gamma_vec)):
        gamma1=gamma_vec[index_g]
        for index_a in range(len(a_vec)):
            a=a_vec[index_a]
            print(count)
            if len(trait_w[count]) == 1:
                axes1[index_g,index_a].plot()
                axes2[index_g,index_a].plot()
                axes3[index_g,index_a].plot()

            else:
                trait = trait_w[count]
                traitvar = trait_v[count]
                pop = pop_w[count]
                #
                # trait = normed_trait[count]
                # traitvar = normed_traitvar[count]
                # pop = normed_pop[count]
                # axes1[index_g, index_a].set_xlim([-0.2,1.2])
                # axes2[index_g, index_a].set_xlim([-0.2,1.2])
                # axes3[index_g, index_a].set_xlim([-0.2,1.2])

                # Generate and plot a random bivariate dataset
                # sns.kdeplot(trait, pop, cmap=cmap, shade=True, cut=5, ax=ax)
                sns.scatterplot(trait, pop,ax=axes1[index_g, index_a])
                sns.scatterplot(trait, traitvar,ax=axes2[index_g, index_a])
                sns.scatterplot(pop, traitvar,ax=axes3[index_g, index_a])

            if count in range(0,6):
                axes1[index_g, index_a].title.set_text(label_a[count])
                axes2[index_g, index_a].title.set_text(label_a[count])
                axes3[index_g, index_a].title.set_text(label_a[count])

            if count in ([5, 11, 17, 23, 29, 35]):
                axes1[index_g, index_a].set_ylabel(label_gamma[int(count/6)])
                axes1[index_g, index_a].yaxis.set_label_position("right")
                axes2[index_g, index_a].set_ylabel(label_gamma[int(count/6)])
                axes2[index_g, index_a].yaxis.set_label_position("right")
                axes3[index_g, index_a].set_ylabel(label_gamma[int(count/6)])
                axes3[index_g, index_a].yaxis.set_label_position("right")
            count += 1
    f1.text(0.5, 0.001, 'Trait mean', ha='center', fontsize=15)
    f1.text(0.001, 0.5, 'Population size', va='center', rotation='vertical', fontsize=15)
    f1.tight_layout()
    f2.text(0.5, 0.001, 'Trait mean', ha='center', fontsize=15)
    f2.text(0.001, 0.5, 'Trait variance', va='center', rotation='vertical', fontsize=15)
    f2.tight_layout()
    f3.text(0.5, 0.001, 'Trait variance', ha='center', fontsize=15)
    f3.text(0.001, 0.5, 'Population size', va='center', rotation='vertical', fontsize=15)
    f3.tight_layout()


    tree = 'tree'+'%d' % no_tree
    dir_fig = 'C:/Liang/Googlebox/Research/Project2/smc_new100replicatesdistribution/'+tree

    f1.savefig(dir_fig+'+meanpopscatter100rep.png')
    plt.close(f1)
    f2.savefig(dir_fig+'+meanvarscatter100rep.png')
    plt.close(f2)
    f3.savefig(dir_fig+'+varpopscatter100rep.png')
    plt.close(f3)
    # f.savefig('C:/Liang/Googlebox/Research/Project2/smc_new100replicatesdistributiontraitvsvar.png')
