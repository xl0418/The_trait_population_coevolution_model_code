import os
import numpy as np
import platform
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(style="white")

gamma_list = []
a_list = []
nv_list = []
fit_list = []

count = 0
timescale_vec = [20000,40000,80000]
heritability_vec = [1,2]
# dividing_vec = [1,4]
data_dir = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/result_cluster/'
for timescaling_index in range(3):
    # for dividing_index in range(2):
        for heritability_index in range(2):
            count += 1
            timescaling = timescale_vec[timescaling_index]
            heritability = heritability_vec[heritability_index]
            # dividing = dividing_vec[dividing_index]
            #
            # data_name = data_dir + 'BWest_t%i_d%i_h%i.npy' % (int(timescaling),int(dividing),int(heritability_index))
            # est_data = np.load(data_name).item()
            # generation = len(est_data['gamma'])
            # population = len(est_data['gamma'][0])
            # gamma_list.append(est_data['gamma'][generation-1])
            # a_list.append(est_data['a'][generation-1])
            # nv_list.append(est_data['nu'][generation-1])

            # random test
            gamma_list.append(np.random.normal(count,1,20000))
            a_list.append(np.random.normal(count,1,20000))
            nv_list.append(np.random.normal(count,1,20000))


gamma_flat_list = [item for sublist in gamma_list for item in sublist]
a_flat_list = [item for sublist in a_list for item in sublist]
nv_flat_list = [item for sublist in nv_list for item in sublist]


timescaling_list = np.repeat(timescale_vec,2*20000)
# dividing_list = np.tile(np.repeat(dividing_vec,2*20000),3)
heritability_list = np.tile(np.repeat(heritability_vec,20000),3)

est_list = {'gamma':gamma_flat_list,'alpha':a_flat_list,'nu':nv_flat_list,'timescale':timescaling_list,
            'heritability':heritability_list}
est_df = pd.DataFrame(est_list)



f, axes = plt.subplots(3, 2,figsize = (9,12),sharex=True,sharey=True)
ax1 = sns.boxplot(x="timescale", y="gamma",
                 hue="heritability", palette=["m", "g"],
                 data=est_df[lambda est_df: est_df.dividing == 1],
                  linewidth=1, ax=axes[0,0], showfliers=False)
ax1.title.set_text('Dividing scalar: 1')
ax2 = sns.boxplot(x="timescale", y="gamma",
                 hue="heritability", palette=["m", "g"],
                 data=est_df[lambda est_df: est_df.dividing == 4],
                  linewidth=1, ax=axes[0,1], showfliers=False)

ax2.title.set_text('Dividing scalar: 4')

ax3 =sns.boxplot(x="timescale", y="alpha",
                 hue="heritability", palette=["m", "g"],
                 data=est_df[lambda est_df: est_df.dividing == 1],
                  linewidth=1, ax=axes[1,0], showfliers=False)


ax4 = sns.boxplot(x="timescale", y="alpha",
                 hue="heritability", palette=["m", "g"],
                 data=est_df[lambda est_df: est_df.dividing == 4],
                  linewidth=1, ax=axes[1,1], showfliers=False)
ax5 = sns.boxplot(x="timescale", y="nu",
                 hue="heritability", palette=["m", "g"],
                 data=est_df[lambda est_df: est_df.dividing == 1],
                  linewidth=1, ax=axes[2,0], showfliers=False)


ax6 = sns.boxplot(x="timescale", y="nu",
                 hue="heritability", palette=["m", "g"],
                 data=est_df[lambda est_df: est_df.dividing == 4],
                  linewidth=1, ax=axes[2,1], showfliers=False)

handles_gamma, labels_gamma = ax1.get_legend_handles_labels()
handles_a, labels_a = ax3.get_legend_handles_labels()
handles_nu, labels_nu = ax5.get_legend_handles_labels()
for ax in [ax1,ax2,ax3,ax4,ax5,ax6]:
    ax.set_ylabel('')
    ax.set_xlabel('')
    ax.legend_.remove()

axes[0, 0].set_ylabel('$\gamma$')
axes[0, 0].yaxis.set_label_position("left")
axes[1, 0].set_ylabel('$\\alpha$')
axes[1, 0].yaxis.set_label_position("left")
axes[2, 0].set_ylabel('$\\nu$')
axes[2, 0].yaxis.set_label_position("left")

handles = handles_gamma #[ item for subhandle in [handles_gamma,handles_a,handles_nu] for item in subhandle]
labels = ['$h^2 = 1$', '$h^2 = 0.5$']
f.text(0.5, 0.04, 'Time scaling parameters', ha='center',fontsize=15)
f.text(0.04, 0.5, 'Estimates', va='center', rotation='vertical',fontsize=15)
l = plt.legend(handles, labels, bbox_to_anchor=(0.85, 3.7), loc=2, borderaxespad=0.)
l.get_frame().set_linewidth(0.0)


f.savefig(dir_fig)
plt.close(f)



