import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
import seaborn as sns


def argsort2D(X):
    return np.arange(len(X))[:, np.newaxis], np.argsort(X, axis=1)


dir_path = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/'

obs_file = dir_path + 'treedata/'
simdata_file = dir_path + 'result_cluster/results_0729_contrast_trait/'
sim_file_TVP = simdata_file + 'predictsimTVP.csv'
sim_file_TV = simdata_file + 'predictsimTV.csv'
sim_file_TVM = simdata_file + 'predictsimTVM.csv'

simN_file_TVP = simdata_file + 'predictsimNTVP.csv'
simN_file_TVM = simdata_file + 'predictsimNTVM.csv'

simV_file_TVP = simdata_file + 'predictsimVTVP.csv'
simV_file_TV = simdata_file + 'predictsimVTV.csv'
simV_file_TVM = simdata_file + 'predictsimVTVM.csv'

with open(obs_file + 'extantspecieslabels.csv') as csv_file:
    csv1_reader = csv.reader(csv_file, delimiter=',')
    extantlabels = list(csv1_reader)

with open(dir_path + 'slater_length_data.csv') as csv_file:
    csv2_reader = csv.reader(csv_file, delimiter=',')
    lengthdata = list(csv2_reader)

extantlabels_array = np.array(extantlabels)[1:, 1]
lengthdata_array = np.array(lengthdata)
length_index = []
for label in extantlabels_array:
    length_index.append(np.where(lengthdata_array[:, 0] == label)[0][0])

logTL = lengthdata_array[length_index, 1].astype(np.float)
length = 10 ** logTL
obsZ = sorted(length)

# load trait data
with open(sim_file_TVP) as csv_file:
    csv2_reader = csv.reader(csv_file, delimiter=',')
    trait_TVP = list(csv2_reader)[1:]

with open(sim_file_TV) as csv_file:
    csv2_reader = csv.reader(csv_file, delimiter=',')
    trait_TV = list(csv2_reader)[1:]

with open(sim_file_TVM) as csv_file:
    csv2_reader = csv.reader(csv_file, delimiter=',')
    trait_TVM = list(csv2_reader)[1:]

# load population data
with open(simN_file_TVP) as csv_file:
    csv2_reader = csv.reader(csv_file, delimiter=',')
    pop_TVP = list(csv2_reader)[1:]

with open(simN_file_TVM) as csv_file:
    csv2_reader = csv.reader(csv_file, delimiter=',')
    pop_TVM = list(csv2_reader)[1:]

# load variance data
with open(simV_file_TVP) as csv_file:
    csv2_reader = csv.reader(csv_file, delimiter=',')
    var_TVP = list(csv2_reader)[1:]

with open(simV_file_TV) as csv_file:
    csv2_reader = csv.reader(csv_file, delimiter=',')
    var_TV = list(csv2_reader)[1:]

with open(simV_file_TVM) as csv_file:
    csv2_reader = csv.reader(csv_file, delimiter=',')
    var_TVM = list(csv2_reader)[1:]

trait_TVP_array = np.array(trait_TVP, dtype=float)
pop_TVP_array = np.array(pop_TVP, dtype=float)
var_TVP_array = np.array(var_TVP, dtype=float)
i, j = argsort2D(trait_TVP_array)
trait_TVP_array = trait_TVP_array[i, j]
pop_TVP_array = pop_TVP_array[i, j]
var_TVP_array = var_TVP_array[i, j]

trait_TV_array = np.array(trait_TV, dtype=float)
var_TV_array = np.array(var_TV, dtype=float)
i, j = argsort2D(trait_TV_array)
trait_TV_array = trait_TV_array[i, j]
var_TV_array = var_TV_array[i, j]

trait_TVM_array = np.array(trait_TVM, dtype=float)
pop_TVM_array = np.array(pop_TVM, dtype=float)
var_TVM_array = np.array(var_TVM, dtype=float)
i, j = argsort2D(trait_TVM_array)
trait_TVM_array = trait_TVM_array[i, j]
pop_TVM_array = pop_TVM_array[i, j]
var_TVM_array = var_TVM_array[i, j]

obsZ_distance = np.linalg.norm(obsZ)

tvp_distance = np.linalg.norm(trait_TVP_array - obsZ, axis=1) / obsZ_distance
tv_distance = np.linalg.norm(trait_TV_array - obsZ, axis=1) / obsZ_distance
tvm_distance = np.linalg.norm(trait_TVM_array - obsZ, axis=1) / obsZ_distance

# Make a separate list for each airline
x1 = list(tvp_distance)
x2 = list(tv_distance)
x3 = list(tvm_distance)

# Assign colors for each airline and the names
colors = ['red', 'green', 'blue']
names = ['TVP', 'TV', 'TVM']

# Make the histogram using a list of lists
# Normalize the flights and assign colors and names
fig = plt.figure(figsize=(9, 9))
plt.hist([x1, x2, x3], bins=20, color=colors, label=names)
plt.legend()
plt.xlabel('Relatvie Distance $\\frac{||Z_{sim}-Z_{obs}||}{||Z_{obs}||}$')
plt.ylabel('Frequence')
plt.title('Goodness of fit distribution')


def gap(targetlist):
    return np.array([targetlist[i + 1] - targetlist[i] for i in range(len(targetlist) - 1)])


gap_tvp = np.apply_along_axis(gap, 1, trait_TVP_array)
gap_tv = np.apply_along_axis(gap, 1, trait_TV_array)
gap_tvm = np.apply_along_axis(gap, 1, trait_TVM_array)

gap_obs = gap(obsZ)

gap_tvp_flatten = gap_tvp.flatten()
gap_tv_flatten = gap_tv.flatten()
gap_tvm_flatten = gap_tvm.flatten()

gap_whole = np.concatenate([gap_tvp_flatten, gap_tv_flatten, gap_tvm_flatten])
gap_label = np.tile(range(14), 3000)
model_label = np.repeat(['TVP', 'TV', 'TVM'], 14000)
df_list = {'gap': gap_whole, 'gap_label': gap_label, 'model': model_label}
df_pd = pd.DataFrame(df_list)
species_gap = range(14)

# Gap plot
fig, axes = plt.subplots(1, 3, sharey=True, sharex=True, figsize=(10, 10))
axes = axes.flatten()

ax_tvp = sns.boxplot(x="gap_label", y="gap", data=df_pd[df_pd['model'] == 'TVP'],
                     showfliers=False, ax=axes[0], color="skyblue")
axes[0].scatter(species_gap, gap_obs, c='red')
ax_tvp.set_title('TVP')

ax_tv = sns.boxplot(x="gap_label", y="gap", data=df_pd[df_pd['model'] == 'TV'],
                    showfliers=False, ax=axes[1], color="skyblue")
axes[1].scatter(species_gap, gap_obs, c='red')
ax_tv.set_title('TV')

ax_tvm = sns.boxplot(x="gap_label", y="gap", data=df_pd[df_pd['model'] == 'TVM'],
                     showfliers=False, ax=axes[2], color="skyblue")
axes[2].scatter(species_gap, gap_obs, c='red')
ax_tvm.set_title('TVM')

# Abundance distribution plot

pop_tvp_flatten = pop_TVP_array.flatten()
pop_tvm_flatten = pop_TVM_array.flatten()

pop_whole = np.concatenate([pop_tvp_flatten, pop_tvm_flatten])
pop_label = np.tile(range(15), 2000)
model_label = np.repeat(['TVP', 'TVM'], 15000)
popdf_list = {'Abundance': pop_whole, 'Species': pop_label, 'model': model_label}
popdf_pd = pd.DataFrame(popdf_list)
species_pop = range(15)

fig_pop, axes_pop = plt.subplots(1, 2, sharey=True, sharex=True, figsize=(10, 10))
axes_pop = axes_pop.flatten()

ax_poptvp = sns.boxplot(x="Species", y="Abundance", data=popdf_pd[popdf_pd['model'] == 'TVP'],
                        showfliers=False, ax=axes_pop[0], color="skyblue")
ax_poptvp.set_title('TVP')

ax_poptvm = sns.boxplot(x="Species", y="Abundance", data=popdf_pd[popdf_pd['model'] == 'TVM'],
                        showfliers=False, ax=axes_pop[1], color="skyblue")
ax_poptvm.set_title('TVM')

axes_pop[0].ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

# Variance distribution plot

var_tvp_flatten = var_TVP_array.flatten()
var_tv_flatten = var_TV_array.flatten()
var_tvm_flatten = var_TVM_array.flatten()

var_whole = np.concatenate([var_tvp_flatten, var_tv_flatten, var_tvm_flatten])
var_label = np.tile(range(15), 3000)
model_label = np.repeat(['TVP', 'TV', 'TVM'], 15000)
vardf_list = {'Variance': var_whole, 'Species': var_label, 'model': model_label}
vardf_pd = pd.DataFrame(vardf_list)
species_var = range(15)

fig_var, axes_var = plt.subplots(1, 3, sharey=True, sharex=True, figsize=(10, 10))
axes_var = axes_var.flatten()

ax_vartvp = sns.boxplot(x="Species", y="Variance", data=vardf_pd[vardf_pd['model'] == 'TVP'],
                        showfliers=False, ax=axes_var[0], color="skyblue")
ax_vartvp.set_title('TVP')

ax_vartv = sns.boxplot(x="Species", y="Variance", data=vardf_pd[vardf_pd['model'] == 'TV'],
                       showfliers=False, ax=axes_var[1], color="skyblue")
ax_vartv.set_title('TV')

ax_vartvm = sns.boxplot(x="Species", y="Variance", data=vardf_pd[vardf_pd['model'] == 'TVM'],
                        showfliers=False, ax=axes_var[2], color="skyblue")
ax_vartvm.set_title('TVM')

# axes_var[0].ticklabel_format(style='sci', axis='y', scilimits=(0, 0))


# Traits distribution plot

trait_TVP_array = trait_TVP_array.flatten()
trait_TV_array = trait_TV_array.flatten()
trait_TVM_array = trait_TVM_array.flatten()

trait_whole = np.concatenate([trait_TVP_array, trait_TV_array, trait_TVM_array])
trait_label = np.tile(range(15), 3000)
model_label = np.repeat(['TVP', 'TV', 'TVM'], 15000)
traitdf_list = {'Trait': trait_whole, 'Species': trait_label, 'model': model_label}
traitdf_pd = pd.DataFrame(traitdf_list)
species_trait = range(15)

fig_trait, axes_trait = plt.subplots(1, 3, sharey=True, sharex=True)
axes_trait = axes_trait.flatten()

ax_traittvp = sns.boxplot(x="Species", y="Trait", data=traitdf_pd[traitdf_pd['model'] == 'TVP'],
                          showfliers=False, ax=axes_trait[0], color="skyblue")
ax_traittvp.set_title('TVP')
axes_trait[0].scatter(species_trait, obsZ, c='red')

ax_traittv = sns.boxplot(x="Species", y="Trait", data=traitdf_pd[traitdf_pd['model'] == 'TV'],
                         showfliers=False, ax=axes_trait[1], color="skyblue")
ax_traittv.set_title('TV')
axes_trait[1].scatter(species_trait, obsZ, c='red')

ax_traittvm = sns.boxplot(x="Species", y="Trait", data=traitdf_pd[traitdf_pd['model'] == 'TVM'],
                          showfliers=False, ax=axes_trait[2], color="skyblue")
ax_traittvm.set_title('TVM')
axes_trait[2].scatter(species_trait, obsZ, c='red')

# Trait vs. population size distribution

trait_TVP_array = trait_TVP_array.flatten()
trait_TVM_array = trait_TVM_array.flatten()

pop_tvp_flatten = pop_TVP_array.flatten()
pop_tvm_flatten = pop_TVM_array.flatten()

trait_whole = np.concatenate([trait_TVP_array, trait_TVM_array])
pop_whole = np.concatenate([pop_tvp_flatten, pop_tvm_flatten])
model_label = np.repeat(['TVP', 'TVM'], 15000)
traitpopdf_list = {'Abundance': pop_whole, 'Trait': trait_whole, 'model': model_label}
traitpopdf_pd = pd.DataFrame(traitpopdf_list)

fig_traitpop, axes_traitpop = plt.subplots(1, 2, sharey=True, sharex=True)
axes_traitpop = axes_traitpop.flatten()

ax_traitpoptvp = sns.scatterplot(x="Trait", y="Abundance", s=10,
                                 data=traitpopdf_pd[traitpopdf_pd['model'] == 'TVP'],
                                 ax=axes_traitpop[0], color="red")
ax_traitpoptvp.set_title('TVP')

ax_traitpoptvm = sns.scatterplot(x="Trait", y="Abundance", s=10,
                                 data=traitpopdf_pd[traitpopdf_pd['model'] == 'TVM'],
                                 ax=axes_traitpop[1], color="red")
ax_traitpoptvm.set_title('TVM')
axes_traitpop[0].ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

# Trait vs. metabolism distribution

trait_TVM_array = trait_TVM_array.flatten()

pop_tvm_flatten = pop_TVM_array.flatten()

meta_TVM_array = pop_tvm_flatten * trait_TVM_array ** (9 / 4)

model_label = np.repeat(['TVM'], 15000)
traitpopdf_list = {'Metabolism': meta_TVM_array, 'Trait': trait_TVM_array, 'model': model_label}
traitpopdf_pd = pd.DataFrame(traitpopdf_list)

fig_meta, axes_meta = plt.subplots(1, 1, sharey=True, sharex=True)

axes_meta = sns.scatterplot(x="Trait", y="Metabolism", s=4,
                            data=traitpopdf_pd[traitpopdf_pd['model'] == 'TVM'],
                            ax=axes_meta, color="red")
axes_meta.set_title('TVM')

# test the energetic equivalence rule
trait_TVM_array = trait_TVM_array.flatten()
log_trait = np.log10(trait_TVM_array)  # .reshape((-1, 1))
pop_tvm_flatten = pop_TVM_array.flatten()

from sklearn.linear_model import LinearRegression

model = LinearRegression().fit(log_trait, pop_tvm_flatten)
print('slope:', model.coef_)

model_label = np.repeat(['TVM'], 15000)
traitpopdf_list = {'Abundance': pop_tvm_flatten, 'Trait': log_trait, 'model': model_label}
traitpopdf_pd = pd.DataFrame(traitpopdf_list)

fig_traitpop, axes_traitpop = plt.subplots(1, 1, sharey=True, sharex=True)

axes_traitpop = sns.scatterplot(x="Trait", y="Abundance", s=10,
                                data=traitpopdf_pd[traitpopdf_pd['model'] == 'TVM'],
                                ax=axes_traitpop, color="red")
ax_traitpoptvp.set_title('TVM')
