import sys
import numpy as np
from dvtraitsim_shared import DVTreeData, DVParamLiang
import dvtraitsim_cpp as dvcpp
import csv
from tp_update import tp_update


#
# argsort of 2-dimensional arrays is awkward
# returns i,j so that X[i,j] == np.sort(X)
def argsort2D(X):
    return np.arange(len(X))[:, np.newaxis], np.argsort(X, axis=1)


def normalized_norm(x, y):
    diff_norm = np.linalg.norm(x - y, axis=1)
    max_err = np.nanmax(diff_norm)
    return diff_norm / max_err


K_TVP=1e6
K_TV = 1e6
K_TVM = 1e12
nu=1e-4


#full tree
dir_path = '/home/p274981/abcpp/'

files = dir_path + 'BaleenWhales/treedata/'

savedir = dir_path + 'BaleenWhales/modelselec2w.npy'

td = DVTreeData(path=files, scalar=20000)


with open(files+'extantspecieslabels.csv') as csv_file:
    csv1_reader = csv.reader(csv_file, delimiter=',')
    extantlabels = list(csv1_reader)

with open(files+'slater_length_data.csv') as csv_file:
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
print('trying to estimate the parameters','...')

s = np.argsort(obsZ)
obsZ = obsZ[s]
obsZ = obsZ.astype(np.float)
meantrait = np.mean(obsZ)
# let's try to find a true simulation:
sampleparam_TVP = DVParamLiang(gamma=1, a=1, K=K_TVP,h=1, nu=nu, r=1, theta=meantrait,V00=.1,V01=.1, Vmax=100, inittrait=meantrait, initpop=1e5,
                                initpop_sigma=10.0, break_on_mu=False)
sampleparam_TV = DVParamLiang(gamma=1, a=1, K=K_TV,h=1, nu=nu, r=1, theta=meantrait,V00=.1,V01=.1, Vmax=100, inittrait=meantrait, initpop=1e5,
                                initpop_sigma=10.0, break_on_mu=False)
sampleparam_TVM = DVParamLiang(gamma=1, a=1, K=K_TVM,h=1, nu=nu, r=1, theta=meantrait,V00=.1,V01=.1, Vmax=100.0, inittrait=meantrait, initpop=1e5,
                                initpop_sigma=10.0, break_on_mu=False)

# pop = dvcpp.DVSim(td, obs_param)

population = 10000
generations = 5
total_population = population*3

prior = [0.0,1e-4,0.0,1e-2,0.0,1e-2,20.0,500.0]


params_TVP = np.tile(sampleparam_TVP, (population, 1))  # duplicate
params_TVP[:, 0] = np.random.uniform(prior[0], prior[1], params_TVP.shape[0])  # randomize 'gamma'
params_TVP[:, 1] = np.random.uniform(prior[2], prior[3], params_TVP.shape[0])  # randomize 'a'
params_TVP[:, 4] = np.random.uniform(prior[4], prior[5], params_TVP.shape[0])  # randomize 'nu'
params_TVP[:, 9] = np.random.uniform(prior[6], prior[7], params_TVP.shape[0])  # randomize 'Vm'

params_TV = np.tile(sampleparam_TV, (population, 1))  # duplicate
params_TV[:, 0] = np.random.uniform(prior[0], prior[1], params_TV.shape[0])  # randomize 'gamma'
params_TV[:, 1] = np.random.uniform(prior[2], prior[3], params_TV.shape[0])  # randomize 'a'
params_TV[:, 4] = np.random.uniform(prior[4], prior[5], params_TV.shape[0])  # randomize 'nu'
params_TV[:, 9] = np.random.uniform(prior[6], prior[7], params_TV.shape[0])  # randomize 'Vm'

params_TVM = np.tile(sampleparam_TVM, (population, 1))  # duplicate
params_TVM[:, 0] = np.random.uniform(prior[0], prior[1], params_TVM.shape[0])  # randomize 'gamma'
params_TVM[:, 1] = np.random.uniform(prior[2], prior[3], params_TVM.shape[0])  # randomize 'a'
params_TVM[:, 4] = np.random.uniform(prior[4], prior[5], params_TVM.shape[0])  # randomize 'nu'
params_TVM[:, 9] = np.random.uniform(prior[6], prior[7], params_TVM.shape[0])  # randomize 'Vm'

# model choice
model_index = np.array([0,1,2])
model_params = np.repeat(model_index,repeats = population)
model_data = np.zeros(shape=(generations+1, total_population))
model_data[0,:] = model_params
propose_model = model_params

# store parameters used
# TVP
gamma_data_TVP = np.zeros(shape=(generations+1, population))
a_data_TVP = np.zeros(shape=(generations+1, population))
nu_data_TVP = np.zeros(shape=(generations+1, population))
vm_data_TVP = np.zeros(shape=(generations+1, population))

# TV
gamma_data_TV = np.zeros(shape=(generations+1, population))
a_data_TV = np.zeros(shape=(generations+1, population))
nu_data_TV = np.zeros(shape=(generations+1, population))
vm_data_TV = np.zeros(shape=(generations+1, population))

# TVM
gamma_data_TVM = np.zeros(shape=(generations+1, population))
a_data_TVM = np.zeros(shape=(generations+1, population))
nu_data_TVM = np.zeros(shape=(generations+1, population))
vm_data_TVM = np.zeros(shape=(generations+1, population))


gamma_data_TVP[0,:] = params_TVP[:,0]
a_data_TVP[0,:] = params_TVP[:,1]
nu_data_TVP[0,:] = params_TVP[:,4]
vm_data_TVP[0,:] = params_TVP[:,9]

gamma_data_TV[0,:] = params_TV[:,0]
a_data_TV[0,:] = params_TV[:,1]
nu_data_TV[0,:] = params_TV[:,4]
vm_data_TV[0,:] = params_TV[:,9]

gamma_data_TV[0,:] = params_TVM[:,0]
a_data_TVM[0,:] = params_TVM[:,1]
nu_data_TVM[0,:] = params_TVM[:,4]
vm_data_TVM[0,:] = params_TVM[:,9]


fitness= np.zeros(shape=(generations, total_population))
# Initialize the weights.
weight_model = np.zeros(total_population)
weight_model.fill(1 / total_population)

# weights for paras of TVP
weight_gamma_TVP = np.zeros(population)
weight_gamma_TVP.fill(1 / population)
weight_a_TVP = np.zeros(population)
weight_a_TVP.fill(1 / population)
weight_nu_TVP = np.zeros(population)
weight_nu_TVP.fill(1 / population)
weight_vm_TVP = np.zeros(population)
weight_vm_TVP.fill(1 / population)

# weights for paras of TV
weight_gamma_TV = np.zeros(population)
weight_gamma_TV.fill(1 / population)
weight_a_TV = np.zeros(population)
weight_a_TV.fill(1 / population)
weight_nu_TV = np.zeros(population)
weight_nu_TV.fill(1 / population)
weight_vm_TV = np.zeros(population)
weight_vm_TV.fill(1 / population)

# weights for paras of TVM
weight_gamma_TVM = np.zeros(population)
weight_gamma_TVM.fill(1 / population)
weight_a_TVM = np.zeros(population)
weight_a_TVM.fill(1 / population)
weight_nu_TVM = np.zeros(population)
weight_nu_TVM.fill(1 / population)
weight_vm_TVM = np.zeros(population)
weight_vm_TVM.fill(1 / population)

for g in range(generations):
    # model 0
    TVP_sample_length = len(np.where(model_data[g,:]==0)[0])
    TV_sample_length = len(np.where(model_data[g,:]==1)[0])
    TVM_sample_length = len(np.where(model_data[g,:]==2)[0])
    Z = np.zeros((1, td.total_species))

    if TVP_sample_length>0:
        print('TVP simulations start...')

        simmodelTVP = dvcpp.DVSimTVP(td, params_TVP)
        valid_TVP = np.where(simmodelTVP['sim_time'] == td.sim_evo_time)[0]
        Z_modelTVP = simmodelTVP['Z'][valid_TVP]
        i, j = argsort2D(Z_modelTVP)
        Z_modelTVP = Z_modelTVP[i, j]
        # V = pop['V'][valid][i, j]
        Z_modelTVP = np.nan_to_num(Z_modelTVP)
        Z = np.vstack([Z,Z_modelTVP])
            # V = np.nan_to_num(V)
            #GOF: Goodness of fit
        if  len(valid_TVP) == 0:
            print('No complete results from TVP model ')


    if TV_sample_length>0:
        print('TV simulations start...')
        # model 1
        # for param_drury in params_DR:
        simmodelTV = dvcpp.DVSimTV(td, params_TV)
        valid_TV = np.where(simmodelTV['sim_time'] == td.sim_evo_time)[0]
        Z_modelTV = simmodelTV['Z'][valid_TV]
        i, j = argsort2D(Z_modelTV)
        Z_modelTV = Z_modelTV[i, j]
        # V = pop['V'][valid][i, j]
        Z_modelTV = np.nan_to_num(Z_modelTV)
        Z = np.vstack([Z, Z_modelTV])


    if TVM_sample_length>0:
        print('TVM simulations start...')
        simmodelTVM = dvcpp.DVSimTVM(td, params_TVM)
        valid_TVM = np.where(simmodelTVM['sim_time'] == td.sim_evo_time)[0]
        Z_modelTVM = simmodelTV['Z'][valid_TVM]
        i, j = argsort2D(Z_modelTVM)
        Z_modelTVM = Z_modelTVM[i, j]
        # V = pop['V'][valid][i, j]
        Z_modelTVM = np.nan_to_num(Z_modelTVM)
        Z = np.vstack([Z, Z_modelTVM])

    Z = Z[1:,]

    valid = np.concatenate([valid_TVP,np.array(valid_TV)+TVP_sample_length,
                            np.array(valid_TVM)+TVP_sample_length+TV_sample_length
                            ]).astype(int)

    eudis = normalized_norm(Z, obsZ)
    # eudis = eudistance(Z, obsZ)


    fitness[g,valid] += 1 - eudis


    # fitness[g,valid] += 1.0 - normalized_norm(np.sqrt(V), np.sqrt(obsV))

    # print something...
    q5 = np.argsort(fitness[g,:])[-int(total_population// 4)]  # best 25%
    fit_index = np.where(fitness[g,:] > fitness[g,q5])[0]

    modelTVPperc = len(np.where(propose_model[fit_index]==0)[0])/len(fit_index)
    modelTVperc = len(np.where(propose_model[fit_index]==1)[0])/len(fit_index)
    modelTVMperc = len(np.where(propose_model[fit_index]==2)[0])/len(fit_index)

    print('Iteration = %d 25th Model TVP: %.1f%% ;  Model TV: %.1f%% ; Model TVM: %.1f%%...'
          % (g,modelTVPperc*100 ,modelTVperc*100,modelTVMperc*100))
    print('Average fitness: %f' % np.mean(fitness[g,fit_index]))
    # reevaluate the weight of the best fitted  models
    weight_model_bestfitted = weight_model[fit_index]*fitness[g,fit_index]/sum(weight_model[fit_index]*fitness[g,fit_index])
    # sample new models from the fitness of previous best fitted models
    sample_model_index = sorted(np.random.choice(fit_index, total_population, p=weight_model_bestfitted))


    propose_model = model_params

    q5_TVP = np.argsort(fitness[g, :population-1])[-int(population // 200)]  # best 5%
    q5_TV = np.argsort(fitness[g, population:2*population-1])[-int(population // 200)]+population  # best 5%
    q5_TVM = np.argsort(fitness[g, 2*population:3*population-1])[-int(population // 200)]+2*population  # best 5%

    fit_index_TVP = np.where(fitness[g, :population-1] > fitness[g, q5_TVP])[0]
    fit_index_TV = np.where(fitness[g, population:2*population-1] > fitness[g, q5_TV])[0]+population
    fit_index_TVM = np.where(fitness[g, 2*population:] > fitness[g, q5_TVM])[0]+2*population

    previous_bestfitted_index_TVP = fit_index_TVP
    previous_bestfitted_index_TV = fit_index_TV - population
    previous_bestfitted_index_TVM = fit_index_TVM- 2*population

    chosengamma_TVP,chosena_TVP,chosennu_TVP,chosenvm_TVP = np.mean(params_TVP[previous_bestfitted_index_TVP,0]),\
                                            np.mean(params_TVP[previous_bestfitted_index_TVP,1]),\
                                                np.mean(params_TVP[previous_bestfitted_index_TVP, 4]),\
                                                np.mean(params_TVP[previous_bestfitted_index_TVP, 9])

    chosengamma_TV,chosena_TV,chosennu_TV,chosenvm_TV = np.mean(params_TV[previous_bestfitted_index_TV,0]),\
                                            np.mean(params_TV[previous_bestfitted_index_TV,1]),\
                                                np.mean(params_TV[previous_bestfitted_index_TV, 4]),\
                                                np.mean(params_TV[previous_bestfitted_index_TV, 9])

    chosengamma_TVM,chosena_TVM,chosennu_TVM,chosenvm_TVM = np.mean(params_TVM[previous_bestfitted_index_TVM,0]),\
                                            np.mean(params_TVM[previous_bestfitted_index_TVM,1]),\
                                                np.mean(params_TVM[previous_bestfitted_index_TVM, 4]),\
                                                np.mean(params_TVM[previous_bestfitted_index_TVM, 9])


    print('Mean estimates: TVP gamma: %.3e ; a: %.3e ; nu: %.3e ; Vm : %f'  % ( chosengamma_TVP,chosena_TVP,chosennu_TVP,chosenvm_TVP))
    print('Mean estimates: TV gamma: %.3e ; a: %.3e ; nu: %.3e ; Vm : %f'  % ( chosengamma_TV,chosena_TV,chosennu_TV,chosenvm_TV))
    print('Mean estimates: TVM gamma: %.3e ; a: %.3e ; nu: %.3e ; Vm : %f'  % ( chosengamma_TVM,chosena_TVM,chosennu_TVM,chosenvm_TVM))




    model_data[g+1,:] = propose_model
    gamma_data_TVP[g+1, :] = params_TVP[:, 0]
    a_data_TVP[g+1, :] = params_TVP[:, 1]
    nu_data_TVP[g+1, :] = params_TVP[:, 4]
    vm_data_TVP[g+1, :] = params_TVP[:, 9]

    gamma_data_TV[g + 1, :] = params_TV[:, 0]
    a_data_TV[g + 1, :] = params_TV[:, 1]
    nu_data_TV[g + 1, :] = params_TV[:, 4]
    vm_data_TV[g + 1, :] = params_TV[:, 9]

    gamma_data_TVM[g+1, :] = params_TVM[:, 0]
    a_data_TVM[g+1, :] = params_TVM[:, 1]
    nu_data_TVM[g+1, :] = params_TVM[:, 4]
    vm_data_TVM[g+1, :] = params_TVM[:, 9]

    if len(np.where(propose_model==0)[0])>0:
        params_TVP_update = params_TVP[:,[0,1,4,9]]
        modelinex = 0
        # update TP paras and weights
        weight_gamma_TVP,weight_a_TVP,weight_nu_TVP,weight_vm_TVP,propose_gamma_TVP,propose_a_TVP,propose_nu_TVP,propose_vm_TVP=\
        tp_update(previous_bestfitted_index_TVP, propose_model, params_TVP_update, weight_gamma_TVP,
                  weight_a_TVP, weight_nu_TVP, weight_vm_TVP,modelinex)
        modelTVP = np.where(propose_model==modelinex)
        params_TVP = np.tile(sampleparam_TVP, (len(modelTVP[0]), 1))
        params_TVP[:, 0] = propose_gamma_TVP
        params_TVP[:, 1] = propose_a_TVP
        params_TVP[:, 4] = propose_nu_TVP
        params_TVP[:, 9] = propose_vm_TVP

    if len(np.where(propose_model==1)[0])>0:
        params_TV_update = params_TV[:,[0,1,4,9]]
        modelinex = 1
        if len(valid_TV)>0:
            # update TP paras and weights
            weight_gamma_TV,weight_a_TV,weight_nu_TV,weight_vm_TV,propose_gamma_TV,propose_a_TV,propose_nu_TV,propose_vm_TV=\
                tp_update(previous_bestfitted_index_TV, propose_model, params_TV_update, weight_gamma_TV,
                          weight_a_TV, weight_nu_TV, weight_vm_TV,modelinex)
            modelTV = np.where(propose_model==modelinex)
            params_TV = np.tile(sampleparam_TV, (len(modelTV[0]), 1))
            params_TV[:, 0] = propose_gamma_TV
            params_TV[:, 1] = propose_a_TV
            params_TV[:, 4] = propose_nu_TV
            params_TV[:, 9] = propose_vm_TV
        else:
            params_TV = np.tile(sampleparam_TV, (population, 1))  # duplicate
            params_TV[:, 0] = np.random.uniform(0, prior[1], params_TV.shape[0])  # randomize 'gamma'
            params_TV[:, 1] = np.random.uniform(0, prior[3], params_TV.shape[0])  # randomize 'a'
            params_TV[:, 4] = np.random.uniform(0, prior[5], params_TV.shape[0])  # randomize 'nu'
            params_TV[:, 9] = np.random.uniform(prior[6], prior[7], params_TV.shape[0])  # randomize 'Vm'

if len(np.where(propose_model==2)[0])>0:
        params_TVM_update = params_TVM[:,[0,1,4,9]]
        modelinex = 2
        # update TP paras and weights
        weight_gamma_TVM,weight_a_TVM,weight_nu_TVM,weight_vm_TVM,propose_gamma_TVM,propose_a_TVM,propose_nu_TVM,propose_vm_TVM=\
            tp_update(previous_bestfitted_index_TVM, propose_model, params_TVM_update, weight_gamma_TVM,
                      weight_a_TVM, weight_nu_TVM, weight_vm_TVM,modelinex)
        modelTVM = np.where(propose_model==modelinex)
        params_TVM = np.tile(sampleparam_TVM, (len(modelTVM[0]), 1))
        params_TVM[:, 0] = propose_gamma_TVM
        params_TVM[:, 1] = propose_a_TVM
        params_TVM[:, 4] = propose_nu_TVM
        params_TVM[:, 9] = propose_vm_TVM


#

para_data = {'model_data': model_data,'fitness': fitness, 'gamma_data_TVP':gamma_data_TVP,
             'a_data_TVP':a_data_TVP,'nu_data_TVP':nu_data_TVP,'vm_data_TVP':vm_data_TVP,
              'gamma_data_TV':gamma_data_TV,
             'a_data_TV':a_data_TV,'nu_data_TV':nu_data_TV,'vm_data_TV':vm_data_TVP,
             'gamma_data_TVM': gamma_data_TVM,
             'a_data_TVM': a_data_TVP, 'nu_data_TVM': nu_data_TVM, 'vm_data_TVM': vm_data_TVM,
             }


np.save(savedir,para_data)
