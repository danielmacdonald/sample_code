"""
The genetic algorithm component of gadrad. The genetic algorithm is used as an
optimisation approach used to determine which simulated radiative transfer
environment  best matches the true physical source.
"""

import random
import numpy as np
import params


def getKey(item):
    return item[1]

def index_find(k,num_pop,fitness):
    lis = []
    for _ in range(len(params.ga_params()['prob_array'])):
        i = random.randint(0,num_pop-1)
        lis.append([i,fitness[i]])
    return sorted(lis, key=getKey)

def parent_prob():
    probrand = random.random()
    neg_test = [probrand-i for i in params.ga_params()['prob_array']]
    return [n for n, i in enumerate(neg_test) if i<0][0]

def mut(chrom,chrom_ind,mut_prob,std_dev):
    for i in range(len(chrom)):
        if random.random() < mut_prob:
            max_index = len(params.var_params().values()[i])
            new_index = int(round(np.random.normal(chrom_ind[i],std_dev[i]+0.5)))
            while new_index < 0:
                new_index = new_index + max_index
            while new_index >= max_index:
                new_index = new_index - max_index
            chrom[i] = params.var_params().values()[i][new_index]
            chrom_ind[i] = new_index
        test = np.append([chrom],[chrom_ind])
    return test


def init():
    values_array = params.var_params().values()
    init_chromo = []
    for i in range(len(values_array)):
        r = random.randint(0,len(values_array[i])-1)
        a = values_array[i][r]
        init_chromo.append([r,a])
    return [[row[1] for row in init_chromo],[row[0] for row in init_chromo]]

def chromo_keys():
    chromo_keys = params.var_params().keys()
    return chromo_keys

def row_in_array(myarray, myrow):
    return (myarray == myrow).all(-1).any()

def roullette_selection(chromes):
    elite_keep = params.ga_params()['elite_keep']
    pop_size = params.ga_params()['pop_size']

    parents_sort = chromes[chromes[:,0].argsort()][0:int(elite_keep*pop_size)]
    np.savetxt('parents_sort',parents_sort)

    fitnesses = [1.0/i for i in parents_sort[:,0]]
    num = len(fitnesses)

    tot_fit = float(np.sum(fitnesses))
    rel_fit = [f/tot_fit for f in fitnesses]

    probs = [sum(rel_fit[:i+1]) for i in range (len(rel_fit))]
    new_population = []
    for i in range(pop_size):
        r = random.random()
        for i in range(pop_size):
            if r <= probs[i]:
                new_population.append(parents_sort[i])
                break
    return new_population

def full_mut(values, indexes,len_chr,num_pop,mut_rate):

    for individual in range(num_pop):
        for i in range(len_chr):
            new_index = int(round(np.random.normal(indexes[individual][i],mut_rate)))
            if new_index >= 0 and new_index < len(params.var_params().values()[i]):
                values[individual][i] = params.var_params().values()[i][new_index]
                indexes[individual][i] = new_index
            else:
                values[individual][i] = values[individual][i]
                indexes[individual][i] = indexes[individual][i]

    full_mut = np.concatenate((np.array(values),np.array(indexes)),axis=1)
#np.array([list(values[ix])+list(indexes[ix]) for ix in range(len(chrom))])

    return np.array(full_mut)

def mutation_convergence(gen, mut_const):
    mut = mut_const/(float(gen)**(1.0/5.0))
    return mut

def std_dev(array):
    std_dev_arr = []
    for i in range(len(array[0])):
    #    xj = array[:,i]
    #    xj = np.append([xj],[[params.var_params().values()[i][0]],[params.var_params().values()[i][1]]])
        std_dev_arr.append(np.std(array[:,i]))
    return std_dev_arr


def children(fit_chromes):
    outfile2 = open('fit_chromes','w')
    print >> outfile2, fit_chromes
    count = 0
    fitness = fit_chromes[:,0]
   #print 'Best', min(fitness)
    num_pop = len(fitness)
    len_chr = int((len(fit_chromes[0])-2)/2.0)
    gen = fit_chromes[:,1][0]
    mut_std_dev = params.ga_params()['mut_std_dev']
    mut_prob = params.ga_params()['prob_mut']
    mut_rate = mutation_convergence(gen,mut_std_dev)
    outfile999 = open('mut_rate','w')
    print >> outfile999, mut_rate
    fittest = fit_chromes[np.argmin(fit_chromes[:,0])]
    fittest = np.array(fittest[2:len_chr+len_chr+2])
    outfile6 = open('fittest','w')
    print >> outfile6, fittest

    std_deviation = std_dev(fit_chromes[:,17:,])
    outfile919 = open('std_deviation','w')
    print >> outfile919, std_deviation

    parents = np.array([[None]*(len_chr*2)])
    parents = np.vstack((parents,fittest))
    while parents.shape[0] <= num_pop:

        ind_1 = index_find(len(params.ga_params()['prob_array']),num_pop,fitness)
        index_1 = [item[0] for item in ind_1]
        par_index_1 = parent_prob()
        p1 = fit_chromes[index_1[par_index_1],:][2:len_chr+2]
        p1_ind= fit_chromes[index_1[par_index_1],:][len_chr+2:len_chr+2+len_chr]

        if random.random() < params.ga_params()['cross_prob']:

            ind_2 = index_find(len(params.ga_params()['prob_array']),num_pop,fitness)
            index_2 = [item[0] for item in ind_2]
            par_index_2 = parent_prob()
                                      #decide which of the parents out of num k becomes parent
            p2 = fit_chromes[index_2[par_index_2],:][2:len_chr+2]
            p2_ind = fit_chromes[index_2[par_index_2],:][len_chr+2:len_chr+2+len_chr]
            # crossover and mutation
            #crossover points
            c = sorted(random.sample(xrange(0,len_chr), 2))
            #crossover
            p1[c[0]:c[1]] = p2[c[0]:c[1]].copy()
            p1_ind[c[0]:c[1]] = p2_ind[c[0]:c[1]].copy()
            #mutation
            p0 = mut(p1,p1_ind,mut_prob,std_deviation)

        else:
            p0 = mut(p1,p1_ind,mut_prob,std_deviation)

        parents = np.vstack((parents,p0))

    parents = np.delete(parents,0,0)
    outfile = open('parents','w')
    print >> outfile, parents
    return np.array(parents)
