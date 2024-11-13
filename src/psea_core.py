import argparse
import matplotlib
matplotlib.use('Agg')
import sys
import math
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import gridspec
from scipy import stats
from scipy.stats import norm 
import random
import random
import pandas as pd

#set seed for testing
np.random.seed(seed=42)

def pcea_score(pc_score, verbose=False):
    '''Calculate the GSEA like enrichment score using the comorbidity 
    occurance in the rank as our set
    
    Parameters
    ----------
    pc_score : list
        with the binary presence or absence of comorbidity
    
    verbose : boolean
        print all the summary statements (default = False)
    
    Returns
    -------
    actual_es : float
        the enrichment score (AUC) between the enrichment curve and the 
        background
        
    normalized_pc_score : list
        normalized patient comorbidity occurance
        
    trend : list
        the trend line for each set
        
    cumscore : list
        cumulative normalized comorbidity hit score
        
    '''
    
    #total patients with a comorbidity
    #and making sure the list is one hot incoded
    pc_score = [1 if i>0 else 0 for i in pc_score]
    total = float(sum(pc_score))
    
    #normalize the positive hits to the total patients with a comorbidity 
    normalized_pc_score = [(float(x)/total) for x in pc_score]

    
    #calculate the cumulative normalized patient comorbidity score
    #the _MAX_ should add to 1
    cumscore = np.cumsum(normalized_pc_score)
    
    #this is the background line
    trend = np.arange(0,1,1.0/float(len(cumscore) - 1))
    
    #enrichment score (ES) is the actual score subtracted from background trend
    actual_es = np.trapz(cumscore) - np.trapz(trend)

    if verbose==True:
        print("Sum of normalized binary scores :"+str(sum(normalized_pc_score)))
        print("Max cumulative score            :"+str(max(cumscore)))
        print("Max of trend                    :"+str(max(trend)))
        print("Len of trend                    :"+str(len(trend)))
        print("Actual PCEA                     :"+str(actual_es))
    
    return [actual_es, normalized_pc_score, trend, cumscore]

def pcea_score_norm(pc_score, verbose=False):
    '''Calculate the GSEA like enrichment score using the comorbidity 
    occurance in the rank as our set
    
    Parameters
    ----------
    pc_score : list
        with the binary presence or absence of comorbidity

    verbose : boolean
        print all the summary statements (default = False)

    Returns
    -------
    actual_es : float
        the enrichment score (AUC) between the enrichment curve and the 
        background
        
    normalized_pc_score : list
        normalized patient comorbidity occurance
        
    trend : list
        the trend line for each set
        
    cumscore : list
        cumulative normalized comorbidity hit score
        
    '''
    
    #number of patients in the dataset
    #num_patients = len(pc_score)
    
    #total patients with a comorbidity
    pc_score = [1 if i>0 else 0 for i in pc_score] #binarize score for normalized scores
    total = float(sum(pc_score))
    
    #get number of bins
    binwidth = 1.0/float(total) 
    
    #normalize the positive hits to the total patients with a comorbidity 
    #normalized_pc_score = [(float(x)/total)*binwidth for x in pc_score]
    normalized_pc_score = np.multiply(np.divide(pc_score, total), binwidth)

    #calculate the cumulative normalized patient comorbidity score
    #the _SUM_ should add to 1
    cumscore = np.cumsum(normalized_pc_score)
    
    #this is the background line
    trend = np.append(np.arange(0,1,1.0/float(len(cumscore) - 1)), 1.0)
    trend = np.multiply(trend, binwidth)
    
    #enrichment score (ES) is the actual score subtracted from background
    actual_es = (np.trapz(cumscore) - np.trapz(trend)) *2 
    
    if verbose==True:
        print("Binwidth                        :"+str(binwidth))
        print("Sum of normalized binary scores :"+str(sum(normalized_pc_score)))
        print("Sum cumulative score            :"+str(sum(cumscore)))
        print("Sum of trend                    :"+str(sum(trend)))
        print("Len of trend                    :"+str(len(trend)))
        print("Actual PCEA                     :"+str(actual_es))
    
    return [actual_es, normalized_pc_score, trend, cumscore]


def pcea_score_norm(pc_score, verbose=False):
    '''Calculate the GSEA like enrichment score using the comorbidity 
    occurance in the rank as our set
    
    Parameters
    ----------
    pc_score : list
        with the binary presence or absence of comorbidity

    verbose : boolean
        print all the summary statements (default = False)

    Returns
    -------
    actual_es : float
        the enrichment score (AUC) between the enrichment curve and the 
        background
        
    normalized_pc_score : list
        normalized patient comorbidity occurance
        
    trend : list
        the trend line for each set
        
    cumscore : list
        cumulative normalized comorbidity hit score
        
    '''
    
    #number of patients in the dataset
    #num_patients = len(pc_score)
    
    #total patients with a comorbidity
    pc_score = [1 if i>0 else 0 for i in pc_score] #binarize score for normalized scores
    total = float(sum(pc_score))
    
    #get number of bins
    binwidth = 1.0/float(len(pc_score)) 
    
    #normalize the positive hits to the total patients with a comorbidity 
    #normalized_pc_score = [(float(x)/total)*binwidth for x in pc_score]
    normalized_pc_score = np.multiply(np.divide(pc_score, total), binwidth)

    #calculate the cumulative normalized patient comorbidity score
    #the _SUM_ should add to 1
    cumscore = np.cumsum(normalized_pc_score)
    
    #this is the background line
    trend = np.append(np.arange(0,1,1.0/float(len(cumscore) - 1)), 1.0)
    trend = np.multiply(trend, binwidth)
    
    #enrichment score (ES) is the actual score subtracted from background
    actual_es = (np.trapz(cumscore) - np.trapz(trend)) *2 
    
    if verbose==True:
        print("Binwidth                        :"+str(binwidth))
        print("Sum of normalized binary scores :"+str(sum(normalized_pc_score)))
        print("Sum cumulative score            :"+str(sum(cumscore)))
        print("Sum of trend                    :"+str(sum(trend)))
        print("Len of trend                    :"+str(len(trend)))
        print("Actual PCEA                     :"+str(actual_es))
    
    return [actual_es, normalized_pc_score, trend, cumscore]

def permute_pcea(ranks, permutations=1000, seed=42):
    '''Generates permutations of the ranks and calculates AUC for each 
        permutation.
        
    Parameters
    ----------
    ranks : list or array
        normalized comorbidity ranks to be shuffled
        
    permutations : int
        number of times to permute (default=1000)
        
    seed : int
        seed for reproducibility (default=42)
        
    Returns
    -------
    es_permute : list 
        list of AUC calculated for N permutations 
       
    '''
    
    #set a seed
    np.random.seed(seed=seed)
    
    es_permute = list()
    for i in range(permutations):

        random_pc_score = np.random.permutation(ranks)
        es_permute_score = pcea_score(random_pc_score)
        es_permute.append(es_permute_score[0])

    return es_permute


def permute_pcea_norm(ranks, permutations=1000, seed=42):
    '''Generates permutations of the ranks and calculates AUC for each 
        permutation.
        
    Parameters
    ----------
    ranks : list or array
        normalized comorbidity ranks to be shuffled
        
    permutations : int
        number of times to permute (default=1000)
        
    seed : int
        seed for reproducibility (default=42)
    Returns
    -------
    es_permute : list 
        list of AUC calculated for N permutations 
       
    '''
    
    #set a seed
    np.random.seed(seed=seed)    
    
    
    
    es_permute = list()
    for i in range(permutations):

        random_pc_score = np.random.permutation(ranks)
        es_permute_score = pcea_score_norm(random_pc_score) #using the normalized PCEA
        es_permute.append(es_permute_score[0])

    return es_permute


def calculateNESpval(actualES, simES):
    if actualES < 0:
            simESsubset = [x for x in simES if x < 0]
            mu = np.mean(simESsubset)
            NES = -(actualES/mu)
            sigma = np.std(simESsubset)
            p = stats.norm.cdf(actualES, mu, sigma)
    else:
            simESsubset = [x for x in simES if x > 0]
            mu = np.mean(simESsubset)
            NES = actualES/mu
            sigma = np.std(simESsubset)
            p = 1-stats.norm.cdf(actualES, mu, sigma)
    return NES, p

def plot_pcea(cumscore, trend, ranks, comorbidity_binary, sim_pcea, pcea, n_bins = 100):
    
    #####################################
    #Begin plotting section             #
    #####################################
    F = plt.figure(figsize=(15.5,8))
    xvals = range(1,len(cumscore)+1)
    limits = [1,len(cumscore)]
    gs = gridspec.GridSpec(4, 1, height_ratios=[2, 0.5, 0.2, 2])

    #This is the enrichment score plot (i.e. line plot)
    ax0 = plt.subplot(gs[0])
    ax0.plot(xvals, cumscore,color='green')
    ax0.plot(trend, '--', alpha=0.5)
    ax0.set_title('Patient Comorbidity Enrichment Score Plot',fontsize=14)
    ax0.set_ylabel('Enrichment Score (ES)', fontsize=12)
    ax0.tick_params(axis='y', which='both', left='on', right='off', labelleft='on')
    ax0.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
    ylims = ax0.get_ylim()
    ymax = math.fabs(max(ylims,key=abs))
    ax0.set_ylim([0,ymax])
    ax0.set_xlim(limits)

    #####################################
    #Heatmap for comorbidity occurance  #
    #####################################
    ax1 = plt.subplot(gs[1])
    ax1.bar(ranks, comorbidity_binary, color ='gray', 
            width = 0.5, alpha=0.5)
    ax1.set_ylabel(' ', fontsize=12)
    ax1.set_xlabel('Patient Ranks', fontsize=12)
    ax1.tick_params(axis='y', which='both', left='off', right='off', labelleft='off')
    ax1.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
    ax1.margins(x=0)

    #####################################
    #Random vs. actual ES               #
    #####################################
    ax2 = plt.subplot(gs[3])
    num_bins = n_bins
    # the histogram of the data
    n, bins, patches = ax2.hist(sim_pcea, num_bins, density=True,edgecolor='black', linewidth=0.75)

    ax2.set_xlabel('PCEA-Score', fontsize=12)
    ax2.set_ylabel('Probability density', fontsize=12)
    ax2.set_title('Distribution of Simulated PCEA-Score', fontsize=14)
    ax2.axvline(x = pcea, color = 'red', label = 'axvline - full height')
    fig.tight_layout()
    plt.show()

def load_file(infilename):
    gene_cormorbid_data = pd.read_csv(infilename)
    return gene_cormorbid_data


def write_psea(infilename, outputfile):
    gene_cormorbid_data = load_file(infilename)
    onegeneNES, onegenep = psea_heart(gene_cormorbid_data, outputfile)
    line = [infilename, onegeneNES, onegenep]
    wf = open(outputfile, "w")
    wf.write("\t".join(map(str,line))+"\n")
    wf.close()

def is_unique(s):
    a = s.to_numpy() # s.values (pandas<0.24)
    return (a[0] == a).all()

def psea_heart(gene_cormorbid_data, outfilename):
    sample = gene_cormorbid_data["sample"].to_list()
    if is_unique(gene_cormorbid_data["value"])==True:
        return "ba_or_value_are_all_identical", "excluded_by_psea"
    elif is_unique(gene_cormorbid_data["binary_attribute"])==True:
        return "ba_or_value_are_all_identical", "excluded_by_psea"
    else:
        value = gene_cormorbid_data["value"].to_list()
        binary_attribute = gene_cormorbid_data["binary_attribute"].to_list()
        df = gene_cormorbid_data[["sample", "value", "binary_attribute"]]
        df = df.sample(frac=1).reset_index(drop=True)
        df = df.sort_values("value")
        thiscomorbidity_binary = df["binary_attribute"].to_list()
        actualES_norm, normalized_pc_score_norm, thistrend_norm, thiscumscore_norm = pcea_score_norm(thiscomorbidity_binary)
        simES_norm = permute_pcea_norm(normalized_pc_score_norm, permutations=1000)
        onegeneNES, onegenep = calculateNESpval(actualES_norm, simES_norm)
        return onegeneNES, onegenep

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run psea on a single set of values and binary attributes')
    parser.add_argument('-i', '--inputfile') 
    parser.add_argument('-o', '--outputfile') 
    args = parser.parse_args()
    inputfile = args.inputfile
    outputfile = args.outputfile
    write_psea(inputfile, outputfile)
