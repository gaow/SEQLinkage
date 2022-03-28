# AUTOGENERATED! DO NOT EDIT! File to edit: nbs/linkage_analysis.ipynb (unless otherwise specified).

__all__ = ['paramlink2', 'pedprobr', 'pedtools', 'get_allele', 'name_haps', 'get_fam_hap', 'format_haps_bunch',
           'calculate_ped_lod', 'parallel_lods', 'sum_variant_lods']

# Cell
import numpy as np
import pandas as pd
import pickle
from itertools import repeat

#Import necessary packages
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
#Must be activated
pandas2ri.activate()
paramlink2=importr('paramlink2')
pedprobr=importr('pedprobr')
pedtools = importr('pedtools')

import time
from concurrent.futures import ProcessPoolExecutor

# Cell
def get_allele(s):
    a = s[1] if s[0].isupper() else s[0]
    return 0 if a=='?' else int(a)

def name_haps(snps):
    name = []
    for i in snps:
        name += [i+'_A0',i+'_A1']
    return name

def get_fam_hap(haps,name=None):
    new_haps,new_iid = [],[]
    iid = haps[:,1]
    haps = haps[:,2:]
    for i in range(0,haps.shape[0],2):
        new_iid.append(iid[i])
        hap_a01 = []
        for a0,a1 in zip(haps[i],haps[i+1]):
            hap_a01 += [get_allele(a0),get_allele(a1)]
        new_haps.append(hap_a01)
    new_haps = pd.DataFrame(new_haps)
    new_haps.index = new_iid
    if name is not None:
        new_haps.columns = name
    return new_haps

# Cell
def format_haps_bunch(dhaps,fam):
    gene_variants = {}
    gene_haps = {}
    for g in dhaps.keys():
        haps = dhaps[g]['predata']
        for f in haps.keys():
            if f not in gene_variants.keys():
                gene_variants[f] = {'genes':[],'variants':[],'freqs':[],'uniq':[]}
                gene_haps[f] = get_fam_hap(haps[f][2],name_haps(haps[f][0]))
            else:
                gene_haps[f] = pd.concat([gene_haps[f],get_fam_hap(haps[f][2],name_haps(haps[f][0]))],axis=1)
            gene_variants[f]['genes'] += [g]*len(haps[f][0])
            gene_variants[f]['variants'] += list(haps[f][0])
            gene_variants[f]['freqs'] += list(haps[f][1])
    for i,j in gene_variants.items():
        redup_idx = ~gene_haps[i].columns.duplicated()
        gene_haps[i] = pd.concat([fam[i],gene_haps[i].iloc[:,redup_idx]],axis=1)
        j['uniq'] = list(redup_idx[range(0,len(redup_idx),2)])
        gene_variants[i] = pd.DataFrame(j)
    return gene_variants,gene_haps

def calculate_ped_lod(ped,rho=0,model = "AD",chrom = "AUTOSOMAL",penetrances = [0.01,0.9,0.9],dfreq=0.001):
    aff=ped.iloc[:,5]
    mped = pedtools.as_ped(ped.drop(ped.columns[5], axis=1),famid_col = 1,id_col = 2,fid_col = 3,mid_col = 4,sex_col = 5)
    modAD = paramlink2.diseaseModel(model,chrom,pd.Series(penetrances),dfreq)
    res = paramlink2.lod(mped, aff = aff, model = modAD,rho=rho)
    try:
        res = pd.DataFrame(res)[['MARKER','LOD']]
    except:
        res = pd.DataFrame([[ped.columns[6],res[0]]],columns=['MARKER','LOD'])
    return res

def parallel_lods(haps,rho=0):
    start = time.perf_counter()
    with ProcessPoolExecutor(max_workers = 10) as executor:
        results = executor.map(calculate_ped_lod,haps,repeat(rho))
    print(time.perf_counter()-start)
    return results

def sum_variant_lods(lods):
    variants = {}
    for lod in lods:
        for m,l in zip(lod['MARKER'],lod['LOD']):
            if m in variants.keys():
                variants[m] += l
            else:
                variants[m] = l
    var_lst = []
    for var,lod in variants.items():
        snp = var[:-3]
        var_lst.append(snp.split(':')+[snp,lod])
    variants=pd.DataFrame(var_lst,columns=['CHR','POS','A0','A1','SNP','LOD'])
    variants.POS = variants.POS.astype(int)
    variants.sort_values('POS')
    return variants