# AUTOGENERATED! DO NOT EDIT! File to edit: nbs/02_Linkage.ipynb (unless otherwise specified).

__all__ = ['base', 'paramlink2', 'pedprobr', 'pedtools', 'update_haps_ped', 'hap2chp', 'generate_marker',
           'recombination_pos', 'recombination_region', 'get_allele', 'name_haps', 'get_fam_hap', 'get_fam_geno',
           'format_haps_bunch', 'hlod_fun', 'calculate_ped_lod', 'parallel_lods', 'linkage_analysis', 'format_fam_lods',
           'rows2one', 'get_lods_batch', 'get_lods_chrom', 'get_hlod_chrom', 'summarize_lods']

# Cell
import sys
import os.path
import glob
import numpy as np
import pandas as pd
import pickle
from itertools import repeat
import numbers

#Import necessary packages
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
base = importr('base')
base.options(expressions = 5e5)
#Must be activated
pandas2ri.activate()
paramlink2=importr('paramlink2')
pedprobr=importr('pedprobr')
pedtools = importr('pedtools')

import time
from concurrent.futures import ProcessPoolExecutor,ThreadPoolExecutor
from scipy.optimize import minimize_scalar

# Cell
def update_haps_ped(genes):
    gpeds={}
    for g,fs in genes.items():
        fs=fs['predata']
        tmp={}
        for f,v in fs.items():
            tmp[f]=hap2chp(v,g)
        if len(tmp)>0:
            gpeds[g]={}
            gpeds[g]['predata']=tmp
    return gpeds

def hap2chp(haps,gene):
    '''input hap is [varnames,freqs,halpotypes]
       output is [chp_varnames,chp_freqs,chps]
    '''
    if len(haps[0])>1:
        hap=haps[2][:,2:]
        rec=recombination_pos(hap)
        new_hap=[]
        for genos in hap:
            new_genos=[]
            for i,j in zip(rec[:-1],rec[1:]):
                new_genos.append(generate_marker(np.array([get_allele(x) for x in genos[i:j]])))
            new_hap.append(new_genos)
        new_hap=np.concatenate([haps[2][:,:2],np.array(new_hap)],axis=1)
        variants=np.array([gene+':'+str(i) for i in range(len(rec)-1)])
        mafs=np.array([1-np.prod(1-haps[1][i:j]) for i,j in zip(rec[:-1],rec[1:])])
    else:
        variants,mafs,new_hap=np.array([gene+':0']),haps[1],haps[2]
    return [variants,mafs,new_hap]

def generate_marker(alleles):
    '''array of 0,1,2. 0 if all 0 -> 2 if any 2 else 1'''
    if np.all(alleles==0):
        return '0:'
    elif np.any(alleles==2):
        return '2:'
    else:
        return '1:'

def recombination_pos(hap):
    rec = [0,hap.shape[1]-1]
    for genos in hap:
        for c,geno in enumerate(genos):
            if geno[-1] not in [':','|']:
                rec.append(c)
    return list(set(rec))
def recombination_region(hap):
    rec=recombination_pos(hap)
    return [(i,j) for i,j in zip(rec[:-1],rec[1:])]

# Cell
def get_allele(s):
    a = s[1] if s[0].isupper() else s[0]
    return 0 if a=='?' else int(a)

def name_haps(snps):
    name = []
    for i in snps:
        name += [i+'_A0',i+'_A1']
    return name

def get_fam_hap(haps,variants,vcf=None):
    new_haps,new_iid = [],[]
    iid = haps[:,1]
    haps = haps[:,2:]
    for i in range(0,haps.shape[0],2):
        cur_iid=iid[i]
        new_iid.append(cur_iid)
        if vcf is None or vcf[cur_iid]:#have vcf
            hap_a01 = []
            for a0,a1 in zip(haps[i],haps[i+1]): #loop through variants
                hap_a01 += [get_allele(a0),get_allele(a1)]
        else:
            hap_a01 = [0,0]*haps.shape[1] #set missing vcf to 0
        new_haps.append(hap_a01)
    new_haps = pd.DataFrame(new_haps)
    new_haps.index = new_iid
    new_haps.columns = name_haps(variants)
    #remove variants with only 1 or 2 as alleles, return None
    idx=[]
    for i in range(0,new_haps.shape[1],2):
        v=set(new_haps.iloc[:,i]).union(set(new_haps.iloc[:,i+1]))
        if 1 not in v or 2 not in v:
            idx.append(False)
        else:
            idx.append(True)
    if sum(idx)==0:
        return None
    return new_haps.loc[:,np.repeat(np.array(idx),2)],idx

def get_fam_geno(haps,variants,vcf=None):
    new_haps,new_iid = [],[]
    iid = haps[:,1]
    haps = haps[:,5:]
    for i in range(haps.shape[0]):
        cur_iid=iid[i]
        new_iid.append(cur_iid)
        if vcf is None or vcf[cur_iid]:#have vcf
            hap_a01 = []
            for a01 in haps[i]: #loop through variants
                hap_a01 += [int(a) for a in a01]
        else:
            hap_a01 = [0,0]*haps.shape[1] #set missing vcf to 0
        new_haps.append(hap_a01)
    new_haps = pd.DataFrame(new_haps)
    new_haps.index = new_iid
    new_haps.columns = name_haps(variants)
    #remove variants with only 1 or 2 as alleles, return None
    idx=[]
    for i in range(0,new_haps.shape[1],2):
        v=set(new_haps.iloc[:,i]).union(set(new_haps.iloc[:,i+1]))
        if 1 not in v or 2 not in v:
            idx.append(False)
        else:
            idx.append(True)
    if sum(idx)==0:
        return None
    return new_haps.loc[:,np.repeat(np.array(idx),2)],idx

def format_haps_bunch(dhaps,fam,vcfs=None,cutoff=None,haplotype=True):
    gene_variants = {}
    gene_haps = {}
    for g in dhaps.keys():
        haps = dhaps[g]['predata']
        with ProcessPoolExecutor(max_workers = 10) as executor:
            if haplotype:
                results = executor.map(get_fam_hap,[haps[k][2] for k in haps.keys()],[haps[k][0] for k in haps.keys()],[vcfs[k] if vcfs else None for k in haps.keys()])
            else:
                results = executor.map(get_fam_geno,[haps[k][2] for k in haps.keys()],[haps[k][0] for k in haps.keys()],[vcfs[k] if vcfs else None for k in haps.keys()])
        for f,hap in  zip(haps.keys(),results):
            if hap is None: #remove only have 1 or 2 variants
                continue
            if f not in gene_variants.keys():
                gene_variants[f] = {'genes':[],'variants':[],'freqs':[]}
                gene_haps[f] = hap[0]
            else:
                gene_haps[f] = pd.concat([gene_haps[f],hap[0]],axis=1)
            idx=hap[1] #False for variants only have 1 or 2.
            gene_variants[f]['genes'] += [g]*sum(idx)
            gene_variants[f]['variants'] += list(haps[f][0][idx])
            gene_variants[f]['freqs'] += list(haps[f][1][idx])
    new_gene_variants,new_gene_haps={},{}
    for i,j in gene_variants.items():
        j=pd.DataFrame(j)
        if cutoff is not None:
            frq_idx=np.array(j['freqs'])<cutoff
            if frq_idx.any()==False:
                continue
            j=j.loc[frq_idx,:]
            gene_haps[i]=gene_haps[i].loc[:,np.repeat(frq_idx,2)]
        redup_idx = ~gene_haps[i].columns.duplicated()
        new_gene_haps[i] = pd.concat([fam[i],gene_haps[i].iloc[:,redup_idx]],axis=1)
        j['uniq'] = list(redup_idx[range(0,len(redup_idx),2)])
        new_gene_variants[i] = j
    return new_gene_variants,new_gene_haps

# Cell
def hlod_fun(Li, sign=1):
    def _fun(alpha):
        return sign * sum(np.log10(alpha*np.power(10, Li) + 1 - alpha))
    return _fun

# Cell
def calculate_ped_lod(ped,afreq=None,rho=0,model = "AD",chrom = "AUTOSOMAL",penetrances = [0.01,0.9,0.9],dfreq=0.001):
    def _calculate_ped_lod(mped, aff, model,rho):
        res = paramlink2.lod(mped, aff, model,rho)
        try:
            res = pd.DataFrame(res)[['MARKER','LOD']]
        except:
            res = pd.DataFrame([[ped.columns[6],res[0]]],columns=['MARKER','LOD'])
        return res
    aff=ped.iloc[:,5]
    mped = pedtools.as_ped(ped.drop(ped.columns[5], axis=1),famid_col = 1,id_col = 2,fid_col = 3,mid_col = 4,sex_col = 5)
    if afreq is not None:
        mped = pedtools.setLocusAttributes(mped,locusAttributes=[base.list(afreq=base.c(1-i,i)) for i in afreq])
    modAD = paramlink2.diseaseModel(model,chrom,pd.Series(penetrances),dfreq)
    if isinstance(rho,numbers.Number):
        res = _calculate_ped_lod(mped, aff = aff, model = modAD,rho=rho)
    else:
        res=None
        for r in rho:
            tmp = _calculate_ped_lod(mped, aff = aff, model = modAD,rho=r)
            if res is None:
                res=tmp
                res.columns = ['MARKER','LOD'+str(round(r,2))]
            else:
                res['LOD'+str(round(r,2))]=tmp.LOD
        res.index=list(res.MARKER)
        res=res.iloc[:,1:]
    return res

def parallel_lods(haps,afreqs=None,rho=0,model = "AD",chrom = "AUTOSOMAL",penetrances = [0.01,0.9,0.9],dfreq=0.001):
    start = time.perf_counter()
    if afreqs is None:
        with ProcessPoolExecutor(max_workers = 10) as executor:
            results = executor.map(calculate_ped_lod,haps.values(),repeat(rho),repeat(model),repeat(chrom),repeat(penetrances),repeat(dfreq))
    else:
        with ProcessPoolExecutor(max_workers = 10) as executor:
            results = executor.map(calculate_ped_lod,haps.values(),afreqs,repeat(rho),repeat(model),repeat(chrom),repeat(penetrances),repeat(dfreq))
    print(time.perf_counter()-start)
    return {k:res for k,res in zip(haps.keys(),results)}

def linkage_analysis(gene_genotype_file,fam,fam_vcf,cutoff,chp=True,rho=np.arange(0,0.5,0.05),model = "AD",chrom = "AUTOSOMAL",penetrances = [0.01,0.9,0.9],dfreq=0.001):
    '''linkage analysis function'''
    linkage_input_file=gene_genotype_file[:-7]+'_AFcutoff'+str(cutoff)+'_linkage.input'
    lod_file=linkage_input_file[:-6]+'.lods'
    hlod_file=linkage_input_file[:-6]+'.hlods'
    besthlod_file=linkage_input_file[:-6]+'.besthlod'
    #preprocess genotypes or phased haplotypes to the format of linkage analysis
    if os.path.isfile(linkage_input_file):
        print('exist! jump',linkage_input_file,file=sys.stderr)
    else:
        print('create',linkage_input_file,file=sys.stdout)
        with open(gene_genotype_file, 'rb') as handle:
            genes = pickle.load(handle)
        if chp: #making CHP markers from phased haplotypes
            genes=update_haps_ped(genes)
            cutoff=None
        gene_variants,gene_fam_haps = format_haps_bunch(genes,fam,fam_vcf,cutoff,chp)
        if len(gene_variants)==0:
            print(gene_geotype_file,'No variants left after filtering')
            return
        with open(linkage_input_file,'wb') as handle:
            pickle.dump([gene_variants,gene_fam_haps], handle, protocol=pickle.HIGHEST_PROTOCOL)

    #linkage analysis and write out results to .lods file
    if os.path.isfile(lod_file):
        print('exist! jump',lod_file,file=sys.stderr)
    else:
        print('create',lod_file,file=sys.stdout)
        with open(linkage_input_file, 'rb') as handle:
            gene_variants,gene_fam_haps = pickle.load(handle)
        afreqs = []
        for k in gene_fam_haps.keys():
            variants= gene_variants[k]
            variants=variants.freqs[variants.uniq]
            afreqs.append(variants)
        res = parallel_lods(gene_fam_haps,afreqs,rho,model,chrom,penetrances,dfreq) #R function for linkage analysis
        with open(lod_file,'wb') as handle:
            pickle.dump(res, handle, protocol=pickle.HIGHEST_PROTOCOL)

    #heterogeneity analysis and write out results to .hlods and .besthlods files
    if os.path.isfile(besthlod_file):
        print('exist! jump',besthlod_file,file=sys.stderr)
    else:
        print('create',besthlod_file,file=sys.stdout)
        with open(lod_file, 'rb') as handle:
            res = pickle.load(handle)
        var_res=format_fam_lods(res.values(),prefix=chp)
        var_sovs,best_sovs=[],[]
        for var,res in var_res.items():
            res=res.fillna(0)
            best_sov=[var,'LOD0.5',0,0]
            for theta in res.index:
                try:
                    sov = minimize_scalar(hlod_fun(list(res.loc[theta]), -1), bounds=(0,1), method='bounded', options={'xatol':1e-8})
                    var_sov=[var,theta,sov.x,-sov.fun]
                except:
                    var_sov=[var,theta,0,0]
                var_sovs.append(var_sov)
                if best_sov[3]<var_sov[3]:
                    best_sov=var_sov
            best_sovs.append(best_sov)
        var_sovs=pd.DataFrame(var_sovs)
        best_sovs=pd.DataFrame(best_sovs)
        with open(hlod_file,'wb') as handle:
            pickle.dump(var_sovs, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(besthlod_file,'wb') as handle:
            pickle.dump(best_sovs, handle, protocol=pickle.HIGHEST_PROTOCOL)

# Cell
def format_fam_lods(res,cutoff=0,prefix=False):
    new_res,variants=[],[]
    for i in res:
        idx=i.index
        if prefix:
            idx=[i.split(':',1)[0] for i in idx]
            i.index=idx
        else:
            idx=[i[:-3] for i in idx]
            i.index=idx
        new_res.append(i)
        variants.append(i.index)
    variants = list(set().union(*variants))
    var_res={}
    for v in variants:
        varlods = [rows2one(r.loc[v]) for r in new_res if v in r.index]
        if len(varlods)>cutoff:
            #v=snp:informative fam
            v=':'.join([v,str(len(varlods))])
            var_res[v]=pd.concat(varlods,axis=1)
    return var_res

def rows2one(x):
    if len(x.shape)>1:
        return x.max(axis=0)
    else:
        return x

def get_lods_batch(path,fams=None,phase=False):
    with open(path, 'rb') as handle:
        res=pickle.load(handle)
        if fams is not None:
            try:
                res=[r for k,r in res.items() if k in fams]
            except:
                with open(path[:-4]+'input', 'rb') as handle:
                    gene_variants,gene_fam_haps=pickle.load(handle)
                res={k:r for k,r in zip(gene_fam_haps.keys(),res)}
                res=[r for k,r in res.items() if k in fams]
        else:
            if type(res) is dict:
                res=res.values()
    res_d=format_fam_lods(res,prefix=phase)
    #sum lods among families
    lods=pd.concat([x.sum(axis=1) for x in res_d.values()],axis=1).T
    lods.index=list(res_d.keys())
    return lods

def get_lods_chrom(prefix,fams=None,phase=False):
    path_ress=glob.glob(prefix)
    ress = []
    for x in path_ress:
        try:
            res=get_lods_batch(x,fams,phase)
        except:
            continue
        if len(res)==0:
            continue
        ress.append(res)
    ress=pd.concat(ress)
    return ress[~ress.index.duplicated(keep=False)]

def get_hlod_chrom(prefix):
    path_ress=glob.glob(prefix)
    ress = []
    for x in path_ress:
        try:
            with open(x,'rb') as handle:
                res=pickle.load(handle)
        except:
            continue
        if len(res)==0:
            continue
        ress.append(res)
    ress=pd.concat(ress)
    ress.columns=['name','theta','alpha','hlod']
    ress.index=list(ress.name)
    return ress[~ress.index.duplicated(keep=False)]

def summarize_lods(input_lod,output_prefix,regions,fams=None,phase=False):
    '''two output files:one is lod scores from rho 0 to 0.5. another is lod at rho=0 and max lod score combined with chr, pos and name'''
    lods_chr = get_lods_chrom(input_lod,fams,phase)
    hlod_chr = get_hlod_chrom(input_lod[:-4]+'besthlod')
    lods_chr.sort_index().to_csv(output_prefix+'_lods.csv',header=True,index=True)

    # get the max lod for each gene
    info_fam=[int(i.rsplit(':',1)[1]) for i in lods_chr.index]
    lods_chr.index=[i.rsplit(':',1)[0] for i in lods_chr.index]
    hlod_chr.index=[i.rsplit(':',1)[0] for i in hlod_chr.index]
    if phase:
        genes=pd.DataFrame(regions).iloc[:,:4]
        genes.columns=['chrom','start','end','name']
        genes.index=list(genes.name)
        genes=genes[genes.index.isin(lods_chr.index)]
    else:
        genes=pd.DataFrame([i.split(':') for i in lods_chr.index])
        genes.columns=['chrom','pos','a0','a1']
        genes.index=list(lods_chr.index)

    lods_chr=lods_chr.loc[genes.index,:]
    hlod_chr=hlod_chr.loc[genes.index,:]
    genes['InfoFam']=info_fam
    genes['LOD0']=list(lods_chr['LOD0.0'])
    genes['LODmax']=list(lods_chr.max(axis=1))
    genes.loc[genes.LODmax<0,'LODmax']=0
    pd.concat([genes,hlod_chr.loc[:,['theta','alpha','hlod']]],axis=1).sort_index().to_csv(output_prefix+'_lod_summary.csv')