# AUTOGENERATED! DO NOT EDIT! File to edit: nbs/00_Core.ipynb (unless otherwise specified).


from __future__ import print_function


__all__ = ['RData', 'RegionExtractor', 'MarkerMaker', 'LinkageWriter', 'EncoderWorker', 'get_family_with_var',
           'phasing_haps', 'run_each_region', 'haplotyper']

# Cell
#nbdev_comment from __future__ import print_function
from .Utils import *
from .Runner import *
from multiprocessing import Process, Queue
from collections import OrderedDict
import itertools
from copy import deepcopy
import sys, faulthandler, platform
import numpy as np
import pandas as pd
import pickle
import os
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from itertools import repeat
if sys.version_info.major == 2:
    from cstatgen import cstatgen_py2 as cstatgen
    from cstatgen.egglib import Align
else:
    from cstatgen import cstatgen_py3 as cstatgen
    import egglib
    from egglib import Align

# Cell
class RData(dict):
    def __init__(self, vcf, tfam,anno_file=None,fam_pop_file=None,ind_sample_file=None,allele_freq_info = None,included_variant_file=None):
        # tfam.samples: a dict of {sid:[fid, pid, mid, sex, trait], ...}
        # tfam.families: a dict of {fid:[s1, s2 ...], ...}
        self.tfam = TFAMParser(tfam)
        # name of allele frequency meta info
        self.af_info = allele_freq_info
        self.vs = self.load_vcf(vcf)
        self.fam_pop = self.load_fam_info(fam_pop_file)
        self.anno = self.load_anno(anno_file,included_variant_file)
        self.samples_vcf = self.vs.GetSampleNames()
        self.samples_not_vcf = checkSamples(self.samples_vcf, self.tfam.samples.keys())[1]
        # samples have to be in both vcf and tfam data
        self.samples = OrderedDict([(k, self.tfam.samples[k]) for k in self.samples_vcf if k in self.tfam.samples])
        # a dict of {fid:[member names], ...}
        self.families = {k : [x for x in self.samples if x in self.tfam.families[k]] for k in self.tfam.families}
        # a dict of {fid:[idx ...], ...}
        self.famsampidx = {}
        # a dict of {fid:[maf1, maf2 ...]}
        self.maf = OrderedDict()
        # finalized sub_regions that are compliant to all families
        self.complied_markers = []
        # finalized sub regions (variants)
        self.combined_regions = []
        self.coordinates_by_region = []
        # RV varnames by family
        self.varnames_by_fam = {}
        self.patterns={}
        self.gnomAD_estimate={'AFR':(1-0.4589)/(2*7652),'AMR':(1-0.4455)/(2*16791),'ASJ':(1-0.2357)/(2*4925),'EAS':(1-0.4735)/(2*8624),'FIN':(1-0.3048)/(2*11150),'NFE':(1-0.5729)/(2*55860),'OTH':(1-0.4386)/(2*2743),'SAS':(1-0.5624)/(2*15391)}
        # reorder family samples based on order of VCF file
        for k in list(self.families.keys()):
            if len(self.families[k]) == 0:
                # skip families having no samples in VCF file
                del self.families[k]
            else:
                self.famsampidx[k] = [i for i, x in enumerate(self.samples_vcf) if x in self.families[k]]
        # a dict of {fid:[idx ...], ...}
        self.famvaridx = {}
        self.famvarmafs = {}
        self.wtvar = {}
        self.freq_by_fam = {}
        self.include_vars = []
        self.total_varnames={}
        self.total_mafs={}
        self.wt_maf={}
        self.freq = []
        self.genotype_all={}
        self.mle_mafs={}
        self.missing_persons=[]
        self.reset()

    def load_vcf(self,vcf):
        # load VCF file header
        return cstatgen.VCFstream(vcf)
    def load_anno(self,anno_file,included_variant_file=None):
        if anno_file is None:
            return None
        anno = pd.read_csv(anno_file)
        anno.index = list(anno.Otherinfo1)
        anno = anno[~anno.index.duplicated()]
        if included_variant_file:
            included_variants= pd.read_csv(included_variant_file)
            anno = anno.loc[included_variants,:]
        tmp = anno[list(set(self.fam_pop.values()))]
        tmp = tmp.replace('.',np.nan)
        tmp = tmp.replace(0,np.nan)
        anno = pd.concat([anno[['Chr','Start']],tmp.astype(np.float64)],axis=1)
        anno = anno.dropna()
        return anno
    def load_fam_info(self,fam_pop_file):
        if fam_pop_file is None:
            return None
        fam_pop = {}
        with open(fam_pop_file) as f:
            for line in f:
                key, value = line.split()
                if value == 'NA':   #Fixme: deal with missing info
                    fam_pop[key]=self.af_info
                else:
                    fam_pop[key] = value
        return fam_pop
    def load_ind_samples(self,ind_sample_file):
        pass

    def get_regions(self,step=1000):
        '''separate chromosome to regions'''
        regions=[]
        chrom=self.anno.Chr.unique()[0]
        for i,s in enumerate(self.anno.Start):
            if i==0:
                pre=None
                cur=s
            elif i%step==0:
                pre=cur
                cur=s
                regions.append([str(chrom),str(pre),str(cur),'R'+str(pre)+'_'+str(cur),'.', '.', '.'])
        if cur!=s:
            pre=cur
            cur=s
            regions.append([str(chrom),str(pre),str(cur),'R'+str(pre)+'_'+str(cur),'.', '.', '.'])
        return regions

    def reset(self):
        for item in self.tfam.samples: #for all samples in fam ( with or without vcfs)
            self[item] = []
            self.genotype_all[item] = []
        self.variants = []
        self.include_vars = []
        self.total_varnames={}
        self.total_mafs={}
        self.chrom = None
        for k in self.families.keys():
            self.famvaridx[k] = []
            self.famvarmafs[k] = []
        self.maf = OrderedDict()
        # superMarkerCount is the max num. of recombinant fragments among all fams
        self.superMarkerCount = 0
        self.complied_markers = []
        self.combined_regions = []
        self.coordinates_by_region = []
        self.patterns={}
        self.missing_persons=[]
        self.gss = {} #test line


    def getMidPosition(self):
        if len(self.variants) == 0:
            return None
        return sum([x[1] for x in self.variants]) / len(self.variants)

    def getFamVariants(self, fam, style = None):
        if style is None:
            return [item for idx, item in enumerate(self.variants) if idx in self.famvaridx[fam]]
        elif style == "map":
            names = []
            pos = []
            for idx in self.famvaridx[fam]:
                names.append(self.variants[idx][0])
                pos.append(self.variants[idx][1])
            mafs = self.famvarmafs[fam]
            return np.array(names), pos, np.array(mafs)  #pos can't be array -> TypeError: in method 'HaplotypingEngine_Execute'
        else:
            raise ValueError("Unknown style '{}'".format(style))

    def getFamSamples(self, fam):
        nvar = len([item for idx, item in enumerate(self.variants) if idx in self.famvaridx[fam]])
        output = [[]] * len(self.tfam.families[fam])
        for idx, item in enumerate(self.tfam.sort_family(fam)):
            # sample info, first 5 columns of ped
            output[idx] = self.tfam.samples[item][:-1]
            # sample genotypes
            if item in self.samples:
                output[idx].extend(self[item])
            else:
                output[idx].extend(["00"] * nvar)
        return output

# Cell
class RegionExtractor:
    '''Extract given genomic region from VCF
    converting genotypes into dictionary of
    genotype list'''
    def __init__(self, filename, build = None, chr_prefix = None):
        self.vcf = cstatgen.VCFstream(filename)
        self.chrom = self.startpos = self.endpos = self.name = None
        self.chr_prefix = chr_prefix
        if build is None:
            build = env.build
        self.xchecker = PseudoAutoRegion('X', build)
        self.ychecker = PseudoAutoRegion('Y', build)

    def apply(self, data):
        # Clean up
        data.reset()
        data.chrom = self.chrom
        self.vcf.Extract(self.chrom, self.startpos, self.endpos)
        if data.anno is None:
            varIdx=self.extract_vcf(data)
        else:
            varIdx=self.extract_vcf_with_anno(data)
        if varIdx == 0:
            return 1
        else:
            with env.variants_counter.get_lock():
                env.variants_counter.value += varIdx
            return 0

    def extract_vcf(self,data):
        varIdx = 0
        # for each variant site
        while (self.vcf.Next()):
            # check if the line's sample number matches the entire VCF sample number
            if not self.vcf.CountSampleGenotypes() == self.vcf.sampleCount:
                raise ValueError('Genotype and sample mismatch for region {}: {:,d} vs {:,d}'.\
                             format(self.name, self.vcf.CountSampleGenotypes(), self.vcf.sampleCount))
            # skip tri-allelic sites
            if not self.vcf.IsBiAllelic():
                with env.triallelic_counter.get_lock():
                    env.triallelic_counter.value += 1
                continue
            # valid line found, get variant info
            try:
                maf = float(self.vcf.GetInfo(data.af_info))
                if maf>0.5: #skip variants with af>0.5
                    continue
            except Exception as e:
                raise ValueError("VCF line {}:{} does not have valid allele frequency field {}!".\
                                 format(self.vcf.GetChrom(), self.vcf.GetPosition(), data.af_info))

            # for each family assign member genotype if the site is non-trivial to the family
            for k in data.families:
                gs = self.vcf.GetGenotypes(data.famsampidx[k])
                if len(set(''.join([x for x in gs if x != "00"]))) <= 1:
                    # skip monomorphic gs
                    continue
                else:
                    # this variant is found in the family
                    data.famvaridx[k].append(varIdx)
                    data.famvarmafs[k].append(maf if maf < 0.5 else 1-maf)
                    for person, g in zip(data.families[k], gs):
                        data[person].append(g if maf<0.5 else self.reverse_genotypes(g))
            data.variants.append([self.vcf.GetVariantID(), self.vcf.GetPosition(), self.name]) #remove maf
            varIdx += 1
        return varIdx

    def extract_vcf_with_anno(self,data):
        '''extract variants and annotation by region'''
        if str(data.anno.Chr[0])!=self.chrom:
            return 0
        anno_idx = (data.anno.Start>=self.startpos) & (data.anno.Start<self.endpos)
        if anno_idx.any()==False:
            return 0
        varmafs = data.anno[anno_idx]
        varIdx = 0
        i = -1
        # for each variant site
        while (self.vcf.Next()):
            i += 1
            # check if the line's sample number matches the entire VCF sample number
            if not self.vcf.CountSampleGenotypes() == self.vcf.sampleCount:
                raise ValueError('Genotype and sample mismatch for region {}: {:,d} vs {:,d}'.\
                             format(self.name, self.vcf.CountSampleGenotypes(), self.vcf.sampleCount))
            # skip tri-allelic sites
            if not self.vcf.IsBiAllelic():
                with env.triallelic_counter.get_lock():
                    env.triallelic_counter.value += 1
                continue
            # valid line found, get variant info
            try:
                mafs=varmafs.loc[self.vcf.GetVariantID()][2:]
            except:
                print(self.vcf.GetVariantID(), 'is not in annotation')
                continue
            if mafs.any()==False:
                continue

            # for each family assign member genotype if the site is non-trivial to the family
            fam_mafs=[]
            for k in data.families:
                gs = self.vcf.GetGenotypes(data.famsampidx[k])
                if len(set(''.join([x for x in gs if x != "00"]))) <= 1:
                    # skip monomorphic gs
                    continue
                else:
                    maf=mafs[data.fam_pop[k]]
                    fam_mafs.append(maf)
                    if maf and maf<0.5: #skip variants with af>0.5
                        # this variant is found in the family
                        data.famvaridx[k].append(varIdx)
                        data.famvarmafs[k].append(maf if maf < 0.5 else 1-maf)
                        for person, g in zip(data.families[k], gs):
                            data[person].append(g if maf<0.5 else self.reverse_genotypes(g))
            data.variants.append([self.vcf.GetVariantID(), self.vcf.GetPosition(), self.name]) #remove maf
            #print(i,varmafs.shape,self.chrom, self.startpos, self.endpos, self.name,self.vcf.GetPosition())
            varIdx += 1
        return varIdx

    def check_gs(self,gs):
        '''skip monomorphic variants and singleton variants in a family'''
        cg={'00':0,'11':0, '12':0, '22':0}
        for i in gs:
            cg[i]+=1
        not00 = cg['11']+cg['12']+cg['22']
        if cg['11']==not00 or cg['12']==not00 or cg['22']==not00:
            #skip monomorphic variants
            return False
        if cg['12']+cg['22']<=1:
            #skip sington variants
            return False
        return True

    def reverse_genotypes(self,g):
        ''' 11->22,12,21,22->11 '''
        if g=='11':
            g='22'
        elif g=='22':
            g='11'
        return g

    def getRegion(self, region):
        self.chrom, self.startpos, self.endpos, self.name = region[:4]
        self.startpos = int(self.startpos)
        self.endpos = int(self.endpos) + 1
        if self.chrom in ['X','23']:
            if self.xchecker.check(self.startpos) or self.xchecker.check(self.endpos):
                self.chrom = 'XY'
        if self.chrom in ['Y','24']:
            if self.ychecker.check(self.startpos) or self.ychecker.check(self.endpos):
                self.chrom = 'XY'
        if self.chr_prefix and not self.chrom.startswith(self.chr_prefix):
            self.chrom = self.chr_prefix + self.chrom


class MarkerMaker:
    def __init__(self, wsize, maf_cutoff = None, recomb=False):
        self.missings = ("0", "0")
        self.gtconv = {'1':0, '2':1}
        self.haplotyper = cstatgen.HaplotypingEngine(verbose = env.debug)
        if wsize == 0 or wsize >= 1:
            self.r2 = None
        else:
            self.r2 = wsize
        self.coder = cstatgen.HaplotypeCoder(wsize)
        self.maf_cutoff = maf_cutoff
        self.rsq = 0.0
        self.recomb = recomb
    def getRegion(self, region):
        self.name = region[3]
        self.dtest = {}
        self.dtest[self.name] = {}
        self.dtest[self.name]['predata']={}

    def apply(self, data):
        #try:
            # haplotyping plus collect found allele counts
            # and computer founder MAFS
        varnames,mafs,haplotypes=self.__Haplotype(data)
        if len(varnames)==0:
            return -1
        if not any([len(varnames[x]) - 1 for x in varnames]):
            # all families have only one variant
            self.__AssignSNVHaplotypes(data, haplotypes, mafs, varnames)
        else:
            # calculate LD clusters using founder haplotypes
            clusters = self.__ClusterByLD(data, haplotypes, varnames)
            # recoding the genotype of the region
            #env.dtest[self.name]['coder']['input'] = [data.copy(), haplotypes, mafs, varnames, clusters]
            self.__CodeHaplotypes(data, haplotypes, mafs, varnames, clusters)
            self.dtest[self.name]['maf']=data.maf
            self.dtest[self.name]['hap']=self.haps
            #env.dtest[self.name]['coder']['output'] = [self.coder.GetHaplotypes(),data.copy(),data.superMarkerCount,deepcopy(data.maf)]
        #except Exception as e:
        #    return -1
        self.__FormatHaplotypes(data)
        #env.dtest[self.name]['format'] = data.copy()
        return 0

    def __Haplotype(self, data):
        '''genetic haplotyping. haplotypes stores per family data'''
        # FIXME: it is SWIG's (2.0.12) fault not to properly destroy the object "Pedigree" in "Execute()"
        # So there is a memory leak here which I tried to partially handle on C++
        #
        # Per family haplotyping
        #
        items = get_family_with_var(data)
        with ProcessPoolExecutor(max_workers = 8) as executor:
            inputs = executor.map(phasing_haps,repeat(data.chrom),items,[data.getFamVariants(item, style = "map") for item in items],[data.getFamSamples(item) for item in items])

        varnames,mafs,haplotypes = OrderedDict(),OrderedDict(),OrderedDict()
        for item,item_varnames,item_mafs,item_haplotypes in inputs:
            if len(item_haplotypes) == 0:
                # C++ haplotyping implementation failed
                with env.chperror_counter.get_lock():
                    env.chperror_counter.value += 1
                    env.log('{} family failed to phase haplotypes.'.format(item))
                for person in data.families[item]:
                    data[person] = self.missings
                    continue
            self.dtest[self.name]['predata'][item]=[item_varnames,item_mafs,item_haplotypes]
            # Drop some variants if maf is greater than given threshold
            if self.maf_cutoff is not None:
                keep_idx = item_mafs<self.maf_cutoff
                if not keep_idx.any():
                    for person in data.families[item]:
                        data[person] = self.missings
                    continue
                item_mafs = item_mafs[keep_idx]
                item_varnames = item_varnames[keep_idx]
                item_haplotypes = item_haplotypes[:,np.concatenate(([True,True],keep_idx))]
            varnames[item],mafs[item],haplotypes[item]= item_varnames,item_mafs,item_haplotypes
        return varnames,mafs,haplotypes


    def __ClusterByLD(self, data, haplotypes, varnames):
        if self.r2 is None:
            return None
        # get founder haplotypes
        founder_haplotypes = []
        markers = sorted(set(itertools.chain(*varnames.values())), key = lambda x: int(x.split("-")[0][1:]))
        for item in haplotypes:
            for ihap, hap in enumerate(haplotypes[item]):
                if not data.tfam.is_founder(hap[1]):
                    continue
                gt = [hap[2 + list(varnames[item]).index(v)] if v in varnames[item] else '?' for v in markers]
                founder_haplotypes.append(("{}-{}".format(hap[1], ihap % 2), "".join([x[1] if x[0].isupper() else x[0] for x in gt])))
        # calculate LD blocks, use r2 measure
        blocks = []
        if sys.version_info.major == 2:
            ld = Align.create(founder_haplotypes).matrixLD(validCharacters="12")["r2"]
            for j in ld:  #upper triangle
                block = [j]
                for k in ld[j]:
                    try:
                        if ld[j][k] > self.r2:
                            block.append(k)
                    except:
                        print('ld value',ld[j][k])
                if len(block) > 1:
                    blocks.append(block)
        else:
            ldi,ld = egglib.stats.matrix_LD(Align.create(founder_haplotypes,egglib.Alphabet(cat='string',expl=['1','2'],miss='?')),('rsq'))
            for j in range(len(ldi)): #lower triangle
                block = [j]
                for k in range(j+1,len(ldi)):
                    try:
                        if ld[k][j] > self.r2:
                            block.append(k)
                    except:
                        print('ld value',ld[k][j])
                if len(block) > 1:
                    blocks.append(block)
        # get LD clusters
        clusters = [[markers[idx] for idx in item] for item in list(connected_components(blocks))]
        #env.dtest[self.name]['ld'] = [ld,blocks,clusters]
        if env.debug:
            with env.lock:
                print("LD blocks: ", blocks, file = sys.stderr)
                print("LD clusters: ", clusters, file = sys.stderr)
        return clusters


    def __CodeHaplotypes(self, data, haplotypes, mafs, varnames, clusters):
        # apply CHP coding
        if clusters is not None:
            clusters_idx = [[[list(varnames[item]).index(x) for x in y if x in varnames[item]] for y in clusters] for item in haplotypes]
        else:
            clusters_idx = [[[]] for item in haplotypes]
        self.coder.Execute(list(haplotypes.values()), [mafs[item] for item in haplotypes], clusters_idx,self.recomb)
        if env.debug:
            with env.lock:
                if clusters:
                    print("Family LD clusters: ", clusters_idx, "\n", file = sys.stderr)
                self.coder.Print()
        # line: [fid, sid, hap1, hap2]
        self.haps = {}
        for line in self.coder.GetHaplotypes():
            #if not line[1] in data:
                # this sample is not in VCF file. Every variant site should be missing
                # they have to be skipped for now
            #    continue
            data[line[1]] = (line[2].split(','), line[4].split(','))
            self.haps[line[0]] = self.haps.get(line[0], line)
            if len(data[line[1]][0]) > data.superMarkerCount:
                data.superMarkerCount = len(data[line[1]][0])
        # get MAF
        for item in haplotypes:
            data.maf[item] = self.coder.GetAlleleFrequencies(item)
            data.maf[item] = tuple(tuple(np.array(v) / np.sum(v)) if np.sum(v) else v
                              for v in data.maf[item])
        if env.debug:
            with env.lock:
                print("marker freqs = ", data.maf, "\n", file = sys.stderr)


    def __AssignSNVHaplotypes(self, data, haplotypes, mafs, varnames):
        print('SNVHap',self.name)
        for item in haplotypes:
            # each person's haplotype
            token = ''
            for idx, line in enumerate(haplotypes[item]):
                if not idx % 2:
                    token = line[2][1] if line[2][0].isupper() else line[2][0]
                else:
                    data[line[1]] = (token, line[2][1] if line[2][0].isupper() else line[2][0])
            # get maf
            data.maf[item] = [(1 - mafs[item][0], mafs[item][0])]
            data.maf[item] = tuple(tuple(np.array(v) / np.sum(v)) if np.sum(v) else v
                              for v in data.maf[item])
        if env.debug:
            with env.lock:
                print("marker freqs = ", data.maf, "\n", file = sys.stderr)


    def __FormatHaplotypes(self, data):
        # Reformat sample genotypes
        for person in data:
            if type(data[person]) is not tuple:
                data[person] = self.missings
                continue
            diff = data.superMarkerCount - len(data[person][0])
            data[person] = list(zip(*data[person]))
            if diff > 0:
                data[person].extend([self.missings] * diff)

    def __PedToHaplotype(self, ped):
        '''convert prephased ped format to haplotype format.
        Input: e.g. [['13346', '5888', '0', '0', '1', '11', '11', '11'], ['13346', '5856', '0', '0', '2', '12', '12', '12'], ['13346', '5920', '5888', '5856', '1', '12', '12', '12'], ['13346', '6589', '5888', '5856', '1', '11', '11', '11']]
        Output: e.g. (('13346', '5856', '1:', '1:', '1:'), ('13346', '5856', '2:', '2:', '2:'), ('13346', '5888', '1:', '1:', '1:'), ('13346', '5888', '1:', '1:', '1:'), ('13346', '6589', '1:', '1|', '1|'), ('13346', '6589', '1:', '1|', '1|'), ('13346', '5920', '2:', '2|', '2|'), ('13346', '5920', '1:', '1|', '1|'))
        '''
        haps = []
        for item in ped:
            entry = [item[0], item[1]] + [x[0] + ':' if x[0] != '0' else '?:' for x in item[5:]]
            haps.append(tuple(entry))
            entry = [item[0], item[1]] + [x[1] + ':' if x[1] != '0' else '?:' for x in item[5:]]
            haps.append(tuple(entry))
        return tuple(haps)


class LinkageWriter:
    def __init__(self, num_missing_append = 0):
        self.chrom = self.prev_chrom = self.name = self.distance = self.distance_avg = self.distance_m = self.distance_f = None
        self.reset()
        self.missings = ["0", "0"]
        self.num_missing = num_missing_append

    def apply(self, data):
        if self.chrom != self.prev_chrom:
            if self.prev_chrom is None:
                self.prev_chrom = self.chrom
            else:
                # new chrom entered,
                # commit whatever is in buffer before accepting new data
                self.commit()
        # write tped output
        position = str(data.getMidPosition())
        if data.superMarkerCount <= 1:
            # genotypes
            gs = [data[s][0] for s in data.tfam.samples]
            if len(set(gs)) == 1:
                # everyone's genotype is the same (most likely missing or monomorphic)
                return 2
            self.tped += env.delimiter.join([self.chrom, self.name, self.distance, position] + \
                list(itertools.chain(*gs)) + self.missings*self.num_missing) + "\n"
            # freqs
            for k in data.maf:
                self.freq += env.delimiter.join([k, self.name] + list(map(str, data.maf[k][0]))) + "\n"
        else:
            # have to expand each region into mutiple chunks to account for different recomb points
            gs = list(zip(*[data[s] for s in data.tfam.samples]))
            # sub-chunk id
            cid = 0
            skipped_chunk = []
            for idx, g in enumerate(gs):
                if len(set(g)) == 1:
                    skipped_chunk.append(idx)
                    continue
                cid += 1
                self.tped += \
                  env.delimiter.join([self.chrom, '{}[{}]'.format(self.name, cid), self.distance, position] + \
                  list(itertools.chain(*g)) + self.missings*self.num_missing) + "\n"
            if cid == 0:
                # everyone's genotype is the same (most likely missing or monomorphic)
                return 2
            # freqs
            for k in data.maf:
                cid = 0
                for idx in range(data.superMarkerCount):
                    if idx in skipped_chunk:
                        continue
                    if idx >= len(data.maf[k]):
                        break
                    cid += 1
                    self.freq += env.delimiter.join([k, '{}[{}]'.format(self.name, cid)] + \
                                                    list(map(str, data.maf[k][idx]))) + "\n"
        if self.counter < env.batch:
            self.counter += data.superMarkerCount
        else:
            self.commit()
        return 0

    def commit(self):
        if self.tped:
            with env.lock:
                with open(os.path.join(env.tmp_cache, '{}.chr{}.tped'.format(env.output, self.prev_chrom)),
                          'a') as f:
                    f.write(self.tped)
        if self.freq:
            with env.lock:
                with open(os.path.join(env.tmp_cache, '{}.chr{}.freq'.format(env.output, self.prev_chrom)),
                          'a') as f:
                    f.write(self.freq)
        self.reset()

    def reset(self):
        self.tped = ''
        self.freq = ''
        self.counter = 0
        self.prev_chrom = self.chrom

    def getRegion(self, region):
        self.chrom = region[0]
        self.name, self.distance_avg, self.distance_m, self.distance_f = region[3:]
        self.distance = ";".join([self.distance_avg, self.distance_m, self.distance_f])

class EncoderWorker(Process):
    def __init__(self, queue, length, data, extractor, coder, writer):
        Process.__init__(self)
        self.queue = queue
        self.numGrps = float(length)
        self.data = data
        self.extractor = extractor
        self.maker = coder
        self.writer = writer

    def report(self):
        env.log('{:,d} units processed {{{:.2%}}} ...'.\
                format(env.success_counter.value, env.total_counter.value / self.numGrps), flush = True)

    def run(self):
        while True:
            try:
                region = self.queue.pop(0) if isinstance(self.queue, list) else self.queue.get()
                if region is None:
                    self.writer.commit()
                    self.report()
                    # total mendelian errors found
                    with env.mendelerror_counter.get_lock():
                        env.mendelerror_counter.value += self.maker.haplotyper.CountMendelianErrors()
                    # total recombination events found
                    with env.recomb_counter.get_lock():
                        env.recomb_counter.value += self.maker.coder.CountRecombinations()
                    break
                else:
                    with env.total_counter.get_lock():
                        env.total_counter.value += 1
                    self.extractor.getRegion(region)
                    self.writer.getRegion(region)
                    self.maker.getRegion(region)
                    isSuccess = True
                    for m in [self.extractor, self.maker, self.writer]:
                        status = m.apply(self.data)
                        if status == -1:
                            with env.chperror_counter.get_lock():
                                # previous module failed
                                env.chperror_counter.value += 1
                        if status == 1:
                            with env.null_counter.get_lock():
                                env.null_counter.value += 1
                        if status == 2:
                            with env.trivial_counter.get_lock():
                                env.trivial_counter.value += 1
                        if status != 0:
                            isSuccess = False
                            break
                    if isSuccess:
                        with env.success_counter.get_lock():
                            env.success_counter.value += 1
                    if env.total_counter.value % (env.batch * env.jobs) == 0:
                        self.report()
            except KeyboardInterrupt:
                break

# Cell
def get_family_with_var(data):
    items = []
    for item,item_vars in data.famvaridx.items():
        if len(item_vars) == 0:  #no variants in the family
            for person in data.families[item]:
                data[person] = ("0", "0")
        else:
            items.append(item)
    return items

haplotyper = cstatgen.HaplotypingEngine(verbose = env.debug)
def phasing_haps(chrom,item,fvar,fgeno):
    item_varnames, positions, item_mafs = fvar
    try:
        item_haplotypes = haplotyper.Execute(chrom, item_varnames, positions, fgeno)[0]
    except:
        env.log("{} fail to phase haplotypes".format(item))
        item_haplotypes = []
    item_haplotypes = np.array(item_haplotypes)
    return item,item_varnames,item_mafs,item_haplotypes

def run_each_region(regions,data,extractor,maker,writer,genotype=True):
    '''get the haplotypes and allele frequency of variants in each region'''
    results = {}
    i=0
    start = time.perf_counter()
    for region in regions:
        extractor.getRegion(region)
        maker.getRegion(region)
        writer.getRegion(region)
        isSuccess = True
        for m in [extractor] if genotype else [extractor, maker, writer]:
            status = m.apply(data)
            if status == -1:
                with env.chperror_counter.get_lock():
                    # previous module failed
                    env.chperror_counter.value += 1
            if status == 1:
                with env.null_counter.get_lock():
                    env.null_counter.value += 1
            if status == 2:
                with env.trivial_counter.get_lock():
                    env.trivial_counter.value += 1
            if status != 0:
                isSuccess = False
                break
        if isSuccess:
            with env.success_counter.get_lock():
                env.success_counter.value += 1
            if genotype:
                #{'gene':{'predata':{'fam':[snp_ids,freq,genos]}}}
                items = get_family_with_var(data)
                predata={}
                for item in items:
                    fvar=data.getFamVariants(item, style = "map")
                    fgeno=np.array(data.getFamSamples(item))
                    predata[item]=[fvar[0],fvar[2],fgeno]
                results[region[3]]={'predata':predata}
            else:
                results[region[3]]=maker.dtest[region[3]]
            if len(results)==env.cache_size:
                env.log('write to pickle: '+os.path.join(env.tmp_cache,env.output+str(i)+'.pickle')+',Gene number:'+str(len(results))+',Time:'+str((time.perf_counter()-start)/3600))
                start = time.perf_counter()
                with open(os.path.join(env.tmp_cache,env.output+str(i)+'.pickle'), 'wb') as handle:
                    pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL)
                results = {}
                i +=1
                #add linkage analysis
    if len(results)>0:
        env.log('write to pickle: '+os.path.join(env.tmp_cache,env.output+str(i)+'.pickle')+',Gene number:'+str(len(results))+',Time:'+str((time.perf_counter()-start)/3600))
        with open(os.path.join(env.tmp_cache,env.output+str(i)+'.pickle'), 'wb') as handle:
            pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL)
        results = {}
        #add linkage analysis
