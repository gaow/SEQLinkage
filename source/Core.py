#!/usr/bin/python2.7
# Copyright (c) 2013, Gao Wang <wang.gao@columbia.edu>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)
from __future__ import print_function
from SEQLinkage.Utils import *
from SEQLinkage.Runner import *
from multiprocessing import Process, Queue
from collections import OrderedDict
import itertools
from copy import deepcopy
import sys, faulthandler, platform
import numpy as np

if sys.version_info.major == 2:
    from cstatgen import cstatgen_py2 as cstatgen
else:
    from cstatgen import cstatgen_py3 as cstatgen
from cstatgen.egglib import Align

def checkParams(args):
    '''set default arguments or make warnings'''
    env.debug = args.debug
    env.quiet = args.quiet
    env.prephased = args.prephased
    args.vcf = os.path.abspath(os.path.expanduser(args.vcf))
    args.tfam = os.path.abspath(os.path.expanduser(args.tfam))
    for item in [args.vcf, args.tfam]:
        if not os.path.exists(item):
            env.error("Cannot find file [{}]!".format(item), exit = True)
    if args.output:
        env.outdir = args.output
        env.output = os.path.split(args.output)[-1]
        env.cache_dir = os.path.join(os.path.dirname(args.output), 'cache')
        env.tmp_log = os.path.join(env.tmp_dir, env.output + ".STDOUT")
    #
    if len([x for x in set(getColumn(args.tfam, 6)) if x.lower() not in env.ped_missing]) > 2:
        env.trait = 'quantitative'
    env.log('{} trait detected in [{}]'.format(env.trait.capitalize(), args.tfam))
    if not args.blueprint:
        args.blueprint = os.path.join(env.resource_dir, 'genemap.{}.txt'.format(args.build))
    args.format = [x.lower() for x in set(args.format)]
    if args.run_linkage and "linkage" not in args.format:
        args.format.append('linkage')
    if None in [args.inherit_mode, args.prevalence, args.wild_pen, args.muta_pen] and "linkage" in args.format:
        env.error('To generate LINKAGE format or run LINKAGE analysis, please specify all options below:\n\t--prevalence, -K\n\t--moi\n\t--wild-pen, -W\n\t--muta-pen, -M', show_help = True, exit = True)
    if not args.tempdir is None:
        env.ResetTempdir(args.tempdir)
    return True

class RData(dict):
    def __init__(self, samples_vcf, tfam):
        # tfam.samples: a dict of {sid:[fid, pid, mid, sex, trait], ...}
        # tfam.families: a dict of {fid:[s1, s2 ...], ...}
        self.tfam = tfam
        # samples have to be in both vcf and tfam data
        self.samples = OrderedDict([(k, tfam.samples[k]) for k in samples_vcf if k in tfam.samples])
        # a dict of {fid:[member names], ...}
        self.families = {k : [x for x in self.samples if x in tfam.families[k]] for k in tfam.families}
        # a dict of {fid:[idx ...], ...}
        self.famsampidx = {}
        # a dict of {fid:[maf1, maf2 ...]}
        self.maf = OrderedDict()
        # finalized sub_regions that are compliant to all families
        self.complied_markers = []
        # finalized sub regions (variants)
        self.combined_regions = []
        # RV varnames by family
        self.varnames_by_fam = {}
        self.patterns={}
        self.gnomAD_estimate={'AFR':(1-0.4589)/(2*7652),'AMR':(1-0.4455)/(2*16791),'ASJ':(1-0.2357)/(2*4925),'EAS':(1-0.4735)/(2*8624),'FIN':(1-0.3048)/(2*11150),'NFE':(1-0.5729)/(2*55860),'OTH':(1-0.4386)/(2*2743),'SAS':(1-0.5624)/(2*15391)}
        # reorder family samples based on order of VCF file
        for k in self.families.keys():
            if len(self.families[k]) == 0:
                # skip families having no samples in VCF file
                del self.families[k]
            else:
                self.famsampidx[k] = [i for i, x in enumerate(samples_vcf) if x in self.families[k]]
        # a dict of {fid:[idx ...], ...}
        self.famvaridx = {}
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

    def reset(self):
        for item in self.samples:
            self[item] = []
            self.genotype_all[item] = []
        self.variants = []
        self.include_vars = []
        self.total_varnames={}
        self.total_mafs={}
        self.wt_maf={}
        self.chrom = None
        for k in self.families:
            self.famvaridx[k] = []
            self.wtvar[k] = []
        self.maf = OrderedDict()
        # superMarkerCount is the max num. of recombinant fragments among all fams
        self.superMarkerCount = 0
        self.complied_markers = []
        self.combined_regions = []
        self.patterns={}
        self.missing_persons=[]

    def getMidPosition(self):
        if len(self.variants) == 0:
            return None
        return sum([x[1] for x in self.variants]) / len(self.variants)

    def getFamVariants(self, fam, style = None, include_wt = False):
        if style is None:
            if include_wt:
                return [item for idx, item in enumerate(self.variants) if idx in self.famvaridx[fam]]
            else:
                return [item for idx, item in enumerate(self.variants) if idx in self.famvaridx[fam] and idx not in self.wtvar[fam]]
        elif style == "map":
            names = []
            pos = []
            mafs = []
            tmp_vars = self.famvaridx[fam]
            if not include_wt:
                tmp_vars=[idx for idx in self.famvaridx[fam] if idx not in self.wtvar[fam]]
            if len(self.freq_by_fam.keys()) != 0:
                pop_idx=self.freq.index(self.freq_by_fam[fam])
            for idx in tmp_vars:
                names.append("V{}-{}".format(idx, self.variants[idx][1]))
                pos.append(self.variants[idx][1])
                tmp_mafs=self.variants[idx][-1]
                if type(tmp_mafs) is list:
                    mafs.append(tmp_mafs[pop_idx])
                else:
                    mafs.append(tmp_mafs)
            return names, pos, mafs
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


class RegionExtractor:
    '''Extract given genomic region from VCF
    converting genotypes into dictionary of
    genotype list'''
    def __init__(self, filename, build = env.build, chr_prefix = None, allele_freq_info = None, include_vars_file=None):
        self.vcf = cstatgen.VCFstream(filename)
        self.chrom = self.startpos = self.endpos = self.name = None
        self.chr_prefix = chr_prefix
        # name of allele frequency meta info
        self.af_info = allele_freq_info
        self.xchecker = PseudoAutoRegion('X', build)
        self.ychecker = PseudoAutoRegion('Y', build)
        self.include_vars_file = include_vars_file

    def apply(self, data):
        # Clean up
        data.reset()
        data.chrom = self.chrom
        self.vcf.Extract(self.chrom, self.startpos, self.endpos)
        varIdx = 0
        # for each variant site
        while (self.vcf.Next()):
            # skip tri-allelic sites
            if not self.vcf.IsBiAllelic():
                with env.triallelic_counter.get_lock():
                    env.triallelic_counter.value += 1
                continue
            if len(data.variants) > 0:
                if self.vcf.GetPosition()==data.variants[-1][1]:
                    continue
            # check if the line's sample number matches the entire VCF sample number
            if not self.vcf.CountSampleGenotypes() == self.vcf.sampleCount:
                raise ValueError('Genotype and sample mismatch for region {}: {:,d} vs {:,d}'.\
                             format(self.name, self.vcf.CountSampleGenotypes(), self.vcf.sampleCount))
            # valid line found, get variant info
            try:
                if type(self.af_info) is list:
                    maf = []
                    large_maf = []
                    for pop_info in self.af_info:
                        large_maf.append(False)
                        try:
                            maf.append(float(self.vcf.GetInfo(pop_info)))
                        except ValueError:
                            maf.append(0.0)
                    for idx in range(len(maf)):
                        if maf[idx] > 0.5:
                            large_maf[idx]=True
                            maf[idx] = 1-maf[idx]
                else:
                    large_maf=False
                    try:
                        maf = float(self.vcf.GetInfo(self.af_info)) if self.af_info else None
                    except ValueError:
                        maf = 0.0
                    if maf > 0.5:
                        large_maf=True
                        maf = 1 - maf
            except Exception:
                raise ValueError("VCF line {}:{} does not have valid allele frequency field {}!".\
                                 format(self.vcf.GetChrom(), self.vcf.GetPosition(), self.af_info))
            data.variants.append([self.vcf.GetChrom(), self.vcf.GetPosition(), self.name, maf])
            # for each family assign member genotype if the site is non-trivial to the family
            for k in data.families:
                gs = self.vcf.GetGenotypes(data.famsampidx[k])
                if len(data.freq_by_fam) > 0:
                    popidx=self.af_info.index(data.freq_by_fam[k])
                    if large_maf[popidx]:
                        tmpgs=[]
                        for tmpg in gs:
                            if tmpg=='00':
                                tmpgs.append(tmpg)
                            else:
                                tmpgs.append(''.join([str(3-int(tg)) for tg in tmpg]))
                        gs=tuple(tmpgs)
                else:
                    if large_maf:
                        tmpgs=[]
                        for tmpg in gs:
                            if tmpg=='00':
                                tmpgs.append(tmpg)
                            else:
                                tmpgs.append(''.join([str(3-int(tg)) for tg in tmpg]))
                        gs=tuple(tmpgs)
                for person, g in zip(data.families[k], gs):
                    data.genotype_all[person].append(g)
                if len(set(''.join(gs))) <= 1:
                    # skip monomorphic gs
                    continue
                else:
                    if len(set(''.join([x for x in gs if x != "00"]))) <= 1:
                        data.wtvar[k].append(varIdx)
                    # this variant is found in the family
                    data.famvaridx[k].append(varIdx)
                    for person, g in zip(data.families[k], gs):
                        data[person].append(g)
            varIdx += 1
        #
        if varIdx == 0:
            return 1
        else:
            if not self.include_vars_file is None:
                with open(self.include_vars_file) as invar_fh:
                    for invar_line in invar_fh:
                        chrom, pos = invar_line.split()
                        for vidx,v in enumerate(data.variants):
                            if v[0] == chrom and v[1] == int(pos):
                                data.include_vars.append("{}".format(pos))
                                break
            else:
                data.include_vars = ["{}".format(item[1]) for item in data.variants]
            with env.variants_counter.get_lock():
                env.variants_counter.value += varIdx
            return 0


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
    def __init__(self, wsize, maf_cutoff = None,single_markers=False,recomb_max = 1,af_info=None,freq_by_fam=False,rsq=0.0,mle=False,rvhaplo=False):
        self.missings = ("0", "0")
        self.gtconv = {'1':0, '2':1}
        self.recomb_max = recomb_max
        self.haplotyper = cstatgen.HaplotypingEngine(verbose = env.debug)
        self.af_info = af_info
        self.freq_by_fam = freq_by_fam
        self.rsq=rsq
        self.mle=mle          #use MLE estimate from families for MAF
        self.count= not mle   #count founder alleles to estimate MAF
        self.rvhaplo=rvhaplo
        if wsize == 0 or wsize >= 1:
            self.r2 = None
        else:
            self.r2 = wsize
        self.coder = cstatgen.HaplotypeCoder(wsize)
        self.maf_cutoff = maf_cutoff
        self.single_markers = single_markers
        self.name = None

    def apply(self, data):
        # temp raw haplotype, maf and variant names data
        haplotypes = OrderedDict()
        mafs = {}   ##Per fam per variant
        uniq_vars = []
        exclude_vars = []
        varnames = {}
        recombPos = {}
        try:
            # haplotyping plus collect found allele counts
            # and computer founder MAFS
            self.__Haplotype(data, haplotypes, mafs, varnames,recombPos,uniq_vars,exclude_vars)
            if len(varnames):
                if not any ([len(varnames[x]) - 1 for x in varnames]):
                    # all families have only one variant
                    self.__AssignSNVHaplotypes(data, haplotypes, mafs, varnames)
                else:
                    # calculate LD clusters using founder haplotypes
                    clusters = self.__ClusterByLD(data, haplotypes, varnames)
                    # recoding the genotype of the region
                    self.__CodeHaplotypes(data, haplotypes, mafs, varnames, clusters)
        except Exception as e:
            if env.debug:
                raise
            return -1
        self.__FormatHaplotypes(data,recombPos,varnames,uniq_vars)
        return 0

    def __getMLEfreq(self,data, markers_to_analyze, pos_all, families, rsq, output_log):
        output_sample=[]
        mle_mafs={}
        if len(markers_to_analyze)==0:
            return mle_mafs
        for fam in families:
            for person in data.tfam.sort_family(fam):
                output_sample.append([])
                last_ele=len(output_sample)-1
                output_sample[last_ele] = data.tfam.samples[person][:-1]
                if person in data.samples:
                    for marker in markers_to_analyze:
                        idx=int(marker.split('-')[0][1:])
                        output_sample[last_ele].append(data.genotype_all[person][idx])
                else:
                    output_sample[last_ele].extend(["00"] * len(markers_to_analyze))
        with stdoutRedirect(to = output_log):
            af = self.haplotyper.Execute(data.chrom, markers_to_analyze, pos_all, output_sample, rsq, output_log,False)
        with open(output_log) as mle_fh:
            for line in mle_fh:
                if line.startswith('V'):
                    tmp_eles = line.split(':')
                    if tmp_eles[0] not in mle_mafs:
                        freqs=tmp_eles[1].split()
                        mle_maf = float(freqs[1])
                        if mle_maf>0.5:
                            mle_mafs[tmp_eles[0]]=float("%.9f"%(1-mle_maf))
                        else:
                            #alt allele is more frequent
                            mle_mafs[tmp_eles[0]]=float("%.9f"%mle_maf)
                            marker_idx=int(tmp_eles[0].split('-')[0][1:])
                            for fam in families:
                                if marker_idx not in data.famvaridx[fam]:
                                    continue
                                tmp_famvaridx=data.famvaridx[fam].index(marker_idx)
                                for person in data.families[fam]:
                                    tmpg=data.genotype_all[person][marker_idx]
                                    tmpg_switch=''.join([str(3-int(tg)) for tg in tmpg]) if tmpg!='00' else tmpg
                                    data.genotype_all[person][marker_idx]=tmpg_switch
                                    tmpg2=data[person][tmp_famvaridx]
                                    tmpg_switch2=''.join([str(3-int(tg)) for tg in tmpg2]) if tmpg2!='00' else tmpg2
                                    data[person][tmp_famvaridx]=tmpg_switch2
        return mle_mafs

    def __computefounderfreq(self,data, families):
        #count founder alleles to estimate MAF
        total_founder_alleles=0
        tmp_haplotypes=OrderedDict()
        tmp_mafs={}
        for item in families:
            tmp_haplotypes[item] = self.__PedToHaplotype(data.getFamSamples(item))
            # count founder alleles
            for hap in tmp_haplotypes[item]:
                if not data.tfam.is_founder(hap[1]):
                    continue
                total_founder_alleles+=1.0
                for idxv, v in enumerate(data.getFamVariants(item,style="map",include_wt=True)[0]):
                    if v not in tmp_mafs:
                        # [#alt, #haplotypes]
                        tmp_mafs[v] = [0, 0]
                    gt = hap[2 + idxv][1] if hap[2 + idxv][0].isupper() else hap[2 + idxv][0]
                    if not gt == "?":
                    #genotyped
                        tmp_mafs[v][0] += self.gtconv[gt]
                    else:
                    #genotype is missing
                        tmp_mafs[v][1] -= 1.0
        #compute MAFs based on counts
        for v in tmp_mafs:
            if type(tmp_mafs[v]) is not list:
                continue
            tmp_mafs[v] = tmp_mafs[v][0] / (tmp_mafs[v][1]+total_founder_alleles) if tmp_mafs[v][1]+total_founder_alleles > 0 else 0.0
        return tmp_mafs

    def __Haplotype(self, data, haplotypes, mafs, varnames,recombPos,uniq_vars,exclude_vars):
        '''genetic haplotyping. haplotypes stores per family data'''
        # FIXME: it is SWIG's (2.0.12) fault not to properly destroy the object "Pedigree" in "Execute()"
        # So there is a memory leak here which I tried to partially handle on C++
        #
        # Per family haplotyping
        #
        self.markers = ["V{}-{}".format(idx, item[1]) for idx, item in enumerate(data.variants)]
        tmp_mafs = {}
        if self.freq_by_fam:
            ## if families are from different populations
            ## estimate MAF by different population
            fam_to_analyze={}
            for fam,pop in data.freq_by_fam.iteritems():
                if pop not in fam_to_analyze:
                    fam_to_analyze[pop]=[fam]
                else:
                    fam_to_analyze[pop].append(fam)
        if self.count:
            ## estimate MAF by counting founder alleles
            if self.freq_by_fam:
                local_count_mafs={}
                for pop in fam_to_analyze:
                    local_count_mafs[pop]=self.__computefounderfreq(data,fam_to_analyze[pop])
            else:
                local_count_mafs=self.__computefounderfreq(data,data.families.keys())
        if self.mle:
            ## estimate MLE allele frequency using all fam
            local_mle_mafs={}
            if self.freq_by_fam:
                for pop in fam_to_analyze:
                    local_mle_mafs[pop]={}
                    markers_to_analyze=[]
                    pos_all=[]
                    markers_analyzed={}
                    if pop not in data.mle_mafs:
                        data.mle_mafs[pop]={}
                    else:
                        for tmpv in data.mle_mafs[pop]:
                            markers_analyzed[tmpv.split('-')[-1]]=data.mle_mafs[pop][tmpv]
                    output_log=env.tmp_log+"AF_{}_{}.log".format(pop,self.name)
                    popidx=self.af_info.index(pop)
                    variants_in_fams=[]
                    for item in fam_to_analyze[pop]:
                        for tmpvar in data.getFamVariants(item):
                            if tmpvar not in variants_in_fams:
                                variants_in_fams.append(tmpvar)
                    variants_in_fams=sorted(variants_in_fams, key=lambda x: x[1])
                    for item in variants_in_fams:
                        idx=data.variants.index(item)
                        if item[-1][popidx]==0:
                            if str(item[1]) in markers_analyzed.keys():
                                #if variant has been analyzed
                                vname="V{}-{}".format(idx,item[1])
                                local_mle_mafs[pop][vname]=markers_analyzed[str(item[1])]
                            else:
                                #variant not analyzed before
                                markers_to_analyze.append("V{}-{}".format(idx,item[1]))
                                pos_all.append(item[1])
                    tmp_mle_mafs=self.__getMLEfreq(data, markers_to_analyze, pos_all, fam_to_analyze[pop], self.rsq, output_log)
                    if len(tmp_mle_mafs) > 0:
                        for vname,vmaf in tmp_mle_mafs.iteritems():
                            data.mle_mafs[pop][vname]=vmaf
                            local_mle_mafs[pop][vname]=vmaf
            else:
                #Homogeneous families
                markers_to_analyze=[]
                pos_all=[]
                markers_analyzed={}
                for tmpv in data.mle_mafs:
                    markers_analyzed[tmpv.split('-')[-1]]=data.mle_mafs[tmpv]
                variants_in_fams=[]
                for item in data.families.keys():
                    var_per_fam=[tuple(tmpvar) for tmpvar in data.getFamVariants(item)]
                    variants_in_fams=list(set(var_per_fam+variants_in_fams))
                variants_in_fams=[list(tmpvar) for tmpvar in sorted(variants_in_fams, key=lambda x: x[1])]
                for item in variants_in_fams:
                    idx=data.variants.index(item)
                    if item[-1]==0 or self.af_info is None:
                        if str(item[1]) in markers_analyzed.keys():
                            #if variant has been analyzed
                            vname="V{}-{}".format(idx,item[1])
                            local_mle_mafs[vname]=markers_analyzed[str(item[1])]
                        else:
                            #variant not analyzed before
                            markers_to_analyze.append("V{}-{}".format(idx,item[1]))
                            pos_all.append(item[1])
                output_log=env.tmp_log+"AF_{}.log".format(self.name)
                tmp_mle_mafs=self.__getMLEfreq(data, markers_to_analyze, pos_all, data.families.keys(), self.rsq, output_log)
                if len(tmp_mle_mafs) > 0:
                    for vname, vmaf in tmp_mle_mafs.iteritems():
                        data.mle_mafs[vname]=vmaf
                        local_mle_mafs[vname]=vmaf
        gnomAD_pop=None
        for item in data.families:
            varnames[item], positions, vcf_mafs = data.getFamVariants(item, style = "map")
            recombPos[item]={}
            var_for_haplotype=[]
            positions_for_haplotype=[]
            output_sample=[]
            if env.debug:
                with env.lock:
                    sys.stderr.write('\n'+repr(varnames[item])+'\n')
                    sys.stderr.write('\n'.join(['\t'.join(x) for x in data.getFamSamples(item)]) + '\n\n')
            # either use privided MAF or compute MAF
            if self.freq_by_fam:
                mafs[item]={}
                tfreq_fam=data.freq_by_fam[item]
                for pop in data.gnomAD_estimate.keys():
                    if pop in tfreq_fam:
                        gnomAD_pop=pop
                        break
            elif gnomAD_pop is None and data.freq is not None:
                for pop in data.gnomAD_estimate.keys():
                    if pop in data.freq:
                        gnomAD_pop=pop
                        break
            for idx, v in enumerate(varnames[item]):
                tmp_maf_var=0
                if self.af_info is None:
                #no vcf freq column specified
                    if v not in tmp_mafs:
                        if self.mle:
                        #use MLE freq for all variants
                            tmp_mafs[v]=local_mle_mafs[v]
                        elif self.count:
                        #estimate MAF based on founder counts if MLE not specified
                            tmp_mafs[v]=local_count_mafs[v]
                        tmp_maf_var=tmp_mafs[v]
                elif not self.af_info is None:
                    #if vcf freq column is specified
                    #use vcf_mafs if possible
                    if vcf_mafs[idx]:
                        tmp_maf_var=vcf_mafs[idx]
                        if self.freq_by_fam:
                            mafs[item][v] = vcf_mafs[idx]
                        else:
                            if v not in tmp_mafs:
                                tmp_mafs[v] = vcf_mafs[idx]
                    else:
                        #if variants do not have valid vcf_mafs values if specified
                        if self.freq_by_fam:
                            if gnomAD_pop is not None:
                                mafs[item][v]=data.gnomAD_estimate[gnomAD_pop]
                            elif self.mle:
                                    mafs[item][v]=local_mle_mafs[data.freq_by_fam[item]][v]
                            elif self.count:
                                    mafs[item][v]=local_count_mafs[data.freq_by_fam[item]][v]
                            tmp_maf_var=mafs[item][v]
                        else:
                            if v not in tmp_mafs:
                                if gnomAD_pop is not None:
                                    tmp_mafs[v]=data.gnomAD_estimate[gnomAD_pop]
                                elif self.mle:
                                    tmp_mafs[v]=local_mle_mafs[v]
                                elif self.count:
                                    tmp_mafs[v]=local_count_mafs[v]
                            tmp_maf_var=tmp_mafs[v]
                if self.rvhaplo:
                    if tmp_maf_var<=self.maf_cutoff:
                        var_for_haplotype.append(v)
                        positions_for_haplotype.append(positions[idx])
            if not self.rvhaplo:
                var_for_haplotype=varnames[item]
                positions_for_haplotype=positions
            #collect sample+genotypes
            for person in data.tfam.sort_family(item):
                output_sample.append([])
                last_ele=len(output_sample)-1
                output_sample[last_ele] = data.tfam.samples[person][:-1]
                if person in data.samples:
                    for marker in var_for_haplotype:
                        idx=int(marker.split('-')[0][1:])
                        output_sample[last_ele].append(data.genotype_all[person][idx])
                else:
                    output_sample[last_ele].extend(["00"] * len(var_for_haplotype))
            # haplotyping
            if len(var_for_haplotype)==0:
                varnames.pop(item,None)
                #for person in data.families[item]:
                #    data[person] = self.missings
                continue
            for person in output_sample:
                if set(person[5:])==set(['00']):
                    data.missing_persons.append(person[1])
            with env.lock:
                if not env.prephased:
                    tmp_log_output=env.tmp_log + str(os.getpid())
                    with stdoutRedirect(to = tmp_log_output + '.log'):
                        haplotypes[item] = self.haplotyper.Execute(data.chrom, var_for_haplotype, positions_for_haplotype, output_sample, self.rsq, tmp_log_output)[0]
                else:
                    haplotypes[item] = self.__PedToHaplotype(data.getFamSamples(item))
            if len(haplotypes[item]) == 0:
                # C++ haplotyping implementation failed
                with env.chperror_counter.get_lock():
                    env.chperror_counter.value += 1
            varnames[item]=var_for_haplotype
        for item in haplotypes:
            for hap_idx,haploid in enumerate(haplotypes[item]):
                for vidx,var in enumerate(haploid[2:]):
                    if not var.endswith(':') and not var.endswith('|') and vidx!=0:
                        postvar_name=varnames[item][vidx]
                        prevar_name=varnames[item][vidx-1]
                        recomb_pair = (prevar_name,postvar_name)
                        try:
                            recombPos[item][recomb_pair].append(hap_idx)
                        except:
                            recombPos[item][recomb_pair]=[hap_idx]
        #
        # Compute founder MAFs
        #
        if len(tmp_mafs) > 0:
            if self.freq_by_fam:
                for pop in tmp_mafs:
                    for v in tmp_mafs[pop]:
                        if type(tmp_mafs[pop][v]) is list:
                            tmp_mafs[pop][v] = tmp_mafs[pop][v][0]/tmp_mafs[pop][v][1] if tmp_mafs[pop][v][1] >0 else 0.0
            else:
                for v in tmp_mafs:
                    if type(tmp_mafs[v]) is list:
                        tmp_mafs[v] = tmp_mafs[v][0]/tmp_mafs[v][1] if tmp_mafs[v][1] > 0 else 0.0
        ## Make mafs consistent in structure regardless of freq_by_fam
        if self.freq_by_fam:
            for item in haplotypes:
                popname=data.freq_by_fam[item]
                if popname not in tmp_mafs:
                    continue
                if item not in mafs:
                    mafs[item]=tmp_mafs[popname]
                else:
                    for v in tmp_mafs[popname]:
                        if v not in mafs[item]:
                            mafs[item][v]=tmp_mafs[popname][v]
        else:
            for item in haplotypes:
                mafs[item]=tmp_mafs
        if env.debug:
            with env.lock:
                print("variant mafs = ", mafs, "\n", file = sys.stderr)
        ##
        #
        # Drop some variants if maf is greater than given threshold
        #
        if not self.maf_cutoff is None or self.single_markers:
            if self.freq_by_fam:
                exclude_vars=[[] for x in range(len(data.freq))]
            for i in haplotypes.keys():
                if self.freq_by_fam:
                    pop_idx=data.freq.index(data.freq_by_fam[i])
                    tmp_exclude_vars=exclude_vars[pop_idx]
                else:
                    tmp_exclude_vars=exclude_vars
                for v in mafs[i].keys():
                    if not self.maf_cutoff is None:
                        if mafs[i][v] > self.maf_cutoff and v not in tmp_exclude_vars or v.split('-')[-1] not in data.include_vars:
                            tmp_exclude_vars.append(v)
                    if self.single_markers:
                        if v.split('-')[-1] not in data.include_vars:
                            tmp_exclude_vars.append(v)
                haplotypes[i] = listit(haplotypes[i])
                tmp_remain_vars=[x for x in varnames[i] if x not in tmp_exclude_vars]
                recomb_remain_vars=[]
                if len(tmp_remain_vars) == 0:
                    recombPos[i]={}
                else:
                    if len(recombPos[i]) > 0:
                        #extend recombination signal to neighbouring RVs
                        #if the original variant is to be excluded
                        #Only allow a maximum of one recombination event between one pair of consecutive markers
                        for pair in recombPos[i].keys():
                            if pair[1] not in tmp_exclude_vars:
                                if tmp_remain_vars.index(pair[1])!=0 and pair[1] not in recomb_remain_vars:
                                    recomb_remain_vars.append(pair[1])
                                else:
                                    del recombPos[i][pair]
                            else:
                                if varnames[i].index(pair[1]) > varnames[i].index(tmp_remain_vars[-1]):
                                    #last variant
                                    del recombPos[i][pair]
                                    continue
                                for tmp_idx in range(varnames[i].index(pair[1])+1,len(varnames[i])):
                                    if varnames[i][tmp_idx] not in tmp_exclude_vars:
                                        if tmp_remain_vars.index(varnames[i][tmp_idx])==0:
                                            #delete recombination pair if the recombination was marked to the first remaining variant
                                            del recombPos[i][pair]
                                            break
                                        for tmp_hap in recombPos[i][pair]:
                                            tmp_var=haplotypes[i][tmp_hap][tmp_idx+2]
                                            if tmp_var.endswith(':') or tmp_var.endswith('|'):
                                                haplotypes[i][tmp_hap][tmp_idx+2]=tmp_var[:-1]+'/'
                                        if varnames[i][tmp_idx] not in recomb_remain_vars:
                                            recomb_remain_vars.append(varnames[i][tmp_idx])
                                        else:
                                            del recombPos[i][pair]
                                        break
                for j in range(len(haplotypes[i])):
                    haplotypes[i][j] = haplotypes[i][j][:2] + \
                      [x for idx, x in enumerate(haplotypes[i][j][2:]) if varnames[i][idx] not in tmp_exclude_vars]
                for tmp_var in varnames[i]:
                    if tmp_var not in uniq_vars:
                             uniq_vars.append(tmp_var)
                varnames[i] = [x for x in varnames[i] if x not in tmp_exclude_vars]
                # handle trivial data
                if len(varnames[i]) == 0:
                    del varnames[i]
                    del haplotypes[i]
                if len(recombPos[i].keys())>self.recomb_max:
                    #treat as missing if recombination events occurred more than speicified times
                    for person in data.families[i]:
                        data[person] = self.missings
                    del varnames[i]
                    del haplotypes[i]
            # count how many variants are removed
            with env.commonvar_counter.get_lock():
                if self.freq_by_fam:
                    tmp_ex_vars=[tmp_var for tmp_vars in exclude_vars for tmp_var in tmp_vars]
                    env.commonvar_counter.value += len(set(tmp_ex_vars))
                else:
                    env.commonvar_counter.value += len(exclude_vars)
            # get total observed variants
            if self.freq_by_fam:
                for item in varnames:
                    pop=data.freq_by_fam[item]
                    if pop not in data.total_mafs:
                        data.total_mafs[pop]={}
                        data.total_varnames[pop]=[]
                    for v in varnames[item]:
                        if v not in data.total_mafs[pop]:
                            data.total_varnames[pop].append(v)
                            data.total_mafs[pop][v]=mafs[item][v]
                for pop in data.total_varnames:
                    data.total_varnames[pop]=sorted(data.total_varnames[pop], key=lambda x: int(x.split("-")[0][1:]))
                    data.wt_maf[pop]=1.0
                    for v,tmaf in data.total_mafs[pop].iteritems():
                        data.wt_maf[pop]*=(1-tmaf)
            else:
                data.total_varnames['pop']=[]
                for item in varnames:
                    for v in varnames[item]:
                        if v not in data.total_mafs:
                            data.total_varnames['pop'].append(v)
                            data.total_mafs[v]=mafs[item][v]
                data.wt_maf['pop']=1.0
                for v,tmaf in data.total_mafs.iteritems():
                    data.wt_maf['pop']*=(1-tmaf)
                data.total_varnames['pop']=sorted(data.total_varnames['pop'], key=lambda x: int(x.split("-")[0][1:]))

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
                gt = [hap[2 + varnames[item].index(v)] if v in varnames[item] else '?' for v in markers]
                founder_haplotypes.append(("{}-{}".format(hap[1], ihap % 2), "".join([x[1] if x[0].isupper() else x[0] for x in gt])))
        # calculate LD blocks, use r2 measure
        ld = Align.create(founder_haplotypes).matrixLD(validCharacters="12")["r2"]
        blocks = []
        for j in ld:
            block = [j]
            for k in ld[j]:
                if ld[j][k] > self.r2:
                    block.append(k)
            if len(block) > 1:
                blocks.append(block)
        # get LD clusters
        clusters = [[markers[idx] for idx in item] for item in list(connected_components(blocks))]
        if env.debug:
            with env.lock:
                print("LD blocks: ", blocks, file = sys.stderr)
                print("LD clusters: ", clusters, file = sys.stderr)
        return clusters


    def __CodeHaplotypes(self, data, haplotypes, mafs, varnames, clusters):
        # apply CHP coding
        ##Add non variant haplotype if not present in a family
        for item in data.famvaridx:
            if item not in haplotypes and data[data.families[item][0]] != ('0','0'):
                if self.freq_by_fam:
                    pop=data.freq_by_fam[item]
                    try:
                        varnames[item]=data.total_varnames[pop]
                        mafs[item]=data.total_mafs[pop]
                    except:
                        continue
                else:
                    varnames[item]=data.total_varnames['pop']
                    mafs[item]=data.total_mafs
                haplotypes[item]=[]
                for person in data.families[item]:
                    tmp_person=[item, person]
                    if '00' in data[person]:
                        tmp_person+=['?:']*len(varnames[item])
                    else:
                        tmp_person+=['1:']*len(varnames[item])
                    haplotypes[item].append(tmp_person)
                    haplotypes[item].append(tmp_person)
            elif item in haplotypes:
                nonvar_hap_flag=False
                for hap in haplotypes[item]:
                    tmp_genes=[]
                    for tmpa in hap[2:]:
                        if 'A' in tmpa or 'B' in tmpa:
                            tmp_genes.append(tmpa[1])
                        else:
                            tmp_genes.append(tmpa[0])
                    if set(tmp_genes)==set(['1']):
                        #non variant haplotype
                        nonvar_hap_flag=True
                        break
                if not nonvar_hap_flag:
                    #if family don't have non variant haplotype
                    var_num=len(varnames[item])
                    fake_person=[item, 'FAKEPERSON']+['1:']*var_num
                    haplotypes[item].append(fake_person)
                for hidx,hap in enumerate(haplotypes[item]):
                    if hap[1] in data.missing_persons:
                        missing_person=[item,hap[1]]+['?:']*len(varnames[item])
                        haplotypes[item][hidx]=missing_person

        if not clusters is None:
            clusters_idx = [[[varnames[item].index(x) for x in y] for y in clusters] for item in haplotypes]
        else:
            clusters_idx = [[[]] for item in haplotypes]
        if env.debug:
            for item in haplotypes:
                with env.lock:
                    print(varnames[item],file=sys.stderr)
                    print("hap{0}\t{1}\n".format(item,haplotypes[item]),file=sys.stderr)
        self.coder.Execute(haplotypes.values(), [[mafs[item][v] for v in varnames[item]] for item in haplotypes], clusters_idx)
        if env.debug:
            with env.lock:
                if clusters:
                    print("Family LD clusters: ", clusters_idx, "\n", file = sys.stderr)
                self.coder.Print()
        # line: [fid, sid, hap1, hap2]
        for line in self.coder.GetHaplotypes():
            if not line[1] in data:
                # this sample is not in VCF file. Every variant site should be missing
                # they have to be skipped for now
                continue
            data[line[1]] = (line[2].split(','), line[4].split(','))
            superMarkerCount=len(data[line[1]][0])
            if line[0] not in data.patterns:
                data.patterns[line[0]]=[[] for x in range(superMarkerCount)]
            for t_Marker in range(superMarkerCount):
                t_pat1=line[3].split(',')[t_Marker]
                t_pat2=line[5].split(',')[t_Marker]
                if t_pat1 not in data.patterns[line[0]][t_Marker]:
                    data.patterns[line[0]][t_Marker].append(t_pat1)
                if t_pat2 not in data.patterns[line[0]][t_Marker]:
                    data.patterns[line[0]][t_Marker].append(t_pat2)
            if len(data[line[1]][0]) > data.superMarkerCount:
                data.superMarkerCount = len(data[line[1]][0])
        # get MAF
        for item in data.famvaridx:
            if item not in haplotypes:
                for person in data.families[item]:
                    data[person]=('0','0')*data.superMarkerCount
        for item in haplotypes:
            data.maf[item] = self.coder.GetAlleleFrequencies(item)
            if not len(data.maf[item][0]):
                continue
            data.varnames_by_fam[item]=varnames[item]
            wt_maf=0
            if self.freq_by_fam:
                try:
                    wt_maf=data.wt_maf[data.freq_by_fam[item]]
                except:
                    pass
            else:
                wt_maf=data.wt_maf['pop']
            tmp_data_maf=[]
            for v in data.maf[item]:
                if len(v)==1:
                    tmp_data_maf.append((v[0],1-v[0]))
                else:
                    if np.sum(v)<1:
                        tmp_ratio=sum(v[1:])/(1-wt_maf)
                        tmp_list=[wt_maf]
                        if tmp_ratio==0:
                            tmp_list.append(1-wt_maf)
                        else:
                            for tmpv in v[1:]:
                                tmp_list.append(tmpv/tmp_ratio)
                        tmp_data_maf.append(tuple(tmp_list))
                    else:
                        tmp_data_maf.append(v)
            data.maf[item]=tuple(tmp_data_maf)
        if env.debug:
            with env.lock:
                print("marker freqs = ", data.maf, "\n", file = sys.stderr)


    def __AssignSNVHaplotypes(self, data, haplotypes, mafs, varnames):
        for item in haplotypes:
            # each person's haplotype
            data.varnames_by_fam[item]=varnames[item]
            token = ''
            for idx,line in enumerate(haplotypes[item]):
                if line[1] in data.missing_persons:
                    data[line[1]]=('0','0')
                else:
                    if not idx % 2:
                        token = line[2][1] if line[2][0].isupper() else line[2][0]
                        if token=='?':
                            token='0'
                    else:
                        tmp_token = line[2][1] if line[2][0].isupper() else line[2][0]
                        if tmp_token=='?':
                            tmp_token='0'
                        data[line[1]] = (token, tmp_token)

            # get MAF
            data.maf[item] = [(1 - mafs[item][varnames[item][0]], mafs[item][varnames[item][0]])]
            data.maf[item] = tuple(tuple(np.array(v) / np.sum(v)) if np.sum(v) else v
                              for v in data.maf[item])
        for item in data.famvaridx:
            if item not in haplotypes and data[data.families[item][0]] != ('0','0'):
                for person in data.families[item]:
                    if '00' in data[person]:
                        data[person]=('0','0')
                    else:
                        data[person]=('1','1')
                t_maf=0
                if self.freq_by_fam:
                    try:
                        t_maf=data.wt_maf[data.freq_by_fam[item]]
                    except:
                        for person in data.families[item]:
                            data[person]=('0','0')
                else:
                    t_maf=data.wt_maf['pop']
                data.maf[item]=((t_maf,1-t_maf),)
        if env.debug:
            with env.lock:
                print("marker freqs = ", data.maf, "\n", file = sys.stderr)


    def __FormatHaplotypes(self, data,recombPos,varnames,uniq_vars):
        # Reformat sample genotypes
        ## Linhai Edit: Reformat to deal with recombination events in families
        tmp_combined_recombPos={}
        sorted_var = sorted(uniq_vars, key=lambda x: int(x.split('-')[0][1:]))
        for fam in data.maf.keys():
            if len(data.maf[fam])>1:
                for pair in sorted(recombPos[fam].keys(), key=lambda x:(sorted_var.index(x[0]),sorted_var.index(x[1]))):
                    if pair[1] == varnames[fam][0]:
                        ##remove recombination event if occurred at 1st RV
                        del recombPos[fam][pair]
                        continue
                    if fam not in tmp_combined_recombPos:
                        tmp_combined_recombPos[fam]=[pair]
                    else:
                            tmp_combined_recombPos[fam].append(pair)
        tmp_all_recombs=[pair for pairs in tmp_combined_recombPos.values() for pair in pairs]
        sorted_combined_recombPos=sorted(list(set(tmp_all_recombs)),key=lambda x:(sorted_var.index(x[0]),sorted_var.index(x[1])))
        recomb_fams=tmp_combined_recombPos.keys()
        ##get sub-regions that applies to all families
        for varidx,variant in enumerate(sorted_var):
            included_fams=len(recomb_fams)
            for recomb_region in sorted_combined_recombPos:
                if varidx > sorted_var.index(recomb_region[0]) and varidx < sorted_var.index(recomb_region[1]):
                    ##if the variant is in a recombination region
                    included_fams-=1
            if included_fams==len(recomb_fams):
                if data.combined_regions==[]:
                    data.combined_regions.append([variant])
                else:
                    if sorted_var.index(data.combined_regions[-1][-1])==varidx-1:
                        neighbour_recomb_flag=False
                        for recomb_region in sorted_combined_recombPos:
                            recomb_idx=sorted_var.index(recomb_region[1])
                            if recomb_idx==varidx:
                                neighbour_recomb_flag=True
                                break
                            elif recomb_idx>varidx:
                                break
                        if neighbour_recomb_flag:
                            data.combined_regions.append([variant])
                        else:
                            data.combined_regions[-1].append(variant)
                    else:
                        data.combined_regions.append([variant])
        ##Get the markers in families compliant with the sub_regions
        for sub_region in data.combined_regions:
            markers={}
            for fam in recomb_fams:
                pidx=0
                for pair in sorted(recombPos[fam].keys(), key=lambda x:(sorted_var.index(x[0]),sorted_var.index(x[1]))):
                    sub_region_start=sorted_var.index(sub_region[0])
                    sub_region_end=sorted_var.index(sub_region[-1])
                    recomb_start=sorted_var.index(pair[0])
                    recomb_end=sorted_var.index(pair[1])
                    if sub_region_end <= recomb_start:
                        markers[fam]=pidx
                        break
                    elif sub_region_end > recomb_start and sub_region_start>recomb_start and sub_region_end<recomb_end:
                        ##within the recombination region
                        markers[fam]=None
                        break
                    pidx+=1
                if fam not in markers:
                    markers[fam]=pidx
            data.complied_markers.append(markers)
        data.superMarkerCount=len(data.combined_regions)
        for person in data:
            if type(data[person]) is not tuple:
                data[person] = self.missings
                continue
            diff = data.superMarkerCount - len(data[person][0])
            data[person] = zip(*data[person])
            if diff > 0:
                if len(data[person]) == 1:
                    ##only one whole region with no recombination
                    data[person].extend(data[person] * diff)
                else:
                    famid=''
                    for fam in data.complied_markers[0].keys():
                        if person in data.families[fam]:
                            famid=fam
                    complied_data=[]
                    for marker in data.complied_markers:
                        complied_data.append(data[person][marker[famid]])
                    data[person]=complied_data

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

    def getRegion(self, region):
        self.name = region[3]

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
            gs = [data[s][0] for s in data.samples]
            if len(set(gs)) == 1:
                # everyone's genotype is the same (most likely missing or monomorphic)
                return 2
            self.tped += env.delimiter.join([self.chrom, self.name, self.distance, position] + \
                list(itertools.chain(*gs)) + self.missings*self.num_missing) + "\n"
            # freqs
            for k in data.maf:
                self.freq += env.delimiter.join([k, self.name] + map(str, data.maf[k][0])) + "\n"
        else:
            # have to expand each region into mutiple chunks to account for different recomb points
            gs = zip(*[data[s] for s in data.samples])
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
                    if len(data.maf[k])>1:
                        matched_idx=data.complied_markers[idx][k]
                        cid += 1
                        self.freq += env.delimiter.join([k, '{}[{}]'.format(self.name, cid)] + \
                                                map(str, data.maf[k][matched_idx])) + "\n"
                    elif len(data.maf[k])==1:
                        cid += 1
                        self.freq += env.delimiter.join([k, '{}[{}]'.format(self.name, cid)] + \
                                                map(str, data.maf[k][0])) + "\n"
        self.chp += "CHP Super Marker positions: "+repr(data.combined_regions)+"\n"
        for item in data.varnames_by_fam:
            try:
                pattern_txt=[tuple(sorted(data.patterns[item][tmarker],key=lambda x:x.count('2') )) for tmarker in range(len(data.patterns[item]))]
            except:
                pattern_txt=''
            self.varfam += "{}\t{}\t{}\n".format(item,data.varnames_by_fam[item],pattern_txt)
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
        if self.chp:
            with env.lock:
                with open(os.path.join(env.tmp_cache, '{}.chr{}.chp'.format(env.output, self.prev_chrom)),
                          'a') as f:
                    f.write(self.chp)
        if self.varfam:
            with env.lock:
                with open(os.path.join(env.tmp_cache, '{}.chr{}.var'.format(env.output, self.prev_chrom)),
                          'a') as f:
                    f.write(self.varfam)
        self.reset()

    def reset(self):
        self.tped = ''
        self.freq = ''
        self.chp = ''
        self.varfam = ''
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


def main(args):
    '''the main encoder function'''
    checkParams(args)
    download_dir = 'http://bioinformatics.org/spower/download/.private'
    downloadResources([('{}/genemap.{}.txt'.format(download_dir, args.build), env.resource_dir),
                       ('{}/{}/mlink'.format(download_dir, platform.system().lower()), env.resource_bin),
                       ('{}/{}/unknown'.format(download_dir, platform.system().lower()), env.resource_bin),
                       ('{}/{}/makeped'.format(download_dir, platform.system().lower()), env.resource_bin),
                       ('{}/{}/pedcheck'.format(download_dir, platform.system().lower()), env.resource_bin)])
    if args.no_save:
        cache = NoCache()
    else:
        cache = Cache(env.cache_dir, env.output, vars(args))
    cache.setID('vcf')
    # STEP 1: write encoded data to TPED format
    if not args.vanilla and cache.check():
        env.log('Loading regional marker data from archive ...')
        cache.load(target_dir = env.tmp_dir, names = ['CACHE'])
        env.success_counter.value = sum(map(fileLinesCount, glob.glob('{}/*.tped'.format(env.tmp_cache))))
        env.batch = 10
    else:
        # load VCF file header
        checkVCFBundle(args.vcf)
        cache.clear()
        try:
            vs = cstatgen.VCFstream(args.vcf)
        except Exception as e:
            env.error("{}".format(e), exit = True)
        samples_vcf = vs.GetSampleNames()
        if len(samples_vcf) == 0:
            env.error("Fail to extract samples from [{}]".format(args.vcf), exit = True)
        env.log('{:,d} samples found in [{}]'.format(len(samples_vcf), args.vcf))
        samples_not_vcf = checkSamples(samples_vcf, getColumn(args.tfam, 2))[1]
        # load sample info
        data = RData(samples_vcf, TFAMParser(args.tfam))
        if len(data.families) == 0:
            env.error('No valid family to process. ' \
                      'Families have to be at least trio with at least one member in VCF file.', exit = True)
        if len(data.samples) == 0:
            env.error('No valid sample to process. ' \
                      'Samples have to be in families, and present in both TFAM and VCF files.', exit = True)
        rewriteFamfile(os.path.join(env.tmp_cache, '{}.tfam'.format(env.output)),
                       data.tfam.samples, data.samples.keys() + samples_not_vcf)
        if args.single_markers:
            regions=[]
            for x in vs.GetGenomeCoordinates():
                region_info = (x[0], x[1], x[1], "{}:{}".format(x[0], x[1]), '.', '.', '.')
                if region_info not in regions:
                    regions.append(region_info)
            args.blueprint = None
        else:
            # load blueprint
            try:
                with open(args.blueprint, 'r') as f:
                    regions = [x.strip().split() for x in f.readlines()]
            except IOError:
                env.error("Cannot load regional marker blueprint [{}]. ".format(args.blueprint), exit = True)
        env.log('{:,d} families with a total of {:,d} samples will be scanned for {:,d} pre-defined units'.\
                format(len(data.families), len(data.samples), len(regions)))
        env.jobs = max(min(args.jobs, len(regions)), 1)
        regions.extend([None] * env.jobs)
        queue = [] if env.jobs == 1 else Queue()
        try:
            faulthandler.enable(file=open(env.tmp_log + '.SEGV', 'w'))
            for i in regions:
                if isinstance(queue, list):
                    queue.append(i)
                else:
                    queue.put(i)
            freq_by_fam_flag = False
            if not args.freq_by_fam is None:
                freq_by_fam_flag = True
                with open(args.freq_by_fam) as freq_fh:
                    for freq_line in freq_fh:
                        tmp_eles=freq_line.split()   #Fam and Population
                        data.freq_by_fam[tmp_eles[0]]=tmp_eles[1]
                data.freq=sorted(list(set(data.freq_by_fam.values())))
            else:
                data.freq=args.freq
            jobs = [EncoderWorker(
                queue, len(regions), deepcopy(data),
                RegionExtractor(args.vcf, chr_prefix = args.chr_prefix, allele_freq_info = data.freq, include_vars_file=args.include_vars),
                MarkerMaker(args.bin, maf_cutoff = args.maf_cutoff,single_markers=args.single_markers,recomb_max=args.recomb_max,af_info=data.freq,freq_by_fam=freq_by_fam_flag,rsq=args.rsq,mle=args.mle, rvhaplo=args.rvhaplo),
                LinkageWriter(len(samples_not_vcf))
                ) for i in range(env.jobs)]
            for j in jobs:
                j.start()
            for j in jobs:
                j.join()
            faulthandler.disable()
        except KeyboardInterrupt:
            # FIXME: need to properly close all jobs
            raise ValueError("Use 'killall {}' to properly terminate all processes!".format(env.prog))
        else:
            env.log('{:,d} units (from {:,d} variants) processed; '\
                '{:,d} Mendelian inconsistencies and {:,d} recombination events handled\n'.\
                format(env.success_counter.value,
                       env.variants_counter.value,
                       env.mendelerror_counter.value,
                       env.recomb_counter.value), flush = True)
            if env.triallelic_counter.value:
                env.log('{:,d} tri-allelic loci were ignored'.format(env.triallelic_counter.value))
            if env.commonvar_counter.value:
                env.log('{:,d} variants ignored due to having MAF > {} and other specified constraints'.\
                        format(env.commonvar_counter.value, args.maf_cutoff))
            if env.null_counter.value:
                env.log('{:,d} units ignored due to absence in VCF file'.format(env.null_counter.value))
            if env.trivial_counter.value:
                env.log('{:,d} units ignored due to absence of variation in samples'.format(env.trivial_counter.value))
            fatal_errors = 0
            try:
                # Error msg from C++ extension
                os.system("cat {}/*.* > {}".format(env.tmp_dir, env.tmp_log))
                fatal_errors = wordCount(env.tmp_log)['fatal']
            except KeyError:
                pass
            if env.chperror_counter.value:
                env.error("{:,d} regional markers failed to be generated due to haplotyping failures!".\
                          format(env.chperror_counter.value))
            if fatal_errors:
                env.error("{:,d} or more regional markers failed to be generated due to runtime errors!".\
                          format(fatal_errors))
            env.log('Archiving regional marker data to directory [{}]'.format(env.cache_dir))
            cache.write(arcroot = 'CACHE', source_dir = env.tmp_cache)
    env.jobs = args.jobs
    # STEP 2: write to PLINK or mega2 format
    tpeds = [os.path.join(env.tmp_cache, item) for item in os.listdir(env.tmp_cache) if item.startswith(env.output) and item.endswith('.tped')]
    for fmt in args.format:
        cache.setID(fmt)
        if not args.vanilla and cache.check():
            env.log('Loading {} data from archive ...'.format(fmt.upper()))
            cache.load(target_dir = env.tmp_dir, names = [fmt.upper()])
        else:
            env.log('{:,d} units will be converted to {} format'.format(env.success_counter.value, fmt.upper()))
            env.format_counter.value = 0
            format(tpeds, os.path.join(env.tmp_cache, "{}.tfam".format(env.output)),
                   args.prevalence, args.wild_pen, args.muta_pen, fmt,
                   args.inherit_mode, args.theta_max, args.theta_inc)
            env.log('{:,d} units successfully converted to {} format\n'.\
                    format(env.format_counter.value, fmt.upper()), flush = True)
            if env.skipped_counter.value:
                # FIXME: perhaps we need to rephrase this message?
                env.log('{} region - family pairs skipped'.\
                        format(env.skipped_counter.value))
            env.log('Archiving {} format to directory [{}]'.format(fmt.upper(), env.cache_dir))
            cache.write(arcroot = fmt.upper(),
                        source_dir = os.path.join(env.tmp_dir, fmt.upper()), mode = 'a')
    mkpath(env.outdir)
    if args.run_linkage:
        cache.setID('analysis')
        if not args.vanilla and cache.check():
            env.log('Loading linkage analysis result from archive ...'.format(fmt.upper()))
            cache.load(target_dir = env.output, names = ['heatmap'])
        else:
            env.log('Running linkage analysis ...'.format(fmt.upper()))
            run_linkage(args.blueprint, args.theta_inc, args.theta_max, args.output_limit)
            env.log('Linkage analysis succesfully performed for {:,d} units\n'.\
                    format(env.run_counter.value, fmt.upper()), flush = True)
            if env.makeped_counter.value:
                env.log('{} "makeped" runtime errors occurred'.format(env.makeped_counter.value))
            if env.pedcheck_counter.value:
                env.log('{} "pedcheck" runtime errors occurred'.format(env.pedcheck_counter.value))
            if env.unknown_counter.value:
                env.log('{} "unknown" runtime errors occurred'.format(env.unknown_counter.value))
            if env.mlink_counter.value:
                env.log('{} "mlink" runtime errors occurred'.format(env.mlink_counter.value))
            cache.write(arcroot = 'heatmap', source_dir = os.path.join(env.output, 'heatmap'), mode = 'a')
        html(args.theta_inc, args.theta_max, args.output_limit)
    else:
        env.log('Saving data to [{}]'.format(os.path.abspath(env.output)))
        cache.load(target_dir = env.output, names = [fmt.upper() for fmt in args.format])
