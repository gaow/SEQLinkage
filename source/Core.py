#!/usr/bin/python2.7
# Copyright (c) 2013, Gao Wang <gaow@uchicago.edu>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)
from __future__ import print_function
from SEQLinkage import HOMEPAGE
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
        env.tmp_log = os.path.join(env.tmp_dir, env.output + ".STDOUT")
    #
    if len([x for x in set(getColumn(args.tfam, 6)) if x.lower() not in env.ped_missing]) > 2:
        env.trait = 'quantitative'
    env.log('{} trait detected in [{}]'.format(env.trait.capitalize(), args.tfam))
    if not args.blueprint:
        args.blueprint = os.path.join(env.resource_dir, 'genemap.txt')
    args.format = [x.lower() for x in set(args.format)]
    if args.run_linkage and "linkage" not in args.format:
        args.format.append('linkage')
    if None in [args.inherit_mode, args.prevalence, args.wild_pen, args.muta_pen] and "linkage" in args.format:
        env.error('To generate LINKAGE format or run LINKAGE analysis, please specify all options below:\n\t--prevalence, -K\n\t--moi\n\t--wild-pen, -W\n\t--muta-pen, -M', show_help = True, exit = True)
    if args.tempdir is not None:
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
        # reorder family samples based on order of VCF file
        for k in self.families.keys():
            if len(self.families[k]) == 0:
                # skip families having no samples in VCF file
                del self.families[k]
            else:
                self.famsampidx[k] = [i for i, x in enumerate(samples_vcf) if x in self.families[k]]
        # a dict of {fid:[idx ...], ...}
        self.famvaridx = {}
        self.reset()

    def reset(self):
        for item in self.samples:
            self[item] = []
        self.variants = []
        self.chrom = None
        for k in self.families:
            self.famvaridx[k] = []
        self.maf = OrderedDict()
        # superMarkerCount is the max num. of recombinant fragments among all fams
        self.superMarkerCount = 0

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
            mafs = []
            for idx in self.famvaridx[fam]:
                names.append("V{}-{}".format(idx, self.variants[idx][1]))
                pos.append(self.variants[idx][1])
                mafs.append(self.variants[idx][-1])
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
    def __init__(self, filename, build = env.build, chr_prefix = None, allele_freq_info = None):
        self.vcf = cstatgen.VCFstream(filename)
        self.chrom = self.startpos = self.endpos = self.name = None
        self.chr_prefix = chr_prefix
        # name of allele frequency meta info
        self.af_info = allele_freq_info
        self.xchecker = PseudoAutoRegion('X', build)
        self.ychecker = PseudoAutoRegion('Y', build)

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
            # check if the line's sample number matches the entire VCF sample number
            if not self.vcf.CountSampleGenotypes() == self.vcf.sampleCount:
                raise ValueError('Genotype and sample mismatch for region {}: {:,d} vs {:,d}'.\
                             format(self.name, self.vcf.CountSampleGenotypes(), self.vcf.sampleCount))
            # valid line found, get variant info
            try:
                maf = float(self.vcf.GetInfo(self.af_info)) if self.af_info else None
                if maf > 0.5:
                    maf = 1 - maf
            except Exception as e:
                raise ValueError("VCF line {}:{} does not have valid allele frequency field {}!".\
                                 format(self.vcf.GetChrom(), self.vcf.GetPosition(), self.af_info))
            data.variants.append([self.vcf.GetChrom(), self.vcf.GetPosition(), self.name, maf])
            # for each family assign member genotype if the site is non-trivial to the family
            for k in data.families:
                gs = self.vcf.GetGenotypes(data.famsampidx[k])
                if len(set(''.join([x for x in gs if x != "00"]))) <= 1:
                    # skip monomorphic gs
                    continue
                else:
                    # this variant is found in the family
                    data.famvaridx[k].append(varIdx)
                    for person, g in zip(data.families[k], gs):
                        data[person].append(g)
            varIdx += 1
        #
        if varIdx == 0:
            return 1
        else:
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
    def __init__(self, wsize, maf_cutoff = None):
        self.missings = ("0", "0")
        self.gtconv = {'1':0, '2':1}
        self.haplotyper = cstatgen.HaplotypingEngine(verbose = env.debug)
        if wsize == 0 or wsize >= 1:
            self.r2 = None
        else:
            self.r2 = wsize
        self.coder = cstatgen.HaplotypeCoder(wsize)
        self.maf_cutoff = maf_cutoff

    def apply(self, data):
        # temp raw haplotype, maf and variant names data
        haplotypes = OrderedDict()
        mafs = {}
        varnames = {}
        try:
            # haplotyping plus collect found allele counts
            # and computer founder MAFS
            self.__Haplotype(data, haplotypes, mafs, varnames)
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
            return -1
        self.__FormatHaplotypes(data)
        return 0

    def __Haplotype(self, data, haplotypes, mafs, varnames):
        '''genetic haplotyping. haplotypes stores per family data'''
        # FIXME: it is SWIG's (2.0.12) fault not to properly destroy the object "Pedigree" in "Execute()"
        # So there is a memory leak here which I tried to partially handle on C++
        #
        # Per family haplotyping
        #
        self.markers = ["V{}-{}".format(idx, item[1]) for idx, item in enumerate(data.variants)]
        for item in data.families:
            varnames[item], positions, vcf_mafs = data.getFamVariants(item, style = "map")
            if len(varnames[item]) == 0:
                for person in data.families[item]:
                    data[person] = self.missings
                continue
            if env.debug:
                with env.lock:
                    sys.stderr.write('\n'.join(['\t'.join(x) for x in data.getFamSamples(item)]) + '\n\n')
            # haplotyping
            with env.lock:
                if not env.prephased:
                    with stdoutRedirect(to = env.tmp_log + str(os.getpid()) + '.log'):
                        haplotypes[item] = self.haplotyper.Execute(data.chrom, varnames[item],
                                                               sorted(positions), data.getFamSamples(item))[0]
                else:
                    haplotypes[item] = self.__PedToHaplotype(data.getFamSamples(item))
            if len(haplotypes[item]) == 0:
                # C++ haplotyping implementation failed
                with env.chperror_counter.get_lock():
                    env.chperror_counter.value += 1
            # either use privided MAF or computer MAF
            if all(vcf_mafs):
                for idx, v in enumerate(varnames[item]):
                    if v not in mafs:
                        mafs[v] = vcf_mafs[idx]
            else:
                # count founder alleles
                for hap in haplotypes[item]:
                    if not data.tfam.is_founder(hap[1]):
                        continue
                    for idxv, v in enumerate(varnames[item]):
                        if v not in mafs:
                            # [#alt, #haplotypes]
                            mafs[v] = [0, 0]
                        gt = hap[2 + idxv][1] if hap[2 + idxv][0].isupper() else hap[2 + idxv][0]
                        if not gt == "?":
                            mafs[v][0] += self.gtconv[gt]
                            mafs[v][1] += 1.0
        #
        # Compute founder MAFs
        #
        for v in mafs:
            if type(mafs[v]) is not list:
                continue
            mafs[v] = mafs[v][0] / mafs[v][1] if mafs[v][1] > 0 else 0.0
        if env.debug:
            with env.lock:
                print("variant mafs = ", mafs, "\n", file = sys.stderr)
        #
        # Drop some variants if maf is greater than given threshold
        #
        if self.maf_cutoff is not None:
            exclude_vars = []
            for v in mafs.keys():
                if mafs[v] > self.maf_cutoff:
                    exclude_vars.append(v)
            for i in haplotypes.keys():
                haplotypes[i] = listit(haplotypes[i])
                for j in range(len(haplotypes[i])):
                    haplotypes[i][j] = haplotypes[i][j][:2] + \
                      [x for idx, x in enumerate(haplotypes[i][j][2:]) if varnames[i][idx] not in exclude_vars]
                varnames[i] = [x for x in varnames[i] if x not in exclude_vars]
                # handle trivial data
                if len(varnames[i]) == 0:
                    for person in data.families[i]:
                        data[person] = self.missings
                    del varnames[i]
                    del haplotypes[i]
            # count how many variants are removed
            with env.commonvar_counter.get_lock():
                env.commonvar_counter.value += len(exclude_vars)


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
        if clusters is not None:
            clusters_idx = [[[varnames[item].index(x) for x in y] for y in clusters] for item in haplotypes]
        else:
            clusters_idx = [[[]] for item in haplotypes]
        self.coder.Execute(haplotypes.values(), [[mafs[v] for v in varnames[item]] for item in haplotypes], clusters_idx)
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
            data[line[1]] = (line[2].split(','), line[3].split(','))
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
        for item in haplotypes:
            # each person's haplotype
            token = ''
            for idx, line in enumerate(haplotypes[item]):
                if not idx % 2:
                    token = line[2][1] if line[2][0].isupper() else line[2][0]
                else:
                    data[line[1]] = (token, line[2][1] if line[2][0].isupper() else line[2][0])
            # get maf
            data.maf[item] = [(1 - mafs[varnames[item][0]], mafs[varnames[item][0]])]
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
            data[person] = zip(*data[person])
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
                    if idx >= len(data.maf[k]):
                        break
                    cid += 1
                    self.freq += env.delimiter.join([k, '{}[{}]'.format(self.name, cid)] + \
                                                    map(str, data.maf[k][idx])) + "\n"
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
                region = self.queue.get()
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
    downloadResources([('{}/uploads/genemap.txt'.format(HOMEPAGE), env.resource_dir),
                       ('{}/uploads/{}/mlink'.format(HOMEPAGE, platform.system().lower()), env.resource_bin),
                       ('{}/uploads/{}/unknown'.format(HOMEPAGE, platform.system().lower()), env.resource_bin),
                       ('{}/uploads/{}/makeped'.format(HOMEPAGE, platform.system().lower()), env.resource_bin),
                       ('{}/uploads/{}/pedcheck'.format(HOMEPAGE, platform.system().lower()), env.resource_bin)])
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
            regions = [(x[0], x[1], x[1], "{}:{}".format(x[0], x[1]), '.', '.', '.')
                       for x in vs.GetGenomeCoordinates()]
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
        queue = Queue()
        try:
            faulthandler.enable(file=open(env.tmp_log + '.SEGV', 'w'))
            for i in regions:
                queue.put(i)
            jobs = [EncoderWorker(
                queue, len(regions), deepcopy(data),
                RegionExtractor(args.vcf, chr_prefix = args.chr_prefix, allele_freq_info = args.freq),
                MarkerMaker(args.bin, maf_cutoff = args.maf_cutoff),
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
                env.log('{:,d} variants ignored due to having MAF > {}'.\
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
