from __future__ import print_function
from SEQLinkage.Utils import *
from SEQLinkage.Runner import *
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
from concurrent.futures import ProcessPoolExecutor

from SEQLinkage.Main import *
from SEQLinkage.Core import *


def main():
    args = Args().parser.parse_args('--bin 1 --fam ../data/new_trim_ped_famless17_no:xx.fam --vcf ../data/first1000snp_full_samples.vcf.gz --anno ../data/first1000_chr1_multianno.csv --pop ../data/full_sample_fam_pop.txt -f MERLIN --blueprint ../data/genemap.hg38.txt --freq AF'.split())

    checkParams(args)

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
        data = RData(args.vcf, args.tfam,args.anno,args.pop,allele_freq_info=args.freq)
        vs = data.vs
        samples_vcf = data.samples_vcf

    if len(samples_vcf) == 0:
        env.error("Fail to extract samples from [{}]".format(args.vcf), exit = True)
    env.log('{:,d} samples found in [{}]'.format(len(samples_vcf), args.vcf))
    samples_not_vcf = data.samples_not_vcf

    if len(data.families) == 0:
        env.error('No valid family to process. ' \
                  'Families have to be at least trio with at least one member in VCF file.', exit = True)
    if len(data.samples) == 0:
        env.error('No valid sample to process. ' \
                  'Samples have to be in families, and present in both TFAM and VCF files.', exit = True)
    rewriteFamfile(os.path.join(env.tmp_cache, '{}.tfam'.format(env.output)),
                   data.tfam.samples, list(data.samples.keys()) + samples_not_vcf)

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

    extractor =RegionExtractor(args.vcf, build=env.build,chr_prefix = args.chr_prefix)
    maker =            MarkerMaker(args.bin, maf_cutoff = args.maf_cutoff,recomb=False)
    writer =             LinkageWriter(len(samples_not_vcf))
    region = regions[2]
    for _ in range(5):
        extractor.getRegion(region)
        maker.getRegion(region)
        writer.getRegion(region)
        extractor.apply(data)
        maker.apply(data)
        #test(regions[2],data,extractor,maker,writer)
    #run_each_region(regions[:10],data,extractor,maker,writer)
    
main()
    
   