# SEQLinkage
> Collapsed Haplotype Pattern Method for Linkage Analysis of Next-Generation Sequencing Data


## Features
- It can do linkage analysis on single variants and CHP markers.
- It can analyze families from different population.
- It can handle large-scale whole-genome linkage analysis.

## Pre-requisites

Make sure you install the pre-requisited before running seqlink:


```
#install cstatgen
conda install -c conda-forge xeus-cling
conda install -c anaconda swig 
conda install -c conda-forge gsl
pip install egglib
git clone https://github.com/statgenetics/cstatgen.git
cd cstatgen
python setup.py install

#install paramlink2
R
install.packages("paramlink2")
```

## Install

`pip install SEQLinkage`

## How to use

```python
!seqlink --help
```

    usage: seqlink [-h] [--single-marker] --fam FILE --vcf FILE [--anno FILE]
                   [--pop FILE] [--included-vars FILE] [-b FILE] [-c P] [-o Name]
                   [--build STRING] [--window INT] [--freq INFO]
                   [--chrom-prefix STRING] [--run-linkage] [-K FLOAT]
                   [--moi STRING] [-W FLOAT] [-M FLOAT] [--theta-max FLOAT]
                   [--theta-inc FLOAT]
    
    SEQLinkage V2, linkage analysis using sequence data
    
    options:
      -h, --help            show this help message and exit
    
    Collapsed haplotype pattern method arguments:
      --single-marker       Use single variant as the marker. Otherwise, use CHP
                            markers.
      --fam FILE            Input pedigree and phenotype information in FAM
                            format.
      --vcf FILE            Input VCF file, bgzipped.
      --anno FILE           Input annotation file from annovar.
      --pop FILE            Input two columns file, first column is family ID,
                            second column population information.
      --included-vars FILE  Variants to be included for linkage analysis, if None,
                            the analysis won't filter any variants. But you can
                            still set AF cutoff by -c argment.
      -b FILE, --blueprint FILE
                            Blueprint file that defines regional marker (format:
                            "chr startpos endpos name avg.distance male.distance
                            female.distance").
      -c P, --maf-cutoff P  MAF cutoff to define variants to be excluded from
                            analyses. this is useful, if you need to analyse
                            multiple population together.
      -o Name, --output Name
                            Output name prefix.
      --build STRING        Reference genome version for VCF file.
      --window INT          If no blueprint, seprate chromosome to pseudogenes
                            with 1000 (as default) variants.
      --freq INFO           Info field name for allele frequency in VCF file.
      --chrom-prefix STRING
                            Prefix to chromosome name in VCF file if applicable,
                            e.g. "chr".
    
    LINKAGE options:
      --run-linkage         Perform Linkage analysis.
      -K FLOAT, --prevalence FLOAT
                            Disease prevalence. Default to 0.001.
      --moi STRING          Mode of inheritance, AD/AR: autosomal
                            dominant/recessive. Default to AD.
      -W FLOAT, --wt-pen FLOAT
                            Penetrance for wild type. Default to 0.01.
      -M FLOAT, --mut-pen FLOAT
                            Penetrance for mutation. Default to 0.9.
      --theta-max FLOAT     Theta upper bound. Default to 0.5.
      --theta-inc FLOAT     Theta increment. Default to 0.05.


### Linkage analysis on specific regions
> Normally, the regions are gene regions. you can also use self-defined regions, such as promoter regions, enhancer regions.

#### 1.run seqlink on CHP marker

```python
!seqlink --fam testdata/test_ped.fam --vcf testdata/test_snps.vcf.gz --anno testdata/test_chr1_anno.csv --pop testdata/test_fam_pop.txt --blueprint testdata/test_blueprint_ext.txt --included-vars testdata/test_chr1_included_vars.txt -o data/test_chp --run-linkage
```

    [1;40;32mMESSAGE: Binary trait detected in [/mnt/vast/hpc/csg/yin/Github/linkage/SEQpy3/nbs/testdata/test_ped.fam][0m
    [1;40;32mMESSAGE: Namespace(single_marker=False, tfam='/mnt/vast/hpc/csg/yin/Github/linkage/SEQpy3/nbs/testdata/test_ped.fam', vcf='/mnt/vast/hpc/csg/yin/Github/linkage/SEQpy3/nbs/testdata/test_snps.vcf.gz', anno='testdata/test_chr1_anno.csv', pop='testdata/test_fam_pop.txt', included_vars='testdata/test_chr1_included_vars.txt', blueprint='testdata/test_blueprint_ext.txt', maf_cutoff=None, output='data/test_chp', build='hg38', window=1000, freq='AF', chr_prefix=None, run_linkage=True, prevalence=0.001, inherit_mode='AD', wild_pen=0.01, muta_pen=0.9, theta_max=0.5, theta_inc=0.05)[0m
    [1;40;32mMESSAGE: 18 samples found in FAM file but not in VCF file:[0m
    
    [1;40;32mMESSAGE: 18 samples found in [/mnt/vast/hpc/csg/yin/Github/linkage/SEQpy3/nbs/testdata/test_snps.vcf.gz][0m
    [1;40;32mMESSAGE: Loading marker map from [testdata/test_blueprint_ext.txt] ...[0m
    [1;40;32mMESSAGE: 6 families with a total of 18 samples will be scanned for 12 pre-defined units[0m
    SNVHap MIR6859-1@1,MIR6859-2@1,MIR6859-3@1,MIR6859-4@1
    [1;40;32mMESSAGE: write to pickle: data/test_chp/chr1result/chr1result0.pickle,Gene number:2,Time:5.62837730265326e-05[0m
    create data/test_chp/chr1result/chr1result0_AFcutoffNone_linkage.input
    create data/test_chp/chr1result/chr1result0_AFcutoffNone_linkage.lods
    0.21258915215730667
    create data/test_chp/chr1result/chr1result0_AFcutoffNone_linkage.besthlod
    [1;40;32mMESSAGE: ============= Finish analysis ==============[0m


#### 2.run seqlink on variants

```python
!seqlink --single-marker --fam testdata/test_ped.fam --vcf testdata/test_snps.vcf.gz --anno testdata/test_chr1_anno.csv --pop testdata/test_fam_pop.txt --blueprint testdata/test_blueprint_ext.txt -c 0.05 -o data/test_var --run-linkage
```

    [1;40;32mMESSAGE: Binary trait detected in [/mnt/vast/hpc/csg/yin/Github/linkage/SEQpy3/nbs/testdata/test_ped.fam][0m
    [1;40;32mMESSAGE: Namespace(single_marker=True, tfam='/mnt/vast/hpc/csg/yin/Github/linkage/SEQpy3/nbs/testdata/test_ped.fam', vcf='/mnt/vast/hpc/csg/yin/Github/linkage/SEQpy3/nbs/testdata/test_snps.vcf.gz', anno='testdata/test_chr1_anno.csv', pop='testdata/test_fam_pop.txt', included_vars=None, blueprint='testdata/test_blueprint_ext.txt', maf_cutoff=0.05, output='data/test_var', build='hg38', window=1000, freq='AF', chr_prefix=None, run_linkage=True, prevalence=0.001, inherit_mode='AD', wild_pen=0.01, muta_pen=0.9, theta_max=0.5, theta_inc=0.05)[0m
    [1;40;32mMESSAGE: 18 samples found in FAM file but not in VCF file:[0m
    
    [1;40;32mMESSAGE: 18 samples found in [/mnt/vast/hpc/csg/yin/Github/linkage/SEQpy3/nbs/testdata/test_snps.vcf.gz][0m
    [1;40;32mMESSAGE: Loading marker map from [testdata/test_blueprint_ext.txt] ...[0m
    [1;40;32mMESSAGE: 6 families with a total of 18 samples will be scanned for 12 pre-defined units[0m
    [1;40;32mMESSAGE: write to pickle: data/test_var/chr1result/chr1result0.pickle,Gene number:4,Time:4.1139241204493574e-05[0m
    create data/test_var/chr1result/chr1result0_AFcutoff0.05_linkage.input
    create data/test_var/chr1result/chr1result0_AFcutoff0.05_linkage.lods
    0.3724569082260132
    create data/test_var/chr1result/chr1result0_AFcutoff0.05_linkage.besthlod
    [1;40;32mMESSAGE: ============= Finish analysis ==============[0m


### No annotation
> If you don't have the annotation file. there is no need to add `--pop`. And `--freq` should be setted based on the `INFO` column in vcf file.

```python
!seqlink --fam testdata/test_ped.fam --vcf testdata/test_snps.vcf.gz --freq='AF' --blueprint testdata/test_blueprint_ext.txt -c 0.05 -o data/test_chp_na --run-linkage
```

    [1;40;32mMESSAGE: Binary trait detected in [/mnt/vast/hpc/csg/yin/Github/linkage/SEQpy3/nbs/testdata/test_ped.fam][0m
    [1;40;32mMESSAGE: Namespace(single_marker=False, tfam='/mnt/vast/hpc/csg/yin/Github/linkage/SEQpy3/nbs/testdata/test_ped.fam', vcf='/mnt/vast/hpc/csg/yin/Github/linkage/SEQpy3/nbs/testdata/test_snps.vcf.gz', anno=None, pop=None, included_vars=None, blueprint='testdata/test_blueprint_ext.txt', maf_cutoff=0.05, output='data/test_chp_na', build='hg38', window=1000, freq='AF', chr_prefix=None, run_linkage=True, prevalence=0.001, inherit_mode='AD', wild_pen=0.01, muta_pen=0.9, theta_max=0.5, theta_inc=0.05)[0m
    [1;40;32mMESSAGE: 18 samples found in FAM file but not in VCF file:[0m
    
    [1;40;32mMESSAGE: 18 samples found in [/mnt/vast/hpc/csg/yin/Github/linkage/SEQpy3/nbs/testdata/test_snps.vcf.gz][0m
    [1;40;32mMESSAGE: Loading marker map from [testdata/test_blueprint_ext.txt] ...[0m
    [1;40;32mMESSAGE: 6 families with a total of 18 samples will be scanned for 12 pre-defined units[0m
    SNVHap MIR6859-1@1,MIR6859-2@1,MIR6859-3@1,MIR6859-4@1
    [1;40;32mMESSAGE: write to pickle: data/test_chp_na/chrallresult/chrallresult0.pickle,Gene number:4,Time:9.55304606921143e-05[0m
    create data/test_chp_na/chrallresult/chrallresult0_AFcutoff0.05_linkage.input
    create data/test_chp_na/chrallresult/chrallresult0_AFcutoff0.05_linkage.lods
    0.3595982789993286
    create data/test_chp_na/chrallresult/chrallresult0_AFcutoff0.05_linkage.besthlod
    [1;40;32mMESSAGE: ============= Finish analysis ==============[0m


### Whole-genome linkage analysis
> if `--blueprint` is not provided, the genomic region will be seperated to pseudogenes with 1000 variants. you can change the variant number per pseudogene by `--window`.

```python
!seqlink --single-marker --fam testdata/test_ped.fam --vcf testdata/test_snps.vcf.gz --anno testdata/test_chr1_anno.csv --pop testdata/test_fam_pop.txt -c 0.05 -o data/test_wg --run-linkage
```

    [1;40;32mMESSAGE: Binary trait detected in [/mnt/vast/hpc/csg/yin/Github/linkage/SEQpy3/nbs/testdata/test_ped.fam][0m
    [1;40;32mMESSAGE: Generate regions by annotation[0m
    [1;40;32mMESSAGE: Namespace(single_marker=True, tfam='/mnt/vast/hpc/csg/yin/Github/linkage/SEQpy3/nbs/testdata/test_ped.fam', vcf='/mnt/vast/hpc/csg/yin/Github/linkage/SEQpy3/nbs/testdata/test_snps.vcf.gz', anno='testdata/test_chr1_anno.csv', pop='testdata/test_fam_pop.txt', included_vars=None, blueprint=None, maf_cutoff=0.05, output='data/test_wg', build='hg38', window=1000, freq='AF', chr_prefix=None, run_linkage=True, prevalence=0.001, inherit_mode='AD', wild_pen=0.01, muta_pen=0.9, theta_max=0.5, theta_inc=0.05)[0m
    [1;40;32mMESSAGE: 18 samples found in FAM file but not in VCF file:[0m
    
    [1;40;32mMESSAGE: 18 samples found in [/mnt/vast/hpc/csg/yin/Github/linkage/SEQpy3/nbs/testdata/test_snps.vcf.gz][0m
    [1;40;32mMESSAGE: separate chromosome to regions[0m
    [1;40;32mMESSAGE: 6 families with a total of 18 samples will be scanned for 1 pre-defined units[0m
    [1;40;32mMESSAGE: write to pickle: data/test_wg/chr1result/chr1result0.pickle,Gene number:1,Time:9.195781416363186e-05[0m
    create data/test_wg/chr1result/chr1result0_AFcutoff0.05_linkage.input
    create data/test_wg/chr1result/chr1result0_AFcutoff0.05_linkage.lods
    0.7846571207046509
    create data/test_wg/chr1result/chr1result0_AFcutoff0.05_linkage.besthlod
    [1;40;32mMESSAGE: ============= Finish analysis ==============[0m


## Input format

- `--fam`, Fam file (required, format: "fid iid fathid mothid sex trait[1 control, 2 case, -9 or 0 missing]")

```python
%%writefile testdata/test_ped.fam
1033    1033_2  0       0       2       -9
1033    1033_1  0       0       1       -9
1033    1033_99 1033_1  1033_2  2       1
1033    1033_9  1033_1  1033_2  2       1
1033    1033_3  1033_1  1033_2  2       2
1036    1036_99 1036_1  1036_2  2       2
1036    1036_6  0       0       1       2
1036    1036_1  0       0       1       -9
1036    1036_3  1036_6  1036_99 2       1
1036    1036_4  1036_6  1036_99 2       1
1036    1036_2  0       0       2       -9
1036    1036_5  1036_6  1036_99 1       1
10J_103 10J_103_10      0       0       1       -9
10J_103 10J_103_4       0       0       1       -9
10J_103 10J_103_3       0       0       2       -9
10J_103 10J_103_2       10J_103_4       10J_103_3       2       2
10J_103 10J_103_1       10J_103_10      10J_103_3       1       2
10J_109 10J_109_2       10J_109_6       10J_109_5       1       2
10J_109 10J_109_3       10J_109_6       10J_109_5       1       2
10J_109 10J_109_4       10J_109_6       10J_109_5       1       2
10J_109 10J_109_6       0       0       1       -9
10J_109 10J_109_1       10J_109_6       10J_109_5       1       2
10J_109 10J_109_5       0       0       2       2
10J_109 10J_109_7       10J_109_6       10J_109_5       1       1
10J_112 10J_112_3       0       0       1       1
10J_112 10J_112_5       10J_112_3       10J_112_2       1       2
10J_112 10J_112_1       10J_112_3       10J_112_2       2       1
10J_112 10J_112_7       10J_112_3       10J_112_2       1       1
10J_112 10J_112_2       0       0       2       2
10J_119 10J_119_2       0       0       1       1
10J_119 10J_119_5       0       0       2       1
10J_119 10J_119_4       0       0       1       1
10J_119 10J_119_6       10J_119_4       10J_119_5       1       2
10J_119 10J_119_7       10J_119_4       10J_119_5       2       2
10J_119 10J_119_1       10J_119_4       10J_119_5       2       2
10J_119 10J_119_3       10J_119_2       10J_119_1       1       1
```

    Overwriting ../testdata/test_ped.fam


- `--vcf`, VCF file (required, vcf.gz with vcf.gz.tbi)
```
bgzip -c file.vcf > file.vcf.gz
tabix -p vcf file.vcf.gz
```


- `--anno`, Annotation file from `ANNOVAR`, It must contains the allele frequency for population in the file of family population information. For example in here, The annotation file must have AF_amr, AF_afr, AF_nfe columns.

```python
anno=pd.read_csv('testdata/test_chr1_anno.csv')
```

```python
anno
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Chr</th>
      <th>Start</th>
      <th>End</th>
      <th>Ref</th>
      <th>Alt</th>
      <th>Func.refGene</th>
      <th>Gene.refGene</th>
      <th>GeneDetail.refGene</th>
      <th>ExonicFunc.refGene</th>
      <th>AAChange.refGene</th>
      <th>...</th>
      <th>CLNDISDB</th>
      <th>CLNREVSTAT</th>
      <th>CLNSIG</th>
      <th>DN ID</th>
      <th>Patient ID</th>
      <th>Phenotype</th>
      <th>Platform</th>
      <th>Study</th>
      <th>Pubmed ID</th>
      <th>Otherinfo1</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>10140</td>
      <td>10147</td>
      <td>ACCCTAAC</td>
      <td>A</td>
      <td>intergenic</td>
      <td>NONE;DDX11L1</td>
      <td>dist=NONE;dist=1727</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>chr1:10140:ACCCTAAC:A</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1</td>
      <td>10146</td>
      <td>10147</td>
      <td>AC</td>
      <td>A</td>
      <td>intergenic</td>
      <td>NONE;DDX11L1</td>
      <td>dist=NONE;dist=1727</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>chr1:10146:AC:A</td>
    </tr>
    <tr>
      <th>2</th>
      <td>1</td>
      <td>10146</td>
      <td>10148</td>
      <td>ACC</td>
      <td>*</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>chr1:10146:ACC:*</td>
    </tr>
    <tr>
      <th>3</th>
      <td>1</td>
      <td>10150</td>
      <td>10151</td>
      <td>CT</td>
      <td>C</td>
      <td>intergenic</td>
      <td>NONE;DDX11L1</td>
      <td>dist=NONE;dist=1723</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>chr1:10150:CT:C</td>
    </tr>
    <tr>
      <th>4</th>
      <td>1</td>
      <td>10172</td>
      <td>10177</td>
      <td>CCCTAA</td>
      <td>C</td>
      <td>intergenic</td>
      <td>NONE;DDX11L1</td>
      <td>dist=NONE;dist=1697</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>chr1:10172:CCCTAA:C</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>995</th>
      <td>1</td>
      <td>66479</td>
      <td>66487</td>
      <td>TATTTATAG</td>
      <td>*</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>chr1:66479:TATTTATAG:*</td>
    </tr>
    <tr>
      <th>996</th>
      <td>1</td>
      <td>66480</td>
      <td>66481</td>
      <td>AT</td>
      <td>A</td>
      <td>intergenic</td>
      <td>FAM138A;OR4F5</td>
      <td>dist=30399;dist=2610</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>chr1:66480:AT:A</td>
    </tr>
    <tr>
      <th>997</th>
      <td>1</td>
      <td>66480</td>
      <td>66483</td>
      <td>ATTT</td>
      <td>*</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>chr1:66480:ATTT:*</td>
    </tr>
    <tr>
      <th>998</th>
      <td>1</td>
      <td>66481</td>
      <td>66488</td>
      <td>TTTATAGA</td>
      <td>T</td>
      <td>intergenic</td>
      <td>FAM138A;OR4F5</td>
      <td>dist=30400;dist=2603</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>chr1:66481:TTTATAGA:T</td>
    </tr>
    <tr>
      <th>999</th>
      <td>1</td>
      <td>66481</td>
      <td>66488</td>
      <td>TTTATAGA</td>
      <td>*</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>chr1:66481:TTTATAGA:*</td>
    </tr>
  </tbody>
</table>
<p>1000 rows × 152 columns</p>
</div>



```python
anno.loc[:,['AF_amr', 'AF_afr', 'AF_nfe']]
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>AF_amr</th>
      <th>AF_afr</th>
      <th>AF_nfe</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0.0006</td>
      <td>0.0008</td>
      <td>0.0007</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0.6380</td>
      <td>0.6300</td>
      <td>0.6413</td>
    </tr>
    <tr>
      <th>2</th>
      <td>.</td>
      <td>.</td>
      <td>.</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0.0357</td>
      <td>0.0426</td>
      <td>0.0370</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0.0086</td>
      <td>0.0097</td>
      <td>0.0084</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>995</th>
      <td>.</td>
      <td>.</td>
      <td>.</td>
    </tr>
    <tr>
      <th>996</th>
      <td>0.0039</td>
      <td>0.0036</td>
      <td>0.0076</td>
    </tr>
    <tr>
      <th>997</th>
      <td>.</td>
      <td>.</td>
      <td>.</td>
    </tr>
    <tr>
      <th>998</th>
      <td>0.0163</td>
      <td>0.0410</td>
      <td>0.0279</td>
    </tr>
    <tr>
      <th>999</th>
      <td>.</td>
      <td>.</td>
      <td>.</td>
    </tr>
  </tbody>
</table>
<p>1000 rows × 3 columns</p>
</div>



Or, create a self-defined annotation file like this:
```
	Chr	Start	AF_amr	AF	AF_nfe	AF_afr
chr1:10140:ACCCTAAC:A	1	10140	0.0006	0.0007	0.0007	0.0008
chr1:10146:AC:A	1	10146	0.638	0.6328	0.6413	0.63
chr1:10150:CT:C	1	10150	0.0357	0.0375	0.037	0.0426
chr1:10172:CCCTAA:C	1	10172	0.0086	0.0082	0.0084	0.0097
chr1:10178:CCTAA:C	1	10178	0.5	0.3333	0.2955	0.4375
chr1:10198:TAACCC:T	1	10198	0.0	0.0	0.0	0.0
chr1:10231:C:A	1	10231	0.2	0.0366	0.0	0.05
chr1:10236:AACCCT:A	1	10236	0.0	0.0	0.0	0.0
chr1:10247:TAAACCCTA:T	1	10247	0.2222	0.2089	0.1429	0.4211
```
The index must match with the ID in vcf file.

- `--pop`, The file of family population information

```python
%%writefile testdata/test_fam_pop.txt
1033 AF_amr
1036 AF_amr
10J_103 AF_afr
10J_109 AF_nfe
10J_112 AF_nfe
10J_119 AF_nfe
```

    Writing ../testdata/test_fam_pop.txt


- `--included-vars`, The file with one column of variants
For example:
```
chr1:10140:ACCCTAAC:A
chr1:10172:CCCTAA:C
chr1:10198:TAACCC:T
chr1:10236:AACCCT:A
chr1:10261:T:TA
chr1:10262:AACCCT:A
```
- `--blueprint`, The blueprint file that defines regional marker (format: "chr startpos endpos name avg.distance male.distance female.distance"). The first four columns are required.

```python
%%writefile testdata/test_blueprint_ext.txt
1       11868   14362   LOC102725121@1  9.177127474362311e-07   1.1657192989882668e-06  6.814189157634088e-07
1       11873   14409   DDX11L1 9.195320788455595e-07   1.1680302941673515e-06  6.82769803434766e-07
1       14361   29370   WASH7P  1.5299877409602128e-06  1.94345806118021e-06    1.136044574393209e-06
1       17368   17436   MIR6859-1@1,MIR6859-2@1,MIR6859-3@1,MIR6859-4@1 1.217692507120495e-06   1.5467668502473368e-06  9.041595098829462e-07
1       30365   30503   MIR1302-10@1,MIR1302-11@1,MIR1302-2@1,MIR1302-9@1       2.1295973889038703e-06  2.705108741548526e-06   1.5812659765416382e-06
1       34610   36081   FAM138A@1,FAM138C@1,FAM138F@1   2.4732411024120156e-06  3.1416201771056266e-06  1.8364278747737466e-06
1       69090   70008   OR4F5   4.866641545668504e-06   6.181823219621424e-06   3.6135725636621673e-06
1       134772  140566  LOC729737       9.633289838108921e-06   1.2236630588823159e-05  7.152898262617822e-06
1       490755  495445  LOC100132062@1,LOC100132287@1   2.2828130832833112e-05  2.8997300893994373e-05  1.6950315013571593e-05
1       450739  451678  OR4F16@1,OR4F29@1,OR4F3@1       2.575942360468604e-05   3.2720758549649544e-05  1.912685483821856e-05
1       627379  629009  LOC101928626    3.943568768003252e-05   5.009295373297854e-05   2.9281737249900675e-05
1       632614  632685  MIR12136        3.974742311959244e-05   5.048893386847169e-05   2.9513206656665908e-05
```

    Overwriting testdata/test_blueprint_ext.txt


Or

```python
%%writefile testdata/test_blueprint.txt
1       11868   14362   LOC102725121@1
1       11873   14409   DDX11L1
1       14361   29370   WASH7P
1       17368   17436   MIR6859-1@1,MIR6859-2@1,MIR6859-3@1,MIR6859-4@1
1       30365   30503   MIR1302-10@1,MIR1302-11@1,MIR1302-2@1,MIR1302-9@1
1       34610   36081   FAM138A@1,FAM138C@1,FAM138F@1
1       69090   70008   OR4F5
1       134772  140566  LOC729737
1       490755  495445  LOC100132062@1,LOC100132287@1
1       450739  451678  OR4F16@1,OR4F29@1,OR4F3@1
1       627379  629009  LOC101928626
1       632614  632685  MIR12136
```

    Overwriting testdata/test_blueprint.txt


## Output format

    - InfoFam: the number of families with the variant or the CHP marker.
### LOD Score. 
It is calculated from 0 to 0.5 with step 0.05 per family per gene. you can change them by `--theta-inc` and `--theta-max`.

    - LOD0: the sum of LOD score at theta=0 among all families
    - LODmax: the max of the sum of LOD score among all families between the range of thetas.
### HLOD Score
    - theta: the theta of best HLOD score.
    - alpha: the alpha of best HLOD score.
    - hlod: the max HLOD of these HLOD between the range of thetas.

### The summary result of CHP markers

```python
result=pd.read_csv('data/test_chp/chr1result_lod_summary.csv',index_col=0)
result
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>chrom</th>
      <th>start</th>
      <th>end</th>
      <th>name</th>
      <th>InfoFam</th>
      <th>LOD0</th>
      <th>LODmax</th>
      <th>theta</th>
      <th>alpha</th>
      <th>hlod</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>MIR6859-1@1,MIR6859-2@1,MIR6859-3@1,MIR6859-4@1</th>
      <td>1</td>
      <td>17368</td>
      <td>17436</td>
      <td>MIR6859-1@1,MIR6859-2@1,MIR6859-3@1,MIR6859-4@1</td>
      <td>3</td>
      <td>-0.864448</td>
      <td>0.000000</td>
      <td>LOD0.5</td>
      <td>0.0</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>WASH7P</th>
      <td>1</td>
      <td>14361</td>
      <td>29370</td>
      <td>WASH7P</td>
      <td>2</td>
      <td>-0.507697</td>
      <td>0.019594</td>
      <td>LOD0.3</td>
      <td>1.0</td>
      <td>0.019594</td>
    </tr>
  </tbody>
</table>
</div>



### The summary result of single variants

```python
result=pd.read_csv('data/test_var/chr1result_lod_summary.csv',index_col=0)
result
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>chrom</th>
      <th>pos</th>
      <th>a0</th>
      <th>a1</th>
      <th>InfoFam</th>
      <th>LOD0</th>
      <th>LODmax</th>
      <th>theta</th>
      <th>alpha</th>
      <th>hlod</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>chr1:13302:C:T</th>
      <td>chr1</td>
      <td>13302</td>
      <td>C</td>
      <td>T</td>
      <td>1</td>
      <td>-0.008503</td>
      <td>0.000000</td>
      <td>LOD0.5</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>chr1:13687:GCCAT:G</th>
      <td>chr1</td>
      <td>13687</td>
      <td>GCCAT</td>
      <td>G</td>
      <td>1</td>
      <td>0.024553</td>
      <td>0.024553</td>
      <td>LOD0.0</td>
      <td>1.000000</td>
      <td>0.024553</td>
    </tr>
    <tr>
      <th>chr1:14464:A:T</th>
      <td>chr1</td>
      <td>14464</td>
      <td>A</td>
      <td>T</td>
      <td>1</td>
      <td>-0.113847</td>
      <td>0.000000</td>
      <td>LOD0.5</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>chr1:14470:G:A</th>
      <td>chr1</td>
      <td>14470</td>
      <td>G</td>
      <td>A</td>
      <td>1</td>
      <td>-0.122592</td>
      <td>0.000000</td>
      <td>LOD0.5</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>chr1:14773:C:T</th>
      <td>chr1</td>
      <td>14773</td>
      <td>C</td>
      <td>T</td>
      <td>1</td>
      <td>-0.007627</td>
      <td>0.000000</td>
      <td>LOD0.5</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>chr1:14843:G:A</th>
      <td>chr1</td>
      <td>14843</td>
      <td>G</td>
      <td>A</td>
      <td>1</td>
      <td>-0.280266</td>
      <td>0.000000</td>
      <td>LOD0.5</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>chr1:14933:G:A</th>
      <td>chr1</td>
      <td>14933</td>
      <td>G</td>
      <td>A</td>
      <td>1</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>LOD0.5</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>chr1:16103:T:G</th>
      <td>chr1</td>
      <td>16103</td>
      <td>T</td>
      <td>G</td>
      <td>4</td>
      <td>-0.414168</td>
      <td>0.079880</td>
      <td>LOD0.0</td>
      <td>0.376008</td>
      <td>0.080043</td>
    </tr>
    <tr>
      <th>chr1:17147:G:A</th>
      <td>chr1</td>
      <td>17147</td>
      <td>G</td>
      <td>A</td>
      <td>1</td>
      <td>-0.000219</td>
      <td>0.000000</td>
      <td>LOD0.5</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>chr1:17358:ACTT:A</th>
      <td>chr1</td>
      <td>17358</td>
      <td>ACTT</td>
      <td>A</td>
      <td>1</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>LOD0.5</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>chr1:17379:G:A</th>
      <td>chr1</td>
      <td>17379</td>
      <td>G</td>
      <td>A</td>
      <td>1</td>
      <td>-0.741666</td>
      <td>0.000000</td>
      <td>LOD0.5</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>chr1:17406:C:T</th>
      <td>chr1</td>
      <td>17406</td>
      <td>C</td>
      <td>T</td>
      <td>1</td>
      <td>-0.122782</td>
      <td>0.000000</td>
      <td>LOD0.5</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>chr1:17407:G:A</th>
      <td>chr1</td>
      <td>17407</td>
      <td>G</td>
      <td>A</td>
      <td>1</td>
      <td>0.016021</td>
      <td>0.016021</td>
      <td>LOD0.0</td>
      <td>1.000000</td>
      <td>0.016021</td>
    </tr>
    <tr>
      <th>chr1:17408:C:G</th>
      <td>chr1</td>
      <td>17408</td>
      <td>C</td>
      <td>G</td>
      <td>1</td>
      <td>0.356048</td>
      <td>0.356048</td>
      <td>LOD0.0</td>
      <td>1.000000</td>
      <td>0.356048</td>
    </tr>
    <tr>
      <th>chr1:17519:G:T</th>
      <td>chr1</td>
      <td>17519</td>
      <td>G</td>
      <td>T</td>
      <td>2</td>
      <td>0.099799</td>
      <td>0.112751</td>
      <td>LOD0.0</td>
      <td>0.728031</td>
      <td>0.114024</td>
    </tr>
    <tr>
      <th>chr1:17594:C:T</th>
      <td>chr1</td>
      <td>17594</td>
      <td>C</td>
      <td>T</td>
      <td>2</td>
      <td>0.100280</td>
      <td>0.113080</td>
      <td>LOD0.0</td>
      <td>0.729654</td>
      <td>0.114300</td>
    </tr>
    <tr>
      <th>chr1:17614:G:A</th>
      <td>chr1</td>
      <td>17614</td>
      <td>G</td>
      <td>A</td>
      <td>1</td>
      <td>-0.120530</td>
      <td>0.000000</td>
      <td>LOD0.5</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>chr1:17716:G:A</th>
      <td>chr1</td>
      <td>17716</td>
      <td>G</td>
      <td>A</td>
      <td>1</td>
      <td>-0.122820</td>
      <td>0.000000</td>
      <td>LOD0.5</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>chr1:17722:A:G</th>
      <td>chr1</td>
      <td>17722</td>
      <td>A</td>
      <td>G</td>
      <td>1</td>
      <td>-0.122079</td>
      <td>0.000000</td>
      <td>LOD0.5</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>chr1:17767:G:A</th>
      <td>chr1</td>
      <td>17767</td>
      <td>G</td>
      <td>A</td>
      <td>1</td>
      <td>-0.122801</td>
      <td>0.000000</td>
      <td>LOD0.5</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>chr1:17928:T:A</th>
      <td>chr1</td>
      <td>17928</td>
      <td>T</td>
      <td>A</td>
      <td>1</td>
      <td>-0.000219</td>
      <td>0.000000</td>
      <td>LOD0.5</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>chr1:17929:C:A</th>
      <td>chr1</td>
      <td>17929</td>
      <td>C</td>
      <td>A</td>
      <td>1</td>
      <td>-0.000219</td>
      <td>0.000000</td>
      <td>LOD0.5</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>chr1:20184:A:G</th>
      <td>chr1</td>
      <td>20184</td>
      <td>A</td>
      <td>G</td>
      <td>1</td>
      <td>-0.000219</td>
      <td>0.000000</td>
      <td>LOD0.5</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>chr1:20231:T:G</th>
      <td>chr1</td>
      <td>20231</td>
      <td>T</td>
      <td>G</td>
      <td>1</td>
      <td>-0.741144</td>
      <td>0.000000</td>
      <td>LOD0.5</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>chr1:20235:G:A</th>
      <td>chr1</td>
      <td>20235</td>
      <td>G</td>
      <td>A</td>
      <td>1</td>
      <td>-0.118623</td>
      <td>0.000000</td>
      <td>LOD0.5</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>chr1:20443:G:A</th>
      <td>chr1</td>
      <td>20443</td>
      <td>G</td>
      <td>A</td>
      <td>1</td>
      <td>-0.280236</td>
      <td>0.000000</td>
      <td>LOD0.5</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>chr1:20485:CA:C</th>
      <td>chr1</td>
      <td>20485</td>
      <td>CA</td>
      <td>C</td>
      <td>1</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>LOD0.5</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>chr1:20522:T:G</th>
      <td>chr1</td>
      <td>20522</td>
      <td>T</td>
      <td>G</td>
      <td>1</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>LOD0.5</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>chr1:29368:G:A</th>
      <td>chr1</td>
      <td>29368</td>
      <td>G</td>
      <td>A</td>
      <td>2</td>
      <td>-0.280978</td>
      <td>0.000000</td>
      <td>LOD0.5</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
  </tbody>
</table>
</div>


