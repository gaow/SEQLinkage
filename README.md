# SEQLinkage
> Collapsed Haplotype Pattern Method for Linkage Analysis of Next-Generation Sequencing Data


### Features
-
-
-
-

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

    Traceback (most recent call last):
      File "/home/yh3455/miniconda3/envs/rpy2/bin/seqlink", line 33, in <module>
        sys.exit(load_entry_point('SEQLinkage', 'console_scripts', 'seqlink')())
      File "/mnt/vast/hpc/csg/yin/Github/linkage/SEQpy3/SEQLinkage/Main.py", line 135, in main
        args = Args().get()
      File "/mnt/vast/hpc/csg/yin/Github/linkage/SEQpy3/SEQLinkage/Main.py", line 33, in __init__
        self.getRuntimeArguments(self.parser)
      File "/mnt/vast/hpc/csg/yin/Github/linkage/SEQpy3/SEQLinkage/Main.py", line 72, in getRuntimeArguments
        vargs.add_argument("-h", "--help", action="help", help="Show help message and exit.")
      File "/home/yh3455/miniconda3/envs/rpy2/lib/python3.10/argparse.py", line 1440, in add_argument
        return self._add_action(action)
      File "/home/yh3455/miniconda3/envs/rpy2/lib/python3.10/argparse.py", line 1642, in _add_action
        action = super(_ArgumentGroup, self)._add_action(action)
      File "/home/yh3455/miniconda3/envs/rpy2/lib/python3.10/argparse.py", line 1454, in _add_action
        self._check_conflict(action)
      File "/home/yh3455/miniconda3/envs/rpy2/lib/python3.10/argparse.py", line 1591, in _check_conflict
        conflict_handler(action, confl_optionals)
      File "/home/yh3455/miniconda3/envs/rpy2/lib/python3.10/argparse.py", line 1600, in _handle_conflict_error
        raise ArgumentError(action, message % conflict_string)
    argparse.ArgumentError: argument -h/--help: conflicting option strings: -h, --help


#### 1.run seqlink on CHP marker

#### 2.run seqlink on variants
