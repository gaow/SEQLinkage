#!/bin/sh
#$ -l h_rt=120:00:00
#$ -l h_vmem=128G
#$ -N jupyter-notebook
#$ -o jupyter-notebook-$JOB_ID.out
#$ -e jupyter-notebook-$JOB_ID.err
#$ -j y
#$ -S /bin/bash

# get tunneling info
XDG_RUNTIME_DIR=""
port=$(shuf -i8000-9999 -n1)
node=$(hostname -s)
user=$(whoami)
cluster=csglogin.neuro.columbia.edu  

# print tunneling instructions jupyter-log
echo -e "

MacOS or linux terminal command to create your ssh tunnel
ssh -N -L ${port}:${node}:${port} ${user}@${cluster}

Windows MobaXterm info
Forwarded port:same as remote port
Remote server: ${node}
Remote port: ${port}
SSH server: ${cluster}
SSH login: $user
SSH port: 22

Use a Browser on your local machine to go to:
https://localhost:${port}  (prefix w/ http:// instead if the browser complains that secured connection cannot be established)
" > jupyter-notebook-$JOB_ID.login_info

# env
export http_proxy=http://bcp3.cumc.columbia.edu:8080
export https_proxy=http://bcp3.cumc.columbia.edu:8080
export PATH=$HOME/miniconda3/envs/rpy2/bin:$PATH
# Start jupyter
cd ~/Github/linkage/SEQpy3
#sos run nbs/seqlink_sos.ipynb makehap --cwd data/wg20220316 --fam_path data/new_trim_ped_famless17_no:xx.fam --chrom 11 -j 1
#sos run nbs/seqlink_sos.ipynb makehap --cwd data/wg20220316 --fam_path data/new_trim_ped_famless17_no:xx.fam --chrom 12 -j 1
#sos run nbs/seqlink_sos.ipynb makehap --cwd data/wg20220316 --fam_path data/new_trim_ped_famless17_no:xx.fam --chrom 13 -j 1
#sos run nbs/seqlink_sos.ipynb makehap --cwd data/wg20220316 --fam_path data/new_trim_ped_famless17_no:xx.fam --chrom 14 15 16 17 -j 4
sos run nbs/seqlink_sos.ipynb lods --cwd data/wg20220316 --fam_path data/new_trim_ped_famless17_no:xx.fam --fam_vcf data/wg20220316/fam17_vcf.pickle --chrom 22 -j 1
