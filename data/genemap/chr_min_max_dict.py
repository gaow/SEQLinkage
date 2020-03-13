import subprocess

chromosome=[]
for i in range(1,23):
	chromosome.append(str(i))
chromosome.append('X')

chr_min_max_dict={}
for item in chromosome:
	command='head -2 RUMap_chr{}.txt |tail -1'.format(item)
	p=subprocess.Popen(command, universal_newlines=True, shell=True, stdout=subprocess.PIPE)
	out=p.stdout.read().split('\t')  #a list
	min_pos=out[6]
	command='tail -1 RUMap_chr{}.txt'.format(item)
	p=subprocess.Popen(command, universal_newlines=True, shell=True, stdout=subprocess.PIPE)
	out=p.stdout.read().split('\t')  #a list
	max_pos=out[6]
	chr_min_max_dict[item]=[min_pos, max_pos]

print chr_min_max_dict
print len(chr_min_max_dict)
