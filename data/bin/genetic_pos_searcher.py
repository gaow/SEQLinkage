#!/usr/bin/python2
# Copyright (c) 2014, Hang Dai <hang.dai@bcm.edu>
import subprocess
import sys

genemap_file=sys.argv[1]
chr_min_max_dict={'20': ['83175', '62911391'], '21': ['14756361', '48077812'], '22': ['17060166', '51162059'], '1': ['797163', '248111078'], '3': ['163913', '197803163'], '2': ['88489', '242926320'], '5': ['90199', '180488260'], '4': ['115926', '190921788'], '7': ['162222', '159018204'], '6': ['188969', '170905977'], '9': ['256875', '141017806'], '8': ['440424', '146143866'], '11': ['214421', '134897304'], '10': ['127299', '135203733'], '13': ['19496434', '115013072'], '12': ['294544', '133777645'], '15': ['22815729', '102364127'], '14': ['20445370', '107286497'], '17': ['65090', '81016049'], '16': ['94314', '90112948'], '19': ['1110829', '59093484'], '18': ['159382', '77864213'], 'X': ['175751', '155150393']}

with open ('cM_'+genemap_file,'w') as write_genemap:
	with open(genemap_file,'r') as read_genemap:
		for line in read_genemap:  #below get the information of each gene
			genemap_info=line.split() #list
			chromosome=genemap_info[0]  #str
			start=genemap_info[1]  #str
			end=genemap_info[2]  #str
			pos=(int(start)+int(end))/2  #int
			gene=genemap_info[3] #str
			
			if chromosome!='X':
				command='tabix RUchr_{}_snp_sorted.txt.gz {}:{}-{}'.format(chromosome, chromosome, str(pos), str(pos))
				p=subprocess.Popen(command, universal_newlines=True, shell=True, stdout=subprocess.PIPE)
				out=p.stdout.read() #str
				if out!='':  #luckily we get it by chance
					cMavg=out.split()[3]  #str
					cMfemale=out.split()[4]  #str
					cMmale=out.split()[5]  #str
					new_line='{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chromosome,start,end,gene,cMavg,cMfemale,cMmale)
				else:  #you cannot get it then you guess
					min_pos=chr_min_max_dict[chromosome][0] #str
					max_pos=chr_min_max_dict[chromosome][1] #str
					if pos<int(min_pos):
						command='tabix RUchr_{}_snp_sorted.txt.gz {}:{}-{}'.format(chromosome, chromosome, min_pos, min_pos)
						p=subprocess.Popen(command, universal_newlines=True, shell=True, stdout=subprocess.PIPE)
						out=p.stdout.read() #str
						cMavg=float(out.split()[3])  #float
						cMfemale=float(out.split()[4])  #float
						cMmale=float(out.split()[5])  #float
						interpolated_cMavg=cMavg*pos/int(min_pos)  #float
						interpolated_cMfemale=cMfemale*pos/int(min_pos)  #float
						interpolated_cMmale=cMmale*pos/int(min_pos)  #float
						new_line='{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chromosome,start,end,gene,str(interpolated_cMavg),str(interpolated_cMfemale),str(interpolated_cMmale))
					elif pos>int(max_pos):
						command='tabix RUchr_{}_snp_sorted.txt.gz {}:{}-{}'.format(chromosome, chromosome, max_pos, max_pos)
						p=subprocess.Popen(command, universal_newlines=True, shell=True, stdout=subprocess.PIPE)
						out=p.stdout.read() #str
						cMavg=float(out.split()[3])  #float
						cMfemale=float(out.split()[4])  #float
						cMmale=float(out.split()[5])  #float
						interpolated_cMavg=cMavg*pos/int(max_pos)  #float
						interpolated_cMfemale=cMfemale*pos/int(max_pos)  #float
						interpolated_cMmale=cMmale*pos/int(max_pos)  #float
						new_line='{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chromosome,start,end,gene,str(interpolated_cMavg),str(interpolated_cMfemale),str(interpolated_cMmale))
					else:  #below are the most common situations
						#since we cannot find it, we go up and down to search; each time the step is 100
						#first go up
						out=''
						loopstart=pos-100  #int
						loopend=pos  #int
						while out=='':  #if not find it go on searching
							command='tabix RUchr_{}_snp_sorted.txt.gz {}:{}-{}'.format(chromosome, chromosome, str(loopstart), str(loopend))
							p=subprocess.Popen(command, universal_newlines=True, shell=True, stdout=subprocess.PIPE)
							out=p.stdout.read() #str
							loopstart=loopstart-100
							loopend=loopend-100
						info=out.split('\n')[-2].split() #list  #get the last line
						small_pos=int(info[2])  #int
						small_cMavg=float(info[3])  #float
						small_cMfemale=float(info[4])  #float
						small_cMmale=float(info[5])  #float
						#then go down
						out=''
						loopstart=pos  #int
						loopend=pos+100  #int
						while out=='':  #if not find it go on searching
							command='tabix RUchr_{}_snp_sorted.txt.gz {}:{}-{}'.format(chromosome, chromosome, str(loopstart), str(loopend))
							p=subprocess.Popen(command, universal_newlines=True, shell=True, stdout=subprocess.PIPE)
							out=p.stdout.read() #str
							loopstart=loopstart+100
							loopend=loopend+100
						info=out.split('\n')[0].split() #list  #get the first line
						big_pos=int(info[2])  #int
						big_cMavg=float(info[3])  #float
						big_cMfemale=float(info[4])  #float
						big_cMmale=float(info[5])  #float
						#then calculate the inperpolated cM
						interpolated_cMavg=small_cMavg+(pos-small_pos)*(big_cMavg-small_cMavg)/(big_pos-small_pos)
						interpolated_cMfemale=small_cMfemale+(pos-small_pos)*(big_cMfemale-small_cMfemale)/(big_pos-small_pos)
						interpolated_cMmale=small_cMmale+(pos-small_pos)*(big_cMmale-small_cMmale)/(big_pos-small_pos)
						new_line='{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chromosome,start,end,gene,str(interpolated_cMavg),str(interpolated_cMfemale),str(interpolated_cMmale))
			else:  #now we tackle with X chromosome
				command='tabix RUchr_{}_snp_sorted.txt.gz {}:{}-{}'.format(chromosome, chromosome, str(pos), str(pos))
				p=subprocess.Popen(command, universal_newlines=True, shell=True, stdout=subprocess.PIPE)
				out=p.stdout.read() #str
				if out!='':  #luckily we get it by chance
					cMavg=out.split()[3]  #str
					cMfemale=out.split()[4]  #str
					cMmale=out.split()[5]  #str
					new_line='{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chromosome,start,end,gene,cMavg,cMfemale,cMmale)
				else:
					if pos in range(175751, 2800678) or pos in range(155045646, 155150394):  #here you can find the cMmale and cMfemale, but all cMavg are 'NA'
						#we go up and down to search; each time the step is 100
						#first go up
						out=''
						loopstart=pos-100  #int
						loopend=pos  #int
						while out=='':  #if not find it go on searching
							command='tabix RUchr_{}_snp_sorted.txt.gz {}:{}-{}'.format(chromosome, chromosome, str(loopstart), str(loopend))
							p=subprocess.Popen(command, universal_newlines=True, shell=True, stdout=subprocess.PIPE)
							out=p.stdout.read() #str
							loopstart=loopstart-100
							loopend=loopend-100
						info=out.split('\n')[-2].split() #list  #get the last line
						small_pos=int(info[2])  #int
						small_cMfemale=float(info[4])  #float
						small_cMmale=float(info[5])  #float
						#then go down
						out=''
						loopstart=pos  #int
						loopend=pos+100  #int
						while out=='':  #if not find it go on searching
							command='tabix RUchr_{}_snp_sorted.txt.gz {}:{}-{}'.format(chromosome, chromosome, str(loopstart), str(loopend))
							p=subprocess.Popen(command, universal_newlines=True, shell=True, stdout=subprocess.PIPE)
							out=p.stdout.read() #str
							loopstart=loopstart+100
							loopend=loopend+100
						info=out.split('\n')[0].split() #list  #get the first line
						big_pos=int(info[2])  #int
						big_cMfemale=float(info[4])  #float
						big_cMmale=float(info[5])  #float
						#then calculate the inperpolated cM
						interpolated_cMavg='NA'
						interpolated_cMfemale=small_cMfemale+(pos-small_pos)*(big_cMfemale-small_cMfemale)/(big_pos-small_pos)
						interpolated_cMmale=small_cMmale+(pos-small_pos)*(big_cMmale-small_cMmale)/(big_pos-small_pos)
						new_line='{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chromosome,start,end,gene,str(interpolated_cMavg),str(interpolated_cMfemale),str(interpolated_cMmale))
					else:  #if not in the two above regions, you can only find cMfemale, cMmale and cMavg are 'NA'
						#we go up and down to search; each time the step is 100
						#first go up
						min_pos=chr_min_max_dict[chromosome][0] #str
						max_pos=chr_min_max_dict[chromosome][1] #str
						if pos<int(min_pos):
							command='tabix RUchr_{}_snp_sorted.txt.gz {}:{}-{}'.format(chromosome, chromosome, min_pos, min_pos)
							p=subprocess.Popen(command, universal_newlines=True, shell=True, stdout=subprocess.PIPE)
							out=p.stdout.read() #str
							cMfemale=float(out.split()[4])  #float
							interpolated_cMavg='NA'
							interpolated_cMfemale=cMfemale*pos/int(min_pos)  #float
							interpolated_cMmale='NA'
							new_line='{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chromosome,start,end,gene,str(interpolated_cMavg),str(interpolated_cMfemale),str(interpolated_cMmale))
						elif pos>int(max_pos):
							command='tabix RUchr_{}_snp_sorted.txt.gz {}:{}-{}'.format(chromosome, chromosome, max_pos, max_pos)
							p=subprocess.Popen(command, universal_newlines=True, shell=True, stdout=subprocess.PIPE)
							out=p.stdout.read() #str
							cMfemale=float(out.split()[4])  #float
							interpolated_cMavg='NA'
							interpolated_cMfemale=cMfemale*pos/int(max_pos)  #float
							interpolated_cMmale='NA'
							new_line='{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chromosome,start,end,gene,str(interpolated_cMavg),str(interpolated_cMfemale),str(interpolated_cMmale))
						else:
							out=''
							loopstart=pos-100  #int
							loopend=pos  #int
							while out=='':  #if not find it go on searching
								command='tabix RUchr_{}_snp_sorted.txt.gz {}:{}-{}'.format(chromosome, chromosome, str(loopstart), str(loopend))
								p=subprocess.Popen(command, universal_newlines=True, shell=True, stdout=subprocess.PIPE)
								out=p.stdout.read() #str
								loopstart=loopstart-100
								loopend=loopend-100
							info=out.split('\n')[-2].split() #list  #get the last line
							small_pos=int(info[2])  #int
							small_cMfemale=float(info[4])  #float
							#then go down
							out=''
							loopstart=pos  #int
							loopend=pos+100  #int
							while out=='':  #if not find it go on searching
								command='tabix RUchr_{}_snp_sorted.txt.gz {}:{}-{}'.format(chromosome, chromosome, str(loopstart), str(loopend))
								p=subprocess.Popen(command, universal_newlines=True, shell=True, stdout=subprocess.PIPE)
								out=p.stdout.read() #str
								loopstart=loopstart+100
								loopend=loopend+100
							info=out.split('\n')[0].split() #list  #get the first line
							big_pos=int(info[2])  #int
							big_cMfemale=float(info[4])  #float
							#then calculate the inperpolated cM
							interpolated_cMavg='NA'
							interpolated_cMfemale=small_cMfemale+(pos-small_pos)*(big_cMfemale-small_cMfemale)/(big_pos-small_pos)
							interpolated_cMmale='NA'
							new_line='{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chromosome,start,end,gene,str(interpolated_cMavg),str(interpolated_cMfemale),str(interpolated_cMmale))
			write_genemap.write(new_line)
