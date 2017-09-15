#!/usr/bin/env python
import sys, re


def global_vars():
	global pos_chr, pos_mut, pos_res, pos_alt, pos_fil, pos_info, pos_rsid
	pos_chr=0
	pos_mut=1
	pos_rsid=2
	pos_res=3
	pos_alt=4
	pos_fil=6
	pos_info=7


def get_options():
	import optparse
        desc = 'Script for parsing the annovar table output with gnomAD frequency'
        parser = optparse.OptionParser("usage: python %prog vcf_file pos_data [-h] [-t var_type] [-f frequency]", description=desc)
        parser.add_option('-f', '--threshold',action='store',type='float',dest='th', default=0.005, help='Maximum allele frequecny cutoff')
        parser.add_option('-t','--snv-type', action='store', type='string', dest='snv_type', help='SNV type')
	parser.add_option('--max-af', action='store_true', dest='af_data', help='Search for AF data')
	parser.add_option('-e', '--exclude', action='store', type='string', dest='excluded', help='List of mutations that have to be excluded')
	(options, args) = parser.parse_args()
	veff=['nonsynonymous_SNV']
	th=0.005
	af_data=False
	excluded=""
	if options.th: th=options.th
	if options.snv_type: veff=options.snv_type.split(',')
	if options.af_data: af_data=True
	if options.excluded: excluded=options.excluded
	return args,th,veff,af_data,excluded
	
	


def parse_annofile_af(filename,pos_id,mf=0.005,veff=['nonsynonymous_SNV']):
	sel_muts=[]
	with open(filename) as fh:
		for line in fh:
			vdata=[]
			if line[0]=='#': continue
			v=line.rstrip().split('\t')
			if v[pos_fil]!='PASS': continue
			valt=v[pos_alt].split(',')
			vaf=re.findall('(?<=;AF=).+?(?=;)|(?<=^AF=).+?(?=;)',v[pos_info])
			val=re.findall('(?<=ANNOVAR_DATE=).+?(?=ALLELE_END)',v[pos_info])
			if len(vaf)==1:
				af=[float(i) for i in vaf[0].split(',')]
			n=len(val)
			if n==0: continue
			for i in range(n):
				mut=re.findall('(?<=AAChange.refGene=).+?(?=;)',val[i])[0].split(',')[0].split(':')
				eff=re.findall('(?<=ExonicFunc.refGene=).+?(?=;)',val[i])[0]
				fa=re.findall('(?<=gnomAD_genome_ALL=).+?(?=;)',val[i])[0]
				if mut[0]=='.' or mut[0]=='UNKNOWN': continue
				try:
					fa=float(fa)
				except:
					fa=0.0
				vdata.append([v[pos_chr],v[pos_mut],v[pos_rsid],v[pos_res],valt[i],mut[0],mut[1],mut[-1].replace('p.',''),eff,max(fa,af[i])])
			if len(vdata)==0: continue
			gt=[int (i) for i in v[pos_id].split(':')[0].replace('|','/').split('/')]
			# Option for Chromosome Y
			if len(gt)==1: gt=[gt[0],gt[0]]
			if (gt[0]!=0 and vdata[gt[0]-1][-1]<=mf and (vdata[int(gt[0]-1)][-2] in veff)):
				sel_muts.append(vdata[gt[0]-1])
				print '\t'.join(str(i) for i in vdata[int(gt[0]-1)])
			if (gt[1]!=0 and gt[1]!=gt[0] and vdata[gt[1]-1][-1]<=mf and (vdata[int(gt[1]-1)][-2] in veff)): 
				sel_muts.append(vdata[gt[1]-1])
				print '\t'.join(str(i) for i in vdata[int(gt[1]-1)])
	return sel_muts


def getFreq(freq):
	try:
		return float(freq)
	except:
		return 0.0

def getMaxFreq(freqList):
	return max(freqList)

def parse_annofile(filename,pos_id,mf=0.005,veff=['nonsynonymous_SNV'],excludedList=[]):
	sel_muts=[]
	with open(filename) as fh:
		for line in fh:
			vdata=[]
			if line[0]=='#': continue
			v=line.rstrip().split('\t')
			if v[pos_fil]!='PASS': continue
			valt=v[pos_alt].split(',')
			val=re.findall('(?<=ANNOVAR_DATE=).+?(?=ALLELE_END)',v[pos_info])
			n=len(val)
			if n==0: continue
			for i in range(n):
				if len(excludedList) != 0 and v[pos_chr]+" "+v[pos_mut]+" "+v[pos_rsid]+" "+v[pos_res]+" "+valt[i] in excludedList: continue
				mut=re.findall('(?<=AAChange.refGene=).+?(?=;)',val[i])[0].split(',')[0].split(':')
				eff=re.findall('(?<=ExonicFunc.refGene=).+?(?=;)',val[i])[0]
				fa=re.findall('(?<=gnomAD_genome_ALL=).+?(?=;)',val[i])[0]
				fa1k=re.findall('(?<=1000g2015aug_all=).+?(?=;)', val[i])[0]
				faExAC=re.findall('(?<=ExAC_ALL=).+?(?=;)', val[i])[0]
				faEsp=re.findall('(?<=esp6500siv2_all=).+?(?=;)', val[i])[0]
				if mut[0]=='.' or mut[0]=='UNKNOWN': continue
				fa=getFreq(fa)				
				fa1k=getFreq(fa1k)
				faExAC=getFreq(faExAC)
				faEsp=getFreq(faEsp)

				maxFreq=getMaxFreq([fa, fa1k, faExAC, faEsp])
				vdata.append([v[pos_chr],v[pos_mut],v[pos_rsid],v[pos_res],valt[i],mut[0],mut[1],mut[-1].replace('p.',''),eff,maxFreq])
			if len(vdata)==0: continue
			gt=[int (i) for i in v[pos_id].split(':')[0].replace('|','/').split('/')]
			# Option for Chromosome Y
			if len(gt)==1: gt=[gt[0],gt[0]]
			if (gt[0]!=0 and vdata[gt[0]-1][-1]<=mf and (vdata[int(gt[0]-1)][-2] in veff)):
				sel_muts.append(vdata[gt[0]-1])
				print '\t'.join(str(i) for i in vdata[int(gt[0]-1)])
			if (gt[1]!=0 and gt[1]!=gt[0] and vdata[gt[1]-1][-1]<=mf and (vdata[int(gt[1]-1)][-2] in veff)): 
				sel_muts.append(vdata[gt[1]-1])
				print '\t'.join(str(i) for i in vdata[int(gt[1]-1)])
	return sel_muts



def print_variants(sel_muts):
	for mut in sel_muts:
		(chr,pos,rsid,ref,alt,gene,tid,mut,eff,fa)=mut
		print '\t'.join(str(i) for i in [gene,tid,mut,eff,fa])

def getExcludedMutations(excludedPath):
	excludedList = []
	if len(excludedPath) != 0:
		excludedList = [i.strip() for i in open(excludedPath).readlines()]
	return set(excludedList)


if __name__ == '__main__':
	global_vars()
	args,th,veff,af_data,excludedPath=get_options()
	filein=args[0]
	pos_ind=int(args[1])-1
	excludedList = getExcludedMutations(excludedPath)
	if af_data:
		sel_muts=parse_annofile_af(filein,pos_ind,th,veff)
	else:
		sel_muts=parse_annofile(filein,pos_ind,th,veff,excludedList)
	#print_variants(sel_muts)
