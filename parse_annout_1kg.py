#!/usr/bin/python
import sys, re
import subprocess
uncompress='zcat'
list_eff=['nonsynonymous_SNV','synonymous_SNV','stopgain','stoploss','frameshift_insertion','frameshift_deletion']
outdir='/home/um99/'



def global_vars():
	global pos_chr, pos_mut, pos_res, pos_alt, pos_fil, pos_info, pos_rsid
	pos_chr=0
	pos_mut=1
	pos_rsid=2
	pos_res=3
	pos_alt=4
	pos_fil=6
	pos_info=7

def getFreq(freq):
        try:
                return float(freq)
        except: 
                return 0.0

def getMaxFreq(freqList):
        return max(freqList)

def get_variants(kgfile,mf=0.005,veff=list_eff):
	vfile=[]
	print "Attempting extraction..."
	proc = subprocess.Popen([uncompress,'-f',kgfile], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = proc.communicate()
	for line in stdout.split('\n'):
			v=line.split('\t')
			if line=='': continue
			if v[0]=='#CHROM':
				nid=len(v)
				for i in range(9,len(v)):
					vfile.append(open(outdir+'/'+v[i]+'.vars','a'))
			if line[0]!='#':
				vdata={}  
				valt=v[pos_alt].split(',')
				vaf=re.findall('(?<=;AF=).+?(?=;)|(?<=^AF=).+?(?=;)',v[pos_info]) 
				val=re.findall('(?<=ANNOVAR_DATE=).+?(?=ALLELE_END)',v[pos_info])
				if len(vaf)==1:
					af=[float(i) for i in vaf[0].split(',')]
				n=len(val)
				if n==0: continue
				for i in range(n):
					#mut=re.findall('(?<=AAChange.refGene=).+?(?=;)',val[i])[0].split(',')[0].split(':')
					mut=re.findall('(?<=Gene.refGene=).+?(?=;)',val[i])[0].split(',')[0].split(':')
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

                	                maxFreq=getMaxFreq([fa, fa1k, faExAC, faEsp, af[i]])
					print v[pos_chr],v[pos_mut],v[pos_rsid],v[pos_res],valt[i],mut[0],mut[1],mut[-1].replace('p.',''),eff,maxFreq
					vdata[i+1]=[v[pos_chr],v[pos_mut],v[pos_rsid],v[pos_res],valt[i],mut[0],mut[1],mut[-1].replace('p.',''),eff,maxFreq]
				if len(vdata)==0: continue
				for pos_id in range(9,nid):
					gt=[int (i) for i in v[pos_id].split(':')[0].replace('|','/').split('/')]
					if len(gt)==1: gt=[gt[0],gt[0]]
					if (vdata.get(gt[0],0)!=0 and gt[0]!=0 and vdata[gt[0]][-1]<=mf and (vdata[gt[0]][-2] in veff)):
						#print '\t'.join(str(i) for i in vdata[gt[0]])+'\n'
						vfile[pos_id-9].write('\t'.join(str(i) for i in vdata[gt[0]])+'\n')
					if (vdata.get(gt[1],0)!=0 and gt[1]!=0 and gt[1]!=gt[0] and vdata[gt[1]][-1]<=mf and (vdata[gt[1]][-2] in veff)):
						#print '\t'.join(str(i) for i in vdata[gt[1]])+'\n'
						vfile[pos_id-9].write('\t'.join(str(i) for i in vdata[gt[1]])+'\n')
	for i in range(len(vfile)):
		vfile[i].close()


if __name__ == '__main__':
	global_vars()
	filename=sys.argv[1]	
	get_variants(filename)
