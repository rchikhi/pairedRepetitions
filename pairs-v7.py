#!/usr/bin/env python
# $Id$

import sys,os,string,progressbar
sys.path.append("./pysarray")
import pysarray

#faster code execution
psyco=0
if (psyco):
	sys.path.append('./psyco-1.6/lib/python2.5/site-packages')
	import psyco
	psyco.full()

#shows memory usage
heapy=0
if (heapy):
  sys.path.append('./')
  from guppy import hpy
  hp=hpy()

if (heapy): hp.setrelheap()

#dnafile = "NC_000913.fna" # ecoli
dnafile = "NC_001416.fna" #lambda phage
#dnafile = "chromo1-unbout.fna"
#dnafile = "../bigdata/drerio_ref_chr1.fa"
#dnafile = "../bigdata/viral/h1n1.fna"
#dnafile= "test_paired_seq_s50_2,4.fna"
#dnafile = "oan_ref_chr10.fa"
#dnafile="randomseq.fna"

# ******** change these ***********
#            parameters!
print_length=range(1,17)
span=300
delta=30
debug=0
figure=1
bar=1*(1-figure)
span_delta_paired_unicity=1 # doesn't give pixel-perfect results if set to 1, probably because of bad range somewhere
#*********************************

import getopt
opts,args= getopt.getopt(sys.argv[1:],"v:d:",["debug=","delta="])
for o,a in opts:
	if o=="--debug":
		debug=int(a)
	if o=="--delta":
		delta=int(a)
def preprocess():
	file = open(dnafile,'r')
	seql = file.readlines()
	for i in xrange(0,len(seql)): 
		if (seql[i][0]=='>'): seql[i]=" "
	seql = [line.rstrip('\n') for line in seql]
	seq = string.join(seql)
	seq = seq.replace(' ','')

	#also consider the reverse strand
	complement=string.translate(seq, string.maketrans('ACTG', 'TGAC'))
	revseq=complement[::-1]
	len_seq=len(seq)
	seq=seq+'Z'+revseq 

	#(sa,lcp,bw) = pysarray.build_sa_lcp_bw(seq);

	(sa,lcp) = pysarray.build_sarray_and_lcp(seq);
	len_sa=len(sa)
	maxl=len_sa # required here

	span_delta_unicity=1-span_delta_paired_unicity


	#define h[[]] and unique[]
	if (figure==0):
		print "defining some arrays.."
	h=dict()
	vmax=0
	for i in range(len(lcp)):
		if (lcp[i]>vmax): 
			vmax=lcp[i]
		if (lcp[i] not in h): h[lcp[i]]=[]
		h[lcp[i]].append(i)
	maxl=vmax
	return (seq,sa,lcp,len_sa,len_seq,maxl,span_delta_unicity,h)


(seq,sa,lcp,len_sa,len_seq,maxl,span_delta_unicity,h)=preprocess()
if (figure==0):
	print 'beginning analysis with '+str(len(seq))+' nt, (0,' + str(maxl) +') nt reads, ' + str(span)+ ' span, ' +str(delta) + ' delta'


def get_range(i,length):
	if delta==0: return (i+length+span,i+length+span)
	left=i+length+span-delta
	if left<0: left=length+span-delta
	elif left<len_seq+1+length+span and left>len_seq: left=len_seq+1+length+span-delta
	right=i+length+span+delta
	if right>len_sa-length: right=len_sa-length
	elif right>len_seq-length and right<len_seq+length+span: right=len_seq-length
	return (left,right)

def build_r(length):
		r=[]
		r2=[]
		invsa=[]
                lastr=0
                for i in xrange(len_sa):
                        if lcp[i]<length: lastr=i
                        r.append(lastr)
			invsa.append(0)
		for i in xrange(len_sa):
			invsa[sa[i]]=i
			r2.append(r[i])
		for i in xrange(len_sa):
			r[i]=r2[invsa[i]]
                return r

def main():

    d=dict()
    v=dict()
    dupes=0
    for length in range(maxl,0,-1): # note the inverted range
	#find new duplicate reads
	for j in h[length]:
		if (j not in v): 
			dupes+=1
			d[sa[j]]=j
			v[j]=1
			if (debug>1): print seq[sa[j]:sa[j]+length]

		if ((j-1) not in v): 
			dupes+=1
			d[sa[j-1]]=j-1
			v[j-1]=1
			if (debug>1): print seq[sa[j-1]:sa[j-1]+length]
		
	if length in print_length:
		left_range=range(len_seq-span+delta-length-length+1)
		left_range.extend(range(len_seq+1,len(seq)-span+delta-length-length+1))

		p=dict()
		paired_dupes=0
		if (bar): pbar=progressbar.ProgressBar(maxval=2*len(left_range))
		cnt=0

		# create an array r that associates dupes within the suffix array, eg. if sa[i] and sa[i+1]
       # have the same length=length prefix, then r[sa[i]]=r[sa[i+1]]
		r=build_r(length)

		# pass 1: fill p where each pair seen either once or more

		for i in xrange(len_sa):
			if i not in d: continue
			if (i>len_seq-span+delta-length-length+1 and i<len_seq+1) or i>len_sa-span+delta-length-length: continue
			if (bar and i%5000==0):	pbar.update(i)
			(left,right)=get_range(i,length)	
			for k in xrange(left,right+1):
				pair=(r[i],r[k]) #instead of hashing seq[i:i+length],seq[k:k+length]
				if pair not in p: p[pair]=0
				else: p[pair]=1


		# pass 2: check if each read is dupe
		for i in xrange(len_sa):
			if (i>len_seq-span+delta-length-length+1 and i<len_seq+1) or i>len_sa-span+delta-length-length: continue
			if (bar and i%5000==0):	pbar.update(len(left_range)+i)
			(left,right)=get_range(i,length)	
			if i not in d: 
				if delta==0: continue
				# treat dupes within the possible right mate
				# rq: it seems that they NEVER happen in real life
				# rq2: actually its very very rare. 
				# so is the symetric dupe form: i1-k, i2-k. but our code takes care of those. 
				# just make sure to never rule out cases where k not in d in the other parts
				rb=dict()
				for k in xrange(left,right+1):
					if k not in d: continue
					rightblock=r[k] #instead of seq[k:k+length]
					if rightblock not in rb: rb[rightblock]=1
					else: 
						if (debug): print i,seq[i:i+length]+"-"+str(rightblock),"(bigblock dupe)"
						paired_dupes+=1
						if not span_delta_paired_unicity:
							break
				del rb
				continue
			for k in xrange(left,right+1):
				pair=(r[i],r[k])
				if (p[pair]!=0): 
					if (debug): print i,pair[0],"-",pair[1]
					paired_dupes+=1
					if not span_delta_paired_unicity:
						break
		del p
		if (heapy):   print hp.heap()
		if (bar): pbar.finish()


	#another formulation of the next formula
	#uniques=100.0*(len_sa-1-(span+length+length-1)*2-paired_dupes)/(len_sa-1-(span+length+length-1)*2)
	
	if length in print_length:
		left_range=range(len_seq-span+delta-length-length+1)
		left_range.extend(range(len_seq+1,len(seq)-span+delta-length-length+1))
		if span_delta_paired_unicity and delta != 0:
			left_range*=(2*delta)+1
		len_range=max(1,len(left_range))
		uniques=100.0*(len_range-paired_dupes)/len_range
		if (figure):
			print length,uniques
		else:
			print "length:",length,"\tpaired dupes:",paired_dupes,"\tdupes:",dupes,"\tpaired-unique:",uniques

main()	
