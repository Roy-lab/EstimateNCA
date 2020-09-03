import sys
import argparse
from math import sqrt

parser = argparse.ArgumentParser(description='Take average of TFA profiles.')
parser.add_argument('--tfalist', metavar='FILE', type=str, nargs=1,help='List of TFA files',required=True)
parser.add_argument('--tfaout', metavar='FILE', type=str, nargs=1,help='Average TFA profile',required=True)
args = parser.parse_args()

def getMean(vec):
	l = float(len(vec))
	m = sum(vec)
	m = m/l
	s = sum([(vec[i]-m)*(vec[i]-m) for i in range(len(vec))  ])
	s = s/(l-1)
	s = sqrt(s)
	return (m,s)

def getCorr(vec1,vec2):
	(m1,s1) = getMean(vec1)
	(m2,s2) = getMean(vec2)
	c = sum([(vec1[i]-m1)*(vec2[i]-m2) for i in range(len(vec1))])
	c = c/(s1*s2*(len(vec1)-1))
	return c

def readData(inname):
	data = {}
	f = open(inname,'r')
	for l in f:
		parts = l.strip().split('\t')
		data[parts[0]] = map(float,parts[1:])
	f.close()
	return data

def getAvg(allprof):
	l = len(allprof[0])
	prof = [allprof[0][i] for i in range(l)]
	for j in range(1,len(allprof)):
		c = getCorr(allprof[0],allprof[j])
		s = 1
		if c<0:
			s = -1
		for i in range(l):
			prof[i] = prof[i]+s*allprof[j][i]
	for i in range(l):
		prof[i] = prof[i]/float(len(allprof))
	return prof

def writeData(outname,data,genes):
	f = open(outname,'w')
	for g in genes:
		f.write('%s\t%s\n' % (g,'\t'.join(map(str,data[g]))))
	f.close()

if __name__ == '__main__':
	alltfa = []
	f = open(args.tfalist[0],'r')
	for l in f:
		tfa = readData(l.strip())
		alltfa.append(tfa)
	f.close()
	genes = set(alltfa[0].keys())
	for i in range(1,len(alltfa)):
		genes = genes.union(set(alltfa[i].keys()))
	genes = list(genes)
	genes.sort()
	tfa = {}
	for g in genes:
		allprof = []
		for i in range(len(alltfa)):
			if g in alltfa[i]:
				allprof.append(alltfa[i][g])
		avg = getAvg(allprof)
		tfa[g] = avg
	writeData(args.tfaout[0],tfa,genes)
