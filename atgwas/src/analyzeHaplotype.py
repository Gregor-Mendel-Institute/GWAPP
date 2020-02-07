#!/usr/bin/env python2.5
"""
11/20/08

Description:
This file analyzes the haplotype structure 

"""
from numpy import *
#from rpy import r
import random, pylab, gc, pdb, math
import dataParsers, snpsdata, phenotypeData

def plot_local_haplotypes(filename, marker_data, focal_start, focal_end, error_tolerance=0, phenotypeData=None, phen_id=None):
	"""
	Plots certain types of haplotype plots...
	"""
	haplotype_ids = range(1000, 0, -1)
	#Fill color matrix up..

	start_i = 0
	cur_pos = 0
	while  start_i < len(marker_data.positions) and cur_pos < focal_start:
		cur_pos = marker_data.positions[start_i]
		start_i += 1
	if start_i == len(marker_data.positions):
		raise Exception("Region is not covered by markers.")
	end_i = start_i
	while  end_i < len(marker_data.positions) and cur_pos < focal_end:
		cur_pos = marker_data.positions[end_i]
		end_i += 1

	center_haplotypes = []
	for a_i in range(len(marker_data.accessions)):
		haplotype = []
		for snp in marker_data.snps[start_i:end_i]:
			haplotype.append(snp[a_i])
		center_haplotypes.append(haplotype)
	haplotype_dict = {}
	hap_pos = (marker_data.positions[end_i - 1] + marker_data.positions[start_i]) / 2
	for a_i, c_h in enumerate(center_haplotypes):
		ch = tuple(c_h)
		if not ch in haplotype_dict:
			haplotype_dict[ch] = [0, 0, [a_i], hap_pos]  #haplotype id number, haplotype frequency count, list of accessions indices, and position.
		else:
			haplotype_dict[ch][2].append(a_i)

	freq_hapl_list = []
	for ch in haplotype_dict:
		hi = haplotype_dict[ch]
		haplotype_dict[ch][1] = len(hi[2])
		freq_hapl_list.append((len(hi[2]), ch))
	freq_hapl_list.sort(reverse=True)
	for (hc, haplotype) in freq_hapl_list:
		if hc == 1:
			haplotype_dict[haplotype][0] = 0
		else:
			haplotype_dict[haplotype][0] = haplotype_ids.pop()

	center_haplotype_dict = haplotype_dict
	left_haplotypes = []
	right_haplotypes = []
	left_haplotypes.append(center_haplotype_dict)
	right_haplotypes = []
	left_positions = [hap_pos]
	right_positions = []

	#Starting with the haplotype structure to the left!
	some_haplotype = True
	i = start_i - 1
	old_hap_dict = center_haplotype_dict
	while old_hap_dict and i >= 0:
		#print i
		#l1  = [len(old_hap_dict[h][2]) for h in old_hap_dict]
		#l2  = [old_hap_dict[h][0] for h in old_hap_dict]
		#print l1,l2, sum(l1)
		haplotype_dict = {}
		hap_pos = marker_data.positions[i]
		left_positions.append(hap_pos)
		for hap in old_hap_dict:
			(h_id, h_count, acc_indices, pos) = old_hap_dict[hap]  #info on the old haplotype
			#print h_id

			temp_hap_dict = {}
			for a_i in acc_indices:
				new_hap = tuple([marker_data.snps[i][a_i]] + list(hap))
				if not new_hap in temp_hap_dict:
					temp_hap_dict[new_hap] = [0, 0, [a_i], hap_pos]  #haplotype id number, haplotype frequency count, list of accessions indices, and position.
				else:
					temp_hap_dict[new_hap][2].append(a_i)

			freq_hapl_list = []
			for h in temp_hap_dict:
				hi = temp_hap_dict[h]
				temp_hap_dict[h][1] = len(hi[2])
				freq_hapl_list.append((len(hi[2]), h))
			freq_hapl_list.sort()

			#print freq_hapl_list
			(hc, h) = freq_hapl_list.pop() #the most frequent haplotype gets colored like the last one.
			if hc == 1:
				del temp_hap_dict[h]
			else:
				temp_hap_dict[h][0] = h_id

			freq_hapl_list.reverse()
			for (hc, h) in freq_hapl_list:
				if hc == 1:
					del temp_hap_dict[h]
				else:
					temp_hap_dict[h][0] = haplotype_ids.pop()
			for h in temp_hap_dict:
				haplotype_dict[h] = temp_hap_dict[h]

		if haplotype_dict:
			left_haplotypes.append(haplotype_dict)
		old_hap_dict = haplotype_dict
		i -= 1



	#Now the same with the haplotype structure to the right!
	i = end_i
	old_hap_dict = center_haplotype_dict
	while old_hap_dict and i < len(marker_data.snps):
		#print i
		#l1  = [len(old_hap_dict[h][2]) for h in old_hap_dict]
		#l2  = [old_hap_dict[h][0] for h in old_hap_dict]
		#print l1,l2, sum(l1)
		haplotype_dict = {}
		hap_pos = marker_data.positions[i]
		right_positions.append(hap_pos)
		for hap in old_hap_dict:
			(h_id, h_count, acc_indices, pos) = old_hap_dict[hap]  #info on the old haplotype

			temp_hap_dict = {}
			for a_i in acc_indices:
				nt = marker_data.snps[i][a_i]
				new_hap = list(hap)
				new_hap.append(nt)
				new_hap = tuple(new_hap)
				#print new_hap
				if not new_hap in temp_hap_dict:
					temp_hap_dict[new_hap] = [0, 0, [a_i], hap_pos]  #haplotype id number, haplotype frequency count, list of accessions indices, and position.
				else:
					temp_hap_dict[new_hap][2].append(a_i)

			freq_hapl_list = []
			for h in temp_hap_dict:
				hi = temp_hap_dict[h]
				temp_hap_dict[h][1] = len(hi[2])
				freq_hapl_list.append((len(hi[2]), h))

			freq_hapl_list.sort()
			(hc, h) = freq_hapl_list.pop() #the most frequent haplotype gets colored like the last one.
			if hc == 1:
				del temp_hap_dict[h]
			else:
				temp_hap_dict[h][0] = h_id

			freq_hapl_list.reverse()
			for (hc, h) in freq_hapl_list:
				if hc == 1:
					del temp_hap_dict[h]
				else:
					temp_hap_dict[h][0] = haplotype_ids.pop()
			for h in temp_hap_dict:
				haplotype_dict[h] = temp_hap_dict[h]

		if haplotype_dict:
			right_haplotypes.append(haplotype_dict)
		old_hap_dict = haplotype_dict
		i += 1


	#Clustering...
	dm = calc_local_dist(marker_data, focal_start, focal_end, error_tolerance=error_tolerance)
	print dm
	import scipy as sp
	import scipy.cluster.hierarchy as hc
	Z = hc.average(dm) #Performing clustering using the dist. matr.
	print Z
	import pylab
	dend_dict = hc.dendrogram(Z, labels=marker_data.accessions)
	new_acc_order = dend_dict['ivl']
	print new_acc_order
	ai_map = [new_acc_order.index(acc) for acc in marker_data.accessions]


	import numpy as np
	#Updating the positions in the figure.
	left_positions.reverse()
	positions = left_positions + right_positions
	x_s = np.zeros((len(positions) + 1, len(marker_data.accessions) + 1))
	start_pos = positions[0] - (0.5 * (positions[1] - positions[0]))
	print len(x_s), len(x_s[0, ])
	for j in range(0, len(x_s[0, ])):
		x_s[0, j] = start_pos
	for j in range(1, len(x_s) - 1): # number of SNPs
		x = positions[j - 1] + 0.5 * (positions[j] - positions[j - 1])
		for k in range(0, len(x_s[j, ])):  # number of NTs
			x_s[j, k] = x
	for j in range(0, len(x_s[0, ])):
		x_s[-1, j] = positions[-1] + (0.5 * (positions[-1] - positions[-2]))


	y_s = np.zeros((len(positions) + 1, len(marker_data.accessions) + 1))
	for j in range(0, len(y_s)): # number of SNPs
		for k in range(0, len(y_s[j, ])):  # number of NTs
			y_s[j, k] = k - 0.5


	#Updating the colors in the figure.
	color_matrix = np.ones((len(positions), len(marker_data.accessions)))
	left_haplotypes.reverse()
	haplotypes = left_haplotypes + right_haplotypes

	max_color = float(haplotype_ids.pop())
	for i, hap_dict in enumerate(haplotypes):
		for h in hap_dict:
			(h_id, h_count, acc_indices, pos) = hap_dict[h]
			for a_i in acc_indices:
				m_ai = ai_map[a_i]
				if h_id == 0:
					color_matrix[i, m_ai] = 1.0
				else:
					color_matrix[i, m_ai] = h_id / max_color

	import phenotypeData as pd
	e_dict = pd._getEcotypeIdInfoDict_()
	accessions = [unicode(e_dict[int(e)][0], 'iso-8859-1') for e in new_acc_order]
	#Plot figure..
	import pylab
	pylab.figure(figsize=(18, 8))
	pylab.axes([0.08, 0.06, 0.9, 0.88])
	pylab.pcolor(x_s, y_s, color_matrix, cmap=pylab.cm.hot)
	#Dealing with the phenotype data
	phenotypeData.removeAccessionsNotInSNPsData(marker_data)
	et_mapping = []
	for i, et in enumerate(new_acc_order):
		et_mapping.append((marker_data.accessions.index(et), i))
	phenotypeData.orderAccessions(et_mapping)
	phen_vals = phenotypeData.getPhenVals(phen_id, noNAs=False)
	acc_strings1 = [accessions[i] + ", " + str(phen_vals[i]) for i in range(len(accessions))]
	acc_strings = [accessions[i] + ", " + str(phen_vals[i]) for i in range(len(accessions))]

	pylab.yticks(range(0, len(marker_data.accessions)), acc_strings, size="small")
	x_range = (x_s[-1, 0] - x_s[0, 0])

	#Retreiving and drawing the genes
	import regionPlotter as rp
	import gwaResults as gr
	genes = gr.get_gene_list(start_pos=x_s[0, 0], end_pos=x_s[-1, 0], chr=5)
	rp.drawGenes(genes, y_shift= -3, rangeVal=40)


	pylab.axis((x_s[0, 0] - 0.05 * x_range, x_s[-1, 0] + 0.05 * x_range, -0.1 * len(marker_data.accessions) - 1, 1.02 * len(marker_data.accessions)))

	pylab.savefig(filename, format='pdf')


def calc_local_dist(marker_data, focal_start, focal_end, error_tolerance=0):
	start_i = 0
	cur_pos = 0
	while  start_i < len(marker_data.positions) and cur_pos < focal_start:
		cur_pos = marker_data.positions[start_i]
		start_i += 1
	if start_i == len(marker_data.positions):
		raise Exception("Region is not covered by markers.")
	end_i = start_i
	while  end_i < len(marker_data.positions) and cur_pos < focal_end:
		cur_pos = marker_data.positions[end_i]
		end_i += 1

	center_haplotypes = []
	for a_i in range(len(marker_data.accessions)):
		haplotype = []
		for snp in marker_data.snps[start_i:end_i]:
			haplotype.append(snp[a_i])
		center_haplotypes.append(haplotype)
	haplotype_dict = {}
	hap_pos = (marker_data.positions[end_i - 1] + marker_data.positions[start_i]) / 2
	for a_i, c_h in enumerate(center_haplotypes):
		ch = tuple(c_h)
		if not ch in haplotype_dict:
			haplotype_dict[ch] = [a_i]  #haplotype id number, haplotype frequency count, list of accessions indices, and position.
		else:
			haplotype_dict[ch].append(a_i)

	import numpy as np
	dm = np.zeros((len(marker_data.accessions), len(marker_data.accessions)))
	for ai1 in range(len(marker_data.accessions)):
		for ai2 in range(0, ai1):
			d = 1
			for h in haplotype_dict:
				acc_indices = haplotype_dict[h]
				if ai1 in acc_indices and ai2 in acc_indices:
					d += 1
			print d
			if d == 2:
				d += 200
				#left search
				i = start_i
				snp = marker_data.snps[i]
				while i > 0 and snp[ai1] == snp[ai2]:
					snp = marker_data.snps[i]
					i -= 1
					d += 1

				#right search
				i = end_i - 1
				snp = marker_data.snps[i]
				while i < len(marker_data.snps) and snp[ai1] == snp[ai2]:
					snp = marker_data.snps[i]
					i += 1
					d += 1

				dm[ai1, ai2] = 1.0 / d
				dm[ai2, ai1] = dm[ai1, ai2]
	return dm




def calc_HS():
	pass

def simulateSNPs(k, n):
	#[pos][acc]
	pass

def sampleSNPs(snps, n, withReplacement=True):
	l = len(snps)
	i = 0
	snpSample = []
	if withReplacement:
		while i < n:
			r = random.randint(0, l - 1)
			snpSample.append(snps[r])
			i += 1
	else:
		while i < n:
			r = random.randint(0, len(snps) - 1)
			snpSample.append(snps[r])
			snps.remove(snps[r])
			i += 1
	return snpSample



def getKinshipDiffs(snpsd, normalGlobalKinship, windowSize=500000, binSize=50000, minSNPs=200):
	kDiffs = []
	binPos = []
	maxPos = snpsd.positions[len(snpsd.positions) - 1]
	numBins = int(maxPos / binSize - windowSize / binSize + 1)
	startPos = 0
	endPos = windowSize
	i = 0
	j = 0
	localKinships = []
	for b_i in range(0, numBins):
		while i < len(snpsd.positions) - 1 and snpsd.positions[i] < startPos:
			i += 1
		while j < len(snpsd.positions) - 1 and snpsd.positions[j] <= endPos:
			j += 1
		if j - 1 - i < minSNPs:
			print "there are two few SNPs"
			kDiffs.append(0)
			localKinships.append(None)
		else:
			print j - 1 - i, "snps found in region ", snpsd.positions[i], "-", snpsd.positions[j - 1]
			if j - 1 - i < minSNPs:
				kDiff = 0
				print "Too few SNPs to estimate K"
			else:
				snps = sampleSNPs(snpsd.snps[i:j], minSNPs, withReplacement=False)
				k = calcKinship(snps)
				localKinships.append(k)
				normal_k = k / mean(k)
				kDiff = mean(abs(normal_k - normalGlobalKinship))
				print "Average k matrix difference:", kDiff
			kDiffs.append(kDiff)
		binPos.append((b_i + (float(windowSize) / binSize) * 0.5) * binSize)
		startPos += binSize
		endPos += binSize
	return (kDiffs, binPos, localKinships)


def getEmmaDiffs(snpsd, phed, p_i, globalKinship, localKinships=None, windowSize=500000, binSize=100000, nSNPs=200, simulate=False, minSNPs=200):
	emmaDiffs = []
	binPos = []
	maxPos = snpsd.positions[len(snpsd.positions) - 1]
	numBins = int(maxPos / binSize - windowSize / binSize + 1)
	startPos = 0
	endPos = windowSize
	i = 0
	j = 0
	for b_i in range(0, numBins):
		while i < len(snpsd.positions) - 1 and snpsd.positions[i] < startPos:
			i += 1
		while j < len(snpsd.positions) - 1 and snpsd.positions[j] <= endPos:
			j += 1
		if j - 1 - i < minSNPs:
			print "there are two few SNPs"
			emmaDiffs.append(0)
		else:
			print j - 1 - i, "snps found in region ", snpsd.positions[i], "-", snpsd.positions[j - 1]
			if localKinships != None and localKinships[b_i] != None:
				k = localKinships[b_i]
			else:
				snps = sampleSNPs(snpsd.snps[i:j], minSNPs)
				k = calcKinship(snps)
			numTries = 0
			unsuccessful = True
			while numTries < 10 and unsuccessful:
				try:
					snps = sampleSNPs(snpsd.snps[i:j], nSNPs)
					emma_res_local = runEmma(phed, p_i, k, snps)
					pvals_local = list(emma_res_local["ps"])
					pvals_local = [-math.log10(pval) for pval in pvals_local]
					emma_res_global = runEmma(phed, p_i, globalKinship, snps)
					pvals_global = list(emma_res_global["ps"])
					pvals_global = [-math.log10(pval) for pval in pvals_global]
					avgPvalDiff = (sum(pvals_global) - sum(pvals_local)) / len(pvals_local)
					print "Average pvalue difference:", avgPvalDiff
					unsuccessful = False
				except Exception:
					print "Emma failed on", (numTries + 1), "try."
					avgPvalDiff = 0
				numTries += 1
			emmaDiffs.append(avgPvalDiff)
		binPos.append((b_i + (float(windowSize) / binSize) * 0.5) * binSize)
		startPos += binSize
		endPos += binSize
	return (emmaDiffs, binPos)





def _locate_max_kinship_acc_(k):
	ai1 = 0
	ai2 = 0
	maxKinship = 0.0
	for i in range(1, len(k)):
		for j in range(0, i):
			if k[i, j] >= maxKinship:
				ai1 = j
				ai2 = i
				maxKinship = k[i, j]
	return (maxKinship, ai1, ai2)

def _merge_accessions_(ai1, ai2, acc_groups):
	acc_groups[ai1] = acc_groups[ai1] + acc_groups[ai2]
	acc_groups.remove(acc_groups[ai2])

def _mergeSNPs_(snpList):
	new_snp = []
	t_snps = zip(*snpList)
	for allele in t_snps:  #For all accessions
		if mean(allele) == 0.5:
			new_snp.append(round(random.random()))
		else:
			new_snp.append(round(mean(allele)))
	return new_snp

def _update_snps_(snps, acc_groups):
	t_snps = map(list, zip(*snps))
	new_snps = []
	for group in acc_groups:
		snpList = []
		for i in group:
			snpList.append(t_snps[i])
		m_snp = _mergeSNPs_(snpList)
		new_snps.append(m_snp)
	new_snps = map(list, zip(*new_snps))
	return new_snps


def _update_phenotype_(phenVals, acc_groups):
	new_phenVals = []
	for group in acc_groups:
		goodCount = 0
		totSum = 0.0
		for a_i in group:
			if phenVals[a_i] != "NA":
				goodCount += 1
				totSum += float(phenVals[a_i])
		#print totSum/goodCount
		phen_val = totSum / goodCount
		new_phenVals.append(phen_val)
	return new_phenVals


def get_KW_pvals(snps, positions, phed, p_i, kinshipThreshold=0.95, method="EMMA"):
	"""
	Takes a set of SNPs, which based on calculates a kinship, and then starts merging strains that look alike, until a threshold is reached.
	
	Methods:
	"Emma"
	"KW"
	"""
	accessions = phed.accessions
	phenVals = phed.getPhenVals(p_i)  #FIXME: make sure that the accessions are ordered correctly.
	pre_count = len(phenVals)
	k = calcKinship(snps)
	acc_groups = []
	for i in range(0, len(accessions)):
		acc_groups.append([i])

	(maxKinship, ai1, ai2) = _locate_max_kinship_acc_(k)
	print "Starting max kinship is", maxKinship
	while maxKinship > kinshipThreshold and len(acc_groups) > 2:
		_merge_accessions_(ai1, ai2, acc_groups) #It updates the acc_groups and map automatically
		merged_snps = _update_snps_(snps, acc_groups)
		print len(merged_snps[0])
		k = calcKinship(merged_snps)
		(maxKinship, ai1, ai2) = _locate_max_kinship_acc_(k)
		print maxKinship

	merged_phenVals = _update_phenotype_(phenVals, acc_groups)
	print "Grouping reduced the # of indiv. by", pre_count - len(merged_phenVals)

	new_snps = []
	new_positions = []
	for i in range(0, len(merged_snps)):
		snp = merged_snps[i]
		if snp.count(0) > 0 and snp.count(1) > 0:
			new_snps.append(snp)
			new_positions.append(positions[i])

	pvals = _run_kw_(new_snps, merged_phenVals)
	#print pvals

	return (pvals, new_positions, acc_groups)

def _plotKinshipDiffs_():

	filterProb = 0.2
	p_i = 1
	res_dir = "/Users/bjarni/tmp/"
	runId = "full_"

	snpsDataFile = "/Network/Data/250k/dataFreeze_080608/250K_f10_080608.csv"
	snpsds = dataParsers.parseCSVData(snpsDataFile, format=1, deliminator=",")#,debug=True)
	phenotypeFile = "/Network/Data/250k/dataFreeze_080608/phenotypes_all_raw_111008.tsv"
	print "Loading phenotype data"
	phed = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter='\t')
	snpsdata.coordinateSnpsAndPhenotypeData(phed, p_i, snpsds)

	for snpsd in snpsds:
		snpsd.filterMinMAF(0.1)
		snpsd.filterMonoMorphicSnps()


	totalSNPs = []
	for i in range(len(snpsds)):
		snpsds[i] = snpsds[i].getSnpsData()
		totalSNPs += snpsds[i].snps

	#For memory, remove random SNPs
	snps = []
	for snp in totalSNPs:
		if random.random() < filterProb:
			snps.append(snp)
	totalSNPs = snps

	print "Calculating the global kinship..."
	globalKinship = calcKinship(totalSNPs)
	print "done."
	normalizedGlobalKinship = globalKinship / mean(globalKinship)
	gc.collect()  #Calling garbage collector, in an attempt to clean up memory..


	for i in range(4, 5):#len(snpsds)):
		chr = i + 1
		snpsd = snpsds[i]
		#pylab.subplot(5,1,chr)
#		pylab.figure(figsize=(18,4))
#		(kinshipDiffs,binPos,local300Kinships) = getKinshipDiffs(snpsd,normalizedGlobalKinship,windowSize=300000)
#		pylab.plot(binPos,kinshipDiffs,"r",label='ws$=300000$')
#		(kinshipDiffs,binPos,local500Kinships) = getKinshipDiffs(snpsd,normalizedGlobalKinship,windowSize=500000)
#		pylab.plot(binPos,kinshipDiffs,"b",label='ws$=500000$')
#		pylab.legend(numpoints=2,handlelen=0.005)
#		pylab.title("Kinship diff. chr. "+str(chr))
#		pylab.savefig(res_dir+runId+"kinshipDiffs_500_300kb_chr"+str(chr)+".pdf",format="pdf")
#		pylab.clf()
		pylab.figure(figsize=(18, 4))
		(emmaDiffs, binPos) = getEmmaDiffs(snpsd, phed, p_i, globalKinship, windowSize=300000)
		pylab.plot(binPos, emmaDiffs, "r", label='ws$=300000$')
		pylab.title("Emma avg. p-value diff. 500kb on chr. " + str(chr))
		(emmaDiffs, binPos) = getEmmaDiffs(snpsd, phed, p_i, globalKinship, windowSize=500000)
		pylab.plot(binPos, emmaDiffs, "b", label='ws$=500000$')
		pylab.title("Emma avg. p-value diff. on chr. " + str(chr))
		pylab.legend(numpoints=2, handlelen=0.005)
		pylab.savefig(res_dir + runId + "EmmaPvalDiffs_500_300kb_chr" + str(chr) + ".pdf", format="pdf")
		pylab.clf()
		gc.collect()  #Calling garbage collector, in an attempt to clean up memory..



def _plotKW_():
	"""
	Analyze how population structure affects KW.
	"""
	filterProb = 0.1
	p_i = 1
	res_dir = "/Users/bjarni/tmp/"
	runId = "_full_quick_"


	snpsDataFile = "/Network/Data/250k/dataFreeze_080608/250K_f10_080608.csv"
	snpsds = dataParsers.parseCSVData(snpsDataFile, format=1, deliminator=",")#,debug=True)
	phenotypeFile = "/Network/Data/250k/dataFreeze_080608/phenotypes_all_raw_111008.tsv"
	print "Loading phenotype data"
	phed = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter='\t')
	snpsdata.coordinateSnpsAndPhenotypeData(phed, p_i, snpsds)

	totalSNPs = []
	for i in range(len(snpsds)):
		snpsds[i] = snpsds[i].getSnpsData()
		totalSNPs += snpsds[i].snps

	#For memory, remove random SNPs
	snps = []
	for snp in totalSNPs:
		if random.random() < filterProb:
			snps.append(snp)
	totalSNPs = snps

	#globalKinship = calcKinship(totalSNPs)
	gc.collect()  #Calling garbage collector, in an attempt to clean up memory..

	#chr = 1
	#for snpsd in snpsds:

	snpsd = snpsds[3]


	k = calcKinship(snpsd.snps[200:1400])
	res = runEmma(phed, p_i, k, snpsd.snps[200:1400]) #runEmma(phed,p_i,k,snps):
	pvals = res["ps"]
	log_pvals = []
	for pval in pvals:
		#print pval
		log_pvals.append(-math.log10(pval))
	pylab.plot(snpsd.positions[200:1400], log_pvals, "c.", label="Emma (local)")

	k = calcKinship(totalSNPs)
	res = runEmma(phed, p_i, k, snpsd.snps[200:1400]) #runEmma(phed,p_i,k,snps):
	pvals = res["ps"]
	log_pvals = []
	for pval in pvals:
		#print pval
		log_pvals.append(-math.log10(pval))
	pylab.plot(snpsd.positions[200:1400], log_pvals, "g.", label="Emma (global)")

	phenVals = phed.getPhenVals(p_i)
	pvals = _run_kw_(snpsd.snps[200:1400], phenVals)
	log_pvals = []
	for pval in pvals:
		#print pval
		log_pvals.append(-math.log10(pval))

	pylab.plot(snpsd.positions[200:1400], log_pvals, "r.", label="KW (full data)")

	(pvals, new_positions, acc_groups) = get_KW_pvals(snpsd.snps[200:1400], snpsd.positions[200:1400], phed, p_i, kinshipThreshold=0.95, method="KW")
	ecot_map = phenotypeData._getEcotypeIdToStockParentDict_()

	for i in range(0, len(acc_groups)):
		acc_list = []
		for a_i in acc_groups[i]:
			e_i = snpsd.accessions[a_i]
			#print e_i
			acc_list.append(ecot_map[int(e_i)][0])
		print "group", i, ":", acc_list

	log_pvals = []
	for pval in pvals:
		#print pval
		log_pvals.append(-math.log10(pval))

	pylab.plot(new_positions, log_pvals, "b.", label="KW (merged data)")

	pylab.legend(numpoints=2, handlelen=0.005)

	pylab.show()


def plotHaplotypes(chr, startPos, endPos):
	snpsd = dataParsers.parseCSVData("/Network/Data/250k/dataFreeze_011209/250K_192_043009.csv")[chr - 1]
	import scipy as sp
	import scipy.cluster.hierarchy as hc
	import Emma
	snpsd = snpsd.getSnpsData()
	newSnps = []
	positions = []
	for i in range(0, len(snpsd.positions)):
		pos = snpsd.positions[i]
		if pos > endPos:
			break
		elif pos >= startPos:
			newSnps.append(snpsd.snps[i])
			positions.append(snpsd.positions[i])

	print "calculating the kinship"
	K = Emma.calcKinship(newSnps)
	#print "K:",K
	Z = hc.average(K)
	#print "Z:",Z
	import pylab
	#hc.leaders(Z)
	dend_dict = hc.dendrogram(Z, labels=snpsd.accessions)
	new_acc_order = dend_dict['ivl']
	print new_acc_order
	print snpsd.accessions
	pylab.savefig("/Users/bjarni/tmp/FRI_tree.pdf", format='pdf')
	#cluster to get ordering??

	acc_mapping = []
	for acc in snpsd.accessions:
		i = new_acc_order.index(acc)
		acc_mapping.append(i)

	snps = []
	for snp in newSnps:
		newSNP = [0] * len(snp)
		for (nt, i) in zip(snp, acc_mapping):
			newSNP[i] = nt
		snps.append(newSNP)

	snps = sp.array(snps)
	pylab.matshow(snps.transpose())
	pylab.savefig("/Users/bjarni/tmp/FRI_haplotype.pdf", format='pdf')


def identity_metric(snps, ignoreSnps=['-', 'N', 'NA']):
	import numpy as np
	seqs = map(list, zip(*snps))
	a = np.zeros((len(seqs), len(seqs)))
	for i, seq1 in enumerate(seqs):
		for j, seq2 in enumerate(seqs):
			error = 0.0
			counts = 0.0
			for nt1, nt2 in zip(seq1, seq2):
				if not (nt1 in ignoreSnps or nt2 in ignoreSnps):
					counts += 1.0
					if nt1 != nt2:
						error += 1.0
			if counts >= 5.0:
				error = error / counts
			else:
				error = 1.0
			a[i, j] = error
	return a

def plot_haplotypes(snps, accessions=None, positions=None, dist_measure="identity", correctScale=False, haplotypeFile=None, treeFile=None, with_genes=False, chr=None):
	"""
	080109 created.
	dist_measure:  type of distance metric used for clustering.
	"""
	import numpy as np
	import scipy as sp
	import scipy.cluster.hierarchy as hc
	if dist_measure == "identity":
		a = identity_metric(snps)
	Z = hc.average(a) #Performing clustering using the dist. matr.
	print Z
	import pylab
	#hc.leaders(Z)
	dend_dict = hc.dendrogram(Z, labels=accessions)
	new_acc_order = dend_dict['ivl']
	print new_acc_order
	#print accessions
	if treeFile:
		pylab.savefig(treeFile, format='pdf')

	acc_mapping = []
	for acc in accessions:
		i = new_acc_order.index(acc)
		acc_mapping.append(i)

	new_snps = []
	for snp in snps:
		new_snp = [0] * len(snp)
		for (nt, i) in zip(snp, acc_mapping):
			new_snp[i] = nt
		new_snps.append(new_snp)
	snps = new_snps

	accessions = new_acc_order

	snps = sp.array(snps)
	class _ntDict_(dict):
		def __missing__(self, key):
			return 0.0
	nt_map = _ntDict_()
	d = {"N":0.0, "NA":0.0, "-":5.0, "A":1.0, "C":2.0, "G":3.0, "T":4.0}
	for key in d:
		#print key,d[key]
		nt_map[key] = d[key]

	for i, nt in enumerate(snps.flat):
		#print nt,nt_map[nt]
		snps.flat[i] = nt_map[nt]
#	for i,nt in enumerate(snps.flat):
#		if nt=='-':
#			snps.flat[i]=5.0
	snps = snps.astype(float)

	#Load genes.. 
	if with_genes:
		import gwaResults as gr
		genes = gr.get_gene_list(chr=chr, start_pos=min(positions), end_pos=max(positions))

	#print snps
	for i in range(0, len(snps), 500):
		pylab.clf()
		if with_genes:
			pylab.figure(figsize=(20, 7))
			pylab.axes([0.04, 0.09, 1.0, 0.87])
		else:
			pylab.figure(figsize=(16, 5))
			pylab.axes([0.06, 0.06, 1.0, 0.9])
		new_snps = snps[i:i + 500]
		if positions: #Prepare SNPs
			x_s = np.zeros((len(new_snps) + 1, len(new_snps[0, ]) + 1))
			start_pos = positions[0] - (0.5 * (positions[1] - positions[0]))
			print len(x_s), len(x_s[0, ])
			for j in range(0, len(x_s[0, ])):
				x_s[0, j] = start_pos
			for j in range(1, len(x_s) - 1): # number of SNPs
				x = positions[j - 1] + 0.5 * (positions[j] - positions[j - 1])
				for k in range(0, len(x_s[j, ])):  # number of NTs
					x_s[j, k] = x
			for j in range(0, len(x_s[0, ])):
				x_s[-1, j] = positions[-1] + (0.5 * (positions[-1] - positions[-2]))


			y_s = np.zeros((len(new_snps) + 1, len(new_snps[0, ]) + 1))
			for j in range(0, len(y_s)): # number of SNPs
				for k in range(0, len(y_s[j, ])):  # number of NTs
					y_s[j, k] = k - 0.5
			#import pdb
			#pdb.set_trace()
#			z_s = []
#			print len(new_snps)
#			for j in range(0,len(new_snps)): # number of SNPs
#				pos = positions[j]
#				for k in range(0,len(new_snps[j,])):  # number of NTs
#					z_s.append(new_snps[j,k])

			pylab.pcolor(x_s, y_s, new_snps)
			x_range = (x_s[-1, 0] - x_s[0, 0])
			pylab.axis((x_s[0, 0] - 0.05 * x_range, x_s[-1, 0] + 0.05 * x_range, -0.1 * len(accessions) - 1, 1.1 * len(accessions)))
		else:
			pylab.pcolor(new_snps.transpose())
		if with_genes:
			import regionPlotter as rp
			rp.drawGenes(genes, y_shift= -1.5)
		pylab.yticks(range(0, len(accessions)), accessions, size="medium")
		cbar = pylab.colorbar(ticks=[0, 1, 2, 3, 4, 5])
		cbar.ax.set_yticklabels(['NA', 'A', 'C', 'G', 'T', '-'])# horizontal colorbar
		#pylab.subplot(2,1,2)

		#Add some more features to the figure...
		if haplotypeFile:
			pylab.savefig(haplotypeFile + "_" + str(i) + ".pdf", format='pdf')
		else:
			pylab.show()
		pylab.clf()

def plot_flc_haplotypes(snps, accessions=None, positions=None, dist_measure="identity", correctScale=False,
		haplotypeFile=None, treeFile=None, acc_250k=None, perlegen_acc=None, flc_250K_positions=None):
	"""
	080609 created.
	dist_measure:  type of distance metric used for clustering.
	"""
	import numpy as np
	import scipy as sp
	import scipy.cluster.hierarchy as hc
	if dist_measure == "identity":
		a = identity_metric(snps)
	elif dist_measure == "":
		pass
	Z = hc.average(a) #Performing clustering using the dist. matr.
	import pylab
	#hc.leaders(Z)
	dend_dict = hc.dendrogram(Z, labels=accessions)
	new_acc_order = dend_dict['ivl']
	#print new_acc_order
	#print accessions
	if treeFile:
		pylab.savefig(treeFile, format='pdf')

	acc_mapping = []
	for acc in accessions:
		i = new_acc_order.index(acc)
		acc_mapping.append(i)

	new_snps = []
	for snp in snps:
		new_snp = [0] * len(snp)
		for (nt, i) in zip(snp, acc_mapping):
			new_snp[i] = nt
		new_snps.append(new_snp)

	snps = sp.array(snps)
	class _ntDict_(dict):
		def __missing__(self, key):
			return 0.0
	nt_map = _ntDict_()
	d = {"N":0.0, "NA":0.0, "-":5.0, "A":1.0, "C":2.0, "G":3.0, "T":4.0}
	for key in d:
		#print key,d[key]
		nt_map[key] = d[key]

	for i, nt in enumerate(snps.flat):
		#print nt,nt_map[nt]
		snps.flat[i] = nt_map[nt]
#	for i,nt in enumerate(snps.flat):
#		if nt=='-':
#			snps.flat[i]=5.0
	snps = snps.astype(float)
	#print snps
	for i in range(0, len(snps), 500):
		pylab.clf()
		fig = pylab.figure(figsize=(16, 5))
		ax1 = pylab.axes([0.06, 0.05, 1.0, 0.9])
		new_snps = snps[i:i + 500]
		if positions: #Prepare SNPs
			x_s = np.zeros((len(new_snps) + 1, len(new_snps[0, ]) + 1))
			start_pos = positions[0] - (0.5 * (positions[1] - positions[0]))
			print len(x_s), len(x_s[0, ])
			for j in range(0, len(x_s[0, ])):
				x_s[0, j] = start_pos
			for j in range(1, len(x_s) - 1): # number of SNPs
				x = positions[j - 1] + 0.5 * (positions[j] - positions[j - 1])
				for k in range(0, len(x_s[j, ])):  # number of NTs
					x_s[j, k] = x
			for j in range(0, len(x_s[0, ])):
				x_s[-1, j] = positions[-1] + (0.5 * (positions[-1] - positions[-2]))


			y_s = np.zeros((len(new_snps) + 1, len(new_snps[0, ]) + 1))
			for j in range(0, len(y_s)): # number of SNPs
				for k in range(0, len(y_s[j, ])):  # number of NTs
					y_s[j, k] = k - 0.5
			#import pdb
			#pdb.set_trace()
#			z_s = []
#			print len(new_snps)
#			for j in range(0,len(new_snps)): # number of SNPs
#				pos = positions[j]
#				for k in range(0,len(new_snps[j,])):  # number of NTs
#					z_s.append(new_snps[j,k])

			pylab.pcolor(x_s, y_s, new_snps)
		else:
			pylab.pcolor(new_snps.transpose())
		yticks, labels = pylab.yticks(range(0, len(accessions)), accessions, size="medium")
		for acc_i in range(0, len(accessions)):
			if acc_i in acc_250k:
				yticks[acc_i].label1.set_color("blue")

		for pos in flc_250K_positions:
			pylab.plot([pos, pos], [-2, -1], color="black")

		x_range = (x_s[-1, 0] - x_s[0, 0])
		y_range = len(accessions)
		pylab.axis((x_s[0, 0] - 0.05 * x_range, x_s[-1, 0] + 0.05 * x_range, -2 - 0.05 * y_range, len(accessions) + 0.05 * y_range))

		cbar = pylab.colorbar(ticks=[0, 1, 2, 3, 4, 5])
		cbar.ax.set_yticklabels(['NA', 'A', 'C', 'G', 'T', '-'])# horizontal colorbar
		#pylab.subplot(2,1,2)
		#ax2 = pylab.axes([0.06,0.04,1.0,0.05],frameon=False)
#		ax2.spines['left'].set_color('none')
#		ax2.spines['right'].set_color('none')
#		ax2.spines['bottom'].set_color('none')
#		ax2.spines['top'].set_color('none')
		#for pos in flc_250K_positions:
		#	ax2.plot([pos,pos],[0,1])
		#ax2.axis((x_s[0,0]-0.05*x_range,x_s[-1,0]+0.05*x_range,-0.05,1.05))


		#Add some more features to the figure...
		if haplotypeFile:
			pylab.savefig(haplotypeFile + "_" + str(i) + ".pdf", format='pdf')
		else:
			pylab.show()
		pylab.clf()


def plot_tree(sd, tree_file, use_emma_kinship=False, verbose=True, kinship_method='ibs'):
	import scipy as sp
	import scipy.cluster.hierarchy as hc
	import pylab
	import phenotypeData
	e_dict = phenotypeData.get_ecotype_id_info_dict()
	#print e_dict
	labels = []
	for acc in snpsd.accessions:
		try:
			s = unicode(e_dict[int(acc)][0], 'iso-8859-1')
		except Exception, err_s:
			print err_s
			print int(acc)
			s = str(acc)
		labels.append(s)
	if verbose:
		print "Plotting tree for SNPs:"
		print "Calculating kinship matrix"
	if kinship_method == 'ibs':
		K = sd.get_ibs_kinship_matrix()
	if kinship_method == 'ibd':
		K = sd.get_ibd_kinship_matrix()
	Z = hc.average(K)
	pylab.figure(figsize=(24, 15))
	pylab.axes([0.03, 0.08, 0.96, 0.91])
	dend_dict = hc.dendrogram(Z, leaf_font_size=7, labels=labels)
	xmin, xmax = pylab.xlim()
	xrange = xmax - xmin
	ymin, ymax = pylab.ylim()
	yrange = ymax - ymin
	pylab.axis([xmin - 0.01 * xrange, xmax + 0.01 * xrange, ymin - 0.02 * yrange, ymax + 0.02 * yrange])
	pylab.savefig(tree_file, format='pdf')
	pylab.clf()
	if verbose:
		print "Done plotting tree, saved in file:", tree_file, "\n"


def plot_250k_Tree(chr=None, startPos=None, endPos=None):
	import scipy as sp
	import scipy.cluster.hierarchy as hc
	import Emma
	import pylab
	import phenotypeData
	e_dict = phenotypeData._getEcotypeIdToStockParentDict_()
	snpsds = dataParsers.parseCSVData("/Network/Data/250k/dataFreeze_011209/250K_192_043009.csv")
	snps = []
	for snpsd in snpsds:
		snps += snpsd.getSnpsData().snps
	snps = sampleSNPs(snps, 100000, False)
	labels = []
	for acc in snpsds[0].accessions:
		try:
			s = unicode(e_dict[int(acc,)][0], 'iso-8859-1')
		except Exception, err_s:
			print err_s
			print e_dict[int(acc)][0]
			s = acc
		labels.append(s)
	print "Calculating kinship matrix"
	K = Emma.calcKinship(snps)
	#print "K:",K
	Z = hc.average(K)
	#print "Z:",Z
	#hc.leaders(Z)
	pylab.figure(figsize=(24, 15))
	pylab.axes([0.03, 0.08, 0.96, 0.91])
	dend_dict = hc.dendrogram(Z, leaf_font_size=7, labels=labels)
	print dend_dict
	pylab.savefig("/Users/bjarni/tmp/250K_192_tree.pdf", format='pdf')
	#cluster to get ordering??



def getLerAndColAccessions(snpsds=None, asFactors=False):
	if not snpsds:
		snpsds = dataParsers.parseCSVData("/Network/Data/250k/dataFreeze_011209/250K_192_043009.csv")
	snpsd = snpsds[3]#.getSnpsData()
	ler_pos = 268809
	col_pos = 269962
	col_accessions = [[], []]
	ler_accessions = [[], []]
	col_factor = []
	ler_factor = []

	for i in range(0, len(snpsd.positions)):
		pos = snpsd.positions[i]
		if pos > col_pos:
			break
		elif pos == ler_pos:
			for j in range(0, len(snpsd.snps[i])):
				if snpsd.snps[i][j] == 0:
					ler_accessions[0].append(snpsd.accessions[j])
					ler_factor.append(0)
				else:
					ler_accessions[1].append(snpsd.accessions[j])
					ler_factor.append(1)

		elif pos == col_pos:
			for j in range(0, len(snpsd.snps[i])):
				if snpsd.snps[i][j] == 0:
					col_accessions[0].append(snpsd.accessions[j])
					col_factor.append(0)
				else:
					col_accessions[1].append(snpsd.accessions[j])
					col_factor.append(1)

	if asFactors:
		return (ler_factor, col_factor)
	else:
		return (ler_accessions[1], col_accessions[1])


def analyze_r_genes(user="bvilhjal", passwd="bjazz32", host="papaya.usc.edu", db="at"):
	"""
	For Nina Riehs.
	"""
	#R-genes of interest
	import gwaResults as gr
	genes = [gr.Gene(chromosome=4, startPos=5962181, endPos=5963791, dbRef="AT4G09420.1"),
		gr.Gene(chromosome=4, startPos=5970932, endPos=5975375, dbRef="AT4G09430.1")]

	#start_pos = 5962181
	#end_pos = 5975375
	start_pos = 5900000
	end_pos = 6000000
	#Figure out what id, Bur and Col have.

	#Retrieve data from DB, regarding Bur and Col
	import MySQLdb
	#from sqlalchemy import select
	print "Connecting to db, host=" + host
	if not user:
		import sys
		sys.stdout.write("Username: ")
		user = sys.stdin.readline().rstrip()
	if not passwd:
		import getpass
		passwd = getpass.getpass()
	try:
		conn = MySQLdb.connect (host=host, user=user, passwd=passwd, db=db)
	except MySQLdb.Error, e:
		print "Error %d: %s" % (e.args[0], e.args[1])
		sys.exit (1)
	cursor = conn.cursor()

	#2010 data:
	s = """
	select distinct l.position, l.offset, g.accession, eva.nativename, al.base
	from at.genotype g, at.allele al, at.locus l, at.alignment an, at.accession2tg_ecotypeid eva 
	where g.allele=al.id and l.id=al.locus and l.alignment=an.id  and an.version>="+dataVersion+locStr+" 
	and g.accession=eva.accession_id and l.chromosome=4 and l.position > 5962180 and l.position < 5975376 
	order by l.position, l.offset, eva.nativename
	"""
	print "Fetcing 2010 data"
	numRows = int(cursor.execute(s))
	print "Number of rows fetched:", numRows
	if numRows:
		alignments = []
		while(1):
			row = cursor.fetchone()
			if not row:
				break;
			print row

	cursor.close ()
	conn.close ()


	#Perlegen data
	db = "chip"
	print "Connecting to db:", db, ", host=" + host
	try:
		conn = MySQLdb.connect (host=host, user=user, passwd=passwd, db=db)
	except MySQLdb.Error, e:
		print "Error %d: %s" % (e.args[0], e.args[1])
		sys.exit (1)
	cursor = conn.cursor()
	s = """
	select distinct d.ecotype 
	from snp_combined_may_9_06_no_van d 
	"""
	s += "where d.chromosome=4 and d.position > " + str(start_pos) + " and d.position < " + str(end_pos) + " and d.mbml98 is not null"
	s += """
	order by d.ecotype	
	"""
	print "Fetcing Perlegen data"
	numRows = int(cursor.execute(s))
	print "Number of rows fetched:", numRows
	accessions = []
	acc_d = dict()
	if numRows:
		i = 0
		alignments = []
		while(1):
			row = cursor.fetchone()
			if not row:
				break;
			acc = row[0]
			accessions.append(acc)
			acc_d[acc] = i
			i += 1
	print accessions

	s = """
	select distinct d.position, d.ecotype, d.mbml98 
	from snp_combined_may_9_06_no_van d 
	"""
	s += "where d.chromosome=4 and d.position > " + str(start_pos) + " and d.position < " + str(end_pos) + " and d.mbml98 is not null"
	s += """
	order by d.position, d.ecotype	
	"""
	numRows = int(cursor.execute(s))
	print "Number of rows fetched:", numRows
	decoder = {'N':'NA', 'A':'A', 'C':'C', 'G':'G', 'T':'T'}
	if numRows:
		alignments = []
		positions = []
		snps = []
		if numRows > 0:
			row = cursor.fetchone()
			newPosition = int(row[0])
			while(1):
				if not row:
					break;
				positions.append(newPosition)
				oldPosition = newPosition
				snp = ['NA'] * 20  #Initialize to missing data.
				while(oldPosition == newPosition):
					snp[acc_d[row[1]]] = decoder[row[2]]
					row = cursor.fetchone()
					if not row:
						break;
					newPosition = int(row[0])
				snps.append(snp)
	print len(snps)
	print positions
	print snps


	cursor.close ()
	conn.close ()


	#Plot haplotypes:
	print len(accessions)
	print len(positions)
	plot_haplotypes(snps, accessions=accessions, positions=positions, haplotypeFile="/Users/bjarnivilhjalmsson/tmp/R-gene_haplotype.pdf", treeFile="/Users/bjarnivilhjalmsson/tmp/R-gene_tree.pdf", with_genes=True, chr=4)




if __name__ == "__main__":
	#plotHaplotypes(4,267000,272000)
	analyze_r_genes()
	#plotTree()
	#getLerAndColAccessions()
	#_plotKinshipDiffs_()
	#_plotKW_()
