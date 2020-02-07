import sys, getopt, traceback, util, pdb, gc
import dataParsers
import phenotypeData
import warnings

import plotResults, gwaResults, pylab

#ft = [1,2,3,4,5,6,7,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,80,81,82]
#ionomics = range(14,32)+range(83,101)
#resDir = "/Network/Data/250k/CEGS_meeting_plots/"

def getTairAnn(startPos, endPos, chr):
	host = "papaya.usc.edu"
	user = "bvilhjal"
	passwd = "bjazz32"
	db = "T8_annotation_TH"

	import MySQLdb
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
	cursor = conn.cursor ()
	#Retrieve the filenames
	print "Fetching data"
	import warnings

	warnings.filterwarnings("ignore")
       	cursor.execute("select distinct pub_locus, start, end from t8_063008 where start > " + str(startPos) + " and end < " + str(endPos) + " and chromosome=" + str(chr) + " and segment_type='gene' order by start")
	genes = []
	currTairID = ""
	while(1):
		try:
			row = cursor.fetchone()
			if not row:
				break;
			gene = gwaResults.Gene()
			gene.startPos = int(row[1])
			gene.endPos = int(row[2])
			gene.tairID = row[0]
			gene.chromosome = chr
			genes.append(gene)
		except Warning:
			pass

	cursor.execute("select distinct t8.tair_id, t8_fd.short_description, t8_fd.description from t8_063008 t8, t8_func_desc t8_fd where t8.pub_locus=t8_fd.tair_id and t8.start > " + str(startPos) + " and t8.end < " + str(endPos) + " and t8.chromosome=" + str(chr) + " order by t8.tair_id")

	functionDescriptions = []
	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		functionDescriptions.append(row)
	cursor.close ()
	conn.close ()

	for gene in genes:
		for fdesc in functionDescriptions:
			if gene.tairID == fdesc[0]:
				gene.shortDescriptions.append(fdesc[1])
				gene.functionDescriptions.append(fdesc[2])


	return genes


def drawGenes(genes, gene_position_cycle=6, rangeVal=15, y_shift=0):
	"""
	2008-09-27 More or less borrowed from Yu Huang..
	"""
	#print "\t Drawing gene model  ..."
	no_of_genes_drawn = 0
	for gene in genes:
		y_value = -y_shift + (no_of_genes_drawn % gene_position_cycle + (rangeVal / 60.0)) * rangeVal / 15.0 #cycling through the y position to avoid clogging
		_plot_one_gene_(gene, y_value=y_value, buffer=0.3 + (rangeVal / 100.0), highlight=gene.highlight)
		no_of_genes_drawn += 1
	#print "Done drawing genes..."



def _plot_one_gene_(gene , y_value=1, buffer=0.4, highlight=False, color=(0.8, 0.4, 0.05)):  #ax: pylab.axis obj.
	"""
	2008-09-29: Code largely borrowed from Yu Huang..		  
	"""
	y_value = buffer + y_value / 2.0
	if highlight:
		pylab.text(gene.startPos, -y_value + 0.08, gene.tairID, fontsize=8, color=color)
		if len(gene.exons) > 0:
			for i in gene.introns:
				pylab.plot([i.startPos, i.endPos], [-y_value, -y_value], color=color, linewidth=2)
			for e in gene.exons:
				pylab.plot([e.startPos, e.endPos], [-y_value, -y_value], color=color, linewidth=4)
		else:
			pylab.plot([gene.startPos, gene.endPos], [-y_value, -y_value], color=color, linewidth=3)
	else:
		pylab.text(gene.startPos, -y_value + 0.08, gene.tairID, size=7)
		if len(gene.exons) > 0:
			for i in gene.introns:
				pylab.plot([i.startPos, i.endPos], [-y_value, -y_value], color=(0.6, 0.6, 0.6), linewidth=1)
			for e in gene.exons:
				pylab.plot([e.startPos, e.endPos], [-y_value, -y_value], color=(0.3, 0.3, 0.3), linewidth=3)
		else:
			pylab.plot([gene.startPos, gene.endPos], [-y_value, -y_value], color=(0.5, 0.5, 0.5), linewidth=2)


class RegionPlotter():
	"""
	Object to plot regions..
	"""


	def __init__(self, phenotypeIndices=None, snpsds=None, results=None, results_map=None, phenotypeFile=None, phed=None):
		if not phenotypeFile and not phed:
			phenotypeFile = "/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/phenotypes_all_raw_042109.tsv"
		elif phed:
			self.phed = phed
		elif phenotypeFile:
			self.phed = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter='\t')
		#if snpsds:
		self.snpsds = snpsds
		#else:
		#	snpsDataFile="/Network/Data/250k/dataFreeze_011209/250K_f13_012509.csv"
		#	self.snpsds = dataParsers.parseCSVData(snpsDataFile, format=1, deliminator=",")
		self.results_map = {}#A map from phenotypes to results (one or more types).
		if results_map:
			self.results_map = results_map
		elif results:
			for result in results:
				if self.results_map.has_key(result.phenotypeID):
					self.results_map[result.phenotypeID].append(result)
				else:
					self.results_map[result.phenotypeID] = [result]
		elif phenotypeIndices:
			self.loadData(phenotypeIndices)

		#print len(results_map[phenotypeIndices[0]])


	def loadData(self, phenotypeIndices):
		"""
		Loads 250K paper data
		
		Deprecated.
		"""
		warnings.warn("This is a deprecated function (regionPlotter.loadData).")
		res_path = "/Network/Data/250k/tmp-bvilhjal/"

#		#SR
#		resultsDirs = [res_path+"kw_results/",res_path+"kw_results/",res_path+"emma_results/",res_path+"emma_results/"]
#		methods=["KW","KW","Emma","Emma"]
#		fileTypes=[".pvals",".sr.pvals",".pvals",".sr.pvals",]
#		datasetNames=["new_raw","new_raw","new_trans","new_trans"]
#		logTransform = [True,True,True,True]
#		mafCutoffs = [0,0,15,15]
#		percentiles = [0.0,0.9,0.0,0.9]
#		interactionResults = [False,True,False,True]

		resultsDirs = [res_path + "kw_results/", res_path + "emma_results/"]
		methods = ["KW", "Emma"]
		fileTypes = [".pvals", ".pvals", ]
		datasetNames = ["raw", "trans"]
		logTransform = [True, True]
		mafCutoffs = [0, 0.1]
		percentiles = [0.0, 0.0]
		interactionResults = [False, False]
		#mrIndex = 1  #The guiding (main) result

		self.results_map = {}
		for i in phenotypeIndices:
			phenName = self.phed.getPhenotypeName(i)

			results = []
			resultTypes = []
			for j in range(0, len(methods)):
				#try:
					resultFile = resultsDirs[j] + methods[j] + "_" + datasetNames[j] + "_" + phenName + fileTypes[j]
					print "Loading result file", resultFile
					rt = gwaResults.ResultType(resultType=methods[j], fileType="pvals", mafCutoff=mafCutoffs[j])
					if self.snpsds:
						result = gwaResults.SNPResult(resultFile, self.snpsds, name=methods[j] + "_" + datasetNames[j] + "_" + phenName, resultType=rt, phenotypeID=i)
					else:
						result = gwaResults.Result(resultFile, name=methods[j] + "_" + datasetNames[j] + "_" + phenName, resultType=rt, phenotypeID=i, interactionResult=interactionResults[j])
					if logTransform[j]:
						print "Log transforming the p-values.."
						result.negLogTransform()

					#result.filterMAF(minMaf=mafCutoffs[j])
					result.filterMARF(minMaf=mafCutoffs[j])
					if percentiles[j] > 0:
						result.filterPercentile(percentiles[j])
					results.append(result)
					resultTypes.append(methods[j])
				#except Exception:
				#	print "Couldn't load",resultFile

			self.results_map[i] = results
			gc.collect()  #Calling garbage collector, in an attempt to clean up memory..

	def plot_small_result(self, results, pdf_file=None, png_file=None, plot_genes=True, tair_file=None, max_val=None, min_val=None,
			bonferroni_level=None, highlight_gene_ids=None, caption=None, ylab=None, xlab=None, results_colors=None):
		"""
		A simple function that plots small regional results.
		
		Only handles 1 result!!!.
		"""
		print results
		result_chr = results[0].chromosomes[0]
		for result in results:
			if len(result.scores) > 1000:
				warnings.warn("This result object contains " + str(len(result.scores)) + " scores.")
			if len(set(result.chromosomes)) > 1:
				raise Exception("Can't plot region that covers more than one chromosome.")
			if len(result.positions) != len(result.scores):
				raise Exception("Number of scores and positions doesn't equal!!!")
			if result_chr != result.chromosomes[0]:
				raise Exception("Can't plot results from different chromosomes.")

		#Preprocess positions, scores, and colors.
		default_colors = ['b', 'r', 'g']

		scores = []
		positions = []
		colors = []
		for i, result in enumerate(results):
			for j in range(len(result.scores)):
				positions.append(result.positions[j])
				scores.append(result.scores[j])
			if results_colors:
				colors += [results_colors[i]] * len(result.scores)
			else:
				colors += [default_colors[i]] * len(result.scores)

		min_x = min(positions)
		max_x = max(positions)
		if max_val:
			max_x = max(max_x, max_val)
		if min_val:
			min_x = min(min_x, min_val)


		chromosome = result.chromosomes[0]
		tair_genes = None
		highlight_genes = None
		if plot_genes or tair_file:
			tair_genes = gwaResults.get_gene_list(min_x - 500, max_x + 500, chromosome)
			if highlight_gene_ids:
				highlight_genes = []
				for g in tair_genes:
					if g.dbRef in highlight_gene_ids:
						highlight_genes.append(g)
			print "Found", len(tair_genes), "genes in region:", min_x, "to", max_x


		if tair_file:
			f = open(tair_file, 'w')
			for gene in tair_genes:
				f.write(str(gene) + "\n")
			f.close()

		if not caption:
			caption = result.name

		self._plot_simple_result_(scores, positions, colors, min_x, max_x, caption=caption, pdf_file=pdf_file,
					png_file=png_file, tair_genes=tair_genes, highlight_genes=highlight_genes,
					bonferroni_level=bonferroni_level)


	def _plot_simple_result_(self, scores, positions, colors, min_x, max_x, caption="", pdf_file=None, png_file=None, tair_genes=None,
				highlight_genes=None, max_val=None, min_val=None, bonferroni_level=None, color=None):
		"""
		A helper function actually drawing the simple plot.
		"""

		print "Plotting", len(scores), "SNPs."


		x_range = max_x - min_x
		y_min = min(scores)
		y_max = max(scores)
		y_range = y_max - y_min
		x_buffer = 0.05
		y_buffer = 0.05


		#plot box boundaries [mix_x,max_x,min_y,max_y]


		if tair_genes:
			pylab.figure(figsize=(8, 6))
			plot_range = [min_x - x_buffer * x_range, max_x + x_buffer * x_range,
				y_min - y_range * (y_buffer + 0.25), y_max + y_range * y_buffer]
		else:
			pylab.figure(figsize=(8, 4))
			plot_range = [min_x - x_buffer * x_range, max_x + x_buffer * x_range,
				y_min - y_range * y_buffer, y_max + y_range * y_buffer]

		pylab.axes([0.06, 0.10, 0.90, 0.82])
		for (position, score, color) in zip(positions, scores, colors):
			pylab.plot(position, score, color=color, marker=".", ls="", alpha=0.8)
		#pylab.plot(positions, scores, marker=".", ls="", alpha=0.8)
		#pylab.text(endPos-0.025*posRange,maxVal-0.15*rangeVal,l,size="x-large")							
		#pylab.legend(numpoints=1,handlelen=0.005,prop=fontProp)
		if bonferroni_level:
			pylab.plot([plot_range[0], plot_range[1]], [bonferroni_level, bonferroni_level], "k-.")
		if tair_genes:
			num_rows = max(3, min(6, 1 + len(tair_genes) / 10))
			print "Drawing TAIR genes"
			drawGenes(tair_genes, gene_position_cycle=num_rows, rangeVal=y_range * 1.5)
		pylab.axis(plot_range)
		pylab.xlabel('Position')
		pylab.ylabel('-log(p-value)')
		if pdf_file:
			pylab.savefig(pdf_file, format="pdf")
		if png_file:
			pylab.savefig(png_file, format="png", dpi=300)
		if not (pdf_file or png_file):
			pylab.show()
		pylab.close()





	def plotReg(self, region, phenotypeID, snpPos=None, pdfFile=None, pngFile=None, tairFile=None, plotGenes=True, printTairInfo=True, binary=False, results=None):
		return self.plotRegion(region.chromosome, region.startPos, region.endPos, phenotypeID, snpPos=snpPos, pdfFile=pdfFile, pngFile=pngFile, tairFile=tairFile, plotGenes=plotGenes, printTairInfo=printTairInfo, binary=binary, results=results)

	def plotRegion(self, chr, startPos, endPos, phenotypeID, snpPos=None, pdfFile=None, pngFile=None, tairFile=None, plotGenes=True, printTairInfo=True, binary=False,
			results=None, maxLogPval=None, plotColors=None, plotOrder=None, bonferoniCorrections=[6.63], labels=None, highlightGenes=None, minVal=0, markers=None,
			alphas=None):
		if not results:
			results = self.results_map[phenotypeID]
		resultTypes = [result.resultType.resultType for result in results]
		phenName_list = self.phed.getPhenotypeName(phenotypeID).split("_")[1:]
		phenName = " ".join(phenName_list)
		tairGenes = gwaResults.get_gene_list(startPos, endPos, chr)

		if tairFile:
			f = open(tairFile, 'w')
			for gene in tairGenes:
				f.write(str(gene) + "\n")
			f.close()
		if printTairInfo:
			print "\nGene annotation:"
			for gene in tairGenes:
				print gene
			print '\n'

		tairInfo = []
		for gene in tairGenes:
			tairInfo.append(str(gene))

		genes = None
		if plotGenes:
			genes = tairGenes
		if not snpPos:
			self.drawSNPPlot(chr, startPos, endPos, phenName, results, resultTypes, genes=genes, pdfFile=pdfFile, pngFile=pngFile, binary=binary,
					maxLogPval=maxLogPval, plotColors=plotColors, plotOrder=plotOrder, bonferoniCorrections=bonferoniCorrections,
					labels=labels, highlightGenes=highlightGenes, minVal=minVal, markers=markers, alphas=alphas)
		else:
			snpsd = self.snpsds[chr - 1]
			pos = snpsd.positions[0]
			posDiff = abs(pos - snpPos)
			i = 0
			while i + 1 < len(snpsd.positions) and abs(snpsd.positions[i + 1] - snpPos) < posDiff:
				i += 1
				posDiff = abs(snpsd.positions[i] - snpPos)

			if posDiff == 0:
				print "SNP was found."
			else:
				print "SNP was not found.  Using the nearest SNP instead."
				print chr, snpsd.positions[i]

			snp = gwaResults.SNP(snpsd.positions[i], chr, alleles=snpsd.snps[i])
			self.drawSNPPlot(chr, startPos, endPos, phenName, results, resultTypes, genes=genes, ldSNP=snp, pdfFile=pdfFile, pngFile=pngFile, binary=binary,
					maxLogPval=maxLogPval, plotColors=plotColors, plotOrder=plotOrder, bonferoniCorrections=bonferoniCorrections, labels=labels,
					highlightGenes=highlightGenes, minVal=minVal, markers=markers, alphas=alphas)

		gc.collect()  #Calling garbage collector, in an attempt to clean up memory..

		return tairInfo



	def drawSNPPlot(self, chr, startPos, endPos, phenName, results, resultTypes, genes=None, ldSNP=None, pdfFile=None,
			pngFile=None, binary=False, maxLogPval=None, plotColors=None, plotOrder=None, bonferoniCorrections=[6.63],
			labels=None, highlightGenes=None, minVal=0, markers=None, alphas=None, separate=True, highlightPosList=[268809, 269962]):
		"""
		Draws a snp-plot
		
		At least on result is required.
		
		requires pylab to be installed.
		"""
		import matplotlib

   		if not markers:
   			markers = ["."] * len(results)

   		if not alphas:
   			alphas = [0.7] * len(results)
   		print "alphas:", alphas

   		if not plotColors:
   			plotColors = []
   			defaultColors = ['b', 'r', 'c', 'g', 'm', 'y']
   			for i in range(0, len(results)):
   				plotColors.append(defaultColors[i % len(defaultColors)])

		if not plotOrder:
			plotOrder = range(1, len(results) + 1)

		if not labels:
   			labels = [""] * len(results)

		import pylab
		size = endPos - startPos

		#throwing out unnecessary highlighting genes
		new_hl_genes = []
		if highlightGenes:
			for gene in highlightGenes:
				if type(gene) == str:
					for g in genes:
						if g.tairID == gene:
							new_hl_genes.append(g)
							g.highlight = True
				else:
					#print gene.chromosome,chr,":",gene.startPos,endPos,":",gene.endPos,endPos
					if gene.chromosome == chr and gene.startPos < endPos and gene.endPos > startPos:
						new_hl_genes.append(gene)
		highlightGenes = new_hl_genes
		print "highlightGenes:", highlightGenes

		hlPositionsList = []
		hlScoreList = []
		positionsList = []
		scoreList = []
		newMaxPos = 0
		for result in results:
			positions = []
			scores = []
			hlPositions = []
			hlScores = []
			i = 0

			currPos = result.positions[0]
			currChr = result.chromosomes[0]

			while currChr < chr:
				i += 1
				currChr = result.chromosomes[i]
			#Found chromsome..

			while currChr == chr and currPos < startPos:
				i += 1
				currPos = result.positions[i]
				currChr = result.chromosomes[i]
			#Found start..

			while currChr == chr and currPos < endPos:
				currPos = result.positions[i]
				pos = currPos
				if result.interactionResult:
					pos = []
					pos.append(currPos)
					pos += result.interactionPositions[i]

				if pos in highlightPosList:
					hlPositions.append(pos)
					if maxLogPval and result.scores[i] > maxLogPval:
						hlScores.append(maxLogPval)
					else:
						hlScores.append(result.scores[i])
				else:
					positions.append(pos)
					#positions.append(pos/1000000.0)			
					if maxLogPval and result.scores[i] > maxLogPval:
						scores.append(maxLogPval)
					else:
						scores.append(result.scores[i])
				i += 1
				currChr = result.chromosomes[i]
				currPos = result.positions[i]

			positionsList.append(positions)
			scoreList.append(scores)
			hlPositionsList.append(hlPositions)
			hlScoreList.append(hlScores)
			if result.interactionResult:
				for p_list in positions:
					if max(p_list) > newMaxPos:
						newMaxPos = max(p_list)

		#Update maxVal for local graph
		if not maxLogPval:
   			maxList = []
   			for scores in scoreList:
   				maxList.append(max(scores))
   			maxVal = max(maxList)
     		else:
     	  		maxVal = maxLogPval
		rangeVal = maxVal - minVal

		endPos = max(endPos, newMaxPos)
		#startPos = startPos/1000000.0
		#endPos = endPos/1000000.0
		posRange = endPos - startPos


		#Now plotting...
		print "Now plotting.."
#		if genes:
#			numGeneLines= int(3+len(genes)/14)
#			#pylab.figure(1,figsize=(8,4+0.4*numGeneLines))
#			pylab.figure(1,figsize=(12,4+0.4*numGeneLines))
#		else:
#			pylab.figure(1,figsize=(8,8))
		#pylab.figure(1,figsize=(8,4))
		pylab.figure(1, figsize=(8, 8))
		fontProp = matplotlib.font_manager.FontProperties(size=12)

		for i in plotOrder:
			print "plotting", i, " result."
			positions = positionsList[i - 1]
			scores = scoreList[i - 1]
			hlPositions = hlPositionsList[i - 1]
			hlScores = hlScoreList[i - 1]
			c = plotColors[i - 1]
			l = labels[i - 1]
			print "interactionResult:", results[i - 1].interactionResult
			if separate:
				pylab.subplot(len(plotOrder) * 100 + 10 + i)
			if results[i - 1].interactionResult and len(positions) > 0:
				for (pos, score) in zip(positions, scores):
					p_list = list(pos)
					s_list = [score] * len(p_list)
					pylab.plot(p_list, s_list, ":", color=c, ms=1.8, lw=0.8)
					for (p, s) in zip(p_list, s_list):
						pylab.plot([p], [s], color=c, ms=1.8, marker=markers[i - 1])
				pylab.plot([p], [s], color=c, label=l, ms=1.8, marker=markers[i - 1])
				#pylab.plot([p],[s],color=c,ms=1.8,marker=markers[i-1])
				#pylab.text(endPos-0.08*posRange,maxVal-0.08*rangeVal,l)							
			else:
				pylab.plot(positions, scores, color=c, marker=markers[i - 1], ls="", alpha=alphas[i - 1])
				pylab.xticks([0, 500000, 1000000, 1500000], ['0', '500,000', '1,000,000', '1,500,000'])
				#pylab.plot(positions,scores,color=c,label=l,marker=markers[i-1],ls="",alpha=alphas[i-1])
				pylab.text(endPos - 0.025 * posRange, maxVal - 0.15 * rangeVal, l, size="x-large")
				#pylab.text(endPos*1.07,minVal+rangeVal*0.40-rangeVal*i*0.065,l,rotation=90,)
				#pylab.legend(numpoints=1,handlelen=0.005,prop=fontProp)
				if separate:
					for bfc in bonferoniCorrections:
						pylab.plot([startPos - 0.05 * posRange, endPos + 0.05 * posRange], [bfc, bfc], "k-.")
					j = 1
					print "almost there"
					for gene in highlightGenes:
						fc = (0.6, 0.6, 0.6)
						ec = (0.6, 0.6, 0.6)
						alpha = 0.6
						print i, j, gene.name, gene.startPos, gene.endPos
						if j == 2: #FRI hack
							print "found FRI"
							fc = (0.8, 0.6, 0.1)#(0.1,0.1,0.9)#
							ec = (0.8, 0.6, 0.1)#(0.1,0.1,0.9)#
							alpha = 0.8
							if i == 1:
								pylab.text(235000, maxVal + 0.08 * rangeVal, "FRI", color=(0.8, 0.6, 0.1), size="large")#(0.1,0.1,0.9))
							pylab.axvspan(gene.startPos, gene.endPos, facecolor=fc, edgecolor=fc, alpha=alpha, linewidth=1.2)
						j += 1
					#pylab.axis([startPos-0.05*posRange,endPos+0.05*posRange,minVal-rangeVal*0.05,maxVal+rangeVal*0.5])

				if i == 1:
					pass
					#pylab.title(phenName+": chromosome "+str(chr)+".")
				#elif i==2:
				#	pylab.ylabel("-log$_{10}$(p-value)")
				elif i == 3:
					pylab.text(-200000, 6.6, "-log$_{10}$(p-value)", rotation='vertical')
				elif i == len(plotOrder):
					pylab.xlabel("Position (bases)")
				hlColors = [(0.8, 0.1, 0.2), (0.2, 0.1, 0.8)]
				for i in range(0, len(hlPositions)):
					if hlPositions[i] == 268809:
						pylab.plot([hlPositions[i]], [hlScores[i]], color=hlColors[0], marker=markers[i - 1], ls="", alpha=alphas[i - 1])
					else:
						pylab.plot([hlPositions[i]], [hlScores[i]], color=hlColors[1], marker=markers[i - 1], ls="", alpha=alphas[i - 1])
				pylab.axis([startPos - 0.05 * posRange, endPos + 0.05 * posRange, minVal - rangeVal * 0.05, maxVal + rangeVal * 0.05])



		for bfc in bonferoniCorrections:
			if not separate:
				pylab.plot([startPos - 0.05 * posRange, endPos + 0.05 * posRange], [bfc, bfc], "k-.")

		#pylab.axis([startPos-0.075*posRange,endPos+0.05*posRange,minVal-rangeVal*0.05,maxVal+rangeVal*0.20])
#		if genes and not separate:
#			self.drawGenes(genes, gene_position_cycle=numGeneLines,rangeVal=rangeVal)
#			geneSpace = 0.49*rangeVal/15.0
#			#pylab.axis([startPos-0.05*posRange,endPos+0.05*posRange,minVal-rangeVal*0.05-0.49*numGeneLines,maxVal+rangeVal*0.05])
#			pylab.axis([startPos-0.05*posRange,endPos+0.05*posRange,minVal-rangeVal*0.05-geneSpace*numGeneLines,maxVal+rangeVal*0.075])
#		else:
#			pylab.axis([startPos-0.05*posRange,endPos+0.05*posRange,minVal-rangeVal*0.05,maxVal+rangeVal*0.05])

		if highlightGenes and not separate:
			i = 1
			for gene in highlightGenes:
				fc = (0.8, 0.5, 0.1)
				ec = (0.8, 0.5, 0.1)
				alpha = 0.6
				print i, gene.name, gene.startPos, gene.endPos
				#pylab.axvspan(gene.startPos,gene.endPos, facecolor=fc,edgecolor=fc, alpha=alpha,linewidth=1.2)
				#pylab.axvspan(gene.startPos/1000000.0,gene.endPos/1000000.0, facecolor=fc,edgecolor=fc, alpha=alpha,linewidth=1.2)
				if i == 2: #FRI hack
					fc = (0.8, 0.6, 0.1)
					ec = (0.8, 0.6, 0.1)
					alpha = 0.8
					pylab.text(265000, maxVal + 0.075 * rangeVal, "FRI", rotation=45, color=(0.8, 0.6, 0.1))
					pylab.axvspan(gene.startPos, gene.endPos, facecolor=fc, edgecolor=ec, alpha=alpha, linewidth=1.2)

				i += 1


		pylab.subplots_adjust(right=0.94)
		pylab.subplots_adjust(left=0.08)
		pylab.subplots_adjust(bottom=0.07)
		pylab.subplots_adjust(top=0.94)
#		pylab.subplots_adjust(right=0.98)
#		pylab.subplots_adjust(left=0.07)
#		pylab.subplots_adjust(bottom=0.1)
#		pylab.subplots_adjust(top=0.985)
		if not separate:
			#pylab.title(phenName+": chromosome "+str(chr)+".")
			pylab.legend(numpoints=1, handlelen=0.005, prop=fontProp, loc=2)
			pylab.xlabel("Position (Mb)")
			pylab.ylabel("-log$_{10}$(p-value)")

		if pdfFile:
			pylab.savefig(pdfFile, format="pdf")
		if pngFile:
			pylab.savefig(pngFile, format="png", dpi=300)
		if not (pdfFile or pngFile):
			pylab.show()

		pylab.close(1)


	def drawSNPPlot_old(self, chr, startPos, endPos, phenName, results, resultTypes, genes=None, ldSNP=None, pdfFile=None,
			pngFile=None, binary=False, maxLogPval=12, interactionResult=[False]):
		"""
		Draws a snp-plot
		
		At least on result is required.
		
		requires pylab to be installed.
		"""

		import pylab
		size = endPos - startPos

		maxVal = maxLogPval
		minVal = 0
		rangeVal = maxVal - minVal
		pvals = [True, True, True, True]

		positionsList = []
		scoreList = []
		snps = []

		result = results[0]
		positions = []
		scores = []
		i = 0
		currPos = result.positions[0]
		currChr = result.chromosomes[0]
		if interactionResults[0]:
			iPositions = []

		while currChr < chr:
			i += 1
			currChr = result.chromosomes[i]
		#Found chromsome..

		while currChr == chr and currPos < startPos:
			i += 1
			currPos = result.positions[i]
			currChr = result.chromosomes[i]
		#Found start..

		while currChr == chr and currPos < endPos:
			currPos = result.positions[i]
			pos = []
			pos.append(currPos)
			if interactionResults[0]:
				pos += result.interactionPositions
			currChr = result.chromosomes[i]
			if ldSNP:
				snps.append(result.snps[i])
			positions.append(pos)
			if result.scores[i] > maxVal and pvals[0]:
				scores.append(maxVal)
			else:
				scores.append(result.scores[i])
			i += 1

		positionsList.append(positions)
		scoreList.append(scores)


		for j in range(1, len(results)):
			result = results[j]
			positions = []
			scores = []
			i = 0
			currPos = result.positions[0]
			currChr = result.chromosomes[0]
			while currChr < chr:
				i += 1
				currChr = result.chromosomes[i]
			#Found chromsome..

			while currChr == chr and currPos < startPos:
				i += 1
				currPos = result.positions[i]
				currChr = result.chromosomes[i]
			#Found start..

			while currChr == chr and currPos < endPos:
				currPos = result.positions[i]
				pos = []
				pos.append(currPos)
				if interactionResults[0]:
					pos.append(result.interactionResults[i])
				currChr = result.chromosomes[i]
				positions.append(pos)
				if result.scores[i] > 10 and pvals[j]:
					scores.append(10.0)
				else:
					scores.append(result.scores[i])
				i += 1

			positionsList.append(positions)
			scoreList.append(scores)
			#Found the end

		startPos = positionsList[0][0]
		endPos = positionsList[0][len(positionsList[0]) - 1]
		for i in range(1, len(positionsList)):
			positions = positionsList[i]
			if positions[0][0] < startPos:
				startPos = positions[0][0]
			if positions[len(positions) - 1][0] > endPos:
				endPos = positions[len(positions) - 1][0]
		posRange = endPos - startPos

		if genes:
			numGeneLines = int(2 + len(genes) / 14)
			pylab.figure(1, figsize=(18, 4 + 0.4 * numGeneLines))
		else:
			pylab.figure(1, figsize=(18, 4))
		if ldSNP:
			r2Values = []
			for i in range(0, len(snps)):
				snp = snps[i]
				pos = positionsList[0][i]
				r2Values.append(self.r2_ld(snp, ldSNP.alleles) * rangeVal + minVal)

			pylab.plot(positionsList[0], r2Values, "y-", label='$r^2$')

		for j in range(2, len(results)):
			minScore = min(results[j].scores)
			maxScore = max(results[j].scores)
			scoreRange = maxScore - minScore
			scores = []
			for score in scoreList[j]:
				scores.append(((score - minScore) / scoreRange) * rangeVal + minVal)
			if resultTypes[j] == "Marg":
				pylab.plot(positionsList[j], scores, "g", label='Marg')
			elif resultTypes[j] == "RF":
				pos = positionsList[j][0]
				score = scores[0]
				pylab.plot([pos, pos], [0, score + 0.03], "c-", label='RF', lw=1.8)
				for k in range(1, len(scores)):
					pos = positionsList[j][k]
					score = scores[k]
					pylab.plot([pos, pos], [0, score + 0.03], "c-", lw=1.8)
			elif resultTypes[j] == "CS":
				pylab.plot(positionsList[j], scores, "b.", label='CS')


		pylab.plot(positionsList[0], scoreList[0], "r.", label='KW')
#		if not binary:
#			pylab.plot(positionsList[1],scoreList[1],"b.",label='Emma')

		pylab.plot([startPos - 0.05 * posRange, endPos + 0.05 * posRange], [6.68, 6.68], "k:")

		if genes:
			numGeneLines = int(2 + len(genes) / 14)
			self.drawGenes(genes, gene_position_cycle=numGeneLines)
			pylab.axis([startPos - 0.05 * posRange, endPos + 0.05 * posRange, minVal - rangeVal * 0.05 - 0.49 * numGeneLines, maxVal + rangeVal * 0.05])
		else:
			pylab.axis([startPos - 0.05 * posRange, endPos + 0.05 * posRange, minVal - rangeVal * 0.05, maxVal + rangeVal * 0.05])

		if ldSNP:
			pylab.title(phenName + ": chromosome " + str(chr) + ", position " + str(ldSNP.position) + ".")
		else:
			pylab.title(phenName + ": chromosome " + str(chr) + ".")

		pylab.subplots_adjust(right=0.98)
		pylab.subplots_adjust(left=0.03)
		pylab.subplots_adjust(bottom=0.15)
		pylab.subplots_adjust(top=0.9)
		pylab.legend(numpoints=2, handlelen=0.005)

		if pdfFile:
			pylab.savefig(pdfFile, format="pdf")
		if pngFile:
			pylab.savefig(pngFile, format="png")
		if not (pdfFile or pngFile):
			pylab.show()

		pylab.close(1)



	def _to01Format_(self, snp):
		all1 = snp[0]
		tSnp = [0] * len(snp)
		for i in range(1, len(snp)):
			allele = snp[i]
			if allele != all1:
				tSnp[i] = 1
		return tSnp


	def r2_ld(self, snp1, snp2):
		tSnp1 = self._to01Format_(snp1)
		tSnp2 = self._to01Format_(snp2)
		delta = 1.0 / float(len(snp1))
		freqs = [0.0] * 4
		for i in xrange(0, len(snp1)):
			val = tSnp1[i] * 2 + tSnp2[i]
			freqs[val] += delta


		f1 = freqs[1] + freqs[3]
		f2 = freqs[2] + freqs[3]
		D = freqs[3] - f1 * f2
		divisor = f1 * f2 * (1 - f1) * (1 - f2)
		if divisor != 0:
			r2 = D * D / divisor
		else:
			r2 = -1
		return r2







	def drawGWPlot(phenotypeID):
		"""
		Draws all the GWA plots for 6 methods.
		"""
		import plotResults
		results = self.results_map[phenotypeID]
		for m_i in range(0, len(results)): #For all methods 
			result = results[m_i]
			print "\nPlotting result", result.name, ":"
			if m_i == 0:
				plotResults.plotResult(result, ylab="KW: -log(p-value)")
			elif m_i == 1:
				plotResults.plotResult(result, ylab="Emma: -log(p-value)")
			elif m_i == 2:
				plotResults.plotResult(result, type="score", ylab="Margarita: ARG score")
			elif m_i == 3:
				plotResults.plotResult(result, type="score", ylab="RF: Importance score")
			elif m_i == 4:
				plotResults.plotResult(result, type="score", ylab="Simple composite rank score")





#def fun1(rp):
#	"""
#	chr5 2.29 - 2.8 M (complex regions)
#	chr5 5.17 - 6.8M (complex regions)
#	chr5 3.5M (myb - emma peak)
#	chr5 6.8 ( Emma peak -HHP1) 
#	chr5 7.78M (Agamous-like)
#	chr5 8.3 (Marg high - other methods not so)
#	chr2 9.58 - 9.62M (SVP and Agamous like)
#	chr1 19.79M (spa4 - nice peak)
#	"""
#	chr = 5
#	startPos = 2250000
#	endPos = 2810000
#	for pi in ft:
#		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
#		if pi != ft[0]:
#			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf")
#		else:
#			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt")
#		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_noGenes"
#		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",plotGenes=False)
#   
# 
#
#def fun3(rp):
#	"""
#	"""
#	chr = 5
#	startPos = 6480000
#	endPos = 7020000
#	for pi in ft:
#		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
#		if pi != ft[0]:
#			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf")
#		else:
#			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt")
#		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_noGenes"
#		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",plotGenes=False)
#
#
#
#def fun4(rp):
#	"""
#	"""
#	chr = 1
#	startPos = 3800000
#	endPos = 4010000
#	for pi in ft:
#		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
#		if pi != ft[0]:
#			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf")
#		else:
#			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt")
#		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_noGenes"
#		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",plotGenes=False)
#
#
#
#def fun5(rp):
#	"""
#	"""
#	chr = 4
#	startPos = 130000
#	endPos = 700000
#	for pi in ft:
#		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
#		if pi != ft[0]:
#			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf")
#		else:
#			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt")
#		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_noGenes"
#		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",plotGenes=False)
#
#
#def fri(rp):
#	"""
#	"""
#	chr = 4
#	startPos = 200000
#	endPos = 350000
#	for pi in ft:
#		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld_ler"
#		if pi != ft[0]:
#			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=268809)
#		else:
#			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt",snpPos=268809)
#		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld_col"
#		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=269962)
#		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
#		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf")
#
#def fri2(rp):
#	"""
#	"""
#	chr = 4
#	startPos = 360000
#	endPos = 510000
#	for pi in ft:
#		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld_fig1"
#		if pi != ft[0]:
#			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=409692)
#		else:
#			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt",snpPos=409692)
#		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld_fig2"
#		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=429928)
#		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld_fig3"
#		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=454542)
#		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
#		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf")
#
def fri_test():
	"""
	"""
	resDir = "/Network/Data/250k/tmp-bvilhjal/"
	chr = 4
	startPos = 0
	#endPos = 700000
	endPos = 1600000

	resFileDir = "/Network/Data/250k/tmp-bvilhjal/lr_results/"
	resultFiles = [("/Users/bjarnivilhjalmsson/Projects/Data/gwas_results/FRI_Emma/LR_logTransform_43_FLC.pvals", "LR", "Full"),
			("/Users/bjarnivilhjalmsson/Projects/Data/gwas_results/FRI_Emma/LR_logTransform_col_43_FLC.pvals", "LR", "Col"),
			("/Users/bjarnivilhjalmsson/Projects/Data/gwas_results/FRI_Emma/LR_logTransform_ler_43_FLC.pvals", "LR", "Ler"),
			("/Users/bjarnivilhjalmsson/Projects/Data/gwas_results/FRI_Emma/LR_logTransform_col_ler_43_FLC.pvals", "LR", "Col_Ler")]
#	resultFiles=[("/Network/Data/250k/tmp-bvilhjal/kw_results/KW_raw_192_43_FLC.pvals","KW","Full"),
#			("/Network/Data/250k/tmp-bvilhjal/kw_results/KW_raw_192_Col_43_FLC.pvals","KW","Col"),
#			("/Network/Data/250k/tmp-bvilhjal/kw_results/KW_raw_192_Ler_43_FLC.pvals","KW","Ler"),
#			("/Network/Data/250k/tmp-bvilhjal/kw_results/KW_raw_192_ColAndLer_43_FLC.pvals","KW","Col_Ler")]
#	resultFiles=[("/Users/bjarnivilhjalmsson/Projects/Data/gwas_results/FRI_Emma/Emma_logTransform_192_43_FLC.pvals","Emma","Full"),
#			("/Users/bjarnivilhjalmsson/Projects/Data/gwas_results/FRI_Emma/Runtype_2_SNP_4_1_4_1750000_Cofactor_4_269962_pheno_43_FLC_transform.tsv","Emma","Col"),
#			("/Users/bjarnivilhjalmsson/Projects/Data/gwas_results/FRI_Emma/Runtype_2_SNP_4_1_4_1750000_Cofactor_4_268809_pheno_43_FLC_transform.tsv","Emma","Ler"),
#			("/Users/bjarnivilhjalmsson/Projects/Data/gwas_results/FRI_Emma/Runtype_2_SNP_4_1_4_1750000_Cofactor_4_268809_4_269962_pheno_43_FLC_transform.tsv","Emma","Col_Ler")]
	results = []

	candGenes = gwaResults.getCandidateGeneList(145) #FIXME

	for (filename, method, name) in resultFiles:
		rt = gwaResults.ResultType(resultType=method, datasetName=name, mafCutoff=15, logTransform=True)
		result = gwaResults.Result(filename, name=method + "_" + name, resultType=rt, phenotypeID=43)
		result.negLogTransform()
		print "Log transformed the p-values"
		result.filterMAF(minMaf=15)
		#result.filterPercentile(0.90)
		results.append(result)


	results_map = {43:results}
	rp = RegionPlotter(phenotypeIndices=[43], results_map=results_map)
	rp.plotRegion(chr, startPos, endPos, 43, pngFile="/Users/bjarnivilhjalmsson/tmp/FRI_large_FLC_LM.png", pdfFile="/Users/bjarnivilhjalmsson/tmp/FRI_large_FLC_LM.pdf",
			plotGenes=False, plotColors=[(0.3, 0.6, 0.6)] * 4, #, [(0.1,0.05,1)]*4 ["#1111bb","#11aa22",(1.0,0.1,0),"#aa11bb"]
			labels=["a", "b", "c", "d"], #["Full data","Col removed","Ler removed", "Col, Ler removed"], #
			markers=[".", ".", ".", "."], highlightGenes=candGenes, alphas=[1] * 4, #[0.8,0.8,0.8,0.8],
			maxLogPval=13) #8.5


def flc_test():
	"""
	"""
	resDir = "/Network/Data/250k/tmp-bvilhjal/"
	chr = 5
	startPos = 2250000
	endPos = 3350000


	resultFiles = [("/Users/bjarni/tmp/FLC_FT_22.pvals", "Emma", "local"), ("/Network/Data/250k/tmp-bvilhjal/emma_results/Emma_newDataset_7_FT_22C.pvals", "Emma", "global")]
	results = []

	for (filename, method, name) in resultFiles:
		rt = gwaResults.ResultType(resultType=method, datasetName=name, mafCutoff=15, logTransform=True)
		result = gwaResults.Result(filename, name=method + "_" + name, resultType=rt, phenotypeID=7)
		result.negLogTransform()
		print "Log transformed the p-values"
		result.filterMAF(minMaf=15)
		results.append(result)

	results_map = {7:results}
	rp = RegionPlotter(phenotypeIndices=[7], results_map=results_map)
	rp.plotRegion(chr, startPos, endPos, 7)


#
#def svp(rp):
#	"""
#	"""
#	chr = 2
#	startPos = 9530000
#	endPos = 9680000
#	for pi in ft:
#		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld_fig1"
#		if pi != ft[0]:
#			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=9588500)
#		else:
#			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt",snpPos=9588500)
#		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld_fig2"
#		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=9611600)
#		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
#		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf")
#
#
#
#def ie1(rp):
#	"""
#	Interesting example
#	"""
#	chr = 3
#	startPos = 18875000
#	endPos = 19025000
#	for pi in ft:
#		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld"
#		if pi != ft[0]:
#			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=18923922)
#		else:
#			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt",snpPos=18923922)
#		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
#		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf")
#
#
#def ie2(rp):
#	"""
#	Interesting example
#	"""
#	chr = 5
#	startPos = 18560000
#	endPos = 18660000
#	for pi in ft:
#		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld_fig1"
#		if pi != ft[0]:
#			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=18625634)
#		else:
#			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt",snpPos=18625634)
#		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld_fig2"
#		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=18582030)
#		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
#		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf")
#
#
#
#
#def ie3(rp):
#	"""
#	Interesting example
#	"""
#	chr = 1
#	startPos = 6340000
#	endPos = 6410000
#	for pi in ft:
#		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld"
#		if pi != ft[0]:
#			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=6369765)
#		else:
#			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt",snpPos=6369765)
#		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
#		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf")
#
#
#def lesioning(rp):
#	"""
#	Interesting example
#	"""
#	chr = 4
#	startPos = 8220000
#	endPos = 8360000
#	pi=77
#	filename = resDir+"Lesioning_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
#	rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt",binary=True)
#	tairFileName = resDir+"Lesioning_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
#	filename = resDir+"Lesioning_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld_fig1"
#	rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=8274503,binary=True)
#	filename = resDir+"Lesioning_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld_fig2"
#	rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=8297531,binary=True)
#
# 
#
#
#def avr(rp):
#	"""
#	Interesting example
#	"""
#	chr = 3
#	startPos = 2150000
#	endPos = 2300000
#	pi=33
#	phenName = "avrB"
#	filename = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
#	rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt",binary=True)
#	tairFileName = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
#	filename = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld"
#	rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=2227823,binary=True)
#
#
#
#
#
#
#def rp2(rp):
#	"""
#	Interesting example
#	"""
#	chr = 4
#	startPos = 13190000
#	endPos = 13270000
#	pi=34
#	phenName = "Rps2"
#	filename = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
#	rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt",binary=True)
#	tairFileName = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
#	filename = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld"
#	rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=13226419,binary=True)
#	
#
#
#
#def rp2_b(rp):
#	"""
#	Interesting example
#	"""
#	chr = 4
#	startPos = 9760000
#	endPos = 9840000
#	pi=34
#	phenName = "Rps2"
#	filename = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
#	rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt",binary=True)
#	tairFileName = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
#	filename = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld"
#	rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=9803559,binary=True)
#	
#
#def hyp(rp):
#	"""
#	Interesting example
#	"""
#	chr = 5
#	startPos = 5015000
#	endPos = 5120000
#	pi=182
#	phenName = "Hypocot_length"
#	filename = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
#	rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt")
#	tairFileName = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
#	filename = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld"
#	rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=5069348)
#	
#
#
#
#
#def ant(rp):
#	"""
#	Interesting example
#	"""
#	chr = 2
#	startPos = 1890000
#	endPos = 2000000
#	pi=170
#	phenName = "Antho_10"
#	filename = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
#	rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt")
#	tairFileName = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
#	filename = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld"
#	rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=1937020)
#	
#
#
#
#def na_23(rp):
#	"""
#	Interesting example
#	"""
#	chr = 4
#	startPos = 6350000
#	endPos = 6470000
#	for pi in [16,85]:
#		phenName = "Na_23_Soil_L"
#		filename = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
#		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt")
#		tairFileName = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
#		filename = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld"
#		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=6391200)
#	
#
#
#def na_23(rp):
#	"""
#	Interesting example
#	"""
#	chr = 4
#	startPos = 6350000
#	endPos = 6470000
#	for pi in [16,85]:
#		phenName = "Na_23_Soil_L"
#		filename = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
#		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt")
#		tairFileName = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
#		filename = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld"
#		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=6391200)
#	
#
#
#def _test_():
#	rp = RegionPlotter([1,2,3,4,5,6,7])
#	rp.plotRegion(5, 18550000, 18670000,1,pngFile='/Users/bjarni/tmp/DOG_region_1.png',labels=["KW","KW SR","Emma","Emma SR"],plotColors=['b','g','r','c'],plotOrder=[2,4,1,3])
#	rp.plotRegion(5, 18550000, 18670000,2,pngFile='/Users/bjarni/tmp/DOG_region_2.png',labels=["KW","KW SR","Emma","Emma SR"],plotColors=['b','g','r','c'],plotOrder=[2,4,1,3])
#	rp.plotRegion(5, 18550000, 18670000,3,pngFile='/Users/bjarni/tmp/DOG_region_3.png',labels=["KW","KW SR","Emma","Emma SR"],plotColors=['b','g','r','c'],plotOrder=[2,4,1,3])
#	rp.plotRegion(5, 18550000, 18670000,4,pngFile='/Users/bjarni/tmp/DOG_region_4.png',labels=["KW","KW SR","Emma","Emma SR"],plotColors=['b','g','r','c'],plotOrder=[2,4,1,3])
#	rp.plotRegion(5, 18550000, 18670000,5,pngFile='/Users/bjarni/tmp/DOG_region_5.png',labels=["KW","KW SR","Emma","Emma SR"],plotColors=['b','g','r','c'],plotOrder=[2,4,1,3])
#	rp.plotRegion(5, 18550000, 18670000,6,pngFile='/Users/bjarni/tmp/DOG_region_6.png',labels=["KW","KW SR","Emma","Emma SR"],plotColors=['b','g','r','c'],plotOrder=[2,4,1,3])
#	rp.plotRegion(5, 18550000, 18670000,7,pngFile='/Users/bjarni/tmp/DOG_region_7.png',labels=["KW","KW SR","Emma","Emma SR"],plotColors=['b','g','r','c'],plotOrder=[2,4,1,3])
#
#
#def _ft_genes_():
#	"""
#	gene name	phenotype	Phenotype label at top of plot	gene id	chr	snp position	 	 
#	2 	ATH1	ft_10	ft_10	At4g32980	4	15918498	 	 
#	3 	SVP	1_LD	LD	At2g22540	2	9588685	 	 
#	4 	FRI	FRI	FRI	At4g00650	4	268809	 	 
#	5 	FPF1	ft_16	FT16	A5g24860	5	8540225	 	 
#	6 	YAP169	ft_16	FT16	At5g07200	5	2241523	 	 
#	7 	RAV-1	3_SD	SD	At1g13260	1	4540454	 	 
#	8 	SPL4	ft_16	FT16	At1g53160	1	19811699	 	 
#	9 	CRP	FLC	FLC	At4g00450	4	200026	 	 
#	10 	SPA2	2_LDV	LDV	At4g11110	4	6776344	 	 
#	11 	SPA4	ft_16	FT16	At1g53090	1	19789536	 	 
#	12 	SIMILAR TO VRN1	2_LDV	LDV	At4g33280	4	16048992	 	 
#	13 	DOG1	LN_10	LN10	At5g45830	5	18607728	 	 
#	14 	DFL2	ft_22	FT22	At4g03400	4	1498770	 	 
#	"""
#	
#	phenotypeFile = "/Network/Data/250k/dataFreeze_011209/phenotypes_all_raw_042109.tsv"
#	phed = phenotypeData.readPhenotypeFile(phenotypeFile)	
#	snpsDataFile="/Network/Data/250k/dataFreeze_011209/250K_192_043009.csv"
#	snpsds = dataParsers.parseCSVData(snpsDataFile, format=1, deliminator=",")
#
#	geneIDs = []
#	genePositions = []
#	geneChrs = []
#	geneNames = []
#	phenIDs = []
#	f = open("/Network/Data/250k/tmp-bvilhjal/Region_plots/Promising_associations.csv","r")
#	lines = f.readlines()
#	for line in lines[1:]:
#		line_lst = line.split(",")
#		geneNames.append(line_lst[0])
#		phenIDs.append(int(line_lst[1]))
#		geneIDs.append(line_lst[3].upper())
#		geneChrs.append(int(line_lst[4]))
#		genePositions.append(int(line_lst[5]))
#	f.close()
#	
#	pi_list = list(set(phenIDs))
#	print "Number of phenotypes:",len(pi_list)
#	rp = RegionPlotter(pi_list,snpsds=snpsds)
#
#	for i in range(0,len(geneNames)):
#		phenName = phed.getPhenotypeName(phenIDs[i])
#		windowSize = 1000000
#		filename = "/Network/Data/250k/tmp-bvilhjal/Region_plots/reg_plot_w"+str(round(windowSize/1000000.0,2))+"_"+geneNames[i]+"_"+phenName
#		pngFile = filename+".png"
#		pdfFile = filename+".pdf"
#		startPos = genePositions[i]-windowSize/2
#		endPos = genePositions[i]+windowSize/2
#		if startPos < 0:
#			startPos = 0
#			endPos = windowSize
#		elif endPos>snpsds[geneChrs[i]-1].positions[-1]:
#			startPos = snpsds[geneChrs[i]-1].positions[-1] - windowSize
#			endPos = snpsds[geneChrs[i]-1].positions[-1]
#			
#		methodLabel = "Wilcoxon"
#		if phed.isBinary(phenIDs[i]):
#			methodLabel = "Fisher's"
#		rp.plotRegion(geneChrs[i], startPos, endPos,phenIDs[i],pngFile=pngFile,pdfFile=pdfFile,labels=[methodLabel,"EMMA",],plotColors=['b','r'],highlightGenes=[geneIDs[i]])
#	
#	#rp = RegionPlotter([1,2,3,5,6,7,43,44,47,59,80,82])
##	rp = RegionPlotter([184,43])
##	rp.plotRegion(2, 12476632, 13476632,184,pdfFile='/Users/bjarni/tmp/ETC2_TCL2_TCL1_184_Trich_avg_JA.pdf',pngFile='/Users/bjarni/tmp/ETC2_TCL2_TCL1_184_Trich_avg_JA.png',labels=["Wilcoxon","EMMA",],plotColors=['b','r'],highlightGenes=["AT2G30420","AT2G30424","AT2G30432"])
##	rp.plotRegion(4, 0, 1000000,43,pdfFile='/Users/bjarni/tmp/CRP_FRI_ETC3_43_FLC.pdf',pngFile='/Users/bjarni/tmp/CRP_FRI_ETC3_43_FLC.png',labels=["Wilcoxon","EMMA",],plotColors=['b','r'],highlightGenes=["AT4G00450","AT4G00650","AT4G01060"])
##	rp.plotRegion(4, 168809, 418809,44,pngFile='/Users/bjarni/tmp/FRI_44_FRI.png',labels=["KW","Emma",],plotColors=['b','r'],highlightGenes=["AT4G00650"])
##	rp.plotRegion(5, 8420225, 8660225,6,pngFile='/Users/bjarni/tmp/FPF1_6_FT16C.png',labels=["KW","Emma",],plotColors=['b','r'],highlightGenes=["AT5G24860"])
##	rp.plotRegion(5, 1941523, 2641523,6,pngFile='/Users/bjarni/tmp/YAP169_6_FT16C.png',labels=["KW","Emma",],plotColors=['b','r'],highlightGenes=["AT5G07200"])
##	rp.plotRegion(1, 4400454, 4680454,3,pngFile='/Users/bjarni/tmp/RAV-1_3_SD.png',labels=["KW","Emma",],plotColors=['b','r'],highlightGenes=["AT1G13260"])
##	rp.plotRegion(1, 19651699, 19961699,6,pngFile='/Users/bjarni/tmp/SPL4_6_FT16C.png',labels=["KW","Emma",],plotColors=['b','r'],highlightGenes=["AT1G53160"])
#
##	rp.plotRegion(4, 15418498, 16418498,5,pngFile='/Users/bjarni/tmp/ATH1_5_FT10C.png',labels=["Wilcoxon","Emma",],plotColors=['b','r'],highlightGenes=["AT4G32980"])
##	rp.plotRegion(2, 9088685, 10088685,1,pngFile='/Users/bjarni/tmp/SVP_1_LD.png',labels=["Wilcoxon","Emma",],plotColors=['b','r'],highlightGenes=["AT2G22540"])
##	rp.plotRegion(4, 0, 1000000,44,pngFile='/Users/bjarni/tmp/FRI_44_FRI.png',labels=["Wilcoxon","Emma",],plotColors=['b','r'],highlightGenes=["AT4G00650"])
##	rp.plotRegion(5, 8040225, 9040225,6,pngFile='/Users/bjarni/tmp/FPF1_6_FT16C.png',labels=["Wilcoxon","Emma",],plotColors=['b','r'],highlightGenes=["AT5G24860"])
##	rp.plotRegion(5, 1741523, 2741523,6,pngFile='/Users/bjarni/tmp/YAP169_6_FT16C.png',labels=["Wilcoxon","Emma",],plotColors=['b','r'],highlightGenes=["AT5G07200"])
##	rp.plotRegion(1, 4040454, 5040454,3,pngFile='/Users/bjarni/tmp/RAV-1_3_SD.png',labels=["Wilcoxon","Emma",],plotColors=['b','r'],highlightGenes=["AT1G13260"])
##	rp.plotRegion(1, 19311699, 20311699,6,pngFile='/Users/bjarni/tmp/SPL4_6_FT16C.png',labels=["Wilcoxon","Emma",],plotColors=['b','r'],highlightGenes=["AT1G53160"])
##	rp.plotRegion(4, 0, 1000000,1,pngFile='/Users/bjarni/tmp/CRP_1_LD.png',labels=["Wilcoxon","Emma",],plotColors=['b','r'],highlightGenes=["AT4G00450"])
##	rp.plotRegion(4, 6276344, 7276344,2,pngFile='/Users/bjarni/tmp/SPA2_2_LDV.png',labels=["Wilcoxon","Emma",],plotColors=['b','r'],highlightGenes=["AT4G11110"])
##	rp.plotRegion(1, 19289536, 20289536,6,pngFile='/Users/bjarni/tmp/SPA4_6_FT16C.png',labels=["Wilcoxon","Emma",],plotColors=['b','r'],highlightGenes=["AT1G53090"])
##	rp.plotRegion(4, 15548992, 16548992,2,pngFile='/Users/bjarni/tmp/SIMILAR_TO_VRN1_2_LDV.png',labels=["Wilcoxon","Emma",],plotColors=['b','r'],highlightGenes=["AT4G33280"])
##	rp.plotRegion(5, 18107728, 19107728,80,pngFile='/Users/bjarni/tmp/DOG1_80_LN10.png',labels=["Wilcoxon","Emma",],plotColors=['b','r'],highlightGenes=["AT5G45830"])
##	rp.plotRegion(4, 998770, 1998770,7,pngFile='/Users/bjarni/tmp/DFL2_7_FT22C.png',labels=["Wilcoxon","Emma",],plotColors=['b','r'],highlightGenes=["AT4G03400"])
##	rp.plotRegion(5, 2678615, 3678615,82,pngFile='/Users/bjarni/tmp/FLC_82_LN22C.png',labels=["Wilcoxon","Emma",],plotColors=['b','r'],highlightGenes=["AT5G10140"])


def plot_simple_region(scores_list, positions_list, chromosome=0, caption='', plot_genes=True,
			png_file=None, pdf_file=None, tair_file=None):
	reg_res_list = []
	for scores, positions in zip(scores_list, positions_list):
		num_snps = len(scores)
		reg_res_list.append(gwaResults.Result(scores=scores, chromosomes=[chromosome] * num_snps, positions=positions))
	plotter = RegionPlotter()
	plotter.plot_small_result(reg_res_list, png_file=png_file, pdf_file=pdf_file, caption=caption,
				tair_file=tair_file, plot_genes=plot_genes)


#def plot_simple_region(snps,positions,scores):
#	def _plot_simple_result_(self,scores,positions,caption="",pdf_file=None,png_file=None,plot_genestair_genes=None,
#				highlight_genes=None,max_val=None,min_val=None,bonferroni_level=None,color=None):
#		"""
#		A helper function actually drawing the simple plot.
#		"""
#		
#		print "Plotting",len(scores),"SNPs."
#		
#		x_max = max(positions)
#		x_min = min(positions)
#		x_range = x_max-x_min
#		y_min = min(scores)
#		y_max = max(scores)
#		y_range = y_max-y_min
#		x_buffer = 0.05
#		y_buffer = 0.05
#					
#		
#		#plot box boundaries [mix_x,max_x,min_y,max_y]
#
#
#		if tair_genes:
#			pylab.figure(figsize=(8,6))			
#			plot_range = [min_x-x_buffer*x_range, max_x+x_buffer*x_range, 
#				y_min-y_range*(y_buffer+0.25), y_max+y_range*y_buffer] 
#		else:
#			pylab.figure(figsize=(8,4))	
#			plot_range = [min_x-x_buffer*x_range, max_x+x_buffer*x_range, 
#				y_min-y_range*y_buffer, y_max+y_range*y_buffer] 
#			
#		pylab.axes([0.06, 0.10, 0.90, 0.82])
#		#for (position,score,color) in zip(positions,scores,colors):
#		#	print "drawing.."
#		#	pylab.plot(position,score,color=color,marker=".",ls="",alpha=0.8)
#		pylab.plot(positions,scores,marker=".",ls="",alpha=0.8)
#		#pylab.text(endPos-0.025*posRange,maxVal-0.15*rangeVal,l,size="x-large")							
#		#pylab.legend(numpoints=1,handlelen=0.005,prop=fontProp)
#		if bonferroni_level:
#			pylab.plot([plot_range[0],plot_range[1]],[bonferroni_level,bonferroni_level],"k-.")
#		if tair_genes:
#			num_rows = max(3,min(6,1+len(tair_genes)/10))
#			print "Drawing TAIR genes"
#			drawGenes(tair_genes,gene_position_cycle=num_rows,rangeVal=y_range*1.5)
#		pylab.axis(plot_range)
#		if pdf_file:
#			pylab.savefig(pdf_file,format="pdf")
#		if png_file:
#			pylab.savefig(png_file,format="png",dpi=300)
#		if not (pdf_file or png_file):
#			pylab.show()
#		pylab.close()
#		
#		
#	



if __name__ == '__main__':
	fri_test()
	#_ft_genes_()
	pass


