"""
Usage: gwa.py [OPTIONS] 

Option:

	-h					Show this help
	-d ...					Filter SNPs randomly (speed-up for debugging)

	-o ...					ID string, used for output files.
	-i ...					The phenotype IDs, to be run. 

	-t ...					What data set is used.  Default is 54.
	-f ...					Load a specific data file, e.g. for heteroplasmy.
	-r ...					Phenotype file, if left out then phenotypes are retireved from the DB 
						(transformed values).
	-k ...					Specify the file containing the kinship matrix.  
						(otherwise default file is used or it's generated.)

	-a ...					Apply specific methods, otherwise all available are applied:
						lm, emmax, rf, etc.
	-b ...				 	Apply a transformation to the data, default is none, other possibilities are 
						log_trans, box_cox_lambda (where lambda is a number)
	-c ...					Should phenotype outliers be removed.  0 (no fence) is the default, 
						else the outlier fence is given in IQRs. (An int is required).
												 
	-u					Use existing results when possible, to speed up.  
						(Otherwise existing files are overwritten)

	-s ...					How many SNPs should be used for the second step.  (Picked by score rank)
						Default is 1000.

	-n					Adds the result(s) to the DB.		
	--comment=...				Comment for DB. (Only applicable if result is added to DB.)	
	--no_phenotype_ids			Phenotypes don't have DB id as an prefix in their names.  
						(The phenotype index should be used instead.)	

	--region_plots=...			Include region plots for the top N (given num) peaks.

	-m ...					MAC threshold which is used for LM, EMMAX, etc.  Default is 15.
	--data_format=...			What type of data should it aim for, binary (default), int, float, etc.
	
	
	#ONLY APPLICABLE FOR CLUSTER RUNS
	-p ...					Run mapping methods on the cluster with standard parameters.  The argument is used for runid 
						as well as output files.  Note that if this option is used then no output file should be specified.
						If phenotype index is missing, then all phenotypes are used.
	-q ...					Request memory (on cluster), otherwise use defaults 4GB for Kruskal-Wallis, 12GB for EMMAX.
	-l ...			 		Request time limit (on cluster), otherwise use defaults
	--only_add_2_db				Does not submit jobs, but only adds available result files to DB. (hack for usc hpc)
	
	
	
Examples:
~/gwas_data$ python gwa.py -o test -i 1,5 -a kw,emmax -r ~/Projects/Data/phenotypes/phen_raw_092910.tsv 
Description:
  Applies various GWA methods to to phenotypes.

  Methods include: EMMAX, LM, RF, etc.
  
  If PHENOTYPE_DATA_FILE is left out, then papaya DB is searched.  
  A PHENOTYPE_ID is required.
    
"""
import sys, getopt, traceback
import os

import phenotypeData
import dataParsers
import snpsdata
import gwaResults
import util
import warnings
import multiprocessing as mp
import time
import cPickle

import linear_models as lm
import scipy as sp
from env import *

import gwa


transformation_method_dict = {
			'none':1,
			'log_trans':2,
			'box_cox':3,
			}


analysis_methods_dict = {"kw":1,
			 "emmax":32,
			 }




def _parse_pids_(pid_arg_str):
	t_pids = pid_arg_str.split(',')
	pids = []
	for s in t_pids:
		if '-' in s:
			pid_range = map(int, s.split('-'))
		        for pid in range(pid_range[0], pid_range[1] + 1):
		        	pids.append(pid)
		else:
			pids.append(int(s))
	return pids


def parse_parameters():
	'Parse the parameters into a dict, etc.'
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)

	long_options_list = ["comment=", 'no_phenotype_ids', 'region_plots=', 'only_add_2_db', 'data_format=']
	try:
		opts, args = getopt.getopt(sys.argv[1:], "o:i:p:a:b:c:d:ef:t:r:s:k:nm:q:l:hu", long_options_list)

	except:
		traceback.print_exc()
		print __doc__
		sys.exit(2)


	p_dict = {'run_id':'donald_duck', 'parallel':None, 'add_to_db':False, 'comment':'', 'mem_req':'1800mb',
		'call_method_id':54, 'walltime_req':'12:00:00',
		'specific_methods':['lm', 'emmax'], 'specific_transformations':['none'],
		'remove_outliers':0, 'kinship_file':None, 'analysis_plots':False, 'use_existing_results':False,
		'region_plots':0, 'debug_filter':1, 'phen_file':None,
		'no_phenotype_ids':False, 'only_add_2_db':False, 'mac_threshold':15, 'data_file':None,
		'data_format':'binary', 'second_step_number':1000}


	for opt, arg in opts:
		if opt in ("-h"):
			print __doc__
			return
		elif opt in ('-i'): p_dict['pids'] = _parse_pids_(arg)
		elif opt in ('-o'): p_dict['run_id'] = arg
		elif opt in ('-p'): p_dict['parallel'] = arg
		elif opt in ('-t'): p_dict['call_method_id'] = int(arg)
		elif opt in ('-n'): p_dict['add_to_db'] = True
		elif opt in ('-m'): p_dict['mac_threshold'] = int(arg)
		elif opt in ('-f'): p_dict['data_file'] = arg
		elif opt in ('-r'): p_dict['phen_file'] = arg
		elif opt in ('-k'): p_dict['kinship_file'] = arg
		elif opt in ('-a'): p_dict['specific_methods'] = arg.split(',')
		elif opt in ("-b"): p_dict['specific_transformations'] = arg.split(',')
		elif opt in ("-c"): p_dict['remove_outliers'] = int(arg)
		elif opt in ("-d"): p_dict['debug_filter'] = float(arg)
		elif opt in ('-s'): p_dict['second_step_number'] = int(arg)
		elif opt in ('-u'): p_dict['use_existing_results'] = True
		elif opt in ('-e'): p_dict['analysis_plots'] = True
		elif opt in ('-q'): p_dict['mem_req'] = arg
		elif opt in ('-l'): p_dict['walltime_req'] = arg
		elif opt in ("--comment"): p_dict['comment'] = arg
		elif opt in ("--no_phenotype_ids"): p_dict['no_phenotype_ids'] = True
		elif opt in ("--region_plots"): p_dict['region_plots'] = int(arg)
		elif opt in ("--only_add_2_db"): p_dict['only_add_2_db'] = True
		elif opt in ("--data_format"): p_dict['data_format'] = arg
		else:
			print "Unkown option!!\n"
			print __doc__
			sys.exit(2)

	return p_dict, args





def _get_file_prefix_(id, p_i, phenotype_name, mapping_method=None, trans_method=None,
		remove_outliers=None, second_step_number=None):
	prefix = env['results_dir'] + id + "_pid" + str(p_i) + "_" + phenotype_name
	if mapping_method:
		prefix += "_" + mapping_method
	if trans_method:
		prefix += "_" + trans_method
	if remove_outliers:
		prefix += '_no' + str(remove_outliers)
	if second_step_number:
		prefix += '_sn' + str(second_step_number)
	return prefix



def run_parallel(p_i, phed, p_dict, mapping_method="analysis", trans_method='none'):
	"""
	If no mapping_method, then analysis run is set up.
	"""

	if not p_i in phed.phenIds:
		print "Phenotype ID not found:%i" % p_i
		return
	if phed.isBinary(p_i):
		if trans_method != 'none':
			return
		elif mapping_method in ["kw"]:
			return
		elif p_dict['analysis_plots']:
			specific_methods = [ 'lm', 'emmax']
	else:
		if mapping_method in ["ft"]:
			return
		elif p_dict['analysis_plots']:
			specific_methods = ['lm', 'emmax']

	phenotype_name = phed.getPhenotypeName(p_i)
	#Cluster specific parameters
	print "Setting up a gwa run for phenotype:%s, pid=%d, using method:%s, with transformation as:%s"\
		% (phenotype_name, p_i, mapping_method, trans_method)
	run_id = p_dict['parallel']

	shstr = "#!/bin/csh\n"
	shstr += "#PBS -l walltime=%s \n" % p_dict['wall_time_req']
	shstr += "#PBS -l mem=%s \n" % p_dict['mem_req']
	shstr += "#PBS -q cmb\n"

	job_id = '%d_%s_%s_%s' % (p_i, mapping_method, trans_method, p_dict['parallel'])
	shstr += "#PBS -N p%s \n" % job_id
	shstr += "(python %s gwa.py -o %s " % (env['script_dir'], run_id)
	shstr += "-i %d -a %s -b %s -c %d -d %n -t %d -m %d -s %d " % (p_i, p_dict['mapping_method'],
					p_dict['trans_method'], p_dict['remove_outliers'], p_dict['debug_filter'],
					p_dict['call_method_id'], p_dict['mac_threshold'], p_dict['second_step_number'])

	shstr += "--region_plots=%d --data_format=%s " % \
		(p_dict['region_plots'], p_dict['data_format'])


	if p_dict['use_existing_results']: shstr += "--use_existing_results  "
	if p_dict['analysis_plots']: shstr += "-e "
	if p_dict['phen_file']: shstr += "-r %s " % p_dict['phen_file']
	if p_dict['kinship_file']: shstr += "-k %s " % p_dict['kinship_file']
	if p_dict['no_phenotype_ids']: shstr += "--no_phenotype_ids "
	if p_dict['add_to_db']: shstr += "-n "
	if p_dict['comment']: shstr += "--comment=%s " % comment

	shstr += "> " + run_id + "_job" + ".out) >& " + run_id + "_job" + ".err\n"
	#print '\n',shstr,'\n'
	script_file_name = p_dict['parallel'] + ".sh"
	f = open(script_file_name, 'w')
	f.write(shstr)
	f.close()

	#Execute qsub script
	os.system("qsub " + script_file_name)




def map_phenotype(p_i, phed, snps_data_file, mapping_method, trans_method, p_dict):
	phenotype_name = phed.getPhenotypeName(p_i)
	phen_is_binary = phed.isBinary(p_i)
	file_prefix = _get_file_prefix_(p_dict['run_id'], p_i, phed.getPhenotypeName(p_i),
				mapping_method, trans_method, p_dict['remove_outliers'])
	result_name = "%s_%s_%s" % (phenotype_name, mapping_method, trans_method)

	res = None
	sd = dataParsers.parse_snp_data(snps_data_file , format=p_dict['data_format'], filter=p_dict['debug_filter'])
	num_outliers = gwa.prepare_data(sd, phed, p_i, trans_method, p_dict['remove_outliers'])
	if p_dict['remove_outliers']:
		assert num_outliers != 0, "No outliers were removed, so it makes no sense to go on and perform GWA."

	phen_vals = phed.getPhenVals(p_i)
	snps = sd.getSnps()
	if mapping_method in ['emmax']:
		#Load genotype file (in binary format)
		sys.stdout.write("Retrieving the Kinship matrix K.\n")
		sys.stdout.flush()
		k_file = env['data_dir'] + "kinship_matrix_cm" + str(p_dict['call_method_id']) + ".pickled"
		kinship_file = p_dict['kinship_file']
		if not kinship_file and os.path.isfile(k_file): #Check if corresponding call_method_file is available
			kinship_file = k_file
		if kinship_file:   #Kinship file was somehow supplied..
			print 'Loading supplied kinship'
			k = lm.load_kinship_from_file(kinship_file, sd.accessions)
		else:
			print "No kinship file was found.  Generating kinship file:", k_file
			sd = dataParsers.parse_snp_data(snps_data_file , format=p_dict['data_format'])
			snps = sd.getSnps()
			k_accessions = sd.accessions[:]
			if p_dict['debug_filter']:
				import random
				snps = random.sample(snps, int(p_dict['debug_filter'] * len(snps)))
			k = lm.calc_kinship(snps)
			f = open(k_file, 'w')
			cPickle.dump([k, sd.accessions], f)
			f.close()
			num_outliers = gwa.prepare_data(sd, phed, p_i, trans_method, p_dict['remove_outliers'])
			k = lm.filter_k_for_accessions(k, k_accessions, sd.accessions)
		sys.stdout.flush()
		sys.stdout.write("Done!\n")

	if p_dict['remove_outliers']:
		assert num_outliers != 0, "No outliers were removed, so it makes no sense to go on and perform GWA."


	#Check whether result already exists.
	if p_dict['use_existing_results']:
		print "\nChecking for existing results."
		result_file = file_prefix + ".pvals"
		if os.path.isfile(result_file):
			res = gwaResults.Result(result_file=result_file, name=result_name, snps=snps)
			pvals = True
		else:
			result_file = file_prefix + ".scores"
			if os.path.isfile(result_file):
				res = gwaResults.Result(result_file=result_file, name=result_name, snps=snps)
				pvals = False
		if res:
			print "Found existing results.. (%s)" % (result_file)
		sys.stdout.flush()


	if not res: #If results weren't found in a file... then do GWA.

		sys.stdout.write("Finished loading and handling data!\n")

		print "FIRST STEP: Applying %s to data. " % (mapping_method)
		sys.stdout.flush()
		kwargs = {}
		additional_columns = []
		if mapping_method in ['emmax']:
			res = lm.emmax(snps, phen_vals, k)
		elif mapping_method in ['lm']:
			res = lm.linear_model(snps, phen_vals)
		else:
			print "Mapping method", mapping_method, 'was not found.'
			sys.exit(2)

		if mapping_method in ['lm', 'emmax']:
			kwargs['genotype_var_perc'] = res['var_perc']
			betas = map(list, zip(*res['betas']))
			kwargs['beta0'] = betas[0]
			kwargs['beta1'] = betas[1]
			additional_columns.append('genotype_var_perc')
			additional_columns.append('beta0')
			additional_columns.append('beta1')
			pvals = res['ps']
			sys.stdout.write("Done!\n")
			sys.stdout.flush()



		kwargs['correlations'] = calc_correlations(snps, phen_vals)
		additional_columns.append('correlations')

		res = gwaResults.Result(scores=pvals, snps_data=sd, name=result_name, **kwargs)

		if mapping_method in ["emmax", 'lm']:
		 	result_file = file_prefix + ".pvals"
		else:
		 	result_file = file_prefix + ".scores"
		res.write_to_file(result_file, additional_columns)

		print "Generating a GW plot."
		sys.stdout.flush()
		png_file = file_prefix + "_gwa_plot.png"
		#png_file_max30 = file_prefix+"_gwa_plot_max30.png"
		if mapping_method in ['lm', "emmax"]:
			res.neg_log_trans()
			if mapping_method in ["kw", "ft"]:# or p_dict['data_format'] != 'binary':
				#res.plot_manhattan(png_file=png_file_max30,percentile=90,type="pvals",ylab="$-$log$_{10}(p)$", 
				#	       plot_bonferroni=True,max_score=30)
				res.plot_manhattan(png_file=png_file, percentile=90, type="pvals", ylab="$-$log$_{10}(p)$",
					       plot_bonferroni=True)
			else:
				if res.filter_attr("mafs", p_dict['mac_threshold']) > 0:
					#res.plot_manhattan(png_file=png_file_max30,percentile=90,type="pvals",ylab="$-$log$_{10}(p)$", 
					#	       plot_bonferroni=True,max_score=30)				
					res.plot_manhattan(png_file=png_file, percentile=90, type="pvals", ylab="$-$log$_{10}(p)$",
						       plot_bonferroni=True)
		else:
			pass

		print "plotting histogram"
		hist_file_prefix = _get_file_prefix_(p_dict['run_id'], p_i, phenotype_name, trans_method, p_dict['remove_outliers'])
		hist_png_file = hist_file_prefix + "_hist.png"
		phed.plot_histogram(p_i, pngFile=hist_png_file)
	else:
		res.neg_log_trans()
		assert res.filter_attr("mafs", p_dict['mac_threshold']), 'All SNPs have MAC smaller than threshold'


	print "SECOND STEP:"
	res.filter_top_snps(p_dict['second_step_number'])
	snps = res.snps
	positions = res.positions
	chromosomes = res.chromosomes
	#Checking res_file exists
	file_prefix = _get_file_prefix_(p_dict['run_id'], p_i, phed.getPhenotypeName(p_i),
				mapping_method, trans_method, p_dict['remove_outliers'], p_dict['second_step_number'])
	res_file = file_prefix + '_res.cpickled'
	if p_dict['use_existing_results'] and os.path.isfile(res_file):
			print 'Found existing results for the second step... loading.'
			f = open(res_file, 'rb')
			second_res = cPickle.load(f)
			f.close()
	else:
		if mapping_method == 'lm':
			second_res = lm.linear_model_two_snps(snps, phen_vals)
		if mapping_method == 'emmax':
			second_res = lm.emmax_two_snps(snps, phen_vals, k)

		#Pickling results..
		print 'Saving results as pickled file:', res_file
		f = open(res_file, 'wb')
		cPickle.dump(second_res, f, protocol=2)
		f.close()



	#Plotting second step plots:
	score_array = -sp.log10(second_res['ps'])
	p3_score_array = -sp.log10(second_res['p3_ps'])
	p4_score_array = -sp.log10(second_res['p4_ps'])
	import plotResults as pr
	pr.plot_snp_pair_result(chromosomes, positions, score_array, file_prefix + '_scatter')
	pr.plot_snp_pair_result(chromosomes, positions, p3_score_array, file_prefix + '_p3_scatter')
	pr.plot_snp_pair_result(chromosomes, positions, p4_score_array, file_prefix + '_p4_scatter')



	if p_dict['region_plots']:
		import regionPlotter as rp
		regions_results = res.get_top_region_results(p_dict['region_plots'])
		plotter = rp.RegionPlotter()
		print "Starting region plots..."
		for reg_res in regions_results:
			chromosome = reg_res.chromosomes[0]
			caption = phenotype_name + "_c" + str(chromosome) + "_" + mapping_method
			png_file = file_prefix + "_reg_plot_c" + str(chromosome) + "_s" + str(reg_res.positions[0]) \
				+ "_e" + str(reg_res.positions[-1]) + ".png"
			tair_file = file_prefix + "_reg_plot_c" + str(chromosome) + "_s" + str(reg_res.positions[0]) \
				+ "_e" + str(reg_res.positions[-1]) + "_tair_info.txt"
			plotter.plot_small_result([reg_res], png_file=png_file, highlight_gene_ids=tair_ids,
						  caption=caption, tair_file=tair_file)

			#Plot Box-plot
			png_file = file_prefix + "_reg_plot_c" + str(chromosome) + "_s" + str(reg_res.positions[0]) \
				+ "_e" + str(reg_res.positions[-1]) + "_box_plot.png"
			(marker, score, chromosome, pos) = reg_res.get_max_snp()
			marker_accessions = sd.accessions
			phed.plot_marker_box_plot(p_i, marker=marker, marker_accessions=marker_accessions, \
						png_file=png_file, title="c" + str(chromosome) + "_p" + str(pos), \
						marker_score=score, marker_missing_val=sd.missing_val)




def _run_():

	p_dict, args = parse_parameters()
	print "GWA runs are being set up with the following parameters:"
	for k, v in p_dict.iteritems(): print k + ': ' + str(v)
	print ''

	#Load phenotype file
	if p_dict['phen_file']:
		print 'Loading phenotypes from file.'
		phed = phenotypeData.readPhenotypeFile(p_dict['phen_file'],
						with_db_ids=(not p_dict['no_phenotype_ids']))  #load phenotype file
	else:
		print 'Retrieving the phenotypes from the DB.'
		phed = phenotypeData.getPhenotypes()

	#If on the cluster, then set up runs..
	if p_dict['parallel']:
		if len(p_dict['pids']) == 0:  #phenotype index arguement is missing, hence all phenotypes are run/analyzed.
			if not p_dict['phen_file']:
				raise Exception('Phenotype file or phenotype ID is missing.')
			p_dict['pids'] = phed.phenIds
		else:
			raise Exception('Too many arguments..')

		if analysis_plots:  #Running on the cluster..
			for p_i in p_dict['pids']:
				run_parallel(p_i, phed, p_dict)
		else:
			for mapping_method in p_dict['specific_methods']:
				for trans_method in p_dict['specific_transformations']:
					for p_i in pids:
						run_parallel(p_i, phed, p_dict, mapping_method, trans_method)
		return #Exiting the program...


	#SNPs data file name
	if not p_dict['data_file']:
		snps_data_file = '%s250K_t%d.csv' % (env['data_dir'], p_dict['call_method_id'])
	else:
		snps_data_file = p_dict['data_file']



	#Plot analysis plots...
	if p_dict['analysis_plots']:
		analysis_plots(snps_data_file, phed, p_dict)
	else:
		#If not analysis plots... then GWAS
		for p_i in p_dict['pids']:
			print '-' * 120, '\n'
			phenotype_name = phed.getPhenotypeName(p_i)
			print "Performing GWAS for phenotype: %s, phenotype_id: %s" % (phenotype_name, p_i)
			for trans_method in p_dict['specific_transformations']:
				print 'Phenotype transformation:', trans_method

				for mapping_method in p_dict['specific_methods']:
					#DO ANALYSIS
					print 'Mapping method:', mapping_method
					map_phenotype(p_i, phed, snps_data_file, mapping_method, trans_method, p_dict)


def calc_correlations(snps, phen_vals):
	import scipy as sp
	corrs = sp.zeros(len(snps))
	for i, snp in enumerate(snps):
		corrs[i] = sp.corrcoef(snp, phen_vals)[1, 0]
	return corrs



if __name__ == '__main__':
	_run_()
	print "Done!"
