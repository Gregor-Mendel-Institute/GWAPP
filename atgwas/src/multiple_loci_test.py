"""
Tests involving stepwise regressions, and model selection.

Option:

	-i ...			The run phenotype index/indices.  
	-o ...			The run_id, used as unique identifier for this run.  
	-s			Collect results and plot things (Does not generate pvalue files...) 

	-t ...			What data set is used. (75 is default)

	-n ...			Number of SNPs (phenotypes) per node, default is 1
	-d ...			Debug filter, random fraction of phenotypes or snps will be used.

	-l ...			Type of latent variable: random_snps (default), random, pc_split, etc..
	-h ...			Heritability in percentages, possible values are 1,10,25,50,75,90,99
	
	-m ...			How to generate the phenotypes: plus, xor, or

	--save_plots           	Plot Manhattan plots
	--save_pvals           	Include p-values into the .pickled files.
	--phen_file=...		File where the phenotypes will be saved, and loaded from.
	--sim_phen		Simulate phenotype, write to phenotype file.
	--parallel		Run parallel on the cluster
	--num_steps=...		Number of steps in the regression. (Default is 10)
	--herit_plots=...	Plot heritabilities
	--var_plots		Plot variation analysis


Examples:
python multiple_loci_test.py run_index test -i 1,5 -a kw,emmax -b most_normal -r ~/Projects/Data/phenotypes/phen_raw_092910.tsv 

"""
import matplotlib
matplotlib.use('Agg')

import pylab
import cPickle
import scipy as sp
import linear_models as lm
import scipy.linalg as linalg
import phenotypeData
import snpsdata
import sys
import os
import env
import random
import dataParsers as dp
import util
import gwaResults as gr
import analyze_gwas_results as agr
import traceback
import getopt
import time
import pdb


mapping_methods = ['LM', 'KW', 'EX', 'Stepw_LM', 'Stepw_EX'] #5 in total

def parse_parameters(min_num_parameters=2):
	'Parse the parameters into a dict, etc.'
	if len(sys.argv) < min_num_parameters:
		print 'Warning: there are not enough parameters..'
		print __doc__
		sys.exit(2)

	long_options_list = ['save_plots', 'phen_file=', 'sim_phen', 'num_steps=', 'parallel', 'herit_plots=', 'var_plots', 'save_pvals']
	try:
		opts, args = getopt.getopt(sys.argv[1:], "i:o:t:k:n:d:l:m:h:s", long_options_list)

	except:
		traceback.print_exc()
		print __doc__
		sys.exit(2)

	p_dict = {'number_per_run':20, 'debug_filter':1.0, 'summarize':False,
		'latent_variable':'random_snp', 'phenotype_model':'plus', 'run_id':'mlt',
		'mapping_method':'emmax', 'heritability':50, 'save_plots':False, 'call_method_id':75,
		'phen_file':env.env['phen_dir'] + 'multi_locus_phen.pickled', 'num_steps':10,
		'phen_index':None, 'sim_phen':False, 'parallel':False, 'herit_plots':None,
		'var_plots':False, 'save_pvals':False}


	for opt, arg in opts:
		if opt in ('-i'): p_dict['phen_index'] = util.parse_ids(arg)
		elif opt in ('-o'): p_dict['run_id'] = arg
		elif opt in ('-t'): p_dict['call_method_id'] = int(arg)
		elif opt in ('-n'): p_dict['number_per_run'] = int(arg)
		elif opt in ('-m'): p_dict['phenotype_model'] = arg
		elif opt in ('-d'): p_dict['debug_filter'] = float(arg)
		elif opt in ('-l'): p_dict['latent_variable'] = arg
		elif opt in ("-s"): p_dict['summarize'] = True
		elif opt in ('-h'): p_dict['heritability'] = int(arg)
		elif opt in ("--phen_file"): p_dict['phen_file'] = arg
		elif opt in ("--save_plots"): p_dict['save_plots'] = True
		elif opt in ("--save_pvals"): p_dict['save_pvals'] = True
		elif opt in ("--sim_phen"): p_dict['sim_phen'] = True
		elif opt in ("--num_steps"): p_dict['num_steps'] = int(arg)
		elif opt in ("--parallel"): p_dict['parallel'] = True
		elif opt in ("--herit_plots"): p_dict['herit_plots'] = util.parse_ids(arg)
		elif opt in ("--var_plots"): p_dict['var_plots'] = True
		else:
			print "Unkown option!!\n"
			print __doc__
			sys.exit(2)

	print p_dict, args
	return p_dict, args




def run_parallel(run_id, start_i, stop_i, latent_var, heritability, phen_model,
		phen_file, summary_run, call_method_id, num_steps, cluster='gmi'):
	"""
	If no mapping_method, then analysis run is set up.
	"""
	phen_ids_str = '%d-%d' % (start_i, stop_i)
	job_id = '%s_%s_lv%s_h%d_m%s_ns%d_t%d' % (run_id, phen_ids_str , latent_var,
						heritability, phen_model, num_steps, call_method_id)
	file_prefix = env.env['results_dir'] + job_id

	#Cluster specific parameters	
	if cluster == 'gmi': #GMI cluster.  
		shstr = '#!/bin/bash\n'
		shstr += '#$ -S /bin/bash\n'
		shstr += '#$ -N %s\n' % job_id
		shstr += '#$ -o %s_job_$JOB_ID.out\n' % file_prefix
		shstr += '#$ -e %s_job_$JOB_ID.err\n' % file_prefix
		shstr += 'source /etc/modules-env.sh\n'
		shstr += 'module load scipy/GotoBLAS2/0.9.0\n'
		shstr += 'module load matplotlib/1.0.0\n'
		shstr += 'module load mysqldb/1.2.3\n'
		shstr += 'module load h5py/2.0.0\n'
		shstr += 'export GOTO_NUM_THREADS=1\n'

	elif cluster == 'usc':  #USC cluster.
		shstr = "#!/bin/csh\n"
		shstr += "#PBS -l walltime=%s \n" % '72:00:00'
		shstr += "#PBS -l mem=%s \n" % '2950mb'
		shstr += "#PBS -q cmb\n"
		shstr += "#PBS -N p%s \n" % job_id

	shstr += "python %smultiple_loci_test.py -i %s -l %s -h %d -m %s --phen_file=%s -t %d --num_steps=%d " % \
		(env.env['script_dir'], phen_ids_str, latent_var, heritability, phen_model, phen_file, call_method_id, num_steps)
	if summary_run:
		shstr += '-s '

	#shstr += "> " + file_prefix + "_job.out) >& " + file_prefix + "_job.err\n"
	print '\n', shstr, '\n'
	script_file_name = run_id + ".sh"
	f = open(script_file_name, 'w')
	f.write(shstr)
	f.close()

	#Execute qsub script
	if cluster == 'gmi':
		os.system("qsub " + script_file_name)
	elif cluster == 'usc':
		os.system("qsub " + script_file_name)




def get_snps_heritabilities(snps, phenotype):
	Y = sp.mat(phenotype).T
	rss_0 = sp.var(Y) * len(phenotype)
	X_0 = sp.mat(sp.ones((len(phenotype), 1)))

	h_expl = []
	for snp in snps:
		rss = linalg.lstsq(sp.hstack([X_0, sp.mat(snp).T]), Y)[1]
		h_expl.append(1 - (rss / rss_0))
	return h_expl


def __get_latent_snps__(ets):
	ecotype_info_dict = phenotypeData.get_ecotype_id_info_dict()
	north_south_split_snp = []
	lats = [ecotype_info_dict[int(et)][2] for et in ets]
	m = sp.median(lats)
	for et in ets:
		latitude = ecotype_info_dict[int(et)][2]
		north_south_split_snp.append(1) if latitude > m else north_south_split_snp.append(0)
	pc_snp = []
	K = dp.load_kinship() #All ecotypes
	(evals, evecs) = linalg.eigh(K)
	pc = (sp.mat(evecs).T[-1]).tolist()[0]
	m = sp.median(pc)
	for v in pc:
		pc_snp.append(1) if v > m else pc_snp.append(0)
	return sp.array(north_south_split_snp, dtype='int8'), sp.array(pc_snp, dtype='int8')


def simulate_phenotypes(phen_file, sd, mac_threshold=0, debug_filter=1.0, num_phens=100):
	"""
	Simulate the phenotypes
	"""
	print 'Generating the phenotypes'
	latent_var_keys = ['random_snp', 'random', 'north_south_split', 'pc_split']
	phenotype_models = ['xor', 'or', 'plus', 'xor2']
	heritabilities = [1, 2, 5, 10, 15, 20, 25, 50] #in percentages

	if mac_threshold > 0:
		sd.filter_mac_snps(mac_threshold)
	num_lines = len(sd.accessions)  #Number of lines
	mafs = sd.get_mafs()["marfs"]
	if debug_filter:
		sd.sample_snps(debug_filter)
	snp_chr_pos_maf_list = sd.get_all_snp_w_info()
	all_indices = range(len(snp_chr_pos_maf_list))
	snp_indices = random.sample(all_indices, num_phens)
	map(all_indices.remove, snp_indices)

	#The first locus..
	snp_chr_pos_maf_list = [snp_chr_pos_maf_list[i] for i in snp_indices]

	#Invert every other SNP (randomize the SNP decoding)
	all_indices = range(len(snp_chr_pos_maf_list))
	invert_indices = random.sample(all_indices, num_phens / 2)
	for i in invert_indices:
		snp, chr, pos, maf = snp_chr_pos_maf_list[i]
		snp_chr_pos_maf_list[i] = (lm.get_anti_snp(snp), chr, pos, maf)

	north_south_split_snp, pc_snp = __get_latent_snps__(sd.accessions)

	phen_dict = {'snp_chr_pos_maf_list': snp_chr_pos_maf_list, 'snp_indices':snp_indices,
			'north_south_split_snp':north_south_split_snp, 'pc_snp':pc_snp}
	for latent_var in latent_var_keys:
		d = {}
		if latent_var == 'random_snp':
			l_snp_indices = random.sample(all_indices, num_phens)
			latent_snps = [snp_chr_pos_maf_list[i][0] for i in l_snp_indices]
			d['latent_chr_pos_maf_list'] = \
				[(snp_chr_pos_maf_list[i][1], snp_chr_pos_maf_list[i][2], \
				snp_chr_pos_maf_list[i][3]) for i in l_snp_indices]
			d['latent_snps'] = latent_snps

		elif latent_var == 'random':
			latent_snps = []
			for i in range(num_phens):
				num_ones = random.randint(1, num_lines - 1)
				l_snp = [0] * num_lines
				one_indices = random.sample(range(num_lines), num_ones)
				for i in one_indices:
					l_snp[i] = 1
				latent_snps.append(sp.array(l_snp, dtype='int8'))
			d['latent_snps'] = latent_snps

		elif latent_var == 'north_south_split':
			latent_snp = north_south_split_snp
			d['latent_snp'] = latent_snp

		elif latent_var == 'pc_snp':
			latent_snp = pc_snp
			d['latent_snp'] = latent_snp

		for h in heritabilities:
			her = h / 100.0
			d2 = {}
			for phen_model in phenotype_models:  #Simulate all phenotype models.
				d3 = {'phenotypes': [], 'h_estimates': [], 'h_loci_est_list': []}
				for i in range(num_phens):
					if latent_var in ['random_snp', 'random']:
						latent_snp = latent_snps[i]
					snp = snp_chr_pos_maf_list[i][0]
					if phen_model == 'xor':
						phenotype = snp ^ latent_snp
					elif phen_model == 'or':
						phenotype = snp | latent_snp
					elif phen_model == 'plus':
						phenotype = snp + latent_snp
					elif phen_model == 'xor2':
						phenotype = (snp ^ latent_snp) + 0.5 * (snp & latent_snp)
					if len(sp.unique(phenotype)) > 1:
						phen_var = sp.var(phenotype, ddof=1)
						error_vector = sp.random.normal(0, 1, size=num_lines)
						error_var = sp.var(error_vector, ddof=1)
						scalar = sp.sqrt((phen_var / error_var) * ((1 - her) / her))
						phenotype = phenotype + error_vector * scalar
						h_est = phen_var / sp.var(phenotype, ddof=1)
						h_est_snp1 = sp.corrcoef(snp, phenotype)[0, 1]
						h_est_snp2 = sp.corrcoef(latent_snp, phenotype)[0, 1]
						#print phen_model, latent_var, her, h_est, h_est_snp1 ** 2, h_est_snp2 ** 2
						d3['h_loci_est_list'].append(h_est)
						d3['h_estimates'].append((h_est_snp1 ** 2, h_est_snp2 ** 2))
					else:
						print 'encountered invalid phenotype for phen_model: %s' % phen_model
						phenotype = None
					d3['phenotypes'].append(phenotype)
				d2[phen_model] = d3
			d[h] = d2
		phen_dict[latent_var] = d


	#phenotype_models for loop ends.
	f = open(phen_file, "wb")
	print "dumping phenotypes to file:", f
	cPickle.dump(phen_dict, f, protocol=2)
	f.close()


def load_phenotypes(phen_file):
	print 'Loading phenotypes and related data'
        f = open(phen_file, "rb")
	phen_dict = cPickle.load(f)
	f.close()
        print 'Loading done..'
	return phen_dict


def summarize_runs(file_prefix, latent_var, heritability, phen_model, phen_d, index_list=None):
	"""
	Summarize runs.. duh
	"""
	pd = phen_d[latent_var][heritability][phen_model]
	if not index_list:
		index_list = range(len(pd['phenotypes']))

	num_pthres = len(pval_thresholds)
	num_winsizes = len(window_sizes)
	summary_dict = {'p_her':[]}
	analysis_methods = ['LM', 'KW', 'EX', 'Stepw_LM_Bonf', 'Stepw_LM_EBIC', 'Stepw_LM_MBIC',
				'Stepw_EX_Bonf', 'Stepw_EX_EBIC', 'Stepw_EX_MBIC']

	#Initializing stuff
	for am in analysis_methods:
		if am in ['LM', 'EX', 'KW']:
			d = {'fdrs':sp.zeros((num_pthres, num_winsizes), dtype='double'),
				'tprs':sp.zeros((num_pthres, num_winsizes), dtype='double')}
			d['ks'] = []
			d['medp'] = []
		else:
			if am in ['Stepw_EX_EBIC', 'Stepw_EX_MBIC', 'Stepw_LM_EBIC', 'Stepw_LM_MBIC']:
				d = {'fdrs':sp.zeros((num_winsizes), dtype='double'),
					'tprs':sp.zeros((num_winsizes), dtype='double')}
			else:
				d = {'fdrs':sp.zeros((num_pthres, num_winsizes), dtype='double'),
					'tprs':sp.zeros((num_pthres, num_winsizes), dtype='double')}
			d['perc_var_expl'] = []
			if 'EX' in am:
				d['rem_p_her'] = []
				d['perc_phen_var_expl'] = []
				d['perc_err_var_expl'] = []

		summary_dict[am] = d

	criteria_map = {'mbonf':('Stepw_EX_Bonf', 'Stepw_LM_Bonf'), 'ebics': ('Stepw_EX_EBIC', 'Stepw_LM_EBIC'),
			'mbics': ('Stepw_EX_MBIC', 'Stepw_LM_MBIC')}

	num_files_found = 0
	print '%s_%d_%s_%s' % (file_prefix, heritability, latent_var, phen_model)
	for i in index_list:
		pickled_file = '%s_%d_%s_%s_%dresults.pickled' % (file_prefix, heritability, latent_var, phen_model, i)
		if os.path.isfile(pickled_file):
			num_files_found += 1
			with open(pickled_file) as f:
				r = cPickle.load(f)
			p_her = r['p_her']
			summary_dict['p_her'].append(p_her)
			for am in mapping_methods:
				if am == 'Stepw_LM':
					for criteria in ['mbonf', 'ebics', 'mbics']:
						for rs in ['fdrs', 'tprs']:
							rs_array = sp.array(r[am][criteria][rs])
							rs_array[rs_array == -1.0] = 0
							summary_dict[criteria_map[criteria][1]][rs] += rs_array

						summary_dict[criteria_map[criteria][1]]['perc_var_expl'].append(r[am][criteria]['perc_var_expl'])
				elif am == 'Stepw_EX':
					for criteria in ['mbonf', 'ebics', 'mbics']:
						for rs in ['fdrs', 'tprs']:
							rs_array = sp.array(r[am][criteria][rs])
							rs_array[rs_array == -1.0] = 0
							summary_dict[criteria_map[criteria][0]][rs] += rs_array
						if 'pseudo_heritability' in r[am][criteria]:
							rem_p_her = r[am][criteria]['pseudo_heritability']
						else:
							rem_p_her = r[am][criteria]['remaining_p_her']
						perc_var_expl = r[am][criteria]['perc_var_expl']
						summary_dict[criteria_map[criteria][0]]['rem_p_her'].append(rem_p_her)
						summary_dict[criteria_map[criteria][0]]['perc_var_expl'].append(perc_var_expl)
						rss_0 = r[am]['step_info_list'][0]['rss']
						rss = r[am]['step_info_list'][r[am][criteria]['opt_i']]['rss']

						if p_her > 0.01:
							summary_dict[criteria_map[criteria][0]]['perc_phen_var_expl'].append(1.0 - ((rss * rem_p_her) / (rss_0 * p_her)))
							summary_dict[criteria_map[criteria][0]]['perc_err_var_expl'].append(1.0 - ((rss * (1 - rem_p_her)) / (rss_0 * (1 - p_her))))
				elif am in ['LM', 'KW', 'EX']:
					summary_dict[am]['fdrs'] += sp.array(r[am]['fdrs'])
					summary_dict[am]['tprs'] += sp.array(r[am]['tprs'])
					summary_dict[am]['ks'].append(r[am]['ks_stat']['D'])
					summary_dict[am]['medp'].append(r[am]['med_pval'])
	print 'Found %d results' % num_files_found
	for am in analysis_methods:
		summary_dict[am]['bonf'] = {}
		for k in ['fdrs', 'tprs']:
			summary_dict[am][k] = summary_dict[am][k] / float(num_files_found)
			if not am in ['Stepw_LM_EBIC', 'Stepw_LM_MBIC', 'Stepw_EX_EBIC', 'Stepw_EX_MBIC']:
				summary_dict[am]['bonf'][k] = summary_dict[am][k][31] * 0.56 + summary_dict[am][k][32] * 0.44
	return summary_dict



def plot_tprs_fdrs(file_prefix, summary_dict):
	"""
	Plot various things relating to run summaries
	"""
	import matplotlib.font_manager
	prop = matplotlib.font_manager.FontProperties(size=10)

	#Heritabilities..
	# - histogram of each category
	# - pseudoheritabilities vs. ks and med pval. of KW and LM

	# TPRs vs. FDRs
	am_list = ['LM', 'KW', 'EX', 'Stepw_LM_Bonf', 'Stepw_EX_Bonf']
	am_colors = ['r', 'g', 'b', 'r', 'b']
	am_ls = ['--', '--', '--', '-', '-']
	am_dot_list = ['Stepw_EX_EBIC', 'Stepw_EX_MBIC', 'Stepw_LM_EBIC', 'Stepw_LM_MBIC']
	am_dot_colors = ['#22DD66', '#22DD66', '#DD2266', '#DD2266']
	am_dot_marker = ['s', '^', 's', '^']

	for w_i, ws in enumerate(window_sizes):
		pylab.figure(figsize=(7, 6))
		pylab.axes([0.09, 0.08, 0.9, 0.85])
		for am, amc, amls in zip(am_list, am_colors, am_ls):
			xs = sp.zeros(len(pval_thresholds))
			ys = sp.zeros(len(pval_thresholds))
			for pt_i, pt in enumerate(pval_thresholds):
				ys[pt_i] = summary_dict[am]['tprs'][pt_i][w_i]
				xs[pt_i] = summary_dict[am]['fdrs'][pt_i][w_i]
			pylab.plot(xs, ys, label=am, color=amc, ls=amls, alpha=0.6, marker='.')
		for am, amc, amm in zip(am_dot_list, am_dot_colors, am_dot_marker):
			pylab.plot(summary_dict[am]['fdrs'][w_i], summary_dict[am]['tprs'][w_i], label=am, marker=amm,
				ls='', color=amc, alpha=0.6)
		png_file = '%s_w%d.png' % (file_prefix, ws)
		pylab.ylabel('Power')
		pylab.xlabel('FDR')
		pylab.legend(loc=4, prop=prop, numpoints=1, scatterpoints=1)
		x_min, x_max = pylab.xlim()
		x_range = x_max - x_min
		y_min, y_max = pylab.ylim()
		y_range = y_max - y_min
		pylab.axis([x_min - 0.025 * x_range, x_max + 0.025 * x_range,
				y_min - 0.025 * y_range, y_max + 0.025 * y_range])

		pylab.savefig(png_file)
		pylab.clf()



def plot_single_tprs_fdrs(summary_dict, ax, ws, w_legend=False, y_label='Power', x_label='FDR', y_lim=None):
	"""
	Plot various things relating to run summaries
	"""
	import matplotlib.font_manager
	prop = matplotlib.font_manager.FontProperties(size=10)

	#Heritabilities..
	# - histogram of each category
	# - pseudoheritabilities vs. ks and med pval. of KW and LM

	# TPRs vs. FDRs
	am_list = ['LM', 'EX', 'Stepw_LM_Bonf', 'Stepw_EX_Bonf']
	am_labels = ['LM', 'MM', 'SWLM', 'MLMM']
	#am_colors = ['#CC9922', '#2299CC', '#FF0022', '#0022FF']
	am_colors = ['#0022FF', '#FF0022', '#0022FF', '#FF0022']
	am_ls = ['--', '--', '-', '-']
	am_dot_list = ['Stepw_LM_EBIC', 'Stepw_EX_EBIC']
	am_dot_colors = ['#0022FF', '#FF0022']
	#am_dot_labels = ['MLML EBIC', 'SWLM EBIC']

	w_i = window_sizes.index(ws)
	for am, amc, amls, am_label in zip(am_list, am_colors, am_ls, am_labels):
		xs = sp.zeros(len(pval_thresholds))
		ys = sp.zeros(len(pval_thresholds))
		for pt_i, pt in enumerate(pval_thresholds):
			ys[pt_i] = summary_dict[am]['tprs'][pt_i][w_i]
			xs[pt_i] = summary_dict[am]['fdrs'][pt_i][w_i]
		ax.plot(xs, ys, label=am_label, color=amc, ls=amls, alpha=0.65, lw=2)
		ax.plot(summary_dict[am]['bonf']['fdrs'][w_i], summary_dict[am]['bonf']['tprs'][w_i], marker='o',
			ls='', color=amc, alpha=0.65)
	for am, amc in zip(am_dot_list, am_dot_colors):
		ax.plot(summary_dict[am]['fdrs'][w_i], summary_dict[am]['tprs'][w_i], marker='v',
			ls='', color=amc, alpha=0.65)
	ax.set_ylabel(y_label)
	ax.set_xlabel(x_label)
	if y_lim:
		ax.set_ylim(y_lim)
	x_min, x_max = ax.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax.get_ylim()
	y_range = y_max - y_min
	ax.plot(x_min - x_range, y_min - y_range, color='k', marker='v', label='EBIC', alpha=0.65)
	ax.plot(x_min - x_range, y_min - y_range, color='k', marker='o', label='Bonferroni', alpha=0.65)
	if w_legend:
		ax.legend(loc=4, prop=prop, numpoints=1, scatterpoints=1)
	ax.axis([x_min - 0.025 * x_range, x_max + 0.025 * x_range,
			y_min - 0.025 * y_range, y_max + 0.025 * y_range])



def plot_single_tprs_fdrs_2(summary_dict, ax, ws, w_legend=False, y_label='Power', x_label='FDR', y_lim=None):
	"""
	Plot various things relating to run summaries
	"""
	import matplotlib.font_manager
	prop = matplotlib.font_manager.FontProperties(size=10)

	#Heritabilities..
	# - histogram of each category
	# - pseudoheritabilities vs. ks and med pval. of KW and LM

	# TPRs vs. FDRs
	am_list = ['SLR', 'EMMAX', 'FBLR-bwd', 'MLMM-bwd']
	am_labels = ['LM', 'EMMAX', 'SWLM', 'MLMM']
	am_colors = ['#0022FF', '#FF0022', '#0022FF', '#FF0022']
	am_ls = ['--', '--', '-', '-']
#	am_colors = ['#CC9922', '#2299CC', '#FF0022', '#0022FF']
#	am_ls = ['-', '-', '-', '-']
	am_dot_list = ['FBLR-bwd', 'MLMM-bwd']
	am_dot_colors = ['#0022FF', '#FF0022']

	for am, amc, amls, am_label in zip(am_list, am_colors, am_ls, am_labels):
		ax.plot(summary_dict[ws][am]['fdr'], summary_dict[ws][am]['power'], label=am_label, color=amc,
			ls=amls, alpha=0.65, lw=2)
		ax.plot(summary_dict[ws][am]['BONF']['fdr'], summary_dict[ws][am]['BONF']['power'], marker='o',
			ls='', color=amc, alpha=0.65)
	for am, amc in zip(am_dot_list, am_dot_colors):
		ax.plot(summary_dict[ws][am]['EBIC']['fdr'], summary_dict[ws][am]['EBIC']['power'], marker='v',
			ls='', color=amc, alpha=0.65)
	ax.set_ylabel(y_label)
	ax.set_xlabel(x_label)
	if y_lim:
		ax.set_ylim(y_lim)
	x_min, x_max = ax.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax.get_ylim()
	y_range = y_max - y_min
	ax.plot(x_min - x_range, y_min - y_range, color='k', marker='v', label='EBIC', alpha=0.65)
	ax.plot(x_min - x_range, y_min - y_range, color='k', marker='o', label='Bonferroni', alpha=0.65)
	if w_legend:
		ax.legend(loc=4, prop=prop, numpoints=1, scatterpoints=1)
	ax.axis([x_min - 0.025 * x_range, x_max + 0.025 * x_range,
			y_min - 0.025 * y_range, y_max + 0.025 * y_range])






def plot_herit_hist(file_prefix, her_dict, latent_var, phen_model):
	import matplotlib.font_manager
	prop = matplotlib.font_manager.FontProperties(size=10)
	file_prefix += '_%s_%s' % (latent_var, phen_model)
	png_file_name = file_prefix + '_h%s_hist.png' % ('_'.join(map(str, her_dict.keys())))
	max_bin_count = 0
	for h in sorted(her_dict):
		bin_counts, bins, patch_list = pylab.hist(her_dict[h]['p_her'], range=(0, 0.8), bins=25, alpha=0.6)
		max_bin_count = max(max_bin_count, max(bin_counts))
		pylab.axvline((h / 100.0), color='k', alpha=0.8, ls='-.')
		pylab.axvline(sp.median(her_dict[h]['p_her']), color='#DD3311', alpha=0.8, ls='-.')
	y_range = max_bin_count - 0
	pylab.axis([-0.8 * 0.025, 0.8 * 1.025, -0.025 * max_bin_count, max_bin_count * 1.025])
	pylab.xlabel('heritability')
	pylab.ylabel('Counts')
	pylab.savefig(png_file_name)

	pylab.figure()
	png_file_name = file_prefix + '_h%s_ks_her_scatter.png' % ('_'.join(map(str, her_dict.keys())))
	for h in sorted(her_dict):
		pylab.plot(her_dict[h]['p_her'], her_dict[h]['LM']['ks'], ls='', marker='.', alpha=0.5, label='herit.=%0.2f' % (h / 100.0))
	pylab.xlabel('pseudo-heritability')
	pylab.ylabel('Kolmogorov-Smirnov statistic')
	pylab.legend(loc=2, prop=prop, numpoints=1, scatterpoints=1)
	pylab.savefig(png_file_name)

	pylab.figure()
	png_file_name = file_prefix + '_h%s_pmed_her_scatter.png' % ('_'.join(map(str, her_dict.keys())))
	for h in sorted(her_dict):
		pylab.plot(her_dict[h]['p_her'], her_dict[h]['LM']['medp'], ls='', marker='.', alpha=0.5, label='herit.=%0.2f' % (h / 100.0))
	pylab.xlabel('pseudo-heritability')
	pylab.ylabel('Median pvalue bias')
	pylab.legend(loc=2, prop=prop, numpoints=1, scatterpoints=1)
	pylab.savefig(png_file_name)


	x = []
	for h in sorted(her_dict):
		pylab.figure()
		png_file_name = file_prefix + '_h%d_hist.png' % (h)
		p_her_bias = (sp.array(her_dict[h]['p_her']) - h / 100.0)
		bin_counts, bins, patch_list = pylab.hist(p_her_bias, range=(-0.25, 0.25), bins=25, alpha=0.6)
		max_bin_count = max(bin_counts)
		pylab.axvline(0.0, color='k', alpha=0.6, ls='-.')
		y_range = max_bin_count - 0
		pylab.axis([-0.25 - 0.025 * 0.5, 0.25 + 0.5 * 0.025, -0.025 * max_bin_count, max_bin_count * 1.025])
		pylab.xlabel('pseudo-heritability bias')
		pylab.ylabel('Counts')
		x.append(p_her_bias)
		pylab.savefig(png_file_name)

	#Box-plot version
	png_file_name = file_prefix + '_h%s_boxplot.png' % ('_'.join(map(str, her_dict.keys())))
	pylab.figure()
	pylab.boxplot(x)
	pylab.ylabel('psuedo-heritability bias')
	pylab.axhline(0.0, color='k', alpha=0.6, ls='-.')
	pylab.xticks(range(1, len(x) + 1), map(str, sorted(her_dict.keys())))
	pylab.xlabel('Heritability')
	pylab.savefig(png_file_name)



def plot_var(file_prefix, d, latent_variable, heritability, phen_models):
	"""
	Plots variance explained by the model for plus,  xor, and or.
	"""

	#Plot remaining heritability
	file_prefix += '_%s_%d' % (latent_variable, heritability)
	for am in ['Stepw_EX_Bonf', 'Stepw_EX_EBIC', 'Stepw_EX_MBIC']:
		png_file_name = file_prefix + '_%s_pher_boxplot.png' % (am)
		pylab.figure()
		rem_pher_list = []
		for m in phen_models:
			rem_pher_list.append(d[m][am]['rem_p_her'])
		pylab.boxplot(rem_pher_list)
		pylab.axhline(0.0, color='k', alpha=0.6, ls='-.')
		pylab.xticks(range(1, len(phen_models) + 1), phen_models)
		pylab.ylabel('Remaining pseudo-heritability in the model')
		pylab.savefig(png_file_name)

	#Plot % explained variance.
	for am in ['Stepw_EX_Bonf', 'Stepw_EX_EBIC', 'Stepw_EX_MBIC', 'Stepw_LM_Bonf', 'Stepw_LM_EBIC', 'Stepw_LM_MBIC']:
		png_file_name = file_prefix + '_%s_var_boxplot.png' % (am)
		pylab.figure()
		perc_var_list = []
		for m in phen_models:
			perc_var_list.append(d[m][am]['perc_var_expl'])
		pylab.boxplot(perc_var_list)
		pylab.axhline(heritability / 100.0, color='k', alpha=0.6, ls='-.')
		pylab.xticks(range(1, len(phen_models) + 1), phen_models)
		pylab.ylabel('Percentage of variance explained in the model')
		pylab.savefig(png_file_name)

	#Plot % explained phenot. variance
	for am in ['Stepw_EX_Bonf', 'Stepw_EX_EBIC', 'Stepw_EX_MBIC']:
		png_file_name = file_prefix + '_%s_phen_var_boxplot.png' % (am)
		pylab.figure()
		perc_var_list = []
		for m in phen_models:
			perc_var_list.append(d[m][am]['perc_phen_var_expl'])
		pylab.boxplot(perc_var_list)
		pylab.axhline(1.0, color='k', alpha=0.6, ls='-.')
		pylab.xticks(range(1, len(phen_models) + 1), phen_models)
		pylab.ylabel('Percentage of phenotypic variance explained in the model')
		pylab.savefig(png_file_name)
	#Plot % explained error variance
	for am in ['Stepw_EX_Bonf', 'Stepw_EX_EBIC', 'Stepw_EX_MBIC']:
		png_file_name = file_prefix + '_%s_err_var_boxplot.png' % (am)
		pylab.figure()
		perc_var_list = []
		for m in phen_models:
			perc_var_list.append(d[m][am]['perc_err_var_expl'])
		pylab.boxplot(perc_var_list)
		pylab.axhline(0.0, color='k', alpha=0.6, ls='-.')
		pylab.xticks(range(1, len(phen_models) + 1), phen_models)
		pylab.ylabel('Percentage of error variance explained in the model')
		pylab.savefig(png_file_name)




#def __get_thresholds(min_thres=10, max_thres=1, num_thres=18):
def __get_thresholds(min_thres=20, max_thres=0.5, num_thres=100):
	thres_step = (min_thres - max_thres) / float(num_thres)
	pval_thresholds = []
	for i in range(num_thres):
		pval_thresholds.append(max_thres + i * thres_step)
	return pval_thresholds

pval_thresholds = __get_thresholds()

window_sizes = [0, 1000, 5000, 10000, 25000, 50000, 100000]



def _update_stats_(gwa_res, c_chr, c_pos, l_chr=None, l_pos=None, significance_threshold=None, sign_res=None):
	"""
	Update result dictionary.
	"""
	res_dict = {}
	cpl = [(c_chr, c_pos)]#Causal chr_pos_list
	if l_chr != None:
		cpl.append((l_chr, l_pos))
	caus_indices = gwa_res.get_indices(cpl)
	gwa_res._rank_scores_()

	#Calculate KS and P-med..
	pvals = gwa_res.snp_results['scores'][:]
	res_dict['ks_stat'] = agr.calc_ks_stats(pvals)
	res_dict['med_pval'] = agr.calc_median(pvals)

	#Get causal p-values, and ranks
	res_dict['causal_pvals'] = [gwa_res.snp_results['scores'][i] for i in caus_indices]
	res_dict['causal_ranks'] = [gwa_res.ranks[i] for i in caus_indices]

	#Get significant chrom_pos_pvals..
	if (not sign_res) and significance_threshold :
		sign_res = gwa_res.filter_attr('scores', significance_threshold, reversed=True, return_clone=True)
		res_dict['sign_chr_pos'] = sign_res.get_chr_pos_score_list()
		res_dict['causal_dist_matrix'] = sign_res.get_distances(cpl)
	elif sign_res:
		res_dict['sign_chr_pos'] = sign_res.get_chr_pos_score_list()
		res_dict['causal_dist_matrix'] = sign_res.get_distances(cpl)


	#Of all SNPs ranked higher than the second causative... which is farthest from a nearest causative.
	dist = gwa_res.get_farthest_w_stronger_association(cpl)
	res_dict['dist_f_w_s_a'] = -1 if dist[0] > 0 else dist[1]

	#Perform power (sensitivity, TPR), FDR, FPR calculations..
	gwa_res.neg_log_trans()
	tprs_list = []
	fdrs_list = []
	causal_tpr_fdr_list = []
	for pval_thres in pval_thresholds:
		#Filter data
		gwa_res.filter_attr('scores', pval_thres)
		pow_dict = gwa_res.get_power_analysis(cpl, window_sizes)
		causal_tpr_fdr_list.append(pow_dict['caus_founds'])
		tprs_list.append(pow_dict['tprs'])
		fdrs_list.append(pow_dict['fdrs'])

	res_dict['tprs'] = tprs_list #[p_valthreshold][window_size]
	res_dict['fdrs'] = fdrs_list #[p_valthreshold][window_size]
	res_dict['causal_tpr_fdr'] = causal_tpr_fdr_list
	return res_dict



def _update_sw_stats_(res_dict, step_info_list, opt_dict, c_chr, c_pos, l_chr=None, l_pos=None,
			significance_threshold=None, type='LM'):
	"""
	Update result dictionary for a stepwise result.
	"""
	res_dict['step_info_list'] = step_info_list
	cpl = [(c_chr, c_pos)]#Causal chr_pos_list
	if l_chr != None:
		cpl.append((l_chr, l_pos))

	for criteria in ['mbonf', 'mbics', 'ebics']:
		opt_i = opt_dict[criteria]
		d = {'opt_i':opt_i}
		si = step_info_list[opt_i]
		if criteria == 'mbonf':
			tprs_list = []
			fdrs_list = []
			causal_tpr_fdr_list = []
			t_opt_i_list = []
			num_steps = len(step_info_list) / 2
			max_cof_pvals = -sp.log10([step_info_list[i]['mbonf'] for i in range(1, 2 * num_steps)])
			for pval_thres in pval_thresholds:
				t_opt_i = 0
				for i in range(num_steps):
					if max_cof_pvals[i] >= pval_thres:
						t_opt_i = i + 1
				for j in range(1, num_steps):
					i = 2 * num_steps - j - 1
					if max_cof_pvals[i] >= pval_thres:
						if j > t_opt_i:
							t_opt_i = i + 1
				if t_opt_i == 0:
					tprs_list.append([-1 for ws in window_sizes])
					fdrs_list.append([-1 for ws in window_sizes])
				else:
					t_si = step_info_list[t_opt_i]
					cpst = map(list, zip(*t_si['cofactors']))
					sign_res = gr.Result(scores=cpst[2], chromosomes=cpst[0], positions=cpst[1])
					pow_dict = sign_res.get_power_analysis(cpl, window_sizes)
					causal_tpr_fdr_list.append(pow_dict['caus_founds'])
					tprs_list.append(pow_dict['tprs'])
					fdrs_list.append(pow_dict['fdrs'])

				t_opt_i_list.append(t_opt_i)
			d['tprs'] = tprs_list #[p_valthreshold][window_size]
			d['fdrs'] = fdrs_list #[p_valthreshold][window_size]
			d['causal_tpr_fdr'] = causal_tpr_fdr_list #[p_valthreshold][window_size]
			d['t_opt_i_list'] = t_opt_i_list


		if opt_i == 0:
			#Set default values (Are these appropriate?)
			d['tprs'] = [-1 for ws in window_sizes]
			d['fdrs'] = [-1 for ws in window_sizes]
			d['sign_chr_pos'] = []
			d['causal_dist_matrix'] = []
		else:
			cpst = map(list, zip(*si['cofactors']))
			#Create a result object..
			sign_res = gr.Result(scores=cpst[2], chromosomes=cpst[0], positions=cpst[1])
			d['sign_chr_pos'] = sign_res.get_chr_pos_score_list()
			d['causal_dist_matrix'] = sign_res.get_distances(cpl)
			pow_dict = sign_res.get_power_analysis(cpl, window_sizes)
			if criteria == 'mbonf':
				d['mbonf_tprs'] = pow_dict['tprs']
				d['mbonf_fdrs'] = pow_dict['fdrs']
				d['mbonf_causal_tpr_fdr'] = pow_dict['caus_founds']
			else:
				d['tprs'] = pow_dict['tprs']
				d['fdrs'] = pow_dict['fdrs']
				d['causal_tpr_fdr'] = pow_dict['caus_founds']
		d['kolmogorov_smirnov'] = si['kolmogorov_smirnov']
		d['pval_median'] = si['pval_median']
		d['perc_var_expl'] = 1.0 - si['rss'] / step_info_list[0]['rss']
		d['num_cofactors'] = len(si['cofactors'])
		if type == 'EX':
			d['remaining_p_her'] = si['pseudo_heritability']
			d['perc_her_var_expl'] = (si['pseudo_heritability'] / step_info_list[0]['pseudo_heritability']) * d['perc_var_expl']
			d['perc_err_var_expl'] = ((1 - si['pseudo_heritability']) / (1 - step_info_list[0]['pseudo_heritability'])) * d['perc_var_expl']
		res_dict[criteria] = d




def run_analysis(sd, K, file_prefix, latent_var, heritability, phen_model, phen_index, phen_d,
		call_method_id, num_steps=10, pickle_results=True, save_plots=False, save_pvals=False):
	"""
	Perform the GWA mapping..
	using the different methods..
	
	Linear model, 
	Kruskal-Wallis
	EMMA
	
	Stepwise Linear Model (bonf. and ext. BIC)
	Stepwise EMMA (bonf. and ext. BIC)
	"""
	file_prefix += '_%d_%s_%s_%d' % (heritability, latent_var, phen_model, phen_index)

	pd = phen_d[latent_var][heritability][phen_model]
	phen_vals = pd['phenotypes'][phen_index]
	if phen_vals == None:
		print 'Found an invalid phenotype... skipping it'
		return

	result_dict = {}
	for mm in mapping_methods:
		result_dict[mm] = {}

	print "Loading SNPS dataset (again)"
	bonferroni_threshold = 1.0 / (20.0 * sd.num_snps())

	snps_list = sd.getSnps()
	(c_snp, c_chr, c_pos, c_maf) = phen_d['snp_chr_pos_maf_list'][phen_index] #Causal SNP
	highlight_loci = [(c_chr, c_pos)]
	if latent_var == 'random_snp':
		(l_chr, l_pos, l_maf) = phen_d[latent_var]['latent_chr_pos_maf_list'][phen_index]
		highlight_loci.append((l_chr, l_pos))
	else:
		l_chr, l_pos = None, None

	print "Running Analysis"
	print 'Running KW'
	p_vals = util.kruskal_wallis(snps_list, phen_vals)['ps'].tolist()
	print len(p_vals)
	kw_res = gr.Result(snps_data=sd, scores=p_vals)
	if save_plots:
		kw_file_prefix = file_prefix + '_kw'
		kw_res.plot_manhattan(png_file=kw_file_prefix + '.png', highlight_loci=highlight_loci, neg_log_transform=True,
					plot_bonferroni=True)
		agr.plot_simple_qqplots(kw_file_prefix, [kw_res], result_labels=['Kruskal-Wallis'])


	print 'Updating stats for KW'
	result_dict['KW'] = _update_stats_(kw_res, c_chr, c_pos, l_chr, l_pos,
					significance_threshold=bonferroni_threshold)
	if save_pvals:
		result_dict['KW']['ps'] = p_vals

	#Finding the interesting plots hack...
#	if result_dict['KW']['dist_f_w_s_a'] >= 0:
#		return

	print 'Running SW LM'
	if save_plots:
		lm_file_prefix = file_prefix + '_lm'
	else:
		lm_file_prefix = None
	ret_dict = lm.lm_step_wise(phen_vals, sd, num_steps=num_steps, file_prefix=lm_file_prefix,
					highlight_loci=highlight_loci, save_pvals=save_pvals)
	lm_step_info = ret_dict['step_info_list']
	lm_pvals = ret_dict['first_lm_res']['ps'].tolist()
	lm_opt_dict = ret_dict['opt_dict']
	lm_res = gr.Result(scores=lm_pvals, snps_data=sd)
	print 'Updating stats for LM'
	result_dict['LM'] = _update_stats_(lm_res, c_chr, c_pos, l_chr, l_pos,
					significance_threshold=bonferroni_threshold)
	if save_pvals:
		result_dict['LM']['ps'] = lm_pvals
	#Finding the interesting plots hack...
#	if result_dict['LM']['dist_f_w_s_a'] >= 0:
#		return
	print 'Updating stats for SW LM'
	_update_sw_stats_(result_dict['Stepw_LM'], lm_step_info, lm_opt_dict, c_chr, c_pos, l_chr, l_pos,
					significance_threshold=bonferroni_threshold)


	print 'Running SW EX'
	if save_plots:
		emmax_file_prefix = file_prefix + '_emmax'
	else:
		emmax_file_prefix = None
	ret_dict = lm.emmax_step_wise(phen_vals, K, sd, num_steps=num_steps, file_prefix=emmax_file_prefix,
					highlight_loci=highlight_loci, save_pvals=save_pvals)
	emmax_step_info = ret_dict['step_info_list']
	emmax_pvals = ret_dict['first_emmax_res']['ps'].tolist()
	emmax_opt_dict = ret_dict['opt_dict']
	emmax_res = gr.Result(scores=emmax_pvals, snps_data=sd)
	print 'Updating stats for EX'
	result_dict['EX'] = _update_stats_(emmax_res, c_chr, c_pos, l_chr, l_pos,
					significance_threshold=bonferroni_threshold)
	if save_pvals:
		result_dict['EX']['ps'] = emmax_pvals
	print 'Updating stats for SW EX'
	_update_sw_stats_(result_dict['Stepw_EX'], emmax_step_info, emmax_opt_dict, c_chr, c_pos, l_chr, l_pos,
					significance_threshold=bonferroni_threshold, type='EX')


	#Record trait pseudo-heritability:
	result_dict['p_her'] = emmax_step_info[0]['pseudo_heritability']

	if pickle_results == True:
		pickled_results_file = file_prefix + 'results.pickled'
		print 'Pickling result dict in file: %s' % pickled_results_file
		with open(pickled_results_file, 'wb') as f:
			cPickle.dump(result_dict, f, protocol=2)

	return result_dict





def _run_():
	p_dict, args = parse_parameters()
	print args
	file_prefix = env.env['results_dir'] + p_dict['run_id']

	if p_dict['sim_phen'] and p_dict['phen_file']: 	#Simulate phenotypes
		print 'Setting up phenotype simulations'
		sd = dp.load_snps_call_method(75, debug_filter=p_dict['debug_filter'])
		simulate_phenotypes(p_dict['phen_file'], sd, debug_filter=p_dict['debug_filter'])

	elif p_dict['parallel']:
		#set up parallel runs
		if p_dict['phen_index'] == None:
			phen_d = load_phenotypes(p_dict['phen_file'])
			phenotypes = phen_d[p_dict['latent_variable']][p_dict['heritability']][p_dict['phenotype_model']]['phenotypes']
			start_i = 0
			end_i = len(phenotypes)
		else:
			start_i = p_dict['phen_index'][0]
			end_i = p_dict['phen_index'][-1]
		num_per_run = p_dict['number_per_run']
		for i in range(start_i, end_i, num_per_run):
			run_parallel(p_dict['run_id'], i, i + num_per_run - 1, p_dict['latent_variable'],
					p_dict['heritability'], p_dict['phenotype_model'],
					p_dict['phen_file'], p_dict['summarize'], p_dict['call_method_id'], p_dict['num_steps'],
					cluster='gmi')

	elif p_dict['phen_index']: #Run things..
		sd = dp.load_snps_call_method(75, debug_filter=p_dict['debug_filter'])
		K = dp.load_kinship(p_dict['call_method_id'])
		phed = load_phenotypes(p_dict['phen_file'])
		results_list = []
		for pid in p_dict['phen_index']:
			try:
				result_dict = run_analysis(sd, K, file_prefix, p_dict['latent_variable'], p_dict['heritability'],
							p_dict['phenotype_model'], pid, phed,
							p_dict['call_method_id'], num_steps=p_dict['num_steps'],
							save_plots=p_dict['save_plots'], save_pvals=p_dict['save_pvals'])
			except Exception, err_str:
				print 'Run for pid=%d failed: %s' % (pid, err_str)
				sys.stderr.write('Run for pid=%d failed: %s\n' % (pid, err_str))
			results_list.append(result_dict)
		#Save as pickled

	else:
		if p_dict['summarize']:
			file_prefix = '/srv/lab/data/mlt_results/' + p_dict['run_id']
			phed = load_phenotypes(p_dict['phen_file'])
			summary_dict = summarize_runs(file_prefix, p_dict['latent_variable'], p_dict['heritability'],
							p_dict['phenotype_model'], phed,
							index_list=p_dict['phen_index'])
			plot_file_prefix = '%s_%d_%s_%s' % (file_prefix, p_dict['heritability'], p_dict['latent_variable'],
								p_dict['phenotype_model'])
			plot_tprs_fdrs(plot_file_prefix, summary_dict)

		if p_dict['herit_plots'] != None:
			d = {}
			file_prefix = '/srv/lab/data/mlt_results/' + p_dict['run_id']
			pd = load_phenotypes(p_dict['phen_file'])
			for her in p_dict['herit_plots']:
				d[her] = summarize_runs(file_prefix, p_dict['latent_variable'], her,
							p_dict['phenotype_model'], pd, index_list=p_dict['phen_index'])
			plot_herit_hist(file_prefix, d, p_dict['latent_variable'], p_dict['phenotype_model'])


		if p_dict['var_plots']:
			d = {}
			file_prefix = '/srv/lab/data/mlt_results/' + p_dict['run_id']
			pd = load_phenotypes(p_dict['phen_file'])
			for mod in ['plus', 'or', 'xor']:
				d[mod] = summarize_runs(file_prefix, p_dict['latent_variable'], p_dict['heritability'], mod, pd, index_list=p_dict['phen_index'])
			plot_var(file_prefix, d, p_dict['latent_variable'], p_dict['heritability'], ['plus', 'or', 'xor'])



#def _run_vincent_scripts_():
#	type = sys.argv[1]
#	start = int(sys.argv[2])
#	end = int(sys.argv[3])
#	if type == 'add_emmax':
#		exec_str = 'sim2loci_add_fwdbwdemmax.sh'
#	elif type == 'add_lm':
#		exec_str = 'sim2loci_add_fwdbwdlm.sh'
#	elif type == 'or_emmax':
#		exec_str = 'sim2loci_or_fwdbwdemmax.sh'
#	elif type == 'or_lm':
#		exec_str = 'sim2loci_or_fwdbwdlm.sh'
#	elif type == 'xor_emmax':
#		exec_str = 'sim2loci_xor_fwdbwdemmax.sh'
#	elif type == 'xor_lm':
#		exec_str = 'sim2loci_xor_fwdbwdlm.sh'
#	for i in range(start, end + 1):
#		exec_st = 'qsub -q cmb -l walltime=6:00:00 -l mem=2950mb /home/cmbpanfs-01/bvilhjal/vincent/jobs/' + exec_str + ' -v VARIABLE=' + str(i)
#		print exec_st
#		os.system(exec_st)



def generate_example_figure_1():
	import gwaResults as gr
	herit = 15
	pid = 13
	i_model = 'or'
	phed = load_phenotypes(env.env['phen_dir'] + 'multi_locus_phen.pickled')
	pickled_file = '%smlt_%d_random_snp_%s_%dresults.pickled' % (env.env['tmp_dir'], herit, i_model, pid)
	if os.path.isfile(pickled_file):
		with open(pickled_file) as f:
			r = cPickle.load(f)

	result_file_prefix = '%smlt_%d_random_snp_%s_%d_' % (env.env['tmp_dir'], herit, i_model, pid)
	result_files = [result_file_prefix + fn for fn in ['lm_step0.pvals', 'emmax_step0.pvals', 'emmax_step1.pvals']]
	#Load pickle file...
	r = cPickle.load(open(result_file_prefix[:-1] + 'results.pickled'))
	results = [gr.Result(result_file=fn) for fn in result_files]
	#Setting up figure
	f = pylab.figure(figsize=(7, 6))
	ax1 = f.add_axes([0.08, 0.07 + (0.667 * 0.94), 0.9, (0.3 * 0.94) ])
	ax1.spines['top'].set_visible(False)
	ax1.xaxis.set_ticks_position('bottom')
	#ax1.xaxis.set_visible(False)
	ax2 = f.add_axes([0.08, 0.07 + (0.333 * 0.94), 0.9, (0.3 * 0.94) ])
	ax2.spines['top'].set_visible(False)
	ax2.xaxis.set_ticks_position('bottom')
	#ax2.xaxis.set_visible(False)
	ax2.set_ylabel(r'$-$log$($p-value$)$')
	ax3 = f.add_axes([0.08, 0.07 , 0.9, (0.3 * 0.94) ])
	ax3.spines['top'].set_visible(False)
	ax3.xaxis.set_ticks_position('bottom')
	ax3.xaxis.set_label_position('bottom')
	res = results[0]
	chrom_ends = res.get_chromosome_ends()
	offset = 0
	tick_positions = []
	tick_labels = []
	for c_end in chrom_ends[:-1]:
		offset += c_end
		tick_positions.append(offset)
		tick_labels.append('')
	ax3.set_xticks(tick_positions)
	ax3.set_xticklabels(tick_labels)

	scpm = phed['snp_chr_pos_maf_list'][pid]
	lcpm = phed['random_snp']['latent_chr_pos_maf_list'][pid]
	highlight_loci = [(lcpm[0], lcpm[1]), (scpm[1], scpm[2])]
	print highlight_loci
	#Fill up the figure..
	cm = {1:'#11BB00', 2:'#1199EE', 3:'#11BB00', 4:'#1199EE', 5:'#11BB00'}
	results[0].plot_manhattan2(ax=ax1, neg_log_transform=True, plot_bonferroni=True, wo_xtick_labels=True,
				chrom_colormap=cm, highlight_loci=highlight_loci, sign_color='#DD1122', plot_xaxis=True,
				percentile=98)
	x_min, x_max = ax1.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax1.get_ylim()
	y_range = y_max - y_min
	ax1.text(0.96 * x_range + x_min, 0.85 * y_range + y_min, 'a')

	results[1].plot_manhattan2(ax=ax2, neg_log_transform=True, plot_bonferroni=True, wo_xtick_labels=True,
				chrom_colormap=cm, highlight_loci=highlight_loci, sign_color='#DD1122', plot_xaxis=True,
				percentile=98)
	x_min, x_max = ax2.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax2.get_ylim()
	y_range = y_max - y_min
	ax2.text(0.96 * x_range + x_min, 0.85 * y_range + y_min, 'b')

	cofactors = r['Stepw_EX']['step_info_list'][1]['cofactors']
	print r['Stepw_EX']['step_info_list'][1]['cofactors']
	results[2].plot_manhattan2(ax=ax3, neg_log_transform=True, plot_bonferroni=True,
				chrom_colormap=cm, highlight_markers=cofactors, highlight_loci=highlight_loci,
				sign_color='#DD1122', plot_xaxis=True,
				percentile=98)
	x_min, x_max = ax3.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax3.get_ylim()
	y_range = y_max - y_min
	ax3.text(0.96 * x_range + x_min, 0.85 * y_range + y_min, 'c')

#	f.text(0.193, 0.04, '1')
#	f.text(0.38, 0.04, '2')
#	f.text(0.542, 0.04, '3')
#	f.text(0.705, 0.04, '4')
#	f.text(0.875, 0.04, '5')
#	f.text(0.43, 0.01, 'Chromosome number')

	#Save the figure?
	pylab.savefig(env.env['tmp_dir'] + 'mlt_figure_1.png', dpi=400)
	pylab.savefig(env.env['tmp_dir'] + 'mlt_figure_1.pdf', format='pdf')
	#pylab.savefig(env.env['tmp_dir'] + 'test.pdf', format='pdf')




def generate_results_figure_2(png_file_name='/tmp/figure2.png', pdf_file_name='/tmp/figure2.pdf', herit=10, window_size=25000):
	#file_prefix = '/srv/lab/data/mlt_results/mlt'
	file_prefix = '/home/GMI/bjarni.vilhjalmsson/results/mlt'
	phed = load_phenotypes(env.env['phen_dir'] + 'multi_locus_phen.pickled')
	f = pylab.figure(figsize=(11, 7))

	summary_dict = summarize_runs(file_prefix, 'random_snp', herit, 'plus', phed, index_list=range(1000))
	ax = f.add_axes([0.06, 0.08 + (0.5 * 0.9), (0.3 * 0.93), (0.46 * 0.9) ])
	ax.set_xticklabels([''] * len(ax.get_xticks()))
	plot_single_tprs_fdrs(summary_dict, ax, window_size, y_lim=(0.3, 1.03), w_legend=True, x_label='')

	summary_dict = summarize_runs(file_prefix, 'random_snp', herit, 'or', phed, index_list=range(1000))
	ax = f.add_axes([0.06 + (0.33 * 0.93), 0.08 + (0.5 * 0.9), (0.3 * 0.93), (0.46 * 0.9) ])
	ax.set_yticklabels([''] * len(ax.get_yticks()))
	ax.set_xticklabels([''] * len(ax.get_xticks()))
	plot_single_tprs_fdrs(summary_dict, ax, window_size, y_lim=(0.3, 1.03), y_label='', x_label='')

	summary_dict = summarize_runs(file_prefix, 'random_snp', herit, 'xor', phed, index_list=range(1000))
	ax = f.add_axes([0.054 + (0.667 * 0.93), 0.08 + (0.5 * 0.9), (0.3 * 0.93), (0.46 * 0.9) ])
	ax.set_yticklabels([''] * len(ax.get_yticks()))
	ax.set_xticklabels([''] * len(ax.get_xticks()))
	plot_single_tprs_fdrs(summary_dict, ax, window_size, y_lim=(0.3, 1.03), y_label='', x_label='')

	summary_dict = summarize_runs(file_prefix, 'north_south_split', herit, 'plus', phed, index_list=range(1000))
	ax = f.add_axes([0.06, 0.08, (0.3 * 0.93), (0.46 * 0.9) ])
	plot_single_tprs_fdrs(summary_dict, ax, window_size, y_lim=(0.0, 1.03))

	summary_dict = summarize_runs(file_prefix, 'north_south_split', herit, 'or', phed, index_list=range(1000))
	ax = f.add_axes([0.06 + (0.33 * 0.93), 0.08, (0.3 * 0.93), (0.46 * 0.9) ])
	ax.set_yticklabels([''] * len(ax.get_yticks()))
	plot_single_tprs_fdrs(summary_dict, ax, window_size, y_lim=(0.0, 1.03), y_label='')

	summary_dict = summarize_runs(file_prefix, 'north_south_split', herit, 'xor', phed, index_list=range(1000))
	ax = f.add_axes([0.054 + (0.667 * 0.93), 0.08, (0.3 * 0.93), (0.46 * 0.9) ])
	ax.set_yticklabels([''] * len(ax.get_yticks()))
	plot_single_tprs_fdrs(summary_dict, ax, window_size, y_lim=(0.0, 1.03), y_label='')

	f.text(0.168, 0.962, 'Additive')
	f.text(0.495, 0.962, "'or'")
	f.text(0.796, 0.962, "'xor'")

	f.text(0.072, 0.91, 'a')
	f.text(0.382, 0.91, "b")
	f.text(0.688, 0.91, "c")
	f.text(0.072, 0.46, 'd')
	f.text(0.382, 0.46, "e")
	f.text(0.688, 0.46, "f")

	f.text(0.97, 0.435, 'North - south latent variable', rotation=90)
	f.text(0.97, 0.822, 'Two random SNPs', rotation=90)
	f.savefig(png_file_name)
	f.savefig(pdf_file_name, format='pdf')





def generate_results_figure_3(png_file_name=' / tmp / test2.png', pdf_file_name=' / tmp / test2.pdf'):
	file_prefix = ' / home / GMI / bjarni.vilhjalmsson / Projects / vincent_plots / plots_sim100loci / '
	caus_in_file = file_prefix + 'ROC_CAUS_IN.csv'
	caus_dropped_file = file_prefix + 'ROC_CAUS_DROP.csv'

	hertis = [25, 50, 75]
	analysis_methods = ['MLMM - fwd', 'MLMM - bwd', 'FBLR - fwd', 'FBLR - bwd', 'SLR', 'EMMAX']
	window_sizes = [0, 5000, 10000, 25000, 50000, 100000]

	caus_in_dict = {25:{}, 50:{}, 75:{}}
	caus_drop_dict = {25:{}, 50:{}, 75:{}}
	for h in hertis:
		for ws in window_sizes:
			caus_in_dict[h][ws] = {}
			caus_drop_dict[h][ws] = {}
			for am in analysis_methods:
				caus_in_dict[h][ws][am] = {'fdr':[], 'power':[]}
				caus_drop_dict[h][ws][am] = {'fdr':[], 'power':[]}

	with open(caus_in_file) as f:
		print f.next()
		for line in f:
			l = line.split(', ')
			h = int(float(l[0]) * 100)
			am = l[1][1:-1]
			ws = int(l[2]) * 1000
			thres = l[3][1:-1]
			if thres == 'ebic':
				if am in ['MLMM - fwd', 'MLMM - bwd', 'FBLR - fwd', 'FBLR - bwd']:
					fdr = float(l[4])
					power = float(l[5])
					caus_in_dict[h][ws][am]['EBIC'] = {'fdr':fdr, 'power':power}
       				continue
			fdr = float(l[4])
                        power = float(l[5])
			if thres == 'bonf':
				caus_in_dict[h][ws][am]['BONF'] = {'fdr':fdr, 'power':power}
			else:
				caus_in_dict[h][ws][am]['fdr'].append(fdr)
				caus_in_dict[h][ws][am]['power'].append(power)

	with open(caus_dropped_file) as f:
		print f.next()
		for line in f:
			l = line.split(', ')
			h = int(float(l[0]) * 100)
			am = l[1][1:-1]
			ws = int(l[2]) * 1000
			thres = l[3][1:-1]
			if thres == 'ebic':
				if am in ['MLMM - fwd', 'MLMM - bwd', 'FBLR - fwd', 'FBLR - bwd']:
					fdr = float(l[4])
					power = float(l[5])
					caus_drop_dict[h][ws][am]['EBIC'] = {'fdr':fdr, 'power':power}
				continue
			fdr = float(l[4])
			power = float(l[5])
			if thres == 'bonf':
				caus_drop_dict[h][ws][am]['BONF'] = {'fdr':fdr, 'power':power}
			else:
				caus_drop_dict[h][ws][am]['fdr'].append(fdr)
				caus_drop_dict[h][ws][am]['power'].append(power)


	f = pylab.figure(figsize=(11, 7))

	ax = f.add_axes([0.06, 0.08 + (0.5 * 0.9), (0.3 * 0.93), (0.46 * 0.9) ])
	ax.set_xticklabels([''] * len(ax.get_xticks()))
	plot_single_tprs_fdrs_2(caus_in_dict[25], ax, 25000, y_lim=(0.16, 1.03), w_legend=True, x_label='')

	ax = f.add_axes([0.06 + (0.33 * 0.93), 0.08 + (0.5 * 0.9), (0.3 * 0.93), (0.46 * 0.9) ])
	ax.set_yticklabels([''] * len(ax.get_yticks()))
	ax.set_xticklabels([''] * len(ax.get_xticks()))
	plot_single_tprs_fdrs_2(caus_in_dict[50], ax, 25000, y_lim=(0.16, 1.03), y_label='', x_label='')

	ax = f.add_axes([0.054 + (0.667 * 0.93), 0.08 + (0.5 * 0.9), (0.3 * 0.93), (0.46 * 0.9) ])
	ax.set_yticklabels([''] * len(ax.get_yticks()))
	ax.set_xticklabels([''] * len(ax.get_xticks()))
	plot_single_tprs_fdrs_2(caus_in_dict[75], ax, 25000, y_lim=(0.16, 1.03), y_label='', x_label='')

	ax = f.add_axes([0.06, 0.08, (0.3 * 0.93), (0.46 * 0.9) ])
	plot_single_tprs_fdrs_2(caus_drop_dict[25], ax, 25000, y_lim=(0.06, 1.03))

	ax = f.add_axes([0.06 + (0.33 * 0.93), 0.08, (0.3 * 0.93), (0.46 * 0.9) ])
	ax.set_yticklabels([''] * len(ax.get_yticks()))
	plot_single_tprs_fdrs_2(caus_drop_dict[50], ax, 25000, y_lim=(0.06, 1.03), y_label='')

	ax = f.add_axes([0.054 + (0.667 * 0.93), 0.08, (0.3 * 0.93), (0.46 * 0.9) ])
	ax.set_yticklabels([''] * len(ax.get_yticks()))
	plot_single_tprs_fdrs_2(caus_drop_dict[75], ax, 25000, y_lim=(0.06, 1.03), y_label='')

	f.text(0.168, 0.962, r'$h ^ 2 = 0.25 $')
	f.text(0.485, 0.962, r'$h ^ 2 = 0.5 $')
	f.text(0.785, 0.962, r'$h ^ 2 = 0.75 $')


	f.text(0.072, 0.91, 'a')
	f.text(0.382, 0.91, "b")
	f.text(0.688, 0.91, "c")
	f.text(0.072, 0.46, 'd')
	f.text(0.382, 0.46, "e")
	f.text(0.688, 0.46, "f")

	f.text(0.97, 0.105, 'Causatives dropped from data', rotation=90)
	f.text(0.97, 0.627, 'Causatives in data', rotation=90)
	f.savefig(png_file_name)
	f.savefig(pdf_file_name, format='pdf')



def _draw_var_plot_():
	pass

def generate_example_figure4():
	"""
	HUMAN FIGURE
	"""
	import math
	#chrom_d, var_d = _parse_map_file_()
	results_dir = '/Users/bjarni.vilhjalmsson/Projects/Data/NFBC_results/'
	result_files = [results_dir + 'NFBC_emmax_step_ibs_ldlres_pid4_step0.pvals',
			results_dir + 'NFBC_emmax_step_ibs_ldlres_pid4_step5.pvals']
	results = []
	for result_file in result_files:
		results.append(gr.Result(result_file))

	with open(results_dir + 'NFBC_emmax_step_ibs_ldlres_pid4_stats.csv') as f:
		lines = f.readlines()
		l5 = map(str.strip, lines[6].split(','))
		c_list = l5[14:]
		cofactors = []
		for c_p_p in c_list:
			l = c_p_p.split('_')
			cofactors.append((int(l[0]), int(l[1]), float(l[2])))
		rss_list = []
		p_her_list = []
		for line in lines[1:]:
			l = map(str.strip, line.split(','))
			rss_list.append(float(l[1]))
			p_her_list.append(float(l[9]))
		num_steps = len(lines) - 1
		rss_list.append(rss_list[0])
		p_her_list.append(p_her_list[0])

	print cofactors

	chrom_col_map = {}
	for i in range(1, 24):
		if i % 2 == 0:
			chrom_col_map[i] = '#1199EE'
		else:
			chrom_col_map[i] = '#11BB00'

	#Draw plot
	f = pylab.figure(figsize=(7, 7))
	ax1 = f.add_axes([0.08, 0.07 + (0.667 * 0.94), 0.9, (0.3 * 0.94) ])
	ax1.spines['top'].set_visible(False)
	ax1.xaxis.set_ticks_position('bottom')
	ax1.xaxis.set_label_position('bottom')
	ax1.set_ylabel(r'$-$log$($p-value$)$')
	ax2 = f.add_axes([0.08, 0.07 + (0.333 * 0.94), 0.9, (0.3 * 0.94) ])
	ax2.spines['top'].set_visible(False)
	ax2.set_ylabel(r'$-$log$($p-value$)$')
	ax2.xaxis.set_ticks_position('bottom')
	ax2.xaxis.set_label_position('bottom')
	ax3 = f.add_axes([0.08, 0.06 , 0.69, (0.24 * 0.94) ]) #Partition of variance axes
	#ax3.spines['top'].set_visible(False)
	ax3.xaxis.set_ticks_position('bottom')
	ax3.xaxis.set_label_position('bottom')
	ax3.set_ylabel('Partition of variance')
	ax3.set_xlabel('Step number')


	#Fill up the figure..
	results[0].plot_manhattan2(ax=ax1, neg_log_transform=True, plot_bonferroni=True, plot_xaxis=True,
				chrom_colormap=chrom_col_map, sign_color='#DD1122', wo_xtick_labels=True,
				percentile=98)
	x_min, x_max = ax1.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax1.get_ylim()
	y_range = y_max - y_min
	ax1.text(0.96 * x_range + x_min, 0.85 * y_range + y_min, 'a')

	results[1].plot_manhattan2(ax=ax2, neg_log_transform=True, plot_bonferroni=True, plot_xaxis=True,
				chrom_colormap=chrom_col_map, highlight_markers=cofactors, sign_color='#DD1122',
				percentile=98)
	x_min, x_max = ax2.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax2.get_ylim()
	y_range = y_max - y_min
	ax2.text(0.96 * x_range + x_min, 0.85 * y_range + y_min, 'b')
	ax2.text(0.061 * x_range + x_min, 0.86 * y_range + y_min, '1', fontsize='small')
	ax2.text(0.112 * x_range + x_min, 0.63 * y_range + y_min, '2', fontsize='small')
	ax2.text(0.845 * x_range + x_min, 0.506 * y_range + y_min, '3', fontsize='small')
	ax2.text(0.897 * x_range + x_min, 0.678 * y_range + y_min, '4', fontsize='small')
	ax2.text(0.856 * x_range + x_min, 0.931 * y_range + y_min, '5', fontsize='small')



	max_rss = max(rss_list)
	rss_array = sp.array(rss_list) / max_rss
	p_her_array = rss_array * sp.array(p_her_list)
	genetic_variance = p_her_array + (1 - rss_array)
	variance_explained = (1 - rss_array)
	ax3.fill_between([0, num_steps], 0, 1.1, color='#FFCC33', alpha=0.8, label='Error')
	ax3.fill_between(sp.arange(num_steps + 1), 0, genetic_variance, color='#22CC44', alpha=0.8, label='Genetic')
	ax3.fill_between(sp.arange(num_steps + 1), -0.1, variance_explained, color='#2255AA', alpha=0.8, label='Explained')
	ax3.axvline(x=num_steps / 2, c='k', linestyle=':')
	ax3.axis([0, num_steps, 0, 1])
	x_min, x_max = ax3.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax3.get_ylim()
	y_range = y_max - y_min
	ax3.text(1.25 * x_range + x_min, 0.95 * y_range + y_min, 'c')
	ax3.bar(-1, 1, color='#FFCC33', alpha=0.8, label='Error')
	ax3.bar(-1, 1, color='#22CC44', alpha=0.8, label='Genetic')
	ax3.bar(-1, 1, color='#2255AA', alpha=0.8, label='Explained')
	ax3.axvline(5, color='#AA0000', alpha=0.8, lw=1.8)
	ax3.text(0.17 * x_range + x_min, 1.04 * y_range + y_min, 'EBIC/MBONF', fontsize='small')
	leg = ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))

	ltext = leg.get_texts()  # all the text.Text instance in the legend
	llines = leg.get_lines()  # all the lines.Line2D instance in the legend
	frame = leg.get_frame()  # the patch.Rectangle instance surrounding the legend

	# see text.Text, lines.Line2D, and patches.Rectangle for more info on
	# the settable properties of lines, text, and rectangles
	frame.set_linewidth(0)
	#frame.set_facecolor('0.80')      # set the frame face color to light gray
	pylab.setp(ltext, fontsize='small')    # the legend text fontsize
	pylab.setp(llines, linewidth=1.5)      # the legend linewidth


	#Save the figure?
	pylab.savefig(env.env['tmp_dir'] + 'mlt_NFBC_example.png', dpi=400)
	pylab.savefig(env.env['tmp_dir'] + 'mlt_NFBC_example.pdf', format='pdf')
	#pylab.savefig(env.env['tmp_dir'] + 'test.pdf', format='pdf')



def generate_example_figure5():
	"""
	SODIUM FIGURE
	"""

	results_dir = env.env['tmp_dir']
	result_files = [results_dir + 'abracadabra_pid226_Na23_Soil_3_emmax_step_log_t75_step0.pvals',
			results_dir + 'abracadabra_pid226_Na23_Soil_3_emmax_step_log_t75_step3.pvals']
	results = []
	for result_file in result_files:
		results.append(gr.Result(result_file))

	with open(results_dir + 'abracadabra_pid226_Na23_Soil_3_emmax_step_log_t75_stats.csv') as f:
		lines = f.readlines()
		l5 = map(str.strip, lines[4].split(','))
		c_list = l5[14:]
		cofactors = []
		for c_p_p in c_list:
			l = c_p_p.split('_')
			cofactors.append((int(l[0]), int(l[1]), float(l[2])))
		rss_list = []
		p_her_list = []
		for line in lines[1:]:
			l = map(str.strip, line.split(','))
			rss_list.append(float(l[1]))
			p_her_list.append(float(l[9]))
		num_steps = len(lines) - 1
		rss_list.append(rss_list[0])
		p_her_list.append(p_her_list[0])

	print cofactors

	chrom_col_map = {}
	for i in range(1, 6):
		if i % 2 == 0:
			chrom_col_map[i] = '#1199EE'
		else:
			chrom_col_map[i] = '#11BB00'

	#Draw plot
	f = pylab.figure(figsize=(7, 7))
	ax1 = f.add_axes([0.08, 0.07 + (0.667 * 0.94), 0.9, (0.3 * 0.94) ])
	ax1.spines['top'].set_visible(False)
	ax1.xaxis.set_ticks_position('bottom')
	ax1.xaxis.set_label_position('bottom')
	ax1.set_ylabel(r'$-$log$($p-value$)$')
	ax2 = f.add_axes([0.08, 0.07 + (0.333 * 0.94), 0.9, (0.3 * 0.94) ])
	ax2.spines['top'].set_visible(False)
	ax2.set_ylabel(r'$-$log$($p-value$)$')
	ax2.xaxis.set_ticks_position('bottom')
	ax2.xaxis.set_label_position('bottom')
	ax3 = f.add_axes([0.08, 0.06 , 0.69, (0.24 * 0.94) ]) #Partition of variance axes
	#ax3.spines['top'].set_visible(False)
	ax3.xaxis.set_ticks_position('bottom')
	ax3.xaxis.set_label_position('bottom')
	ax3.set_ylabel('Partition of variance')
	ax3.set_xlabel('Step number')


	#Fill up the figure..
	results[0].plot_manhattan2(ax=ax1, neg_log_transform=True, plot_bonferroni=True, plot_xaxis=True,
				chrom_colormap=chrom_col_map, sign_color='#DD1122', wo_xtick_labels=True,
				percentile=98)
	x_min, x_max = ax1.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax1.get_ylim()
	y_range = y_max - y_min
	ax1.text(0.96 * x_range + x_min, 0.85 * y_range + y_min, 'a')

	results[1].plot_manhattan2(ax=ax2, neg_log_transform=True, plot_bonferroni=True, plot_xaxis=True,
				chrom_colormap=chrom_col_map, highlight_markers=cofactors, sign_color='#DD1122',
				percentile=98)
	x_min, x_max = ax2.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax2.get_ylim()
	y_range = y_max - y_min
	ax2.text(0.96 * x_range + x_min, 0.85 * y_range + y_min, 'b')
	ax2.text(0.68 * x_range + x_min, 0.93 * y_range + y_min, '1', fontsize='small')
	ax2.text(0.642 * x_range + x_min, 0.619 * y_range + y_min, '2', fontsize='small')
	ax2.text(0.682 * x_range + x_min, 0.483 * y_range + y_min, '3', fontsize='small')


	max_rss = max(rss_list)
	rss_array = sp.array(rss_list) / max_rss
	p_her_array = rss_array * sp.array(p_her_list)
	genetic_variance = p_her_array + (1 - rss_array)
	variance_explained = (1 - rss_array)
	ax3.fill_between([0, num_steps], 0, 1.1, color='#FFCC33', alpha=0.8, label='Error')
	ax3.fill_between(sp.arange(num_steps + 1), 0, genetic_variance, color='#22CC44', alpha=0.8, label='Genetic')
	ax3.fill_between(sp.arange(num_steps + 1), -0.1, variance_explained, color='#2255AA', alpha=0.8, label='Explained')
	ax3.axvline(3, color='#AA0000', alpha=0.8, lw=1.8)
	ax3.axvline(x=num_steps / 2, c='k', linestyle=':')
	ax3.axis([0, num_steps, 0, 1])
	x_min, x_max = ax3.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax3.get_ylim()
	y_range = y_max - y_min
	ax3.text(0.1 * x_range + x_min, 1.04 * y_range + y_min, 'EBIC/MBONF', fontsize='small')
	ax3.text(1.25 * x_range + x_min, 0.95 * y_range + y_min, 'c')
	ax3.bar(-1, 1, color='#FFCC33', alpha=0.8, label='Error')
	ax3.bar(-1, 1, color='#22CC44', alpha=0.8, label='Genetic')
	ax3.bar(-1, 1, color='#2255AA', alpha=0.8, label='Explained')
	leg = ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))

	ltext = leg.get_texts()  # all the text.Text instance in the legend
	llines = leg.get_lines()  # all the lines.Line2D instance in the legend
	frame = leg.get_frame()  # the patch.Rectangle instance surrounding the legend

	# see text.Text, lines.Line2D, and patches.Rectangle for more info on
	# the settable properties of lines, text, and rectangles
	frame.set_linewidth(0)
	#frame.set_facecolor('0.80')      # set the frame face color to light gray
	pylab.setp(ltext, fontsize='small')    # the legend text fontsize
	pylab.setp(llines, linewidth=1.5)      # the legend linewidth


	#Save the figure?
	pylab.savefig(env.env['tmp_dir'] + 'mlt_sodium_example.png', dpi=400)
	pylab.savefig(env.env['tmp_dir'] + 'mlt_sodium_example.pdf', format='pdf')
	#pylab.savefig(env.env['tmp_dir'] + 'test.pdf', format='pdf')


def generate_example_bayesian_figures():
	import gwaResults as gr
	ex_res = gr.Result(result_file=env.env['tmp_dir'] + 'bf_qtl_priors_DTF2ndSwAverage2009_317_step0.pvals')
	ex_ppa_res = gr.Result(result_file=env.env['tmp_dir'] + 'bf_qtl_priors_DTF2ndSwAverage2009_317_opt_mbonf_step15.ppas')
	c_ends = ex_res.get_chromosome_ends()

	cg_tair_ids = ['AT1G04400', 'AT1G65480', 'AT1G77080', 'AT2G26330', 'AT4G00650', 'AT5G10140', 'AT5G65050', \
			'AT5G65060', 'AT5G65070', 'AT5G65080'] \
			#FT, MAF, ER, FRI, FLC, MAF2-MAF5 (Salome et al. 2011)
	cgs = gr.get_genes_w_tair_id(cg_tair_ids)

	chromosomes = []
	positions = []
	scores = []
	last_chr = 1
	with open('/Users/bjarni.vilhjalmsson/Projects/Data/DTF1.scan.tsv') as f:
		print f.next()
		for l in f:
			line = l.split()
			chrom = int(line[1])
			if last_chr != chrom:
				chromosomes.append(last_chr)
				positions.append(c_ends[last_chr - 1])
				scores.append(scores[-1])
				last_chr = chrom
			chromosomes.append(chrom)
			positions.append(int(line[2]))
			scores.append(float(line[4]))


	linkage_res = gr.Result(positions=positions, chromosomes=chromosomes, scores=scores)


	with open(env.env['tmp_dir'] + 'bf_qtl_priors_DTF2ndSwAverage2009_317_stats.csv') as f:
		lines = f.readlines()
		l15 = map(str.strip, lines[16].split(','))
		c_list = l15[14:]
		cofactors = []
		for c_p_p in c_list:
			l = c_p_p.split('_')
			cofactors.append((int(l[0]), int(l[1]), float(l[2])))
		rss_list = []
		p_her_list = []
		for line in lines[1:]:
			l = map(str.strip, line.split(','))
			rss_list.append(float(l[1]))
			p_her_list.append(float(l[9]))
		num_steps = len(lines) - 1
		rss_list.append(rss_list[0])
		p_her_list.append(p_her_list[0])

	print cofactors
	ppa_cofactors = [[5, 3188327, 0.99389222891327844], [4, 177047, 0.99613567887355359], [4, 1330863, 0.95803471534031548], [4, 341989, 0.98581862771623219], [5, 21512983, 0.96576709966410257]]

	chrom_col_map = {}
	for i in range(1, 6):
		if i % 2 == 0:
			chrom_col_map[i] = '#1199EE'
		else:
			chrom_col_map[i] = '#11BB00'

	#Draw plot
	f = pylab.figure(figsize=(9, 9))
	ax1 = f.add_axes([0.08, 0.07 + (0.84 * 0.94), 0.9, (0.13 * 0.94) ])
	ax1.spines['top'].set_visible(False)
	ax1.xaxis.set_ticks_position('bottom')
	ax1.xaxis.set_label_position('bottom')
	ax1.set_ylabel(r'$-$log$($p-value$)$')
	ax2 = f.add_axes([0.08, 0.07 + (0.54 * 0.94), 0.9, (0.25 * 0.94) ])
	ax2.spines['top'].set_visible(False)
	ax2.set_ylabel(r'$-$log$($p-value$)$')
	ax2.xaxis.set_ticks_position('bottom')
	ax2.xaxis.set_label_position('bottom')
	ax3 = f.add_axes([0.08, 0.07 + (0.25 * 0.94), 0.9, (0.25 * 0.94) ])
	ax3.spines['top'].set_visible(False)
	ax3.set_ylabel('Posterior prob. of assoc.')
	ax3.xaxis.set_ticks_position('bottom')
	ax3.xaxis.set_label_position('bottom')
	ax4 = f.add_axes([0.08, 0.06 , 0.73, (0.18 * 0.94) ]) #Partition of variance axes
	#ax4.spines['top'].set_visible(False)
	ax4.xaxis.set_ticks_position('bottom')
	ax4.xaxis.set_label_position('bottom')
	ax4.set_ylabel('Partition of variance')
	ax4.set_xlabel('Step number')


	linkage_res.plot_manhattan2(ax=ax1, neg_log_transform=False, plot_bonferroni=False, plot_xaxis=True,
				chrom_colormap=chrom_col_map, wo_xtick_labels=True, percentile=0.0)
	x_min, x_max = ax1.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax1.get_ylim()
	y_range = y_max - y_min
	ax1.text(0.96 * x_range + x_min, 0.93 * y_range + y_min, 'a')

	ex_res.plot_manhattan2(ax=ax2, neg_log_transform=True, plot_bonferroni=True, plot_xaxis=True, wo_xtick_labels=True,
				chrom_colormap=chrom_col_map, sign_color='#DD1122', cand_genes=cgs)
	x_min, x_max = ax2.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax2.get_ylim()
	y_range = y_max - y_min
	ax2.text(0.96 * x_range + x_min, 0.94 * y_range + y_min, 'b')

	ex_ppa_res.plot_manhattan2(ax=ax3, neg_log_transform=False, plot_bonferroni=False, plot_xaxis=True,
				chrom_colormap=chrom_col_map, highlight_markers=ppa_cofactors, sign_color=None,
				cand_genes=cgs)
	x_min, x_max = ax3.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax3.get_ylim()
	y_range = y_max - y_min
	ax3.text(0.96 * x_range + x_min, 0.94 * y_range + y_min, 'c')

	max_rss = max(rss_list)
	rss_array = sp.array(rss_list) / max_rss
	p_her_array = rss_array * sp.array(p_her_list)
	genetic_variance = p_her_array + (1 - rss_array)
	variance_explained = (1 - rss_array)
	ax4.fill_between([0, num_steps], 0, 1.1, color='#FFCC33', alpha=0.8, label='Error')
	ax4.fill_between(sp.arange(num_steps + 1), 0, genetic_variance, color='#22CC44', alpha=0.8, label='Genetic')
	ax4.fill_between(sp.arange(num_steps + 1), -0.1, variance_explained, color='#2255AA', alpha=0.8, label='Explained')
	ax4.axvline(x=num_steps / 2, c='k', linestyle=':')
	ax4.axis([0, num_steps, 0, 1])
	x_min, x_max = ax4.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax4.get_ylim()
	y_range = y_max - y_min
	ax4.text(1.185 * x_range + x_min, 0.94 * y_range + y_min, 'd')
	ax4.bar(-1, 1, color='#FFCC33', alpha=0.8, label='Error')
	ax4.bar(-1, 1, color='#22CC44', alpha=0.8, label='Genetic')
	ax4.bar(-1, 1, color='#2255AA', alpha=0.8, label='Explained')
	ax4.axvline(15, color='#AA0000', alpha=0.8, lw=1.8)
	ax4.text(0.712 * x_range + x_min, 1.04 * y_range + y_min, 'MBONF', fontsize='small')
	ax4.axvline(2, color='#AA0000', alpha=0.8, lw=1.8)
	ax4.text(0.079 * x_range + x_min, 1.04 * y_range + y_min, 'EBIC', fontsize='small')
	leg = ax4.legend(loc='center left', bbox_to_anchor=(1, 0.5))

	ltext = leg.get_texts()  # all the text.Text instance in the legend
	llines = leg.get_lines()  # all the lines.Line2D instance in the legend
	frame = leg.get_frame()  # the patch.Rectangle instance surrounding the legend

	# see text.Text, lines.Line2D, and patches.Rectangle for more info on
	# the settable properties of lines, text, and rectangles
	frame.set_facecolor('0.80')      # set the frame face color to light gray
	pylab.setp(ltext, fontsize='small')    # the legend text fontsize
	pylab.setp(llines, linewidth=1.5)      # the legend linewidth



	pylab.savefig(env.env['tmp_dir'] + 'mlt_bayesian_example.png', dpi=400)




def generate_example_bayesian_figures_flc():
	import gwaResults as gr
	ex_res = gr.Result(result_file=env.env['tmp_dir'] + 'bf_flc_FLC_43_step0.pvals')
	ex_ppa_res = gr.Result(result_file=env.env['tmp_dir'] + 'bf_flc_FLC_43_step2.ppas')
	c_ends = ex_res.get_chromosome_ends()

	cg_tair_ids = ['AT4G00650']#, 'AT5G10140'] \
			#FRI, FLC
	cgs = gr.get_genes_w_tair_id(cg_tair_ids)

	with open(env.env['tmp_dir'] + 'bf_flc_FLC_43_stats.csv') as f:
		lines = f.readlines()
		l15 = map(str.strip, lines[16].split(','))
		c_list = l15[14:]
		cofactors = []
		for c_p_p in c_list:
			l = c_p_p.split('_')
			cofactors.append((int(l[0]), int(l[1]), float(l[2])))
		rss_list = []
		p_her_list = []
		for line in lines[1:]:
			l = map(str.strip, line.split(','))
			rss_list.append(float(l[1]))
			p_her_list.append(float(l[9]))
		num_steps = len(lines) - 1
		rss_list.append(rss_list[0])
		p_her_list.append(p_her_list[0])

	print cofactors
	ppa_cofactors = [[4, 268809, 0.81607629105215196], [4, 269962, 0.44641420123071723]]

	chrom_col_map = {}
	for i in range(1, 6):
		if i % 2 == 0:
			chrom_col_map[i] = '#1199EE'
		else:
			chrom_col_map[i] = '#11BB00'

	#Draw plot
	f = pylab.figure(figsize=(7, 7))
	ax1 = f.add_axes([0.08, 0.07 + (0.667 * 0.94), 0.9, (0.3 * 0.94) ])
	ax1.spines['top'].set_visible(False)
	ax1.xaxis.set_ticks_position('bottom')
	ax1.xaxis.set_label_position('bottom')
	ax1.set_ylabel(r'$-$log$($p-value$)$')
	ax2 = f.add_axes([0.08, 0.07 + (0.333 * 0.94), 0.9, (0.3 * 0.94)  ])
	ax2.spines['top'].set_visible(False)
	ax2.set_ylabel('Posterior prob. of assoc.')
	ax2.xaxis.set_ticks_position('bottom')
	ax2.xaxis.set_label_position('bottom')
	ax4 = f.add_axes([0.08, 0.06 , 0.69, (0.24 * 0.94) ]) #Partition of variance axes
	#ax4.spines['top'].set_visible(False)
	ax4.xaxis.set_ticks_position('bottom')
	ax4.xaxis.set_label_position('bottom')
	ax4.set_ylabel('Partition of variance')
	ax4.set_xlabel('Step number')


	ex_res.plot_manhattan2(ax=ax1, neg_log_transform=True, plot_bonferroni=True, plot_xaxis=True, wo_xtick_labels=True,
				chrom_colormap=chrom_col_map, sign_color='#DD1122', cand_genes=cgs, percentile=98)
	x_min, x_max = ax1.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax1.get_ylim()
	y_range = y_max - y_min
	ax1.text(0.96 * x_range + x_min, 0.93 * y_range + y_min, 'a')


	ex_ppa_res.plot_manhattan2(ax=ax2, neg_log_transform=False, plot_bonferroni=False, plot_xaxis=True,
				chrom_colormap=chrom_col_map, highlight_markers=ppa_cofactors, sign_color=None,
				cand_genes=cgs, max_score=1, percentile=98)
	x_min, x_max = ax2.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax2.get_ylim()
	y_range = y_max - y_min
	ax2.text(0.96 * x_range + x_min, 0.94 * y_range + y_min, 'b')
	ax2.text(0.629 * x_range + x_min, 0.767 * y_range + y_min, '1', fontsize='small')
	ax2.text(0.629 * x_range + x_min, 0.43 * y_range + y_min, '2', fontsize='small')

	max_rss = max(rss_list)
	rss_array = sp.array(rss_list) / max_rss
	p_her_array = rss_array * sp.array(p_her_list)
	genetic_variance = p_her_array + (1 - rss_array)
	variance_explained = (1 - rss_array)
	ax4.fill_between([0, num_steps], 0, 1.1, color='#FFCC33', alpha=0.8, label='Error')
	ax4.fill_between(sp.arange(num_steps + 1), 0, genetic_variance, color='#22CC44', alpha=0.8, label='Genetic')
	ax4.fill_between(sp.arange(num_steps + 1), -0.1, variance_explained, color='#2255AA', alpha=0.8, label='Explained')
	ax4.axvline(x=num_steps / 2, c='k', linestyle=':')
	ax4.axis([0, num_steps, 0, 1])
	x_min, x_max = ax4.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax4.get_ylim()
	y_range = y_max - y_min
	ax4.text(1.25 * x_range + x_min, 0.98 * y_range + y_min, 'c')
	ax4.bar(-1, 1, color='#FFCC33', alpha=0.8, label='Error')
	ax4.bar(-1, 1, color='#22CC44', alpha=0.8, label='Genetic')
	ax4.bar(-1, 1, color='#2255AA', alpha=0.8, label='Explained')
	ax4.axvline(2, color='#AA0000', alpha=0.8, lw=1.8)
	ax4.text(0.04 * x_range + x_min, 1.04 * y_range + y_min, 'FRI-indels', fontsize='small')
	leg = ax4.legend(loc='center left', bbox_to_anchor=(1, 0.5))

	ltext = leg.get_texts()  # all the text.Text instance in the legend
	llines = leg.get_lines()  # all the lines.Line2D instance in the legend
	frame = leg.get_frame()  # the patch.Rectangle instance surrounding the legend
	# see text.Text, lines.Line2D, and patches.Rectangle for more info on
	# the settable properties of lines, text, and rectangles
	frame.set_linewidth(0)
	#frame.set_facecolor('0.80')      # set the frame face color to light gray
	pylab.setp(ltext, fontsize='small')    # the legend text fontsize
	pylab.setp(llines, linewidth=1.5)      # the legend linewidth



	pylab.savefig(env.env['tmp_dir'] + 'mlt_bayesian_example.png', dpi=400)
	pylab.savefig(env.env['tmp_dir'] + 'mlt_bayesian_example.pdf', format='pdf')



def plot_remaining_bayesian_figures():
	with open(env.env['tmp_dir'] + 'bf_qtl_priors_DTF1stSpAverage2009_314_stats.csv') as f:
		lines = f.readlines()
		rss_list1 = []
		p_her_list1 = []
		for line in lines[1:]:
			l = map(str.strip, line.split(','))
			rss_list1.append(float(l[1]))
			p_her_list1.append(float(l[9]))
		num_steps = len(lines) - 1
		rss_list1.append(rss_list1[0])
		p_her_list1.append(p_her_list1[0])
	ppa_cofactors1 = [[5, 3188327, 0.99631213930493612]]

	with open(env.env['tmp_dir'] + 'bf_qtl_priors_DTF2ndSpAverage2009_315_stats.csv') as f:
		lines = f.readlines()
		rss_list2 = []
		p_her_list2 = []
		for line in lines[1:]:
			l = map(str.strip, line.split(','))
			rss_list2.append(float(l[1]))
			p_her_list2.append(float(l[9]))
		num_steps2 = len(lines) - 1
		rss_list2.append(rss_list2[0])
		p_her_list2.append(p_her_list2[0])
	ppa_cofactors2 = [[5, 3188327, 0.99846328919628757]]

	with open(env.env['tmp_dir'] + 'bf_qtl_priors_DTF1stSwAverage2009_316_stats.csv') as f:
		lines = f.readlines()
		rss_list3 = []
		p_her_list3 = []
		for line in lines[1:]:
			l = map(str.strip, line.split(','))
			rss_list3.append(float(l[1]))
			p_her_list3.append(float(l[9]))
		num_steps3 = len(lines) - 1
		rss_list3.append(rss_list3[0])
		p_her_list3.append(p_her_list3[0])
	ppa_cofactors3 = [[4, 161496, 0.99931725487952472], [4, 269260, 0.9914874431415277]]

	f = pylab.figure(figsize=(9, 7))
	ax1 = f.add_axes([0.08, 0.07 + (0.667 * 0.94), 0.9, (0.3 * 0.94) ])
	ax1.spines['top'].set_visible(False)
	ax1.xaxis.set_ticks_position('bottom')
	ax1.xaxis.set_label_position('bottom')
	ax1.set_ylabel('Posterior prob. of assoc.')
	ax2 = f.add_axes([0.08, 0.07 + (0.333 * 0.94), 0.9, (0.3 * 0.94) ])
	ax2.spines['top'].set_visible(False)
	ax2.set_ylabel('Posterior prob. of assoc.')
	ax2.xaxis.set_ticks_position('bottom')
	ax2.xaxis.set_label_position('bottom')
	ax3 = f.add_axes([0.08, 0.07, 0.9, (0.3 * 0.94) ])
	ax3.spines['top'].set_visible(False)
	ax3.set_ylabel('Posterior prob. of assoc.')
	ax3.xaxis.set_ticks_position('bottom')
	ax3.xaxis.set_label_position('bottom')

	ex_ppa_res1 = gr.Result(result_file=env.env['tmp_dir'] + 'bf_qtl_priors_DTF1stSpAverage2009_314_step1.ppas')
	ex_ppa_res2 = gr.Result(result_file=env.env['tmp_dir'] + 'bf_qtl_priors_DTF2ndSpAverage2009_315_step1.ppas')
	ex_ppa_res3 = gr.Result(result_file=env.env['tmp_dir'] + 'bf_qtl_priors_DTF1stSwAverage2009_316_opt_mbonf_step18.ppas')

	chrom_col_map = {}
	for i in range(1, 6):
		if i % 2 == 0:
			chrom_col_map[i] = '#1199EE'
		else:
			chrom_col_map[i] = '#11BB00'


	ex_ppa_res1.plot_manhattan2(ax=ax1, neg_log_transform=False, plot_bonferroni=False, plot_xaxis=True,
				chrom_colormap=chrom_col_map, highlight_markers=ppa_cofactors1, sign_color=None)
	x_min, x_max = ax1.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax1.get_ylim()
	y_range = y_max - y_min
	ax1.text(0.96 * x_range + x_min, 0.94 * y_range + y_min, 'a')

	ex_ppa_res2.plot_manhattan2(ax=ax2, neg_log_transform=False, plot_bonferroni=False, plot_xaxis=True,
				chrom_colormap=chrom_col_map, highlight_markers=ppa_cofactors2, sign_color=None)
	x_min, x_max = ax2.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax2.get_ylim()
	y_range = y_max - y_min
	ax3.text(0.96 * x_range + x_min, 0.94 * y_range + y_min, 'b')

	ex_ppa_res3.plot_manhattan2(ax=ax3, neg_log_transform=False, plot_bonferroni=False, plot_xaxis=True,
				chrom_colormap=chrom_col_map, highlight_markers=ppa_cofactors3, sign_color=None)
	x_min, x_max = ax3.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax3.get_ylim()
	y_range = y_max - y_min
	ax3.text(0.96 * x_range + x_min, 0.94 * y_range + y_min, 'c')

	pylab.savefig(env.env['tmp_dir'] + 'mlt_remaining_bayesian_examples.png', dpi=400)



def plot_all_human_figures():
	"""Plots EBIC plot
	"""
	import gwaResults as gr
	chrom_d, var_d = dp._parse_map_file_()
	human_dir = env.env['data_dir'] + 'NFBC_results/multiple_loci/'
	file_list = [#('BMI', 'bmires_pval_steps.csv', 'bmires_steptable_ok.csv', 0, 'a'),
			('CRP', 'crpres_pval_steps.csv', 'crpres_steptable_ok.csv', 3, 'a'),
			#('DBP', 'diares_pval_steps.csv', 'diares_steptable_ok.csv', 0, 'c'),
			('Glucose', 'glures_pval_steps.csv', 'glures_steptable_ok.csv', 3, 'b'),
			('HDL', 'hdlres_pval_steps.csv', 'hdlres_steptable_ok.csv', 3, 'c'),
			#('Height', 'heightres_pval_steps.csv', 'heightres_steptable_ok.csv', 0, 'a'),
			#('Insulin', 'insres_pval_steps.csv', 'insres_steptable_ok.csv', 0, 'b'),
			#('LDL', 'ldlres_pval_steps.csv', 'ldlres_steptable_ok.csv', 5, 'c'),
			('SBP', 'sysres_pval_steps.csv', 'sysres_steptable_ok.csv', 0, 'd'),
			('TG', 'tgres_pval_steps.csv', 'tgres_steptable_ok.csv', 2, 'e'),
			]

	chrom_col_map = {}
	for i in range(1, 24):
		if i % 2 == 0:
			chrom_col_map[i] = '#1199EE'
		else:
			chrom_col_map[i] = '#11BB00'

	for name, pval_file, step_info_file, step_nr, fig_id in file_list:
		#One plot per file
		pvals = []
		chromosomes = []
		positions = []

		with open(human_dir + pval_file) as f:
			print f.next()
			for line in f:
				l = line.split(',')
				var_id = l[0].replace('"', '')
				d = var_d[var_id]
				chromosomes.append(d['chrom'])
				positions.append(d['position'])
				pvals.append(float(l[step_nr + 1]))
		l = zip(chromosomes, positions, pvals)
		l.sort()
		l = map(list, zip(*l))
		res = gr.Result(scores=l[2], chromosomes=l[0], positions=l[1])
		res.neg_log_trans()
		cofactors = []
		if step_nr > 0:
			with open(human_dir + step_info_file) as f:
				print f.next()
				print f.next()
				for i in range(step_nr):
					line = f.next()
					l = line.split(',')
					var_id = (l[2].replace('"', ''))[1:]
					d = var_d[var_id]
					cofactors.append((d['chrom'], d['position']))

		f = pylab.figure(figsize=(10, 2.5))
		ax = f.add_axes([0.06, 0.14, 0.93, 0.83])
		ax.spines['top'].set_visible(False)
		ax.xaxis.set_ticks_position('bottom')
		ax.xaxis.set_label_position('bottom')
		ax.set_ylabel(r'$-$log$($p-value$)$')

		res.plot_manhattan2(ax=ax, neg_log_transform=False, plot_bonferroni=True, plot_xaxis=True,
					chrom_colormap=chrom_col_map, highlight_markers=cofactors,
					sign_color='#DD1122', max_score=32, percentile=98)
		x_min, x_max = ax.get_xlim()
		x_range = x_max - x_min
		y_min, y_max = ax.get_ylim()
		y_range = y_max - y_min
		ax.text(0.96 * x_range + x_min, 0.94 * y_range + y_min, fig_id)
		pylab.savefig(env.env['tmp_dir'] + 'mlt_HUMAN_%s.png' % name, dpi=400)
		pylab.savefig(env.env['tmp_dir'] + 'mlt_HUMAN_%s.pdf' % name, format='pdf')




def plot_local_sodium_figure():
	cg_tair_id = 'AT4G10310' #6391984 - 6395877 bp
	cgs = [cg_tair_id]#gr.get_genes_w_tair_id([cg_tair_id])
	results_dir = env.env['tmp_dir']
	result_files = [results_dir + 'abracadabra_pid226_Na23_Soil_3_emmax_step_log_t75_4_6291984_6495877_step0.pvals',
			results_dir + 'abracadabra_pid226_Na23_Soil_3_emmax_step_log_t75_4_6291984_6495877_step6.pvals']
	results = []
	for result_file in result_files:
		res = gr.Result(result_file)
		res.neg_log_trans()
		results.append(res)

	with open(results_dir + 'abracadabra_pid226_Na23_Soil_3_emmax_step_log_t75_4_6291984_6495877_stats.csv') as f:
		lines = f.readlines()
		l6 = map(str.strip, lines[7].split(','))
		c_list = l6[14:]
		cofactors = []
		for c_p_p in c_list:
			l = c_p_p.split('_')
			cofactors.append((int(l[0]), int(l[1]), float(l[2])))
		rss_list = []
		p_her_list = []
		for line in lines[1:]:
			l = map(str.strip, line.split(','))
			rss_list.append(float(l[1]))
			p_her_list.append(float(l[9]))
		num_steps = len(lines) - 1
		rss_list.append(rss_list[0])
		p_her_list.append(p_her_list[0])

	print cofactors

	chrom_col_map = {}
	for i in range(1, 6):
		if i % 2 == 0:
			chrom_col_map[i] = '#1199EE'
		else:
			chrom_col_map[i] = '#11BB00'

	#Draw plot
	f = pylab.figure(figsize=(9, 10))
	ax1 = f.add_axes([0.08, 0.07 + (0.65 * 0.94), 0.9, (0.28 * 0.94) ])
	ax1.spines['top'].set_visible(False)
	ax1.xaxis.set_ticks_position('bottom')
	ax1.xaxis.set_label_position('bottom')
	ax1.set_ylabel(r'$-$log$($p-value$)$')
	ax2 = f.add_axes([0.08, 0.07 + (0.33 * 0.94), 0.9, (0.29 * 0.94) ])
	ax2.spines['top'].set_visible(False)
	ax2.set_ylabel(r'$-$log$($p-value$)$')
	ax2.xaxis.set_ticks_position('bottom')
	ax2.xaxis.set_label_position('bottom')
	ax3 = f.add_axes([0.08, 0.07 + (0.19 * 0.94), 0.9, (0.09 * 0.94) ]) #Gene model axes
	ax4 = f.add_axes([0.08, 0.05 , 0.73, (0.17 * 0.94) ]) #Partition of variance axes
	#ax3.spines['top'].set_visible(False)
	ax4.xaxis.set_ticks_position('bottom')
	ax4.xaxis.set_label_position('bottom')
	ax4.set_ylabel('Partition of variance')
	ax4.set_xlabel('Step number')


	#Fill up the figure..
	results[0]._plot_small_manhattan_2_(plot_ax=ax1, plot_bonferroni=True, chromosome=4,
				color_map=chrom_col_map, sign_color='#DD1122', xlab='',
				cand_genes=cgs, max_score=15.5)
	x_min, x_max = ax1.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax1.get_ylim()
	y_range = y_max - y_min
	ax1.text(0.96 * x_range + x_min, 0.85 * y_range + y_min, 'a')

	results[1]._plot_small_manhattan_2_(plot_ax=ax2, gene_model_ax=ax3, plot_bonferroni=True, chromosome=4,
					color_map=chrom_col_map, highlight_markers=cofactors,
					sign_color='#DD1122', title='', cand_genes=cgs, max_score=15.5)
	x_min, x_max = ax2.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax2.get_ylim()
	y_range = y_max - y_min
	ax2.text(0.96 * x_range + x_min, 0.85 * y_range + y_min, 'b')


	max_rss = max(rss_list)
	rss_array = sp.array(rss_list) / max_rss
	p_her_array = rss_array * sp.array(p_her_list)
	genetic_variance = p_her_array + (1 - rss_array)
	variance_explained = (1 - rss_array)
	ax4.fill_between([0, num_steps], 0, 1.1, color='#FFCC33', alpha=0.8, label='Error')
	ax4.fill_between(sp.arange(num_steps + 1), 0, genetic_variance, color='#22CC44', alpha=0.8, label='Genetic')
	ax4.fill_between(sp.arange(num_steps + 1), -0.1, variance_explained, color='#2255AA', alpha=0.8, label='Explained')
	ax4.axvline(6, color='#AA0000', alpha=0.8, lw=1.8)
	ax4.axvline(5, color='#AA0000', alpha=0.8, lw=1.8)
	ax4.axvline(x=num_steps / 2, c='k', linestyle=':')
	ax4.axis([0, num_steps, 0, 1])
	x_min, x_max = ax4.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax4.get_ylim()
	y_range = y_max - y_min
	ax4.text(0.2 * x_range + x_min, 1.04 * y_range + y_min, 'MBONF', fontsize='small')
	ax4.text(0.29 * x_range + x_min, 1.04 * y_range + y_min, 'EBIC', fontsize='small')
	ax4.text(1.185 * x_range + x_min, 0.9 * y_range + y_min, 'c')
	ax4.bar(-1, 1, color='#FFCC33', alpha=0.8, label='Error')
	ax4.bar(-1, 1, color='#22CC44', alpha=0.8, label='Genetic')
	ax4.bar(-1, 1, color='#2255AA', alpha=0.8, label='Explained')
	leg = ax4.legend(loc='center left', bbox_to_anchor=(1, 0.5))

	ltext = leg.get_texts()  # all the text.Text instance in the legend
	llines = leg.get_lines()  # all the lines.Line2D instance in the legend
	frame = leg.get_frame()  # the patch.Rectangle instance surrounding the legend

	# see text.Text, lines.Line2D, and patches.Rectangle for more info on
	# the settable properties of lines, text, and rectangles
	frame.set_linewidth(0)
	#frame.set_facecolor('0.80')      # set the frame face color to light gray
	pylab.setp(ltext, fontsize='small')    # the legend text fontsize
	pylab.setp(llines, linewidth=1.5)      # the legend linewidth


	#Save the figure?
	pylab.savefig(env.env['tmp_dir'] + 'mlt_local_sodium_example.png', dpi=400)
	pylab.savefig(env.env['tmp_dir'] + 'mlt_local_sodium_example.pdf', format='pdf')


def plot_first_second_power_fdr(heritability, fdr_window_size=25000, pdf_file_name='/tmp/first_second_power.pdf',
				png_file_name='/tmp/first_second_power.png'):
	"""
	Plot the power-FDR plot for the first and second markers.
	"""
	#file_prefix = '/srv/lab/data/mlt_results/mlt'  #Because that's where the results are
	file_prefix = '/home/GMI/bjarni.vilhjalmsson/results/mlt'  #Because that's where the results are
	p_dict, args = parse_parameters(min_num_parameters=0)
	phen_d = load_phenotypes(p_dict['phen_file']) #Loading simulated phenotypes
	latent_var = 'random_snp'  #This is only intersting for this case
	ws_i = window_sizes.index(fdr_window_size)  #This only works for window size 0 at the moment.  (Which is not fair)

	fig = pylab.figure(figsize=(11, 3.8))
	ax1 = fig.add_axes([0.06, 0.11, (0.3 * 0.97), 0.81 ])
	ax2 = fig.add_axes([0.06 + (0.33 * 0.97), 0.11 , (0.3 * 0.97), 0.81 ])
	ax2.set_yticklabels([''] * len(ax2.get_yticks()))
	ax3 = fig.add_axes([0.054 + (0.667 * 0.97), 0.11, (0.3 * 0.97), 0.81 ])
	ax3.set_yticklabels([''] * len(ax3.get_yticks()))

	efig = pylab.figure(figsize=(11, 3.8))
	eax1 = efig.add_axes([0.06, 0.11, (0.3 * 0.97), 0.81 ])

	eax2 = efig.add_axes([0.06 + (0.33 * 0.97), 0.11 , (0.3 * 0.97), 0.81 ])
	eax2.set_yticklabels([''] * len(eax2.get_yticks()))

	eax3 = efig.add_axes([0.06 + (0.667 * 0.97), 0.11, (0.3 * 0.97), 0.81 ])
	eax3.set_yticklabels([''] * len(eax3.get_yticks()))

	for phen_model, ax, eax in zip(['plus', 'or', 'xor'], [ax1, ax2, ax3], [eax1, eax2, eax3]):
		pd = phen_d[latent_var][heritability][phen_model]
		index_list = range(len(pd['phenotypes']))

		#Load results (for EX and SW)
		print '%s_%d_%s_%s' % (file_prefix, heritability, latent_var, phen_model)

		#Set up FDR and Power datastructures
		plot_dict = {'EX':{'fdrs':[[] for pt in pval_thresholds], 'f_pow':[[] for pt in pval_thresholds],
				's_pow':[[] for pt in pval_thresholds]},
			'SW_MBONF':{'fdrs':[[] for pt in pval_thresholds], 'f_pow':[[] for pt in pval_thresholds],
				's_pow':[[] for pt in pval_thresholds]},
			'SW_EBIC':{'fdrs':[], 'f_pow':[], 's_pow':[]},
			'SW_BONF':{'fdrs':[], 'f_pow':[], 's_pow':[]},
			'EX_BONF':{'fdrs':[], 'f_pow':[], 's_pow':[]},
			'efs':[], 'ess':[]}

		num_files_found = 0
		#For each SW,EX result
		for i in index_list:
			pickled_file = '%s_%d_%s_%s_%dresults.pickled' % (file_prefix, heritability, latent_var, phen_model, i)
			if not os.path.isfile(pickled_file):
				#print "Didn't find pickle file: %s" % pickled_file
				continue

			num_files_found += 1
			with open(pickled_file) as f:
				r = cPickle.load(f)

			ex_dict = r['EX']
			sw_dict = r['Stepw_EX']
			t_opt_i_list = sw_dict['mbonf']['t_opt_i_list']
			max_i = len(sw_dict['step_info_list']) - 1

			#Decide the order of the SNPs
			f_i = 1 if ex_dict['causal_ranks'][0] > ex_dict['causal_ranks'][1] else 0
			s_i = abs(f_i - 1)

			effect_sizes = phen_d['random_snp'][heritability][phen_model]['h_estimates']
			plot_dict['efs'].append(effect_sizes[f_i])
			plot_dict['ess'].append(effect_sizes[s_i])
			#Power for first and second causal over total FDR
			#For each threshold and window	
			for pt_i, pt in enumerate(pval_thresholds):
				#For EX
				plot_dict['EX']['f_pow'][pt_i].append(ex_dict['causal_tpr_fdr'][pt_i][ws_i][f_i])
				plot_dict['EX']['s_pow'][pt_i].append(ex_dict['causal_tpr_fdr'][pt_i][ws_i][s_i])
				plot_dict['EX']['fdrs'][pt_i].append(ex_dict['fdrs'][pt_i][ws_i])

				#Now for SW
				#MBONF
				opt_i = t_opt_i_list[pt_i]
				if opt_i == 0 or opt_i == max_i:
					plot_dict['SW_MBONF']['fdrs'][pt_i].append(0.0)
					plot_dict['SW_MBONF']['f_pow'][pt_i].append(0.0)
					plot_dict['SW_MBONF']['s_pow'][pt_i].append(0.0)
					continue
				cofactors = sw_dict['step_info_list'][opt_i]['cofactors']
				cofactors = sp.array([(c[0], c[1]) for c in cofactors])
				if sw_dict['mbonf']['fdrs'] == [-1] * len(window_sizes):
					plot_dict['SW_MBONF']['fdrs'][pt_i].append(0)
				else:
					plot_dict['SW_MBONF']['fdrs'][pt_i].append(sw_dict['mbonf']['fdrs'][pt_i][ws_i])
				plot_dict['SW_MBONF']['f_pow'][pt_i].append(sw_dict['mbonf']['causal_tpr_fdr'][pt_i][ws_i][f_i])
				plot_dict['SW_MBONF']['s_pow'][pt_i].append(sw_dict['mbonf']['causal_tpr_fdr'][pt_i][ws_i][s_i])

			#Now for SW
			#EBIC
			sw_ebic_dict = r['Stepw_EX']['ebics']
			opt_i = sw_ebic_dict['opt_i']
			if opt_i == 0 or opt_i == max_i:
				plot_dict['SW_EBIC']['fdrs'].append(0.0)
				plot_dict['SW_EBIC']['f_pow'].append(0.0)
				plot_dict['SW_EBIC']['s_pow'].append(0.0)
				continue
			cofactors = sw_dict['step_info_list'][opt_i]['cofactors']
			cofactors = sp.array([(c[0], c[1]) for c in cofactors])
			plot_dict['SW_EBIC']['fdrs'].append(sw_ebic_dict['fdrs'][ws_i])
			plot_dict['SW_EBIC']['f_pow'].append(sw_ebic_dict['causal_tpr_fdr'][ws_i][f_i])
			plot_dict['SW_EBIC']['s_pow'].append(sw_ebic_dict['causal_tpr_fdr'][ws_i][s_i])


		#Now bin the effects
#			'efs':[[] for i in range(10)], 'efs_counts':[[] for i in range(10)],
#			'ess':[[] for i in range(10)], 'ess_counts':[[] for i in range(10)]}


		#Now take the averages
		for pt_i, pt in enumerate(pval_thresholds):
			plot_dict['EX']['f_pow'][pt_i] = sp.mean(plot_dict['EX']['f_pow'][pt_i])
			plot_dict['EX']['s_pow'][pt_i] = sp.mean(plot_dict['EX']['s_pow'][pt_i])
			plot_dict['EX']['fdrs'][pt_i] = sp.mean(plot_dict['EX']['fdrs'][pt_i])
			plot_dict['SW_MBONF']['f_pow'][pt_i] = sp.mean(plot_dict['SW_MBONF']['f_pow'][pt_i])
			plot_dict['SW_MBONF']['s_pow'][pt_i] = sp.mean(plot_dict['SW_MBONF']['s_pow'][pt_i])
			plot_dict['SW_MBONF']['fdrs'][pt_i] = sp.mean(plot_dict['SW_MBONF']['fdrs'][pt_i])
		plot_dict['SW_EBIC']['f_pow'] = sp.mean(plot_dict['SW_EBIC']['f_pow'])
		plot_dict['SW_EBIC']['s_pow'] = sp.mean(plot_dict['SW_EBIC']['s_pow'])
		plot_dict['SW_EBIC']['fdrs'] = sp.mean(plot_dict['SW_EBIC']['fdrs'])
		plot_dict['SW_BONF']['f_pow'] = plot_dict['SW_MBONF']['f_pow'][31] * 0.56 + plot_dict['SW_MBONF']['f_pow'][32] * 0.44
		plot_dict['SW_BONF']['s_pow'] = plot_dict['SW_MBONF']['s_pow'][31] * 0.56 + plot_dict['SW_MBONF']['s_pow'][32] * 0.44
		plot_dict['SW_BONF']['fdrs'] = plot_dict['SW_MBONF']['fdrs'][31] * 0.56 + plot_dict['SW_MBONF']['fdrs'][32] * 0.44
		plot_dict['EX_BONF']['f_pow'] = plot_dict['EX']['f_pow'][31] * 0.56 + plot_dict['EX']['f_pow'][32] * 0.44
		plot_dict['EX_BONF']['s_pow'] = plot_dict['EX']['s_pow'][31] * 0.56 + plot_dict['EX']['s_pow'][32] * 0.44
		plot_dict['EX_BONF']['fdrs'] = plot_dict['EX']['fdrs'][31] * 0.56 + plot_dict['EX']['fdrs'][32] * 0.44

		#Now plot the plot 
		y_label = ''
		w_legend = False
		if phen_model == 'plus':
			y_label = 'Power (for the SNPs)'
			w_legend = True
		plot_single_tprs_fdrs_3(plot_dict, ax, y_lim=(0, 1.11), y_label=y_label, w_legend=w_legend)

		#Now plot the effects plot 



	fig.text(0.174, 0.95, 'Additive')
	fig.text(0.511, 0.95, "'or'")
	fig.text(0.828, 0.95, "'xor'")

	fig.text(0.072, 0.861, 'a')
	fig.text(0.392, 0.861, "b")
	fig.text(0.713, 0.861, "c")

	fig.savefig(png_file_name)
	fig.savefig(pdf_file_name, format='pdf')


def plot_single_tprs_fdrs_3(dict, ax, w_legend=False, y_label='Power (for the SNPs)', x_label='FDR', y_lim=None):
	"""
	Plots the first-second power plots..
	"""
	import matplotlib.font_manager
	prop = matplotlib.font_manager.FontProperties(size=10)

	#Heritabilities..
	# - histogram of each category
	# - pseudoheritabilities vs. ks and med pval. of KW and LM

	# TPRs vs. FDRs
	xss = [dict['EX']['fdrs'], dict['SW_MBONF']['fdrs'], dict['EX']['fdrs'], dict['SW_MBONF']['fdrs']]
	yss = [dict['EX']['f_pow'], dict['SW_MBONF']['f_pow'], dict['EX']['s_pow'], dict['SW_MBONF']['s_pow']]
	am_labels = ['MM first SNP', 'MLMM first SNP', 'MM second SNP', 'MLMM second SNP']
	am_colors = ['#0022FF', '#0022FF', '#FF0022', '#FF0022']
	am_ls = ['--', '-', '--', '-', ]
	am_dot_list = ['Stepw_EX_EBIC']
	am_dot_colors = ['#FF0022']
	#am_dot_labels = ['MLML EBIC', 'SWLM EBIC']

	#First the plot
	for xs, ys, amc, amls, am_label in zip(xss, yss, am_colors, am_ls, am_labels):
		ax.plot(xs, ys, ls=amls, color=amc, label=am_label, alpha=0.65, lw=2)
	#Now the dots
	ax.plot([dict['SW_EBIC']['fdrs']], [dict['SW_EBIC']['f_pow']], marker='v',
		ls='', color='#0022FF', alpha=0.65)
	ax.plot([dict['SW_EBIC']['fdrs']], [dict['SW_EBIC']['s_pow']], marker='v',
		ls='', color='#FF0022', alpha=0.65)
	ax.plot([dict['SW_BONF']['fdrs']], [dict['SW_BONF']['f_pow']], marker='o',
		ls='', color='#0022FF', alpha=0.65)
	ax.plot([dict['SW_BONF']['fdrs']], [dict['SW_BONF']['s_pow']], marker='o',
		ls='', color='#FF0022', alpha=0.65)
	ax.plot([dict['EX_BONF']['fdrs']], [dict['EX_BONF']['f_pow']], marker='o',
		ls='', color='#0022FF', alpha=0.65)
	ax.plot([dict['EX_BONF']['fdrs']], [dict['EX_BONF']['s_pow']], marker='o',
		ls='', color='#FF0022', alpha=0.65)
	ax.set_ylabel(y_label)
	ax.set_xlabel(x_label)
	if y_lim:
		ax.set_ylim(y_lim)
	x_min, x_max = ax.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax.get_ylim()
	y_range = y_max - y_min
	ax.plot(x_min - x_range, y_min - y_range, color='k', marker='v', label='EBIC', alpha=0.65)
	ax.plot(x_min - x_range, y_min - y_range, color='k', marker='o', label='Bonferroni', alpha=0.65)
	if w_legend:
		ax.legend(loc=4, prop=prop, numpoints=1, scatterpoints=1)
	ax.axis([x_min - 0.025 * x_range, x_max + 0.025 * x_range,
			y_min - 0.025 * y_range, y_max + 0.025 * y_range])






def generate_example_figure_7_old():
	import gwaResults as gr

	result_file_prefix = env.env['results_dir'] + 'bf_most_norm_'
	cg_tair_ids = ['AT1G04400', 'AT1G65480', 'AT1G77080', 'AT2G26330', 'AT4G00650', 'AT5G10140', 'AT5G65050', \
			'AT5G65060', 'AT5G65070', 'AT5G65080'] \
			#FT, MAF, ER, FRI, FLC, MAF2-MAF5 (Salome et al. 2011)
	cgs = gr.get_genes_w_tair_id(cg_tair_ids)
	result_files = [result_file_prefix + fn for fn in ['DTF1stSpAverage2009_314_step1.ppas', \
							'DTF2ndSpAverage2009_315_opt_min_cof_ppa_step13.ppas', \
							'DTF1stSwAverage2009_316_step3.ppas', \
							'DTF2ndSwAverage2009_317_step4.ppas']]
	#Load pickle file...
	results = [gr.Result(result_file=fn) for fn in result_files]
	#Setting up figure
	f = pylab.figure(figsize=(9, 8))
	ax1 = f.add_axes([0.08, 0.07 + (0.75 * 0.94), 0.9, (0.225 * 0.94) ])
	ax1.spines['top'].set_visible(False)
	ax1.xaxis.set_visible(False)
	ax2 = f.add_axes([0.08, 0.07 + (0.5 * 0.94), 0.9, (0.225 * 0.94) ])
	ax2.spines['top'].set_visible(False)
	ax2.xaxis.set_visible(False)
	ax3 = f.add_axes([0.08, 0.07 + (0.25 * 0.94), 0.9, (0.225 * 0.94) ])
	ax3.spines['top'].set_visible(False)
	ax3.xaxis.set_visible(False)
	ax4 = f.add_axes([0.08, 0.07 , 0.9, (0.225 * 0.94) ])
	ax4.spines['top'].set_visible(False)
	ax4.xaxis.set_ticks_position('bottom')
	ax4.xaxis.set_label_position('bottom')
	res = results[0]
	chrom_ends = res.get_chromosome_ends()
	offset = 0
	tick_positions = []
	tick_labels = []
	for c_end in chrom_ends[:-1]:
		offset += c_end
		tick_positions.append(offset)
		tick_labels.append('')
	ax4.set_xticks(tick_positions)
	ax4.set_xticklabels(tick_labels)
	y_tick_positions = [0, 0.5, 1.0]
	y_tick_labels = ['0.0', '0.5', '1.0']

	max_y = offset + chrom_ends[-1]

	for ax in [ax1, ax2, ax3, ax4]:
		ax.set_ylim((-0.05, 1.05))
		ax.set_xlim((-0.02 * max_y, 1.02 * max_y))
		ax.set_yticks(y_tick_positions)
		ax.set_yticklabels(y_tick_labels)


	#Fill up the figure..
	cm = {1:'#1199EE', 2:'#11BB00', 3:'#1199EE', 4:'#11BB00', 5:'#1199EE'}
	results[0].plot_manhattan2(ax=ax1, neg_log_transform=False, plot_bonferroni=False,
				chrom_colormap=cm, highlight_markers=[(5, 3188327, 0.99957464909496818)],
				cand_genes=cgs)
	x_min, x_max = ax1.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax1.get_ylim()
	y_range = y_max - y_min
	ax1.text(0.95 * x_range + x_min, 0.85 * y_range + y_min, 'A')

	results[1].plot_manhattan2(ax=ax2, neg_log_transform=False, plot_bonferroni=False,
				chrom_colormap=cm, highlight_markers=[(5, 3188327, 0.99987154081342355), (4, 161496, 0.99792660595760496), (1, 24341345, 0.99816162117265428), (2, 8509438, 0.99185074339847878), (4, 633609, 0.87069974648101967), (4, 7232003, 0.70642884690410201), (1, 1630085, 0.62628047124132957)],
				cand_genes=cgs)
	x_min, x_max = ax2.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax2.get_ylim()
	y_range = y_max - y_min
	ax2.text(0.95 * x_range + x_min, 0.85 * y_range + y_min, 'B')

	results[2].plot_manhattan2(ax=ax3, neg_log_transform=False, plot_bonferroni=False,
				chrom_colormap=cm, highlight_markers=[(4, 1356197, 0.86299217462489453), (4, 493905, 0.95488774939952203), (5, 3188327, 0.5687697089878565)],
				cand_genes=cgs)
	x_min, x_max = ax3.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax3.get_ylim()
	y_range = y_max - y_min
	ax3.text(0.95 * x_range + x_min, 0.85 * y_range + y_min, 'C')

	results[3].plot_manhattan2(ax=ax4, neg_log_transform=False, plot_bonferroni=False,
				chrom_colormap=cm, highlight_markers=[(5, 3188327, 0.99975387199431631), (4, 493905, 0.98796470842089601), (5, 7476090, 0.78288774499946212), (1, 24341345, 0.62528073836407139)],
				cand_genes=cgs)
	x_min, x_max = ax4.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax4.get_ylim()
	y_range = y_max - y_min
	ax4.text(0.95 * x_range + x_min, 0.85 * y_range + y_min, 'D')

	f.text(0.195, 0.04, '1')
	f.text(0.25, 0.04, 'FT')
	f.text(0.381, 0.04, '2')
	f.text(0.542, 0.04, '3')
	f.text(0.65, 0.04, 'FRI')
	f.text(0.704, 0.04, '4')
	f.text(0.82, 0.04, 'FLC')
	f.text(0.873, 0.04, '5')
	f.text(0.43, 0.01, 'Chromosome number')

	#Save the figure?
	pylab.savefig(env.env['tmp_dir'] + 'test.png')




def benchmark_human_emmax():
	import dataParsers as dp
	import phenotypeData as pd
	import env
	import gwaResults as gr
	import random
	import sys
	#What to benchmark:
	plink_prefix = env.env['phen_dir'] + 'NFBC_20091001/NFBC_20091001'
	#Calculation of kinship.
	sd = dp.parse_plink_tped_file(plink_prefix)
	individs = sd.accessions[:]
	phed = pd.parse_phenotype_file(env.env['phen_dir'] + 'NFBC_20091001/phenotype.csv')
	sd.coordinate_w_phenotype_data(phed, 1)
#	sd.filter_mac_snps(50)
	s0 = time.time()
	K = sd.get_ibd_kinship_matrix(num_dots=10)
	secs = time.time() - s0
	if secs > 60:
		mins = int(secs) / 60
		secs = secs - mins * 60
		print 'Took %d mins and %f seconds to calculate the kinship.' % (mins, secs)
	else:
		print 'Took %f seconds to calculate the kinship..' % (secs)

        snps = sd.getSnps()
        y = phed.get_values(1)
#        emmax_res = lm.emmax(snps, y, K)
#        res = gr.Result(scores=emmax_res['ps'].tolist(), snps_data=sd)

	s0 = time.time()
	emmax_res = lm.emmax_step_wise(y, K, sd=sd, num_steps=10)
	secs = time.time() - s0
	print 'EX stepwise took %0.2f seconds' % (secs)

#	s0 = time.time()
#	y = phed.get_values(1)
#	res = lm.get_emma_reml_estimates(y, K)
#	secs = time.time() - s0
#	if secs > 60:
#		mins = int(secs) / 60
#		secs = secs - mins * 60
#		print 'Took %d mins and %f seconds to calculate the variance components.' % (mins, secs)
#	else:
#		print 'Took %f seconds to calculate the variance components..' % (secs)
#	print res
	#GWAS scan.


def perform_human_emmax(pid=1):
	import dataParsers as dp
	import phenotypeData as pd
	import env
	import gwaResults as gr
	import random
	import sys
	s1 = time.time()
	plink_prefix = env.env['phen_dir'] + 'NFBC_20091001/NFBC_20091001'
	sd = dp.parse_plink_tped_file(plink_prefix)
	#sd.sample_snps(0.05)
	individs = sd.accessions[:]
	phed = pd.parse_phenotype_file(env.env['phen_dir'] + 'NFBC_20091001/phenotype.csv')
	#phed.filter_ecotypes(pid, random_fraction=0.2)
	sd.coordinate_w_phenotype_data(phed, pid)
	sd.filter_mac_snps(50)
	#K = sd.get_ibd_kinship_matrix(num_dots=10)
	K = lm.load_kinship_from_file('/home/GMI/bjarni.vilhjalmsson/Projects/data/NFBC_20091001/NFBC_20091001_kinship_diploid.ibs.pickled', sd.accessions)
	#K = prepare_k(K, individs, sd.accessions)
	phen_vals = phed.get_values(pid)
	phen_name = phed.get_name(pid)
	print 'Working on %s' % phen_name
	sys.stdout.flush()
	file_prefix = env.env['results_dir'] + 'NFBC_emmax_step_ibs_%s_pid%d' % (phen_name, pid)
	secs = time.time() - s1
	if secs > 60:
		mins = int(secs) / 60
		secs = secs - mins * 60
		print 'Took %d mins and %f seconds to load and preprocess the data.' % (mins, secs)
	else:
		print 'Took %f seconds to load and preprocess the data..' % (secs)
	chrom_col_map = {}
	for i in range(1, 24):
		if i % 2 == 0:
			chrom_col_map[i] = '#1199EE'
		else:
			chrom_col_map[i] = '#11BB00'
	emmax_res = lm.emmax_step_wise(phen_vals, K, sd=sd, num_steps=10, file_prefix=file_prefix, markersize=5,
				chrom_col_map=chrom_col_map, save_pvals=True)
#	snps = sd.getSnps()
#	emmax_res = emmax(snps, phen_vals, K)
#	res = gr.Result(scores=emmax_res['ps'].tolist(), snps_data=sd)
	#res.write_to_file(env.env['results_dir'] + 'NFBC_emmax_pid%d.pvals' % pid)
	#res.neg_log_trans()
	#res.plot_manhattan(png_file=file_prefix + '.png', plot_xaxis=False, plot_bonferroni=True)
	secs = time.time() - s1
	if secs > 60:
		mins = int(secs) / 60
		secs = secs - mins * 60
		print 'Took %d mins and %f seconds in total.' % (mins, secs)
	else:
		print 'Took %f seconds in total.' % (secs)



def perform_a_thal_emmax(pid=226):
	import dataParsers as dp
	import phenotypeData as pd
	import env
	import gwaResults as gr
	import random
	import sys
	s1 = time.time()
	sd = dp.load_snps_call_method(76)
	#sd.sample_snps(0.05)
	individs = sd.accessions[:]
	phed = pd.get_phenotypes_from_db([pid])
	phed.log_transform(pid)
	#phed.filter_ecotypes(pid, random_fraction=0.2)
	sd.coordinate_w_phenotype_data(phed, pid)
	sd.filter_mac_snps(15)
	K = sd.get_ibs_kinship_matrix(num_dots=10)
	#K = prepare_k(K, individs, sd.accessions)
	phen_vals = phed.get_values(pid)
	phen_name = phed.get_name(pid)
	print 'Working on %s' % phen_name
	sys.stdout.flush()
	file_prefix = env.env['results_dir'] + 'MLMM_t76_ibs_%s_pid%d' % (phen_name, pid)
	secs = time.time() - s1
	if secs > 60:
		mins = int(secs) / 60
		secs = secs - mins * 60
		print 'Took %d mins and %f seconds to load and preprocess the data.' % (mins, secs)
	else:
		print 'Took %f seconds to load and preprocess the data..' % (secs)
	chrom_col_map = {}
	for i in range(1, 6):
		if i % 2 == 0:
			chrom_col_map[i] = '#1199EE'
		else:
			chrom_col_map[i] = '#11BB00'
#	emmax_res = lm.emmax_step_wise(phen_vals, K, sd=sd, num_steps=10, file_prefix=file_prefix, markersize=5,
#				chrom_col_map=chrom_col_map)
	#snps = sd.getSnps()
	#emmax_res = emmax(snps, phen_vals, K)
	#res = gr.Result(scores=emmax_res['ps'].tolist(), snps_data=sd)
	#res.write_to_file(env.env['results_dir'] + 'NFBC_emmax_pid%d.pvals' % pid)
	#res.neg_log_trans()
	#res.plot_manhattan(png_file=file_prefix + '.png', plot_xaxis=False, plot_bonferroni=True)
	secs = time.time() - s1
	if secs > 60:
		mins = int(secs) / 60
		secs = secs - mins * 60
		print 'Took %d mins and %f seconds in total.' % (mins, secs)
	else:
		print 'Took %f seconds in total.' % (secs)

	cgs = gr.get_genes_w_tair_id(['AT4G10310'])
	s1 = time.time()
	sd = snpsdata.SNPsDataSet([sd.get_region_snpsd(4, 6360000, 6460000)],
				[4], data_format=sd.data_format)
	emmax_res = lm.emmax_step_wise(phen_vals, K, sd=sd, num_steps=0, file_prefix=file_prefix, markersize=5,
				chrom_col_map=chrom_col_map, local=True, cand_gene_list=cgs)
	#snps = sd.getSnps()
	#emmax_res = emmax(snps, phen_vals, K)
	#res = gr.Result(scores=emmax_res['ps'].tolist(), snps_data=sd)
	#res.write_to_file(env.env['results_dir'] + 'NFBC_emmax_pid%d.pvals' % pid)
	#res.neg_log_trans()
	#res.plot_manhattan(png_file=file_prefix + '.png', plot_xaxis=False, plot_bonferroni=True)
	secs = time.time() - s1
	if secs > 60:
		mins = int(secs) / 60
		secs = secs - mins * 60
		print 'Took %d mins and %f seconds in total.' % (mins, secs)
	else:
		print 'Took %f seconds in total.' % (secs)



def test_bayes_factor_enrichment():
	import phenotypeData as pd
	import env
	cpp_list = []
	with open('/Users/bjarni.vilhjalmsson/Projects/Data/DTF1.scan.tsv') as f:
		print f.next()
		for l in f:
			line = l.split()
			cpp_list.append((int(line[1]), int(line[2]), float(line[4])))

	for pid in [314, 315, 316, 317]:
		sd = dp.load_snps_call_method(75)
		#phed = pd.parse_phenotype_file(env.env['phen_dir'] + 'phen_raw_112210.csv')
		phed = pd.get_phenotypes_from_db([pid])
		phed.convert_to_averages()
		phed.transform(pid, 'most_normal')
		sd.coordinate_w_phenotype_data(phed, pid)
		sd.filter_mac_snps(15)
		phen_vals = phed.get_values(pid)
		phen_name = phed.get_name(pid)
		K = lm.scale_k(sd.get_ibs_kinship_matrix())

		#Candidate genes (TAIR IDs)
		cg_tair_ids = ['AT1G04400', 'AT1G65480', 'AT1G77080', 'AT2G26330', 'AT4G00650', 'AT5G10140', 'AT5G65050', \
				'AT5G65060', 'AT5G65070', 'AT5G65080'] \
				#FT, MAF, ER, FRI, FLC, MAF2-MAF5 (Salome et al. 2011)
		cgs = gr.get_genes_w_tair_id(cg_tair_ids)
		#snp_priors = sd.get_cand_genes_snp_priors(cgs, radius=25000, cg_prior_fold_incr=20)
		snp_priors = sd.get_snp_priors(cpp_list)#, cand_genes=cgs)
		lm.emmax_step_wise(phen_vals, K, sd, 10, env.env['tmp_dir'] + 'bf_qtl_priors_%s_%d' % (phen_name, pid),
				snp_priors=snp_priors, cand_gene_list=cgs, save_pvals=True, snp_choose_criteria='ppas')



def test_bayes_factor_enrichment_flc():
	"""
	FLC expression... etc.
	"""
	import phenotypeData as pd
	import env
	import pdb
	for pid in [43]:
		sd = dp.load_snps_call_method(32)
		#phed = pd.parse_phenotype_file(env.env['phen_dir'] + 'phen_raw_112210.csv')
		phed = pd.get_phenotypes_from_db([pid])
		phed.convert_to_averages()
		phed.transform(pid, 'most_normal')
		sd.coordinate_w_phenotype_data(phed, pid)
		sd.filter_mac_snps(15)
		phen_vals = phed.get_values(pid)
		phen_name = phed.get_name(pid)
		K = lm.scale_k(sd.get_ibs_kinship_matrix())

		#Candidate genes (TAIR IDs)
#		cg_tair_ids = ['AT1G04400', 'AT1G65480', 'AT1G77080', 'AT2G26330', 'AT4G00650', 'AT5G10140', 'AT5G65050', \
#				'AT5G65060', 'AT5G65070', 'AT5G65080'] \
#				#FT, MAF, ER, FRI, FLC, MAF2-MAF5 (Salome et al. 2011)
		cg_tair_ids = [ 'AT4G00650', 'AT5G10140'] #FRI and FLC
		cgs = gr.get_genes_w_tair_id(cg_tair_ids)
		snp_priors = sd.get_cand_genes_snp_priors(cand_genes=cgs, radius=10000, cg_prior_fold_incr=100.0)
		#pdb.set_trace()
		#snp_priors = sd.get_snp_priors(cpp_list)#, cand_genes=cgs)
		lm.emmax_step_wise(phen_vals, K, sd, 10, env.env['tmp_dir'] + 'bf_flc_%s_%d' % (phen_name, pid),
				snp_priors=snp_priors, cand_gene_list=cgs, save_pvals=True, snp_choose_criteria='ppas',
				emma_num=1000)







if __name__ == '__main__':
	if len(sys.argv) > 2:
		_run_()
	else:
		#generate_example_figure_1()
		generate_example_figure4()
		generate_example_figure5()
		#generate_example_bayesian_figures_flc()
		plot_local_sodium_figure()
		#plot_all_human_figures()
		#generate_results_figure_2()
		#plot_first_second_power_fdr(10)
	    #benchmark_human_emmax()
#	for herit in [5, 10, 20, 25]:
#		for ws in [0, 1000, 5000, 10000, 25000]:
#			file_name = '/tmp/fig2_h%d_ws%d.png' % (herit, ws)
#			generate_results_figure_2(file_name=file_name, herit=herit, window_size=ws)
	#generate_example_figure_1()
#	generate_results_figure_2(png_file_name='/tmp/mlt_two_loci_sim_h10_d1000.png', pdf_file_name='/tmp/mlt_two_loci_sim_h10_d1000.pdf', herit=10, window_size=1000)
#	generate_results_figure_2(png_file_name='/tmp/mlt_two_loci_sim_h10_d25000.png', pdf_file_name='/tmp/mlt_two_loci_sim_h10_d25000.pdf', herit=10, window_size=25000)
#	generate_results_figure_2(png_file_name='/tmp/mlt_two_loci_sim_h25_d1000.png', pdf_file_name='/tmp/mlt_two_loci_sim_h25_d1000.pdf', herit=25, window_size=1000)
#	generate_results_figure_2(png_file_name='/tmp/mlt_two_loci_sim_h25_d25000.png', pdf_file_name='/tmp/mlt_two_loci_sim_h25_d25000.pdf', herit=25, window_size=25000)
#	generate_results_figure_3(png_file_name='/tmp/100_loci_sim.png', pdf_file_name='/tmp/100_loci_sim.pdf')
	#perform_human_emmax(1)
	#plot_all_human_figures()
	#generate_example_bayesian_figures()
#	sd = dp.load_250K_snps()
#	simulate_phenotypes(env.env['tmp_dir'] + 'simulated_phenotypes.pickled', sd)
	#perform_human_emmax(4)
	#test_bayes_factor_enrichment_flc()
	#generate_example_bayesian_figures_flc()
	print "Done!!\n"
