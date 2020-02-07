"""
Simulations for joint analysis of correlated traits.
"""
import scipy as sp
from scipy import stats
import random
import pdb
import sys
import dataParsers as dp
import phenotypeData as pd
import analyze_gwas_results as agr
import env
import os
import cPickle
import linear_models as lm
import math
import pylab
import time
import gwaResults as gr
import bisect


max_num = 200

#Simulate phenotypes..
def simulate_traits(sd, num_traits=max_num, heritabilities=['0.01', '0.1', '0.25', '0.5', '0.75', '0.9', '0.99'],
			num_effects=100):
	"""
	Return a dictionary object of traits
	
	100 causal SNPs, with exponentially distributed effect sizes.
	"""

	snps = sd.getSnps()
	cpl = sd.getChrPosList()
	num_snps = len(snps)
	num_accessions = len(sd.accessions)

	sim_traits_dict = {}
	for her in heritabilities:
		h = float(her)
		trait_pairs = []
		sample_indices_list = []
		snp_effects_list = []
		num_non_overlapping_list = []
		rho_est_list = []
		trait1_perc_var_list = []
		trait2_perc_var_list = []
		for i in range(num_traits):
			if (i + 1) % int(num_traits / 100) == 0:
				sys.stdout.write('.')
				sys.stdout.flush()
				#print (i + 1) / int(num_traits / 100)
			#Simulate trait pair..

			num_non_overlapping = int(round(stats.beta.rvs(4, 2) * num_effects))
			num_non_overlapping_list.append(num_non_overlapping)
			num_causatives = num_effects + num_non_overlapping
			sample_indices = random.sample(range(num_snps), num_causatives)
			chosen_snps = sp.array([snps[i] for i in sample_indices])
			c = sp.random.random_integers(0, 1, (num_causatives, 1))
			chosen_snps = sp.absolute(c - chosen_snps)
			exp_effects = stats.expon.rvs(scale=1, size=(num_causatives, 1))
			#exp_effects = stats.norm.rvs(scale=1, size=(num_causatives, 1))
			snp_effects = chosen_snps * exp_effects
			snp_effects1 = snp_effects[:num_effects]
			snp_effects2 = snp_effects[-num_effects:]
			trait1 = sp.sum(snp_effects1, 0)
			trait2 = sp.sum(snp_effects2, 0)

			gv = sp.var(trait1, ddof=1)
			error = stats.norm.rvs(0, 1, size=num_accessions)
			ev = sp.var(error, ddof=1)
			n_trait1 = trait1 + error * sp.sqrt(((1.0 - h) / h) * (gv / ev))
			trait1_perc_var_list.append(sp.var(snp_effects1, 1) / sp.var(n_trait1))
			n_trait1 = (n_trait1 - sp.mean(n_trait1)) / sp.std(n_trait1)
			gv = sp.var(trait2, ddof=1)
			error = stats.norm.rvs(0, 1, size=num_accessions)
			ev = sp.var(error, ddof=1)
			n_trait2 = trait2 + error * sp.sqrt(((1.0 - h) / h) * (gv / ev))
			trait2_perc_var_list.append(sp.var(snp_effects2, 1) / sp.var(n_trait2))
			n_trait2 = (n_trait2 - sp.mean(n_trait2)) / sp.std(n_trait2)
			trait_pairs.append((n_trait1, n_trait2))
			sample_indices_list.append(sample_indices)
			snp_effects_list.append(exp_effects)
			rho_est = sp.corrcoef(trait1, trait2)
			rho_est_list.append(rho_est)

			#Variance contributions.


		sim_traits_dict[her] = {'trait_pairs':trait_pairs, 'sample_indices_list':sample_indices_list,
				'num_non_overlapping_list':num_non_overlapping_list, 'snp_effects_list':snp_effects_list,
				'rho_est_list':rho_est_list, 'trait1_perc_var_list':trait1_perc_var_list,
				'trait2_perc_var_list':trait2_perc_var_list}
	return sim_traits_dict



def _load_sim_data_(use_pickle=True):
	sim_data_file = env.env['data_dir'] + 'corr_trait_sim_data.pickled'
	if use_pickle and os.path.isfile(sim_data_file):
		with open(sim_data_file) as f:
			sim_traits_dict = cPickle.load(f)
	else:
		import dataParsers as dp
		sd = dp.load_250K_snps()
		sim_traits_dict = simulate_traits(sd)
		if use_pickle:
			with open(sim_data_file, 'wb') as f:
				cPickle.dump(sim_traits_dict, f)
	return sim_traits_dict




def run_joint_analysis(start_i, stop_i, heritability, debug_filter=1.0, run_id='joint_trait'):

	file_prefix = env.env['results_dir'] + '%s_corr_sim_h%0.2f_si%d_%d' % (run_id, float(heritability), start_i, stop_i)
	pickle_file = file_prefix + '_res_dict.pickled'
	if os.path.isfile(pickle_file):
		with open(pickle_file) as f:
			res_dict = cPickle.load(f)
	else:
		sim_traits_dict = _load_sim_data_()[heritability]

		if debug_filter < 1.0:
			sd = dp.load_250K_snps(debug_filter=1.0)
			cpl = sd.getChrPosList() #For simulated SNPs lookup.
			sd = dp.load_250K_snps(debug_filter=debug_filter)
		else:
			sd = dp.load_250K_snps(debug_filter=1.0)
			cpl = sd.getChrPosList() #For simulated SNPs lookup.

		snps = sd.getSnps()
	 	f_cpl = sd.getChrPosList()  #Filtered chrom, pos list.
	 	positions = sd.getPositions()
		chromosomes = sd.get_chr_list()
		K = lm.load_kinship_from_file(env.env['data_dir'] + 'kinship_matrix_cm72.pickled', sd.accessions)
		num_ecotypes = len(sd.accessions)

		res_dict = {'ks_stats_her':{'h01':[], 'h02':[], 'h12':[]},
				'ks_stats_pher':{'h01':[], 'h02':[], 'h12':[]},
				'ks_stats_indep':{'h01':[], 'h02':[], 'h12':[]},
				'ks_stats_marginal':[],
				'pval_medians_her':{'h01':[], 'h02':[], 'h12':[]},
				'pval_medians_pher':{'h01':[], 'h02':[], 'h12':[]},
				'pval_medians_indep':{'h01':[], 'h02':[], 'h12':[]},
				'pval_medians_marginal':[],
				'rho':[], 'rho_her':[], 'rho_pher':[],
				'p_hers':[], 'vgs':[], 'trait_corr':[],
				'causal_snps_pvals_pher':{'h01':[], 'h02':[], 'h12':[]},
				'causal_snps_pvals_her':{'h01':[], 'h02':[], 'h12':[]},
				'causal_snps_pvals_indep':{'h01':[], 'h02':[], 'h12':[]},
				'causal_snps_pvals_marginal':{'trait1':[], 'trait2':[]},
				'res_distances_pher':{'h01':[], 'h02':[], 'h12':[]},
				'res_distances_her':{'h01':[], 'h02':[], 'h12':[]},
				'res_distances_indep':{'h01':[], 'h02':[], 'h12':[]},
				'res_distances_marginal':{'trait1':[], 'trait2':[]},
				'num_non_overlapping': [], 'snp_effects':[],
				'sample_indices':[], 'trait_pair':[]}

		#Add some summary statistics on the results.
		for i in range(start_i, stop_i):
			print "Working on the %d't trait:" % i
			trait_pair = sim_traits_dict['trait_pairs'][i]
			res_dict['trait_pair'].append(trait_pair)
			sample_indices = sim_traits_dict['sample_indices_list'][i]
			res_dict['sample_indices'].append(sample_indices)
			num_non_overlapping = sim_traits_dict['num_non_overlapping_list'][i]
			res_dict['num_non_overlapping'].append(num_non_overlapping)
			snp_effects = sim_traits_dict['snp_effects_list'][i]
			res_dict['snp_effects'].append(snp_effects)
			rho = sim_traits_dict['rho_est_list'][i][0, 1]


			joint_phen_vals = trait_pair[0].tolist() + trait_pair[1].tolist()
			joint_ecotypes = sd.accessions + sd.accessions

			trait1_res = lm.emmax(snps, trait_pair[0], K)
			trait2_res = lm.emmax(snps, trait_pair[1], K)
			f_prefix = file_prefix + '_i%d' % i
			qq_file_name = f_prefix + '_qq_plot.png'
			log_qq_file_name = f_prefix + '_log_qq_plot.png'

			marginal_res1 = gr.Result(scores=trait1_res['ps'].tolist(), positions=positions, chromosomes=chromosomes)
			marginal_res2 = gr.Result(scores=trait2_res['ps'].tolist(), positions=positions, chromosomes=chromosomes)
			(areas, medians) = agr.qq_plot({'EMMAX_trait1':marginal_res1, 'EMMAX_trait2':marginal_res2},
					1000, method_types=['emma', 'emma'], mapping_labels=['EMMAX_trait1', 'EMMAX_trait2'],
					phenName='marginal EMMAX', pngFile=qq_file_name)
			(ds, areas, slopes) = agr.log_qq_plot({'EMMAX_trait1':marginal_res1, 'EMMAX_trait2':marginal_res2},
					1000, 7, method_types=['emma', 'emma'], mapping_labels=['EMMAX_trait1', 'EMMAX_trait2'],
					phenName='marginal EMMAX', pngFile=log_qq_file_name)
			res_dict['ks_stats_marginal'].append(ds)
			res_dict['pval_medians_marginal'].append(medians)


			gen_var_list = [trait1_res['vg'], trait2_res['vg']]
			err_var_list = [trait1_res['ve'], trait2_res['ve']]
			her_list = [trait1_res['pseudo_heritability'], trait2_res['pseudo_heritability']]
			res_dict['p_hers'].append(her_list)
			res_dict['vgs'].append(gen_var_list)

			min_ps = sp.minimum(trait1_res['ps'], trait2_res['ps'])




			#FINISH EMMA stuff..  from here on.. for comparison..

			#E in the model.
			joint_X = [[1] * (2 * num_ecotypes), [0] * num_ecotypes + [1] * num_ecotypes]
			E = sp.array([0] * num_ecotypes + [1] * num_ecotypes)

			#Get correlations between traits..
			corr_mat = sp.corrcoef(trait_pair[0], trait_pair[1])
			print her_list, gen_var_list, err_var_list, corr_mat[0, 1]
			res_dict['trait_corr'].append(corr_mat[0, 1])



			#Construct the full variance matrix
			V = sp.zeros((2 * num_ecotypes, 2 * num_ecotypes))
			V_2 = sp.zeros((2 * num_ecotypes, 2 * num_ecotypes))
			V_3 = sp.zeros((2 * num_ecotypes, 2 * num_ecotypes))
			V[0:num_ecotypes, 0:num_ecotypes] = gen_var_list[0] * K + err_var_list[0] * sp.eye(num_ecotypes)
			V[num_ecotypes:2 * num_ecotypes, num_ecotypes:2 * num_ecotypes] = \
					gen_var_list[1] * K + err_var_list[1] * sp.eye(num_ecotypes)
			V_2[0:num_ecotypes, 0:num_ecotypes] = V[0:num_ecotypes, 0:num_ecotypes]
			V_2[num_ecotypes:2 * num_ecotypes, num_ecotypes:2 * num_ecotypes] = \
					V[num_ecotypes:2 * num_ecotypes, num_ecotypes:2 * num_ecotypes]
			V_3[0:num_ecotypes, 0:num_ecotypes] = V[0:num_ecotypes, 0:num_ecotypes]
			V_3[num_ecotypes:2 * num_ecotypes, num_ecotypes:2 * num_ecotypes] = \
					V[num_ecotypes:2 * num_ecotypes, num_ecotypes:2 * num_ecotypes]

			rho_est = corr_mat[0, 1] / math.sqrt(her_list[0] * her_list[1])
			if rho_est > 1: rho_est = 1.0
			if rho_est < -1: rho_est = -1.0
			rho_est_2 = corr_mat[0, 1] / float(heritability)
			if rho_est_2 > 1: rho_est_2 = 1.0
			if rho_est_2 < -1: rho_est_2 = -1.0

			print rho, rho_est, rho_est_2
			res_dict['rho'].append(rho)
			res_dict['rho_pher'].append(rho_est)
			res_dict['rho_her'].append(rho_est_2)
			V[0:num_ecotypes, num_ecotypes:2 * num_ecotypes] = \
					rho_est * math.sqrt(gen_var_list[0] * gen_var_list[1]) * K
			V[num_ecotypes:2 * num_ecotypes, 0:num_ecotypes] = V[0:num_ecotypes, num_ecotypes:2 * num_ecotypes]
			V_2[0:num_ecotypes, num_ecotypes:2 * num_ecotypes] = \
					rho_est_2 * math.sqrt(gen_var_list[0] * gen_var_list[1]) * K
			V_2[num_ecotypes:2 * num_ecotypes, 0:num_ecotypes] = V_2[0:num_ecotypes, num_ecotypes:2 * num_ecotypes]

			print 'Performing Cholesky decompositions'
			H_sqrt = lm.cholesky(V)
			H_sqrt_inv = sp.mat(H_sqrt.T).I
			H_sqrt_2 = lm.cholesky(V_2)
			H_sqrt_inv_2 = sp.mat(H_sqrt_2.T).I

			H_sqrt_3 = lm.cholesky(V_3) #Independence..
			H_sqrt_inv_3 = sp.mat(H_sqrt_3.T).I

			#pdb.set_trace()

			#Set up analysis..
			lmm = lm.LinearMixedModel(joint_phen_vals)
			lmm.set_factors(joint_X, False)

			#Doubling the SNPs!!!
			Z = sp.int16(sp.mat(joint_ecotypes).T == sp.mat(sd.accessions))

			#Running EMMAX(s)
			t0 = time.time()
			print 'Running EMMAX full model'
			res1 = lmm.emmax_full_model_gxe(snps, E, H_sqrt_inv, Z)
			t = time.time() - t0
			print 'EMMAX took %d minutes and %0.2f seconds.' % (int(t / 60), t % 60)
			t0 = time.time()
			print 'Running EMMAX full model'
			res2 = lmm.emmax_full_model_gxe(snps, E, H_sqrt_inv_2, Z)
			t = time.time() - t0
			print 'Second EMMAX took %d minutes and %0.2f seconds.' % (int(t / 60), t % 60)
			t0 = time.time()
			print 'Running EMMAX full model'
			res3 = lmm.emmax_full_model_gxe(snps, E, H_sqrt_inv_3, Z)
			t = time.time() - t0
			print 'Indep EMMAX took %d minutes and %0.2f seconds.' % (int(t / 60), t % 60)

			res1_h01 = gr.Result(scores=res1['ps_h01'].tolist(), positions=positions, chromosomes=chromosomes)
			res1_h02 = gr.Result(scores=res1['ps_h02'].tolist(), positions=positions, chromosomes=chromosomes)
			res1_h12 = gr.Result(scores=res1['ps_h12'].tolist(), positions=positions, chromosomes=chromosomes)

			res2_h01 = gr.Result(scores=res2['ps_h01'].tolist(), positions=positions, chromosomes=chromosomes)
			res2_h02 = gr.Result(scores=res2['ps_h02'].tolist(), positions=positions, chromosomes=chromosomes)
			res2_h12 = gr.Result(scores=res2['ps_h12'].tolist(), positions=positions, chromosomes=chromosomes)

			res3_h01 = gr.Result(scores=res3['ps_h01'].tolist(), positions=positions, chromosomes=chromosomes)
			res3_h02 = gr.Result(scores=res3['ps_h02'].tolist(), positions=positions, chromosomes=chromosomes)
			res3_h12 = gr.Result(scores=res3['ps_h12'].tolist(), positions=positions, chromosomes=chromosomes)

			qq_file_name = f_prefix + '_qq_plot_h01.png'
			log_qq_file_name = f_prefix + '_log_qq_plot_h01.png'
			(areas, medians_h01) = agr.qq_plot({'EMMAX_joint_her':res2_h01, 'EMMAX_joint_pher':res1_h01, 'EMMAX_joint_indep':res3_h01},
					1000, method_types=['emma', 'emma', 'emma'], mapping_labels=['EMMAX_joint_her', 'EMMAX_joint_pher', 'EMMAX_joint_indep'],
					phenName='G', pngFile=qq_file_name)
			(ds_h01, areas, slopes) = agr.log_qq_plot({'EMMAX_joint_her':res2_h01, 'EMMAX_joint_pher':res1_h01, 'EMMAX_joint_indep':res3_h01},
					1000, 7, method_types=['emma', 'emma', 'emma'], mapping_labels=['EMMAX_joint_her', 'EMMAX_joint_pher', 'EMMAX_joint_indep'],
					phenName='G', pngFile=log_qq_file_name)
			qq_file_name = f_prefix + '_qq_plot_h02.png'
			log_qq_file_name = f_prefix + '_log_qq_plot_h02.png'
			(areas, medians_h02) = agr.qq_plot({'EMMAX_joint_her':res2_h02, 'EMMAX_joint_pher':res1_h02, 'EMMAX_joint_indep':res3_h02},
					1000, method_types=['emma', 'emma', 'emma'], mapping_labels=['EMMAX_joint_her', 'EMMAX_joint_pher', 'EMMAX_joint_indep'],
					phenName='G+GxE', pngFile=qq_file_name)
			(ds_h02, areas, slopes) = agr.log_qq_plot({'EMMAX_joint_her':res2_h02, 'EMMAX_joint_pher':res1_h02, 'EMMAX_joint_indep':res3_h02},
					1000, 7, method_types=['emma', 'emma', 'emma'], mapping_labels=['EMMAX_joint_her', 'EMMAX_joint_pher', 'EMMAX_joint_indep'],
					phenName='G+GxE', pngFile=log_qq_file_name)
			qq_file_name = f_prefix + '_qq_plot_h12.png'
			log_qq_file_name = f_prefix + '_log_qq_plot_h12.png'
			(areas, medians_h12) = agr.qq_plot({'EMMAX_joint_her':res2_h12, 'EMMAX_joint_pher':res1_h12, 'EMMAX_joint_indep':res3_h12},
					1000, method_types=['emma', 'emma', 'emma'], mapping_labels=['EMMAX_joint_her', 'EMMAX_joint_pher', 'EMMAX_joint_indep'],
					phenName='GxE', pngFile=qq_file_name)
			(ds_h12, areas, slopes) = agr.log_qq_plot({'EMMAX_joint_her':res2_h12, 'EMMAX_joint_pher':res1_h12, 'EMMAX_joint_indep':res3_h12},
					1000, 7, method_types=['emma', 'emma', 'emma'], mapping_labels=['EMMAX_joint_her', 'EMMAX_joint_pher', 'EMMAX_joint_indep'],
					phenName='GxE', pngFile=log_qq_file_name)


			res_dict['ks_stats_her']['h01'].append(ds_h01[0])
			res_dict['ks_stats_her']['h02'].append(ds_h02[0])
			res_dict['ks_stats_her']['h12'].append(ds_h12[0])
			res_dict['ks_stats_pher']['h01'].append(ds_h01[1])
			res_dict['ks_stats_pher']['h02'].append(ds_h02[1])
			res_dict['ks_stats_pher']['h12'].append(ds_h12[1])
			res_dict['ks_stats_indep']['h01'].append(ds_h01[2])
			res_dict['ks_stats_indep']['h02'].append(ds_h02[2])
			res_dict['ks_stats_indep']['h12'].append(ds_h12[2])
			res_dict['pval_medians_her']['h01'].append(medians_h01[0])
			res_dict['pval_medians_her']['h02'].append(medians_h02[0])
			res_dict['pval_medians_her']['h12'].append(medians_h12[0])
			res_dict['pval_medians_pher']['h01'].append(medians_h01[1])
			res_dict['pval_medians_pher']['h02'].append(medians_h02[1])
			res_dict['pval_medians_pher']['h12'].append(medians_h12[1])
			res_dict['pval_medians_indep']['h01'].append(medians_h01[2])
			res_dict['pval_medians_indep']['h02'].append(medians_h02[2])
			res_dict['pval_medians_indep']['h12'].append(medians_h12[2])

			#Process results.
			#Do we detect the causal loci?
			sample_cpl = []
			for j in sample_indices:
				sample_cpl.append(cpl[j])
			for t in sample_cpl:
				pi = bisect.bisect(f_cpl, t) - 1
				if f_cpl[pi] == t:
					d = res_dict['causal_snps_pvals_pher']
					d['h01'].append(res1_h01.snp_results['scores'][pi])
					d['h02'].append(res1_h02.snp_results['scores'][pi])
					d['h12'].append(res1_h12.snp_results['scores'][pi])
					d = res_dict['causal_snps_pvals_her']
					d['h01'].append(res2_h01.snp_results['scores'][pi])
					d['h02'].append(res2_h02.snp_results['scores'][pi])
					d['h12'].append(res2_h12.snp_results['scores'][pi])
					d = res_dict['causal_snps_pvals_indep']
					d['h01'].append(res3_h01.snp_results['scores'][pi])
					d['h02'].append(res3_h02.snp_results['scores'][pi])
					d['h12'].append(res3_h12.snp_results['scores'][pi])
					d = res_dict['causal_snps_pvals_marginal']
					d['trait1'].append(marginal_res1.snp_results['scores'][pi])
					d['trait2'].append(marginal_res2.snp_results['scores'][pi])
				else:
					for h1 in ['causal_snps_pvals_pher', 'causal_snps_pvals_her', 'causal_snps_pvals_indep']:
						d = res_dict[h1]
						for h2 in ['h01', 'h02', 'h12']:
							d[h2].append(None)
					d = res_dict['causal_snps_pvals_marginal']
					d['trait1'].append(None)
					d['trait2'].append(None)


			trait1_only_causal = sample_cpl[:100]
			trait2_only_causal = sample_cpl[-100:]
			trait1_perc_var = sim_traits_dict['trait1_perc_var_list'][i]
			trait2_perc_var = sim_traits_dict['trait2_perc_var_list'][i]
			print 'Trait 1 specific'
			for j in range(num_non_overlapping):
				pv = trait1_perc_var[j]
				pi = j
				if pv > 0.01 and res_dict['causal_snps_pvals_pher']['h01'][pi]:
					print pv, trait1_only_causal[j]
					print res_dict['causal_snps_pvals_marginal']['trait1'][pi], res_dict['causal_snps_pvals_marginal']['trait2'][pi]
					for h in ['h01', 'h02', 'h12']:
						print h, res_dict['causal_snps_pvals_pher'][h][pi], res_dict['causal_snps_pvals_her'][h][pi], res_dict['causal_snps_pvals_indep'][h][pi]
			print 'Trait 2 specific'
			for j in range(num_non_overlapping):
				pv = trait2_perc_var[(100 - num_non_overlapping) + j]
				pi = 100 + j
				if pv > 0.01 and res_dict['causal_snps_pvals_pher']['h01'][pi]:
					print pv, trait2_only_causal[(100 - num_non_overlapping) + j]
					print res_dict['causal_snps_pvals_marginal']['trait1'][pi], res_dict['causal_snps_pvals_marginal']['trait2'][pi]
					for h in ['h01', 'h02', 'h12']:
						print h, res_dict['causal_snps_pvals_pher'][h][pi], res_dict['causal_snps_pvals_her'][h][pi], res_dict['causal_snps_pvals_indep'][h][pi]
			print 'Common SNPs'
			for j in range(num_non_overlapping, 100):
				pv1 = trait1_perc_var[j]
				pv2 = trait2_perc_var[j - num_non_overlapping]
				if (pv1 > 0.01 or pv2 > 0.01) and res_dict['causal_snps_pvals_pher']['h01'][j]:
					print pv1, pv2, sample_cpl[j]
					print res_dict['causal_snps_pvals_marginal']['trait1'][j], res_dict['causal_snps_pvals_marginal']['trait2'][j]
					for h in ['h01', 'h02', 'h12']:
						print h, res_dict['causal_snps_pvals_pher'][h][j], res_dict['causal_snps_pvals_her'][h][j], res_dict['causal_snps_pvals_indep'][h][j]



			#Generate manhattan plots..
			png_file_name = f_prefix + '_manhattan_res1.png'
			marginal_res1.neg_log_trans()
			marginal_res1.plot_manhattan(png_file=png_file_name, plot_bonferroni=True)
			png_file_name = f_prefix + '_manhattan_res2.png'
			marginal_res2.neg_log_trans()
			marginal_res2.plot_manhattan(png_file=png_file_name, plot_bonferroni=True)


			f_prefix = file_prefix + '_pher_i%d' % i
			png_file_name = f_prefix + '_manhattan_h01.png'
			res1_h01.neg_log_trans()
			res1_h01.plot_manhattan(png_file=png_file_name, plot_bonferroni=True)
			png_file_name = f_prefix + '_manhattan_h02.png'
			res1_h02.neg_log_trans()
			res1_h02.plot_manhattan(png_file=png_file_name, plot_bonferroni=True)
			png_file_name = f_prefix + '_manhattan_h12.png'
			res1_h12.neg_log_trans()
			res1_h12.plot_manhattan(png_file=png_file_name, plot_bonferroni=True)

			f_prefix_2 = file_prefix + '_her_i%d' % i
			png_file_name_2 = f_prefix_2 + '_manhattan_h01.png'
			res2_h01.neg_log_trans()
			res2_h01.plot_manhattan(png_file=png_file_name_2, plot_bonferroni=True)
			png_file_name_2 = f_prefix_2 + '_manhattan_h02.png'
			res2_h02.neg_log_trans()
			res2_h02.plot_manhattan(png_file=png_file_name_2, plot_bonferroni=True)
			png_file_name_2 = f_prefix_2 + '_manhattan_h12.png'
			res2_h12.neg_log_trans()
			res2_h12.plot_manhattan(png_file=png_file_name_2, plot_bonferroni=True)

			f_prefix_3 = file_prefix + '_indep_i%d' % i
			png_file_name_3 = f_prefix_3 + '_manhattan_h01.png'
			res3_h01.neg_log_trans()
			res3_h01.plot_manhattan(png_file=png_file_name_3, plot_bonferroni=True)
			png_file_name_3 = f_prefix_3 + '_manhattan_h02.png'
			res3_h02.neg_log_trans()
			res3_h02.plot_manhattan(png_file=png_file_name_3, plot_bonferroni=True)
			png_file_name_3 = f_prefix_3 + '_manhattan_h12.png'
			res3_h12.neg_log_trans()
			res3_h12.plot_manhattan(png_file=png_file_name_3, plot_bonferroni=True)


			#Where are the peaks?
			#Filter at Bonferroni thresholds..
			print 'Filtering results, leaving only SNPs above Bonferroni threshold.'
			bonf_threshold = -math.log10(1.0 / (20 * len(snps)))

			dl = marginal_res1.get_min_distances(trait1_only_causal) \
				if marginal_res1.filter_attr('scores', bonf_threshold) else None
			res_dict['res_distances_marginal']['trait1'].append(dl)
			dl = marginal_res2.get_min_distances(trait2_only_causal) \
				if marginal_res2.filter_attr('scores', bonf_threshold) else None
			res_dict['res_distances_marginal']['trait2'].append(dl)
			dl = res1_h01.get_min_distances(sample_cpl) if res1_h01.filter_attr('scores', bonf_threshold) else None
			res_dict['res_distances_pher']['h01'].append(dl)
			dl = res1_h02.get_min_distances(sample_cpl) if res1_h02.filter_attr('scores', bonf_threshold) else None
			res_dict['res_distances_pher']['h02'].append(dl)
			dl = res1_h12.get_min_distances(sample_cpl) if res1_h12.filter_attr('scores', bonf_threshold) else None
			res_dict['res_distances_pher']['h12'].append(dl)
			dl = res2_h01.get_min_distances(sample_cpl) if res2_h01.filter_attr('scores', bonf_threshold) else None
			res_dict['res_distances_her']['h01'].append(dl)
			dl = res2_h02.get_min_distances(sample_cpl) if res2_h02.filter_attr('scores', bonf_threshold) else None
			res_dict['res_distances_her']['h02'].append(dl)
			dl = res2_h12.get_min_distances(sample_cpl) if res2_h12.filter_attr('scores', bonf_threshold) else None
			res_dict['res_distances_her']['h12'].append(dl)
			dl = res3_h01.get_min_distances(sample_cpl) if res3_h01.filter_attr('scores', bonf_threshold) else None
			res_dict['res_distances_indep']['h01'].append(dl)
			dl = res3_h02.get_min_distances(sample_cpl) if res3_h02.filter_attr('scores', bonf_threshold) else None
			res_dict['res_distances_indep']['h02'].append(dl)
			dl = res3_h12.get_min_distances(sample_cpl) if res3_h12.filter_attr('scores', bonf_threshold) else None
			res_dict['res_distances_indep']['h12'].append(dl)



	#Pickle these results.
	with open(pickle_file, 'wb') as f:
		cPickle.dump(res_dict, f)
	return res_dict


def plot_results(run_id='joint_trait'):

	rd_map = {}
	kl = ['num_non_overlapping', 'snp_effects', 'sample_indices', 'trait_pair', 'p_hers', 'vgs', 'trait_corr',
		'rho', 'rho_her', 'rho_pher', 'ks_stats_marginal', 'pval_medians_marginal']
	kl2 = ['ks_stats_her', 'ks_stats_pher', 'ks_stats_indep', 'pval_medians_her', 'pval_medians_pher',
		'pval_medians_indep', 'causal_snps_pvals_pher', 'causal_snps_pvals_her', 'causal_snps_pvals_indep',
		'res_distances_pher', 'res_distances_her', 'res_distances_indep']
	kl3 = ['causal_snps_pvals_marginal', 'res_distances_marginal']
	sd = dp.load_250K_snps(debug_filter=1.0)
	snps = sd.getSnps()
	neg_log_bonf_threshold = -math.log10(1.0 / (20.0 * len(snps)))
	bonf_threshold = 1.0 / (20.0 * len(snps))
	print bonf_threshold, neg_log_bonf_threshold
	disjoint_pher_h01_pvals = {'0.25':[], '0.50':[], '0.75':[]}
	disjoint_her_h01_pvals = {'0.25':[], '0.50':[], '0.75':[]}
	disjoint_indep_h01_pvals = {'0.25':[], '0.50':[], '0.75':[]}
	disjoint_pher_h02_pvals = {'0.25':[], '0.50':[], '0.75':[]}
	disjoint_her_h02_pvals = {'0.25':[], '0.50':[], '0.75':[]}
	disjoint_indep_h02_pvals = {'0.25':[], '0.50':[], '0.75':[]}
	disjoint_pher_h12_pvals = {'0.25':[], '0.50':[], '0.75':[]}
	disjoint_her_h12_pvals = {'0.25':[], '0.50':[], '0.75':[]}
	disjoint_indep_h12_pvals = {'0.25':[], '0.50':[], '0.75':[]}
	disjoint_min_marg_pvals = {'0.25':[], '0.50':[], '0.75':[]}
	common_pher_h01_pvals = {'0.25':[], '0.50':[], '0.75':[]}
	common_her_h01_pvals = {'0.25':[], '0.50':[], '0.75':[]}
	common_indep_h01_pvals = {'0.25':[], '0.50':[], '0.75':[]}
	common_pher_h02_pvals = {'0.25':[], '0.50':[], '0.75':[]}
	common_her_h02_pvals = {'0.25':[], '0.50':[], '0.75':[]}
	common_indep_h02_pvals = {'0.25':[], '0.50':[], '0.75':[]}
	common_pher_h12_pvals = {'0.25':[], '0.50':[], '0.75':[]}
	common_her_h12_pvals = {'0.25':[], '0.50':[], '0.75':[]}
	common_indep_h12_pvals = {'0.25':[], '0.50':[], '0.75':[]}
	common_min_marginal_pvals = {'0.25':[], '0.50':[], '0.75':[]}
	disjoint_fraction_found = {'0.25':{'her':{'h01':0, 'h02':0, 'h12':0}, 'pher':{'h01':0, 'h02':0, 'h12':0},
					'indep':{'h01':0, 'h02':0, 'h12':0}, 'marginal':0},
				'0.50':{'her':{'h01':0, 'h02':0, 'h12':0}, 'pher':{'h01':0, 'h02':0, 'h12':0},
				'indep':{'h01':0, 'h02':0, 'h12':0}, 'marginal':0},
				'0.75':{'her':{'h01':0, 'h02':0, 'h12':0}, 'pher':{'h01':0, 'h02':0, 'h12':0},
				'indep':{'h01':0, 'h02':0, 'h12':0}, 'marginal':0}}
	common_fraction_found = {'0.25':{'her':{'h01':0, 'h02':0, 'h12':0}, 'pher':{'h01':0, 'h02':0, 'h12':0},
					'indep':{'h01':0, 'h02':0, 'h12':0}, 'marginal':0},
				'0.50':{'her':{'h01':0, 'h02':0, 'h12':0}, 'pher':{'h01':0, 'h02':0, 'h12':0},
				'indep':{'h01':0, 'h02':0, 'h12':0}, 'marginal':0},
				'0.75':{'her':{'h01':0, 'h02':0, 'h12':0}, 'pher':{'h01':0, 'h02':0, 'h12':0},
				'indep':{'h01':0, 'h02':0, 'h12':0}, 'marginal':0}}
	for (heritability, chunk_size) in [('0.25', 1), ('0.50', 2), ('0.75', 1)]:
		#Combine res_dicts..
		rd = {'ks_stats_her':{'h01':[], 'h02':[], 'h12':[]},
				'ks_stats_pher':{'h01':[], 'h02':[], 'h12':[]},
				'ks_stats_indep':{'h01':[], 'h02':[], 'h12':[]},
				'ks_stats_marginal':[],
				'pval_medians_her':{'h01':[], 'h02':[], 'h12':[]},
				'pval_medians_pher':{'h01':[], 'h02':[], 'h12':[]},
				'pval_medians_indep':{'h01':[], 'h02':[], 'h12':[]},
				'pval_medians_marginal':[],
				'rho':[], 'rho_her':[], 'rho_pher':[],
				'p_hers':[], 'vgs':[], 'trait_corr':[],
				'causal_snps_pvals_pher':{'h01':[], 'h02':[], 'h12':[]},
				'causal_snps_pvals_her':{'h01':[], 'h02':[], 'h12':[]},
				'causal_snps_pvals_indep':{'h01':[], 'h02':[], 'h12':[]},
				'causal_snps_pvals_marginal':{'trait1':[], 'trait2':[]},
				'res_distances_pher':{'h01':[], 'h02':[], 'h12':[]},
				'res_distances_her':{'h01':[], 'h02':[], 'h12':[]},
				'res_distances_indep':{'h01':[], 'h02':[], 'h12':[]},
				'res_distances_marginal':{'trait1':[], 'trait2':[]},
				'num_non_overlapping': [], 'snp_effects':[],
				'sample_indices':[], 'trait_pair':[]}
		for i in range(0, max_num, chunk_size):
			file_prefix = env.env['results_dir'] + '%s_corr_sim_h%s_si%d_%d' % (run_id, heritability, i, i + chunk_size)
			pickle_file = file_prefix + '_res_dict.pickled'
			if os.path.isfile(pickle_file):
				with open(pickle_file) as f:
					res_dict = cPickle.load(f)
				for k in kl:
					rd[k].extend(res_dict[k])
				for k in kl2:
					for k2 in ['h01', 'h02', 'h12']:
						rd[k][k2].extend(res_dict[k][k2])
				for k in kl3:
					for k2 in ['trait1', 'trait2']:
						rd[k][k2].extend(res_dict[k][k2])
		rd_map[heritability] = rd

		print 'Plotting.'
		#Basic plots
		pylab.figure(figsize=(5, 5))
		pylab.axes([0.15, 0.1, 0.83, 0.84])

		pylab.plot(rd['rho'], rd['rho_pher'], 'g.', label='pseudo-herit. est.', alpha=0.6)
		pylab.plot(rd['rho'], rd['rho_her'], 'b.', label='herit. est.', alpha=0.6)
		min_rho = min(rd['rho'])
		max_rho = 1.0#max(rd['rho'])
		range_rho = max_rho - min_rho
		max_y = 1.0
		min_y = min(min(rd['rho_pher']), min(rd['rho_her']))
		range_y = max_y - min_y
		pylab.xlabel(r'True $\rho$', fontsize=14)
		if heritability == '0.25':
			pylab.ylabel(r'$\rho$ estimated using herit.', fontsize=14)

		#pylab.rcParams['legend.fontsize'] = 'small'
		if heritability == '0.75':
			pylab.legend(loc=2, numpoints=1)
		pylab.axis([min_rho - 0.05 * range_rho, max_rho + 0.05 * range_rho,
				min_y - 0.05 * range_y, max_y + 0.05 * range_y])
		lim = [min(pylab.xlim()[0], pylab.ylim()[0]), max(pylab.xlim()[1], pylab.ylim()[1])]
		pylab.plot(lim, lim, alpha=0.6, color='k', ls='--')
		pylab.title(r'$h^2= %s$' % heritability)
		pylab.savefig(env.env['tmp_dir'] + 'rho_est_plot_h%s.png' % heritability)
		pylab.clf()

		rho_dist_pher = sp.absolute(sp.array(rd['rho']) - sp.array(rd['rho_pher']))
		rho_dist_her = sp.absolute(sp.array(rd['rho']) - sp.array(rd['rho_her']))
		min_rho = min(min(rho_dist_pher), min(rho_dist_her))
		max_rho = max(max(rho_dist_pher), max(rho_dist_her))
		range_rho = max_rho - min_rho
		pylab.axes([0.17, 0.1, 0.8, 0.87])
		for h in ['h01', 'h02', 'h12']:
			png_file_name = env.env['tmp_dir'] + 'rho_median_%s_pval_plot_h%s.png' % (h, heritability)
			pylab.plot(rho_dist_pher, rd['pval_medians_pher'][h], 'g.', label='pseudo-herit. est.', alpha=0.6)
			pylab.plot(rho_dist_her, rd['pval_medians_her'][h], 'b.', label='true herit. est.', alpha=0.6)
			max_y = max(map(max, [rd['pval_medians_pher'][h], rd['pval_medians_her'][h]]))
			min_y = min(map(min, [rd['pval_medians_pher'][h], rd['pval_medians_her'][h]]))
			range_y = max_y - min_y
			pylab.xlabel(r'$\rho - \rho_{\mathrm{est}}$', fontsize=14)
			pylab.ylabel('Median p-value deviation.', fontsize=14)
			#pylab.rcParams['legend.fontsize'] = 'small'
			pylab.legend(loc=2, numpoints=1)
			pylab.axis([min_rho - 0.05 * range_rho, max_rho + 0.05 * range_rho,
					min_y - 0.05 * range_y, max_y + 0.1 * range_y])
			pylab.plot(pylab.xlim(), [0, 0], alpha=0.6, color='k', ls='--')
			pylab.savefig(png_file_name)
			pylab.clf()

		for h in ['h01', 'h02', 'h12']:
			png_file_name = env.env['tmp_dir'] + 'rho_ks_stat_%s_pval_plot_h%s.png' % (h, heritability)
			pylab.plot(rho_dist_pher, rd['ks_stats_pher'][h], 'g.', label='pseudo-herit. est.', alpha=0.6)
			pylab.plot(rho_dist_her, rd['ks_stats_her'][h], 'b.', label='true herit. est.', alpha=0.6)
			max_y = max(map(max, [rd['ks_stats_pher'][h], rd['ks_stats_her'][h]]))
			min_y = min(map(min, [rd['ks_stats_pher'][h], rd['ks_stats_her'][h]]))
			range_y = max_y - min_y
			pylab.xlabel(r'$\rho - \rho_{\mathrm{est}}$')
			pylab.ylabel('Kolmogorov-Smirnov statistic.')
			#pylab.rcParams['legend.fontsize'] = 'small'
			pylab.legend(loc=2, numpoints=1)
			pylab.axis([min_rho - 0.05 * range_rho, max_rho + 0.05 * range_rho,
					min_y - 0.05 * range_y, max_y + 0.1 * range_y])
			pylab.savefig(png_file_name)
			pylab.clf()

		trait1_perc_var = []
		trait2_perc_var = []
		shift = 0
		tot_overlapping = 0
		tot_non_overlapping = 0
		for i in range(len(rd['num_non_overlapping'])):
			cpl = sd.getChrPosList() #For simulated SNPs lookup.
			sample_indices = rd['sample_indices'][i]
			chosen_snps = sp.array([snps[j] for j in sample_indices])
			snp_effects = rd['snp_effects'][i]
			snp_effects = chosen_snps * snp_effects
			num_non_overlapping = rd['num_non_overlapping'][i]
			trait1 = rd['trait_pair'][i][0]
			trait2 = rd['trait_pair'][i][1]
			trait1_pv = []
			for snp in chosen_snps[:100]:
				trait1_pv.append(sp.corrcoef(snp, trait1)[0, 1] ** 2)
			trait1_perc_var.append(trait1_pv)
			trait2_pv = []
			for snp in chosen_snps[-100:]:
				trait2_pv.append(sp.corrcoef(snp, trait2)[0, 1] ** 2)
			trait2_perc_var.append(trait2_pv)

			for j in range(num_non_overlapping) + range(100, 100 + num_non_overlapping):
				pval1 = rd['causal_snps_pvals_marginal']['trait1'][shift + j]
				pval2 = rd['causal_snps_pvals_marginal']['trait2'][shift + j]
				p = min(pval1, pval2)
				if p < bonf_threshold: disjoint_fraction_found[heritability]['marginal'] += 1
				disjoint_min_marg_pvals[heritability].append(p)
				p = rd['causal_snps_pvals_pher']['h01'][shift + j]
				if p < bonf_threshold: disjoint_fraction_found[heritability]['pher']['h01'] += 1
				disjoint_pher_h01_pvals[heritability].append(p)
				p = rd['causal_snps_pvals_pher']['h02'][shift + j]
				if p < bonf_threshold: disjoint_fraction_found[heritability]['pher']['h02'] += 1
				disjoint_pher_h02_pvals[heritability].append(p)
				p = rd['causal_snps_pvals_pher']['h12'][shift + j]
				if p < bonf_threshold: disjoint_fraction_found[heritability]['pher']['h12'] += 1
				disjoint_pher_h12_pvals[heritability].append(p)
				p = rd['causal_snps_pvals_her']['h01'][shift + j]
				if p < bonf_threshold: disjoint_fraction_found[heritability]['her']['h01'] += 1
				disjoint_her_h01_pvals[heritability].append(p)
				p = rd['causal_snps_pvals_her']['h02'][shift + j]
				if p < bonf_threshold: disjoint_fraction_found[heritability]['her']['h02'] += 1
				disjoint_her_h02_pvals[heritability].append(p)
				p = rd['causal_snps_pvals_her']['h12'][shift + j]
				if p < bonf_threshold: disjoint_fraction_found[heritability]['her']['h12'] += 1
				disjoint_her_h12_pvals[heritability].append(p)
				p = rd['causal_snps_pvals_indep']['h01'][shift + j]
				if p < bonf_threshold: disjoint_fraction_found[heritability]['indep']['h01'] += 1
				disjoint_indep_h01_pvals[heritability].append(p)
				p = rd['causal_snps_pvals_indep']['h02'][shift + j]
				if p < bonf_threshold: disjoint_fraction_found[heritability]['indep']['h02'] += 1
				disjoint_indep_h02_pvals[heritability].append(p)
				p = rd['causal_snps_pvals_indep']['h12'][shift + j]
				if p < bonf_threshold: disjoint_fraction_found[heritability]['indep']['h12'] += 1
				disjoint_indep_h12_pvals[heritability].append(p)

			tot_overlapping += 100 - num_non_overlapping
			tot_non_overlapping += num_non_overlapping
			for j in range(num_non_overlapping, 100):

				pval1 = rd['causal_snps_pvals_marginal']['trait1'][shift + j]
				pval2 = rd['causal_snps_pvals_marginal']['trait2'][shift + j]
				p = min(pval1, pval2)
				if p < bonf_threshold: common_fraction_found[heritability]['marginal'] += 1
				common_min_marginal_pvals[heritability].append(p)
				p = rd['causal_snps_pvals_her']['h01'][shift + j]
				if p < bonf_threshold: common_fraction_found[heritability]['her']['h01'] += 1
				common_her_h01_pvals[heritability].append(p)
				p = rd['causal_snps_pvals_her']['h02'][shift + j]
				if p < bonf_threshold: common_fraction_found[heritability]['her']['h02'] += 1
				common_her_h02_pvals[heritability].append(p)
				p = rd['causal_snps_pvals_her']['h12'][shift + j]
				if p < bonf_threshold: common_fraction_found[heritability]['her']['h12'] += 1
				common_her_h12_pvals[heritability].append(p)
				p = rd['causal_snps_pvals_pher']['h01'][shift + j]
				if p < bonf_threshold: common_fraction_found[heritability]['pher']['h01'] += 1
				common_pher_h01_pvals[heritability].append(p)
				p = rd['causal_snps_pvals_pher']['h02'][shift + j]
				if p < bonf_threshold: common_fraction_found[heritability]['pher']['h02'] += 1
				common_pher_h02_pvals[heritability].append(p)
				p = rd['causal_snps_pvals_pher']['h12'][shift + j]
				if p < bonf_threshold: common_fraction_found[heritability]['pher']['h12'] += 1
				common_pher_h12_pvals[heritability].append(p)
				p = rd['causal_snps_pvals_indep']['h01'][shift + j]
				if p < bonf_threshold: common_fraction_found[heritability]['indep']['h01'] += 1
				common_indep_h01_pvals[heritability].append(p)
				p = rd['causal_snps_pvals_indep']['h02'][shift + j]
				if p < bonf_threshold: common_fraction_found[heritability]['indep']['h02'] += 1
				common_indep_h02_pvals[heritability].append(p)
				p = rd['causal_snps_pvals_indep']['h12'][shift + j]
				if p < bonf_threshold: common_fraction_found[heritability]['indep']['h12'] += 1
				common_indep_h12_pvals[heritability].append(p)
			shift = shift + 100 + num_non_overlapping
		for k1 in ['her', 'pher', 'indep']:
			for k2 in ['h01', 'h02', 'h12']:
				common_fraction_found[heritability][k1][k2] = \
					common_fraction_found[heritability][k1][k2] / float(tot_overlapping)
				disjoint_fraction_found[heritability][k1][k2] = \
					disjoint_fraction_found[heritability][k1][k2] / float(tot_non_overlapping)
		common_fraction_found[heritability]['marginal'] = \
			common_fraction_found[heritability]['marginal'] / float(tot_overlapping)
		disjoint_fraction_found[heritability]['marginal'] = \
			disjoint_fraction_found[heritability]['marginal'] / float(tot_non_overlapping)


#		a1 = -sp.log10(rd['causal_snps_pvals_her']['h02'])
#		a2 = sp.minimum(sp.array(rd['causal_snps_pvals_marginal']['trait1']),
#				sp.array(rd['causal_snps_pvals_marginal']['trait2']))
#		a2 = -sp.log10(a2)
#		pylab.plot(a1, a2, '.', alpha=0.5)
#		max_x = a1.max()
#		max_y = a2.max()
#		max_val = max(max_x, max_y)
#		pylab.axis([-0.025 * max_x, 1.025 * max_x, -0.025 * max_y, 1.025 * max_y])
#		pylab.plot([-0.025 * max_val, 1.025 * max_val], [-0.025 * max_val, 1.025 * max_val],
#			ls=' - -', color='k', alpha=0.5)
#		pylab.plot([neg_log_bonf_threshold, neg_log_bonf_threshold], [-0.025 * max_val, 1.025 * max_val], color='m', ls='..')
#		pylab.plot([-0.025 * max_val, 1.025 * max_val], [neg_log_bonf_threshold, neg_log_bonf_threshold], color='m', ls='..')
#		pylab.savefig(env.env['tmp_dir'] + 'total_pvals_h % s.png' % heritability)

	#Common SNPs vs disjoint ones, power to detect them..
	pylab.clf()
	pylab.figure(figsize=(8, 6))
	l1 = disjoint_min_marg_pvals['0.25'] + disjoint_min_marg_pvals['0.50'] + disjoint_min_marg_pvals['0.75']
	a1 = -sp.log10(l1)
	l2 = disjoint_pher_h01_pvals['0.25'] + disjoint_pher_h01_pvals['0.50'] + disjoint_pher_h01_pvals['0.75']
	a2 = -sp.log10(l2)
	max_x = a1.max()
	max_y = a2.max()
	pylab.plot(a1, a2, '.', alpha=0.5, label='Disjoint causal SNPs')
	l1 = common_min_marginal_pvals['0.25'] + common_min_marginal_pvals['0.50'] + common_min_marginal_pvals['0.75']
	a1 = -sp.log10(l1)
	l2 = common_her_h01_pvals['0.25'] + common_her_h01_pvals['0.50'] + common_her_h01_pvals['0.75']
	a2 = -sp.log10(l2)
	pylab.plot(a1, a2, '.', alpha=0.5, label='Common causal SNPs')
	pylab.legend(loc=2)
	pylab.xlabel('Minimum marginal EMMAX p - value')
	pylab.ylabel('Joint EMMAX full model p - value')
	max_x = max(max_x, a1.max())
	max_y = max(max_y, a2.max())
	max_val = max(max_x, max_y)
	pylab.axis([-0.025 * max_x, 1.025 * max_x, -0.025 * max_y, 1.025 * max_y])
	pylab.plot([-0.025 * max_val, 1.025 * max_val], [-0.025 * max_val, 1.025 * max_val],
		ls=' - -', color='k', alpha=0.5)
	pylab.plot([neg_log_bonf_threshold, neg_log_bonf_threshold], [-0.025 * max_val, 1.025 * max_val], color='m', ls='-.')
	pylab.plot([-0.025 * max_val, 1.025 * max_val], [neg_log_bonf_threshold, neg_log_bonf_threshold], color='m', ls='-.')
	pylab.savefig(env.env['tmp_dir'] + 'common_pvals.png')

	print common_fraction_found
	print disjoint_fraction_found

	pylab.clf()


	#Common fraction barchart.
	bars_xs = [[], [], []]
	bars_ys = [[], [], []]
	vals = []
	bar_labels = ['Minimum marginal', 'Joint pseudo-herit.', 'Joint herit.']
	color_map = ['#1155BB', '#BB5511', '#11BB55']
	curr_x = 0
	pylab.figure(figsize=(6, 5))
	pylab.axes([0.16, 0.12, 0.82, 0.84])

	for h in ['0.25', '0.50', '0.75']:
		for j, m in enumerate(['marginal', 'pher', 'her']):
			if m != 'marginal':
				bars_ys[j].append(common_fraction_found[h][m]['h01'])
			else:
				bars_ys[j].append(common_fraction_found[h][m])
			bars_xs[j].append(curr_x)
			curr_x += 1
		curr_x += 1
	curr_x = curr_x - 1.2
	for i, lab in enumerate(bar_labels):
		vals.extend(bars_ys[i])
		pylab.bar(bars_xs[i], bars_ys[i], color=color_map[i], alpha=0.7, label=lab)
	tick_positions = [1.4, 5.4, 9.4]
	tick_labels = [r'$h^2=0.25$', r'$h^2=0.5$', r'$h^2=0.75$']
	pylab.xticks(tick_positions, tick_labels)#, rotation=45, ha='right')
	pylab.ylabel('Fraction of causative SNPs above Bonferroni')
	pylab.axis([-0.02 * curr_x, 1.02 * curr_x, 0.0, 1.15 * max(vals)])
#	pylab.text(0.8, 1.05 * max(bars_ys), r'$h^2=0.25$')
#	pylab.text(4.8, 1.05 * max(bars_ys), r'$h^2=0.5$')
#	pylab.text(8.8, 1.05 * max(bars_ys), r'$h^2=0.75$')
	pylab.legend(loc=2)
	pylab.savefig(env.env['tmp_dir'] + 'power_common.png')

	pylab.clf()
	#Disjoint fraction barchart.
	bars_xs = []
	bars_ys = []
	bar_labels = []
	curr_x = 0
	tick_positions = []

	pylab.axes([0.16, 0.22, 0.82, 0.73])

	for h in ['0.25', '0.50', '0.75']:
		for m, lab in [('marginal', 'Minimum marginal'), ('pher', 'Joint pseudo-herit.'), \
				('her', 'Joint herit.')]:
			if m != 'marginal':
				bars_ys.append(disjoint_fraction_found[h][m]['h01'])
			else:
				bars_ys.append(disjoint_fraction_found[h][m])
			bars_xs.append(curr_x)
			bar_labels.append(lab)
			tick_positions.append(curr_x + 0.4)
			curr_x += 1
		curr_x += 1
	curr_x = curr_x - 1.2
	pylab.bar(bars_xs, bars_ys, color="#1155BB", alpha=0.7)
	pylab.xticks(tick_positions, bar_labels, rotation=45, ha='right')
	pylab.ylabel('Fraction of sigmificant causative SNPs')
	pylab.axis([-0.02 * curr_x, 1.02 * curr_x, 0.0, 1.15 * max(bars_ys)])
	pylab.text(0.8, 1.05 * max(bars_ys), r'$h^2=0.25$')
	pylab.text(4.8, 1.05 * max(bars_ys), r'$h^2=0.5$')
	pylab.text(8.8, 1.05 * max(bars_ys), r'$h^2=0.75$')
	pylab.savefig(env.env['tmp_dir'] + 'power_disjoint.png')



	#Plot box-plots of median p-values..
	pher_ks = []
	her_ks = []
	marg_ks = []
	indep_ks = []

	pher_med = []
	her_med = []
	marg_med = []
	indep_med = []
	for her in ['0.25', '0.50', '0.75']:
		pher_ks.append(rd_map[her]['ks_stats_pher']['h01'])
		her_ks.append(rd_map[her]['ks_stats_her']['h01'])
		indep_ks.append(rd_map[her]['ks_stats_indep']['h01'])
		marg_ks.append(rd_map[her]['ks_stats_marginal'])

		pher_med.append(rd_map[her]['pval_medians_pher']['h01'])
		her_med.append(rd_map[her]['pval_medians_her']['h01'])
		indep_med.append(rd_map[her]['pval_medians_indep']['h01'])
		marg_med.append(rd_map[her]['pval_medians_marginal'])

	pylab.clf()
	pylab.figure(figsize=(5, 5))
	pylab.axes([0.16, 0.26, 0.82, 0.7])

	pylab.boxplot([marg_ks, pher_ks, her_ks])
	tick_positions = [1, 2, 3]
	tick_labels = ['Minimum marginal', 'Joint pseudo-herit.', 'Joint herit.']
	pylab.xticks(tick_positions, tick_labels, rotation=45, ha='right')
	pylab.ylabel('Kolmogorov-Smirnov statistic')
	pylab.savefig(env.env['tmp_dir'] + 'ks.png')

	pylab.clf()
	pylab.axes([0.16, 0.26, 0.82, 0.7])
	pylab.boxplot([marg_med, pher_med, her_med])
	tick_labels = ['Minimum marginal', 'Joint pseudo-herit.', 'Joint herit.']
	pylab.xticks(tick_positions, tick_labels, rotation=45, ha='right')
	pylab.ylabel('0.5-(median pvalue)')
	pylab.savefig(env.env['tmp_dir'] + 'median.png')





def run_parallel(heritability, x_start_i, x_stop_i, cluster='usc'):
	"""
	If no mapping_method, then analysis run is set up.
	"""
	run_id = 'corr_trait_sim'
	job_id = ' % s_ % d_ % d' % (run_id, x_start_i, x_stop_i)
	file_prefix = env.env['results_dir'] + run_id + '_' + str(x_start_i) + '_' + str(x_stop_i)

	#Cluster specific parameters	
	if cluster == 'gmi': #GMI cluster.
		shstr = '#!/bin/sh\n'
		shstr += '#$ -N %s\n' % job_id
		shstr += "#$ -q q.norm@blade*\n"
		shstr += '#$ -o %s.log\n' % job_id
		#shstr += '#$ -cwd /home/GMI/$HOME\n'
		#shstr += '#$ -M bjarni.vilhjalmsson@gmi.oeaw.ac.at\n\n'

	elif cluster == 'usc':  #USC cluster.
		shstr = "#!/bin/csh\n"
		shstr += "#PBS -l walltime=%s \n" % '72:00:00'
		shstr += "#PBS -l mem=%s \n" % '1950mb'
		shstr += "#PBS -q cmb\n"
		shstr += "#PBS -N p%s \n" % job_id

	shstr += "(python %scorr_trait_sim.py %s %d %d " % (env.env['script_dir'], heritability, x_start_i, x_stop_i)

	shstr += "> " + file_prefix + "_job.out) >& " + file_prefix + "_job.err\n"
	print '\n', shstr, '\n'
	script_file_name = run_id + ".sh"
	f = open(script_file_name, 'w')
	f.write(shstr)
	f.close()

	#Execute qsub script
	os.system("qsub " + script_file_name)



def run_parallel_gwas():
	if len(sys.argv) > 3:
		run_joint_analysis(int(sys.argv[2]), int(sys.argv[3]), sys.argv[1], debug_filter=1.0, run_id='joint_trait')
	else:
		chunck_size = int(sys.argv[2])
		for i in range(0, max_num, chunck_size):
			run_parallel(sys.argv[1], i, i + chunck_size)



def plot_nfbc_traits():
	import gwaResults as gr
	chrom_d, var_d = dp._parse_map_file_()
	results_dir = '/Users/bjarni.vilhjalmsson/Projects/Data/NFBC_results/'
	file_list = [#('HDL', 'CRP', 'multiple_traits/hdlres_crpres_EMMAX.csv', 'multiple_loci/hdlres_pval_steps.csv', 'multiple_loci/crpres_pval_steps.csv'),
			#('HDL', 'LDL', 'multiple_traits/hdlres_ldlres_EMMAX.csv', 'multiple_loci/hdlres_pval_steps.csv', 'multiple_loci/ldlres_pval_steps.csv'),
			#('LDL', 'CRP', 'multiple_traits/ldlres_crpres_EMMAX.csv', 'multiple_loci/ldlres_pval_steps.csv', 'multiple_loci/crpres_pval_steps.csv'),
			#('TG', 'CRP', 'multiple_traits/tgres_crpres_EMMAX.csv', 'multiple_loci/tgres_pval_steps.csv', 'multiple_loci/crpres_pval_steps.csv'),
			#('TG', 'HDL', 'multiple_traits/tgres_hdlres_EMMAX.csv', 'multiple_loci/tgres_pval_steps.csv', 'multiple_loci/hdlres_pval_steps.csv'),
			('TG', 'LDL', 'multiple_traits/tgres_ldlres_EMMAX.csv', 'multiple_loci/tgres_pval_steps.csv', 'multiple_loci/ldlres_pval_steps.csv', 15), ]

	chrom_col_map = {}
	for i in range(1, 24):
		if i % 2 == 0:
			chrom_col_map[i] = '#1199EE'
		else:
			chrom_col_map[i] = '#11BB00'

	for t1, t2, mt_pval_file, t1_pval_file, t2_pval_file, max_score in file_list:
		#One plot per file
		full_pvals = [] #6
		g_pvals = [] #8
		ge_pvals = [] #7
		chromosomes = []
		positions = []

		with open(results_dir + mt_pval_file) as f:
			print f.next()
			for line in f:
				l = line.split(',')
				var_id = l[1].replace('"', '')
				d = var_d[var_id]
				chromosomes.append(d['chrom'])
				positions.append(d['position'])
				full_pvals.append(float(l[6]))
				ge_pvals.append(float(l[7]))
				g_pvals.append(float(l[8]))
		l = zip(chromosomes, positions, full_pvals, g_pvals, ge_pvals)
		l.sort()
		l = map(list, zip(*l))
		full_res = gr.Result(scores=l[2], chromosomes=l[0][:], positions=l[1][:])
		g_res = gr.Result(scores=l[3], chromosomes=l[0][:], positions=l[1][:])
		ge_res = gr.Result(scores=l[4], chromosomes=l[0], positions=l[1])
		quantiles = []
		log_quantiles = []
		labels = []
		for r, n in [(full_res, 'Full model'), (g_res, 'Common effects model'), (ge_res, 'Differing effects model')]:
			s = r.get_scores()
			quantiles.append(agr.get_quantiles(s, 1000))
			log_quantiles.append(agr.get_log_quantiles(s, 1000, max_val=5.5))
			labels.append(n)
			r.neg_log_trans()

		#full_res.plot_qq(env.env['tmp_dir'] + 'mlt_HUMAN_full_%s' % name)
		#g_res.plot_qq(env.env['tmp_dir'] + 'mlt_HUMAN_g_%s' % name)
		#ge_res.plot_qq(env.env['tmp_dir'] + 'mlt_HUMAN_ge_%s' % name)

		marg_results = []
		for pval_file in [t1_pval_file, t2_pval_file]:
			chromosomes = []
			positions = []
			pvals = []
			with open(results_dir + pval_file) as f:
				print f.next()
				for line in f:
					l = line.split(',')
					var_id = l[0].replace('"', '')
					d = var_d[var_id]
					chromosomes.append(d['chrom'])
					positions.append(d['position'])
					pvals.append(float(l[1]))
			l = zip(chromosomes, positions, pvals)
			l.sort()
			l = map(list, zip(*l))
			res = gr.Result(scores=l[2], chromosomes=l[0], positions=l[1])
			marg_results.append(res)
		for r, n in [(marg_results[0], 'TG EMMAX'), (marg_results[1], 'LDL EMMAX')]:
			s = r.get_scores()
			quantiles.append(agr.get_quantiles(s, 1000))
			log_quantiles.append(agr.get_log_quantiles(s, 1000, max_val=5.5))
			labels.append(n)
			r.neg_log_trans()

		agr.simple_qqplot(quantiles, pdf_file=env.env['tmp_dir'] + '%s_%s_qq_plot.pdf' % (t1, t2),
				quantile_labels=labels, line_colors=['g', 'r', 'b', 'c', 'm'], num_dots=1000,
				plot_label='a')
		agr.simple_log_qqplot(log_quantiles, pdf_file=env.env['tmp_dir'] + '%s_%s_log_qq_plot.pdf' % (t1, t2),
				line_colors=['g', 'r', 'b', 'c', 'm'], num_dots=1000, max_val=5.5, plot_label='b')



		f = pylab.figure(figsize=(10, 10))
		ax1 = f.add_axes([0.05, 0.04 + (0.8 * 0.97), 0.945, (0.18 * 0.97) ])
		ax1.spines['top'].set_visible(False)
		ax1.xaxis.set_ticks_position('bottom')
		ax1.xaxis.set_label_position('bottom')
		#ax1.set_ylabel(r'$-$log$($p-value$)$')
		ax2 = f.add_axes([0.05, 0.04 + (0.6 * 0.97), 0.945, (0.18 * 0.97) ])
		ax2.spines['top'].set_visible(False)
		ax2.xaxis.set_ticks_position('bottom')
		ax2.xaxis.set_label_position('bottom')
		#ax2.set_ylabel(r'$-$log$($p-value$)$')
		ax3 = f.add_axes([0.05, 0.04 + (0.4 * 0.97), 0.945, (0.18 * 0.97) ])
		ax3.spines['top'].set_visible(False)
		ax3.xaxis.set_ticks_position('bottom')
		ax3.xaxis.set_label_position('bottom')
		ax3.set_ylabel(r'$-$log$($p-value$)$')
		ax4 = f.add_axes([0.05, 0.04 + (0.2 * 0.97), 0.945, (0.18 * 0.97) ])
		ax4.spines['top'].set_visible(False)
		ax4.xaxis.set_ticks_position('bottom')
		ax4.xaxis.set_label_position('bottom')
		#ax4.set_ylabel(r'$-$log$($p-value$)$')
		ax5 = f.add_axes([0.05, 0.04, 0.945, (0.18 * 0.97) ])
		ax5.spines['top'].set_visible(False)
		ax5.xaxis.set_ticks_position('bottom')
		ax5.xaxis.set_label_position('bottom')
		#ax25.set_ylabel(r'$-$log$($p-value$)$')


		marg_results[0].plot_manhattan2(ax=ax1, neg_log_transform=False, plot_bonferroni=True, plot_xaxis=True,
					chrom_colormap=chrom_col_map, sign_color='#DD1122', wo_xtick_labels=True,
					max_score=max_score)
		x_min, x_max = ax1.get_xlim()
		x_range = x_max - x_min
		y_min, y_max = ax1.get_ylim()
		y_range = y_max - y_min
		ax1.text(0.96 * x_range + x_min, 0.94 * y_range + y_min, 'a')

		marg_results[1].plot_manhattan2(ax=ax2, neg_log_transform=False, plot_bonferroni=True, plot_xaxis=True,
					chrom_colormap=chrom_col_map, sign_color='#DD1122', wo_xtick_labels=True,
					max_score=max_score)
		x_min, x_max = ax2.get_xlim()
		x_range = x_max - x_min
		y_min, y_max = ax2.get_ylim()
		y_range = y_max - y_min
		ax2.text(0.96 * x_range + x_min, 0.94 * y_range + y_min, 'b')

		full_res.plot_manhattan2(ax=ax3, neg_log_transform=False, plot_bonferroni=True, plot_xaxis=True,
					chrom_colormap=chrom_col_map, sign_color='#DD1122', wo_xtick_labels=True,
					max_score=max_score)
		x_min, x_max = ax3.get_xlim()
		x_range = x_max - x_min
		y_min, y_max = ax3.get_ylim()
		y_range = y_max - y_min
		ax3.text(0.96 * x_range + x_min, 0.94 * y_range + y_min, 'c')

		g_res.plot_manhattan2(ax=ax4, neg_log_transform=False, plot_bonferroni=True, plot_xaxis=True,
					chrom_colormap=chrom_col_map, sign_color='#DD1122', wo_xtick_labels=True,
					max_score=max_score)
		x_min, x_max = ax4.get_xlim()
		x_range = x_max - x_min
		y_min, y_max = ax4.get_ylim()
		y_range = y_max - y_min
		ax4.text(0.96 * x_range + x_min, 0.94 * y_range + y_min, 'd')

		ge_res.plot_manhattan2(ax=ax5, neg_log_transform=False, plot_bonferroni=True, plot_xaxis=True,
					chrom_colormap=chrom_col_map, sign_color='#DD1122', max_score=max_score)
		x_min, x_max = ax5.get_xlim()
		x_range = x_max - x_min
		y_min, y_max = ax5.get_ylim()
		y_range = y_max - y_min
		ax5.text(0.96 * x_range + x_min, 0.94 * y_range + y_min, 'e')

		pylab.savefig(env.env['tmp_dir'] + 'mlt_HUMAN_%s_%s.png' % (t1, t2), dpi=400)



def plot_yan_li_traits():
	data_file = '/Users/bjarni.vilhjalmsson/Projects/Data/Yan_Li_2009_log.csv'
	chromosomes = []
	positions = []
	comb_list = []
#	pvals = [[], [], [], []]
#	combined_pvals = [[], [], [], [], []]
	with open(data_file) as f:
		print f.next()
		for line in f:
			l = map(str.strip, line.split(','))
			res_list = [int(l[2]), int(l[3]), int(l[6])] + [float(l[i]) for i in range(8, 17)]
			comb_list.append(tuple(res_list))
	comb_list.sort()
	comb_list = map(list, zip(*comb_list))
	chromosomes = comb_list[0]
	positions = comb_list[1]
	macs = comb_list[2]
	labels = ['Spain spring', 'Spain summer', 'Sweden spring', 'Sweden summer', 'Full model', 'Differing effects',
			'Genotype x country', 'Genotype x season', 'Common effects']

	quantiles = []
	log_quantiles = []
	results = []
	for ps in comb_list[3:]:
		r = gr.Result(scores=ps, chromosomes=chromosomes, positions=positions, macs=macs)
		r.filter_attr('macs', 20)
		s = r.get_scores()
		quantiles.append(agr.get_quantiles(s, 1000))
		log_quantiles.append(agr.get_log_quantiles(s, 1000, max_val=5))
		r.neg_log_trans()
		results.append(r)

	real_files = ['bf_qtl_priors_DTF2ndSwAverage2009_317_step0.pvals', 'bf_qtl_priors_DTF1stSwAverage2009_316_step0.pvals',
			'bf_qtl_priors_DTF2ndSpAverage2009_315_step0.pvals', 'bf_qtl_priors_DTF1stSpAverage2009_314_step0.pvals']
	for fn , n in zip(real_files, ['Sweden summer', 'Sweden spring', 'Spain summer', 'Spain spring' ]):
		r = gr.Result(env.env['tmp_dir'] + fn)
		r.filter_attr('macs', 20)
		s = r.get_scores()
		quantiles.append(agr.get_quantiles(s, 1000))
		log_quantiles.append(agr.get_log_quantiles(s, 1000, max_val=5))
		labels.append(n)
		r.neg_log_trans()
		results.append(r)


	qs = quantiles[9:] + [quantiles[7], quantiles[6]]
	log_qs = log_quantiles[9:] + [log_quantiles[7], log_quantiles[6]]
	labs = labels[9:] + [labels[7], labels[6]]
	agr.simple_qqplot(qs, pdf_file=env.env['tmp_dir'] + 'mlt_FT_li_qq_plot.pdf',
			quantile_labels=labs, line_colors=['#6a5acd', '#0000cd', '#ff0000', '#ffa500', '#228b22', '#66cdaa', ], num_dots=1000,
			plot_label='a')
	agr.simple_log_qqplot(log_qs, pdf_file=env.env['tmp_dir'] + 'mlt_FT_li_log_qq_plot.pdf',
			line_colors=['#6a5acd', '#0000cd', '#ff0000', '#ffa500', '#228b22', '#66cdaa', ],
			num_dots=1000, max_val=5, plot_label='b')


	chrom_col_map = {}
	for i in range(1, 24):
		if i % 2 == 0:
			chrom_col_map[i] = '#1199EE'
		else:
			chrom_col_map[i] = '#11BB00'

	f = pylab.figure(figsize=(10, 10))
	ax1 = f.add_axes([0.05, 0.04 + (0.8 * 0.97), 0.945, (0.18 * 0.97) ])
	ax1.spines['top'].set_visible(False)
	ax1.xaxis.set_ticks_position('bottom')
	ax1.xaxis.set_label_position('bottom')
	#ax1.set_ylabel(r'$-$log$($p-value$)$')
	ax2 = f.add_axes([0.05, 0.04 + (0.6 * 0.97), 0.945, (0.18 * 0.97) ])
	ax2.spines['top'].set_visible(False)
	ax2.xaxis.set_ticks_position('bottom')
	ax2.xaxis.set_label_position('bottom')
	#ax2.set_ylabel(r'$-$log$($p-value$)$')
	ax3 = f.add_axes([0.05, 0.04 + (0.4 * 0.97), 0.945, (0.18 * 0.97) ])
	ax3.spines['top'].set_visible(False)
	ax3.xaxis.set_ticks_position('bottom')
	ax3.xaxis.set_label_position('bottom')
	ax3.set_ylabel(r'$-$log$($p-value$)$')
	ax4 = f.add_axes([0.05, 0.04 + (0.2 * 0.97), 0.945, (0.18 * 0.97) ])
	ax4.spines['top'].set_visible(False)
	ax4.xaxis.set_ticks_position('bottom')
	ax4.xaxis.set_label_position('bottom')
	#ax4.set_ylabel(r'$-$log$($p-value$)$')
	ax5 = f.add_axes([0.05, 0.04, 0.945, (0.18 * 0.97) ])
	ax5.spines['top'].set_visible(False)
	ax5.xaxis.set_ticks_position('bottom')
	ax5.xaxis.set_label_position('bottom')
	#ax25.set_ylabel(r'$-$log$($p-value$)$')

	max_score = 9

	results[9].plot_manhattan2(ax=ax1, neg_log_transform=False, plot_bonferroni=True, plot_xaxis=True,
				chrom_colormap=chrom_col_map, sign_color='#DD1122', wo_xtick_labels=True,
				max_score=max_score)
	x_min, x_max = ax1.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax1.get_ylim()
	y_range = y_max - y_min
	ax1.text(0.96 * x_range + x_min, 0.94 * y_range + y_min, 'a')

	results[10].plot_manhattan2(ax=ax2, neg_log_transform=False, plot_bonferroni=True, plot_xaxis=True,
				chrom_colormap=chrom_col_map, sign_color='#DD1122', wo_xtick_labels=True,
				max_score=max_score)
	x_min, x_max = ax2.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax2.get_ylim()
	y_range = y_max - y_min
	ax2.text(0.96 * x_range + x_min, 0.94 * y_range + y_min, 'b')

	results[11].plot_manhattan2(ax=ax3, neg_log_transform=False, plot_bonferroni=True, plot_xaxis=True,
				chrom_colormap=chrom_col_map, sign_color='#DD1122', wo_xtick_labels=True,
				max_score=max_score)
	x_min, x_max = ax3.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax3.get_ylim()
	y_range = y_max - y_min
	ax3.text(0.96 * x_range + x_min, 0.94 * y_range + y_min, 'c')

	results[12].plot_manhattan2(ax=ax4, neg_log_transform=False, plot_bonferroni=True, plot_xaxis=True,
				chrom_colormap=chrom_col_map, sign_color='#DD1122', wo_xtick_labels=True,
				max_score=max_score)
	x_min, x_max = ax4.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax4.get_ylim()
	y_range = y_max - y_min
	ax4.text(0.96 * x_range + x_min, 0.94 * y_range + y_min, 'd')

	results[7].plot_manhattan2(ax=ax5, neg_log_transform=False, plot_bonferroni=True, plot_xaxis=True,
				chrom_colormap=chrom_col_map, sign_color='#DD1122', max_score=max_score)
	x_min, x_max = ax5.get_xlim()
	x_range = x_max - x_min
	y_min, y_max = ax5.get_ylim()
	y_range = y_max - y_min
	ax5.text(0.96 * x_range + x_min, 0.94 * y_range + y_min, 'e')

	print results[7].get_top_snps(10, True)

	print results[12].get_top_snps(10, True)

	pylab.savefig(env.env['tmp_dir'] + 'mlt_FT_li.png', dpi=400)


def _test_():
	#sim_traits_dict = _load_sim_data_()
	#chunk_size = 2
#	for i in range(0, max_num, chunk_size):
#		run_joint_analysis(i, i + chunk_size, 0.75, debug_filter=0.1)
	plot_yan_li_traits()


if __name__ == '__main__':
	_test_()
	#run_parallel_gwas()


