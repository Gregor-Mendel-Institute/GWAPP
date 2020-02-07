"""
GWA project using full sequence data.

Quan, Bjarni, Dazhe, et al.
"""
import sys
import os
import env
import dataParsers as dp
import phenotypeData as pd
import linear_models as lm


def run_parallel(pid, call_method_id, run_id='gwas', kinship_method='ibd'):
        """
        If no mapping_method, then analysis run is set up.
        """
        job_id = '%s_%s_%d_%d' % (run_id, kinship_method, call_method_id, pid)
        file_prefix = env.env['results_dir'] + job_id

        #Cluster specific parameters        
        shstr = '#!/bin/bash\n'
        shstr += '#$ -S /bin/bash\n'
        shstr += '#$ -N %s\n' % job_id
        #shstr += '#$ -o %s_job_$JOB_ID.out\n' % file_prefix
        #shstr += '#$ -e %s_job_$JOB_ID.err\n' % file_prefix
        shstr += '#$ -o %s_job.out\n' % file_prefix
        shstr += '#$ -e %s_job.err\n' % file_prefix
        shstr += 'source /etc/modules-env.sh\n'
        shstr += 'module load scipy/GotoBLAS2/0.9.0\n'
        shstr += 'module load matplotlib/1.0.0\n'
        shstr += 'module load mysqldb/1.2.3\n'
	shstr += 'module load h5py/2.0.0\n'
        shstr += 'export GOTO_NUM_THREADS=1\n'


        shstr += "python %sfullseq_gwas_project.py %s %s %d %d" % \
                        (env.env['script_dir'], run_id, kinship_method, call_method_id, pid)

        #shstr += "> " + file_prefix + "_job.out) >& " + file_prefix + "_job.err\n"
        print '\n', shstr, '\n'
        script_file_name = run_id + ".sh"
        f = open(script_file_name, 'w')
        f.write(shstr)
        f.close()

        #Execute qsub script
        os.system("qsub " + script_file_name)



def run_gwas(pid, call_method_id, run_id, kinship_method, debug_filter=1):
        #import snpsdata

        #LOAD DATA
	sd = dp.load_snps_call_method(call_method_id)
	if debug_filter < 1:
		sd.sample_snps(debug_filter)
	phenotype_file = env.env['phen_dir'] + 'phen_with_swedish_082211.csv'
	phed = pd.parse_phenotype_file(phenotype_file)
	phed.convert_to_averages()
	phen_name = phed.get_name(pid)
	sd.coordinate_w_phenotype_data(phed, pid)
        phed.transform(pid, 'most_normal')
	phen_vals = phed.get_values(pid)

	if kinship_method == 'ibd':
		global_k = sd.get_ibd_kinship_matrix()
	elif kinship_method == 'ibs':
		global_k = sd.get_ibs_kinship_matrix()

        p_her = phed.get_pseudo_heritability(pid, global_k)
        hist_file = env.env['results_dir'] + '%s_%s_%d_%d_%s_hist.png' % \
                                                (run_id, kinship_method, call_method_id, pid, phen_name)

        phed.plot_histogram(pid, p_her=p_her, png_file=hist_file)

        #Set up GWAS

	#Chromosomes.
	res_dict = lm.chrom_vs_rest_mm(phen_vals, sd, kinship_method, global_k)
	print res_dict
	file_prefix = env.env['results_dir'] + '%s_loc_v_glob_chrom_%s_%d_%d_%s' % \
						(run_id, kinship_method, call_method_id, pid, phen_name)
	res_file_name = file_prefix + '.csv'
	_write_res_dict_to_file_2_(res_file_name, res_dict)

	#Now 'normal' window sizes
	for ws in [3000000, 1000000, 500000, 200000, 100000, 50000, 20000]:
		file_prefix = env.env['results_dir'] + '%s_loc_v_glob_%s_%d_%d_%d_%s' % \
							(run_id, kinship_method, call_method_id, ws, pid, phen_name)
		res_dict = lm.local_vs_global_mm_scan(phen_vals, sd, file_prefix, ws, ws / 2, kinship_method, global_k)
		res_file_name = file_prefix + '.csv'
		_write_res_dict_to_file_(res_file_name, res_dict)

	#Now gene-centralized.
	for radius in [20000, 10000, 5000]:
		file_prefix = env.env['results_dir'] + '%s_loc_v_glob_gene_%s_%d_%d_%d_%s' % \
							(run_id, kinship_method, call_method_id, radius, pid, phen_name)
		res_dict = lm.local_vs_global_gene_mm_scan(phen_vals, sd, file_prefix, radius, kinship_method, global_k)
		res_file_name = file_prefix + '.csv'
		_write_res_dict_to_file_3_(res_file_name, res_dict)

        sd.filter_mac_snps(15)
        file_prefix = env.env['results_dir'] + '%s_emmax_stepwise_%s_%d_%d_%s' % \
                                                (run_id, kinship_method, call_method_id, pid, phen_name)
        lm.emmax_step_wise(phen_vals, global_k, sd=sd, num_steps=10, file_prefix=file_prefix, save_pvals=True)



def _write_res_dict_to_file_(filename, rd):
	with open(filename, 'w') as f:
		f.write('h0_heritability: %f\n' % rd['h0_heritability'])
		f.write('chromosomes, positions, pvalues, perc_variance_local, perc_variance_global, h1_heritabilities\n')
		num_res = len(rd['chromosomes'])
		for i in range(num_res):
			f.write('%d, %d, %s, %f, %f, %f\n' % (rd['chromosomes'][i], rd['positions'][i], str(rd['pvals'][i]),
						rd['perc_variances2'][i], rd['perc_variances1'][i], rd['h1_heritabilities'][i]))

def _write_res_dict_to_file_2_(filename, rd):
	with open(filename, 'w') as f:
		f.write('h0_heritability: %f\n' % rd['h0_heritability'])
		f.write('chromosomes, pvalues, perc_variance_local, perc_variance_global, h1_heritabilities\n')
		num_res = len(rd['chromosomes'])
		for i in range(num_res):
			f.write('%d, %s, %f, %f, %f\n' % (rd['chromosomes'][i], str(rd['pvals'][i]), rd['perc_variances2'][i],
							rd['perc_variances1'][i], rd['h1_heritabilities'][i]))

def _write_res_dict_to_file_3_(filename, rd):
	with open(filename, 'w') as f:
		f.write('h0_heritability: %f\n' % rd['h0_heritability'])
		f.write('tair_ids, chromosomes, positions, pvalues, perc_variance_local, perc_variance_global, h1_heritabilities\n')
		num_res = len(rd['chromosomes'])
		for i in range(num_res):
			f.write('%s, %d, %d, %s, %f, %f, %f\n' % (rd['tair_ids'][i], rd['chromosomes'][i], rd['positions'][i],
								str(rd['pvals'][i]), rd['perc_variances2'][i],
								rd['perc_variances1'][i], rd['h1_heritabilities'][i]))


def telomere_example_plots(debug_filter=1.0, pid=1365, call_method_id=78, radius=20000, kinship_method='ibs'):
        genes_of_interest = ['AT1G21390', 'AT1G21400', 'AT1G21410', 'AT1G21420', 'AT1G21430', 'AT1G21440', 'AT1G21450',
                             'AT1G21460', 'AT1G21470', 'AT1G21480', 'AT1G21490']
        sd = dp.load_snps_call_method(call_method_id)
        if debug_filter < 1:
                sd.sample_snps(debug_filter)
        phenotype_file = env.env['phen_dir'] + 'phen_with_swedish_082211.csv'
        phed = pd.parse_phenotype_file(phenotype_file)
        phed.convert_to_averages()
        phen_name = phed.get_name(pid)
        sd.coordinate_w_phenotype_data(phed, pid)
        phed.transform(pid, 'most_normal')
        png_file = env.env['results_dir'] + 'histogram_%s_hist.png' % phed.get_name(pid)
        phed.plot_histogram(pid, png_file=png_file)
        phen_vals = phed.get_values(pid)
        file_prefix = env.env['results_dir'] + 'loc_v_glob_gene_%d_%d_%d_%s' % \
                                                (call_method_id, radius, pid, phen_name)
        res_dict = lm.local_vs_global_gene_mm_scan(phen_vals, sd, file_prefix, radius, kinship_method,
                                                   tair_ids=genes_of_interest, plot_gene_trees=True, ets=sd.accessions)

def run():
	phenotype_file = env.env['phen_dir'] + 'phen_with_swedish_082211.csv'
	run_id = sys.argv[1]
        kinship_method = sys.argv[2]
	call_method_id = int(sys.argv[3])
        if len(sys.argv) < 5:
                print 'Setting up a cluster run'
        	phed = pd.parse_phenotype_file(phenotype_file)
		pids = phed.phen_ids
		for pid in pids:
			run_parallel(pid, call_method_id, run_id, kinship_method)


	else:
                print 'Setting up a test run'
		pid = int(sys.argv[4])
		run_gwas(pid, call_method_id, run_id, kinship_method)



if __name__ == '__main__':
        #telomere_example_plots()
        run()
