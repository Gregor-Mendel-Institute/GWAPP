"""
Basic analysis of RNA data.

Bjarni J. Vilhjalmsson (bjarni.vilhjalmsson@gmi.oeaw.ac.at)
"""

#TODO
# - Run calculations on cluster (using full sequence data):
# - Draw cool plots..
# - Using pylons, create interactive web based plots.

import phenotypeData as pd
import dataParsers as dp
import gwaResults as gr
import linear_models as lm
import scipy as sp
from env import *
import sys
import cPickle
import util
import gwaResults
import pylab
import pdb
import itertools as it
import math
import analyze_gwas_results as agr
import h5py
import random
import snpsdata
import bisect
import time
import kinship

#Annoyingly ad hoc constants
near_const_filter = 20
phen_file_prefix = env['phen_dir'] + 'rna_seq_081411'
#phen_file_prefix = env['phen_dir'] + 'rna_seq_061611'


def run_parallel(x_start_i, x_stop_i, temperature, call_method_id, conf_matrix_type, cluster='gmi', run_id='rna_seq'):
    """
    If no mapping_method, then analysis run is set up.
    """
    job_id = '%s_%d_%d_%s' % (run_id, x_start_i, x_stop_i, temperature)
    file_prefix = env['results_dir'] + job_id

    #Cluster specific parameters    
    if cluster == 'gmi': #GMI cluster.  
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

    elif cluster == 'usc':  #USC cluster.
        shstr = "#!/bin/csh\n"
        shstr += "#PBS -l walltime=%s \n" % '72:00:00'
        shstr += "#PBS -l mem=%s \n" % '2950mb'
        shstr += "#PBS -q cmb\n"
        shstr += "#PBS -N p%s \n" % job_id

    shstr += "python %sanalyze_rna_seq.py %d %d %s %d %s %s" % \
            (env['script_dir'], x_start_i, x_stop_i, temperature, call_method_id, run_id, conf_matrix_type)

    #shstr += "> " + file_prefix + "_job.out) >& " + file_prefix + "_job.err\n"
    print '\n', shstr, '\n'
    script_file_name = '/tmp/' + run_id + '.sh'
    f = open(script_file_name, 'w')
    f.write(shstr)
    f.close()

    #Execute qsub script
    os.system("qsub " + script_file_name)







def summarize_stepwise(summary_dict, gene, step_info_list, opt_dict):
    #Store results for PPAs, MBONF, EBIC
    sw_d = {}
    for criteria in ['ebics', 'mbonf', 'min_cof_ppa']:
        i_opt = opt_dict[criteria]
        step_info = step_info_list[i_opt]
        cof_list = step_info['cofactors']
        ppa_cof_list = step_info['ppa_cofactors']
        cofactors = [(chrom, pos) for chrom, pos, pval in cof_list]
        cof_pvals = [pval for chrom, pos, pval in cof_list]
        cof_ppas = [ppa for chrom, pos, ppa in ppa_cof_list]

        # How close to the gene are the selected loci.
        da = []
        for chrom, pos in cofactors:
            if gene.chromosome == chrom:
                if gene.startPos < pos < gene.endPos:
                    da.append(0)
                else:
                    da.append(min(abs(pos - gene.startPos), abs(pos - gene.endPos)))
            else:
                da.append(-1)

        # Number of selected loci near the gene (bins)
        bin_distances = [100000, 50000, 25000, 10000, 5000, 1000, 0]
        bin_counts = [da.count(-1)]
        bin_count = 0
        for d in da:
            if d > bin_distances[0]: bin_count += 1
        bin_counts.append(bin_count)
        for bin_dist in bin_distances:
            bin_count = 0
            for d in da:
                if d == -1: continue
                elif d <= bin_dist: bin_count += 1
            bin_counts.append(bin_count)

        # Percentage of variance (error and genetic) explained by the selected loci (bin it as well)
        pass

        d = {'cofactors':cofactors, 'cof_pvals':cof_pvals, 'cof_ppas':cof_ppas, 'cof_gene_dist':da,
        'bin_counts':bin_counts, 'i_opt':i_opt}
        #print d
        sw_d[criteria] = d
        sw_d['step_info_list'] = step_info_list
    summary_dict['SW'] = sw_d



def summarize_st_results(temperature, call_method_id=79, debug_filter=1.0, conf_matrix_type='none',
        run_id='rna_seq', mms=['EX']):
    """
    Summarize the single trait results
    
    (works on the cluster)
    """
    #Figure out prefix
    file_prefix = env['rna_seq_results_dir'] + '%s/cm_%d/%s/' % (temperature, call_method_id, conf_matrix_type)

    f = h5py.File(env['rna_seq_data_dir'] + 'rna_seq_gwas_data_cm%d.h5py' % call_method_id)
    g = f[temperature]
    num_genes = g['gene_ids'].len()
    if debug_filter < 1.0:
        gene_filter = random.sample(range(num_genes), int(round(num_genes * debug_filter)))
        gene_filter.sort()
        gene_ids = g['gene_ids'][gene_filter].tolist()
    else:
        gene_ids = g['gene_ids'][...].tolist()
    num_snps = g['snps'].len()
    bonf_threshold = 1.0 / (20.0 * num_snps)
    #phen_vals_list = g['phen_values'][gene_filter]
    #ets = g['ecotypes'][...].tolist()

#    gene_dict = dp.parse_tair_gff_file()#_load_genes_list_('rna_seq_031311_%sC' % temperature)

    radii = [500000, 100000, 50000, 25000, 10000, 5000, 1000, 0]
    tss_dists = [200000, 100000, 50000, 25000, 10000, 5000, 1000]
    num_radii = len(radii)
    num_tss_dists = len(tss_dists)
    if conf_matrix_type == 'none':
        cvt_summary_dict = {'radius':{'avg_cis_trans_var_ratio':[0.0] * num_radii,
                        'avg_cis_herit':[0.0] * num_radii,
                        'avg_trans_herit':[0.0] * num_radii,
                        'counts':[0.0] * num_radii},
                    'radius_herit':{'avg_cis_trans_var_ratio':[0.0] * num_radii,
                        'avg_cis_herit':[0.0] * num_radii,
                        'avg_trans_herit':[0.0] * num_radii,
                        'counts':[0.0] * num_radii},
                    'tss_dist':{'avg_cis_trans_var_ratio':[0.0] * num_tss_dists,
                        'avg_cis_herit':[0.0] * num_tss_dists,
                        'avg_trans_herit':[0.0] * num_tss_dists,
                        'counts':[0.0] * num_tss_dists}}

    heritabilities = []
    pval_infl_dict = {}
    distance_bins = [(0, 5000), (0, 10000), (0, 25000), (0, 50000), (0, 100000), (1, -1), (6, -1)]
    radius_bins = [0, 1000, 5000, 10000, 25000, 50000, 100000]
    bonf_sign_bin_dict = {}


    sign_count_ex = 0
    min_pval_dist_dict = {}
    min_pvals_ex = [] #Only EX
    for mm in mms:
        min_pval_dist_dict[mm] = []
        pval_infl_dict[mm] = {'kolmogorov_smirnov':[], 'median_pvals':[]}
        bonf_sign_bin_dict[mm] = {}
        for bin_radius in radius_bins:
            bonf_sign_bin_dict[mm][bin_radius] = {'count':0.0, 'total':0.0}


    cofactor_count_dict = {}
    for criteria in ['ebics', 'mbonf', 'min_cof_ppa']:
        cofactor_count_dict[criteria] = {'num_cofactor_list':[], 'bin_counts':sp.zeros(9),
                        'num_cis_cofactor_list':[], 'num_found':0}


    num_cis_peaks = 0
    num_trans_peaks = 0

    #Iterate over result files
    for i, gene_tair_id in enumerate(gene_ids):

        chrom = int(gene_tair_id[2])
        chr_prefix = '%schr_%d/rna_seq_%s_both_pid%d_%s' % (file_prefix, chrom, run_id, i, gene_tair_id)
        summary_file = chr_prefix + '_summary.pickled'
        #Check if the file is already there.
        print 'Checking %s' % (summary_file)
        if not os.path.isfile(summary_file):
            print 'File is missing... moving on :%s' % summary_file
            continue

        #Load the summary dictionary
        f = open(summary_file, 'r')
        sd = cPickle.load(f)
        f.close()

        #Summarize interesting statistics

        #Heritabilities
        heritabilities.append(sd['pseudo_heritability'])


        #cis vs. trans stuff
        if conf_matrix_type == 'none':
            cvt_dict = sd['CVT']
            for r_i, r in enumerate(radii):
                if cvt_dict['radius'][r] != None:
                    pvg = cvt_dict['radius'][r]['perc_var1']
                    pvl = cvt_dict['radius'][r]['perc_var2']
                    herit = cvt_dict['radius'][r]['pseudo_heritability1']
                    cvt_summary_dict['radius']['avg_cis_trans_var_ratio'][r_i] += pvl / (pvl + pvg)
                    cvt_summary_dict['radius']['avg_cis_herit'][r_i] += pvl * herit
                    cvt_summary_dict['radius']['avg_trans_herit'][r_i] += pvg * herit
                    cvt_summary_dict['radius']['counts'][r_i] += 1.0

            for r_i, r in enumerate(radii):
                if cvt_dict['radius'][r] != None:
                    herit = cvt_dict['radius'][r]['pseudo_heritability1']
                    if herit > 0.05:
                        pvg = cvt_dict['radius'][r]['perc_var1']
                        pvl = cvt_dict['radius'][r]['perc_var2']
                        cvt_summary_dict['radius_herit']['avg_cis_trans_var_ratio'][r_i] += pvl / (pvl + pvg)
                        cvt_summary_dict['radius_herit']['avg_cis_herit'][r_i] += pvl * herit
                        cvt_summary_dict['radius_herit']['avg_trans_herit'][r_i] += pvg * herit
                        cvt_summary_dict['radius_herit']['counts'][r_i] += 1.0

            for td_i, td in enumerate(tss_dists):
                if cvt_dict['tss_upstream'][td] != None:
                    pvg = cvt_dict['tss_upstream'][td]['perc_var1']
                    pvl = cvt_dict['tss_upstream'][td]['perc_var2']
                    herit = cvt_dict['tss_upstream'][td]['pseudo_heritability1']
                    cvt_summary_dict['tss_dist']['avg_cis_trans_var_ratio'][td_i] += pvl / (pvl + pvg)
                    cvt_summary_dict['tss_dist']['avg_cis_herit'][td_i] += pvl * herit
                    cvt_summary_dict['tss_dist']['avg_trans_herit'][td_i] += pvg * herit
                    cvt_summary_dict['tss_dist']['counts'][td_i] += 1.0

        #Specific results for the mapping methods
        for mm in mms:
            pval_infl_dict[mm]['kolmogorov_smirnov'].append(sd[mm]['kolmogorov_smirnov']['D'])
            pval_infl_dict[mm]['median_pvals'].append(sd[mm]['pval_median'])
            min_pval_dist_dict[mm].append(tuple(sd[mm]['dist_to_min_pval']))
            for bin_radius in radius_bins:
                pval = sd[mm]['bin_dict'][bin_radius]['min_pval']
                num_snps = sd[mm]['bin_dict'][bin_radius]['num_snps']
                if num_snps > 0:
                    bonf_sign_bin_dict[mm][bin_radius]['total'] += 1
                    if pval < 1.0 / (20 * num_snps):
                        bonf_sign_bin_dict[mm][bin_radius]['count'] += 1

        #Stepwise stuff 
        for criteria in ['ebics', 'mbonf', 'min_cof_ppa']:
            num_cofactors = len(sd['SW'][criteria]['cofactors'])
            cofactor_count_dict[criteria]['num_cofactor_list'].append(num_cofactors)
            if num_cofactors > 0:
                cofactor_count_dict[criteria]['num_found'] += 1
                cofactor_count_dict[criteria]['bin_counts'] += sp.array(sd['SW'][criteria]['bin_counts'])
                cofactor_count_dict[criteria]['num_cis_cofactor_list'].append(sd['SW'][criteria]['bin_counts'][2])

        min_pval_ex = sd['SW']['step_info_list'][0]['min_pval']
        min_pvals_ex.append(min_pval_ex)
        if min_pval_ex < bonf_threshold:
            sign_count_ex += 1

        #Count cis vs. trans, when looking at the smallest p-value
        min_pval_dist = min_pval_dist_dict['EX'][-1]
        chrom_dist = min_pval_dist[0]
        pos_dist = min_pval_dist[1]
        if chrom_dist == 0 and pos_dist < 100000:
            num_cis_peaks += 1
        else:
            num_trans_peaks += 1


    #pre-process the cis vs. trans stuff
    if conf_matrix_type == 'none':
        for r_i, r in enumerate(radii):
            r_counts = cvt_summary_dict['radius']['counts'][r_i]
            cvt_summary_dict['radius']['avg_cis_trans_var_ratio'][r_i] = \
                cvt_summary_dict['radius']['avg_cis_trans_var_ratio'][r_i] / r_counts
            cvt_summary_dict['radius']['avg_cis_herit'][r_i] = \
                cvt_summary_dict['radius']['avg_cis_herit'][r_i] / r_counts
            cvt_summary_dict['radius']['avg_trans_herit'][r_i] = \
                cvt_summary_dict['radius']['avg_trans_herit'][r_i] / r_counts


        for r_i, r in enumerate(radii):
            r_counts = cvt_summary_dict['radius_herit']['counts'][r_i]
            cvt_summary_dict['radius_herit']['avg_cis_trans_var_ratio'][r_i] = \
                cvt_summary_dict['radius_herit']['avg_cis_trans_var_ratio'][r_i] / r_counts
            cvt_summary_dict['radius_herit']['avg_cis_herit'][r_i] = \
                cvt_summary_dict['radius_herit']['avg_cis_herit'][r_i] / r_counts
            cvt_summary_dict['radius_herit']['avg_trans_herit'][r_i] = \
                cvt_summary_dict['radius_herit']['avg_trans_herit'][r_i] / r_counts


        for td_i, td in enumerate(tss_dists):
            td_counts = cvt_summary_dict['tss_dist']['counts'][td_i]
            cvt_summary_dict['tss_dist']['avg_cis_trans_var_ratio'][td_i] = \
                cvt_summary_dict['tss_dist']['avg_cis_trans_var_ratio'][td_i] / td_counts
            cvt_summary_dict['tss_dist']['avg_cis_herit'][td_i] = \
                cvt_summary_dict['tss_dist']['avg_cis_herit'][td_i] / td_counts
            cvt_summary_dict['tss_dist']['avg_trans_herit'][td_i] = \
                cvt_summary_dict['tss_dist']['avg_trans_herit'][td_i] / td_counts

    #Now plot stuff (better than before)
    results_prefix = env['results_dir'] + 'RNAseq_summary_%s_cm%d_%s' % (temperature, call_method_id, conf_matrix_type)

    if conf_matrix_type == 'none':
        #cis vs. trans
        pylab.figure()
        pylab.plot(cvt_summary_dict['radius']['avg_cis_trans_var_ratio'])
        pylab.ylabel('Avg. perc. of cis genetic var.')
        pylab.xlabel('Dist. from gene (kb)')
        pylab.xticks(range(8), [500, 100, 50, 25, 10, 5, 1, 0])
        pylab.savefig(results_prefix + '_avg_perc_cis_gen_var_rad.png')
        pylab.clf()

        pylab.figure()
        pylab.plot(cvt_summary_dict['tss_dist']['avg_cis_trans_var_ratio'])
        pylab.ylabel('Avg. perc. of cis genetic var.')
        pylab.xlabel('Dist. upstream from gene TSS (kb)')
        pylab.xticks(range(7), [200, 100, 50, 25, 10, 5, 1])
        pylab.savefig(results_prefix + '_avg_perc_cis_gen_var_td.png')
        pylab.clf()

        tot_herit = sp.array(cvt_summary_dict['radius']['avg_cis_herit']) + \
            sp.array(cvt_summary_dict['radius']['avg_trans_herit'])
        cis_herit = sp.array(cvt_summary_dict['radius']['avg_cis_herit'])
        pylab.figure(figsize=(10, 6))
        pylab.axes([0.06, 0.08, 0.92, 0.90])
        pylab.fill_between([0, 7], 0, 1, color='#DD3333', alpha=0.8, label='Error')
        pylab.fill_between(sp.arange(8), 0, tot_herit, color='#22CC44', alpha=0.8, label='Heritable variance')
        pylab.fill_between(sp.arange(8), 0, cis_herit, color='#2255AA', \
                    alpha=0.8, label='Heritable variance (cis)')
        pylab.ylabel('Average partition of variance')
        pylab.xlabel('Dist. from gene (kb)')
        pylab.xticks(range(8), [500, 100, 50, 25, 10, 5, 1, 0])
        pylab.legend(loc=1, ncol=3, shadow=True)
        pylab.axis([0, 7, 0, 1])
        pylab.savefig(results_prefix + 'avg_herit_rad.png')

        tot_herit = sp.array(cvt_summary_dict['radius_herit']['avg_cis_herit']) + \
            sp.array(cvt_summary_dict['radius_herit']['avg_trans_herit'])
        cis_herit = sp.array(cvt_summary_dict['radius_herit']['avg_cis_herit'])
        pylab.figure(figsize=(10, 6))
        pylab.axes([0.06, 0.08, 0.92, 0.90])
        pylab.fill_between([0, 7], 0, 1, color='#DD3333', alpha=0.8, label='Error')
        pylab.fill_between(sp.arange(8), 0, tot_herit, color='#22CC44', alpha=0.8, label='Heritable variance')
        pylab.fill_between(sp.arange(8), 0, cis_herit, color='#2255AA', \
                    alpha=0.8, label='Heritable variance (cis)')
        pylab.ylabel('Average partition of variance')
        pylab.xlabel('Dist. from gene (kb)')
        pylab.xticks(range(8), [500, 100, 50, 25, 10, 5, 1, 0])
        pylab.legend(loc=1, ncol=3, shadow=True)
        pylab.axis([0, 7, 0, 1])
        pylab.savefig(results_prefix + 'avg_herit_2_rad.png')



        tot_herit = sp.array(cvt_summary_dict['tss_dist']['avg_cis_herit']) + \
            sp.array(cvt_summary_dict['tss_dist']['avg_trans_herit'])
        cis_herit = sp.array(cvt_summary_dict['tss_dist']['avg_cis_herit'])
        pylab.figure(figsize=(10, 6))
        pylab.axes([0.06, 0.08, 0.92, 0.90])
        pylab.fill_between([0, 6], 0, 1, color='#DD3333', alpha=0.8, label='Error')
        pylab.fill_between(sp.arange(7), 0, tot_herit, color='#22CC44', alpha=0.8, label='Heritable variance')
        pylab.fill_between(sp.arange(7), 0, cis_herit, color='#2255AA', \
                    alpha=0.8, label='Heritable variance (cis)')
        pylab.ylabel('Average partition of variance')
        pylab.xlabel('Dist. upstream from gene TSS (kb)')
        pylab.xticks(range(7), [200, 100, 50, 25, 10, 5, 1])
        pylab.legend(loc=1, ncol=3, shadow=True)
        pylab.axis([0, 6, 0, 1])
        pylab.savefig(results_prefix + 'avg_herit_td.png')



    pylab.figure()
    pylab.hist(heritabilities, bins=20, alpha=0.7)
    pylab.xlabel('Pseudo-heritability')
    pylab.xlim((-0.025, 1.025))
    pylab.savefig(results_prefix + '_herits_hist.png')
    pylab.clf()

    ks_list = []
    pm_list = []
    for mm in mms:
        ks_list.append(pval_infl_dict[mm]['kolmogorov_smirnov'])
        pm_list.append(pval_infl_dict[mm]['median_pvals'])

    png_file_name = results_prefix + '_kolmogorov_smirnov_boxplot.png'
    pylab.figure()
    pylab.boxplot(ks_list)
    pylab.axhline(0, color='k', alpha=0.6, ls='-.')
    pylab.xticks(range(1, len(mms) + 1), mms)
    pylab.ylabel('Kolmogorov-Smirnov statistic D.')
    pylab.savefig(png_file_name)
    pylab.clf()


    png_file_name = results_prefix + '_median_pvals_boxplot.png'
    pylab.figure()
    pylab.boxplot(pm_list)
    pylab.axhline(0, color='k', alpha=0.6, ls='-.')
    pylab.xticks(range(1, len(mms) + 1), mms)
    pylab.ylabel('Median p-value bias')
    pylab.savefig(png_file_name)
    pylab.clf()


#    x_positions = sp.arange(len(distance_bins), dtype='d64')
#    width = 0.25
#    png_file_name = results_prefix + '_dist_min_pval_hist.png'
#    pylab.axes([0.08, 0.2, 0.91, 0.75])
#    for mm, color in zip(['EX', 'LM'], ['b', 'c']):
#        l = [dist_min_pval_dict[mm][dist_bin] for dist_bin in distance_bins]
#        tot_sum = sum(l)
#        l = map(lambda x: x / float(tot_sum), l)
#        pylab.bar(x_positions, l, width, color=color, alpha=0.7, label=mm)
#        x_positions += width
#
#
#    pylab.ylabel('Frequency')
#    pylab.xticks(x_positions - 3 * width / 2.0, (r'$d \leq 5$', r'$5< d \leq 10$', r'$10< d \leq 25$', \
#                        r'$25< d \leq 50$', r'$50< d \leq 100$', r'$d>100$', \
#                        'Other chrom.'), rotation='45')
#    pylab.xlabel('Distance $d$ (kb) to the smallest p-value from the gene.')
#    pylab.xlim((-0.25, len(distance_bins)))
#    pylab.legend(loc=2)
#    pylab.savefig(png_file_name)
#    pylab.clf()


    x_positions = sp.arange(len(radius_bins), dtype='d64')
    width = 0.25
    png_file_name = results_prefix + 'bonf_sign_bin_hist.png'
    pylab.axes([0.08, 0.22, 0.91, 0.73])
    for mm, color in zip(mms, ['b', 'c']):
        l = [bonf_sign_bin_dict[mm][b]['count'] / bonf_sign_bin_dict[mm][b]['total'] for b in radius_bins]
#        l.append(sign_count[mm] / float(num_genes))
        pylab.bar(x_positions, l, width, color=color, alpha=0.7, label=mm)
        x_positions += width


    pylab.ylabel('Fraction of sign. results')
    pylab.xticks(x_positions - 3 * width / 2.0, ('Within gene', r'$d \leq 1$', r'$d \leq 5$', \
                        r'$d \leq 10$', r'$d \leq 25$', r'$d \leq 50$', \
                        r'$d \leq 100$'), rotation='45')
    pylab.xlabel(r'Among SNPs with distance $d$ (kb) from gene.')
    pylab.xlim((-0.25, len(radius_bins)))
    pylab.legend(loc=2)
    pylab.savefig(png_file_name)
    pylab.clf()


    png_file_name = results_prefix + 'cofactor_count_hist.png'
    x_positions = sp.arange(6, dtype='d64')
    width = 0.25
    for criteria, color in zip(['ebics', 'mbonf', 'min_cof_ppa'], ['b', 'c', 'g']):
        bin_counts = list(sp.bincount(cofactor_count_dict[criteria]['num_cofactor_list']))
        while len(bin_counts) < 6:
            bin_counts.append(0)
        pylab.bar(x_positions, bin_counts, width, color=color, alpha=0.7, label=criteria)
        x_positions += width
    pylab.xlabel('Number of cofactor SNPs')
    pylab.ylabel('Number of genes')
    pylab.xticks(x_positions - 3 * width / 2.0, ('0', '1', '2', '3', '4', '5'))
    pylab.legend(loc=1)
    pylab.xlim((-0.2, 6))
    pylab.savefig(png_file_name)
    pylab.clf()


    png_file_name = results_prefix + 'cis_cofactor_count_hist.png'
    x_positions = sp.arange(6, dtype='d64')
    for criteria, color in zip(['ebics', 'mbonf', 'min_cof_ppa'], ['b', 'c', 'g']):
        bin_counts = list(sp.bincount(cofactor_count_dict[criteria]['num_cis_cofactor_list']))
        while len(bin_counts) < 6:
            bin_counts.append(0)
        pylab.bar(x_positions, bin_counts, width, color=color, alpha=0.7, label=criteria)
        x_positions += width
    pylab.xlabel('Number of cis cofactor SNPs')
    pylab.ylabel('Number of genes')
    pylab.xticks(x_positions - 3 * width / 2.0, ('0', '1', '2', '3', '4', '5'))
    pylab.legend(loc=1)
    pylab.xlim((-0.2, 6))
    pylab.savefig(png_file_name)
    pylab.clf()


    png_file_name = results_prefix + 'cofactor_bin_count_hist.png'
    x_positions = sp.arange(9, dtype='d64')
    width = 0.25
    pylab.axes([0.08, 0.2, 0.91, 0.75])
    for criteria, color in zip(['ebics', 'mbonf', 'min_cof_ppa'], ['b', 'c', 'g']):
        cofactor_count_dict[criteria]['bin_counts'] = \
            cofactor_count_dict[criteria]['bin_counts'] / cofactor_count_dict[criteria]['num_found']
        l = list(cofactor_count_dict[criteria]['bin_counts'])
        l.reverse()
        pylab.bar(x_positions, l, width, color=color, alpha=0.7, label=criteria)
        x_positions += width
    pylab.ylabel('Fraction all genes with cofactors.')
    pylab.xlabel(r'Distance $d$ (kb) to cofactor from gene.')
    pylab.xticks(x_positions - 3 * width / 2.0, ('Within gene', r'$1\geq d$', r'$5\geq d$', r'$10\geq d$', \
                        r'$25\geq d$', r'$50\geq d$', r'$100\geq d$', \
                        r'$d>100$', 'Other chrom.'), rotation='45')
    pylab.xlim((-0.2, 9))
    pylab.legend(loc=2)
    pylab.savefig(png_file_name)
    pylab.clf()

    print 'cis: %d' % num_cis_peaks
    print 'trans: %d' % num_trans_peaks
    #What else???
    #Distance to min p-value...


    #Plot across temperatures
    #    - Number of cis peaks and trans peaks
    #    - cis vs. trans using different things... 
    #    - cis/trans ratios for different...





def run_mtmm_gwas(file_prefix, start_i, stop_i, call_method_id=79, debug_filter=1.0, w_conf_matrix=False,
        w_marginal_analysis=False, run_id='blurb'):
    sys.path.append(env['script_dir'] + '.. / .. / pygp / ')
    sys.path.append(env['script_dir'] + '.. / .. / cegs / code / ')
    import fit_covariance_models as fcm

    f = h5py.File(env['rna_seq_data_dir'] + 'rna_seq_gwas_data_cm % d.h5py' % call_method_id)
    if w_conf_matrix:
        conf_f = h5py.File(' / srv / lab / data / rna_seq / standard / confounders_cm % d.hdf5' % call_method_id)
        Kexpr = conf_f['both']['Kexpr'][...]
    else:
        Kexpr = None
    g = f['both']

    if debug_filter < 1.0:
        num_snps = len(g['positions'])
        filter = random.sample(range(num_snps), int(round(num_snps * debug_filter)))
        filter.sort()
        snps = g['snps'][...][filter]
        mafs = g['mafs'][...][filter]
        macs = g['macs'][...][filter]
        positions = g['positions'][...][filter]
        chromosomes = g['chromosomes'][...][filter]
    else:
        snps = g['snps'][...]
        mafs = g['mafs'][...]
        macs = g['macs'][...]
        positions = g['positions'][...]
        chromosomes = g['chromosomes'][...]
    snps = [snp for snp in snps]
    positions = positions.tolist()
    chromosomes = chromosomes.tolist()
    macs = mafs.tolist()
    mafs = mafs.tolist()
    phen_vals_list = g['phen_values'][start_i:stop_i]
    gene_ids = g['gene_ids'][...].tolist()[start_i:stop_i]
    ets = g['ecotypes'][...].tolist()
    K = sp.matrix(g['kinship'][...])
    Z = sp.matrix(g['Z'][...])
    K = Z * K * Z.T
    correlations = g['correlations'][start_i:stop_i]
    p_herits = g['p_herits'][start_i:stop_i]
    C = g['env_vector'][...]

    print 'In total there are % d SNPs.' % len(snps)
    b_threshold = -sp.log10(1.0 / (len(snps) * 20.0))
    gene_dict = dp.parse_tair_gff_file()#_load_genes_list_('rna_seq_031311_ % sC' % temperature)

    if w_marginal_analysis:
        g10 = f['10C']
        g16 = f['16C']
        if debug_filter < 1.0:
            num_snps = len(g10['positions'])
            filter = random.sample(range(num_snps), int(round(num_snps * debug_filter)))
            filter.sort()
            snps10 = g10['snps'][...][filter]
            mafs10 = g10['mafs'][...][filter]
            macs10 = g10['macs'][...][filter]
            positions10 = g10['positions'][...][filter]
            chromosomes10 = g10['chromosomes'][...][filter]
            num_snps = len(g16['positions'])
            filter = random.sample(range(num_snps), int(round(num_snps * debug_filter)))
            filter.sort()
            snps16 = g16['snps'][...][filter]
            mafs16 = g16['mafs'][...][filter]
            macs16 = g16['macs'][...][filter]
            positions16 = g16['positions'][...][filter]
            chromosomes16 = g16['chromosomes'][...][filter]
        else:
            snps10 = g10['snps'][...]
            mafs10 = g10['mafs'][...]
            macs10 = g10['macs'][...]
            positions10 = g10['positions'][...]
            chromosomes10 = g10['chromosomes'][...]
            snps16 = g16['snps'][...]
            mafs16 = g16['mafs'][...]
            macs16 = g16['macs'][...]
            positions16 = g16['positions'][...]
            chromosomes16 = g16['chromosomes'][...]

        K10 = g10['kinship'][...]
        phen_vals_list10 = g10['phen_values']
        gene_ids10 = g10['gene_ids'][...].tolist()
        K16 = g16['kinship'][...]
        phen_vals_list16 = g16['phen_values']
        gene_ids16 = g16['gene_ids'][...].tolist()


    for i, gene_tair_id in enumerate(gene_ids):

        chrom = int(gene_tair_id[2])
        chr_prefix = ' % schr_ % d / rna_seq_ % s_both_pid % d_ % s' % (file_prefix, chrom, run_id, i + start_i, gene_tair_id)
        hist_prefix = ' % shistograms / rna_seq_ % s_both_pid % d_ % s' % (file_prefix, run_id, i + start_i, gene_tair_id)
        manhattan_prefix = ' % smanhattan_plots / rna_seq_ % s_both_pid % d_ % s' % (file_prefix, run_id, i + start_i, gene_tair_id)
        qq_prefix = ' % sqq_plots / rna_seq_ % s_both_pid % d_ % s' % (file_prefix, run_id, i + start_i, gene_tair_id)
        summary_file = chr_prefix + '_summary.pickled'
        #Check if the file is already there.
        print 'Checking % s' % (summary_file)
        if os.path.isfile(summary_file):
            print 'File already exists, skipping this file'
            continue
        print 'File not found, running analysis..'

        summary_dict = {}
        d = gene_dict[gene_tair_id]
        gene_strand = d['strand']
        try:
            chrom = int(d['chromosome'])
        except Exception:
            raise
        gene = gwaResults.Gene(chromosome=int(d['chromosome']), startPos=d['start_pos'],
                endPos=d['end_pos'], name=gene_tair_id, description=None, dbRef=gene_tair_id,
                tairID=gene_tair_id)
        print i + start_i, gene
        print 'Heritabilties: % 0.4f and % 0.4f' % (p_herits[i, 0], p_herits[i, 1])

        #Check if marginal analysis should be performed
        if w_marginal_analysis:
            i10 = bisect.bisect(gene_ids10, gene_tair_id) - 1
            r10 = lm.emmax(snps10, phen_vals_list10[i10], K10)
            r10 = gr.Result(scores=r10['ps'], macs=macs10, mafs=mafs10, positions=positions10,
                    chromosomes=chromosomes10, perc_var_expl=r10['var_perc'])
            i16 = bisect.bisect(gene_ids16, gene_tair_id) - 1
            r16 = lm.emmax(snps16, phen_vals_list16[i16], K16)
            r16 = gr.Result(scores=r16['ps'], macs=macs16, mafs=mafs16, positions=positions16,
                    chromosomes=chromosomes16, perc_var_expl=r16['var_perc'])

        #Standardizing trait
        y = sp.copy(phen_vals_list[i])
        y0 = y[(C == 0)[:, 0]]
        y1 = y[(C == 1)[:, 0]]
        mean0 = y0.mean(axis=0)
        mean1 = y1.mean(axis=0)
        y0 -= mean0
        y1 -= mean1
        std0 = y0.std(axis=0)
        std1 = y1.std(axis=0)
        y0 /= std0
        y1 /= std1
        y[(C == 0)[:, 0]] = y0
        y[(C == 1)[:, 0]] = y1
        y = sp.matrix(y).T
        E = sp.matrix(C)

        s1 = time.time()
        KV = fcm.fit_gXe_model(sp.array(y), sp.array(E), [sp.array(K)], corr=correlations[i],
                her=p_herits[i], Nk=3, Nd=50)
        secs = time.time() - s1
        if secs > 60:
            mins = int(secs) / 60
            secs = secs - mins * 60
            print 'Took % d mins and % f seconds to obtain variance components.' % (mins, secs)
        else:
            print 'Took % f seconds to obtain variance components.' % (secs)
        print KV.keys()
        print KV['p_herits'], KV['variances'], KV['LML'], KV['grid_lml'], KV['hyperparams_o']

        s1 = time.time()
        lmm = lm.LinearMixedModel(sp.array(y).flatten())
        H = sp.matrix(KV['K'])
        lmm.add_random_effect(H)
        lmm.add_factor(C.T)
        H_sqrt = lm.cholesky(H).T
        H_sqrt_inv = H_sqrt.I
        res = lmm._emmax_GxT_f_test_(snps, H_sqrt_inv, E, Z)
        secs = time.time() - s1
        if secs > 60:
            mins = int(secs) / 60
            secs = secs - mins * 60
            print 'Took % d mins and % f seconds to perform genome - wide association scan.' % (mins, secs)
        else:
            print 'Took % f seconds to perform genome - wide association scan.' % (secs)

        ex_g_res = gr.Result(scores=res['g_res']['ps'], macs=macs, mafs=mafs, positions=positions, chromosomes=chromosomes,
                perc_var_expl=res['g_res']['var_perc'], betas=res['g_res']['betas'])
        g_ks = agr.calc_ks_stats(res['g_res']['ps'])
        g_med = agr.calc_median(res['g_res']['ps'])
        ex_gt_res = gr.Result(scores=res['gt_res']['ps'], macs=macs, mafs=mafs, positions=positions, chromosomes=chromosomes,
                perc_var_expl=res['gt_res']['var_perc'], betas=res['gt_res']['betas'])

        gt_ks = agr.calc_ks_stats(res['gt_res']['ps'])
        gt_med = agr.calc_median(res['gt_res']['ps'])
        ex_gt_g_res = gr.Result(scores=res['gt_g_res']['ps'], macs=macs, mafs=mafs, positions=positions, chromosomes=chromosomes,
                perc_var_expl=res['gt_g_res']['var_perc'])
        gt_g_ks = agr.calc_ks_stats(res['gt_g_res']['ps'])
        gt_g_med = agr.calc_median(res['gt_g_res']['ps'])

        #Save statistics!!!
        summary_dict['g_test'] = {'KS':g_ks, 'med_pval_bias':g_med,
                        'analysis_dict': ex_g_res.get_gene_analysis(gene)}
        summary_dict['gt_test'] = {'KS':gt_ks, 'med_pval_bias':gt_med,
                        'analysis_dict': ex_gt_res.get_gene_analysis(gene)}
        summary_dict['gt_g_test'] = {'KS':gt_g_ks, 'med_pval_bias':gt_g_med,
                        'analysis_dict': ex_gt_g_res.get_gene_analysis(gene)}
        summary_dict['phen_info'] = {'correlations':correlations[i], 'pseudo_heritabilities':p_herits[i],
                        'SD':sp.array([std0, std1]), 'mean':sp.array([mean0, mean1])}
        #Heritatbilities, correlations, etc.
        summary_dict['var_components'] = {'sigma_g':KV['variances']['Sk'], 'sigma_e':KV['variances']['Se'],
                        'rho_g':KV['variances']['rho'],
                        'heritabilities':KV['p_herits'], 'mle':KV['LML']}


        quantiles = []
        log_quantiles = []
        labels = []
        ms = []
        l = [(ex_g_res, 'Common effects', 'g'), (ex_gt_res, 'Full model', 'gt'),
            (ex_gt_g_res, 'GxT', 'gt_g')]
        line_colors = ['g', 'r', 'b']
        if w_marginal_analysis:
            l.extend([(r10, '10C', '10C'), (r16, '16C', '16C')])
            line_colors.extend(['m', 'c'])
        for r, n, fl in l:
            s = r.get_scores()
            quantiles.append(agr.get_quantiles(s, 1000))
            log_quantiles.append(agr.get_log_quantiles(s, 1000, max_val=5.5))
            labels.append(n)
            r.neg_log_trans()
            ms.append(r.max_score())
            r.filter_percentile(0.98) #Filter lowest  98% 
            #Dump results into pickle files. 
            additional_columns = ['perc_var_expl']
            if fl in  ['g', 'gt']:
                additional_columns.append('betas')
            r.write_to_file(' % s_ % s_.pvals' % (chr_prefix, fl), only_pickled=True,
                    additional_columns=additional_columns)


        dq = {}
        dlq = {}
        for k_i, k in enumerate(['g', 'gt', 'gt_g']):
            dq[k] = quantiles[k_i]
            dlq[k] = log_quantiles[k_i]
        summary_dict['quantiles'] = dq
        dlq['max_val'] = 5.5
        summary_dict['log_quantiles'] = d

        #Dump summary_dict into a pickle file.
        cPickle.dump(summary_dict, open(summary_file, 'w'), protocol=2)


        max_scores = [r.max_score() for r in [ex_g_res, ex_gt_res, ex_gt_g_res]]
        #If interesting results, then plot stuff...
        if max(max_scores) >= 10:
            #Plot histogram
            y = phen_vals_list[i]
            y0 = y[(C == 0)[:, 0]]
            y1 = y[(C == 1)[:, 0]]
            pylab.figure(figsize=(6, 4))
            n, bins, patches = pylab.hist([y0, y1], 10, histtype='bar', alpha=0.7,
                        label=['10C', '16C'])
            pylab.legend()
            pylab.xlabel('Anscombe transformed RNAseq counts')
            hist_file = hist_prefix + '_histogram.pdf'
            pylab.savefig(hist_file, format='pdf')
            pylab.clf()

            agr.simple_qqplot(quantiles, pdf_file=qq_prefix + '_qq_plot.pdf' ,
                    quantile_labels=labels, line_colors=line_colors, plot_label='a')
            agr.simple_log_qqplot(log_quantiles, pdf_file=qq_prefix + '_log_qq_plot.pdf',
                    line_colors=line_colors, max_val=5.5, plot_label='b')

            chrom_col_map = {}
            for i in range(1, 6):
                chrom_col_map[i] = '#1199EE' if i % 2 == 0 else '#11BB00'

            if w_marginal_analysis:
                f = pylab.figure(figsize=(10, 10))
                ax1 = f.add_axes([0.05, 0.05 + (0.8 * 0.94), 0.945, (0.18 * 0.94) ])
                ax2 = f.add_axes([0.05, 0.05 + (0.6 * 0.94), 0.945, (0.18 * 0.94) ])
                ax3 = f.add_axes([0.05, 0.05 + (0.4 * 0.94), 0.945, (0.18 * 0.94) ])
                ax3.set_ylabel(r'$ - $log$($p - value$)$')
                ax4 = f.add_axes([0.05, 0.05 + (0.2 * 0.94), 0.945, (0.18 * 0.94) ])
                ax5 = f.add_axes([0.05, 0.05, 0.945, (0.18 * 0.94) ])
                for ax in [ax1, ax2, ax3, ax4, ax5]:
                    ax.spines['top'].set_visible(False)
                    ax.xaxis.set_ticks_position('bottom')
                    ax.xaxis.set_label_position('bottom')

                max_score = max(ms)
                for r, ax, lab in [(ex_gt_res, ax1, 'a'), (ex_g_res, ax2, 'b'),
                            (ex_gt_g_res, ax3, 'c'), (r10, ax4, 'd'), (r16, ax5, 'e')]:
                    wo_xtick_labels = lab != 'e'
                    r.plot_manhattan2(ax=ax, neg_log_transform=False, plot_bonferroni=True,
                            plot_xaxis=True, chrom_colormap=chrom_col_map,
                            sign_color='#DD1122', wo_xtick_labels=wo_xtick_labels,
                            max_score=max_score, cand_genes=[gene],
                            b_threshold=b_threshold)
                    x_min, x_max = ax.get_xlim()
                    x_range = x_max - x_min
                    y_min, y_max = ax.get_ylim()
                    y_range = y_max - y_min
                    ax.text(0.96 * x_range + x_min, 0.94 * y_range + y_min, lab)
            else:
                f = pylab.figure(figsize=(10, 6))
                ax1 = f.add_axes([0.05, 0.06 + (0.667 * 0.94), 0.945, (0.3 * 0.94) ])
                ax2 = f.add_axes([0.05, 0.06 + (0.333 * 0.94), 0.945, (0.3 * 0.94) ])
                ax2.set_ylabel(r'$ - $log$($p - value$)$')
                ax3 = f.add_axes([0.05, 0.06, 0.945, (0.3 * 0.94) ])
                for ax in [ax1, ax2, ax3]:
                    ax.spines['top'].set_visible(False)
                    ax.xaxis.set_ticks_position('bottom')
                    ax.xaxis.set_label_position('bottom')

                max_score = max(ms)
                for r, ax, lab in [(ex_gt_res, ax1, 'a'), (ex_g_res, ax2, 'b'),
                            (ex_gt_g_res, ax3, 'c')]:
                    wo_xtick_labels = lab != 'c'
                    r.plot_manhattan2(ax=ax, neg_log_transform=False, plot_bonferroni=True,
                            plot_xaxis=True, chrom_colormap=chrom_col_map,
                            sign_color='#DD1122', wo_xtick_labels=wo_xtick_labels,
                            max_score=max_score, cand_genes=[gene],
                            b_threshold=b_threshold)
                    x_min, x_max = ax.get_xlim()
                    x_range = x_max - x_min
                    y_min, y_max = ax.get_ylim()
                    y_range = y_max - y_min
                    ax.text(0.96 * x_range + x_min, 0.94 * y_range + y_min, lab)

            pylab.savefig(manhattan_prefix + '_EX_manhattan_all.png', dpi=400)




def run_gwas(file_prefix, start_i, stop_i, temperature, filter_threshold=0.02, call_method_id=79,
        data_format='diploid_int', debug_filter=1.0, test_run=False, conf_matrix_type='none',
        run_id='blurb'):
    """
    GWAS
    """
    f = h5py.File(env['rna_seq_data_dir'] + 'rna_seq_gwas_data_cm%d.h5py' % call_method_id)
    Kexpr = None
    if conf_matrix_type != 'none':
        if conf_matrix_type == 'k_expr_standard':
            conf_file_name = env['rna_seq_data_dir'] + 'standard/confounders_cm%d.hdf5' % (call_method_id)
            print 'Attempting to load %s' % conf_file_name
            conf_f = h5py.File(conf_file_name)
            Kexpr = conf_f[temperature]['Kexpr'][...]
        elif conf_matrix_type == 'k_expr_blup':
            conf_file_name = env['rna_seq_data_dir'] + 'blubPop/confounders_cm%d.hdf5' % (call_method_id)
            print 'Attempting to load %s' % conf_file_name
            conf_f = h5py.File(conf_file_name)
            Kexpr = conf_f[temperature]['Kexpr'][...]
        elif conf_matrix_type == 'k_joint_standard':
            conf_file_name = env['rna_seq_data_dir'] + 'standard/confounders_cm%d.hdf5' % (call_method_id)
            print 'Attempting to load %s' % conf_file_name
            conf_f = h5py.File(conf_file_name)
            K = conf_f[temperature]['K_joint'][...]
        elif conf_matrix_type == 'k_joint_blup':
            conf_file_name = env['rna_seq_data_dir'] + 'blubPop/confounders_cm%d.hdf5' % (call_method_id)
            print 'Attempting to load %s' % conf_file_name
            conf_f = h5py.File(conf_file_name)
            K = conf_f[temperature]['K_joint'][...]
    file_prefix = file_prefix + '%s/' % (conf_matrix_type)

    g = f[temperature]
    if debug_filter < 1.0:
        num_snps = len(g['positions'])
        filter = random.sample(range(num_snps), int(round(num_snps * debug_filter)))
        filter.sort()
        snps = g['snps'][...][filter]
        mafs = g['mafs'][...][filter]
        macs = g['macs'][...][filter]
        positions = g['positions'][...][filter]
        chromosomes = g['chromosomes'][...][filter]
    else:
        snps = g['snps'][...]
        mafs = g['mafs'][...]
        macs = g['macs'][...]
        positions = g['positions'][...]
        chromosomes = g['chromosomes'][...]
    snps = [snp for snp in snps]
    positions = positions.tolist()
    chromosomes = chromosomes.tolist()
    macs = mafs.tolist()
    mafs = mafs.tolist()
    phen_vals_list = g['phen_values'][start_i:stop_i]
    gene_ids = g['gene_ids'][...].tolist()[start_i:stop_i]
    ets = g['ecotypes'][...].tolist()
    if conf_matrix_type in ['none', 'k_expr_blup', 'k_expr_standard']:
        K = g['kinship'][...]
    if conf_matrix_type == 'none':
        sd = snpsdata.construct_snps_data_set(snps, positions, chromosomes, ets)
        filtered_K = kinship.scale_k(sd.get_ibs_kinship_matrix())

    print 'In total there are %d SNPs.' % len(snps)
    gene_dict = dp.parse_tair_gff_file()#_load_genes_list_('rna_seq_031311_%sC' % temperature)

    if test_run:
        heritabilities = []
        heritabilities_with_e = []
        perc_var_e = []
        med_pvals = []
        KS_stats = []
        num_cis_peaks = 0
        num_trans_peaks = 0

    for i, gene_tair_id in enumerate(gene_ids):

        chrom = int(gene_tair_id[2])
        chr_prefix = '%schr_%d/rna_seq_%s_both_pid%d_%s' % (file_prefix, chrom, run_id, i + start_i, gene_tair_id)
        hist_prefix = '%shistograms/rna_seq_%s_both_pid%d_%s' % (file_prefix, run_id, i + start_i, gene_tair_id)
        manhattan_prefix = '%smanhattan_plots/rna_seq_%s_both_pid%d_%s' % (file_prefix, run_id, i + start_i, gene_tair_id)
        summary_file = chr_prefix + '_summary.pickled'
        #Check if the file is already there.
        print 'Checking %s' % (summary_file)
        if os.path.isfile(summary_file):
            print 'File already exists, skipping this file'
            continue
        print 'File not found, running analysis..'

        phen_vals = phen_vals_list[i]
        summary_dict = {}
        d = gene_dict[gene_tair_id]
        gene_strand = d['strand']
        try:
            chrom = int(d['chromosome'])
        except Exception:
            raise
        gene = gwaResults.Gene(chromosome=int(d['chromosome']), startPos=d['start_pos'],
                endPos=d['end_pos'], name=gene_tair_id, description=None, dbRef=gene_tair_id,
                tairID=gene_tair_id)
        print i + start_i, gene

        print 'Calculating heritabilities, etc.'
        res = lm.get_emma_reml_estimates(phen_vals, K)
        pherit = res['pseudo_heritability']
        if conf_matrix_type in ['k_expr_standard', 'k_expr_blup']:
            res = lm.get_emma_reml_estimates(phen_vals, K, K2=Kexpr)
            gen_perc_var = res['perc_var1']
            expr_perc_var = res['perc_var2']
            gen_expr_perc_var = res['pseudo_heritability']
            print 'pseudo-heritability (only kinship): %0.4f, pseudo-heritability (with expr. het.): %0.4f, ' % (pherit, gen_perc_var)
            print expr_perc_var, gen_expr_perc_var
        else:
            print 'pseudo-heritability (only kinship): %0.4f' % pherit
        if test_run:
            heritabilities.append(pherit)
            if conf_matrix_type not in ['none', 'k_joint_standard', 'k_joint_blup']:
                heritabilities_with_e.append(gen_perc_var)
                perc_var_e.append(expr_perc_var)

        elif conf_matrix_type == 'none':
            print 'Applying LM'
            res = lm.linear_model(snps, phen_vals)
            pvals = res['ps'].tolist()
            perc_var_expl = res['var_perc'].tolist()
            lm_res = gr.Result(scores=pvals, macs=macs, mafs=mafs, positions=positions, chromosomes=chromosomes,
                    perc_var_expl=perc_var_expl)
            print 'Summarizing LM'
            summary_dict['LM'] = lm_res.get_gene_analysis(gene)
            summary_dict['LM']['kolmogorov_smirnov'] = agr.calc_ks_stats(res['ps'])
            summary_dict['LM']['pval_median'] = agr.calc_median(res['ps'])

        print 'Applying EX Stepwise'
        ex_sw_res = lm.emmax_step_wise(phen_vals, K, macs=macs, mafs=mafs, positions=positions,
                    chromosomes=chromosomes, snps=snps, num_steps=5, cand_gene_list=[gene],
                    with_qq_plots=False, log_qq_max_val=6.0,
                    pval_file_prefix=chr_prefix + '_ex', K2=Kexpr)
        print 'Summarizing the step-wise mixed model'
        pvals = ex_sw_res['first_emmax_res']['ps'].tolist()
        perc_var_expl = ex_sw_res['first_emmax_res']['var_perc'].tolist()
        ex_res = gr.Result(scores=pvals, macs=macs, mafs=mafs, positions=positions, chromosomes=chromosomes,
                perc_var_expl=perc_var_expl)
        summary_dict['EX'] = ex_res.get_gene_analysis(gene)
        summary_dict['pseudo_heritability'] = ex_sw_res['step_info_list'][0]['pseudo_heritability']
        summary_dict['EX']['kolmogorov_smirnov'] = agr.calc_ks_stats(ex_sw_res['first_emmax_res']['ps'])
        summary_dict['EX']['pval_median'] = agr.calc_median(ex_sw_res['first_emmax_res']['ps'])
        if test_run:
            med_pvals.append(summary_dict['EX']['pval_median'])
            KS_stats.append(summary_dict['EX']['kolmogorov_smirnov']['D'])

        #Does the linear mixed model fit the data better?
        summary_dict['MM_LRT'] = lm.mm_lrt_test(phen_vals, K)

        #FINISH summarizing the stepwise!!!
        summarize_stepwise(summary_dict, gene, ex_sw_res['step_info_list'], ex_sw_res['opt_dict'])


        if test_run:
            cof_gene_dists = summary_dict['SW']['mbonf']['cof_gene_dist']
            for cof_gene_dist in cof_gene_dists:
                if cof_gene_dist == -1 or cof_gene_dist > 100000:
                    num_trans_peaks += 1
                else:
                    num_cis_peaks += 1


        #Run local vs. global kinship analysis only if using traditional kinship matrix.
        if conf_matrix_type == 'none':
            cvt_dict = {'radius':{}, 'tss_upstream':{}}
            print 'Comparing cis vs. trans kinship'

            #Check 1 mb, 200kb, 100kb, 50kb, 20kb, 10kb, 2kb, 0kb
            for radius in [500000, 100000, 50000, 25000, 10000, 5000, 1000, 0]:
                print radius
                r_start_pos = max(gene.startPos - radius, 0)
                r_end_pos = gene.endPos + radius
                d = sd.get_region_split_kinships([(chrom, r_start_pos, r_end_pos)],
                                kinship_method='ibs', global_kinship=filtered_K)
                reg_k = d['regional_k']
                glob_k = d['global_k']
                if reg_k != None:
                    cvt_dict['radius'][radius] = lm.local_vs_global_mm(phen_vals, reg_k, glob_k, filtered_K)
                else:
                    cvt_dict['radius'][radius] = None
                print cvt_dict['radius'][radius]

            #Check TSS, 100kb, 50kb,25kb, 10kb,5kb,0kb, (all upstream)
            for dist in [200000, 100000, 50000, 25000, 10000, 5000, 1000]:
                print dist, gene_strand
                if gene_strand == '+':
                    r_start_pos = max(gene.startPos - dist, 0)
                    r_end_pos = gene.startPos
                else:
                    r_start_pos = gene.endPos
                    r_end_pos = gene.endPos + dist
                d = sd.get_region_split_kinships([(chrom, r_start_pos, r_end_pos)],
                                kinship_method='ibs', global_kinship=filtered_K)
                reg_k = d['regional_k']
                glob_k = d['global_k']
                if reg_k != None:
                    cvt_dict['tss_upstream'][dist] = lm.local_vs_global_mm(phen_vals, reg_k, glob_k, filtered_K)
                else:
                    cvt_dict['tss_upstream'][dist] = None
                print cvt_dict['tss_upstream'][dist]

            summary_dict['CVT'] = cvt_dict

        #Write info to file..
        cPickle.dump(summary_dict, open(summary_file, 'w'), protocol=2)

        f_prefix = hist_prefix + '_hist'
        phed = pd.phenotype_data(phen_dict={1:{'name':gene_tair_id, 'ecotypes':ets,
                                'values':phen_vals, 'transformation':None,
                                'raw_values':phen_vals}})
        phed.plot_histogram(1, title='Gene expressions for %s' % gene_tair_id,
                png_file=f_prefix + '.png', p_her=summary_dict['pseudo_heritability'],
                x_label='RNA seq expression levels (Anscombe transformed)')
        #Plot GWAs...
        res_list = [(ex_res, 'EX')]
        if conf_matrix_type == 'none':
            res_list.append((lm_res, 'LM'))
        for res, method_name in res_list:
            res.filter_percentile(filter_threshold, reversed=True)
            res.write_to_file('%s_%s_.pvals' % (chr_prefix, method_name), only_pickled=True,
                        additional_columns=['perc_var_expl'])
            if ex_res.min_score() < 10e-10:
                #print [cg.tairID for cg in cgs]
                f_prefix = '%s_%s_manhattan' % (manhattan_prefix, method_name)
                res.plot_manhattan(png_file=f_prefix + '.png', percentile=0, cand_genes=[gene],
                        plot_bonferroni=True, neg_log_transform=True)

        if test_run:
            print 'Cis peaks per gene:%0.4f , Trans peaks per gene:%0.4f' % \
                (num_cis_peaks / float(i + 1), num_trans_peaks / float(i + 1))


    if test_run:
        #plot stuff
        num_genes = len(gene_ids)
        print 'Cis peaks per gene:%0.4f , Trans peaks per gene:%0.4f' % (num_cis_peaks / float(num_genes), num_trans_peaks / float(num_genes))
        pylab.hist(heritabilities, bins=20, alpha=0.6, range=(0, 1))
        pylab.xlabel('herit. (wo expr. het. matrix)')
        pylab.savefig(file_prefix + '_herit_hist.png')
        pylab.figure()
        pylab.hist(heritabilities_with_e, bins=20, alpha=0.6, range=(0, 1))
        pylab.xlabel('herit. (w expr. het. matrix)')
        pylab.savefig(file_prefix + '_herit_with_expr_hist.png')
        pylab.figure()
        pylab.hist(perc_var_e, bins=20, alpha=0.6, range=(0, 1))
        pylab.xlabel('expr. het. perc. of var.')
        pylab.savefig(file_prefix + '_perc_var_expr_hist.png')
        pylab.figure()
        pylab.boxplot([med_pvals, KS_stats])
        pylab.axhline(0, color='k', alpha=0.6, ls='-.')
        pylab.xticks(range(1, 3), ['med. p-val', 'KS'])
        pylab.ylabel('Distance (p-value)')
        pylab.savefig(file_prefix + '_pval_dist_stats.png')
        pylab.figure()
        pylab.plot(heritabilities, heritabilities_with_e, marker='.', color='k', alpha=0.7, ls='None')
        pylab.xlabel('herit. (wo expr. het. matrix)')
        pylab.ylabel('herit. (w expr. het. matrix)')
        pylab.axis([-0.02, 1.02, -0.02, 1.02])
        pylab.savefig(file_prefix + '_herits_scatter.png')
        pylab.figure()
        pylab.plot(perc_var_e, heritabilities_with_e, marker='.', color='k', alpha=0.7, ls='None')
        pylab.xlabel('expr. het. perc. of var.')
        pylab.ylabel('herit. (w expr. het. matrix)')
        pylab.axis([-0.02, 1.02, -0.02, 1.02])
        pylab.savefig(file_prefix + '_herits_expr_perc_var_e_scatter.png')
        pylab.figure()
        pylab.plot(perc_var_e, heritabilities, marker='.', color='k', alpha=0.7, ls='None')
        pylab.xlabel('expr. het. perc. of var.')
        pylab.ylabel('herit. (wo expr. het. matrix)')
        pylab.axis([-0.02, 1.02, -0.02, 1.02])
        pylab.savefig(file_prefix + '_herits_perc_var_e_scatter.png')
        pylab.figure()
        pylab.plot(perc_var_e, med_pvals, marker='.', color='k', alpha=0.7, ls='None')
        pylab.xlabel('expr. het. perc. of var.')
        pylab.ylabel('median p-val bias')
        pylab.axis([-0.02, 1.02, -0.02, 1.02])
        pylab.savefig(file_prefix + '_med_pval_perc_var_e_scatter.png')



def plot(temperature=10, call_method_id=75, mapping_method='EX', mac_threshold=15, min_score=5,
        near_const_filter=20, data_format='binary', plot_data=True):
    #Load in chromosome dict..

    #file_prefix = '/srv/lab/data/rna_seq_062911/%dC/cm_%d/' % (temperature, call_method_id)
    file_prefix = '/srv/lab/data/rna_seq_083011/%dC/cm_%d/' % (temperature, call_method_id)

    results_dict_file = '%sresults_%s_mac%d.pickled' % (file_prefix, mapping_method, mac_threshold)
    res_dict = cPickle.load(open(results_dict_file))

    phen_file = '%s_%dC.csv' % (phen_file_prefix, temperature)
    phed = pd.parse_phenotype_file(phen_file, with_db_ids=False)  #load phenotype file
    phed.filter_near_const_phens(near_const_filter)
    phed.convert_to_averages()
    num_traits = phed.num_traits()
    pids = phed.phen_ids
    sd = dp.load_snps_call_method(call_method_id=call_method_id, data_format=data_format, debug_filter=0.01)
    indices_to_keep = sd.coordinate_w_phenotype_data(phed, 1, coord_phen=False)  #All phenotypes are ordered the same way, so we pick the first one.
    phed.filter_ecotypes(indices_to_keep, pids=pids)

    chrom_dict = {}
    for x_chrom in [1, 2, 3, 4, 5]:
        for y_chrom in [1, 2, 3, 4, 5]:
            chrom_dict[(x_chrom, y_chrom)] = {'scores':[], 'x_positions':[], 'y_positions':[],
                                'tair_ids':[], 'r2':[], 'mac':[]}
    scores = []
    for x_chrom, x_pos in res_dict:
        d = res_dict[(x_chrom, x_pos)]
        tair_id = d['tair_id']
        for y_chrom in [1, 2, 3, 4, 5]:
            cps_d = d['chrom_pos_score'][y_chrom]
            for i in range(len(cps_d['scores'])):
                s = cps_d['scores'][i]
                if s > min_score:
                    if s > 25:
                        s = 25
                    scores.append(s)
                    chrom_dict[(x_chrom, y_chrom)]['scores'].append(s)
                    chrom_dict[(x_chrom, y_chrom)]['tair_ids'].append(tair_id)
                    chrom_dict[(x_chrom, y_chrom)]['x_positions'].append(x_pos)
                    chrom_dict[(x_chrom, y_chrom)]['y_positions'].append(cps_d['positions'][i])

    #Write chrom_dict to file..
    if not plot_data:
        for x_chrom in [1, 2, 3, 4, 5]:
            for y_chrom in [1, 2, 3, 4, 5]:
                file_name = file_prefix + 'result_plots/pvalues_chrom%d_chrom%d_%s_min%d.txt' % (x_chrom, y_chrom, mapping_method, min_score)
                print 'Writing to file:', file_name
                with open(file_name, 'w') as f:
                    d = chrom_dict[(x_chrom, y_chrom)]
                    f.write('x_position, y_position, score, tair_id\n')
                    l = zip(d['x_positions'], d['y_positions'], d['scores'], d['tair_ids'])
                    l.sort()
                    for t in l:
                        f.write('%d,%d,%f,%s\n' % t)



    chrom_sizes = [30425061, 19694800, 23456476, 18578714, 26974904]
    cum_chrom_sizes = [sum(chrom_sizes[:i]) for i in range(5)]
    tot_num_bases = float(sum(chrom_sizes))
    rel_chrom_sizes = map(lambda x: 0.925 * (x / tot_num_bases), chrom_sizes)
    rel_cum_chrom_sizes = map(lambda x: 0.925 * (x / tot_num_bases), cum_chrom_sizes)
    for i in range(5):
        rel_cum_chrom_sizes[i] = rel_cum_chrom_sizes[i] + 0.02 + 0.01 * i

    chromosome_ends = {1:30.425061, 2:19.694800, 3:23.456476, 4:18.578714, 5:26.974904}
    print rel_chrom_sizes, rel_cum_chrom_sizes

    #Filter data..

    #Now plot data!!
    if plot_data:
        alpha = 0.8
        linewidths = 0
        vmin = min_score
        f = pylab.figure(figsize=(40, 35))
        chromosomes = [1, 2, 3, 4, 5]
        plot_file_name = file_prefix + 'result_plots/pvalues_%s_min%d.png' % (mapping_method, min_score)
        label = '$-log_{10}$(p-value)'
        vmax = max(scores)

        for yi, chr2 in enumerate(chromosomes):
            for xi, chr1 in enumerate(chromosomes):

                l = chrom_dict[(chr1, chr2)]['scores']
                if len(l) == 0:
                    continue
                ax = f.add_axes([0.96 * (rel_cum_chrom_sizes[xi] + 0.01), rel_cum_chrom_sizes[yi] - 0.02,
                        0.96 * (rel_chrom_sizes[xi]), rel_chrom_sizes[yi] ])
                ax.spines['right'].set_visible(False)
                ax.spines['bottom'].set_visible(False)
                #ax.tick_params(fontsize='x-large')
                if xi > 0:
                    ax.spines['left'].set_visible(False)
                    ax.yaxis.set_visible(False)
                else:
                    ax.yaxis.set_ticks_position('left')
                    ax.set_ylabel('Chromosome %d (Mb)' % chr2, fontsize='x-large')
                if yi < 4:
                    ax.spines['top'].set_visible(False)
                    ax.xaxis.set_visible(False)
                else:
                    ax.xaxis.set_ticks_position('top')
                    ax.xaxis.set_label_position('top')
                    ax.set_xlabel('Chromosome %d (Mb)' % chr1, fontsize='x-large')
                    #ax.set_xlabel('Chromosome %d' % chr1)

                #l = -sp.log10(l)
                #l = l.tolist()
                l_zxy = zip(l, chrom_dict[(chr1, chr2)]['x_positions'],
                    chrom_dict[(chr1, chr2)]['y_positions'])
                l_zxy.sort()
                l = map(list, zip(*l_zxy))
                zs = l[0]
                xs = map(lambda x: x / 1000000.0, l[1])
                ys = map(lambda x: x / 1000000.0, l[2])

                scatter_plot = ax.scatter(xs, ys, c=zs, alpha=alpha, linewidths=linewidths, vmin=vmin,
                            vmax=vmax)
                ax.axis([-0.025 * chromosome_ends[chr1], 1.025 * chromosome_ends[chr1],
                    - 0.025 * chromosome_ends[chr2], 1.025 * chromosome_ends[chr2]])

        cax = f.add_axes([0.965, 0.7, 0.01, 0.2])
        cb = pylab.colorbar(scatter_plot, cax=cax)
        cb.set_label(label, fontsize='xx-large')
        #cb.set_tick_params(fontsize='x-large')
        f.text(0.005, 0.47, 'Associated SNP position', size='xx-large', rotation='vertical')
        f.text(0.47, 0.988, 'Expressed gene position', size='xx-large')
        print 'Saving figure:', plot_file_name
        f.savefig(plot_file_name, format='png')








def load_and_plot_info_files(call_method_id=75, temperature=10, mac_threshold=15, debug_filter=1,
            near_const_filter=20, data_format='binary'):
    import random

    phen_file = '%s_%dC.csv' % (phen_file_prefix, temperature)
    phed = pd.parse_phenotype_file(phen_file, with_db_ids=False)  #load phenotype file
    phed.filter_near_const_phens(near_const_filter)
    phed.convert_to_averages()
    num_traits = phed.num_traits()
    pids = phed.phen_ids
    sd = dp.load_snps_call_method(call_method_id=call_method_id, data_format=data_format, debug_filter=0.01)
    indices_to_keep = sd.coordinate_w_phenotype_data(phed, 1, coord_phen=False)  #All phenotypes are ordered the same way, so we pick the first one.
    phed.filter_ecotypes(indices_to_keep, pids=pids)

    print 'Loading the gene annotation dictionary'
    gene_dict = dp.parse_tair_gff_file()
    run_id = 'd081511'
    #run_id = 'rs_%d' % call_method_id


    file_prefix = '/srv/lab/data/rna_seq_083011/%dC/cm_%d/' % (temperature, call_method_id)


    num_genes = 0

    radii = [500000, 100000, 50000, 25000, 10000, 5000, 1000, 0]
    tss_dists = [200000, 100000, 50000, 25000, 10000, 5000, 1000]
    cvt_summary_dict = {'radius':{'avg_cis_trans_var_ratio':[0.0 for r in radii],
                    'avg_cis_herit':[0.0 for r in radii],
                    'avg_trans_herit':[0.0 for r in radii],
                    'counts':[0.0 for td in radii]},
                'radius_herit':{'avg_cis_trans_var_ratio':[0.0 for r in radii],
                    'avg_cis_herit':[0.0 for r in radii],
                    'avg_trans_herit':[0.0 for r in radii],
                    'counts':[0.0 for td in radii]},
                'tss_dist':{'avg_cis_trans_var_ratio':[0.0 for td in tss_dists],
                    'avg_cis_herit':[0.0 for td in tss_dists],
                    'avg_trans_herit':[0.0 for td in tss_dists],
                    'counts':[0.0 for td in tss_dists]}}

    heritabilities = []
    transformations = []
    shapiro_wilk_pvals = []
    tair_ids = []
    pval_infl_dict = {}
    dist_min_pval_dict = {}
    distance_bins = [(0, 5000), (0, 10000), (0, 25000), (0, 50000), (0, 100000), (1, -1), (6, -1)]
    radius_bins = [0, 1000, 5000, 10000, 25000, 50000, 100000]
    bonf_sign_bin_dict = {}
    res_dict = {}
    sign_count = {}
    for mm in ['EX', 'LM', 'KW']:
        pval_infl_dict[mm] = {'kolmogorov_smirnov':[], 'median_pvals':[]}
        dist_min_pval_dict[mm] = {}
        for bin in distance_bins:
            dist_min_pval_dict[mm][bin] = 0
        bonf_sign_bin_dict[mm] = {}
        for bin in radius_bins:
            bonf_sign_bin_dict[mm][bin] = {'count':0.0, 'total':0.0}
        sign_count[mm] = 0

    cofactor_count_dict = {}
    for criteria in ['ebics', 'mbonf', 'min_cof_ppa']:
        cofactor_count_dict[criteria] = {'num_cofactor_list':[], 'bin_counts':sp.zeros(9),
                        'num_cis_cofactor_list':[], 'num_found':0}

    pickle_file_dict = {}
    for mm in ['EX', 'LM', 'KW']:
        pickle_file_dict[mm] = {}
        pickle_file_dict[mm]['file_name'] = '%sresults_%s_mac%d.pickled' % (file_prefix, mm, mac_threshold)
        pickle_file_dict[mm]['res_dict'] = {}

    pids = phed.get_pids()
    for i, pid in enumerate(pids):
        tair_id = phed.get_name(pid)
        chrom = int(tair_id[2])
        curr_file_prefix = '%schr_%d/rna_seq_%s_%dC_mac%d_pid%d_%s' % \
                    (file_prefix, chrom, run_id, temperature, mac_threshold, pid, tair_id)
        info_file_name = '%s_info.pickled' % curr_file_prefix
        for mm in ['EX', 'LM', 'KW']:
            res_dict[mm] = '%s_%s_.pvals' % (curr_file_prefix, mm)
        if random.random() > debug_filter:
            continue
        if os.path.isfile(info_file_name) and os.path.isfile(res_dict['EX'] + ".pickled") \
                and os.path.isfile(res_dict['LM'] + ".pickled") and os.path.isfile(res_dict['KW'] + ".pickled"):
            print 'Loading info file: %s' % info_file_name
            num_genes += 1
            info_dict = cPickle.load(open(info_file_name)) #Loading the info dict
            for mm in ['EX', 'LM', 'KW']:
                res_dict[mm] = gr.Result(res_dict[mm]) #Loading the result

            #Saving some basic statistics
            transformations.append(info_dict['transformation_type'])
            shapiro_wilk_pvals.append(info_dict['transformation_shapiro_pval'])
            heritabilities.append(info_dict['pseudo_heritability'])

            #cis vs. trans stuff
            cvt_dict = info_dict['CVT']
            for r_i, r in enumerate(radii):
                if cvt_dict['radius'][r] != None:
                    pvg = cvt_dict['radius'][r]['perc_var1']
                    pvl = cvt_dict['radius'][r]['perc_var2']
                    herit = cvt_dict['radius'][r]['pseudo_heritability1']
                    cvt_summary_dict['radius']['avg_cis_trans_var_ratio'][r_i] += pvl / (pvl + pvg)
                    cvt_summary_dict['radius']['avg_cis_herit'][r_i] += pvl * herit
                    cvt_summary_dict['radius']['avg_trans_herit'][r_i] += pvg * herit
                    cvt_summary_dict['radius']['counts'][r_i] += 1.0

            for r_i, r in enumerate(radii):
                if cvt_dict['radius'][r] != None:
                    herit = cvt_dict['radius'][r]['pseudo_heritability1']
                    if herit > 0.05:
                        pvg = cvt_dict['radius'][r]['perc_var1']
                        pvl = cvt_dict['radius'][r]['perc_var2']
                        cvt_summary_dict['radius_herit']['avg_cis_trans_var_ratio'][r_i] += pvl / (pvl + pvg)
                        cvt_summary_dict['radius_herit']['avg_cis_herit'][r_i] += pvl * herit
                        cvt_summary_dict['radius_herit']['avg_trans_herit'][r_i] += pvg * herit
                        cvt_summary_dict['radius_herit']['counts'][r_i] += 1.0

            for td_i, td in enumerate(tss_dists):
                if cvt_dict['tss_upstream'][td] != None:
                    pvg = cvt_dict['tss_upstream'][td]['perc_var1']
                    pvl = cvt_dict['tss_upstream'][td]['perc_var2']
                    herit = cvt_dict['tss_upstream'][td]['pseudo_heritability1']
                    cvt_summary_dict['tss_dist']['avg_cis_trans_var_ratio'][td_i] += pvl / (pvl + pvg)
                    cvt_summary_dict['tss_dist']['avg_cis_herit'][td_i] += pvl * herit
                    cvt_summary_dict['tss_dist']['avg_trans_herit'][td_i] += pvg * herit
                    cvt_summary_dict['tss_dist']['counts'][td_i] += 1.0



            tair_ids.append(tair_id)
            for mm in ['EX', 'LM', 'KW']:
                pval_infl_dict[mm]['kolmogorov_smirnov'].append(info_dict[mm]['kolmogorov_smirnov']['D'])
                pval_infl_dict[mm]['median_pvals'].append(info_dict[mm]['pval_median'])
                dist_min_pval = tuple(info_dict[mm]['dist_to_min_pval'])
                if res_dict[mm].min_score() < 1 / (20.0 * res_dict[mm].num_scores()):
                    sign_count[mm] += 1
                    for bin in distance_bins:
                        if dist_min_pval <= bin:
                            dist_min_pval_dict[mm][bin] += 1
                            break

                for bin in radius_bins:
                    pval = info_dict[mm]['bin_dict'][bin]['min_pval']
                    num_snps = info_dict[mm]['bin_dict'][bin]['num_snps']
                    if num_snps > 0:
                        bonf_sign_bin_dict[mm][bin]['total'] += 1
                        if pval < 1.0 / (20 * num_snps):
                            bonf_sign_bin_dict[mm][bin]['count'] += 1

            #Stepwise stuff 
            for criteria in ['ebics', 'mbonf', 'min_cof_ppa']:
                num_cofactors = len(info_dict['SW'][criteria]['cofactors'])
                cofactor_count_dict[criteria]['num_cofactor_list'].append(num_cofactors)
                if num_cofactors > 0:
                    cofactor_count_dict[criteria]['num_found'] += 1
                    cofactor_count_dict[criteria]['bin_counts'] += sp.array(info_dict['SW'][criteria]['bin_counts'])
                    cofactor_count_dict[criteria]['num_cis_cofactor_list'].append(info_dict['SW'][criteria]['bin_counts'][2])


            #Pre-process the results..
            for mm in ['EX', 'LM', 'KW']:
                res = res_dict[mm]
                #Trim results
                res.neg_log_trans()
                if mm == 'EX':
                    res.filter_attr('scores', 3) #Filter everything below 10^-2.5
                else:
                    res.filter_attr('scores', 4) #Filter everything below 10^-4
                if res.num_scores() == 0:
                    print "Skipping file since nothing is below 10^-5"
                    continue

                gene_d = gene_dict[tair_id]
                avg_g_pos = (gene_d['start_pos'] + gene_d['end_pos']) / 2.0
                chrom = int(gene_d['chromosome']) #Current gene chromosome


                #Prepare for plotting results.. x,y style, where gene is x, and y is p-values
                chrom_pos_score_dict = res.get_chrom_score_pos_dict()

                dist_dict = {}
                for score_threshold in [5, 6, 7]: #negative log10 thresholds.
                    if len(res.snp_results['scores']) == 0:
                        dist_dict[score_threshold] = -2 #No results
                    else:
                        res.filter_attr('scores', score_threshold)
                        if len(res.snp_results['scores']) == 0:
                            dist_dict[score_threshold] = -2 #No results
                        else:
                            cps_dict = res.get_chrom_score_pos_dict()
                            pos_list = cps_dict[chrom]['positions']
                            if len(pos_list) > 0:
                                distances = sp.absolute(sp.array(pos_list) - avg_g_pos)
                                d_i = sp.argmin(distances)
                                dist_dict[score_threshold] = distances[d_i] #Min distance.
                            else:
                                dist_dict[score_threshold] = -1 #Different chromosome

                pickle_file_dict[mm]['res_dict'][(chrom, avg_g_pos)] = {'tair_id':tair_id,
                            'chrom_pos_score':chrom_pos_score_dict, 'dist_dict':dist_dict,
                            'pid':pid}
                print dist_dict
        else:
            print "Didn't find file: %s or %s" % (info_file_name, res_dict['EX'] + ".pickled")

    for mm in ['EX', 'LM', 'KW']:
        cPickle.dump(pickle_file_dict[mm]['res_dict'], open(pickle_file_dict[mm]['file_name'], 'wb'), protocol=2)


    for r_i, r in enumerate(radii):
        r_counts = cvt_summary_dict['radius']['counts'][r_i]
        cvt_summary_dict['radius']['avg_cis_trans_var_ratio'][r_i] = \
            cvt_summary_dict['radius']['avg_cis_trans_var_ratio'][r_i] / r_counts
        cvt_summary_dict['radius']['avg_cis_herit'][r_i] = \
            cvt_summary_dict['radius']['avg_cis_herit'][r_i] / r_counts
        cvt_summary_dict['radius']['avg_trans_herit'][r_i] = \
            cvt_summary_dict['radius']['avg_trans_herit'][r_i] / r_counts


    for r_i, r in enumerate(radii):
        r_counts = cvt_summary_dict['radius_herit']['counts'][r_i]
        cvt_summary_dict['radius_herit']['avg_cis_trans_var_ratio'][r_i] = \
            cvt_summary_dict['radius_herit']['avg_cis_trans_var_ratio'][r_i] / r_counts
        cvt_summary_dict['radius_herit']['avg_cis_herit'][r_i] = \
            cvt_summary_dict['radius_herit']['avg_cis_herit'][r_i] / r_counts
        cvt_summary_dict['radius_herit']['avg_trans_herit'][r_i] = \
            cvt_summary_dict['radius_herit']['avg_trans_herit'][r_i] / r_counts


    for td_i, td in enumerate(tss_dists):
        td_counts = cvt_summary_dict['tss_dist']['counts'][td_i]
        cvt_summary_dict['tss_dist']['avg_cis_trans_var_ratio'][td_i] = \
            cvt_summary_dict['tss_dist']['avg_cis_trans_var_ratio'][td_i] / td_counts
        cvt_summary_dict['tss_dist']['avg_cis_herit'][td_i] = \
            cvt_summary_dict['tss_dist']['avg_cis_herit'][td_i] / td_counts
        cvt_summary_dict['tss_dist']['avg_trans_herit'][td_i] = \
            cvt_summary_dict['tss_dist']['avg_trans_herit'][td_i] / td_counts



    results_prefix = env['results_dir'] + 'RNAseq_summary_%dC_cm%d' % (temperature, call_method_id)

    pylab.figure()
    pylab.plot(cvt_summary_dict['radius']['avg_cis_trans_var_ratio'])
    pylab.ylabel('Avg. perc. of cis genetic var.')
    pylab.xlabel('Dist. from gene (kb)')
    pylab.xticks(range(8), [500, 100, 50, 25, 10, 5, 1, 0])
    pylab.savefig(results_prefix + '_avg_perc_cis_gen_var_rad.png')
    pylab.clf()

    pylab.figure()
    pylab.plot(cvt_summary_dict['tss_dist']['avg_cis_trans_var_ratio'])
    pylab.ylabel('Avg. perc. of cis genetic var.')
    pylab.xlabel('Dist. upstream from gene TSS (kb)')
    pylab.xticks(range(7), [200, 100, 50, 25, 10, 5, 1])
    pylab.savefig(results_prefix + '_avg_perc_cis_gen_var_td.png')
    pylab.clf()

#    pylab.figure()
#    pylab.plot(cvt_summary_dict['tss_dist']['avg_cis_herit'])
#    pylab.ylabel('Avg. cis heritability')
#    pylab.xlabel('Dist. upstream from gene TSS (kb)')
#    pylab.xticks(range(7), [200, 100, 50, 25, 10, 5, 1])
#    pylab.savefig(results_prefix + 'avg_cis_herit_td.png')
#    pylab.clf()
#
#
#    pylab.figure()
#    pylab.plot(cvt_summary_dict['tss_dist']['avg_trans_herit'])
#    pylab.ylabel('Avg. remaining heritability')
#    pylab.xlabel('Dist. upstream from gene TSS (kb)')
#    pylab.xticks(range(7), [200, 100, 50, 25, 10, 5, 1])
#    pylab.savefig(results_prefix + 'avg_trans_herit_td.png')
#    pylab.clf()


#    pylab.figure()
#    pylab.plot(cvt_summary_dict['radius']['avg_trans_herit'])
#    pylab.ylabel('Avg. remaining heritability')
#    pylab.xlabel('Dist. from gene (kb)')
#    pylab.xticks(range(8), [500, 100, 50, 25, 10, 5, 1, 0])
#    pylab.savefig(results_prefix + 'avg_trans_herit_rad.png')
#    pylab.clf()
#
#    pylab.figure()
#    pylab.plot(cvt_summary_dict['radius']['avg_cis_herit'])
#    pylab.ylabel('Avg. cis heritability')
#    pylab.xlabel('Dist. from gene (kb)')
#    pylab.xticks(range(8), [500, 100, 50, 25, 10, 5, 1, 0])
#    pylab.savefig(results_prefix + 'avg_cis_herit_rad.png')
#    pylab.clf()

    tot_herit = sp.array(cvt_summary_dict['radius']['avg_cis_herit']) + \
        sp.array(cvt_summary_dict['radius']['avg_trans_herit'])
    cis_herit = sp.array(cvt_summary_dict['radius']['avg_cis_herit'])
    pylab.figure(figsize=(10, 6))
    pylab.axes([0.06, 0.08, 0.92, 0.90])
    pylab.fill_between([0, 7], 0, 1, color='#DD3333', alpha=0.8, label='Error')
    pylab.fill_between(sp.arange(8), 0, tot_herit, color='#22CC44', alpha=0.8, label='Heritable variance')
    pylab.fill_between(sp.arange(8), 0, cis_herit, color='#2255AA', \
                alpha=0.8, label='Heritable variance (cis)')
    pylab.ylabel('Average partition of variance')
    pylab.xlabel('Dist. from gene (kb)')
    pylab.xticks(range(8), [500, 100, 50, 25, 10, 5, 1, 0])
    pylab.legend(loc=1, ncol=3, shadow=True)
    pylab.axis([0, 7, 0, 1])
    pylab.savefig(results_prefix + 'avg_herit_rad.png')

    tot_herit = sp.array(cvt_summary_dict['radius_herit']['avg_cis_herit']) + \
        sp.array(cvt_summary_dict['radius_herit']['avg_trans_herit'])
    cis_herit = sp.array(cvt_summary_dict['radius_herit']['avg_cis_herit'])
    pylab.figure(figsize=(10, 6))
    pylab.axes([0.06, 0.08, 0.92, 0.90])
    pylab.fill_between([0, 7], 0, 1, color='#DD3333', alpha=0.8, label='Error')
    pylab.fill_between(sp.arange(8), 0, tot_herit, color='#22CC44', alpha=0.8, label='Heritable variance')
    pylab.fill_between(sp.arange(8), 0, cis_herit, color='#2255AA', \
                alpha=0.8, label='Heritable variance (cis)')
    pylab.ylabel('Average partition of variance')
    pylab.xlabel('Dist. from gene (kb)')
    pylab.xticks(range(8), [500, 100, 50, 25, 10, 5, 1, 0])
    pylab.legend(loc=1, ncol=3, shadow=True)
    pylab.axis([0, 7, 0, 1])
    pylab.savefig(results_prefix + 'avg_herit_2_rad.png')



    tot_herit = sp.array(cvt_summary_dict['tss_dist']['avg_cis_herit']) + \
        sp.array(cvt_summary_dict['tss_dist']['avg_trans_herit'])
    cis_herit = sp.array(cvt_summary_dict['tss_dist']['avg_cis_herit'])
    pylab.figure(figsize=(10, 6))
    pylab.axes([0.06, 0.08, 0.92, 0.90])
    pylab.fill_between([0, 6], 0, 1, color='#DD3333', alpha=0.8, label='Error')
    pylab.fill_between(sp.arange(7), 0, tot_herit, color='#22CC44', alpha=0.8, label='Heritable variance')
    pylab.fill_between(sp.arange(7), 0, cis_herit, color='#2255AA', \
                alpha=0.8, label='Heritable variance (cis)')
    pylab.ylabel('Average partition of variance')
    pylab.xlabel('Dist. upstream from gene TSS (kb)')
    pylab.xticks(range(7), [200, 100, 50, 25, 10, 5, 1])
    pylab.legend(loc=1, ncol=3, shadow=True)
    pylab.axis([0, 6, 0, 1])
    pylab.savefig(results_prefix + 'avg_herit_td.png')



    pylab.figure()
    pylab.hist(heritabilities, bins=20, alpha=0.7)
    pylab.xlabel('Pseudo-heritability')
    pylab.xlim((-0.025, 1.025))
    pylab.savefig(results_prefix + '_herits_hist.png')
    pylab.clf()

    ks_list = []
    pm_list = []
    for mm in ['EX', 'LM', 'KW']:
        ks_list.append(pval_infl_dict[mm]['kolmogorov_smirnov'])
        pm_list.append(pval_infl_dict[mm]['median_pvals'])

    png_file_name = results_prefix + '_kolmogorov_smirnov_boxplot.png'
    pylab.figure()
    pylab.boxplot(ks_list)
    pylab.axhline(0, color='k', alpha=0.6, ls='-.')
    pylab.xticks(range(1, 4), ['EX', 'LM', 'KW'])
    pylab.ylabel('Kolmogorov-Smirnov statistic D.')
    pylab.savefig(png_file_name)
    pylab.clf()


    png_file_name = results_prefix + '_median_pvals_boxplot.png'
    pylab.figure()
    pylab.boxplot(pm_list)
    pylab.axhline(0, color='k', alpha=0.6, ls='-.')
    pylab.xticks(range(1, 4), ['EX', 'LM', 'KW'])
    pylab.ylabel('Median p-value bias')
    pylab.savefig(png_file_name)
    pylab.clf()


    x_positions = sp.arange(len(distance_bins), dtype='d64')
    width = 0.25
    png_file_name = results_prefix + '_dist_min_pval_hist.png'
    pylab.axes([0.08, 0.2, 0.91, 0.75])
    for mm, color in zip(['EX', 'LM', 'KW'], ['b', 'c', 'g']):
        l = [dist_min_pval_dict[mm][bin] for bin in distance_bins]
        tot_sum = sum(l)
        l = map(lambda x: x / float(tot_sum), l)
        pylab.bar(x_positions, l, width, color=color, alpha=0.7, label=mm)
        x_positions += width


    pylab.ylabel('Frequency')
    pylab.xticks(x_positions - 3 * width / 2.0, (r'$d \leq 5$', r'$5< d \leq 10$', r'$10< d \leq 25$', \
                        r'$25< d \leq 50$', r'$50< d \leq 100$', r'$d>100$', \
                        'Other chrom.'), rotation='45')
    pylab.xlabel('Distance $d$ (kb) to the smallest p-value from the gene.')
    pylab.xlim((-0.25, len(distance_bins)))
    pylab.legend(loc=2)
    pylab.savefig(png_file_name)
    pylab.clf()


    x_positions = sp.arange(len(radius_bins) + 1, dtype='d64')
    width = 0.25
    png_file_name = results_prefix + 'bonf_sign_bin_hist.png'
    pylab.axes([0.08, 0.22, 0.91, 0.73])
    for mm, color in zip(['EX', 'LM', 'KW'], ['b', 'c', 'g']):
        l = [bonf_sign_bin_dict[mm][bin]['count'] / bonf_sign_bin_dict[mm][bin]['total'] for bin in radius_bins]
        l.append(sign_count[mm] / float(num_genes))
        pylab.bar(x_positions, l, width, color=color, alpha=0.7, label=mm)
        x_positions += width


    pylab.ylabel('Fraction of sign. results')
    pylab.xticks(x_positions - 3 * width / 2.0, ('Within gene', r'$d \leq 1$', r'$d \leq 5$', \
                        r'$d \leq 10$', r'$d \leq 25$', r'$d \leq 50$', \
                        r'$d \leq 100$', 'Whole genome'), rotation='45')
    pylab.xlabel(r'Among SNPs with distance $d$ (kb) from gene.')
    pylab.xlim((-0.25, len(radius_bins) + 1))
    pylab.legend(loc=2)
    pylab.savefig(png_file_name)
    pylab.clf()


    png_file_name = results_prefix + 'cofactor_count_hist.png'
    x_positions = sp.arange(6, dtype='d64')
    width = 0.25
    for criteria, color in zip(['ebics', 'mbonf', 'min_cof_ppa'], ['b', 'c', 'g']):
        bin_counts = list(sp.bincount(cofactor_count_dict[criteria]['num_cofactor_list']))
        while len(bin_counts) < 6:
            bin_counts.append(0)
        pylab.bar(x_positions, bin_counts, width, color=color, alpha=0.7, label=criteria)
        x_positions += width
    pylab.xlabel('Number of cofactor SNPs')
    pylab.ylabel('Number of genes')
    pylab.xticks(x_positions - 3 * width / 2.0, ('0', '1', '2', '3', '4', '5'))
    pylab.legend(loc=1)
    pylab.xlim((-0.2, 6))
    pylab.savefig(png_file_name)
    pylab.clf()


    png_file_name = results_prefix + 'cis_cofactor_count_hist.png'
    x_positions = sp.arange(6, dtype='d64')
    for criteria, color in zip(['ebics', 'mbonf', 'min_cof_ppa'], ['b', 'c', 'g']):
        bin_counts = list(sp.bincount(cofactor_count_dict[criteria]['num_cis_cofactor_list']))
        while len(bin_counts) < 6:
            bin_counts.append(0)
        pylab.bar(x_positions, bin_counts, width, color=color, alpha=0.7, label=criteria)
        x_positions += width
    pylab.xlabel('Number of cis cofactor SNPs')
    pylab.ylabel('Number of genes')
    pylab.xticks(x_positions - 3 * width / 2.0, ('0', '1', '2', '3', '4', '5'))
    pylab.legend(loc=1)
    pylab.xlim((-0.2, 6))
    pylab.savefig(png_file_name)
    pylab.clf()


    png_file_name = results_prefix + 'cofactor_bin_count_hist.png'
    x_positions = sp.arange(9, dtype='d64')
    width = 0.25
    pylab.axes([0.08, 0.2, 0.91, 0.75])
    for criteria, color in zip(['ebics', 'mbonf', 'min_cof_ppa'], ['b', 'c', 'g']):
        cofactor_count_dict[criteria]['bin_counts'] = \
            cofactor_count_dict[criteria]['bin_counts'] / cofactor_count_dict[criteria]['num_found']
        l = list(cofactor_count_dict[criteria]['bin_counts'])
        l.reverse()
        pylab.bar(x_positions, l, width, color=color, alpha=0.7, label=criteria)
        x_positions += width
    pylab.ylabel('Fraction all genes with cofactors.')
    pylab.xlabel(r'Distance $d$ (kb) to cofactor from gene.')
    pylab.xticks(x_positions - 3 * width / 2.0, ('Within gene', r'$1\geq d$', r'$5\geq d$', r'$10\geq d$', \
                        r'$25\geq d$', r'$50\geq d$', r'$100\geq d$', \
                        r'$d>100$', 'Other chrom.'), rotation='45')
    pylab.xlim((-0.2, 9))
    pylab.legend(loc=2)
    pylab.savefig(png_file_name)
    pylab.clf()





def run_parallel_rna_seq_gwas():
    if len(sys.argv) > 5:
        run_id = sys.argv[5]
        call_method_id = int(sys.argv[4])
        temperature = sys.argv[3]
        file_prefix = env['rna_seq_results_dir'] + '%s/cm_%d/' % (temperature, call_method_id)
        if temperature == 'both':
            run_mtmm_gwas(file_prefix, int(sys.argv[1]), int(sys.argv[2]), call_method_id=call_method_id,
                    w_marginal_analysis=False, run_id=run_id)
        else:
            conf_matrix_type = sys.argv[6]
            run_gwas(file_prefix, int(sys.argv[1]), int(sys.argv[2]), temperature,
                call_method_id=call_method_id, conf_matrix_type=conf_matrix_type, run_id=run_id)
    else:
        if len(sys.argv) == 5:
            conf_matrix_type = sys.argv[4]
        call_method_id = int(sys.argv[3])
        temperature = sys.argv[2]
        f = h5py.File(env['rna_seq_data_dir'] + 'rna_seq_gwas_data_cm%d.h5py' % call_method_id)
        g = f[temperature]
        num_genes = len(g['gene_ids'][...])
        f.close()
        print 'Found %d gene expression values' % num_genes
        chunck_size = int(sys.argv[1])
        for i in range(0, num_genes, chunck_size):
            run_parallel(i, i + chunck_size, temperature, call_method_id, conf_matrix_type)



def generate_filtered_csv_files(file_prefix, temperature, call_method_id, run_id):
    pass


def _plot_both_summary_results_(call_method_id, conf_matrix_type, run_id, debug_filter=1.0):

    radius_bins = [0, 1000, 5000, 10000, 25000, 50000, 100000]
    distance_bins = [(0, 5000), (0, 10000), (0, 25000), (0, 50000), (0, 100000), (1, -1), (6, -1)]

    #Initializing the summary dictionary
    summary_dict = {'var_components':{'rho_g':[], 'heritabilities_10':[], 'heritabilities_16':[]}, #, 'sigma_g':[], 'sigma_e':[], 'mle':[]},
            'phen_info':{'correlations':[], 'pseudo_heritabilities_10':[], 'pseudo_heritabilities_16':[]}, #, 'SD':[], 'mean':[]}}
            'dist_min_pval_dict':{}}
    for k in ['g_test', 'gt_test', 'gt_g_test']:
        summary_dict[k] = {'ks':[], 'med_pval_bias':[], 'dist_to_min_pval':[], 'sign_count':0, 'dist_min_pval_dict':{}}
        bonf_sign_dict = {}
        for r in radius_bins:
            bonf_sign_dict[r] = {'total':0, 'count':0}
        summary_dict[k]['bonf_sign_dict'] = bonf_sign_dict
        for bin in distance_bins:
            summary_dict[k]['dist_min_pval_dict'][bin] = 0

    #results = []


    f = h5py.File(env['rna_seq_data_dir'] + 'rna_seq_gwas_data_cm%d.h5py' % call_method_id)
    g = f['both']
    gene_ids = g['gene_ids'][...].tolist()
    num_snps = len(g['positions'])
    bonf_threshold = 1 / (20.0 * num_snps)
    file_prefix = env['rna_seq_results_dir'] + 'both/cm_%d/' % (call_method_id)
    gene_dict = dp.parse_tair_gff_file()#_load_genes_list_('rna_seq_031311_%sC' % temperature)
    num_genes = len(gene_ids)

    for i, gene_tair_id in enumerate(gene_ids):
        if random.random() > debug_filter:
            continue
        chrom = int(gene_tair_id[2])
        chr_prefix = '%schr_%d/rna_seq_%s_both_pid%d_%s' % (file_prefix, chrom, run_id, i, gene_tair_id)
        hist_prefix = '%shistograms/rna_seq_%s_both_pid%d_%s' % (file_prefix, run_id, i, gene_tair_id)
        manhattan_prefix = '%smanhattan_plots/rna_seq_%s_both_pid%d_%s' % (file_prefix, run_id, i, gene_tair_id)
        qq_prefix = '%sqq_plots/rna_seq_%s_both_pid%d_%s' % (file_prefix, run_id, i, gene_tair_id)
        summary_file = chr_prefix + '_summary.pickled'
        #Check if the file is already there.
        print 'Checking %s' % (summary_file)
        if not os.path.isfile(summary_file):
            print "File doesn't exists, hence, skipping this file."
            continue

        print 'File found, loading and summarizing results.'
        sd = cPickle.load(open(summary_file))

#        print 'Loading filtered results'
#        res_dict = {}
#        for k, rk in zip(['g_test', 'gt_test', 'gt_g_test'], ['g', 'gt', 'gt_g']):
#            res_file = '%s_%s_.pvals' % (chr_prefix, rk)
#            r = gr.Result(res_file)
#            res_dict[k] = r

        d = gene_dict[gene_tair_id]
        gene_strand = d['strand']
        try:
            chrom = int(d['chromosome'])
        except Exception:
            raise
        gene = gwaResults.Gene(chromosome=int(d['chromosome']), startPos=d['start_pos'],
                endPos=d['end_pos'], name=gene_tair_id, description=None, dbRef=gene_tair_id,
                tairID=gene_tair_id)
        #print i, gene

        #Aggregate statistics!!!
        for k in ['g_test', 'gt_test', 'gt_g_test']:
            d = sd[k]
            for bin in radius_bins:
                pval = d['analysis_dict']['bin_dict'][bin]['min_pval']
                n_snps = d['analysis_dict']['bin_dict'][bin]['num_snps']
                if n_snps > 0:
                    summary_dict[k]['bonf_sign_dict'][bin]['total'] += 1.0
                    if pval < 1.0 / (20 * n_snps):
                        summary_dict[k]['bonf_sign_dict'][bin]['count'] += 1.0

            dist_min_pval = tuple(d['analysis_dict']['dist_to_min_pval'])
            #if res_dict[k].min_score() < bonf_threshold:
#            summary_dict[k]['sign_count'] += 1
            for bin in distance_bins:
                if dist_min_pval <= bin:
                    summary_dict[k]['dist_min_pval_dict'][bin] += 1
                    break

            summary_dict[k]['dist_to_min_pval'].append(d['analysis_dict']['dist_to_min_pval'])
            summary_dict[k]['ks'].append(d['KS']['D'])
            summary_dict[k]['med_pval_bias'].append(d['med_pval_bias'])

        summary_dict['phen_info']['correlations'].append(sd['phen_info']['correlations'])
        summary_dict['phen_info']['pseudo_heritabilities_10'].append(sd['phen_info']['pseudo_heritabilities'][0])
        summary_dict['phen_info']['pseudo_heritabilities_16'].append(sd['phen_info']['pseudo_heritabilities'][1])
#        summary_dict['phen_info']['SD'].append(sd['phen_info']['SD'])
#        summary_dict['phen_info']['mean'].append(sd['phen_info']['mean'])

#        {'sigma_g':[], 'sigma_e':[], 'rho_g':[], 'heritabilities':[], 'mle':[]}
        summary_dict['var_components']['rho_g'].append(sd['var_components']['rho_g'])
        summary_dict['var_components']['heritabilities_10'].append(sd['var_components']['heritabilities'][0])
        summary_dict['var_components']['heritabilities_16'].append(sd['var_components']['heritabilities'][1])
        summary_dict['var_components']['mle'].append(sd['var_components']['mle'])
        summary_dict['var_components']['sigma_g'].append(sd['var_components']['sigma_g'])
        summary_dict['var_components']['sigma_e'].append(sd['var_components']['sigma_e'])

        #Load GWAS results? cis vs. trans (the old way)

    #Now plot stuff...
    #Specifying the directory for the figures.
    results_prefix = '%ssummary_plots/%s' % (file_prefix, run_id)

    #First heritiabilities

    #Histograms
    pylab.figure()
    p_herits = sp.array([summary_dict['phen_info']['pseudo_heritabilities_10'], summary_dict['phen_info']['pseudo_heritabilities_16']])
    pylab.hist(sp.transpose(p_herits), bins=20, alpha=0.6, color=['g', 'c'], label=['10C', '16C'], range=(0, 1))
    pylab.xlabel('Pseudo-heritability')
    pylab.xlim((-0.025, 1.025))
    pylab.legend()
    pylab.savefig(results_prefix + '_p_herits_hist.png')
    pylab.clf()

    pylab.figure()
    vc_p_herits = sp.array([summary_dict['var_components']['heritabilities_10'], summary_dict['var_components']['heritabilities_16']])
    pylab.hist(sp.transpose(vc_p_herits), bins=20, alpha=0.6, color=['g', 'c'], label=['10C', '16C'], range=(0, 1))
    pylab.xlabel('VC pseudo-heritability')
    pylab.xlim((-0.025, 1.025))
    pylab.legend()
    pylab.savefig(results_prefix + '_vc_p_herits_hist.png')
    pylab.clf()

    #Correlation plots
    pylab.figure()
    xs = summary_dict['phen_info']['pseudo_heritabilities_10']
    ys = summary_dict['var_components']['heritabilities_10']
    print len(xs), len(ys)
    pylab.plot(xs, ys, ls='', marker='.', alpha=0.5, color='k')
    pylab.xlabel('Marginally estimated heritability')
    pylab.ylabel('Jointly estimated heritability')
    pylab.title('10C')
    pylab.savefig(results_prefix + 'herits_10C_corr_plot.png')
    pylab.clf()

    pylab.figure()
    xs = summary_dict['phen_info']['pseudo_heritabilities_16']
    ys = summary_dict['var_components']['heritabilities_16']
    print len(xs), len(ys)
    pylab.plot(xs, ys, ls='', marker='.', alpha=0.5, color='k')
    pylab.xlabel('Marginally estimated heritability')
    pylab.ylabel('Jointly estimated heritability')
    pylab.title('16C')
    pylab.savefig(results_prefix + 'herits_16C_corr_plot.png')
    pylab.clf()

    pylab.figure()
    xs = summary_dict['phen_info']['pseudo_heritabilities_10'] + summary_dict['phen_info']['pseudo_heritabilities_16']
    ys = summary_dict['var_components']['heritabilities_10'] + summary_dict['var_components']['heritabilities_16']
    print len(xs), len(ys)
    pylab.plot(xs, ys, ls='', marker='.', alpha=0.5, color='k')
    pylab.xlabel('Marginally estimated heritability')
    pylab.ylabel('Jointly estimated heritability')
    pylab.title('All')
    pylab.savefig(results_prefix + 'herits_all_corr_plot.png')
    pylab.clf()


    #Population structure confounding quantification
    ks_list = []
    pm_list = []
    for mm in ['g_test', 'gt_test', 'gt_g_test']:
        ks_list.append(summary_dict[mm]['ks'])
        pm_list.append(summary_dict[mm]['med_pval_bias'])

    png_file_name = results_prefix + '_kolmogorov_smirnov_boxplot.png'
    pylab.figure()
    pylab.boxplot(ks_list)
    pylab.axhline(0, color='k', alpha=0.6, ls='-.')
    pylab.xticks(range(1, 4), ['g_test', 'gt_test', 'gt_g_test'])
    pylab.ylabel('Kolmogorov-Smirnov statistic D.')
    pylab.savefig(png_file_name)
    pylab.clf()


    png_file_name = results_prefix + '_median_pvals_boxplot.png'
    pylab.figure()
    pylab.boxplot(pm_list)
    pylab.axhline(0, color='k', alpha=0.6, ls='-.')
    pylab.xticks(range(1, 4), ['g_test', 'gt_test', 'gt_g_test'])
    pylab.ylabel('Median p-value bias')
    pylab.savefig(png_file_name)
    pylab.clf()


#    Where is the strongest signal?
    x_positions = sp.arange(len(distance_bins), dtype='d64')
    width = 0.25
    png_file_name = results_prefix + '_dist_min_pval_hist.png'
    pylab.axes([0.08, 0.2, 0.91, 0.75])
    for k, color in zip(['g_test', 'gt_test', 'gt_g_test'], ['b', 'c', 'g']):
        l = [summary_dict[k]['dist_min_pval_dict'][bin] for bin in distance_bins]
        tot_sum = sum(l)
        l = map(lambda x: x / float(tot_sum), l)
        pylab.bar(x_positions, l, width, color=color, alpha=0.7, label=k)
        x_positions += width


    pylab.ylabel('Frequency')
    pylab.xticks(x_positions - 3 * width / 2.0, (r'$d \leq 5$', r'$5< d \leq 10$', r'$10< d \leq 25$', \
                        r'$25< d \leq 50$', r'$50< d \leq 100$', r'$d>100$', \
                        'Other chrom.'), rotation='45')
    pylab.xlabel('Distance $d$ (kb) to the smallest p-value from the gene.')
    pylab.xlim((-0.25, len(distance_bins)))
    pylab.legend(loc=2)
    pylab.savefig(png_file_name)
    pylab.clf()


    #x_positions = sp.arange(len(radius_bins) + 1, dtype='d64')
    x_positions = sp.arange(len(radius_bins), dtype='d64')
    width = 0.25
    png_file_name = results_prefix + 'bonf_sign_bin_hist.png'
    pylab.axes([0.08, 0.22, 0.91, 0.73])
    for mm, color in zip(['g_test', 'gt_test', 'gt_g_test'], ['b', 'c', 'g']):
        l = [summary_dict[mm]['bonf_sign_dict'][bin]['count'] / summary_dict[mm]['bonf_sign_dict'][bin]['total'] for bin in radius_bins]
#        l.append(summary_dict[mm]['sign_count'] / float(num_genes))
        pylab.bar(x_positions, l, width, color=color, alpha=0.7, label=mm)
        x_positions += width


    pylab.ylabel('Fraction of sign. results')
#    pylab.xticks(x_positions - 3 * width / 2.0, ('Within gene', r'$d \leq 1$', r'$d \leq 5$', \
#                        r'$d \leq 10$', r'$d \leq 25$', r'$d \leq 50$', \
#                        r'$d \leq 100$', 'Whole genome'), rotation='45')
    pylab.xticks(x_positions - 3 * width / 2.0, ('Within gene', r'$d \leq 1$', r'$d \leq 5$', \
                        r'$d \leq 10$', r'$d \leq 25$', r'$d \leq 50$', \
                        r'$d \leq 100$'), rotation='45')
    pylab.xlabel(r'Among SNPs with distance $d$ (kb) from gene.')
    #pylab.xlim((-0.25, len(radius_bins) + 1))
    pylab.xlim((-0.25, len(radius_bins)))
    pylab.legend(loc=2)
    pylab.savefig(png_file_name)
    pylab.clf()


    #Analyze correlations
    pylab.figure()
    genetic_corrs = sp.array(summary_dict['var_components']['rho_g'])
    pylab.hist(genetic_corrs, bins=20, alpha=0.6, color='b')
    pylab.xlabel('Genetic correlations')
    #pylab.xlim((-0.025, 1.025))
    pylab.legend()
    pylab.savefig(results_prefix + '_genetic_correlations.png')
    pylab.clf()

    pylab.figure()
    phen_corrs = sp.array(summary_dict['phen_info']['correlations'])
    pylab.hist(phen_corrs, bins=20, alpha=0.6, color='b')
    pylab.xlabel('Phenotypic correlations')
    #pylab.xlim((-0.025, 1.025))
    pylab.legend()
    pylab.savefig(results_prefix + '_phen_correlations.png')
    pylab.clf()

    pylab.figure()
    pylab.plot(phen_corrs, genetic_corrs, ls='', marker='.', alpha=0.5, color='k')
    pylab.xlabel('Phenotype correlations')
    pylab.ylabel('Genotype correlations')
    pylab.savefig(results_prefix + '_phen_gen_correlations.png')
    pylab.clf()

    pylab.figure()
    pylab.plot(phen_corrs, summary_dict['phen_info']['pseudo_heritabilities_10'], ls='', marker='.', alpha=0.5, color='k')
    pylab.xlabel('Phenotype correlations')
    pylab.ylabel('Heritabilities (10C)')
    pylab.savefig(results_prefix + '_phen_corrs_herits_10C.png')
    pylab.clf()

    pylab.figure()
    pylab.plot(phen_corrs, summary_dict['phen_info']['pseudo_heritabilities_16'], ls='', marker='.', alpha=0.5, color='k')
    pylab.xlabel('Phenotype correlations')
    pylab.ylabel('Heritabilities (16C)')
    pylab.savefig(results_prefix + '_phen_corrs_herits_16C.png')
    pylab.clf()


    pylab.figure()
    pylab.plot(phen_corrs, summary_dict['var_components']['heritabilities_10'], ls='', marker='.', alpha=0.5, color='k')
    pylab.xlabel('Phenotype correlations')
    pylab.ylabel('Jointly est. heritabilities (10C)')
    pylab.savefig(results_prefix + '_phen_corrs_vc_herits_10C.png')
    pylab.clf()

    pylab.figure()
    pylab.plot(phen_corrs, summary_dict['var_components']['heritabilities_16'], ls='', marker='.', alpha=0.5, color='k')
    pylab.xlabel('Phenotype correlations')
    pylab.ylabel('Jointly est. heritabilities (16C)')
    pylab.savefig(results_prefix + '_phen_corrs_vc_herits_16C.png')
    pylab.clf()


def create_lmm_test_data():
    phen_file = '%s_%dC.csv' % (phen_file_prefix, 10)
    phed = pd.parse_phenotype_file(phen_file, with_db_ids=False)  #load phenotype file
    phed.convert_to_averages()
    sd = dp.load_snps_call_method(79)
    indices_to_keep = sd.coordinate_w_phenotype_data(phed, 1, coord_phen=False)  #All phenotypes are ordered the same way, so we pick the first one.
    phed.filter_ecotypes(indices_to_keep)
    phed.filter_near_const_phens(20)
    raw_phen_vals = phed.get_values_list(clone=True)
    trans_types = phed.transform_pids()
    K = kinship.scale_k(sd.get_ibs_kinship_matrix())
    sd.filter_mac_snps(20)
    snps = sd.get_snps()
    h5f = h5py.File('/tmp/lmm_test_10C.h5py')
    h5f.create_dataset('gene_ids', data=phed.get_names(), compression='lzf')
    h5f.create_dataset('snps', data=snps, compression='lzf')
    h5f.create_dataset('ecotypes', data=sd.accessions, compression='lzf')
    h5f.create_dataset('positions', data=sd.get_positions(), compression='lzf')
    h5f.create_dataset('chromosomes', data=sd.get_chr_list(), compression='lzf')
    h5f.create_dataset('kinship', data=K, compression='lzf')
    h5f.create_dataset('phen_values', data=phed.get_values_list(), compression='lzf')
    h5f.create_dataset('raw_phen_values', data=raw_phen_vals, compression='lzf')
    h5f.create_dataset('trans_types', data=trans_types, compression='lzf')
    g = h5f.create_group('emmax_results')
    pids = phed.get_pids()
    for i, pid in enumerate(pids[:100]):
        r = lm.emmax(snps, phed.get_values(pid), K)
        g.create_dataset('r%d' % i, data=r['ps'], compression='lzf')
    h5f.close()


def create_rna_seq_data_file(call_method_id=79):
    h5f = h5py.File('/tmp/rna_seq_gwas_data_cm%d.h5py' % call_method_id)
    h5f.create_group('10C')
    h5f.create_group('16C')
    h5f.create_group('both')
    h5f.create_group('raw_data')


    #10C
    phen_file_10 = '%s/rna_seq_%dC.csv' % (env['phen_dir'], 10)
    phed_10 = pd.parse_phenotype_file(phen_file_10, with_db_ids=False)  #load phenotype file
    phed_10.convert_to_averages()
    sd_10 = dp.load_snps_call_method(call_method_id)
    indices_to_keep_10 = sd_10.coordinate_w_phenotype_data(phed_10, 1, coord_phen=False)  #All phenotypes are ordered the same way, so we pick the first one.
    K_10 = kinship.scale_k(sd_10.get_ibs_kinship_matrix())
    phed_10.filter_ecotypes(indices_to_keep_10)
    phed_10.filter_near_const_phens(20)
    phen_vals_10 = phed_10.get_values_list(clone=True)
    d = phed_10.get_blups(K_10)
#    norm_trans_types_10 = phed_10.transform_pids()
#    norm_trans_phen_vals_10 = phed_10.get_values_list(clone=True)
#    norm_d = phed_10.get_blups(K_10)
    sd_10.filter_mac_snps(20)
    ets_10 = sd_10.accessions
    names_10 = phed_10.get_names()
    r = sd_10.get_mafs()
    macs_10 = r['mafs']
    mafs_10 = r['marfs']


    h5f['10C'].create_dataset('gene_ids', data=names_10, compression='lzf')
    h5f['10C'].create_dataset('phen_values', data=phen_vals_10, compression='lzf')
#    h5f['10C'].create_dataset('most_norm_phen_values', data=norm_trans_phen_vals_10, compression='lzf')
    h5f['10C'].create_dataset('p_herits', data=d['phers'], compression='lzf')
    h5f['10C'].create_dataset('p_herit_pvals', data=d['pvals'], compression='lzf')
    h5f['10C'].create_dataset('residuals', data=d['blup_residuals'], compression='lzf')
#    h5f['10C'].create_dataset('norm_p_herits', data=norm_d['phers'], compression='lzf')
#    h5f['10C'].create_dataset('norm_p_herit_pvals', data=norm_d['pvals'], compression='lzf')
#    h5f['10C'].create_dataset('norm_residuals', data=norm_d['blup_residuals'], compression='lzf')
    h5f['10C'].create_dataset('snps', data=sd_10.get_snps(), compression='lzf')
    h5f['10C'].create_dataset('ecotypes', data=ets_10, compression='lzf')
    h5f['10C'].create_dataset('positions', data=sd_10.get_positions(), compression='lzf')
    h5f['10C'].create_dataset('chromosomes', data=sd_10.get_chr_list(), compression='lzf')
    h5f['10C'].create_dataset('macs', data=macs_10, compression='lzf')
    h5f['10C'].create_dataset('mafs', data=mafs_10, compression='lzf')
    h5f['10C'].create_dataset('kinship', data=K_10, compression='lzf')
#    h5f['10C'].create_dataset('trans_types', data=norm_trans_types_10, compression='lzf')

    #16C
    phen_file_16 = '%s/rna_seq_%dC.csv' % (env['phen_dir'], 16)
    phed_16 = pd.parse_phenotype_file(phen_file_16, with_db_ids=False)  #load phenotype file
    phed_16.convert_to_averages()
    sd_16 = dp.load_snps_call_method(call_method_id)
    indices_to_keep_16 = sd_16.coordinate_w_phenotype_data(phed_16, 1, coord_phen=False)  #All phenotypes are ordered the same way, so we pick the first one.
    K_16 = kinship.scale_k(sd_16.get_ibs_kinship_matrix())
    phed_16.filter_ecotypes(indices_to_keep_16)
    phed_16.filter_near_const_phens(20)
    phed_16.transform_pids(trans_type='ascombe')
    phen_vals_16 = phed_16.get_values_list(clone=True)
    d = phed_16.get_blups(K_16)
#    norm_trans_types_16 = phed_16.transform_pids()
#    norm_trans_phen_vals_16 = phed_16.get_values_list(clone=True)
#    norm_d = phed_16.get_blups(K_16)
    sd_16.filter_mac_snps(20)
    names_16 = phed_16.get_names()
    ets_16 = sd_16.accessions
    r = sd_16.get_mafs()
    macs_16 = r['mafs']
    mafs_16 = r['marfs']

    h5f['16C'].create_dataset('gene_ids', data=names_16, compression='lzf')
#    h5f['16C'].create_dataset('most_norm_phen_values', data=norm_trans_phen_vals_16, compression='lzf')
    h5f['16C'].create_dataset('phen_values', data=phen_vals_16, compression='lzf')
    h5f['16C'].create_dataset('p_herits', data=d['phers'], compression='lzf')
    h5f['16C'].create_dataset('p_herit_pvals', data=d['pvals'], compression='lzf')
    h5f['16C'].create_dataset('residuals', data=d['blup_residuals'], compression='lzf')
#    h5f['16C'].create_dataset('norm_p_herits', data=norm_d['phers'], compression='lzf')
#    h5f['16C'].create_dataset('norm_p_herit_pvals', data=norm_d['pvals'], compression='lzf')
#    h5f['16C'].create_dataset('norm_residuals', data=norm_d['blup_residuals'], compression='lzf')
    h5f['16C'].create_dataset('snps', data=sd_16.get_snps(), compression='lzf')
    h5f['16C'].create_dataset('ecotypes', data=ets_16, compression='lzf')
    h5f['16C'].create_dataset('positions', data=sd_16.get_positions(), compression='lzf')
    h5f['16C'].create_dataset('chromosomes', data=sd_16.get_chr_list(), compression='lzf')
    h5f['16C'].create_dataset('macs', data=macs_16, compression='lzf')
    h5f['16C'].create_dataset('mafs', data=mafs_16, compression='lzf')
    h5f['16C'].create_dataset('kinship', data=K_16, compression='lzf')
#    h5f['16C'].create_dataset('trans_types', data=norm_trans_types_16, compression='lzf')

    #10C and 16C
    common_ets = list(set(sd_10.accessions).union(set(sd_16.accessions)))
    joint_ets = ets_10 + ets_16
    sd_common = dp.load_snps_call_method(call_method_id)
    sd_common.filter_accessions(common_ets)
    pids_to_keep_10 = []
    pids_to_keep_16 = []
    pids_10 = phed_10.get_pids()
    pids_16 = phed_16.get_pids()
    p_herits = []
    for i10, pid_10 in enumerate(pids_10):
        name_10 = names_10[i10]
        i16 = bisect.bisect(names_16, name_10) - 1
        if i16 < len(names_16) and name_10 == names_16[i16]:
            pids_to_keep_10.append(pid_10)
            pids_to_keep_16.append(pids_16[i16])
            p_herits.append([h5f['10C']['p_herits'][i10], h5f['16C']['p_herits'][i16]])
    print "Keeping %d RNAseq expressions" % len(pids_to_keep_10)
    phed_10.filter_phenotypes(pids_to_keep_10)
    phed_16.filter_phenotypes(pids_to_keep_16)
    correlations = []
    for pid_10, pid_16 in it.izip(pids_to_keep_10, pids_to_keep_16):
        correlations.append(phed_10.get_correlation(pid_10, phed_16, pid_16))
    env_vector = sp.matrix([0.0] * len(ets_10) + [1.0] * len(ets_16)).T
    Z = sp.int8(sp.mat(joint_ets).T == sp.mat(common_ets))
    phen_vals_10 = phed_10.get_values_list(clone=True)
    phen_vals_16 = phed_16.get_values_list(clone=True)
    phen_vals = [pvs_10 + pvs_16 for pvs_10, pvs_16 in zip(phen_vals_10, phen_vals_16)]
    std_phen_vals = [((sp.array(pvs_10) - sp.mean(pvs_10)) / sp.std(pvs_10)).tolist() +
            ((sp.array(pvs_16) - sp.mean(pvs_16)) / sp.std(pvs_16)).tolist()
            for pvs_10, pvs_16 in zip(phen_vals_10, phen_vals_16)]
#    norm_trans_phen_vals = [pvs_10 + pvs_16 for pvs_10, pvs_16 in zip(norm_trans_phen_vals_10, norm_trans_phen_vals_16)]
#    std_norm_trans_phen_vals = [((sp.array(pvs_10) - sp.mean(pvs_10)) / sp.std(pvs_10)).tolist() +
#                ((sp.array(pvs_16) - sp.mean(pvs_16)) / sp.std(pvs_16)).tolist()
#                for pvs_10, pvs_16 in zip(norm_trans_phen_vals_10, norm_trans_phen_vals_16)]
    K_common = kinship.scale_k(sd_common.get_ibs_kinship_matrix())
    sd_common.filter_mac_snps(20)
    r = sd_common.get_mafs()
    macs = r['mafs']
    mafs = r['marfs']


    h5f['both'].create_dataset('gene_ids', data=phed_10.get_names(), compression='lzf')
#    h5f['both'].create_dataset('most_norm_phen_values', data=norm_trans_phen_vals, compression='lzf')
    h5f['both'].create_dataset('phen_values', data=phen_vals, compression='lzf')
#    h5f['both'].create_dataset('std_most_norm_phen_values', data=std_norm_trans_phen_vals, compression='lzf')
    h5f['both'].create_dataset('std_phen_values', data=std_phen_vals, compression='lzf')
    h5f['both'].create_dataset('correlations', data=correlations, compression='lzf')
    h5f['both'].create_dataset('p_herits', data=p_herits, compression='lzf')
    h5f['both'].create_dataset('snps', data=sd_common.get_snps(), compression='lzf')
    h5f['both'].create_dataset('ecotypes', data=common_ets, compression='lzf')
    h5f['both'].create_dataset('env_vector', data=env_vector, compression='lzf')
    h5f['both'].create_dataset('Z', data=Z, compression='lzf')
    h5f['both'].create_dataset('positions', data=sd_common.get_positions(), compression='lzf')
    h5f['both'].create_dataset('chromosomes', data=sd_common.get_chr_list(), compression='lzf')
    h5f['both'].create_dataset('macs', data=macs, compression='lzf')
    h5f['both'].create_dataset('mafs', data=mafs, compression='lzf')
    h5f['both'].create_dataset('kinship', data=K_common, compression='lzf')
    h5f.close()


#def generate_test_data(call_method_id=75):
#    """
#    Generate test data for Arthur and Oliver.
#    """
#    import random
#    import h5py
#    #10C
#    phen_file_10 = '%s/rna_seq_%dC.csv' % (env['phen_dir'], 10)
#    phed_10 = pd.parse_phenotype_file(phen_file_10, with_db_ids=False)  #load phenotype file
#    phed_10.filter_near_const_phens(20)
#    pids_10 = random.sample(phed_10.get_pids(), 100)
#    pids_10.sort()
#    phed_10.filter_phenotypes(pids_10)
#    names_10 = phed_10.get_names()
#    phen_file_16 = '%s/rna_seq_%dC.csv' % (env['phen_dir'], 16)
#    phed_16 = pd.parse_phenotype_file(phen_file_16, with_db_ids=False)  #load phenotype file
#    phed_16.filter_near_const_phens(20)
#    #Find the corresponding pids
#    pids_16 = phed_16.get_pids()
#    pids_to_keep_16 = []
#    pids_to_keep_10 = []
#    names_16 = phed_16.get_names()
#    for n, pid in zip(names_10, pids_10):
#        i_16 = bisect.bisect(names_16, n) - 1
#        if names_16[i_16] == n:
#            pids_to_keep_16.append(pids_16[i_16])
#            pids_to_keep_10.append(pid)
#    pids_to_keep_16.sort()
#    phed_10.filter_phenotypes(pids_to_keep_10)
#    phed_16.filter_phenotypes(pids_to_keep_16)
#    pids_10 = pids_to_keep_10
#    pids_16 = pids_to_keep_16
#    common_ets = list(set(phed_10.get_ecotypes(pids_10[0])).intersection(set(phed_16.get_ecotypes(pids_16[0]))))
#    phed_10.filter_ecotypes_2(common_ets)
#    phed_16.filter_ecotypes_2(common_ets)
#
#    sd = dp.load_snps_call_method(call_method_id)
#    indices_to_keep = sd.coordinate_w_phenotype_data(phed_10, pids_10[0], coord_phen=False)  #All phenotypes are ordered the same way, so we pick the first one.
#    K = lm.scale_k(sd.get_ibs_kinship_matrix())
#    sd.filter_mac_snps(20)
#    r = sd.get_mafs()
#    macs = r['mafs']
#    mafs = r['marfs']
#    ets = sd.accessions
#
#    phed_10.filter_ecotypes(indices_to_keep)
#    phed_16.filter_ecotypes(indices_to_keep)
#    phen_vals_10 = phed_10.get_values_list(clone=True)
#    phen_vals_16 = phed_16.get_values_list(clone=True)
#    phen_vals = [pvs_10 + pvs_16 for pvs_10, pvs_16 in zip(phen_vals_10, phen_vals_16)]
#
#    #Calculate BLUPs and heritabilites
#    d10 = phed_10.get_blups(K)
#    d16 = phed_16.get_blups(K)
#    p_herits = zip(d10['phers'], d16['phers'])
#    #Calculate correlations
#    correlations = [phed_10.get_correlation(pid_10, phed_16, pid_16) for pid_10, pid_16 in zip(pids_10, pids_16)]
#    #Calculate ENV vector
#    ets = sd.accessions
#    joint_ets = ets + ets
#    env_vector = sp.matrix([0.0] * len(ets) + [1.0] * len(ets)).T
#    #incidence matrix
#    Z = sp.int8(sp.mat(joint_ets).T == sp.mat(common_ets))
#
#    #Save Oliver's data in HDF5 file.
#    h5f = h5py.File('/tmp/mtmm_test.h5py')
#    h5f.create_dataset('phen_values', data=phen_vals, compression='lzf')
#    h5f.create_dataset('gene_ids', data=phed_10.get_names(), compression='lzf')
#    h5f.create_dataset('correlations', data=correlations, compression='lzf')
#    h5f.create_dataset('p_herits', data=p_herits, compression='lzf')
#    h5f.create_dataset('snps', data=sd.get_snps(), compression='lzf')
#    h5f.create_dataset('ecotypes', data=common_ets, compression='lzf')
#    h5f.create_dataset('env_vector', data=env_vector, compression='lzf')
#    h5f.create_dataset('Z', data=Z, compression='lzf')
#    h5f.create_dataset('positions', data=sd.get_positions(), compression='lzf')
#    h5f.create_dataset('chromosomes', data=sd.get_chr_list(), compression='lzf')
#    h5f.create_dataset('macs', data=macs, compression='lzf')
#    h5f.create_dataset('mafs', data=mafs, compression='lzf')
#    h5f.create_dataset('kinship', data=K, compression='lzf')
#    h5f.close()
#
#    #Save Arthur's data in csv files.
#    phed_10.write_to_file('/tmp/mtmm_test_10C.csv')
#    phed_16.write_to_file('/tmp/mtmm_test_16C.csv')
#    lm.save_kinship_in_text_format('/tmp/mtmm_test_kinship.csv', K, sd.accessions)
#    sd.writeToFile('/tmp/mtmm_test_genotypes.csv')



#def update_blups(call_method_id=79):
#    """
#    """
#    h5f = h5py.File('/tmp/rna_seq_gwas_data_cm%d.h5py' % call_method_id)
#    g_10 = h5f['10C']
##    resids = []
##    for i, p_herit in enumerate(g_10['p_herits']):
##        pvs = g_10['phen_values'][i]
##        K = sp.matrix(g_10['kinship'])
##        pvs -= sp.mean(pvs)
##        Y = sp.matrix(pvs).T
##        delta = (1 - p_herit) / p_herit
##        M = (K + delta * sp.eye(K.shape[0]))
##        u_blup = K * (M.I * Y)
##        resids.append(Y - u_blup)
##    h5f['10C'].create_dataset('residuals', data=resids, compression='lzf')
##
##
##
#    g_16 = h5f['16C']
##    resids = []
##    for i, p_herit in enumerate(g_16['p_herits']):
##        pvs = g_16['phen_values'][i]
##        K = sp.matrix(g_16['kinship'])
##        pvs -= sp.mean(pvs)
##        Y = sp.matrix(pvs).T
##        delta = (1 - p_herit) / p_herit
##        M = (K + delta * sp.eye(K.shape[0]))
##        u_blup = K * (M.I * Y)
##        resids.append(Y - u_blup)
##    h5f['16C'].create_dataset('residuals', data=resids, compression='lzf')
#    g_both = h5f['both']
#    gene_ids_10 = g_10['gene_ids'][...]
#    resids_10 = g_10['residuals'][...]
#    gene_ids_16 = g_16['gene_ids'][...]
#    resids_16 = g_16['residuals'][...]
#    gene_ids = g_both['gene_ids'][...]
#    residuals = []
#    for gi in gene_ids:
#        i10 = bisect.bisect(gene_ids_10, gi) - 1
#        i16 = bisect.bisect(gene_ids_16, gi) - 1
#        r10 = resids_10[i10].flatten().tolist()
#        r16 = resids_16[i16].flatten().tolist()
#        residuals.append(r10 + r16)
#    h5f['both'].create_dataset('residuals', data=residuals, compression='lzf')


if __name__ == '__main__':
#    summarize_st_results('10C', 79, conf_matrix_type='none')
    summarize_st_results('10C', 79, conf_matrix_type='k_expr_standard')
#    print sys.argv
#    if  len(sys.argv) > 4:
#        run_parallel_rna_seq_gwas()
#    else: #plotting stuff
#        temperature = sys.argv[1]
#        call_method_id = int(sys.argv[2])
#        conf_matrix_type = 'none'
#        if len(sys.argv) == 4:
#            conf_matrix_type = sys.argv[3]
#        plot_summary_results(temperature=temperature, call_method_id=call_method_id,
#                    conf_matrix_type=conf_matrix_type)
    print  'Done'
