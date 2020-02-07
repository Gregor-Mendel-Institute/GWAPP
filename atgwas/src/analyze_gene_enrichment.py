"""
Identifying gene enrichment
"""
import gwaResults as gr
import pylab as plt
#import matplotlib
import scipy as sp
import random
from scipy import stats as st
import env
import sys


def plot_gene_distances(cgl_id=145):
    cand_genes = gr.getCandidateGeneList(cgl_id, passwd='bjazz32')
    cand_gene_distances = []
    for i in range(len(cand_genes)):
        g1 = cand_genes[i]
        for j in range(i + 1, len(cand_genes)):
            g2 = cand_genes[j]
            if g1.chromosome == g2.chromosome:
                dist = max(max(g1.startPos - g2.endPos, 0), max(g2.startPos - g1.endPos, 0))
                if dist < 500000:
                    cand_gene_distances.append(dist)
    print cand_gene_distances
    plt.figure(figsize=(8, 5))
    plt.hist(cand_gene_distances, bins=100)
    plt.savefig("/Users/bjarnivilhjalmsson/tmp/Gene_distances_cgl" + str(cgl_id) + ".pdf")
    plt.clf()
    all_genes = gr.get_gene_list(passwd='bjazz32')
    all_gene_distances = []
    for i in range(len(all_genes)):
        g1 = all_genes[i]
        for j in range(i + 1, len(all_genes)):
            g2 = all_genes[j]
            if g1.chromosome == g2.chromosome:
                dist = max(max(g1.startPos - g2.endPos, 0), max(g2.startPos - g1.endPos, 0))
                if dist < 500000:
                    all_gene_distances.append(dist)
        print i
    plt.hist(all_gene_distances, bins=100)
    plt.savefig("/Users/bjarnivilhjalmsson/tmp/Gene_distances.pdf")
    (D_stat, pval) = stats.ks_2samp(all_gene_distances, cand_gene_distances)
    print pval




def calc_enrichment(all_genes, cg_indices, regions):
    """
    Calculate the enrichment (efficiently iterating over all genes.)
    """
    num_genes = len(all_genes)
    num_cand_genes = len(cg_indices)
    g_iter = enumerate(all_genes)
    g_i, g = g_iter.next()
    g_end_chr_pos = sp.array([g.chromosome, g.endPos])
    g_start_chr_pos = sp.array([g.chromosome, g.startPos])
    r_iter = enumerate(regions)
    r_i, r = r_iter.next()
    r_start_chr_pos = r[0]
    r_end_chr_pos = r[1]
    num_close_genes = 0
    num_close_cand_genes = 0
    cg_iter = enumerate(cg_indices)
    cg_ii, cg_i = cg_iter.next()
    obs_cg_indices = []
    while g_i < num_genes - 1 and r_i < len(regions) - 1:
        #count until overlap
        #if g_i % 100 == 0: print g_i
        while g_i < num_genes - 1 and tuple(g_end_chr_pos) < tuple(r_start_chr_pos):
            g_i, g = g_iter.next()
            g_end_chr_pos = sp.array([g.chromosome, g.endPos])
            g_start_chr_pos = sp.array([g.chromosome, g.startPos])

        while r_i < len(regions) - 1 and  tuple(r_end_chr_pos) < tuple(g_start_chr_pos):
            r_i, r = r_iter.next()
            r_start_chr_pos = r[0]
            r_end_chr_pos = r[1]

        while g_i < num_genes - 1 and tuple(g_end_chr_pos) >= tuple(r_start_chr_pos) \
                and tuple(r_end_chr_pos) >= tuple(g_start_chr_pos):
            #there is an overlap
            num_close_genes += 1
            while cg_ii < num_cand_genes - 1 and cg_i < g_i:
                cg_ii, cg_i = cg_iter.next()
            if g_i == cg_i:
                num_close_cand_genes += 1
                obs_cg_indices.append(g_i)
            g_i, g = g_iter.next()
            g_end_chr_pos = sp.array([g.chromosome, g.endPos])
            g_start_chr_pos = sp.array([g.chromosome, g.startPos])
        #print g_i, r_i
        #print num_close_genes, num_close_cand_genes
        #pdb.set_trace()

    #r1 = (num_close_cand_genes / float(num_close_genes))
    #r2 = (num_cand_genes / float(num_genes))
    return (num_close_cand_genes, num_close_genes, obs_cg_indices)


def get_chi_square_pval(num_close_cand_genes, num_close_genes, num_cand_genes, num_genes):
    num_exp_ccg = (num_close_genes) * num_cand_genes / float(num_genes)
    num_exp_cg = (num_close_genes) * (1 - num_cand_genes / float(num_genes))
    num_exp_dcg = (float(num_genes) - num_close_genes) * (num_cand_genes / float(num_genes))
    num_exp_dg = (float(num_genes) - num_close_genes) * (1 - num_cand_genes / float(num_genes))
    num_obs_ccg = num_close_cand_genes
    num_obs_cg = num_close_genes - num_close_cand_genes
    num_obs_dcg = num_cand_genes - num_close_cand_genes
    num_obs_dg = num_genes - num_close_genes
    f_obs = sp.array([num_obs_ccg, num_obs_cg, num_obs_dcg, num_obs_dg])
    f_exp = sp.array([num_exp_ccg, num_exp_cg, num_exp_dcg, num_exp_dg])
    chi_sq_stat, chi_sq_pval = st.chisquare(f_obs, f_exp, 2)
    if num_obs_ccg >= num_exp_ccg:
        chi_sq_pval = chi_sq_pval / 2
    else:
        chi_sq_pval = 1 - chi_sq_pval / 2
    print 'Cand. gene enrichment p-value from Chi-square test:', chi_sq_pval
    return chi_sq_pval, chi_sq_stat


def get_gene_perm_pval(obs_stat, regions, all_genes, cand_gene_indices, num_perm=500, early_stop_threshold=25):
    print "Doing gene permutations"
    num_genes = len(all_genes)
    num_cand_genes = len(cand_gene_indices)
    r2 = num_cand_genes / float(num_genes)
    perm_stats = []
    if obs_stat != float('-inf'):
        sign_count = 0
        for perm_i in range(num_perm):
            perm_cand_gene_indices = random.sample(range(num_genes), num_cand_genes)
            perm_cand_gene_indices.sort()
            perm_enrichments = calc_enrichment(all_genes, perm_cand_gene_indices, regions)
            r1 = perm_enrichments[0] / float(perm_enrichments[1])
            if r1 == 0.0:
                perm_stat = float('-inf')
            else:
                perm_stat = sp.log(r1 / r2)
            perm_stats.append(perm_stat)
            if perm_stat >= obs_stat:
                sign_count += 1
                if sign_count == early_stop_threshold:
                    break
            sys.stdout.write('.')
            sys.stdout.flush()
        sys.stdout.write('\n')
        sys.stdout.flush()
        perm_stats.sort()
        p_val = sign_count / float(perm_i + 1)
    else:
        p_val = 1.0
    print 'Permutation p-value estimate: % f' % p_val
    return p_val, perm_stats



def get_snps_perm_pval(obs_stat, regions, all_genes, cand_gene_indices, chrom_ends, num_perm=500,
        early_stop_threshold=25):
    print "Doing SNPs permutations"
    num_genes = len(all_genes)
    num_cand_genes = len(cand_gene_indices)
    r2 = num_cand_genes / float(num_genes)
    def _get_perm_regions_(region_dict):
        new_regions = []
        for chrom in region_dict:
            chrom_end = chrom_ends[chrom - 1]
            regions = region_dict[chrom]
            while True:
                shift = int(random.random() * chrom_end)
                for region in regions:
                    if (region[0][1] + shift) < chrom_end \
                            and (region[1][1] + shift) > chrom_end:
                        break #for loop
                else:
                    break #while loop
            start_pos_list = []
            for i, region in enumerate(regions):
                start_pos = (region[0][1] + shift) % chrom_end
                region[0][1] = start_pos
                region[1][1] = (region[1][1] + shift) % chrom_end
                start_pos_list.append((start_pos, i))
            start_pos_list.sort()
            for sp, i in start_pos_list:
                new_regions.append(regions[i])
        return new_regions


    region_dict = {1:[], 2:[], 3:[], 4:[], 5:[]}
    for region in regions:
        region_dict[region[0][0]].append(region)

    perm_stats = []
    if obs_stat != float('-inf'):
        sign_count = 0
        for perm_i in range(num_perm):
            perm_regions = _get_perm_regions_(region_dict)
            perm_enrichments = calc_enrichment(all_genes, cand_gene_indices, perm_regions)
            r1 = perm_enrichments[0] / float(perm_enrichments[1])
            if r1 == 0.0:
                perm_stat = float('-inf')
            else:
                perm_stat = sp.log(r1 / r2)
            perm_stats.append(perm_stat)
            if perm_stat >= obs_stat:
                sign_count += 1
                if sign_count == early_stop_threshold:
                    break
            sys.stdout.write('.')
            sys.stdout.flush()
        sys.stdout.write('\n')
        sys.stdout.flush()
        perm_stats.sort()
        p_val = sign_count / float(perm_i + 1)
    else:
        p_val = 1.0
    print 'Permutation p-value estimate: % f' % p_val
    return p_val, perm_stats


def _test_cand_gene_enrichment_():
    #load result
    res_files = []
    res_files.append((env.env['tmp_dir'] + 'ss_pid1_seed_size_2n_emmax_none.pvals', 'seed_size_2n EMMAX'))
#    res_files.append((env.env['tmp_dir'] + 'ss_pid2_seed_size_3n_2x4_emmax_none.pvals', 'seed_size_3n_2x4 EMMAX'))
#    res_files.append((env.env['tmp_dir'] + 'ss_pid3_seed_size_3n_4x2_emmax_none.pvals', 'seed_size_3n_4x2 EMMAX'))
#    res_files.append((env.env['tmp_dir'] + 'ss_pid4_seed_size_spss_emmax_none.pvals', 'seed_size_spss EMMAX'))
#    res_files.append((env.env['tmp_dir'] + 'ss_pid1_seed_size_2n_kw_none.pvals', 'seed_size_2n KW'))
#    res_files.append((env.env['tmp_dir'] + 'ss_pid2_seed_size_3n_2x4_kw_none.pvals', 'seed_size_3n_2x4 KW'))
#    res_files.append((env.env['tmp_dir'] + 'ss_pid3_seed_size_3n_4x2_kw_none.pvals', 'seed_size_3n_4x2 KW'))
#    res_files.append((env.env['tmp_dir'] + 'ss_pid4_seed_size_spss_kw_none.pvals', 'seed_size_spss KW'))
#    res_files.append((env.env['results_dir'] + 'donald_duck_pid0_25-DHBAG_emmax_log_trans.pvals', '25-DHBAG'))
#    res_files.append((env.env['results_dir'] + 'dilkes_pid1_23-DHBAG_emmax_log_trans.pvals', '23-DHBAG'))
#    res_files.append((env.env['results_dir'] + 'dilkes_pid2_25-DHBAP_emmax_log_trans.pvals', '25-DHBAP'))
#    res_files.append((env.env['results_dir'] + 'dilkes_pid3_23-DHBAP_emmax_log_trans.pvals', '23-DHBAP'))
#    res_files.append((env.env['results_dir'] + 'dilkes_pid4_DHBAG_emmax_log_trans.pvals', 'DHBAG'))
#    res_files.append((env.env['results_dir'] + 'dilkes_pid5_DHBAP_emmax_none.pvals', 'DHBAP'))
#    res_files.append((env.env['results_dir'] + 'dilkes_pid6_total_emmax_log_trans.pvals', 'total'))
#    res_files.append((env.env['results_dir'] + 'dilkes_pid7_25-DHBAg_emmax_sqrt_trans.pvals', '25-DHBA_g'))
#    res_files.append((env.env['results_dir'] + 'dilkes_pid8_23-DHBAg_emmax_log_trans.pvals', '23-DHBA_g'))
#    res_files.append((env.env['results_dir'] + 'dilkes_pid9_Gpct_emmax_none.pvals', 'Gpct'))
#    res_files.append((env.env['results_dir'] + 'dilkes_pid10_25pct_emmax_none.pvals', '25pct'))
#    res_files.append((env.env['results_dir'] + 'dilkes_pid12_SM_emmax_none.pvals', 'SM'))
#    res_files.append((env.env['results_dir'] + 'dilkes_pid14_unkown1_emmax_log_trans.pvals', 'unknown1'))
#    res_files.append((env.env['results_dir'] + 'donald_duck_pid1_telomere_length_emmax_none.pvals', 'telomere_emmax'))
#    res_files.append((env.env['results_dir'] + 'donald_duck_pid1_telomere_length_emmax_sqrt_trans.pvals', 'telomere_emmax_sqrt_trans'))
#    res_files.append((env.env['results_dir'] + 'donald_duck_pid1_telomere_length_kw_none.pvals', 'telomere_kw'))
    cgl_file = env.env['phen_dir'] + 'seeds_cgl.csv'
    #cgl_file = env.env['phen_dir'] + 'teleomere_cgl.csv'
    #cgl_file = env.env['phen_dir'] + 'Dilkes_data_candidate_gene_list.csv'
    for res_file, name in reversed(res_files):

        r = gr.Result(res_file)
        r.neg_log_trans()
        r.filter_attr('mafs', 10)
        res = r.candidate_gene_enrichments(cgl_file=cgl_file, methods=['snps_perm', 'gene_perm'],
                    pval_thresholds=[ 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0005, 0.0002, 0.0001, 0.00005],
                    num_perm=500, obs_genes_file=env.env['tmp_dir'] + 'enrichments_' + name + '_obs_genes.csv',
                    file_prefix=env.env['tmp_dir'] + 'enrichments_' + name, early_stop_threshold=25)



if __name__ == '__main__':
    pass
    #cand_genes_enrichment(5,43)
    #plot_gene_distances()
    #plot_gene_pvalue_count()
    #plot_enrichment_stat([1,2,3,4,5,6,7],145,gene_window_size=20000)
    #plot_enrichment_histogram([1,2,3,4,5,6,7],145)
    #plot_enrichment_ratios([1,2,3,4,5,6,7],145,use_fixed_pvalue_limits=False,gene_window_size=10000,max_num_pvals=160000,num_steps=160)
    _test_cand_gene_enrichment_()

