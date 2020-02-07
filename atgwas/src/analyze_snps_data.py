"""
Contains basic functions to analyze SNP data, e.g. LD levels, polymorphism rates, Fst, recombination, etc.
"""
import sys
from env import *
import dataParsers as dp
import scipy as sp
import scipy.stats as st
from scipy import linalg
import linear_models as lm
import pylab
import os
import cPickle
import math
import random
import pdb
import time
import h5py
from itertools import *
from bisect import bisect
min_float = 5e-324


cm_num_snps_dict = {75:214051, 76:4988387, 78:1031827} #After MAC filtering
chromosome_ends = {1:30425061, 2:19694800, 3:23456476, 4:18578714, 5:26974904}
def get_chromsome_grid(tick_size=50000):
    chrom_grid_dict = {}
    chromosomes = [1, 2, 3, 4, 5]
    for chrom1 in chromosomes:
        clen1 = chromosome_ends[chrom1]
        for chrom2 in chromosomes[:chrom1]:
            clen2 = chromosome_ends[chrom2]
            n1 = 1 + clen1 / tick_size
            n2 = 1 + clen2 / tick_size
            chrom_grid_dict[(chrom1, chrom2)] = sp.empty((n1, n2))

#def update_grid(chrom_grid,):


cm_data_format_dict = {75:'binary', 76:'binary', 78:'diploid_int'}

def test_correlation(sample_num=4000, mac_filter=15, debug_filter=0.05):
    dtype = 'single' #To increase matrix multiplication speed... using 32 bits.
    sd = dp.parse_numerical_snp_data(env['data_dir'] + '250K_t72.csv.binary',
                    filter=debug_filter)
    sd.filter_mac_snps(mac_filter)
    snps = sd.getSnps()[:]
    kinship_matrix_file = env['data_dir'] + 'snp_corr_kinship_cm72.pickled'
    if not os.path.isfile(kinship_matrix_file):
        K = sd.get_snp_cov_matrix()
        lm.save_kinship_to_file(kinship_matrix_file, K, sd.accessions)
    else:
        K = lm.load_kinship_from_file(kinship_matrix_file, dtype='single')
    h_inverse_matrix_file = env['data_dir'] + 'snp_corr_kinship_h_inv_cm72.pickled'
    if not os.path.isfile(h_inverse_matrix_file):
        H_sqrt = lm.cholesky(K)
        H_sqrt_inv = (H_sqrt).I
        with file(h_inverse_matrix_file, 'wb') as f:
            cPickle.dump(H_sqrt_inv.tolist(), f)
    else:
        with file(h_inverse_matrix_file) as f:
            H_sqrt_inv = sp.mat(cPickle.load(f))

    t_snps_file = env['data_dir'] + 'snps_trans_kinship_cm72.pickled'
    if debug_filter < 1 or not os.path.isfile(t_snps_file):
        if os.path.isfile(kinship_matrix_file):
            n_snps_mat = sd.get_normalized_snps().T
        t_snps_mat = n_snps_mat * H_sqrt_inv
        t_snps = t_snps_mat.tolist()
        with open(t_snps_file, 'wb') as f:
            cPickle.dump(t_snps, f)
    else:
        with open(t_snps_file) as f:
            t_snps = cPickle.load(f)

    cp_list = sd.getChrPosList()
    l = range(len(cp_list))
    y_sample = random.sample(l, sample_num)
    for i in reversed(sorted(y_sample)):
        del l[i]
    x_sample = random.sample(l, sample_num)

    res_d = {}
    for i, xi in enumerate(x_sample):
        (x_c, x_p) = cp_list[xi]
        x_snp = snps[xi]
        x_snp_t = t_snps[xi]
        print 'SNP %d: chromosome=%d, position=%d' % (i, x_c, x_p)
        for j, yi in enumerate(y_sample):
            #pdb.set_trace()
            (y_c, y_p) = cp_list[yi]
            y_snp = snps[yi]
            y_snp_t = t_snps[yi]
            (r, pearson_pval) = st.pearsonr(x_snp, y_snp)
            r2 = r * r
            if r2 > 0.05:
                (r_t, pearson_pval_t) = st.pearsonr(x_snp_t, y_snp_t)
                r2_t = r_t * r_t
                if pearson_pval == 0:
                    pearson_pval = min_float
                if pearson_pval_t == 0:
                    pearson_pval_t = min_float
                res_d[(yi, xi)] = [r2_t, r2, pearson_pval, pearson_pval_t]

    xs = []
    for yi, xi in res_d:
        xs.append(res_d[(yi, xi)][0])
    pylab.hist(xs, bins=20, log=True)
    pylab.xlabel('Transformed $r^2$')
    pylab.savefig(env['results_dir'] + 'r2_t_hist.png')


    pylab.clf()

#    #xs = -sp.log10(xs)
#    #ys = -sp.log10(ys)
#    pylab.plot(xs, ys, '.', markersize=2)
#    pylab.savefig(env['results_dir'] + 'norm_emma_log_pvals.png')
#    pylab.clf()

    xs = []
    for yi, xi in res_d:
        xs.append(res_d[(yi, xi)][1])
    pylab.hist(xs, bins=20, log=True)
    pylab.xlabel('Standard $r^2$')
    pylab.savefig(env['results_dir'] + 'r2_hist.png')
    pylab.clf()


    xs = []
    ys = []
    for yi, xi in res_d:
        ys.append(res_d[(yi, xi)][0])
        xs.append(res_d[(yi, xi)][1])
    x_max = max(xs)
    y_max = max(ys)
    lim = min(x_max, y_max)
    pylab.plot(xs, ys, '.', markersize=2)
    pylab.plot([0, lim], [0, lim], color='g', alpha=0.5)
    pylab.xlabel('Standard $r^2$')
    pylab.ylabel('Transformed $r^2$')
    pylab.savefig(env['results_dir'] + 'r2_r2_t.png')
    pylab.clf()

    pylab.hexbin(xs, ys)
    pylab.xlabel('Standard $r^2$')
    pylab.ylabel('Transformed $r^2$')
    cb = pylab.colorbar()
    cb.set_label('Count')
    pylab.savefig(env['results_dir'] + 'r2_r2_t_3Dhist.png')
    pylab.clf()

    xs = []
    ys = []
    for yi, xi in res_d:
        ys.append(res_d[(yi, xi)][3])
        xs.append(res_d[(yi, xi)][2])
    xs = -sp.log10(xs)
    ys = -sp.log10(ys)
    x_max = max(xs)
    y_max = max(ys)
    lim = min(x_max, y_max)
    pylab.plot(xs, ys, '.', markersize=2)
    pylab.plot([0, lim], [0, lim], color='g', alpha=0.5)
    pylab.xlabel('Standard $r^2$ -log(p-value)')
    pylab.ylabel('Transformed $r^2$ -log(p-value)')
    pylab.savefig(env['results_dir'] + 'pval_pval_t.png')
    pylab.clf()

    pylab.hexbin(xs, ys)
    pylab.xlabel('Standard $r^2$ -log(p-value)')
    pylab.ylabel('Transformed $r^2$ -log(p-value)')
    cb = pylab.colorbar()
    cb.set_label('Count')
    pylab.savefig(env['results_dir'] + 'pval_pval_t_3Dhist.png')


#    ys = -sp.log10(ys)
#    pylab.plot(xs, ys, '.', markersize=2)
#    pylab.savefig(env['results_dir'] + 'norm_r2_emma_pval.png')



def calc_r2_levels(file_prefix, x_start_i, x_stop_i, call_method_id=78, data_format='diploid_int',
        mac_filter=15, save_threshold=0.2, save_threshold2=0.3, debug_filter=1):
    """
    Returns statistics on LD levels, and plot them.
    """

    dtype = 'single' #To increase matrix multiplication speed... using 32 bits.
    sd = dp.load_snps_call_method(call_method_id=call_method_id, data_format=data_format, debug_filter=debug_filter, min_mac=mac_filter)
    #sd.filter_mac_snps(mac_filter)
    h_inverse_matrix_file = env['data_dir'] + 'snp_cov_mat_h_inv_cm%d.pickled' % (call_method_id)
    if not os.path.isfile(h_inverse_matrix_file):
        K = sd.get_snp_cov_matrix()
        H_sqrt = lm.cholesky(K)
        H_sqrt_inv = (H_sqrt).I
        with file(h_inverse_matrix_file, 'wb') as f:
            cPickle.dump(H_sqrt_inv, f, protocol=2)
    else:
        with file(h_inverse_matrix_file) as f:
            H_sqrt_inv = cPickle.load(f)

    cps_list = sd.getChrPosSNPList()
    x_cps = cps_list[x_start_i:x_stop_i]
    y_cps = cps_list
    result_dict = {}
    n = len(sd.accessions)
    print 'Starting calculation'
    sys.stdout.flush()
    hdf5_file_name = file_prefix + '_x_' + str(x_start_i) + '_' + str(x_stop_i) + ".hdf5"
    h5_file = h5py.File(hdf5_file_name, 'w')
    for i, (x_c, x_p, x_snp) in enumerate(x_cps):
        print '%d: chromosome=%d, position=%d' % (i, x_c, x_p)
        #Normalize SNP..
        xs = sp.array(x_snp)
        t_x_snp = sp.dot(((xs - sp.mean(xs)) / sp.std(xs)), H_sqrt_inv).T
        s1 = time.time()
        y_cs = []
        y_ps = []
        r2s = []
        t_r2s = []
        ps = []
        t_ps = []
        n_saved = 0
        for (y_c, y_p, y_snp) in reversed(y_cps):
            if (x_c, x_p) < (y_c, y_p):
                ys = sp.array(y_snp)
                mac = ys.sum()
                (r, pearson_pval) = st.pearsonr(xs, ys)
                r2 = r * r
                if x_c == y_c and y_p - x_p <= 50000 and r2 > save_threshold2 :
                    t_y_snp = sp.dot(((ys - sp.mean(ys)) / sp.std(ys)), H_sqrt_inv).T
                    (t_r, t_pearson_pval) = st.pearsonr(t_x_snp, t_y_snp) #Done twice, but this is fast..
                    t_r, t_pearson_pval = float(t_r), float(t_pearson_pval)
                    t_r2 = t_r * t_r
                    y_cs.append(y_c)
                    y_ps.append(y_p)
                    r2s.append(r2)
                    t_r2s.append(t_r2)
                    ps.append(pearson_pval)
                    t_ps.append(t_pearson_pval)
                    n_saved += 1

                elif ((x_c == y_c and y_p - x_p > 50000) or x_c != y_c) and r2 > save_threshold:
                    t_y_snp = sp.dot(((ys - sp.mean(ys)) / sp.std(ys)), H_sqrt_inv).T
                    (t_r, t_pearson_pval) = st.pearsonr(t_x_snp, t_y_snp) #Done twice, but this is fast..
                    t_r, t_pearson_pval = float(t_r), float(t_pearson_pval)
                    t_r2 = t_r * t_r
                    y_cs.append(y_c)
                    y_ps.append(y_p)
                    r2s.append(r2)
                    t_r2s.append(t_r2)
                    ps.append(pearson_pval)
                    t_ps.append(t_pearson_pval)
                    n_saved += 1
            else:
                break
    if n_saved > 0:
        grp = h5_file.create_group('x%d' % i)
        grp.create_dataset("n_saved", data=n_saved)
        grp.create_dataset("x_c", data=x_c)
        grp.create_dataset("x_p", data=x_p)
        grp.create_dataset("x_snp", compression='gzip', data=x_snp)
        grp.create_dataset("y_cs", compression='gzip', data=y_cs)
        grp.create_dataset("y_ps", compression='gzip', data=y_ps)
        grp.create_dataset("r2s", compression='gzip', data=r2s)
        grp.create_dataset("t_r2s", compression='gzip', data=t_r2s)
        grp.create_dataset("ps", compression='gzip', data=ps)
        grp.create_dataset("t_ps", compression='gzip', data=t_ps)

    time_secs = time.time() - s1
    print 'It took %d minutes and %d seconds to finish.' % (time_secs / 60, time_secs % 60)
    print '%d values were saved.' % n_saved
    sys.stdout.flush()

    h5_file.close()




def calc_r2_levels_w_mixed_model(file_prefix, x_start_i, x_stop_i, mac_filter=10, emma_r2_threshold=0.15,
                min_emma_dist=0, save_threshold=0.15, debug_filter=1):
    """
    Returns statistics on LD levels, and plot them.
    
    Calculates emma_
    """

    dtype = 'single' #To increase matrix multiplication speed... using 32 bits.
    sd = dp.parse_numerical_snp_data(env['data_dir'] + '250K_t72.csv.binary',
                    filter=debug_filter)
    sd.filter_mac_snps(mac_filter)
    K = lm.load_kinship_from_file(env['data_dir'] + 'kinship_matrix_cm72.pickled',
                sd.accessions)
    cps_list = sd.getChrPosSNPList()
    x_cps = cps_list[x_start_i:x_stop_i]
    y_cps = cps_list
    result_list = []
    q = 1  # Single SNP is being tested
    p = 2
    n = len(sd.accessions)
    n_p = n - p
    for (x_c, x_p, x_snp) in x_cps:
        #print '%d: chromosome=%d, position=%d' % (i, x_c, x_p)
        #get the EMMA REML
        res = lm.get_emma_reml_estimates(x_snp, K)
        print 'Pseudo-heritability:', res['pseudo_heritability']
        h_sqrt_inv = res['H_sqrt_inv']
        h0_X = res['X_t']
        Y = res['Y_t']    #The transformed outputs.
        std_Y = sp.std(Y)
        (h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X, Y)
        num_emmas = 0
        for (y_c, y_p, y_snp) in y_cps:
            #if (x_c, x_p) < (y_c, y_p):
            (r, pearson_pval) = st.pearsonr(x_snp, y_snp) #Done twice, but this is fast..
            r2 = r * r
            if r2 > save_threshold:
                emma_pval = 2
                f_stat = 0
                beta = 0
                emmax_r2 = 0
                l = [x_c, x_p, y_c, y_p, r2, pearson_pval, f_stat, emma_pval, beta, emmax_r2]
                if (x_c != y_c or abs(x_p - y_p) > min_emma_dist) and r2 > emma_r2_threshold:
                    num_emmas += 1
                    #Do EMMAX
                    yt = y_snp * h_sqrt_inv
                    (b, rss, p, s) = linalg.lstsq(sp.hstack([h0_X, sp.matrix(yt).T]), Y)

                    if rss:
                        f_stat = (h0_rss / rss[0] - 1) * n_p / float(q)
                        emma_pval = st.f.sf(f_stat, q, n_p)[0]
                        beta = b[1, 0]
                        emmax_r = beta * (std_Y / sp.std(yt))
                        emmax_r2 = emmax_r * emmax_r
                if emma_pval < 0.1:
                    l = [x_c, x_p, y_c, y_p, r2, pearson_pval, f_stat[0], emma_pval, beta, emmax_r2]
                    print l
                    print std_Y, sp.std(yt)
                result_list.append(l)
        print '%d EMMAX run' % num_emmas
    file_name = file_prefix + '_x_' + str(x_start_i) + '_' + str(x_stop_i) + ".csv"
    f = open(file_name, 'w')
    for r in result_list:
        st = ','.join(map(str, r))
        f.write(st + '\n')
    f.close()
    return result_list




def run_parallel(call_method_id, x_start_i, x_stop_i, cluster='gmi'):
    """
    If no mapping_method, then analysis run is set up.
    """

    job_id = 'ld_%d_%d_%d' % (call_method_id, x_start_i, x_stop_i)
    file_prefix = '/projects/long_range_LD/raw_results/long_range_ld_min02_mac15'
    job_output_file_prefix = file_prefix + job_id

    #Cluster specific parameters    
    if cluster == 'gmi': #GMI cluster.  
        shstr = '#!/bin/bash\n'
        shstr += '#$ -S /bin/bash\n'
        shstr += '#$ -N %s\n' % job_id
        #shstr += '#$ -o %s_job_$JOB_ID.out\n' % file_prefix
        #shstr += '#$ -e %s_job_$JOB_ID.err\n' % file_prefix
        shstr += '#$ -o %s_job.out\n' % job_output_file_prefix
        shstr += '#$ -e %s_job.err\n' % job_output_file_prefix
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

    shstr += "python %sanalyze_snps_data.py %d %d %d %s" % \
            (env['script_dir'], call_method_id, x_start_i, x_stop_i, file_prefix)

    print '\n', shstr, '\n'
    script_file_name = "long_range_ld.sh"
    f = open(script_file_name, 'w')
    f.write(shstr)
    f.close()

    #Execute qsub script
    os.system("qsub " + script_file_name)



def run_r2_calc():

    call_method_id = int(sys.argv[1])
    if len(sys.argv) > 3:
        x_start_i = int(sys.argv[2])
        x_stop_i = int(sys.argv[3])
        if len(sys.argv) > 4:
            file_prefix = sys.argv[4]
        else:
            file_prefix = env['results_dir'] + 'long_range_ld_min02_mac15'
        data_format = cm_data_format_dict[call_method_id]
        calc_r2_levels(file_prefix, x_start_i, x_stop_i, call_method_id=call_method_id,
                data_format=data_format,)
    else:
        num_snps = cm_num_snps_dict[call_method_id]
        chunck_size = int(sys.argv[2])
        for i in range(0, num_snps, chunck_size):
            run_parallel(call_method_id, i, i + chunck_size)


#def _load_r2_results_(file_prefix='/storage/r2_results/250K_r2_min015'): #/Users/bjarni.vilhjalmsson/Projects/250K_r2/results/
#    headers = ['x_chr', 'x_pos', 'y_chr', 'y_pos', 'r2', 'pval', 'f_stat', 'emmax_pval', 'beta', 'emmax_r2']
#    if os.path.isfile(file_prefix + '.pickled'):
#        print 'Loading pickled data..'
#        f = open(file_prefix + '.pickled', 'rb')
#        res_dict = cPickle.load(f)
#        f.close()
#        print 'Done'
#    else:
#        sd = dp.parse_numerical_snp_data(env['data_dir'] + '250K_t72.csv.binary')
#        num_snps = len(sd.getSnps())
#        chunck_size = int(sys.argv[1])
#        res_dict = {}
#        for h in headers:
#            res_dict[h] = []
#        delim = ','
#        for i in range(0, num_snps, chunck_size):
#            file_name = file_prefix + '_x_' + str(i) + '_' + str(i + chunck_size) + ".csv"
#            print i
#            try:
#                f = open(file_name)
#                for line in f:
#                    l = map(str.strip, line.split(delim))
#                    for j, st in enumerate(l):
#                        h = headers[j]
#                        if h in ['x_chr', 'x_pos', 'y_chr', 'y_pos']:
#                            res_dict[h].append(int(st))
#                        elif h in ['pval', 'emmax_pval']:
#                            v = float(st)
#                            res_dict[h].append(v if v != 0.0 else min_float)
#                        elif h in ['r2', 'beta', 'emmax_r2']:
#                            res_dict[h].append(float(st))
#                        elif h == 'f_stat':
#                            v = float(st)
#                            res_dict[h].append(v if v != sp.nan else 0)
#                        else:
#                            raise Exception()
#            except Exception, err_str:
#                print "Problems with file %s: %s" % (file_name, err_str)
#        f = open(file_prefix + '.pickled', 'wb')
#        cPickle.dump(res_dict, f, 2)
#        f.close()
#    return res_dict

#def _load_r2_res_file_(file_name, res_dict, headers):
#    delim = ','
#    try:
#        with open(file_name) as f:
#            for line in f:
#                l = map(str.strip, line.split(delim))
#                for j, st in enumerate(l):
#                    h = headers[j]
#                    if h in ['x_chr', 'x_pos', 'y_chr', 'y_pos']:
#                        res_dict[h].append(int(st))
#                    elif h in ['pval', 't_pval']:
#                        v = float(st)
#                        res_dict[h].append(v if v != 0.0 else min_float)
#                    elif h in ['r2', 't_r2']:
#                        res_dict[h].append(float(st))
#                    else:
#                        raise Exception('Unknown value')
#    except Exception, err_str:
#        print "Problems with file %s: %s" % (file_name, err_str)


#def _load_r2_results_(file_prefix='/storage/r2_results/250K_r2_min01_mac15'):#_mac15'): #/Users/bjarni.vilhjalmsson/Projects/250K_r2/results/
#    headers = ['x_chr', 'x_pos', 'y_chr', 'y_pos', 'r2', 'pval', 't_r2', 't_pval']
#    if os.path.isfile(file_prefix + '.pickled'):
#        print 'Loading pickled data..'
#        f = open(file_prefix + '.pickled', 'rb')
#        res_dict = cPickle.load(f)
#        f.close()
#        print 'Done'
#    else:
#        sd = dp.parse_numerical_snp_data(env['data_dir'] + '250K_t72.csv.binary')
#        num_snps = len(sd.getSnps())
#        chunck_size = int(sys.argv[1])
#        res_dict = {}
#        for h in headers:
#            res_dict[h] = []
#        for i in range(0, num_snps, chunck_size):
#            file_name = file_prefix + '_x_' + str(i) + '_' + str(i + chunck_size) + ".csv"
#            _load_r2_res_file_(file_name, res_dict, headers)
#            print i
#        f = open(file_prefix + '.pickled', 'wb')
#        cPickle.dump(res_dict, f, 2)
#        f.close()
#    return res_dict


#
#def load_chr_res_dict(r2_thresholds=[(0.4, 25000), (0.2, 50000), (0.1, 100000), (0.1, 400000), (0.1, 1000000)], final_r2_thres=0.1):
#    headers = ['x_chr', 'x_pos', 'y_chr', 'y_pos', 'r2', 'pval', 't_r2', 't_pval']
#    res_dict = _load_r2_results_()
#    num_res = len(res_dict['x_chr'])
#    chromosomes = [1, 2, 3, 4, 5]
#    chr_res_dict = {}
#    for chr2 in chromosomes:
#        for chr1 in chromosomes[:chr2]:
#            d = {}
#            for h in headers:
#                d[h] = []
#            chr_res_dict[(chr1, chr2)] = d
#    num_retained = 0
#    chr_pos_set = set()
#    for i in range(num_res):
#        x_chr = res_dict['x_chr'][i]
#        y_chr = res_dict['y_chr'][i]
#        x_pos = res_dict['x_pos'][i]
#        y_pos = res_dict['y_pos'][i]
#        r2 = res_dict['t_r2'][i]
#        x_chr_pos = (x_chr, x_pos)
#        y_chr_pos = (y_chr, y_pos)
#        if x_chr <= y_chr:
#            if x_chr == y_chr and x_pos < y_pos:
#                for r2_thres, window in r2_thresholds:
#                    if y_pos - x_pos < window:
#                        if r2 > r2_thres:
#                            for h in headers:
#                                chr_res_dict[(x_chr, y_chr)][h].append(res_dict[h][i])
#                            num_retained += 1
#                            chr_pos_set.add((x_chr, x_pos))
#                            chr_pos_set.add((y_chr, y_pos))
#                        break
#                else:
#                    if r2 > final_r2_thres:
#                        for h in headers:
#                                chr_res_dict[(x_chr, y_chr)][h].append(res_dict[h][i])
#                        num_retained += 1
#                        chr_pos_set.add((x_chr, x_pos))
#                        chr_pos_set.add((y_chr, y_pos))
#            elif x_chr < y_chr:
#                if r2 > final_r2_thres:
#                    for h in headers:
#                            chr_res_dict[(x_chr, y_chr)][h].append(res_dict[h][i])
#                    num_retained += 1
#                    chr_pos_set.add((x_chr, x_pos))
#                    chr_pos_set.add((y_chr, y_pos))
#
#    print 'Number of results which were retained:', num_retained
#    print len(chr_pos_set)
#    return chr_res_dict


def load_chr_res_dict(results_prefix='/srv/lab/data/long_range_r2/swedish_seq/long_range_ld_min02_mac15_',
              final_r2_thres=0.8, final_t_r2_thres=0.2, chunk_size=500, call_method_id=78, ignore_dist=0,
              save_types=['r2s', 't_r2']):

    r2_thresholds = [(1.1, 1000000)]
#        chrom_res_dict_pickled_file = '%sfr2%0.2f_chunk%d_cm_%d.pickled' % \
#                                                (results_prefix, final_r2_thres, chunk_size, call_method_id)
#    if os.path.isfile(chrom_res_dict_pickled_file):
#                print 'Found pickled chr_res_dict, now loading it from file: %s' % chrom_res_dict_pickled_file
#                chr_res_dict = cPickle.load(open(chrom_res_dict_pickled_file))
#                print 'Pickled chr_res_dict loaded'
#                return chr_res_dict

    chrom_res_hdf5_file = '%sfr%0.2f_chunk%d_cm_%d.hdf5' % \
                    (results_prefix, final_r2_thres, chunk_size, call_method_id)

    print chrom_res_hdf5_file
    headers = ['x_pos', 'y_pos', 'r2', 'pval', 't_r2', 't_pval']
    chromosomes = [1, 2, 3, 4, 5]
    chr_res_dict = {}
    for chr2 in chromosomes:
        for chr1 in chromosomes[:chr2]:
            d = {}
            for h in headers:
                d[h] = []
            chr_res_dict[(chr1, chr2)] = d

    if os.path.isfile(chrom_res_hdf5_file):
        print 'Found hdf5 file now loading it from file: %s' % chrom_res_hdf5_file
        res_h5f = h5py.File(chrom_res_hdf5_file, 'r')

        for chr2 in chromosomes:
            for chr1 in chromosomes[:chr2]:
                h5d = res_h5f['c%d_c%d' % (chr1, chr2)]
                tup = (chr1, chr2)
                chr_res_dict[tup]['x_pos'] = h5d['x_pos'][...]
                chr_res_dict[tup]['y_pos'] = h5d['y_pos'][...]
                chr_res_dict[tup]['r2'] = h5d['r2'][...]
                chr_res_dict[tup]['t_r2'] = h5d['t_r2'][...]
                chr_res_dict[tup]['pval'] = h5d['pval'][...]
                chr_res_dict[tup]['t_pval'] = h5d['t_pval'][...]
        res_h5f.close()
        return chr_res_dict


    num_retained = 0
    chr_pos_set = set()
    num_snps = cm_num_snps_dict[call_method_id]
    for i in range(0, num_snps, chunk_size):
        result_file = '%sx_%d_%d.hdf5' % (results_prefix, i, i + chunk_size)
        if os.path.isfile(result_file):
            print 'Plowing through hdf5 file:', result_file
            h5f = h5py.File(result_file, 'r')
            #res_dict = cPickle.load(open(result_file))
        else:
            print 'Could not find hdf5 file:', result_file
            continue
        for xi in range(len(h5f.items())):
            si = 'x%d' % xi
            if si in h5f:
                d = h5f[si]
                x_pos = int(d['x_p'][...])
                x_chr = int(d['x_c'][...])
                r2s = d['r2s'][...]
                t_r2s = d['t_r2s'][...]
                af = (r2s > final_r2_thres) * (t_r2s > final_t_r2_thres) #array filter
                r2s = r2s[af]
                t_r2s = t_r2s[af]
                y_positions = map(int, d['y_ps'][af])
                y_chromosomes = map(int, d['y_cs'][af])
                pvals = d['ps'][af]
                t_pvals = d['t_ps'][af]
                n_saved = len(r2s)
                for yi in range(n_saved):
                    y_pos = y_positions[yi]
                    y_chr = y_chromosomes[yi]
                    tup = (x_chr, y_chr)
                    r2 = r2s[yi]
                    t_r2 = t_r2s[yi]
                    dist = abs(y_pos - x_pos)
                    if x_chr == y_chr and dist <= ignore_dist:
                        continue
    #                                    elif x_chr == y_chr and dist <= 5000000:
    #                                            for r2_thres, window in r2_thresholds:
    #                                if dist < window:
    #                                    if t_r2 > r2_thres:
    #                                    pval = pvals[yi]
    #                                    t_pval = t_pvals[yi]
    #                                    t_r2 = t_r2s[yi]
    #                                                                    chr_res_dict[tup]['x_pos'].append(x_pos)
    #                                                                    chr_res_dict[tup]['y_pos'].append(y_pos)
    #                                                                    chr_res_dict[tup]['r2'].append(r2)
    #                                                                    chr_res_dict[tup]['pval'].append(pval)
    #                                                                    chr_res_dict[tup]['t_r2'].append(t_r2)
    #                                                                    chr_res_dict[tup]['t_pval'].append(t_pval)
    #                                        num_retained += 1
    #                                        chr_pos_set.add(tup)
    #                                     break

                    elif x_chr != y_chr or dist > ignore_dist:
                        pval = pvals[yi]
                        t_pval = t_pvals[yi]
                        t_r2 = t_r2s[yi]
                        chr_res_dict[tup]['x_pos'].append(x_pos)
                        chr_res_dict[tup]['y_pos'].append(y_pos)
                        chr_res_dict[tup]['r2'].append(r2)
                        chr_res_dict[tup]['pval'].append(pval)
                        chr_res_dict[tup]['t_r2'].append(t_r2)
                        chr_res_dict[tup]['t_pval'].append(t_pval)
                        num_retained += 1
                        chr_pos_set.add(tup)
                del r2s
                del t_r2s
                del pvals
                del t_pvals
                del y_positions
                del y_chromosomes

        h5f.close()
        print 'Number of results which have been retained:', num_retained
    print 'Number of results which were retained:', num_retained
    print 'Number of positions involved: %d' % len(chr_pos_set)
#        print 'Pickling chr_res_dict'
#        cPickle.dump(chr_res_dict, open(chrom_res_dict_pickled_file, 'wb'), protocol=2)
#        print 'Done pickling'
    print 'Saving chr_res_dict to a HDF5 file'
    res_h5f = h5py.File(chrom_res_hdf5_file)
    for chr2 in chromosomes:
        for chr1 in chromosomes[:chr2]:
            h5d = res_h5f.create_group('c%d_c%d' % (chr1, chr2))
            tup = (chr1, chr2)
            print tup
            if len(chr_res_dict[tup]['x_pos']) > 0:
                h5d.create_dataset('x_pos', data=chr_res_dict[tup]['x_pos'])
                h5d.create_dataset('y_pos', data=chr_res_dict[tup]['y_pos'])
                h5d.create_dataset('r2', data=chr_res_dict[tup]['r2'])
                h5d.create_dataset('t_r2', data=chr_res_dict[tup]['t_r2'])
                h5d.create_dataset('pval', data=chr_res_dict[tup]['pval'])
                h5d.create_dataset('t_pval', data=chr_res_dict[tup]['t_pval'])
    res_h5f.close()
    print 'Done saving as HDF5 file'
    return chr_res_dict



#def plot_pval_emmax_correlations(filter=1.0, file_prefix='/storage/r2_results/250K_r2_min015'):
#    pickled_file = file_prefix + '_corr_info.pickled'
#    x_pvals = []
#    y_pvals = []
#    if os.path.isfile(pickled_file):
#        print 'Loading pickled data..'
#        f = open(pickled_file, 'rb')
#        d = cPickle.load(f)
#        f.close()
#        print 'Done'
#        for t in d:
#            x_pvals.append(d[t]['x'][3])
#            y_pvals.append(d[t]['y'][3])
#    else:
#        res_dict = _load_r2_results_()
#        #find pairs...
#        d = {}
#        num_res = len(res_dict['x_chr'])
#        print 'Plowing through %i results..' % num_res
#        for i in xrange(num_res):
#            if (i + 1) % (num_res / 100) == 0:
#                sys.stdout.write('.')
#                sys.stdout.flush()
#            if sp.rand() > filter: continue
#            x_chr = res_dict['x_chr'][i]
#            y_chr = res_dict['y_chr'][i]
#            x_pos = res_dict['x_pos'][i]
#            y_pos = res_dict['y_pos'][i]
#            #headers = ['x_chr', 'x_pos', 'y_chr', 'y_pos', 'r2', 'pval', 'f_stat', 'emmax_pval', 'beta', 'emmax_r2']
#            r2 = res_dict['r2'][i]
#            pval = res_dict['pval'][i]
#            f_stat = res_dict['f_stat'][i]
#            emmax_pval = res_dict['emmax_pval'][i]
#            if emmax_pval > 0.1:
#                continue
#            beta = res_dict['beta'][i]
#            emmax_r2 = res_dict['emmax_r2'][i]
#            y_t = (y_chr, y_pos)
#            x_t = (x_chr, x_pos)
#            if y_t < x_t:
#                t = (x_chr, x_pos, y_chr, y_pos)
#                flipped = True
#            else:
#                t = (y_chr, y_pos, x_chr, x_pos)
#                flipped = False
#
#            if not t in d:  #Slow as hell!!
#                d[t] = {}
#
#            if flipped:
#                d[t]['x'] = [r2, pval, f_stat, emmax_pval, beta, emmax_r2]
#            else:
#                d[t]['y'] = [r2, pval, f_stat, emmax_pval, beta, emmax_r2]
#
#        l = d.keys()[:]
#        for t in l:
#            if 'x' in d[t] and 'y' in d[t]:
#                x_emmax_pval = d[t]['x'][3]
#                y_emmax_pval = d[t]['y'][3]
#                if x_emmax_pval < 1 and y_emmax_pval < 1:
#                    x_pvals.append(x_emmax_pval)
#                    y_pvals.append(y_emmax_pval)
#                else:
#                    del d[t]
#            else:
#                del d[t]
#        f = open(pickled_file, 'wb')
#        cPickle.dump(d, f, 2)
#        f.close()
#    sp.corrcoef(x_pvals, y_pvals)[0, 1]
#    pylab.plot(x_pvals, y_pvals, '.')
#    pylab.xlabel('p-value')
#    pylab.ylabel('p-value')
#    pval_corr = sp.corrcoef(x_pvals, y_pvals)[0, 1]
#    pylab.title('Pval. corr.: %0.2f' % pval_corr)
#    pylab.savefig(env['results_dir'] + 'pval_corr_plot.png')
#    pylab.clf()
#    pylab.hexbin(x_pvals, y_pvals, gridsize=1000)
#    pylab.xlabel('p-value')
#    pylab.ylabel('p-value')
#    pylab.title('Pval. corr.: %0.2f' % pval_corr)
#    pylab.colorbar()
#    pylab.savefig(env['results_dir'] + 'pval_corr_2d_hist.png')
#
#    x_log_pvals = map(lambda x:-sp.log10(x), x_pvals)
#    y_log_pvals = map(lambda x:-sp.log10(x), y_pvals)
#    pylab.plot(x_log_pvals, y_log_pvals, '.')
#    pylab.xlabel('p-value')
#    pylab.ylabel('p-value')
#    log_pval_corr = sp.corrcoef(x_pvals, y_pvals)[0, 1]
#    pylab.title('Neg. log. pval. corr.: %0.2f' % log_pval_corr)
#    pylab.savefig(env['results_dir'] + 'log_pval_corr_plot.png')
#    pylab.clf()
#    pylab.hexbin(x_log_pvals, y_log_pvals, gridsize=1000)
#    pylab.xlabel('p-value')
#    pylab.ylabel('p-value')
#    pylab.title('Pval. corr.: %0.2f' % pval_corr)
#    pylab.colorbar()
#    pylab.savefig(env['results_dir'] + 'log_pval_corr_2d_hist.png')



def plot_r2_results(file_prefix='/srv/lab/data/long_range_r2/swedish_r2_min01_mac15_filtered', save_to_file=False):
    chrom_sizes = [30425061, 19694800, 23456476, 18578714, 26974904]
    cum_chrom_sizes = [sum(chrom_sizes[:i]) for i in range(5)]
    tot_num_bases = float(sum(chrom_sizes))
    rel_chrom_sizes = map(lambda x: 0.93 * (x / tot_num_bases), chrom_sizes)
    rel_cum_chrom_sizes = map(lambda x: 0.93 * (x / tot_num_bases), cum_chrom_sizes)
    for i in range(5):
        rel_cum_chrom_sizes[i] = rel_cum_chrom_sizes[i] + 0.01 + 0.01 * i

    chromosome_ends = {1:30.425061, 2:19.694800, 3:23.456476, 4:18.578714, 5:26.974904}
    print rel_chrom_sizes, rel_cum_chrom_sizes

    chr_res_dict = load_chr_res_dict()
    max_pval = -math.log10(min_float)
    #Filter data..
    #Now plot data!!
    alpha = 0.8
    linewidths = 0
    vmin = 0
    f = pylab.figure(figsize=(50, 46))
    chromosomes = [1, 2, 3, 4, 5]
    plot_info = [#('r2', file_prefix + '_r2s', 'Pairwise correlation ($r^2$)', 1.0),
            #('pval', file_prefix + '_pvals', 'Correlation p-value', -math.log10(min_float)),
            ('t_r2', file_prefix + '_t_r2', 'Pairwise correlation between transformed SNPs', 1.0),
            ('t_pval', file_prefix + '_t_pvals', 'Correlation p-value for pairs of transformed SNPs',
            - math.log10(min_float))]


    for h, plot_file_name, label, vmax in plot_info:
        for yi, chr2 in enumerate(chromosomes):
            for xi, chr1 in enumerate(chromosomes[:chr2]):

                ax = f.add_axes([rel_cum_chrom_sizes[xi] + 0.01, rel_cum_chrom_sizes[yi],
                        rel_chrom_sizes[xi], rel_chrom_sizes[yi] ])
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

                l = chr_res_dict[(chr1, chr2)][h]
                if h in ['pval', 't_pval']:
                    l[l[:] == 0.0] = min_float
                    l = -sp.log10(l)
                    l = l.tolist()
                l_zxy = zip(l, chr_res_dict[(chr1, chr2)]['x_pos'],
                    chr_res_dict[(chr1, chr2)]['y_pos'])
                l_zxy.sort()
                l = map(list, zip(*l_zxy))
                zs = l[0]
                xs = map(lambda x: x / 1000000.0, l[1])
                ys = map(lambda x: x / 1000000.0, l[2])

                scatter_plot = ax.scatter(xs, ys, c=zs, alpha=alpha, linewidths=linewidths, vmin=vmin,
                            vmax=vmax)
                ax.axis([-0.025 * chromosome_ends[chr1], 1.025 * chromosome_ends[chr1],
                    - 0.025 * chromosome_ends[chr2], 1.025 * chromosome_ends[chr2]])

        cax = f.add_axes([0.62, 0.3, 0.01, 0.2])
        cb = pylab.colorbar(scatter_plot, cax=cax)
        cb.set_label(label, fontsize='x-large')
        #cb.set_tick_params(fontsize='x-large')
        f.savefig(plot_file_name + '.png', format='png')

    if save_to_file:
        chromosomes = [1, 2, 3, 4, 5]
        for y_chrom in chromosomes:
            for x_chrom in chromosomes[0:y_chrom]:
                file_name = file_prefix + '_chrom%d_chrom%d_values.csv' % (x_chrom, y_chrom)
                print 'Writing to file:', file_name
                with open(file_name, 'w') as f:
                    d = chr_res_dict[(x_chrom, y_chrom)]
                    f.write('x_position, y_position, r2, t_r2\n')
                    l = zip(d['x_pos'], d['y_pos'], d['r2'], d['t_r2'])
                    l.sort()
                    for t in l:
                        f.write('%d,%d,%f,%f\n' % t)




def plot_r2_results_small(plot_file_name='/srv/lab/data/long_range_r2/r2.png',
            results_prefix='/srv/lab/data/long_range_r2/swedish_seq/long_range_ld_min02_mac15_',
            final_r2_thres=0.8, final_t_r2_thres=0.2, chunk_size=500, call_method_id=78, vmin=0.5,
            ignore_dist=1000000, symmetric=True):
    x_square_width = 0.8
    y_square_width = 0.9
    x_chrom_padding = 0.008
    y_chrom_padding = 0.0088
    chrom_sizes = [30425061, 19694800, 23456476, 18578714, 26974904]
    cum_chrom_sizes = [sum(chrom_sizes[:i]) for i in range(5)]
    tot_num_bases = float(sum(chrom_sizes))
    x_rel_chrom_sizes = map(lambda x: x_square_width * (x / tot_num_bases), chrom_sizes)
    x_rel_cum_chrom_sizes = map(lambda x: x_square_width * (x / tot_num_bases), cum_chrom_sizes)
    y_rel_chrom_sizes = map(lambda x: y_square_width * (x / tot_num_bases), chrom_sizes)
    y_rel_cum_chrom_sizes = map(lambda x: y_square_width * (x / tot_num_bases), cum_chrom_sizes)
    for i in range(5):
        x_rel_cum_chrom_sizes[i] = x_rel_cum_chrom_sizes[i] + x_chrom_padding * (i + 1)
        y_rel_cum_chrom_sizes[i] = y_rel_cum_chrom_sizes[i] + y_chrom_padding * (i + 1)

    chromosome_ends = {1:30.425061, 2:19.694800, 3:23.456476, 4:18.578714, 5:26.974904}
    print x_rel_chrom_sizes, x_rel_cum_chrom_sizes, y_rel_chrom_sizes, y_rel_cum_chrom_sizes

    chr_res_dict = load_chr_res_dict(results_prefix=results_prefix, final_r2_thres=final_r2_thres,
                    final_t_r2_thres=final_t_r2_thres, chunk_size=chunk_size,
                    call_method_id=call_method_id, save_types=['r2s', 't_r2'],
                    ignore_dist=ignore_dist)
    max_pval = -math.log10(min_float)
    #Filter data..
    #Now plot data!!
    alpha = 0.8
    linewidths = 0
    vmin = 0
    f = pylab.figure(figsize=(10, 8))
    chromosomes = [1, 2, 3, 4, 5]
#    plot_info = [('r2', file_prefix + '_r2s', 'Pairwise correlation ($r^2$)', 1.0),
#            #('pval', file_prefix + '_pvals', 'Correlation p-value', -math.log10(min_float)),
#            ('t_r2', file_prefix + '_t_r2', 'Pairwise correlation between transformed SNPs', 1.0),
#            #('t_pval', file_prefix + '_t_pvals', 'Correlation p-value for pairs of transformed SNPs', -math.log10(min_float)),
#            ]


    print 'Plotting'
    for yi, chr2 in enumerate(chromosomes):
        for xi, chr1 in enumerate(chromosomes[:chr2]):
            print 'Chromosome-pair:', (chr2, chr1)


            if chr1 == chr2:
                ax = f.add_axes([x_rel_cum_chrom_sizes[xi] + 0.05 * x_square_width,
                        y_rel_cum_chrom_sizes[yi] + 0.005 * y_square_width,
                        x_rel_chrom_sizes[xi], y_rel_chrom_sizes[yi] ])
                ax.spines['right'].set_visible(False)
                ax.spines['bottom'].set_visible(False)
                if xi > 0:
                    ax.yaxis.set_visible(False)
                    ax.spines['left'].set_visible(False)
                else:
                    ax.yaxis.set_visible(True)
                    ax.yaxis.set_ticks_position('left')
                    ax.set_ylabel('Chr. %d' % chr2)
                if yi < 4:
                    ax.xaxis.set_visible(False)
                    ax.spines['top'].set_visible(False)
                else:
                    ax.xaxis.set_visible(True)
                    ax.xaxis.set_ticks_position('top')
                    ax.xaxis.set_label_position('top')
                    ax.set_xlabel('Chr. %d' % chr1)
                    #ax.set_xlabel('Chromosome %d' % chr1)
                l = chr_res_dict[(chr1, chr2)]['r2']
                l_zxy = zip(l, chr_res_dict[(chr1, chr2)]['x_pos'],
                    chr_res_dict[(chr1, chr2)]['y_pos'])
                l_zxy.sort()
                l = map(list, zip(*l_zxy))
                zs = l[0]
                xs = map(lambda x: x / 1000000.0, l[1])
                ys = map(lambda x: x / 1000000.0, l[2])
                if symmetric:
                    l = chr_res_dict[(chr1, chr2)]['t_r2']
                    l_zxy = zip(l, chr_res_dict[(chr1, chr2)]['x_pos'],
                        chr_res_dict[(chr1, chr2)]['y_pos'])
                    l_zxy.sort()
                    l = map(list, zip(*l_zxy))
                    zs = zs + l[0]
                    xs = xs + map(lambda x: x / 1000000.0, l[2])
                    ys = ys + map(lambda x: x / 1000000.0, l[1])
                scatter_plot = ax.scatter(xs, ys, c=zs, alpha=alpha, linewidths=linewidths, vmin=vmin,
                            vmax=1.0, s=1.0, edgecolors='none')
                ax.axis([-0.025 * chromosome_ends[chr1], 1.025 * chromosome_ends[chr1],
                    - 0.025 * chromosome_ends[chr2], 1.025 * chromosome_ends[chr2]])
            else:
                ax1 = f.add_axes([x_rel_cum_chrom_sizes[xi] + 0.05 * x_square_width,
                        y_rel_cum_chrom_sizes[yi] + 0.005 * y_square_width,
                        x_rel_chrom_sizes[xi], y_rel_chrom_sizes[yi] ])
                ax1.spines['right'].set_visible(False)
                ax1.spines['bottom'].set_visible(False)
                if xi > 0:
                    ax1.yaxis.set_visible(False)
                    ax1.spines['left'].set_visible(False)
                else:
                    ax1.yaxis.set_visible(True)
                    ax1.yaxis.set_ticks_position('left')
                    ax1.set_ylabel('Chr. %d' % chr2)
                if yi < 4:
                    ax1.xaxis.set_visible(False)
                    ax1.spines['top'].set_visible(False)
                else:
                    ax1.xaxis.set_visible(True)
                    ax1.xaxis.set_ticks_position('top')
                    ax1.xaxis.set_label_position('top')
                    ax1.set_xlabel('Chr. %d' % chr1)
                l = chr_res_dict[(chr1, chr2)]['r2']
                l_zxy = zip(l, chr_res_dict[(chr1, chr2)]['x_pos'],
                    chr_res_dict[(chr1, chr2)]['y_pos'])
                l_zxy.sort()
                l = map(list, zip(*l_zxy))
                zs = l[0]
                xs = map(lambda x: x / 1000000.0, l[1])
                ys = map(lambda x: x / 1000000.0, l[2])
                scatter_plot = ax1.scatter(xs, ys, c=zs, alpha=alpha, linewidths=linewidths, vmin=vmin,
                            vmax=1.0, s=1.0, edgecolors='none')
                ax1.axis([-0.025 * chromosome_ends[chr1], 1.025 * chromosome_ends[chr1],
                    - 0.025 * chromosome_ends[chr2], 1.025 * chromosome_ends[chr2]])

                #Now the opposite side
                if symmetric:
                    ax2 = f.add_axes([x_rel_cum_chrom_sizes[yi] + 0.05 * x_square_width,
                            y_rel_cum_chrom_sizes[xi] + 0.005 * y_square_width,
                            x_rel_chrom_sizes[yi], y_rel_chrom_sizes[xi] ])
                    ax2.spines['right'].set_visible(False)
                    ax2.spines['bottom'].set_visible(False)
                    ax2.spines['left'].set_visible(False)
                    ax2.yaxis.set_visible(False)
                    ax2.spines['top'].set_visible(False)
                    ax2.xaxis.set_visible(False)
                    l = chr_res_dict[(chr1, chr2)]['t_r2']
                    l_zxy = zip(l, chr_res_dict[(chr1, chr2)]['x_pos'],
                        chr_res_dict[(chr1, chr2)]['y_pos'])
                    l_zxy.sort()
                    l = map(list, zip(*l_zxy))
                    zs = l[0]
                    xs = map(lambda x: x / 1000000.0, l[2])
                    ys = map(lambda x: x / 1000000.0, l[1])
                    scatter_plot = ax2.scatter(xs, ys, c=zs, alpha=alpha, linewidths=linewidths, vmin=vmin,
                                vmax=1.0, s=1.0, edgecolors='none')
                    ax2.axis([-0.025 * chromosome_ends[chr2], 1.025 * chromosome_ends[chr2],
                        - 0.025 * chromosome_ends[chr1], 1.025 * chromosome_ends[chr1]])

    cax = f.add_axes([0.91, 0.15, 0.03, 0.7])
    cb = pylab.colorbar(scatter_plot, cax=cax)
    cb.set_label(r'Pairwise $r^2$', fontsize='x-large')
    #cb.set_tick_params(fontsize='x-large')
    f.savefig(plot_file_name)



def plot_gw_r2_decay(file_prefix, num_random_xs=20000, max_dist=1000000, call_method_id=78, mac_filter=18, debug_filter=0.5, n_bins=40):
    """
    Plots r2 decay on the genome-wide scale
    """
    dtype = 'single' #To increase matrix multiplication speed... using 32 bits.
    sd = dp.load_snps_call_method(call_method_id=call_method_id, debug_filter=debug_filter, min_mac=mac_filter)
    sd.filter_mac_snps(mac_filter)
    h_inverse_matrix_file = env['data_dir'] + 'snp_cov_mat_h_inv_cm%d.pickled' % (call_method_id)
    if not os.path.isfile(h_inverse_matrix_file):
        K = sd.get_snp_cov_matrix()
        H_sqrt = lm.cholesky(K)
        H_sqrt_inv = (H_sqrt).I
        with file(h_inverse_matrix_file, 'wb') as f:
            cPickle.dump(H_sqrt_inv, f, protocol=2)
    else:
        with file(h_inverse_matrix_file) as f:
            H_sqrt_inv = cPickle.load(f)

    cps_list = sd.getChrPosSNPList()
    x_cps = random.sample(cps_list, num_random_xs)
    y_cps = cps_list
    result_dict = {}
    n = float(len(sd.accessions))
    print 'Starting calculation'
    sys.stdout.flush()
    dists = []
    r2s = []
    t_r2s = []
    x_mafs = []
    y_mafs = []
    n_saved = 0
    s1 = time.time()
    for i, (x_c, x_p, x_snp) in enumerate(x_cps):
        print '%d: chromosome=%d, position=%d' % (i, x_c, x_p)
        #Normalize SNP..
        xs = sp.array(x_snp)
        x_mac = sum(xs)
        x_maf = x_mac / n
        if x_maf < 0.1:
            continue
        t_x_snp = sp.dot(((xs - sp.mean(xs)) / sp.std(xs)), H_sqrt_inv).T
        for (y_c, y_p, y_snp) in reversed(y_cps):
            if x_c != y_c:
                continue
            if abs(x_p - y_p) > max_dist:
                continue
            ys = sp.array(y_snp)
            y_maf = sum(ys) / n
            if y_maf < 0.1:
                continue
            x_i = bisect([0.0, 0.1, 0.2, 0.3, 0.4], x_maf)
            y_i = bisect([0.0, 0.1, 0.2, 0.3, 0.4], y_maf)
            if x_i != y_i:
                continue
            x_mafs.append(x_maf)
            y_mafs.append(y_maf)
            (r, pearson_pval) = st.pearsonr(xs, ys)
            r2 = r * r
            t_y_snp = sp.dot(((ys - sp.mean(ys)) / sp.std(ys)), H_sqrt_inv).T
            (t_r, t_pearson_pval) = st.pearsonr(t_x_snp, t_y_snp) #Done twice, but this is fast..
            t_r, t_pearson_pval = float(t_r), float(t_pearson_pval)
            t_r2 = t_r * t_r
            dists.append(abs(x_p - y_p))
            r2s.append(r2)
            t_r2s.append(t_r2)
            n_saved += 1


    time_secs = time.time() - s1
    print 'It took %d minutes and %d seconds to finish.' % (time_secs / 60, time_secs % 60)
    print '%d values were saved.' % n_saved
    sys.stdout.flush()

    #Now binning by MAC..
    mac_dict = {}
    for f in [1, 2, 3, 4, 5]: #For 0.1,0.2,0.3,0.4,0.5 MAF thresholds
        mac_dict[f] = {'t_r2s':[], 'r2s':[], 'dists':[]}
    n = float(n)
    left_overs = 0
    for dist, r2, t_r2, x_maf, y_maf in izip(dists, r2s, t_r2s, x_mafs, y_mafs):
        x_i = bisect([0.0, 0.1, 0.2, 0.3, 0.4], x_maf)
        y_i = bisect([0.0, 0.1, 0.2, 0.3, 0.4], y_maf)
        if x_i == y_i:
            d = mac_dict[x_i]
            d['r2s'].append(r2)
            d['t_r2s'].append(t_r2)
            d['dists'].append(dist)
            left_overs += 1
    print 'After MAF level filtering %d pairs were left.' % left_overs

    #Now plotting and binning..
    for m_dist in [50000, 100000, 200000, 500000, 1000000]:
        kbs = m_dist / 1000
        filtered_r2s = []
        filtered_t_r2s = []
        filtered_dists = []
        bin_dict = {}
        for bid in range(n_bins):
            bin_dict[bid] = {'r2s': [], 't_r2s': []}
        for f in [2, 3, 4, 5]:
            d = mac_dict[f]
            bin_ids = sp.digitize(d['dists'], sp.arange(0, m_dist, m_dist / n_bins)) - 1
            for bid, r2, t_r2, dist in izip(bin_ids, d['r2s'], d['t_r2s'], d['dists']):
                if dist > m_dist:
                    continue
                bin_dict[bid]['r2s'].append(r2)
                filtered_r2s.append(r2)
                bin_dict[bid]['t_r2s'].append(t_r2)
                filtered_t_r2s.append(t_r2)
                filtered_dists.append(dist)


        pylab.figure()
        pylab.plot(filtered_dists, filtered_r2s, alpha=0.3, color='k', marker='.', ls='None')
        pylab.xlabel('Distance (bases)')
        pylab.ylabel(r'$r^2$')
        pylab.savefig(file_prefix + '_%dkb_r2s.png' % (kbs))
        pylab.figure()
        pylab.plot(filtered_dists, filtered_t_r2s, alpha=0.3, color='k', marker='.', ls='None')
        pylab.xlabel('Distance (bases)')
        pylab.ylabel(r'$r^2$')
        pylab.savefig(file_prefix + '_%dkb_t_r2s.png' % (kbs))


        r2_avgs = []
        t_r2_avgs = []
        xs = []
        l = sp.arange(0, m_dist, m_dist / n_bins) + (m_dist / (2 * n_bins))
        for bid in range(n_bins):
            n = float(len(bin_dict[bid]['r2s']))
            if n > 0:
                r2_avgs.append(sp.sum(bin_dict[bid]['r2s']) / n)
                t_r2_avgs.append(sp.sum(bin_dict[bid]['t_r2s']) / n)
                xs.append(l[bid])

        pylab.figure()
        pylab.plot(xs, r2_avgs, alpha=0.7, color='b', lw=1.8, label=r'standard $r^2$')
        pylab.plot(xs, t_r2_avgs, alpha=0.7, color='m', lw=1.8, label=r'transformed $r^2$')
        pylab.legend(loc=1)
        pylab.xlabel('Distance (bases)')
        pylab.ylabel(r'$r^2$')
        pylab.savefig(file_prefix + '_%dkb_r2s_avgs.png' % (kbs))



#
#def plot_r2_results():
#    chr_res_dict = load_chr_res_dict()
#    max_pval = -math.log10(min_float)
#    max_emmax_pval = -math.log10(min(res_dict['emmax_pval']))
#    #Filter data..
#    #Now plot data!!
#    alpha = 0.8
#    linewidths = 0
#    vmin = 0.0
#    for chr2 in chromosomes:
#        for chr1 in chromosomes[:chr2]:
#            x_min = min(chr_res_dict[(chr1, chr2)]['x_pos'])
#            x_max = max(chr_res_dict[(chr1, chr2)]['x_pos'])
#            y_min = min(chr_res_dict[(chr1, chr2)]['y_pos'])
#            y_max = max(chr_res_dict[(chr1, chr2)]['y_pos'])
#            x_range = x_max - x_min
#            y_range = y_max - y_min
#            x_lab = 'Chromosome %d' % chr1
#            y_lab = 'Chromosome %d' % chr2
#            left = 0.05
#            bottom = 0.04
#            width = 0.94
#            height = 0.94
#            r2_plot_file_name = file_prefix + 'c_' + str(chr1) + 'x' + str(chr2) + '_r2s.png'
#            pval_file_name = file_prefix + 'c_' + str(chr1) + 'x' + str(chr2) + '_pvals.png'
#            emma_pval_file_name = file_prefix + 'c_' + str(chr1) + 'x' + str(chr2) + '_emmax_pvals.png'
##            print r2_plot_file_name, pval_file_name , emma_pval_file_name
#            pylab.figure(figsize=(18, 16))
#            pylab.axes([left, bottom, width, height])
#
#            l_zxy = zip(chr_res_dict[(chr1, chr2)]['r2'], chr_res_dict[(chr1, chr2)]['x_pos'],
#                chr_res_dict[(chr1, chr2)]['y_pos'])
#            l_zxy.sort()
#            l = map(list, zip(*l_zxy))
#            zs = l[0]
#            xs = l[1]
#            ys = l[2]
##            print len(chr_res_dict[(chr1, chr2)]['x_pos']), len(chr_res_dict[(chr1, chr2)]['y_pos']), \
##                len(chr_res_dict[(chr1, chr2)]['r2'])
#            pylab.scatter(xs, ys, c=zs, alpha=alpha, linewidths=linewidths, vmin=vmin, vmax=1.0)
#            pylab.axis([x_min - 0.025 * x_range, x_max + 0.025 * x_range, y_min - 0.025 * y_range,
#                y_max + 0.025 * y_range])
#            pylab.colorbar()
#            pylab.xlabel(x_lab)
#            pylab.ylabel(y_lab)
#            pylab.savefig(r2_plot_file_name, format='png')
#            pylab.clf()
#            pylab.figure(figsize=(18, 16))
#            pylab.axes([left, bottom, width, height])
#            log_pvals = map(lambda x:-math.log10(x), chr_res_dict[(chr1, chr2)]['pval'])
#            l_zxy = zip(log_pvals, chr_res_dict[(chr1, chr2)]['x_pos'], chr_res_dict[(chr1, chr2)]['y_pos'])
#            l_zxy.sort()
#            l = map(list, zip(*l_zxy))
#            zs = l[0]
#            xs = l[1]
#            ys = l[2]
#            pylab.scatter(xs, ys, c=zs, alpha=alpha, linewidths=linewidths, vmin=vmin, vmax=max_pval)
#            pylab.axis([x_min - 0.025 * x_range, x_max + 0.025 * x_range, y_min - 0.025 * y_range,
#                y_max + 0.025 * y_range])
#            pylab.colorbar()
#            pylab.xlabel(x_lab)
#            pylab.ylabel(y_lab)
#            pylab.savefig(pval_file_name, format='png')
#            pylab.clf()
#            pylab.figure(figsize=(18, 16))
#            pylab.axes([left, bottom, width, height])
#            log_pvals = map(lambda x:-math.log10(x), chr_res_dict[(chr1, chr2)]['emmax_pval'])
#            l_zxy = zip(log_pvals, chr_res_dict[(chr1, chr2)]['x_pos'], chr_res_dict[(chr1, chr2)]['y_pos'])
#            l_zxy.sort()
#            i = 0
#            while l_zxy[i] < 0:
#                i += 1
#            l_zxy = l_zxy[i:]
#            l = map(list, zip(*l_zxy))
#            zs = l[0]
#            xs = l[1]
#            ys = l[2]
#            pylab.scatter(xs, ys, c=zs, alpha=alpha, linewidths=linewidths, vmin=vmin, vmax=max_emmax_pval)
#            pylab.axis([x_min - 0.025 * x_range, x_max + 0.025 * x_range, y_min - 0.025 * y_range,
#                y_max + 0.025 * y_range])
#            pylab.colorbar()
#            pylab.xlabel(x_lab)
#            pylab.ylabel(y_lab)
#            pylab.savefig(emma_pval_file_name, format='png')
#plt.subplot(121)
#plt.scatter(xyc[:13], xyc[:13], c=xyc[:13], s=35, vmin=0, vmax=20)
#plt.colorbar()
#plt.xlim(0, 20)
#plt.ylim(0, 20)
#
#plt.subplot(122)
#plt.scatter(xyc[8:20], xyc[8:20], c=xyc[8:20], s=35, vmin=0, vmax=20)   
#plt.colorbar()
#plt.xlim(0, 20)
#plt.ylim(0, 20)









if __name__ == "__main__":
    #run_r2_calc()
    #plot_gw_r2_decay(env['results_dir'] + 'ld_cm75', num_random_xs=5000, call_method_id=75,
    #            debug_filter=0.2)
    #plot_r2_results(save_to_file=True)
    #plot_r2_results_small(final_r2_thres=0.6, final_t_r2_thres=0.0, vmin=0.0)
    #plot_pval_emmax_correlations()
    #test_correlation()
    plot_gw_r2_decay('/tmp/r2_test')
