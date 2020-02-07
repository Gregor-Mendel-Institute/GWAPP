"""
Basic analysis of Metagenomic data.

Matt Horton (mhorton@uchicago.edu)
(hack/port of code from linear_models.py to account for reps.).
"""

#TODO
# - run on cluster (using imputed SNP data):
# - 1.) synchronize phenotype data/snp data on ecotype id.
# - 2.) thin snp data.
# - 3.) perform emma-x on transformed phenotype data
# - 4.) make manhattan plots; write p-values out to file.

from env import *
import sys
import time

import phenotypeData as pd
import gwaResults as gr
import genotypes as gt
import kinship as kin

import linear_models as lm
import dataParsers as dp

import math
import scipy as sp
from scipy import linalg
from scipy import stats

class PermuteGwas:

        def _emmax_permutations(self, snps, phenotypes, num_perm, K=None, Z=None, method='REML'):
                """
                EMMAX permutation test
                Single SNPs
                
                Returns the list of max_pvals and max_fstats 
                """
                lmm = lm.LinearMixedModel(phenotypes)
                lmm.add_random_effect( Z * K * Z.T )

                eig_L = lmm._get_eigen_L_()

                print 'Getting variance estimates'
                res = lmm.get_estimates(eig_L, method=method)

                q = 1  # Single SNP is being tested
                p = len(lmm.X.T) + q
                n = lmm.n
                n_p = n - p
                H_sqrt_inv = res['H_sqrt_inv']

                Y = H_sqrt_inv * lmm.Y    #The transformed outputs.
                h0_X = H_sqrt_inv * lmm.X
                (h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X, Y)
                Y = Y - h0_X * h0_betas

                num_snps = len(snps)
                max_fstat_list = []
                min_pval_list = []
                chunk_size = len(Y)
                print "Working with chunk size: " + str(chunk_size)
                print "and " + str(num_snps) + " SNPs."
                Ys = sp.mat(sp.zeros((chunk_size, num_perm)))

                for perm_i in range(num_perm):
                    #print 'Permutation nr. % d' % perm_i
                    sp.random.shuffle(Y)
                    Ys[:, perm_i] = Y

                min_rss_list = sp.repeat(h0_rss, num_perm)
                for i in range(0, num_snps, chunk_size): #Do the dot-product in chunks!
                    snps_chunk = sp.matrix(snps[i:(i + chunk_size)])
                    snps_chunk = snps_chunk * Z.T
                    Xs = snps_chunk * (H_sqrt_inv.T)
                    Xs = Xs - sp.mat(sp.mean(Xs, axis=1))
                    for j in range(len(Xs)): # for each snp
                        (betas, rss_list, p, sigma) = linalg.lstsq(Xs[j].T, Ys, overwrite_a=True) # read the lstsq lit
                        for k, rss in enumerate(rss_list):
                            if not rss:
                                print 'No predictability in the marker, moving on...'
                                continue
                            if min_rss_list[k] > rss:
                                min_rss_list[k] = rss
                        if num_snps >= 10 and (i + j + 1) % (num_snps / num_perm) == 0: #Print dots
                            sys.stdout.write('.')
                            sys.stdout.flush()

                if num_snps >= 10:
                    sys.stdout.write('\n')

                #min_rss = min(rss_list)
                max_f_stats = ((h0_rss / min_rss_list) - 1.0) * n_p / float(q)
                min_pvals = (stats.f.sf(max_f_stats, q, n_p))

                res_d = {'min_ps':min_pvals, 'max_f_stats':max_f_stats}
                print "There are: " + str(len(min_pvals))
                return res_d

        def _run_otu_wperm(self,
                        file_prefix, 
                        phenotype_file, 
                        delimiter=',', 
                        covariate_file=None, 
                        phenotype_id=1,
                        call_method_id=1307, 
                        maf_threshold=5,
                        number_of_permutations=10):
##
#                phenotype_file = "/home/GMI/matt.horton/meta/metagenomics/gwas/leaf/16S/min800_cca/phenotypes/leaf.16S.800.2sampPerOTU.rare.cca.abd.2reps.n100.cca.txt"
#                call_method_id = 1308
#                maf_threshold = 5
#                phenotype_id = 1
#                delimiter = ','

                print "Opening snp and phenotype files."
                sys.stdout.flush()
        
                if '/' in phenotype_file:
                        print "Opening phenotype-file: " + phenotype_file
                        phenotype = pd.parse_phenotype_file(phenotype_file, delim=delimiter)  #load phenotype file
                        results_directory = phenotype_file.partition("phenotypes") # parse this off of the phenotypeFileName and sub the phenotypes dir for the results dir (which needs to be at the same level!!!)
                        results_directory = results_directory[0] + 'results/'
                        print "Outputing results to: " + results_directory
                else:
                        phenotype = pd.parse_phenotype_file(env['phen_dir'] + phenotype_file, delim=delimiter)  #load phenotype file
                        results_directory = env['results_dir']
        
                sd = dp.load_snps_call_method(call_method_id=call_method_id, data_format='binary')
                indices_to_keep = sd.coordinate_w_phenotype_data(phenotype, phenotype_id) #truncate to the phenotype of interest.
                indices_to_keep = indices_to_keep.get('pd_indices_to_keep')

                # determine whether to use mac or maf (I might have to use the mac code after determining what the mac should be from the maf)
                if maf_threshold > 0:
                        sd.filter_mac_snps(10);
#                        mac_threshold = int(math.ceil(len(sd.accessions) * (float(maf_threshold) / 100)))
#                        print "Applying maf threshold: " + str(maf_threshold) + "% to " + str(len(sd.accessions)) + " accessions (mac < " + str(mac_threshold) + ")"
#                        sd.filter_mac_snps(mac_threshold)

                phenotype_name = phenotype.get_name(phenotype_id)
                phenotype_values = phenotype.get_values(phenotype_id)
                Z = phenotype.get_incidence_matrix( phenotype_id )

                print "There are: " + str(sd.num_snps()) + " SNPs."
                print "in: " + str(len(sd.accessions)) + " accessions"
                print "and " + str(len(indices_to_keep)) + " observations."
                print "The average number of observations per genotype is " + str(float(len(indices_to_keep)) / float(len(sd.accessions)))
                sys.stdout.flush()

                K = sd.get_ibs_kinship_matrix()
                K = sp.matrix(K)
                Z = sp.matrix(Z)

                print "Examining phenotype: '" + phenotype_name + "' (phenotype_id: " + str(phenotype_id) + ")."
                print 'Applying Permutation tests.'

                snps = sd.get_snps()

                print "Running %d EMMAX-permutations (writes %d dots)" % (number_of_permutations, number_of_permutations)
                s1 = time.time()
                res_perm = self._emmax_permutations(snps, phenotype_values, number_of_permutations, K=K, Z=Z)
                p_f_list = zip(res_perm['min_ps'], res_perm['max_f_stats'])
                p_f_list.sort()
                print p_f_list[:10]
                threshold = p_f_list[len(p_f_list) / 20]
                res_perm['threshold_05'] = threshold
                print 'Threshold should be:', threshold
                secs = time.time() - s1
                if secs > 60:
                    mins = int(secs) / 60
                    secs = secs - mins * 60
                    print 'Took %d mins and %f seconds.' % (mins, secs)
                else:
                    print 'Took %f seconds.' % (secs)

                print "Permutation tests done for phenotype: " + phenotype_name

                results = {}
                results['perm_pval'] = res_perm['min_ps'].tolist()
                results['perm_fstat'] = res_perm['max_f_stats'].tolist()

                output_file = '%s/%s_perm.pvals_pid_%d_%s' % ( results_directory, file_prefix, phenotype_id, phenotype_name)
                columns = ['perm_pval', 'perm_fstat']
                with open(output_file, "w") as f:
                        f.write(','.join(columns) + "\n")
                        for i in range(1, (number_of_permutations + 1)):
                                l = [results[c][i-1] for c in columns]
                                l = map(str, l)
                                f.write(",".join(l) + "\n")
                print "Permutation p-values written."

if __name__ == '__main__': 
        pm = PermuteGwas()
        pm._run_otu_wperm(  file_prefix = sys.argv[1], 
                        phenotype_file = sys.argv[2], 
                        phenotype_id = int(sys.argv[3]),
                        call_method_id = int(sys.argv[4]),
                        maf_threshold = int(sys.argv[5]),
                        number_of_permutations = int(sys.argv[6]));

        print "We are now finished running permutation tests for: " + str(sys.argv[3])