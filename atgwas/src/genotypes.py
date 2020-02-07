"""
A module which contains a updated version of genotype data classes
"""
import h5py
import scipy as sp
from scipy import linalg
from scipy import stats
import sys
import bisect
#import kinship
import linear_models
import time

class genotype_data:
    """
    A class which uses HDF5 to store SNPs, but otherwise implements a 
    similar interface as the the older SNPsDataSet.    
    """
    def __init__(self, hdf5_file_name, sd=None, hdf5_chunk_size=100, filling_chunk_size=2000, file_access_mode='a'):

        """
        If sd is provided, then the hdf5_file_name is filled with that data. 
        """
        self.hdf5_file_name = hdf5_file_name
        if sd is not None:
            self.h5file = h5py.File(hdf5_file_name, file_access_mode)

            #Fill file
            print 'Filling file'
            self.h5file.create_dataset('indiv_ids', data=sd.accessions) #i.e. ecotypes
            self.h5file.create_dataset('num_indivs', data=len(sd.accessions))
            self.h5file.create_dataset('chromosomes', data=sd.get_chr_list(), compression='lzf')
            self.h5file.create_dataset('positions', data=sd.get_positions(), compression='lzf')
            num_snps = sd.num_snps()
            self.h5file.create_dataset('num_snps', data=num_snps)
            self.h5file.create_dataset('data_format', data=sp.array(sd.data_format))
            if sd.data_format in ['binary', 'diploid_int']:
                self.h5file.create_dataset('snps', shape=(sd.num_snps(), len(sd.accessions)),
                                dtype='int8', compression='lzf', chunks=((hdf5_chunk_size, len(sd.accessions))))
                offset = 0
                for snpsd in sd.snpsDataList:
                    n_snps = len(snpsd.snps)
                    for i in range(0, n_snps, filling_chunk_size):
                        sys.stdout.write('\b\b\b\b\b%.1f%%' % (100 * (float(offset + (1 + i))) / \
                                    num_snps))
                        sys.stdout.flush()
                        stop_i = min(i + filling_chunk_size, n_snps)
                        snps_chunk = sp.array(snpsd.snps[i:stop_i], dtype='int8')
                        self.h5file['snps'][offset + i: offset + stop_i ] = snps_chunk
                    offset += n_snps
                sys.stdout.write('\b\b\b\b\b100.0%\n')
                sys.stdout.flush()
            else:
                raise NotImplementedError

            self.h5file.create_group('filters')
            self.h5file.close()
            self.h5file = h5py.File(hdf5_file_name, 'r')

        else:
            self.h5file = h5py.File(hdf5_file_name, 'r')

        self.data_format = str(self.h5file['data_format'][...])
        self.indiv_filter = None
        self.snps_filter = None
        self.cached_snps = None


    def num_individs(self):
        if self.indiv_filter is None:
            return int(self.h5file['num_indivs'][...])
        else:
            return len(self.indiv_filter)


    def num_snps(self):
        if self.snps_filter is None:
            return int(self.h5file['num_snps'][...])
        else:
            return int(sp.sum(self.snps_filter))

    def _get_cached_group_(self):

        if self.indiv_filter is None:
            cache_tuple = 'full_data'
        else:
            cache_tuple = str(tuple(self.h5file['indiv_ids'][self.indiv_filter]))
        if cache_tuple in self.h5file['filters'].keys():
            g = self.h5file['filters'][cache_tuple]
            return g, True
        else:
            g = self.h5file['filters'].create_group(cache_tuple)
            return g, False



#    def _update_macs_(self, g, chunk_size=1024):
#        temp_filter = self.snps_filter
#        self.snps_filter = None
#        n_snps = self.num_snps()
#        g.create_dataset('macs', shape=(n_snps,), compression='gzip')
#        print 'Calculating MACs'
#        offset = 0
#        snps_chunks = self.snps_chunks(chunk_size)
#        for chunk_i, snps_chunk in enumerate(snps_chunks):
#            a = sp.empty(len(snps_chunk))
#            if self.data_format == 'binary':
#                for j, snp in enumerate(snps_chunk):
#                    bc = sp.bincount(snp)
#                    a[j] = 0 if len(bc) < 2 else bc.min()
#            else:
#                for j, snp in enumerate(snps_chunk):
#                    bc = sp.bincount(snp)
#                    if len(bc) < 2:
#                        a[j] = 0
#                    else:
#                        l = sp.array([bc[0], bc[2]]) + bc[1] / 2.0
#                        a[j] = l.min()
#            g['macs'][offset:offset + len(snps_chunk)] = a
#            offset += len(snps_chunk)
#            sys.stdout.write('\b\b\b\b\b\b%0.2f%%' % (100.0 * (min(1, \
#                                    ((chunk_i + 1.0) * chunk_size) / n_snps))))
#            sys.stdout.flush()
#        self.snps_filter = temp_filter
#        print '\nFinished calculating MACs'



#    def filter_mac(self, min_mac=15):
#        """
#        Sets a mac filter which is applied at runtime.
#        """
#        #Check wether cached! otherwise..
#        g, already_exists = self._get_cached_group_()
#        if not already_exists:
#            self._update_macs_(g)
#        if self.snps_filter != None:
#            self.snps_filter = self.snps_filter * (g['macs'][...] >= min_mac)
#        else:
#            self.snps_filter = g['macs'][...] >= min_mac


#    def _get_macs_(self, chunk_size=1024):
#        temp_filter = self.snps_filter
#        self.snps_filter = None
#        n_snps = self.num_snps()
#        print 'Calculating MACs'
#        offset = 0
#        snps_chunks = self.snps_chunks(chunk_size)
#        for chunk_i, snps_chunk in enumerate(snps_chunks):
#            a = sp.empty(len(snps_chunk))
#            if self.data_format == 'binary':
#                for j, snp in enumerate(snps_chunk):
#                    bc = sp.bincount(snp)
#                    a[j] = 0 if len(bc) < 2 else bc.min()
#            else:
#                for j, snp in enumerate(snps_chunk):
#                    bc = sp.bincount(snp)
#                    if len(bc) < 2:
#                        a[j] = 0
#                    else:
#                        l = sp.array([bc[0], bc[2]]) + bc[1] / 2.0
#                        a[j] = l.min()
#            g['macs'][offset:offset + len(snps_chunk)] = a
#            offset += len(snps_chunk)
#            sys.stdout.write('\b\b\b\b\b\b%0.2f%%' % (100.0 * (min(1, \
#                                    ((chunk_i + 1.0) * chunk_size) / n_snps))))
#            sys.stdout.flush()
#        self.snps_filter = temp_filter
#        print '\nFinished calculating MACs'
#
#
#    def filter_mac(self, min_mac=15):
#        """
#        Sets a mac filter which is applied at runtime.
#        """
#        #Check wether cached! otherwise..
#        macs = self._get_macs_()
#        if self.snps_filter != None:
#            self.snps_filter = self.snps_filter * (macs >= min_mac)
#        else:
#            self.snps_filter = macs >= min_mac


    def filter_indivs(self, ids):
        """
        Filter individuals, leaving the ones given as an argument.
        """
        raise NotImplementedError
#        indiv_ids = self.h5file['indiv_ids'][...]
#        indiv_filter = sp.zeros(len(indiv_ids))
#        for i in ids:
#            pass




    def filter_random_snps(self, fraction):
        """
        Filters random SNPs, which may be helpful when debugging
        """
        if fraction < 1.0:
            num_snps = self.h5file['num_snps'][...]
            rand_floats = sp.random.random(num_snps)
            snps_filter = sp.zeros(num_snps)
            snps_filter[rand_floats < fraction] = 1
            if self.snps_filter is None:
                self.snps_filter = snps_filter
            else:
                self.snps_filter = self.snps_filter * snps_filter

    def sync_phenotype_data_w_snps(self, pheno_i, pid):
        """
        Onus is on the user to use genotypes.genotype_data to reduce data to snp-set and ecotype-set of interest.
        This function reduces the phenotype set to the desired snp-set.
        """

        print "Synchronizing Phenotype data with SNP data."
        ets = pheno_i.phen_dict[pid]['ecotypes']

        #Checking which accessions to keep and which to remove from the phenotype file.
        pd_indices_to_keep = []

        #Filter accessions which do not have phenotype values (from the genotype data).
        for i, iid in enumerate(self.h5file['indiv_ids']):
            for j, et in enumerate(ets):
                if et == iid:
                    pd_indices_to_keep.append(j)

        num_values = len(pheno_i.phen_dict[pid]['ecotypes'])
        print "Filtering phenotype data."
        pheno_i.filter_ecotypes(pd_indices_to_keep, pids=[pid]) #Removing accessions that don't have genotypes or phenotype values
        ets = pheno_i.phen_dict[pid]['ecotypes']
        print "Out of %d, leaving %d values." % (num_values, len(ets))

        return pd_indices_to_keep

#    def coordinate_w_phenotype_data(self, phend, pid, coord_phen=True):
#
#        """
#        Deletes accessions which are not common, and sorts the accessions, removes monomorphic SNPs, etc.
#        """
#        print "Coordinating SNP and Phenotype data."
#        ets = phend.phen_dict[pid]['ecotypes']
#
#        #Checking which accessions to keep and which to remove.
#        sd_indices_to_keep = set()#[]
#        pd_indices_to_keep = []
#
#        for i, iid in enumerate(self.h5file['indiv_ids']):
#            for j, et in enumerate(ets):
#                if et == iid:
#                    sd_indices_to_keep.add(i)
#                    pd_indices_to_keep.append(j)
#
#        sd_indices_to_keep = list(sd_indices_to_keep)
#        sd_indices_to_keep.sort()
#
#
#        #Filter accessions which do not have phenotype values (from the genotype data).
#        print "Filtering genotype data"
#        #if len(sd_indices_to_keep) != len(self.accessions):
#        self.indiv_filter = sd_indices_to_keep
#        if coord_phen:
#            num_values = len(phend.phen_dict[pid]['ecotypes'])
#            print "Filtering phenotype data."
#            phend.filter_ecotypes(pd_indices_to_keep, pids=[pid]) #Removing accessions that don't have genotypes or phenotype values
#            ets = phend.phen_dict[pid]['ecotypes']
#            print "Out of %d, leaving %d values." % (num_values, len(ets))
#
#        return pd_indices_to_keep


    def get_snps(self, chunk_size=1000):
        n_snps = self.num_snps()
        n_indivs = self.num_individs()
        if self.snps_filter is None and self.indiv_filter is None:
            snps = self.h5file['snps'][...]
        else:
            if self.data_format in ['binary', 'diploid_int']:
                print 'Allocating memory'
                snps = sp.empty((n_snps, n_indivs), dtype='int8')
                print 'done allocating.'
            else:
                raise NotImplementedError
            offset = 0
            print 'Extracting the SNPs'
            for chunk_i, snps_chunk in enumerate(self.snps_chunks(chunk_size)):
                snps[offset:offset + len(snps_chunk)] = snps_chunk
                offset += len(snps_chunk)
                sys.stdout.write('\b\b\b\b\b\b%0.2f%%' % (100.0 * (min(1, \
                                    ((chunk_i + 1.0) * chunk_size) / n_snps))))
                sys.stdout.flush()

            print '\nDone extracting the SNPs.'
        return snps



    def get_positions(self):
        if self.snps_filter is None:
            return self.h5file['positions'][...]
        else:
            return self.h5file['positions'][self.snps_filter == 1]


    def get_chromosomes(self):
        if self.snps_filter is None:
            return self.h5file['chromosomes'][...]
        else:
            return self.h5file['chromosomes'][self.snps_filter == 1]


    def get_chr_pos_list(self):
        return zip(self.get_chromosomes(), self.get_positions())

    def get_macs(self):
        g, already_exists = self._get_cached_group_()
        if not already_exists:
            self._update_macs_(g)
        if self.snps_filter is None:
            return g['macs'][...]
        else:
            return g['macs'][self.snps_filter == 1]

    def get_indivs(self):
        if self.indiv_filter is None:
            return self.h5file['indiv_ids'][...]
        else:
            return self.h5file['indiv_ids'][self.indiv_filter]


    def get_mafs(self):
        return self.get_macs() / self.num_individs()

    def get_old_snps_data_set(self):
        """
        """
        raise NotImplementedError

    def close(self):
        self.h5file.close()


    def snps_chunks(self, chunk_size=1000):
        """
        An generator/iterator for SNP chunks.
        """
        n_snps = int(self.h5file['num_snps'][...])
        if self.snps_filter is None:
            if self.indiv_filter is None:
                for i in range(0, n_snps, chunk_size):
                    stop_i = min(i + chunk_size, n_snps)
                    yield self.h5file['snps'][i:stop_i]
            else:
                for i in range(0, n_snps, chunk_size):
                    stop_i = min(i + chunk_size, n_snps)
                    yield self.h5file['snps'][i:stop_i, self.indiv_filter]
        else:
            if self.indiv_filter is None:
                for i in range(0, n_snps, chunk_size):
                    stop_i = min(i + chunk_size, n_snps)
                    filter_chunk = self.snps_filter[i:stop_i]
                    snps_chunk = self.h5file['snps'][i:stop_i]
                    yield snps_chunk[filter_chunk == 1]
            else:
                for i in range(0, n_snps, chunk_size):
                    stop_i = min(i + chunk_size, n_snps)
                    filter_chunk = self.snps_filter[i:stop_i]
                    snps_chunk = self.h5file['snps'][i:stop_i, self.indiv_filter]
                    yield snps_chunk[filter_chunk == 1]


    def snps(self, chunk_size=1000):
        """
        An generator/iterator for the SNPs.
        """
        for snps_chunk in self.snps_chunks(chunk_size):
            for snp in snps_chunk:
                yield snp


    def _calc_ibs_kinship_(self, dtype='single', chunk_size=None):
        n_snps = self.num_snps()
        n_indivs = self.num_individs()
        if chunk_size is None:
            chunk_size = n_indivs
        #print 'Allocating K matrix'
        k_mat = sp.zeros((n_indivs, n_indivs), dtype=dtype)
        #print 'Starting calculation'
        i = 0
        snps_chunks = self.snps_chunks(chunk_size)
        for snps_chunk in snps_chunks: #FINISH!!!
            i += len(snps_chunk)
            snps_array = snps_chunk.T
            if self.data_format == 'diploid_int':
                for i in range(n_indivs):
                    for j in range(i):
                        bin_counts = sp.bincount(sp.absolute(snps_array[j] - snps_array[i]))
                        if len(bin_counts) > 1:
                            k_mat[i, j] += (bin_counts[0] + 0.5 * bin_counts[1])
                        else:
                            k_mat[i, j] += bin_counts[0]
                        k_mat[j, i] = k_mat[i, j]
            elif self.data_format == 'binary':
                sm = sp.mat(snps_array * 2.0 - 1.0)
                k_mat = k_mat + sm * sm.T
            sys.stdout.write('\b\b\b\b\b\b%0.1f%%' % (100.0 * i / n_snps))
            sys.stdout.flush()
        if self.data_format == 'diploid_int':
            k_mat = k_mat / float(n_snps) + sp.eye(n_indivs)
        elif self.data_format == 'binary':
            k_mat = k_mat / (2 * float(n_snps)) + 0.5
        return k_mat



    def _calc_ibd_kinship_(self, dtype='single', chunk_size=None):
        n_snps = self.num_snps()
        n_indivs = self.num_individs()
        if chunk_size is None:
            chunk_size = n_indivs
        k_mat = sp.zeros((n_indivs, n_indivs), dtype=dtype)
        snps_chunks = self.snps_chunks(chunk_size)
        i = 0
        for snps_chunk in snps_chunks:
            i += len(snps_chunk)
            snps_array = snps_chunk.T
            norm_snps_array = (snps_array - sp.mean(snps_array, 0)) / sp.std(snps_array, 0)
            x = sp.mat(norm_snps_array.T)
            k_mat += x.T * x
            sys.stdout.write('\b\b\b\b\b\b%0.1f%%' % (100.0 * i / n_snps))
            sys.stdout.flush()
        k_mat = k_mat / float(n_snps)
        return k_mat



    def get_kinship(self, method='ibs', dtype='single', chunk_size=1024):
        """
        Returns kinship
        """
        print method
        if method == 'ibd':
            print 'Starting IBD calculation'
            k_mat = self._calc_ibd_kinship_(dtype=dtype, chunk_size=chunk_size)
            print '\nFinished calculating IBD kinship matrix'
            return k_mat
        elif method == 'ibs':
            print 'Starting IBS calculation'
            k_mat = self._calc_ibs_kinship_(dtype=dtype, chunk_size=chunk_size)
            print '\nFinished calculating IBS kinship matrix'
            return k_mat




    def get_region_split_snps(self, chrom, start_pos, stop_pos):
        """
        Returns two SNP sets, one with the SNPs within the given region, 
        and the other with the remaining SNPs.
        """
        chr_pos_l = self.get_chr_pos_list()
        start_i = bisect.bisect(chr_pos_l, (chrom, start_pos))
        stop_i = bisect.bisect(chr_pos_l, (chrom, stop_pos))
        if self.cached_snps is None:
            self.cached_snps = self.get_snps()
        snps = self.cached_snps
        local_snps = snps[start_i:stop_i]
        global_snps = snps[:start_i] + snps[stop_i:]
        return local_snps, global_snps


#    def get_local_n_global_kinships(self, focal_chrom_pos=None, window_size=25000, chrom=None, start_pos=None,
#                    stop_pos=None, kinship_method='ibd', global_kinship=None, verbose=False):
#        """
#        Returns local and global kinship matrices.
#        """
#        if focal_chrom_pos != None:
#            chrom, pos = focal_chrom_pos
#            start_pos = pos - window_size
#            stop_pos = pos + window_size
#
#        local_snps, global_snps = self.get_region_split_snps(chrom, start_pos, stop_pos)
#        if verbose:
#            print 'Found %d local SNPs' % len(local_snps)
#            print 'and %d global SNPs' % len(global_snps)
#        if kinship_method == 'ibd':
#            local_k = self._calc_ibd_kinship_2_(local_snps, num_dots=0) if len(local_snps) else None
#            if global_kinship is None:
#                global_k = self._calc_ibd_kinship_2_(global_snps, num_dots=0) if len(global_snps) else None
#        elif kinship_method == 'ibs':
#            local_k = self._calc_ibs_kinship_2_(local_snps, num_dots=0) if len(local_snps) else None
#            if global_kinship is None:
#                global_k = self._calc_ibs_kinship_2_(global_snps, num_dots=0) if len(global_snps) else None
#        else:
#            raise NotImplementedError
#        if global_kinship != None:
#            global_k = (global_kinship * self.num_snps() - local_k * len(local_snps)) / len(global_snps)
#        return {'local_k':local_k, 'global_k':global_k, 'num_local_snps':len(local_snps),
#            'num_global_snps':len(global_snps)}



class mixed_model():
    """
    A wrapper class for conducting mixed models GWAS using h5py datasets
    """

    def __init__(self, y, fixed_effects=None, K=None, Z=None):
        self.lmm = linear_models.LinearMixedModel(Y=y)
        if Z is not None:
            self.lmm.add_random_effect(Z * K * Z.T)
            if fixed_effects is not None:
                for cofactor in fixed_effects:
                    self.lmm.add_factor(Z * cofactor)
        else:
            self.lmm.add_random_effect(K)
            if fixed_effects:
                for cofactor in fixed_effects:
                    self.lmm.add_factor(cofactor)

    def _emmax_(self, gd, Z=None, with_betas=False, inference_method='REML', emma_num=100):
        """
        EMMAX implementation (in python)
        Single SNPs
        """
        print 'Calculating the eigenvalues of K'
        s0 = time.time()
        eig_L = self.lmm._get_eigen_L_()
        print 'Done.'
        print 'Took %0.2f seconds' % (time.time() - s0)

        print "Calculating the eigenvalues of S(K+I)S where S = I-X(X'X)^-1X'"
        s0 = time.time()
        eig_R = self.lmm._get_eigen_R_(X=self.lmm.X)
        print 'Done'
        print 'Took %0.2f seconds' % (time.time() - s0)

        print 'Getting variance estimates'
        s0 = time.time()
        res = self.lmm.get_estimates(eig_L, method=inference_method, eig_R=eig_R) #Get the variance estimates..
        print 'Done.'
        print 'Took %0.2f seconds' % (time.time() - s0)
        print 'pseudo_heritability:', res['pseudo_heritability']
        H_sqrt_inv = res['H_sqrt_inv']

        s0 = time.time()

        dtype = 'single'
        q = 1  # Single SNP is being tested
        p = len(self.lmm.X.T) + q
        n = self.lmm.n
        n_p = n - p
        num_snps = gd.num_snps()

        h0_X = sp.mat(H_sqrt_inv * self.lmm.X, dtype=dtype)
        Y = H_sqrt_inv * self.lmm.Y    #The transformed outputs.
        (h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X, Y)
        Y = sp.mat(Y - h0_X * h0_betas, dtype=dtype)
        h0_betas = map(float, list(h0_betas))

        if Z is not None:
            H_sqrt_inv = H_sqrt_inv * Z

        if not with_betas:
            (Q, R) = linear_models.qr_decomp(h0_X)  #Do the QR-decomposition for the Gram-Schmidt process.
            Q = sp.mat(Q)
            Q2 = Q * Q.T
            M = sp.mat(H_sqrt_inv.T * (sp.eye(n) - Q2), dtype=dtype)
        else:
            betas_list = [h0_betas] * num_snps
            M = H_sqrt_inv.T

        rss_list = sp.repeat(h0_rss, num_snps)
        chunk_size = len(Y)
        i = 0
        snps_chunks = gd.snps_chunks(chunk_size)
        for snps_chunk in snps_chunks: #Do the dot-product in chuncks!
            if len(snps_chunk) == 0:
                continue

            Xs = sp.matrix(snps_chunk) * M
            for X_j in Xs:
                if with_betas:
                    (betas, rss, p, sigma) = linalg.lstsq(sp.hstack([h0_X, X_j.T]), Y, \
                                    overwrite_a=True)
                    if rss:
                        betas_list[i] = map(float, list(betas))
                else:
                    (betas, rss, p, sigma) = linalg.lstsq(X_j.T, Y, overwrite_a=True)
                if rss:
                    rss_list[i] = rss[0]

                sys.stdout.write('\b\b\b\b\b\b%0.1f%%' % (100.0 * i / num_snps))
                sys.stdout.flush()
                i += 1

        sys.stdout.write('\b\b\b\b\b\b\b100.0%\n')
        sys.stdout.flush()

        rss_ratio = h0_rss / rss_list
        var_perc = 1 - 1 / rss_ratio
        #assert sp.all(var_perc < 1.01), '%f\n%s\n%s' % (h0_rss, str(var_perc[var_perc < 1.01]), str(rss_list[var_perc < 1.01]))
        f_stats = (rss_ratio - 1) * n_p / float(q)
        p_vals = stats.f.sf(f_stats, q, n_p)

        res_d = {'ps':p_vals, 'f_stats':f_stats, 'rss':rss_list, 'var_perc':var_perc,
            'h0_rss':h0_rss, 'h0_betas':h0_betas}
        if with_betas:
            res_d['betas'] = betas_list

        return res_d


    def emmax(self, gd, test_method='f_test', Z=None, with_betas=False):
        """
        gd is a h5py genotype_data (coordinated)
        """
        assert test_method == 'f_test', 'Only F-test is implemented at the moment'
        print "Running EMMAX"
        s1 = time.time()
        res = self._emmax_(gd, Z=Z, with_betas=with_betas)
        secs = time.time() - s1
        if secs > 60:
            mins = int(secs) / 60
            secs = secs - mins * 60
            print 'Took %d mins and %f seconds.' % (mins, secs)
        else:
            print 'Took %f seconds.' % (secs)
        return res




def test_genotype_data():
    import dataParsers as dp
    sd = dp.load_snps_call_method(75)
    gd = genotype_data('/tmp/test.h5py')



