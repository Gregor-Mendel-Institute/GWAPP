'''
Created on Nov 6, 2012
@author: dazhe.meng

Script to transform a phenotype/snps according to emmax assumption
effectively, on the transformed set, simple least squares can be used to estimate betas

borrows and uses Bjarni's scripts heavily
'''

import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from itertools import izip

import math
import numpy as np
import scipy as sp
from scipy import linalg
from scipy import stats

import phenotypeData as pd
import kinship as kin
import linear_models as lm
import dataParsers as dp

def sys_write(message):
    sys.stdout.write(message)
    sys.stdout.flush()

parser = ArgumentParser(description="Script to transform phenotypes/snps according to emmax assumptipn", epilog="", formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("-s","--snpfile", help="The SNP file to use (support for call method coming)", default="/projects/1001genomes/dazhe/SNPdata/250K_m75.csv")
parser.add_argument("-p","--phenfile", help="Phenotype file to load", default="/projects/1001genomes/dazhe/phen_data/phen_raw_071910.tsv")
parser.add_argument("-i","--phenid", help="Specifies which phenotype to actually run", type=int, default=1)
parser.add_argument("-o","--outprefix", help="File prefix for output files, should contain directory as well", default="test/test")
parser.add_argument("--minmaf", help="Minimum minor allele frequency to filter out", type=float, default=0.05)
args = parser.parse_args()

class c_transform_data():
    ' using a class to hold all globals '
    def __init__(self):
        self._load_data_()
    
    def run(self):
        self._prepare_data_()
        self._get_estimates_()
        self._transform_()
        self.output()
    
    def _load_data_(self, delimiter=','):
        print "Loading phenotype and genotype data..."
        self.pd = pd.parse_phenotype_file(args.phenfile, delim=delimiter)
        (snpsds, chromosomes) = dp.parse_raw_snps_data(args.snpfile, target_format='binary', missing_val='?', debug_filter=0.01, return_chromosomes=True) # retain only 1% data for now 
        self.sd = dp.SNPsDataSet(snpsds, chromosomes, data_format='binary')
        
    def _prepare_data_(self):
        print 'Converting replicates of phenotypes to averages...'
        self.pd.convert_to_averages([args.phenid])
        
        print "Coordinating phenotype and genotype data..."
        d = self.sd.coordinate_w_phenotype_data(self.pd, args.phenid)
        
        if args.minmaf > 0:
            mac_threshold = int(math.ceil(len(self.sd.accessions) * (args.minmaf)))
            print "Applying mac filtering of %s..."%mac_threshold
            self.sd.filter_mac_snps(mac_threshold)
    
        print "Calculating kinship..."
        self.k = self.sd.get_ibs_kinship_matrix()
        
        print "Extracting SNP matrix X and phenotype vector Y..."
        self.snps = self.sd.getSnps()
        self.pvls = self.pd.get_values(args.phenid)
    
    def _get_estimates_(self):
        print "Initializing mixed model..."
        self.lmm = lm.LinearMixedModel(self.pvls)
        self.lmm.add_random_effect(self.k)
        eig_L = self.lmm._get_eigen_L_()
        print "Estimating variance components..."
        self.est = self.lmm.get_estimates(eig_L, self.k)
    
    def _transform_(self, dtype='single'):
        print "Tranforming data..."
        H_sqrt_inv = self.est['H_sqrt_inv']
        H_sqrt_inv_t = H_sqrt_inv.T
        self.Yt = H_sqrt_inv * self.lmm.Y # transformed phenotypes
        self.h0_X = H_sqrt_inv * self.lmm.X # transformed intercept
        chunk_size = len(self.pvls)
        num_snps = len(self.snps)
        self.t_snps = []
        for i in range(0, num_snps, chunk_size): #Do the dot-product in chuncks!
            snps_chunk = sp.matrix(self.snps[i:i + chunk_size], dtype=dtype)
            Xs = snps_chunk * H_sqrt_inv_t
            for X_j in Xs:
                self.t_snps.append(sp.array(X_j).flatten())
        
    def output(self):
        print "Outputting data..."
        f_y = open(args.outprefix+'.y', 'w')
        f_y.write("\n".join([str(yt) for yt in self.Yt.reshape(-1,).tolist()[0]]))
        f_i = open(args.outprefix+'.i','w')
        f_i.write("\n".join([str(xt) for xt in self.h0_X.reshape(-1,).tolist()[0]]))
        f_s = open(args.outprefix+'.bsnps', 'w')
        f_ts = open(args.outprefix+'.tsnps', 'w')
        f_s.write("Chr,Pos,"+",".join(self.sd.snpsDataList[0].accessions)+'\n')
        f_ts.write("Chr,Pos,"+",".join(self.sd.snpsDataList[0].accessions)+'\n')
        l_pos = []
        for chromosome, snpsd in izip(self.sd.chromosomes, self.sd.snpsDataList):
            for pos in snpsd.positions:
                l_pos.append((chromosome,pos))
        if len(l_pos)!=len(self.snps) or len(self.snps)!=len(self.t_snps):
            raise ValueError("Number of SNPs does not match! pos: %s ; snps: %s ; tsnps: %s"%(len(l_pos), len(self.snps), len(self.t_snps)))
        for snpi in xrange(len(l_pos)):
            f_s.write("%s,%s,"%l_pos[snpi]+",".join(map(str,list(self.snps[snpi])))+'\n')
            f_ts.write("%s,%s,"%l_pos[snpi]+",".join(map(str,list(self.t_snps[snpi])))+'\n')
        self._output_hdf5_()
    
    def _output_hdf5_(self):
        ' does nothing now but should add later! '
        pass
    
TD = c_transform_data()
TD.run()
