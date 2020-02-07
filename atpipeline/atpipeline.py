#!/usr/local/bin/python2.7
# encoding: utf-8
'''
atpipeline -- shortdesc

atpipeline is a wrapper script for unning GWAS

It defines classes_and_methods

@author:     Ümit Seren
        
@copyright:  2012 Gregor Mendel Institute. All rights reserved.
        
@license:    MIT

@contact:    uemit.seren@gmail.com
@deffield    updated: Updated
'''

import sys
import os

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from progress import ProgressMessenger, StdoutMessenger
#from atgwas import phenotypeData
#from atgwas import dataParsers
#from atgwas import kinship
#from atgwas import linear_models as lm
#from atgwas import  util
#from atgwas import mtcorr
#from atgwas import analyze_gwas_results as agr

import phenotypeData
import dataParsers
import kinship
import linear_models as lm
import  util
import mtcorr
import analyze_gwas_results as agr

import scipy
import itertools
import math
import h5py
import numpy
import pdb

__all__ = []
__version__ = 0.1
__date__ = '2012-12-06'
__updated__ = '2012-12-06'




def main(argv=None): 
    '''Command line options.'''
    
    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''%s

  Created by Ümit Seren on %s.
  Copyright 2012 Gregor Mendel Institute. All rights reserved.
  
  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0
  
  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
''' % (program_shortdesc, str(__date__))

    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument("-l", "--list_genotypes", dest="list_genotypes", help="display available genotype dataset", action='store_true')
        parser.add_argument("-t", "--transformation", dest="transformation", help="Apply a transformation to the data. Default[None]", choices=["log", "sqrt", "exp", "sqr", "arcsin_sqrt", "box_cox"])
        parser.add_argument("-a", "--analysis_method", dest="analysis_method", help="analyis method to use",required=True,choices=["lm", "emma", "emmax", "kw", "ft", "emmax_anova", "lm_anova", "emmax_step", "lm_step","loc_glob_mm","amm"])
        parser.add_argument("-g", "--genotype", dest="genotype", help="genotype dataset to be used in the GWAS analysis (run with option -l to display list of available genotype datasets)", required=True, type=int,metavar="INTEGER" )
        parser.add_argument("-k", "--kinship", dest="kinship", help="Specify the file containing the kinship matrix. (otherwise default file is used or it's generated.)", metavar="FILE" )
        parser.add_argument("-s", "--kinship_type", dest="kinship_type", help="Type of kinship calculated. Possible types are ibs (default) or ibd ", choices=["ibs", "ibd"],default="ibs")
        parser.add_argument("-q", "--queue", dest="queue", help="Send status updates to Message Broker", action='store_true')
        parser.add_argument("-z", "--queue_host", dest="queue_host", help="Host of the Message Broker")
        parser.add_argument("-o", "--output_file", dest="outputfile", help="Name of the output file")
        parser.add_argument('-V', '--version', action='version', version=program_version_message)
        parser.add_argument(dest="file", help="csv file containing phenotype values", metavar="FILE")
        
        # Process arguments
        args = parser.parse_args()
        messenger = StdoutMessenger()
        if args.queue: 
            messenger = ProgressMessenger(args.queue_host,5672,'admin','eastern')
        
        messenger.update_status(progress=0.0, task_status='Loading phenotype data')
        phenData = phenotypeData.parse_phenotype_file(args.file,False)  #load phenotype file
        phen_ids = phenData.phen_dict.keys()  # get phenotype ids
        #If not analysis plots... then GWAS
        for phen_id in phen_ids:
            phenotype_name = phenData.get_name(phen_id)
            messenger.update_status(progress=0.0, task_status='Loading phenotype data')
            print "Performing GWAS for phenotype: %s, phenotype_id: %s" % (phenotype_name, phen_id)
            _perform_gwas_(phen_id, phenData, args.analysis_method, args.transformation,args.genotype,args.kinship_type,args.kinship,messenger,args.outputfile)
        return 0
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception, e:
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2
    

def _perform_gwas_(phen_id,phenData,analysis_method,transformation,genotype,kinship_type,kinshipFile=None,messenger=None,outputfile=None):
    additional_columns = {}
    messenger.update_status(progress=0.0, task_status='Loading genotype data')
    genotypeData = dataParsers.load_snps_call_method(genotype)
    #genotypeData = dataParsers.load_hdf5_snps_call_method(genotype)
    K = None
    messenger.update_status(step=0.05, task_status='Preparing data')
    n_filtered_snps = _prepare_data_(genotypeData,phenData,phen_id)
    phen_vals = phenData.get_values(phen_id)
    if analysis_method in ['emma', 'emmax', 'emmax_anova', 'emmax_step', 'loc_glob_mm','amm']:
        #Load genotype file (in binary format)
        sys.stdout.write("Retrieving the Kinship matrix K.\n")
        sys.stdout.flush()
        if kinshipFile:   #Kinship file was supplied..
            messenger.update_status(progress=0.15, task_status='Loading supplied kinship file: %s' % kinshipFile)
            print 'Loading supplied kinship file: %s' % kinshipFile
            K = kinship.load_kinship_from_file(kinshipFile, genotypeData.accessions)
        else:
            messenger.update_status(progress=0.15, task_status='Loading kinship file')
            print 'Loading kinship file.'
            K = kinship.get_kinship(call_method_id=genotype,
                                            method=kinship_type, n_removed_snps=n_filtered_snps,
                                            remain_accessions=genotypeData.accessions)
            sys.stdout.flush()
            sys.stdout.write("Done!\n")

    snps = genotypeData.getSnps()
    positions = genotypeData.getPositions()
    chromosomes = []
    for i, (s, c) in enumerate(itertools.izip(genotypeData.snpsDataList, genotypeData.chromosomes)):
        chromosomes.extend([c] * len(s.snps))
        maf_dict = genotypeData.get_mafs()
    
    if analysis_method in ['kw']:
        messenger.update_status(progress=0.7, task_status='Performing KW')
        res = util.kruskal_wallis(snps, phen_vals)
        
    elif analysis_method in ['loc_glob_mm']:
        raise NotImplementedError
    elif analysis_method in ['emma']:
        res = lm.emma(snps, phen_vals, K)
    elif analysis_method in ['emmax','amm']:
        d = lm.emmax_step(phen_vals, genotypeData, K, [], emma_num=100)
        res = d['res']
        #additional_columns['stats'] = d['stats']
    elif analysis_method in ['lm']:
        d = lm.lin_reg_step(phen_vals, genotypeData, [])
        res = d['res']
        #additional_columns['stats'] = d['stats']
    else:
        raise Exception('analysis method %s not supported' % analysis_method)
    
    pvals = res['ps']
    
    #Calculate Benjamini-Hochberg threshold
    bh_thres_d = mtcorr.get_bhy_thres(res['ps'], fdr_thres=0.05)
    #Calculate Median p-value
    med_pval = agr.calc_median(res['ps'])
    #Calculate the Kolmogorov-Smirnov statistic
    ks_res = agr.calc_ks_stats(res['ps'])
    
    quantiles_dict = _calculate_qqplot_data_(pvals)
    scores = map(lambda x:-math.log10(x), pvals)
    
    if analysis_method in ['lm', 'emma', 'emmax','amm']:
        additional_columns['genotype_var_perc'] = res['var_perc']
        if 'betas' in res:
            betas = map(list, zip(*res['betas']))
            additional_columns['beta0'] = betas[0]
            if len(betas) > 1:
                additional_columns['beta1'] = betas[1]
    
    #calculate ld
    if outputfile is None:
         outputfile = "%s.hdf5" % phen_id
    messenger.update_status(progress=0.8, task_status='Processing and saving results')
    _save_hdf5_pval_file(outputfile, analysis_method, transformation,chromosomes, positions, scores, maf_dict['marfs'], maf_dict['mafs'], 
                         quantiles_dict,ks_res,bh_thres_d['thes_pval'],med_pval,additional_columns)
    
def _prepare_data_(genotypeData, phenData, phen_id,with_replicates=False):
    """
    Coordinates phenotype and snps data for different mapping methods.
    """
    if not with_replicates:
        print 'Converting replicates of phenotypes to averages'
        phenData.convert_to_averages([phen_id])
    d = genotypeData.coordinate_w_phenotype_data(phenData, phen_id)
    return d['n_filtered_snps']

def _calculate_qqplot_data_(pvals,num_dots=1000,max_val=6):
    quantiles = agr.get_quantiles(pvals, num_dots=num_dots)
    exp_quantiles = agr._getExpectedPvalueQuantiles_(num_dots)
    log_quantiles = agr.get_log_quantiles(pvals, num_dots=num_dots, max_val=max_val)
    exp_log_quantiles = scipy.arange(1, num_dots + 1, dtype='single') / (num_dots + 1) * max_val

    quantiles_dict = {'quantiles':quantiles, 'exp_quantiles':exp_quantiles,
            'log_quantiles':log_quantiles, 'exp_log_quantiles':exp_log_quantiles}
    return quantiles_dict

def _save_hdf5_pval_file(outputfile,analysis_method,transformation,chromosomes,positions,scores,mafs,macs,
                         quantiles_dict,ks_stats,bh_thresh,med_pval,additional_columns):
    f = h5py.File(outputfile,'w') 
    
    # store quantiles
    quant_group = f.create_group('quantiles')
    quantiles_array = zip(quantiles_dict['exp_quantiles'],quantiles_dict['quantiles'])
    log_quantiles_array = zip(quantiles_dict['exp_log_quantiles'],quantiles_dict['log_quantiles'])
    quant_group.create_dataset('quantiles',(len(quantiles_dict['quantiles']), 2),'f8',data=quantiles_array)
    quant_group.create_dataset('log_quantiles',(len(quantiles_dict['log_quantiles']), 2),'f8',data=log_quantiles_array)
    
    #store pvalues
    pvals_group = f.create_group('pvalues')
    pvals_group.attrs['numberOfSNPs'] = len(scores)
    pvals_group.attrs['max_score'] = max(scores)
    pvals_group.attrs['analysis_method'] = analysis_method
    if transformation is None:
        transformation = "raw"
    pvals_group.attrs['transformation'] = transformation
    pvals_group.attrs['bonferroni_threshold'] = -math.log10(0.05 / len(scores))
    pvals_group.attrs['ks_stat'] = ks_stats['D']
    pvals_group.attrs['ks_pval'] = ks_stats['p_val']
    pvals_group.attrs['med_pval'] = med_pval
    pvals_group.attrs['bh_thres'] =-math.log10(bh_thresh)
    data = numpy.array(zip(chromosomes, positions, scores, mafs, macs,*additional_columns.values()))
    
    for chr in range(1,6):
        chr_group = pvals_group.create_group('chr%s' % chr)
        chr_data = data[numpy.where(data[:,0] == chr)]
        chr_data =chr_data[chr_data[:,2].argsort()[::-1]]
        positions = chr_data[:,1]
        chr_group.create_dataset('positions',(len(positions),),'i4',data=positions)
        scores = chr_data[:,2]
        chr_group.create_dataset('scores',(len(scores),),'f8',data=scores)
        mafs = chr_data[:,3]
        chr_group.create_dataset('mafs',(len(mafs),),'f8',data=mafs)
        macs = chr_data[:,4]
        chr_group.create_dataset('macs',(len(macs),),'i4',data=macs)

        if chr_data.shape[1] > 5: 
            for i,key in enumerate(additional_columns.keys()):
                values = chr_data[:,5+i]
                chr_group.create_dataset(key,values.shape,values.dtype,data=values) 
    f.close()

if __name__ == "__main__":
    sys.exit(main())
