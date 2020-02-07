'''
Created on Oct 13, 2011

@author: uemit.seren
'''

import tables
import numpy
import numpy.lib.recfunctions
import csv
import time
import re
import cPickle
import pdb
import os
import itertools
import scipy as sp
from optparse import OptionParser
from datetime import datetime
import math

chromosome_ends = [30429953, 19701870, 23467451, 18578708, 26992130]
cvt_types = tables.Enum(['radius', 'tss_upstream'])
result_types = tables.Enum(['EX', 'LM', 'KW', 'full', 'genetic', 'environ'])
environ_types = tables.Enum(['10', '16', 'GxE'])

class Accessions(tables.IsDescription):

    accession_id = tables.Int32Col()
    name = tables.StringCol(16)
    latitude = tables.Float32Col()
    longitude = tables.Float32Col()
    country = tables.StringCol(256)
    country_ISO = tables.StringCol(3)
    collector = tables.StringCol(256)
    collection_date = tables.Time32Col()
    dataset = tables.EnumCol(['10', '16', 'both'], '10', base='uint8')


class PhenotypeValue(tables.IsDescription):
    """
    Phenotype value class
    """
    ecotype = tables.Int32Col()
    value = tables.Float32Col()

class ResultRecord(tables.IsDescription):
    """
    Kruskal Wallis
    """
    chromosome = tables.Int32Col()
    position = tables.Int32Col()
    score = tables.Float64Col()
    maf = tables.Float32Col()
    mac = tables.Int32Col()



class ResultRecordLM(ResultRecord):
    """
    Linear model, mixed models, etc.
    """
    genotype_var_perc = tables.Float32Col()
    beta0 = tables.Float32Col()
    beta1 = tables.Float32Col()
    correlation = tables.Float32Col()

def loadPhenotypeFile(phenotype_file):
    reader = csv.reader(open(phenotype_file, 'rb'), delimiter=',')
    header = reader.next()
    last_pid = -1
    phenotypes = []
    d = {}
    for row in reader:
        pid = int(row[0])
        if last_pid != pid:
            if last_pid != -1:
                d['values'].sort()
                phenotypes.append(d)
            phen_name = row[1]
            d = {'name':phen_name, 'values':[]}
        if row[2] == 'ACap':
            ecotype = -1
        elif row[2] == 'ALyr':
            ecotype = -2
        else:
            ecotype = int(row[2])
        d['values'].append((ecotype, float(row[3])))
        last_pid = pid
    d['values'].sort()
    phenotypes.append(d)
    return phenotypes



class RNASeqRecords:

    def __init__(self, hdf5_file_name, mode='r'):
        self.h5_f = tables.open_file(hdf5_file_name, mode)

    def close(self):
        self.h5_f.close()

    def getAccessions(self):
        table = self.h5_f.root.accessions.infos
        accessions = []
        datasets = table.get_enum('dataset')
        accessions_ix = {}
        i = 0
        for row in table:
            collection_date = datetime.fromtimestamp(int(row['collection_date']))
            collection_date = datetime.strftime(collection_date, '%d.%m.%Y')
            accession = {'accession_id':row['accession_id'], 'collection_date':collection_date, 'collector':unicode(row['collector'], 'latin1'), 'country':row['country_ISO'], 'dataset':datasets.__call__(row['dataset']), 'latitude':row['latitude'], 'longitude':row['longitude'], 'name':unicode(row['name'], 'utf8')}
            accessions.append(accession)
            accessions_ix[row['accession_id']] = i
            i = i + 1
        return accessions, accessions_ix

    def getTopResults(self, environ_type , result_type , range_start, range_length, gene='', snp_chr='', snp_pos='', gene_chr='', min_score=''):
        keys_to_return = ['gene', 'snp_pos', 'snp_chr', 'gene_mid_pos', 'gene_chr', 'mac', 'maf', 'perc_var_expl', 'score', 'gene_start', 'gene_end']
        if result_type not in ['EX', 'LM', 'KW']:
            result_type = result_type.lower()
        results_to_filter = []
        results_to_return = []
        stop = int(range_start) + int(range_length)
        try:
            table = self.h5_f.get_node('/top_results/%s/%s' % (environ_type, result_type))
        except Exception, err:
            return [], 0, 0
        isCriteria = False
        if (gene != '')  \
            or (snp_pos != '') or (snp_chr != '') \
            or (gene_chr != '') or (min_score != '') :
            isCriteria = True
        if isCriteria == False:
            results_to_filter = table
        else:
            op = ' and '
            filter = ''
            if snp_chr != '':
                filter = filter + '%s(snp_chr == %s)' % (op if filter != '' else '', snp_chr)
            if snp_pos != '':
                filter = filter + '%s(snp_pos == %s)' % (op if filter != '' else '', snp_pos)
            if gene_chr != '':
                filter = filter + '%s(gene_chr == %s)' % (op if filter != '' else '', gene_chr)
            if min_score != '':
                filter = filter + '%s(score >= %s)' % (op if filter != '' else '', float(min_score))
            if gene == '':
                results_to_filter = table.read_where(filter)
            else:
                ix_to_fetch = []
                results_to_filter = []
                if filter != '':
                    iterator = table.where(filter)
                else:
                    iterator = table
                for result in iterator:
                    if gene.lower() in unicode(result['gene'], 'utf8').lower():
                        ix_to_fetch.append(result.nrow)
                results_to_filter = table.read_coordinates(ix_to_fetch)
        count = len(results_to_filter)
        if stop == -1:
            stop = count
        if stop > count:
            range_start = range_start - stop
            stop = count
        if range_start < 0:
            range_start = 0
        for i in range(range_start, stop):
            result = results_to_filter[i]
            result_to_return = dict((k, result[k].tolist()) for k in keys_to_return)
            results_to_return.append(result_to_return)
        return results_to_return, count, range_start

    def getPhenotypes(self, range_start, range_length, name='', chr='', start='', end=''):
        phenotypes_to_filter = []
        phenotypes_to_return = []
        stop = int(range_start) + int(range_length)
        isCriteria = False
        if (name != '') or (chr != '') or (start != '') or (end != ''):
            isCriteria = True
        if isCriteria == False:
            phenotypes_to_filter = self.h5_f.root.phenotype_list
        else:
            op = ' and '
            filter = ''
            if chr != '':
                filter = filter + '%s(chr == %s)' % (op if filter != '' else '', chr)
            if start != '':
                filter = filter + '%s(start >= %s)' % (op if filter != '' else '', start)
            if end != '':
                filter = filter + '%s(end <= %s)' % (op if filter != '' else '', end)
            if name == '':
                phenotypes_to_filter = self.h5_f.root.phenotype_list.read_where(filter)
            else:
                phenotypes_to_filter = []
                if filter != '':
                    iterator = self.h5_f.root.phenotype_list.where(filter)
                else:
                    iterator = self.h5_f.root.phenotype_list
                for phenotype_to_filter in iterator:
                    if name.lower() in unicode(phenotype_to_filter['name'], 'utf8').lower():
                        phenotype_to_add_filter = {'name':phenotype_to_filter['name'], 'chr':phenotype_to_filter['chr'],
                                                   'start':phenotype_to_filter['start'], 'end':phenotype_to_filter['end'],
                                                   'id':phenotype_to_filter['id'],
                                                   'min_score_10C':phenotype_to_filter['min_score_10C'], 'min_score_16C':phenotype_to_filter['min_score_16C'],
                                                   'min_score_full':phenotype_to_filter['min_score_full'],
                                                   'pseudo_heritability_10C':phenotype_to_filter['pseudo_heritability_10C'], 'pseudo_heritability_16C':phenotype_to_filter['pseudo_heritability_16C']}
                        phenotypes_to_filter.append(phenotype_to_add_filter)

        count = len(phenotypes_to_filter)
        if stop > count:
            range_start = range_start - stop
            stop = count
        if range_start < 0:
            range_start = 0
        for i in range(range_start, stop):
            phenotype = phenotypes_to_filter[i]
            phenotypes_to_return.append({'name':phenotype['name'], 'chr':int(phenotype['chr']), 'start':int(phenotype['start'])
                                         , 'end':int(phenotype['end']), 'phenotype_id':int(phenotype['id'])
                                         , 'maxScore10C':round(float(phenotype['min_score_10C']), 2), 'maxScore16C':round(float(phenotype['min_score_16C']), 2)
                                         , 'maxScoreFull':round(float(phenotype['min_score_full']), 2)
                                         , 'pseudoHeritability10C':round(float(phenotype['pseudo_heritability_10C']), 4), 'pseudoHeritability16C':round(float(phenotype['pseudo_heritability_16C']), 4)})
        return phenotypes_to_return, count, range_start

    def get_phenotype_values(self, phen_name, environment, dataset, transformation='raw'):
        import bisect
        """
        Returns the phenotype values
        """
        #self.h5file = self._open(mode="r")

        table = self.h5_f.get_node('/phenotypes/%s/%s/Fullset/%s/values' % (phen_name, environment, transformation))
        d = {'ecotype' : [], 'value' : []}
        if dataset == "Fullset":
            for x in table.iterrows():
                for k in d:
                    d[k].append(x[k])
        else:
            group = self.h5_f.get_node("/phenotypes/%s/%s/%s" % (phen_name, environment, dataset))
            ids = table.read(field='ecotype')
            indices = [bisect.bisect(ids, val) - 1 for val in group._v_attrs.ecotypes]
            for x in table.itersequence(indices):
                for k in d:
                    d[k].append(x[k])
        #self._close()
        return d

    def get_phenotype_bins(self, phen_name, environment, dataset='Fullset', transformation='raw', bin_number=20):
        if environment == 'both':
            phen_vals = self.get_phenotype_values(phen_name, "T10C", dataset, transformation)['value']
            phen_vals2 = self.get_phenotype_values(phen_name, "T16C", dataset, transformation)['value']
        else:
            phen_vals = self.get_phenotype_bins(phen_name, environment, transformation, dataset, bin_number)
        return self._get_phenotype_bins(phen_vals, bin_number, phen_vals2)



    def get_phenotype_info(self, id):
        phenotype_group = self.h5_f.get_node('/phenotypes', id)
        phenotype_info = {'id':0, 'name': id, 'environments':self._get_phenotype_environments_(id), 'gxeResults':self._get_gxe_results(id)}
        return phenotype_info

    def _get_gxe_results(self, phenotype):
        phenotype_group = self.h5_f.get_node('/phenotypes/', phenotype)
        gxeResults = []
        g = self.h5_f.get_node("/phenotypes/%s/" % phenotype)
        if 'GxE' in g:
            for gxeResult in self.h5_f.iter_nodes(g.GxE, 'Group'):
                gxeResults.append({'phenotype':phenotype, 'type':gxeResult._v_name.title()})
        return gxeResults

    def _get_phenotype_environments_(self, phenotype):
        phenotype_group = self.h5_f.get_node('/phenotypes/', phenotype)
        environments = []
        for environment in self.h5_f.iter_nodes("/phenotypes/%s" % phenotype, 'Group'):
            if environment._v_name != 'GxE':
                environments.append({'phenotype':phenotype, 'name':environment._v_name, 'datasets':self._get_environment_datasets_(phenotype, environment._v_name)})
        return environments

    def _get_environment_datasets_(self, phenotype, environment):
        datasets = []
        for dataset in self.h5_f.iter_nodes('/phenotypes/%s/%s' % (phenotype, environment), 'Group'):
            datasets.append({'phenotype':phenotype, 'environment':environment, 'name':dataset._v_name, 'transformations':self._get_dataset_transformations_(phenotype, environment, dataset._v_name)})
        return datasets

    def _get_cissvstrans_stats(self, transformation, type):
        cvt_stats = []
        if 'cvt' not in transformation:
            return cvt_stats
        table = self.h5_f.get_node(transformation, 'cvt')
        e_type = cvt_types[type]
        result = table.read_where('(cvt_type==e_type)')
        for cvt in result:
            pval = float(cvt['pval'])
            if pval > 0.0:
                pval = round(float(-math.log10(cvt['pval'])), 2)
            cvt_stats.append({'type':type, 'percVarGlobal':round(float(cvt['perc_var_global']), 4)
                               , 'percVarLocal':round(float(cvt['perc_var_local']), 4), 'distance':int(cvt['position']),
                               'pseudoHeritabilityGlobal':round(float(cvt['pseudo_heritability_global']), 4), 'pseudoHeritabilityGlobalLocal':round(float(cvt['pseudo_heritability_local']), 4)
                               , 'score':pval})
        return cvt_stats

    def _get_dataset_transformations_(self, phenotype, environment, dataset):
        transformations = []
        for transformation in self.h5_f.iter_nodes('/phenotypes/%s/%s/%s' % (phenotype, environment, dataset), 'Group'):
            tssUpstream = self._get_cissvstrans_stats(transformation, "tss_upstream")
            radius = self._get_cissvstrans_stats(transformation, "radius")
            cofactors = []
            if "EX" in transformation:
                cofactors = self._get_step_wise_statistics(self.h5_f.get_node(transformation, 'EX'))
            transformations.append({'phenotype':phenotype, 'environment':environment, 'dataset':dataset, 'name':transformation._v_name, 'results':self._get_transformation_results_(phenotype, environment, dataset, transformation._v_name), 'cofactors':cofactors, 'tssUpstream':tssUpstream, 'radius':radius})
        return transformations

    def _get_transformation_results_(self, phenotype, environment, dataset, transformation):
        results = []
        stats = []
        for result in self.h5_f.iter_nodes('/phenotypes/%s/%s/%s/%s' % (phenotype, environment, dataset, transformation), 'Group'):
            results.append({'phenotype':phenotype, 'environment':environment, 'dataset':dataset, 'transformation':transformation, 'name':result._v_name})
            if (result._v_name == 'EX'):
                results.extend(self._get_step_wise_results(phenotype, environment, dataset, transformation, result))
        return results

    def _get_step_wise_results(self, phenotype, environment, dataset, transformation, result):
        results = []
        step = 0
        for ex_step in self.h5_f.iter_nodes(result, 'Table'):
            if (ex_step._v_name != 'results'):
                step = int(ex_step._v_name[-1])
                results.append({'phenotype':phenotype, 'environment':environment, 'dataset':dataset, 'transformation':transformation, 'name':ex_step._v_name , 'step':step})
        return results

    def _get_step_wise_statistics(self, result):
        stats_to_return = ['chr', 'pos', 'bic', 'ebic', 'mbic', 'step', 'pseudo_heritability', 'max_cof_pval']
        stat_headers = result._v_attrs['stepwise_header']
        stats = []
        for stat in result._v_attrs['stepwise_stats']:
            cofactors = stat[stat_headers.index('cofactors')]
            if len(cofactors) == 0:
                chr = None
                pos = None
            else:
                cofactor = cofactors[-1]
                chr = int(cofactor[0])
                pos = int(cofactor[1])
            mbonf = stat[stat_headers.index('mbonf')]
            if mbonf > 0 :
                mbonf = -math.log10(mbonf)
            stats.append({'chr':chr, 'pos':pos, 'step':len(cofactors), 'mbonf':round(mbonf, 4), 'bic':round(stat[stat_headers.index('bic')], 4), 'mbic':round(stat[stat_headers.index('m_bic')], 4), 'ebic':round(stat[stat_headers.index('e_bic')], 4), 'pseudo_heritability':round(stat[stat_headers.index('pseudo_heritability')], 4)})
        return stats

    def _get_phenotype_bins(self, phen_vals, bin_number, phen2_vals):
        min_phen_val = min(phen_vals)
        max_phen_val = max(phen_vals)
        if phen2_vals is not None:
            min_phen2_val = min(phen2_vals)
            max_phen2_val = max(phen2_vals)
            if min_phen2_val < min_phen_val:
                min_phen_val = min_phen2_val
            if max_phen2_val > max_phen_val:
                max_phen_val = max_phen2_val
        bin_counts, edges = sp.histogram(phen_vals, bins=bin_number, range=(min_phen_val, max_phen_val))
        if phen2_vals is not None:
            bin2_counts, edges2 = sp.histogram(phen2_vals, bins=bin_number, range=(min_phen_val, max_phen_val))
        bin_list = []
        for i, bin_count in enumerate(bin_counts):
            x1 = edges[i]
            x2 = edges[i + 1]
            d = {'x_axis':'%0.2f-%0.2f' % (x1, x2), 'frequency':bin_count}
            if phen2_vals is not None:
                d['frequency_alt'] = bin2_counts[i]
            bin_list.append(d)
        return bin_list


    def get_gxe_results_by_chromosome(self, phen_name, result, min_mac=0, min_score=0.0, \
                top_fraction=1, chromosomes=[1, 2, 3, 4, 5]):
        """
        Return results..
        """
        cd = {}
        #h5file = self._open(mode="r")

        try:
            if result != 'combined':
                d_keys = ['position', 'score', 'maf', 'mac']
                table = self.h5_f.get_node('/phenotypes/%s/GxE/%s/results' % (phen_name, result))
                for chromosome in chromosomes:
                    res = table.read_where('(chromosome==%s) & (score>=%f) & (mac >=%d)' % (chromosome, min_score, min_mac))
                    d = {}
                    for k in d_keys:
                        d[k] = res[k]
                    cd[chromosome] = d
                cd['chromosome_ends'] = chromosome_ends
                cd['max_score'] = table._v_attrs['max_score']
            else:
                start = 0
                stop = 3333
                full_table = self.h5_f.get_node('/phenotypes/%s/GxE/full/results' % (phen_name))
                environ_table = self.h5_f.get_node('/phenotypes/%s/GxE/environ/results' % (phen_name))
                genetic_table = self.h5_f.get_node('/phenotypes/%s/GxE/genetic/results' % (phen_name))
                for chromosome in chromosomes:
                    full = full_table.read_where('(chromosome==%s) & (score>=%f) & (mac >=%d)' % (chromosome, min_score, min_mac))[start:stop]
                    genetic = genetic_table.read_where('(chromosome==%s) & (score>=%f) & (mac >=%d)' % (chromosome, min_score, min_mac))[start:stop]
                    environ = environ_table.read_where('(chromosome==%s) & (score>=%f) & (mac >=%d)' % (chromosome, min_score, min_mac))[start:stop]
                    data = numpy.lib.recfunctions.join_by('position', numpy.lib.recfunctions.join_by('position', full[['position', 'score']], genetic[['position', 'score']], jointype='outer'), environ[['position', 'score']], jointype='outer')
                    d = {}
                    d['position'] = data['position'].tolist()
                    d['full'] = data['score1'].tolist()
                    d['genetic'] = data['score2'].tolist()
                    d['environ'] = data['score'].tolist()
                    cd[chromosome] = d
                max_score = [full_table._v_attrs['max_score'], genetic_table._v_attrs['max_score'], environ_table._v_attrs['max_score']]
                max_score = max(max_score)

                cd['chromosome_ends'] = chromosome_ends
                cd['max_score'] = max_score
        except Exception, err:
            raise(err)
        finally:
            test = ''
            #self._close(h5file)
        return cd

    def get_results_by_chromosome(self, phen_name, environment, dataset, transformation, result, min_mac=0, min_score=0.0, \
                top_fraction=1, chromosomes=[1, 2, 3, 4, 5], log_transform=True):
        """
        Return results..
        """
        cd = {}
        #h5file = self._open(mode="r")
        try:
            cofactors_grouped = {}
            if result in ['EX', 'LM', 'KW']:
                table = self.h5_f.get_node('/phenotypes/%s/%s/%s/%s/%s/results' % (phen_name, environment, dataset, transformation, result))
            else:
                group = self.h5_f.get_node('/phenotypes/%s/%s/%s/%s/%s/' % (phen_name, environment, dataset, transformation, 'EX'))
                step = int(result[-1])
                table = self.h5_f.get_node(group, result)
                stepwise_header = group._v_attrs['stepwise_header']
                cofactors = group._v_attrs['stepwise_stats'][step][stepwise_header.index('cofactors')]
                cofactors.sort()
                from itertools import groupby
                for k, g in groupby(cofactors, lambda x: x[0]):
                    cofactors_grouped[k] = ([(item[1], item[2]) for item in g])

            d_keys = ['score', 'position', 'maf', 'mac']
            """ if analysis_method == 'kw':
                d_keys.append('statistic')
            else:
                d_keys.extend(['beta0', 'beta1', 'correlation', 'genotype_var_perc'])
            """
            for chromosome in chromosomes:
                cofactors_chr = []
                if chromosome in cofactors_grouped:
                    cofactors_chr = cofactors_grouped[chromosome]
                sort_list = []
                #for x in table.where('(chromosome==%d) &(score>=%f) & (mac>=%d)' % (chromosome, min_score, min_mac)):
                res = table.read_where('(chromosome==%s) & (score> %f) & (mac >=%d)' % (chromosome, min_score, min_mac))
                d = {}
                for k in d_keys:
                    d[k] = res[k]
                d['score'] = numpy.round(numpy.negative(numpy.log10(d['score'])), decimals=2)
                for cofactor in cofactors_chr:
                    d['score'] = numpy.append(d['score'], [cofactor[1]])
                    d['position'] = numpy.append(d['position'], [cofactor[0]])
                    d['mac'] = numpy.append(d['mac'], [0])
                    d['maf'] = numpy.append(d['maf'], [0])
                cd[chromosome] = d
            cd['chromosome_ends'] = chromosome_ends
            cd['max_score'] = -math.log10(table._v_attrs['max_score'])
        except Exception, err:
            raise(err)
        finally:
            test = ''
            #self._close(h5file)
        return cd



    def get_results_for_csv(self, phen_name, environment, dataset=None, transformation=None, result=None, min_mac=0):
        """
        Return results..
        """
        header = ['chromosome', 'position', 'score', 'mac', 'maf', 'perc_var_expl']
        if environment != 'GxE':
            if result in ['LM', 'KW', 'EX']:
                table = self.h5_f.get_node('/phenotypes/%s/%s/%s/%s/%s/results' % (phen_name, environment, dataset, transformation, result))
            else:
                group = self.h5_f.get_node('/phenotypes/%s/%s/%s/%s/%s/' % (phen_name, environment, dataset, transformation, 'EX'))
                step = int(result[-1])
                table = self.h5_f.get_node(group, result)
        else:
            table = self.h5_f.get_node('/phenotypes/%s/%s/%s/results' % (phen_name, environment, result.lower()))
        try:
            data = []
            data.append(header)
            res = table.read_where('mac >=%f' % min_mac)
            for row in res:
                data.append([row[k] for k in header])
        except Exception, err:
            raise(err)
        finally:
            test = ''
            #self._close(h5file)
        return data

    def importAccessions(self, csv_file, delimiter='\t'):
        csv_f = None
        try:
            if 'accessions' not in self.h5_f.root:
                group = self.h5_f.create_group("/", "accessions", "Accessions")
            else:
                group = self.h5_f.get_node('/', 'accessions')
                self.h5_f.remove_node(group, 'infos', True)
            csv_f = csv.reader(open(csv_file, 'rb'), delimiter=delimiter)
            table = self.h5_f.create_table(group, 'infos', Accessions, 'Accessions Infos')
            datasets = table.get_enum('dataset')
            for accession in csv_f:
                row = table.row
                row['accession_id'] = accession[0]
                row['name'] = accession[1]
                if accession[2] != "\\N":
                    row['longitude'] = accession[2]
                if accession[3] != "\\N":
                    row['latitude'] = accession[3]
                row['country'] = accession[4]
                row['collector'] = accession[5]
                if accession[6] != "\\N":
                    unix_time = int(time.mktime(datetime.strptime(accession[6], "%Y-%m-%d %H:%M:%S").timetuple()))
                    print unix_time
                row['collection_date'] = unix_time
                row['dataset'] = datasets[accession[7]]
                row.append()
            table.flush()
            print "%s accessions successuflly imported" % table.nrows
        except Exception, err:
            print str(err)


    def importPhenotypes(self, phenotype_file, environment):
        phenotypes = loadPhenotypeFile(phenotype_file)
        if environment not in ['T10C', 'T16C']:
            raise Exception('% environment unknown' % environment)
        filters = tables.Filters(complevel=5, complib='blosc')
        if 'phenotypes' not in self.h5_f.root:
            group = self.h5_f.create_group('/', 'phenotypes', 'Phenotype Group')
        else:
            group = self.h5_f.root.phenotypes
        #pdb.set_trace()
        for phenotype in phenotypes:
            phen_name = phenotype['name']
            if phen_name not in group:
                phen_group = self.h5_f.create_group(group, phen_name, 'Phenotype')
            else:
                phen_group = self.h5_f.get_node(self.h5_f.root.phenotypes, phen_name)
            if environment not in phen_group:
                env_group = self.h5_f.create_group(phen_group, environment, 'Environment %s ' % environment)
            else:
                env_group = self.h5_f.get_node(phen_group, environment)
            if 'Fullset' not in env_group:
                dataset_group = self.h5_f.create_group(env_group, 'Fullset', 'Fullset')
            else:
                dataset_group = self.h5_f.get_node(env_group, 'Fullset')
            if 'raw' not in dataset_group:
                trans_group = self.h5_f.create_group(dataset_group, 'raw', 'Raw transformation')
            else:
                trans_group = self.h5_f.get_node(dataset_group, 'raw')
            if 'values' in trans_group:
                table = self.h5_f.get_node(trans_group, 'values')
                if table.nrows != len(phenotype['values']):
                    table.remove_rows(0, table.nrows)
                else:
                    continue
            else:
                table = self.h5_f.create_table(trans_group, 'values', PhenotypeValue, filters=filters)
            phen_value = table.row
            for phenotype_value in phenotype['values']:
                phen_value['ecotype'] = phenotype_value[0]
                phen_value['value'] = phenotype_value[1]
                phen_value.append()
            table.flush()
        self.h5_f.flush()



if __name__ == '__main__':
    usage = "usage: %prog [options] action"
    parser = OptionParser(usage=usage)
    parser.add_option("-f", "--file", dest="hdf5_filename", help="hdf5 file", metavar="FILE")
    parser.add_option("-p", "--phenotype_file", dest="phenotype_filename", help="Phenotype file name (use only with --phenotype action)", metavar="FILE")
    parser.add_option("-r", "--results_folder", dest="results_folder", help="Results folder (use only with --results action)", metavar="FOLDER")
    parser.add_option("-a", "--accessions_file", dest="accessions_file", help="Accessions file (use only with --accessions action)", metavar="FILE")
    parser.add_option("-e", "--environment", dest="environment", help="Environment (T10C or T16C)")
    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("specify the action (phenotype or results)")
    elif not options.hdf5_filename:
        parser.error("you have to specify the hdf5_filename")
    elif args[0] == "phenotype" and not options.phenotype_filename:
        parser.error("you have to specify the phenotype_filename for the phenotype action")
    elif args[0] == "results" and not options.results_folder:
        parser.error("you have to specify the results_folder for the results action")
    elif args[0] == "accessions" and not options.accessions_file:#
        parser.error("you have to specify the accessions_file for the accessions action")
    elif options.environment != "T10C"  and options.environment != "T16C" and args[0] != 'accessions':
        parser.error("Environment has to be T10C or T16C")
    else:
        rnaseq_records = RNASeqRecords(options.hdf5_filename, "a")
        if args[0] == 'phenotype':
            rnaseq_records.importPhenotypes(options.phenotype_filename, options.environment)
        elif args[0] == 'results':
            rnaseq_records.importGWASResults(options.results_folder, options.environment)
        else:
            rnaseq_records.importAccessions(options.accessions_file)
        rnaseq_records.close()
    parser.destroy()

