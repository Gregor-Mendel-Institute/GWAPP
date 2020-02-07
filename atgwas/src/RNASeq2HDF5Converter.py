#!/usr/bin/env python

import tables
import csv
import re
import cPickle
import pdb
import os
import itertools
from optparse import OptionParser


gene_symbol = re.compile(r'^AT[1-5]{1}[GE]\d+$')

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
    chromosome = tables.UInt8Col()
    position = tables.Int32Col()
    score = tables.Float32Col()
    maf = tables.Float32Col()
    mac = tables.Int32Col()
    perc_var_expl = tables.Float32Col()


class ResultRecordLM(ResultRecord):
    """
    Linear model, mixed models, etc.
    """
    genotype_var_perc = tables.Float32Col()
    beta0 = tables.Float32Col()
    beta1 = tables.Float32Col()
    correlation = tables.Float32Col()

cvt_types = tables.Enum(['radius', 'tss_upstream'])

class CVT(tables.IsDescription):
    cvt_type = tables.EnumCol(cvt_types, 'radius', base='uint8')
    position = tables.Int32Col()
    perc_var_global = tables.Float32Col()
    perc_var_local = tables.Float32Col()
    pseudo_heritability_global = tables.Float32Col()
    pseudo_heritability_local = tables.Float32Col()
    pval = tables.Float32Col()

class GxERecord(tables.IsDescription):
    chromosome = tables.UInt8Col()
    position = tables.Int32Col()
    score = tables.Float32Col()
    maf = tables.Float32Col()
    mac = tables.Int32Col()
    perc_var_expl = tables.Float32Col()
    betas = tables.Float32Col(shape=(3,))

class GxERecordFull(tables.IsDescription):
    chromosome = tables.UInt8Col()
    position = tables.Int32Col()
    score = tables.Float32Col()
    maf = tables.Float32Col()
    mac = tables.Int32Col()
    perc_var_expl = tables.Float32Col()
    betas = tables.Float32Col(shape=(4,))

class Phenotype(tables.IsDescription):
    id = tables.Int32Col()
    chr = tables.UInt8Col()
    start = tables.Int32Col()
    end = tables.Int32Col()
    name = tables.StringCol(itemsize=20)
    min_score_10C = tables.Float32Col()
    min_score_16C = tables.Float32Col()
    min_score_full = tables.Float32Col()
    pseudo_heritability_10C = tables.Float32Col()
    pseudo_heritability_16C = tables.Float32Col()


result_types = tables.Enum(['EX', 'LM', 'KW', 'full', 'genetic', 'environ'])
environ_types = tables.Enum(['10', '16', 'GxE'])

METADATA_CACHE_SIZE = 2 * 1024 * 1024
NODE_CACHE_SLOTS = 1024


class TopSNPResults(tables.IsDescription):
    snp_chr = tables.UInt8Col()
    snp_pos = tables.UInt32Col()
    score = tables.Float32Col()
    maf = tables.Float32Col()
    mac = tables.Int32Col()
    perc_var_expl = tables.Float32Col()
    gene = tables.StringCol(itemsize=20)
    gene_chr = tables.UInt8Col()
    gene_mid_pos = tables.UInt32Col()
    gene_start = tables.UInt32Col()
    gene_end = tables.UInt32Col()
#    result_type = tables.EnumCol(result_types,'EX',base='uint8')
#    environ_type = tables.EnumCol(environ_types,'10',base='uint8')


def loadPhenotypeFile(phenotype_file, gene_file):
    reader = csv.reader(open(phenotype_file, 'rb'), delimiter=',')
    header = reader.next()
    last_pid = -1
    phenotypes = []
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

def _getResultsFromFile(result_file):
    f = open(result_file, 'r')
    try:
        result_data = cPickle.load(f)
    except Exception, err:
        result_data = None
    finally:
        f.close()
    return result_data


def generateTopResultsTable(hdf5_file, min_score=5, environ_types=['T10C', 'T16C', 'GxE'],
                result_types=['EX', 'LM', 'KW', 'full', 'genetic', 'environ'], replace=False):
    import math
    errors = []
    filters = tables.Filters(complevel=4, complib='blosc')
    h5_f = tables.openFile(hdf5_file, 'r+', NODE_CACHE_SLOTS=0, METADATA_CACHE_SIZE=METADATA_CACHE_SIZE)
    group_name = 'top_results'
#    pdb.set_trace()
    if 'phenotypes' not in h5_f.root:
        raise Exception('First run importPhenotypes')
    else:
        group = h5_f.root
    if 'phenotype_list' not in group:
        raise Exception('First run generatePhenotypeList method')

    if 'top_results' in h5_f.root:
        group = h5_f.getNode('/', group_name)
    else:
        group = h5_f.createGroup('/', group_name)
    def _getRow(gene_name, gene_chr, gene_mid_pos, result_type, environ_type, snp, convert=True):
        score = snp['score']
        if convert and score > 0:
            score = -math.log10(score)

        row = (gene_name, gene_chr, gene_end, gene_mid_pos, gene_start, snp['mac'], snp['maf'], \
              snp['perc_var_expl'], score, snp['chromosome'], snp['position'])
        return row

    environ_ix = 0
    environ_count = len(environ_types)
    result_count = len(result_types)
    phen_count = len(h5_f.root.phenotype_list)
    for environ_type in environ_types:
        environ_ix = environ_ix + 1
        convert = False
        if environ_type != 'GxE':
            convert = True
        if environ_type in group:
            environ_group = h5_f.getNode(group, environ_type)
        else:
            environ_group = h5_f.createGroup(group, environ_type)
        result_ix = 0
        for result_type in result_types:
            result_ix = result_ix + 1
            if environ_type == 'GxE' and result_type not in ['full', 'genetic', 'environ']:
                continue
            elif environ_type != 'GxE' and result_type not in ['EX', 'LM', 'KW']:
                continue
            table = None
            if result_type in environ_group:
                table = h5_f.getNode(environ_group, result_type)
                if replace:
                    h5_f.removeNode(environ_group, result_type)
                    table = None
                elif table.nrows > 0:
                    continue

            rows_to_add = []
            phen_ix = 0
            for phenotype in h5_f.root.phenotype_list:
                try:
 #                   if phenotype['name'] == 'AT1G01040':
  #                      pdb.set_trace()
                #    if phen_ix == 3: 
                 #       break
                    phen_ix = phen_ix + 1
                    gene_name = phenotype['name']
                    gene_mid_pos = int(round((phenotype['end'] - phenotype['start']) / 2))
                    gene_start = phenotype['start']
                    gene_end = phenotype['end']
                    gene_chr = phenotype['chr']
                    if environ_type == 'GxE':
                        res_table = h5_f.getNode('/phenotypes/%s/GxE/%s/results' % (gene_name, result_type))
                    else:
                        res_table = h5_f.getNode('/phenotypes/%s/%s/Fullset/raw/%s/results' % (gene_name, environ_type, result_type))
                    results = _getTopSNPs(res_table, min_score, convert)
                    for snp in results:
                        rows_to_add.append(_getRow(gene_name, gene_chr, gene_mid_pos, gene_start, gene_end, snp, convert))
                except Exception, err:
                    errors.append({'environ':environ_type, 'result':result_type, 'phenotype':gene_name, 'error':str(err)})
                print('Environment %s (%s/%s) | Result %s (%s/%s) | Phenotype %s (%s/%s) loaded' % (environ_type, environ_ix, environ_count, result_type, result_ix, result_count, gene_name, phen_ix, phen_count))
          #  pdb.set_trace()
            rows_to_add.sort(key=lambda x : x[8] * -1)
            if table == None:
                table = h5_f.createTable(environ_group, result_type, TopSNPResults, 'Top SNPs from %s/%s' % (environ_type, result_type), filters=filters, expectedrows=len(rows_to_add))
            table.append(rows_to_add)
            print('Environment %s (%s/%s) | Result %s (%s/%s) saved' % (environ_type, environ_ix, environ_count, result_type, result_ix, result_count))
            table._v_attrs.min_score = min_score
            table.flush()
            table.cols.snp_chr.createIndex()
            table.cols.snp_pos.createIndex()
            table.cols.score.createIndex()
            table.cols.maf.createIndex()
            table.cols.mac.createIndex()
            table.cols.perc_var_expl.createIndex()
            table.cols.gene.createIndex()
            table.cols.gene_chr.createIndex()
            table.cols.gene_mid_pos.createIndex()
            table.cols.gene_start.createIndex()
            table.cols.gene_end.createIndex()
            table = None
    h5_f.flush()
    h5_f.close()
    for error in errors:
        print('Environment %s | Result  %s | Phenotype %s : %s' % (error['environ'], error['result'], error['phenotype'], error['error']))



def _getTopSNPs(table, min_score=5, convert=True):
    import math
    if convert:
        min_score = math.pow(10, -min_score)
        results = table.readWhere('(score <= min_score)')
    else:
        results = table.readWhere('(score >= min_score)')
    return results


def generatePhenotypeList(hdf5_file, gene_file):
    import math
    errors = []
    table_name = 'phenotype_list'
    h5_f = tables.openFile(hdf5_file, 'r+', NODE_CACHE_SLOTS=NODE_CACHE_SLOTS, METADATA_CACHE_SIZE=METADATA_CACHE_SIZE)
    gene_lookup = _getResultsFromFile(gene_file)
    if gene_lookup is None:
        raise Exception('Could not load gene_file %s' % gene_file)
    if 'phenotypes' not in h5_f.root:
        raise Exception('First run importPhenotypes')
    else:
        group = h5_f.root.phenotypes

    if table_name in h5_f.root:
        table = h5_f.getNode("/", table_name)
        h5_f.removeNode(table)
        table = None
    pdb.set_trace()
    table = h5_f.createTable("/", table_name, Phenotype, 'Phenotype list')
    phen_value = table.row
    i = 1
    for phenotype in h5_f.iterNodes(group, 'Group'):
        try:
            phen_name = phenotype._v_name
            if phen_name not in gene_lookup:
                raise Exception('phenotype %s not found in lookup' % phen_name)
            phen_group = h5_f.getNode('/phenotypes/%s' % phen_name);
            min_score_10C = 0.0
            min_score_16C = 0.0
            pseudo_heritability_10C = 0.0
            pseudo_heritability_16C = 0.0
            min_score_full = 0.0
            try:
                min_score_10C = -math.log10(min(phen_group.T10C.Fullset.raw.EX.results[:]['score']))
                min_score_16C = -math.log10(min(phen_group.T16C.Fullset.raw.EX.results[:]['score']))
                pseudo_heritability_10C = phen_group.T10C.Fullset.raw.EX._v_attrs['stepwise_stats'][0][4]
                pseudo_heritability_16C = phen_group.T16C.Fullset.raw.EX._v_attrs['stepwise_stats'][0][4]
                min_score_full = max(phen_group.GxE.full.results[:]['score'])
            except:
                test = ''

            gene = gene_lookup[phen_name]
            phen_value['id'] = i
            phen_value['chr'] = gene['chromosome']
            phen_value['name'] = phen_name
            phen_value['start'] = gene['start_pos']
            phen_value['end'] = gene['end_pos']
            phen_value['min_score_10C'] = min_score_10C
            phen_value['min_score_16C'] = min_score_16C
            phen_value['min_score_full'] = min_score_full
            phen_value['pseudo_heritability_10C'] = pseudo_heritability_10C
            phen_value['pseudo_heritability_16C'] = pseudo_heritability_16C
            phen_value.append()
            i = i + 1
        except Exception, err:
            errors.append({'phenotype':phen_name, 'error':str(err)})

    table.flush()
    table.cols.chr.createIndex()
    table.cols.id.createIndex()
    table.cols.start.createIndex()
    table.cols.end.createIndex()
    table.cols.name.createIndex()

    table.cols.min_score_10C.createIndex()
    table.cols.min_score_16C.createIndex()
    table.cols.min_score_full.createIndex()
    table.cols.pseudo_heritability_10C.createIndex()
    table.cols.pseudo_heritability_16C.createIndex()
    h5_f.flush()
    h5_f.close()
    for error in errors:
        print('%s : %s' % (error['phenotype'], error['error']))

def importPhenotypes(hdf5_file, phenotype_file, environment):
    h5_f = tables.openFile(hdf5_file, 'a', NODE_CACHE_SLOTS=NODE_CACHE_SLOTS, METADATA_CACHE_SIZE=METADATA_CACHE_SIZE)
    p_f = open(phenotype_file)
    header = p_f.next()
    phenotypes = loadPhenotypeFile(phenotype_file)
    if environment not in ['T10C', 'T16C']:
        raise Exception('% environment unknown' % environment)
    filters = tables.Filters(complevel=5, complib='blosc')
    if 'phenotypes' not in h5_f.root:
        group = h5_f.createGroup('/', 'phenotypes', 'Phenotype Group')
    else:
        group = h5_f.root.phenotypes
    #pdb.set_trace()
    for phenotype in phenotypes:
        phen_name = phenotype['name']
        if phen_name not in group:
            phen_group = h5_f.createGroup(group, phen_name, 'Phenotype')
        else:
            phen_group = h5_f.getNode(h5_f.root.phenotypes, phen_name)
        if environment not in phen_group:
            env_group = h5_f.createGroup(phen_group, environment, 'Environment %s ' % environment)
        else:
            env_group = h5_f.getNode(phen_group, environment)
        if 'Fullset' not in env_group:
            dataset_group = h5_f.createGroup(env_group, 'Fullset', 'Fullset')
        else:
            dataset_group = h5_f.getNode(env_group, 'Fullset')
        if 'raw' not in dataset_group:
            trans_group = h5_f.createGroup(dataset_group, 'raw', 'Raw transformation')
        else:
            trans_group = h5_f.getNode(dataset_group, 'raw')

        if 'values' in trans_group:
            table = h5_f.getNode(trans_group, 'values')
            if table.nrows != len(phenotype['values']):
                table.removeRows(0, table.nrows)
            else:
                continue
        else:
            table = h5_f.createTable(trans_group, 'values', PhenotypeValue, filters=filters)
        phen_value = table.row
        for phenotype_value in phenotype['values']:
            phen_value['ecotype'] = phenotype_value[0]
            phen_value['value'] = phenotype_value[1]
            phen_value.append()
        table.flush()
    h5_f.flush()
    h5_f.close()






def importGWASResults(hdf5_file, gwas_results_folder, env, type='all', replace=False):
    errors = []
    filters = tables.Filters(complevel=4, complib='blosc')
    h5_f = tables.openFile(hdf5_file, 'r+', NODE_CACHE_SLOTS=NODE_CACHE_SLOTS, METADATA_CACHE_SIZE=METADATA_CACHE_SIZE)
    if 'phenotypes' not in h5_f.root:
        raise Exception('First run importPhenotypes')
    if env not in ['T10C', 'T16C']:
        raise Exception('%s environment unknown' % env)
    phen_group = h5_f.root.phenotypes
    chr_folder = gwas_results_folder + '/chr_%s'
    #pdb.set_trace()
    for chr_index in range(1, 6):
        files = os.listdir(chr_folder % chr_index)
        files.sort()
        for file in files:
            try:
                ex_step = 0
                parts = file.split('_')
                phen_name = parts[-3]
                analysis = parts[-2]
                if len(parts) == 8:
                    phen_name = parts[-2]
                    analysis = 'stats'
                if analysis not in ['EX', 'LM', 'KW', 'ex', 'stats']:
                    raise Exception('%s analysis not supported' % analysis)
                if analysis == 'ex':
                    ex_step = int(parts[-1].split('.')[0][-1])
                    analysis = 'EX'
                if gene_symbol.search(phen_name) is  None:
                    raise Exception('%s has unknown phenotype name' % file)
                if type != 'all'  and type != analysis:
                    continue
                trans_group = h5_f.getNode('/phenotypes/%s/%s/%s/%s' % (phen_name, env, 'Fullset', 'raw'))
                if analysis in ['EX', 'LM', 'KW']:
                    if analysis not in trans_group:
                        analysis_group = h5_f.createGroup(trans_group, analysis, '%s Analysis' % analysis)
                    else:
                        analysis_group = h5_f.getNode(trans_group, analysis)

                    dataset_name = 'results'
                    if analysis == 'EX' and ex_step > 0:
                        dataset_name = '%s%i' % (dataset_name, ex_step)
                    if dataset_name not in analysis_group:
                        table = h5_f.createTable(analysis_group, dataset_name, ResultRecord, dataset_name, filters=filters, expectedrows=35000)
                    else:
                        table = h5_f.getNode(analysis_group, dataset_name)
                    if table.nrows > 0:
                        continue
                    result = table.row
                    result_data = _getResultsFromFile('%s/chr_%s/%s' % (gwas_results_folder, chr_index, file))
                    if result_data == None:
                        h5_f.removeNode(analysis_group, dataset_name)
                        raise Exception('%s failed to retrieve result data' % file)
                    snp_results = result_data['snp_results']
                    max_score = min(snp_results['scores'])
                    table._v_attrs.max_score = max_score
                    data = zip(snp_results['chromosomes'], snp_results['positions'], snp_results['scores'], snp_results['mafs'], snp_results['macs'], snp_results['perc_var_expl'])
                    data.sort(key=lambda k: (k[0], k[2]))
                    for cpsmm in data:
                        (result['chromosome'], result['position'], result['score'], result['maf'], result['mac'], result['perc_var_expl']) = cpsmm
                        result.append()
                    table.flush()
                    table.cols.chromosome.createIndex()
                    table.cols.score.createIndex()
                    table.cols.mac.createIndex()
                    print("%s results (%s) imported" % (phen_name, analysis))
                elif analysis in ['stats']:
                    #i analysis_group._v_nchildren > 0:
                     #  continue
                    result_data = _getResultsFromFile('%s/chr_%s/%s' % (gwas_results_folder, chr_index, file))
                    if result_data == None:
                        raise Exception('%s failed to retrieve result data' % file)
                    #trans_group._v_attrs.transformation_type = result_data['transformation_type']
                    trans_group._v_attrs.pseudo_heritability = result_data['pseudo_heritability']
                    #trans_group._v_attrs.transformation_shapiro_pval = result_data['transformation_shapiro_pval']

                    # store GWAS statistics
                    statistics = ['LM', 'KW', 'EX', 'SW']
                    for stats in statistics:
                        if stats not in result_data:
                            continue
                        if stats != 'SW':
                            if  stats not in trans_group:
                                analysis_stats_group = h5_f.createGroup(trans_group, stats, '%s Analysis' % stats)
                            else:
                                analysis_stats_group = h5_f.getNode(trans_group, stats)
                            analysis_stats_group._v_attrs.statistics = result_data[stats]
                        else:
                            if 'EX' not in trans_group:
                                analysis_stats_group = h5_f.createGroup(trans_group, 'EX', '%s Analysis' % stats)
                            else:
                                analysis_stats_group = h5_f.getNode(trans_group, 'EX')
                            step_wise_header = ['m_bic', 'e_bic', 'bic', 'mbonf', 'pseudo_heritability', 'cofactors']
                            step_wise_stats = []
                            for step_stats in result_data[stats]['step_info_list']:
                                stat_to_add = [step_stats[key] for key in step_wise_header]
                                step_wise_stats.append(stat_to_add)
                            analysis_stats_group._v_attrs.stepwise_header = step_wise_header
                            analysis_stats_group._v_attrs.stepwise_stats = step_wise_stats
                            analysis_stats_group._v_attrs.ebics = result_data[stats]['ebics']
                    # store CVT statistics
                    cvt_stats = result_data['CVT']
                    if 'cvt' in trans_group:
                        cvt_table = h5_f.getNode(trans_group, 'cvt')
                        if replace:
                            h5_f.removeNode(trans_group, 'cvt')
                            cvt_table = None
                        elif cvt_table.nrows > 0:
                            continue
                    if 'cvt' not in trans_group:
                        #pdb.set_trace()
                        cvt_table = h5_f.createTable(trans_group, 'cvt', CVT, 'cis_vs_trans', filters=filters)
                        #pdb.set_trace()
                        cvt_stat_row = cvt_table.row
                        for cvt_stat_key in cvt_stats.keys():
                            cvt_stat = cvt_stats[cvt_stat_key]
                            positions = cvt_stat.keys()
                            positions.sort()
                            cvt_type = cvt_types[cvt_stat_key]
                            for position  in positions:
                                cvt_stat_values = cvt_stat[position]
                                if cvt_stat_values is None:
                                    continue
                                cvt_stat_row['cvt_type'] = cvt_type
                                cvt_stat_row['position'] = position
                                cvt_stat_row['perc_var_global'] = cvt_stat_values['perc_var1']
                                cvt_stat_row['perc_var_local'] = cvt_stat_values['perc_var2']
                                cvt_stat_row['pseudo_heritability_global'] = cvt_stat_values['pseudo_heritability0']
                                cvt_stat_row['pseudo_heritability_local'] = cvt_stat_values['pseudo_heritability1']
                                cvt_stat_row['pval'] = cvt_stat_values['pval']
                                cvt_stat_row.append()
                        cvt_table.flush()
                    print("%s  Stats imported" % phen_name)
            except Exception, err:
                errors.append({'file':file, 'error':str(err)})
            #break
        #break
    h5_f.flush()
    h5_f.close()
    for error in errors:
        print('%s : %s\n' % (error['file'], error['error']))


def importGxEResults(hdf5_file, folder, gxe_type_to_import):
    errors = []
    tb_descs = {'full':GxERecordFull, 'genetic':GxERecord, 'environ':ResultRecord}
    filters = tables.Filters(complevel=4, complib='blosc')
    h5_f = tables.openFile(hdf5_file, 'r+', NODE_CACHE_SLOTS=NODE_CACHE_SLOTS, METADATA_CACHE_SIZE=METADATA_CACHE_SIZE)
    if 'phenotypes' not in h5_f.root:
        raise Exception('First run importPhenotypes')
    phenotypes_group = h5_f.root.phenotypes
    chr_folder = folder + '/chr_%s'
    for chr_index in range(1, 6):
        files = os.listdir(chr_folder % chr_index)
        files.sort()
        for file in files:
#            pdb.set_trace()
            try:
                parts = file.split('_')
                if len(parts) == 10:
                    phen_name = parts[-4]
                    gxe_type = 'environ'
                elif len(parts) == 8:
                    phen_name = parts[-2]
                    gxe_type = 'stats'
                else:
                    phen_name = parts[-3]
                    if parts[-2] == 'g':
                        gxe_type = 'genetic'
                    else:
                        gxe_type = 'full'
                #if phen_name != 'AT1G01010':
                #    exit()
                if gxe_type_to_import != 'all' and gxe_type_to_import != gxe_type:
                    continue

                if gxe_type not in ['full', 'genetic', 'environ', 'stats']:
                    raise Exception('%s analysis not supported' % analysis)
                if gene_symbol.search(phen_name) is  None:
                    raise Exception('%s has unknown phenotype name' % file)
                phen_group = h5_f.getNode(phenotypes_group, phen_name)

                if 'GxE' not in phen_group:
                    gxe_group = h5_f.createGroup(phen_group, 'GxE', 'GxE Analysis')
                else:
                    gxe_group = h5_f.getNode(phen_group, 'GxE')
                if gxe_type in ['genetic', 'full', 'environ']:
                    if gxe_type not in gxe_group:
                        gxe_type_group = h5_f.createGroup(gxe_group, gxe_type, gxe_type)
                    else:
                        gxe_type_group = h5_f.getNode(gxe_group, gxe_type)
                    tb_desc = tb_descs[gxe_type]
                    if 'results' not in gxe_type_group:
                        table = h5_f.createTable(gxe_type_group, 'results', tb_desc, 'pValues', filters=filters, expectedrows=35000)
                    else:
                        table = h5_f.getNode(gxe_type_group, 'results')

                    if table.nrows > 0:
                        continue
                    else:
                        h5_f.removeNode(gxe_type_group, 'results')
                        table = h5_f.createTable(gxe_type_group, 'results', tb_desc, 'pValues', filters=filters, expectedrows=35000)
                    result = table.row
                    result_data = _getResultsFromFile('%s/chr_%s/%s' % (folder, chr_index, file))
                    if result_data == None:
                        h5_f.removeNode(gxe_group, gxe_type)
                        raise Exception('%s failed to retrieve result data' % file)
                    snp_results = result_data['snp_results']
                    max_score = max(snp_results['scores'])
                    table._v_attrs.max_score = max_score
                    if gxe_type != 'environ':
                        data = zip(snp_results['chromosomes'], snp_results['positions'], snp_results['scores'], snp_results['mafs'], snp_results['macs'], snp_results['betas'], snp_results['perc_var_expl'])
                    else:
                        data = zip(snp_results['chromosomes'], snp_results['positions'], snp_results['scores'], snp_results['mafs'], snp_results['macs'], snp_results['perc_var_expl'])
                    data.sort(key=lambda k: (k[0], k[2] * -1))
                    top_fraction = int(0.04 * len(data))
                    curr_count = 0
                    for cpsmm in data:
                        if gxe_type != 'environ':
                            (result['chromosome'], result['position'], result['score'], result['maf'], result['mac'], result['betas'], result['perc_var_expl']) = cpsmm
                        else:
                            (result['chromosome'], result['position'], result['score'], result['maf'], result['mac'], result['perc_var_expl']) = cpsmm
                        if curr_count <= top_fraction * cpsmm[0]:
                            result.append()
                            curr_count = curr_count + 1
                    table.flush()
                    table.cols.chromosome.createIndex()
                    table.cols.score.createIndex()
                    table.cols.mac.createIndex()
                    print ("%s Result (%s) finished" % (phen_name, gxe_type))
                elif gxe_type in ['stats']:
                    result_data = _getResultsFromFile('%s/chr_%s/%s' % (folder, chr_index, file))
                    if result_data == None:
                        raise Exception('%s failed to retrieve result data' % file)
                    statistics = {'g_test':'genetic', 'gt_g_test':'environ', 'gt_test':'full'}

                    # store phenotype relevant stats
                    gxe_group._v_attrs.phen_info = result_data['phen_info']
                    gxe_group._v_attrs.log_quantiles = result_data['log_quantiles']
                    gxe_group._v_attrs.var_components = result_data['var_components']

                    quantiles = result_data['quantiles']
                    quantile_keys = {'g_test':'g', 'gt_test':'gt', 'gt_g_test':'gt_g'}
                    for stats in statistics.keys():
                        stat_k = statistics[stats]
                        if  stat_k not in gxe_group:
                            gxe_type_group = h5_f.createGroup(gxe_group, stat_k, '%s Stats' % stat_k)
                        else:
                            gxe_type_group = h5_f.getNode(gxe_group, stat_k)
                        gxe_type_group._v_attrs.statistics = result_data[stats]
                        quantile = quantiles[quantile_keys[stats]]
                        # quantiles
                        shape = (len(quantile),)
                        atom = tables.Float32Atom()
                        quantiles_array = h5_f.createCArray(gxe_type_group, 'quantiles', atom, shape, filters=filters)
                        quantiles_array[:] = quantile
                    print ("%s stats finished" % (phen_name))

                else:
                    raise Exception('Result file %s unknown' % file)
            except Exception, err:
                errors.append({'file':file, 'error':str(err)})
            #break
        #break
    h5_f.flush()
    h5_f.close()
    for error in errors:
        print('%s : %s\n' % (error['file'], error['error']))

def exportTopResultsToCSV(hdf5_file, environ_types=['T10C', 'T16C', 'GxE'], result_types=['EX', 'LM', 'KW', 'full', 'genetic', 'environ'], outputFolder="."):
    import csv
    pdb.set_trace()
    h5_f = tables.openFile(hdf5_file, 'r')
    fieldnames = ['gene', 'gene_mid_pos', 'gene_chr', 'snp_chr', 'snp_pos', 'score', 'mac', 'maf', 'perc_var_expl']

    if 'top_results' not in h5_f.root:
        raise Exception('Run top_results action first')
    for environ_type in environ_types:
        for result_type in result_types:
            if environ_type == 'GxE' and result_type not in ['full', 'genetic', 'environ']:
                continue
            elif environ_type != 'GxE' and result_type not in ['EX', 'LM', 'KW']:
                continue
            try:
                table = h5_f.getNode("/top_results/%s/%s" % (environ_type, result_type))
                writer = csv.writer(open("%s/%s_%s.csv" % (outputFolder, environ_type, result_type), 'a'), delimiter=',')
                writer.writerow(fieldnames)
                for result in table:
                    row = [result[key] for key  in fieldnames]
                    writer.writerow(row)
            except Exception, err:
                print "Error %s" % str(err)
    h5_f.close()


if __name__ == '__main__':
    usage = "usage: %prog [options] action"
    parser = OptionParser(usage=usage)
    parser.add_option("-f", "--file", dest="hdf5_filename", help="hdf5 file", metavar="FILE")
    parser.add_option("-p", "--phenotype_file", dest="phenotype_filename", help="Phenotype file name (use only with --phenotype action)", metavar="FILE")
    parser.add_option("-g", "--gene_file", dest="gene_filename", help="Gene annotation file name (use only with --phenotype action)", metavar="FILE")
    parser.add_option("-t", "--gxe_type", dest="gxe_type", help="Type of GxE to import (full, genetic, environ or stats. Default: all)", default='all')
    parser.add_option("-r", "--results_folder", dest="results_folder", help="Results folder (use only with --results action)", metavar="FOLDER")
    parser.add_option("-e", "--environment", dest="environment", help="Environment (T10C or T16C)")
    parser.add_option("-s", "--min_score", dest="min_score", help="Minimum score (-log10(x)) for generating top_results table", default=5)
    parser.add_option("-u", "--result_type", dest="result_type", help="For which result should the TOP-SNP Caching table generated (EX,LM,KW,full,genetic,environ)", default="EX,LM,KW,full,genetic,environ,genetic")
    parser.add_option("-c", "--replace", dest="replace", help="Replace existing records", action="store_true")
    parser.add_option("-o", "--output_folder", dest="output_folder", help="Output folder for CSV export")
    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("specify the action (phenotype or results or GxE or phenotypelist or top_results)")
    elif not options.hdf5_filename:
        parser.error("you have to specify the hdf5_filename")
    elif args[0] == "phenotype" and not options.phenotype_filename:
        parser.error("you have to specify the phenotype_filename for the phenotype action")
    elif args[0] in ["results", 'GxE']  and not options.results_folder and options.gxe_type not in ['environ', 'full', 'genetic', 'stats', 'all']:
        parser.error("you have to specify the results_folder for the results action and the gxe_type")
    elif args[0] == 'phenotypelist' and not options.gene_filename:
        parser.error("you have to specify the gene_filename  for the phenotypelist action")
    elif options.environment != "T10C"  and options.environment != "T16C" and args[0] == 'results':
        parser.error("Environment has to be T10C or T16C")
    else:
        if args[0] == 'phenotype':
            importPhenotypes(options.hdf5_filename, options.phenotype_filename, options.environment)
        elif args[0] == 'results':
            importGWASResults(options.hdf5_filename, options.results_folder, options.environment, options.gxe_type, options.replace)
        elif args[0] == 'GxE':
            importGxEResults(options.hdf5_filename, options.results_folder, options.gxe_type)
        elif args[0] == 'phenotypelist':
            generatePhenotypeList(options.hdf5_filename, options.gene_filename)
        elif args[0] == 'top_results':
            generateTopResultsTable(options.hdf5_filename, float(options.min_score), options.environment.split(','), options.result_type.split(','), options.replace)
        elif args[0] == 'top_results_csv':
            exportTopResultsToCSV(options.hdf5_filename, options.environment.split(','), options.result_type.split(','), options.output_folder)

    parser.destroy()
