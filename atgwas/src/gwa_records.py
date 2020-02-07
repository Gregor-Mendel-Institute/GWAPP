"""
A file which contains hdf5 interface for the phenotypes and results.

Overall data structure is: 

One hdf5 file per phenotype.
Three types of tables.
    - A info table, one for each transformation/accession subset.
    - phenotype table, one for each transformation/accession subset
    - result table, one for each transformation/accession subset and one for each analysis method.
    
    The minimum contains an info record, and a raw phenotype table.
"""
import tables
import phenotypeData as pd
import itertools
import scipy as sp
from scipy import stats
import util
import linear_models as lm
import dataParsers as dp
import math
import time
import pdb
import mtcorr

import analyze_gwas_results as agr



chromosome_ends = [30429953, 19701870, 23467451, 18578708, 26992130]


class ProgressFileWriter():
    """
    A small object to facilitate a progress bar update.
    """
    def __init__(self, file_name, file_mode='w'):
        self.file_mode = file_mode
        self.f = open(file_name, file_mode, 0)
        self.progress = 0.0
        self.task_status = ''
        self.step = 0.01

    def update_progress_bar(self, progress=None, step=None, task_status=''):
        assert self.file_mode == 'w', 'This requires file in write mode'
        if progress is not None:
            self.progress = progress
        elif step is not None:
            self.progress = self.progress + step
        else:
            self.progress = self.progress + self.step
        if self.progress > 1.0:
            self.progress = 1.0
        self.task_status = task_status
        self.f.write('%f\n%s\n' % (self.progress, self.task_status))

    def refresh_from_file(self):
        assert self.file_mode == 'r', 'This requires file in read mode'
        lines = self.f.readlines()
        if len(lines) >= 1:
            self.progress = lines[-2].strip()
        if len(lines) > 1:
            self.task_status = lines[-1].strip()

    def set_step(self, step):
        self.step = step

    def get_progress(self, from_file=True):
        if from_file:
            assert self.file_mode == 'r', 'This requires file in read mode'
            self.refresh_from_file()
        return self.progress

    def get_task_status(self, from_file=True):
        if from_file:
            assert self.file_mode == 'r', 'This requires file in read mode'
            self.refresh_from_file()
        return self.task_status

    def close_file(self):
        self.f.close()

class GenomeStats(tables.IsDescription):
    chromosome = tables.Int32Col()
    position = tables.Int32Col()
    value = tables.Float64Col()

class PhenotypeValue(tables.IsDescription):
    """
    Phenotype value wrapper
    """
    ecotype = tables.Int32Col()
    accession_name = tables.StringCol(16)
    mean_value = tables.Float32Col()
    std_dev = tables.Float32Col()
    comment = tables.StringCol(256)



class ResultRecordLM(tables.IsDescription):
    """
    Linear model, mixed models, etc.
    """
    chromosome = tables.Int32Col()
    position = tables.Int32Col()
    score = tables.Float64Col()
    maf = tables.Float32Col()
    mac = tables.Int32Col()
    genotype_var_perc = tables.Float32Col()
    beta0 = tables.Float32Col()
    beta1 = tables.Float32Col()
    correlation = tables.Float32Col()


class ResultRecordKW(tables.IsDescription):
    """
    Kruskal Wallis
    """
    chromosome = tables.Int32Col()
    position = tables.Int32Col()
    score = tables.Float64Col()
    maf = tables.Float32Col()
    mac = tables.Int32Col()
    statistic = tables.Float32Col()


class ResultRecordFT(tables.IsDescription):
    """
    Fisher's exact test
    """
    chromosome = tables.Int32Col()
    position = tables.Int32Col()
    score = tables.Float64Col()
    maf = tables.Float32Col()
    mac = tables.Int32Col()
    odds_ratio = tables.Float32Col()





_file_counter = {}

class GWASRecord():

    def __init__(self, hdf5_file_name):
        self.filename = hdf5_file_name

    def open(self, mode, **args):
        if self.is_file_exists(self.filename):
            self.h5file = tables.open_file(self.filename, mode, **args)
        else:
            self.h5file = tables.open_file(self.filename, "w", **args)

    def close(self):
        self.h5file.close()

    def _open(self, mode, **args):
        h5file = tables.open_file(self.filename, mode, **args)
        if self.filename not in _file_counter:
            _file_counter[self.filename] = 1
        else:
            _file_counter[self.filename] += 1
        return h5file

    def _close(self, h5file=None):
        if self.filename in _file_counter:
            _file_counter[self.filename] -= 1
            if _file_counter[self.filename] == 0:
                print "closing hdf5 file"
                if h5file is not None:
                    h5file.close()
                elif self.h5file is not None:
                    self.h5file.close()

    def is_file_exists(self, file_name=None):
        import os.path
        if file_name is None:
            file_name = self.filename
        return os.path.isfile(file_name)


    def init_file(self):
        #self.h5file = self._open(mode="w", title="Phenotype_results_file")
        g = self.h5file.create_group("/", 'phenotypes', 'Basic phenotype folder')
        #h5file.create_table(g, 'info', PhenotypeInfo, "Phenotyping information")
        self.h5file.flush()
        #self._close()

    def add_new_phenotype(self, phen_name, phenotype_values, ecotypes, accession_names=None, growth_conditions='',
                phenotype_scoring='', method_description='', measurement_scale='', is_binary=False):
        """
        Initializes the phenotype group for this phenotype and inserts it into the file object.
        """
        #Now parsing the phenotype file
        #self.h5file = self._open(mode="r+")
        self._init_phenotype_(phen_name, num_vals=len(phenotype_values), std_dev=sp.std(phenotype_values),
                growth_conditions=growth_conditions, phenotype_scoring=phenotype_scoring,
                method_description=method_description, measurement_scale=measurement_scale,
                is_binary=is_binary)

        self._add_phenotype_values_(phen_name, ecotypes, phenotype_values, 'Fullset', transformation='raw',
                accessions=accession_names, std_dev_values=None, value_comments=None)
        self.h5file.flush()
        #self._close(self.h5file)


    def _init_phenotype_(self, phen_name, num_vals=0.0, std_dev=0.0, growth_conditions='', phenotype_scoring='',
                method_description='', measurement_scale='', bs_herits=None, bs_avg_herits=None, bs_herit_pvals=None, is_binary=False):
        """
        Insert a new phenotype into the DB
        """
        group = self.h5file.create_group("/phenotypes", phen_name, 'Phenotype folder for ' + phen_name)
        #table = self.h5file.create_table(group, 'transformation_info', TransformationInfo, "Transformation information")
        #table = self.h5file.get_node('/phenotypes/info')
        #info = table.row
        group._v_attrs.name = phen_name
        group._v_attrs.num_vals = num_vals
        group._v_attrs.std_dev = std_dev
        group._v_attrs.growth_conditions = growth_conditions
        group._v_attrs.phenotype_scoring = phenotype_scoring
        group._v_attrs.method_description = method_description
        group._v_attrs.measurement_scale = measurement_scale
        group._v_attrs.std_dev = is_binary
        group._v_attrs.bs_herits = bs_herits
        group._v_attrs.bs_avg_herits = bs_avg_herits
        group._v_attrs.bs_herit_pvals = bs_herit_pvals
        group._v_attrs.current_dataset_id = 1
        fullset_group = self.h5file.create_group(group, 'Fullset', 'Fullset phenotype values')
        fullset_group._v_attrs.name = 'Fullset'
        fullset_group._v_attrs.id = 'Fullset';
        fullset_group._v_attrs.description = 'Fullset of all phenotypes'
        """info['name'] = phen_name
        info['num_values'] = num_vals
        info['std_dev'] = std_dev
        info['growth_conditions'] = growth_conditions
        info['phenotype_scoring'] = phenotype_scoring
        info['method_description'] = method_description
        info['measurement_scale'] = measurement_scale
        info['is_binary'] = is_binary
        info.append()
        table.flush()"""

    def _check_phenotype_exists_(self, phen_name):
        if "phenotypes" not in self.h5file.root:
            self.init_file()
        if phen_name in self.h5file.root.phenotypes:
            return True
        return False


    def add_phenotype_file(self, ecotype_ids, phen_file_name=None, file_object=None, transformation='raw',
                           transformation_description=None):
        """
        Adds phenotype values, to an existing phenotype, e.g. when applying different transformations.
        """
        retval = {'status':'OK', 'statustext':''}
        phed = pd.parse_phenotype_file(phen_file_name, file_object, with_db_ids=False)
        phed.convert_to_averages()
        phed.filter_ecotypes_2(ecotype_ids)
        
        #self.h5file = self._open(mode="r+")
        growth_conditions = None
        phenotype_scoring = None
        method_description = None
        measurement_scale = None
        is_binary = None
        try:
            for pid in phed.phen_dict:
                phen_vals = phed.get_values(pid)
                ecotypes = phed.get_ecotypes(pid)
                phen_name = phed.get_name(pid)
                try:
                    if self._check_phenotype_exists_(phen_name):
                        raise Exception("Phenotype %s already exists, delete it first to upload it again.<br>" % phen_name)
                    (bs_herits, bs_pids, bs_avg_herits, bs_herit_pvals) = phed.get_broad_sense_heritability()
                    self._init_phenotype_(phen_name, num_vals=len(phen_vals), std_dev=sp.std(phen_vals),
                          growth_conditions=growth_conditions, phenotype_scoring=phenotype_scoring,
                          method_description=method_description, measurement_scale=measurement_scale, bs_herit_pvals=bs_herit_pvals, bs_avg_herits=bs_avg_herits, bs_herits=bs_herits,
                          is_binary=is_binary)
                    self._add_phenotype_values_(phen_name, ecotypes, phen_vals, dataset='Fullset', transformation=transformation, transformation_description=transformation_description)
                except Exception, err:
                    retval['status'] = "ERROR"
                    retval['statustext'] += str(err)
        except Exception, err:
            raise(err)
        finally:
            test = ''
            #self._close()
        return retval

    def add_phenotype_values(self, phen_name, ecotypes, values, dataset='Fullset', transformation='raw', transformation_description=None,
                accessions=None, std_dev_values=None, value_comments=None, sp_pval=None):
        """
        Adds phenotype values, to an existing phenotype, e.g. when applying different transformations.
        """
        #self.h5file = self._open(mode="r+")

        #Calculating Shapiro-Wilks p-value
        if sp_pval is None:
            sp_pval = self.calculate_sp_pval(values)

        try:
            self._add_phenotype_values_(phen_name, ecotypes, values, dataset, transformation, transformation_description,
                                        accessions, std_dev_values, value_comments, sp_pval)
            self.h5file.flush()
        except Exception, err:
            raise(err)
        finally:
            test = ''
            #self._close()

    def delete_phenotype(self, phen_name):
        #self.h5file = self._open(mode="r+")
        try:
            try:
                self._delete_phenotype(phen_name)
            except Exception, err:
                raise err
        finally:
            self.h5file.flush()
            #self._close()

    def _delete_phenotype(self, phen_name):
        self.h5file.remove_node('/phenotypes' , phen_name, True)


    def delete_analysis(self, phen_name, dataset, transformation, analysis):
        #self.h5file = self._open(mode="r+")
        try:
            try:
                self._delete_analysis(phen_name, dataset, transformation, analysis)
            except Exception, err:
                raise err
        finally:
            self.h5file.flush()
            #self._close()

    def delete_result(self, phen_name, dataset, transformation, analysis, result_name):
        try:
            try:
                self._delete_result(phen_name, dataset, transformation, analysis, result_name)
            except Exception, err:
                raise err
        finally:
            self.h5file.flush()
            
    def delete_custom_genomestat(self, phen_name, dataset, stat):
        try:
            try:
                self._delete_custom_genomestat(phen_name, dataset, stat)
            except Exception, err:
                raise err
        finally:
            self.h5file.flush()
            #self._close()
    def _delete_custom_genomestat(self, phen_name, dataset, stat):
        self.h5file.remove_node('/phenotypes/%s/%s/stats' % (phen_name, dataset), stat)

    def add_genome_stats(self, phen_name, dataset, genomestat_name, genomestat_file):
        try:
            retval = self._add_genome_stats(phen_name, dataset, genomestat_name, genomestat_file)
            self.h5file.flush()
        except Exception, err:
            raise(err)
        finally:
            test = ''
        return retval
    
    def _add_genome_stats(self, phen_name, dataset, genomestat_name, genomestat_file):
        dataset_group = self.h5file.get_node('/phenotypes/%s/%s' % (phen_name, dataset))
        if not 'stats' in dataset_group:
            g = self.h5file.create_group(dataset_group, 'stats', 'Genomestats')
        else:
            g = self.h5file.get_node(dataset_group, 'stats')
        if genomestat_name in g:
            raise Exception('%s already exists' % genomestat_name)
        import csv 
        csvfile = open(genomestat_file, 'rb')
        dialect = csv.Sniffer().sniff(csvfile.read(1024))
        csvfile.seek(0)
        reader = csv.DictReader(csvfile, dialect=dialect)
        first_row = reader.next()
        for fieldname in ['chromosome', 'position', 'value']:
            if fieldname not in first_row:
                raise Exception('%s fieldname must appear in the csv' % fieldname)
        chr = int(first_row['chromosome'])
        position = int(first_row['position'])
        value = float(first_row['value'])
        csvfile.seek(0)
        reader.next()
        table = self.h5file.create_table(g, genomestat_name, GenomeStats, "Genomestats %s" % genomestat_name)
        value = table.row
        chr_regions = []
        chr = 1
        start_ix = 0
        stop_ix = 0
        for row in reader:
            if int(row['chromosome']) != chr:
                chr_regions.append((start_ix, stop_ix)) 
                start_ix = stop_ix
                chr = int(row['chromosome'])
            stop_ix = stop_ix + 1 
            value['chromosome'] = int(row['chromosome'])
            value['position'] = int(row['position'])
            value['value'] = float(row['value'])
            value.append()
        chr_regions.append((start_ix, stop_ix))
        table._v_attrs.chr_regions = chr_regions
        table.flush()
        return {'status':'OK', 'statustext':'', 'data':{'name':genomestat_name, 'isCustom':True, 'isStackable':False}}
            
            

    def _delete_analysis(self, phen_name, dataset, transformation, analysis):
        self.h5file.remove_node('/phenotypes/%s/%s/%s' % (phen_name, dataset, transformation), analysis, True)

    def _delete_result(self, phen_name, dataset, transformation, analysis, result_name):
        g = self.h5file.get_node('/phenotypes/%s/%s/%s/%s' % (phen_name, dataset, transformation, analysis))
        self.h5file.remove_node(g, result_name, True)
        quantiles_table_name = 'quantiles'
        if result_name != 'results':
            quantiles_table_name = 'quantiles_%s' % result_name
        log_quantiles_table_name = "log_%s" % quantiles_table_name
        if quantiles_table_name in g:
            self.h5file.remove_node(g, quantiles_table_name)
        if log_quantiles_table_name in g:
            self.h5file.remove_node(g, log_quantiles_table_name)
        ld_value_node = '%s_LD' % result_name 
        ld_snp_node = '%s_LD_SNPS' % result_name
        if ld_value_node in g:
            self.h5file.remove_node(g, ld_value_node)
        if ld_snp_node in g:
            self.h5file.remove_node(g, ld_snp_node)

    def delete_transformation(self, phen_name, dataset, transformation='raw'):
        """
        Deletes transformation from an existing phenotype.
        """
        #self.h5file = self._open(mode="r+")
        try:
            try:
                self._delete_transformation_(phen_name, dataset, transformation)
            except Exception, err:
                raise err
        finally:
            self.h5file.flush()
            #self._close()

    def delete_dataset(self, phen_name, dataset):
        try:
            try:
                self._delete_dataset_(phen_name, dataset)
            except Exception, err:
                raise err
        finally:
            self.h5file.flush()

    def _delete_transformation_(self, phen_name, dataset, transformation='raw'):
        self.h5file.remove_node('/phenotypes/%s/%s' % (phen_name, dataset), transformation, True)

    def _is_delete_transformation_allowed_(self, phen_name, transformation):
        for dataset_group in self.h5file.iter_nodes("/phenotypes/%s" % phen_name, 'Group'):
            if dataset_group._v_attrs.id != 'Fullset' and transformation in dataset_group:
                return False
        return True

    def _delete_dataset_(self, phen_name, dataset):
        self.h5file.remove_node("/phenotypes/%s" % phen_name, dataset, True)

    def _add_dataset_(self, phen_name):
        phen_group = self.h5file.get_node("/phenotypes", phen_name)
        phen_group._v_attrs.current_dataset_id = phen_group._v_attrs.current_dataset_id + 1
        dataset_id = 'S%s' % phen_group._v_attrs.current_dataset_id
        dataset_group = self.h5file.create_group(phen_group, str(dataset_id))
        dataset_group._v_attrs.id = dataset_id
        return dataset_group

    def _get_subset_phenotype_values(self, phen_name, ecotypes):
        import bisect
        table = self.h5file.get_node('/phenotypes/%s/Fullset/raw/values' % phen_name)
        ids = table.read(field='ecotype')
        indices = [bisect.bisect(ids, val) - 1 for val in ecotypes]
        phen_vals = table.read_coordinates(indices, field='mean_value')
        return phen_vals


    def _add_phenotype_values_(self, phen_name, ecotypes, values, dataset='Fullset', transformation='raw', transformation_description=None,
                accessions=None, std_dev_values=None, value_comments=None, sp_pval=None):
        """
        """

        phen_group = self.h5file.get_node('/phenotypes/%s/%s' % (phen_name, dataset))
        #table = self.h5file.get_node('/phenotypes/%s/transformation_info' % phen_name)
        #info = table.row
        #info['name'] = transformation
        #if transformation_description: info['description'] = transformation_description
        #info.append()
        #table.flush()

        trans_group = self.h5file.create_group(phen_group, transformation, 'Transformation: ' + transformation)
        trans_group._v_attrs.name = transformation
        trans_group._v_attrs.description = transformation_description
        trans_group._v_attrs.sp_pval = sp_pval
        table = self.h5file.create_table(trans_group, 'values', PhenotypeValue, "Phenotype values")
        acc_values = zip(ecotypes, values)
        acc_values.sort()
        value = table.row
        for i, (ei, v) in enumerate(acc_values):
            value['ecotype'] = ei
            ecotypes[i] = int(ei)
            value['mean_value'] = v
            if accessions: value['accession_name'] = accessions[i]
            if std_dev_values: value['std_dev'] = std_dev_values[i]
            if value_comments: value['comment'] = value_comments[i]
            value.append()
        phen_group._v_attrs.ecotypes = ecotypes
        table.flush()




    def get_phenotype_values(self, phen_name, dataset, transformation='raw'):
        import bisect
        """
        Returns the phenotype values
        """
        #self.h5file = self._open(mode="r")
        sp_pval = None
        table = self.h5file.get_node('/phenotypes/%s/%s/%s/values' % (phen_name, dataset, transformation))
        d = {'ecotype' : [], 'mean_value' : [], 'accession_name': [], 'std_dev': [], 'comment':[]}
        if 'sp_pval' in table._v_attrs:
            sp_pval = table._v_attrs.sp_pval
        for x in table.iterrows():
            for k in d:
                d[k].append(x[k])
        if sp_pval is None:
            sp_pval = self.calculate_sp_pval(d['mean_value'])
            table._v_attrs.sp_pval = sp_pval
        d['sp_pval'] = sp_pval
        #self._close()

        return d



    def get_phenotype_info(self, phen_name=None):
        """
        Returns the phenotype meta data in a dict.
        """
        dict_list = []
        #self.h5file = self._open(mode="r")
        try:
            if "phenotypes" not in self.h5file.root:
                return dict_list
            if not phen_name:
                for phenotype_table in self.h5file.iter_nodes("/phenotypes", 'Group'):
                    d = {'id':phenotype_table._v_attrs.name, 'name': phenotype_table._v_attrs.name, 'num_values': phenotype_table._v_attrs.num_vals, 'std_dev': phenotype_table._v_attrs.std_dev, 'growth_conditions': phenotype_table._v_attrs.growth_conditions,
                        'phenotype_scoring': phenotype_table._v_attrs.phenotype_scoring, 'method_description': phenotype_table._v_attrs.method_description, 'measurement_scale': phenotype_table._v_attrs.measurement_scale,
                        'is_binary': False}
                    d['datasets'] = self._get_phenotype_datasets_(d['name'])
                    dict_list.append(d)
            else:
                x = self.h5file.get_node("/phenotypes/%s" % phen_name)
                dict_list = [{'id':x._v_attrs.name, 'name': x._v_attrs.name, 'num_values': x._v_attrs.num_vals, 'std_dev': x._v_attrs.std_dev, 'growth_conditions': x._v_attrs.growth_conditions,
                            'phenotype_scoring': x._v_attrs.phenotype_scoring, 'method_description': x._v_attrs.method_description, 'measurement_scale': x._v_attrs.measurement_scale,
                            'is_binary': False, 'datasets':self._get_phenotype_datasets_(x._v_attrs.name)}]
        except Exception, err:
            raise(err)
        finally:
            test = ''
            #self._close()
        return dict_list



    def get_qq_plots(self, phen_name, dataset, transformation, analysis=None, result_name='results'):
            """
            Returns the qq-plot  for one or more analysis methods
            """
        #try:
            import matplotlib.pyplot as plt
            import cStringIO
            colors = {'emmax':'red', 'amm':'red', 'kw':'blue', 'lm':'green'}
            output = cStringIO.StringIO()
            quantile_ls = []
            quantile_ls_log = []
            quantile_labels = []
            line_colors = []
            if analysis is not None: # qqplot for a single GWAS result
                label = analysis
                quantiles_table_name = 'quantiles'
                if result_name != 'results':
                    label = "%s_%s" % (label, result_name)
                    quantiles_table_name = 'quantiles_%s' % result_name
                log_quantiles_table_name = "log_%s" % quantiles_table_name

                data_log = self.h5file.get_node('/phenotypes/%s/%s/%s/%s/%s' % (phen_name, dataset, transformation, analysis.lower(), log_quantiles_table_name))
                data = self.h5file.get_node('/phenotypes/%s/%s/%s/%s/%s' % (phen_name, dataset, transformation, analysis.lower(), quantiles_table_name))
                quantile_ls.append(data[:, 1])
                quantile_ls_log.append(data_log[:, 1])
                quantile_labels.append(label)
            else:
                group = self.h5file.get_node('/phenotypes/%s/%s/%s/' % (phen_name, dataset, transformation))
                for analysis in self.h5file.iter_nodes(group, classname='Group'):
                    for result in self.h5file.iter_nodes(analysis, classname='Table'):
                        result_name = result.name
                        quantiles_table_name = 'quantiles'
                        if result_name != 'results':
                            quantiles_table_name = 'quantiles_%s' % result_name
                        log_quantiles_table_name = "log_%s" % quantiles_table_name
                        data_log = self.h5file.get_node(analysis, log_quantiles_table_name)
                        data = self.h5file.get_node(analysis, quantiles_table_name)
                        quantile_ls.append(data[:, 1])
                        quantile_ls_log.append(data_log[:, 1])
                        if result_name == 'results':
                            line_colors.append(colors[result._v_attrs.name])
                        else:
                            line_colors.append(colors[analysis._v_attrs.type])
                        if result_name == 'results':
                            quantile_labels.append(result._v_attrs.name)
                        else:
                            quantile_labels.append(result_name)
            image_data = None
            f = plt.figure(figsize=(10.8, 5))
            ax = f.add_axes([0.07, 0.09, 0.4, 0.86])
            agr.simple_log_qqplot(quantile_ls_log, ax=ax, format='png', quantile_labels=quantile_labels, line_colors=line_colors, max_val=6)
            ax2 = f.add_axes([0.57, 0.09, 0.4, 0.86])
            agr.simple_qqplot(quantile_ls, ax=ax2, format='png', quantile_labels=quantile_labels, line_colors=line_colors)
            f.savefig(output)
            output.seek(0)
            image_data = output.readlines()
            output.close()
        #except Exception, err:
           # raise(err)
        #finally:
            
            return image_data



    def get_qq_plot_data(self, phen_name, dataset, transformation, analysis=None, results_name='results'):
        """
        Returns the qq-plot data for one or more analysis methods
        """
        try:

            if analysis is not None: # qqplot for a single GWAS result
                dataset = self.h5file.get_node('/phenotypes/%s/%s/%s/%s/results' % (phen_name, dataset, transformation, analysis))
            else:
                group = self.h5file.get_node('/phenotypes/%s/%s/%s/' % (phen_name, dataset, transformation))
                # loop through all results and retreive qqplot data
            result = {'expected':[], 'observed':{'RESULT1':[], 'RESULT2':[]}}
            """
            result must have following structure:
            result = {'expected':[int,int,int,int],'observered':{'EMMAX':[int,int,int,int],'KW':[int,int,int],...}}
            """

        except Exception, err:
            raise(err)
        finally:
            test = ''
        return result



    def _get_phenotype_datasets_(self, phen_name):
        dict_list = []
        #self.h5file.flush()
        for dataset in self.h5file.iter_nodes("/phenotypes/%s" % phen_name, 'Group'):
            s = {'id':str(dataset._v_attrs.id), 'name':dataset._v_attrs.name, 'accessionCount':len(dataset._v_attrs.ecotypes), 'phenotype':phen_name, 'accession_ids':dataset._v_attrs.ecotypes, 'description':dataset._v_attrs.description, 'transformations':self._get_dataset_transformations_(phen_name, dataset._v_attrs.id)}
            dict_list.append(s)
        return dict_list

    def _get_dataset_transformations_(self, phen_name, dataset_name):
        dict_list = []
        #table = self.h5file.get_node('/phenotypes/%s/transformation_info' % phen_name)
        for x in self.h5file.iter_nodes('/phenotypes/%s/%s' % (phen_name, dataset_name), 'Group'):
            if x._v_name == 'stats':
                continue
            d = {'id':x._v_attrs.name, 'name': x._v_attrs.name, 'description': x._v_attrs.description}
            d['phenotype'] = phen_name
            d['dataset'] = str(dataset_name)
            d['analysis_methods'] = self._get_analysis_methods_(phen_name, dataset_name, d['name'])
            dict_list.append(d)
        return dict_list



    def get_dataset_transformations(self, phen_name, dataset):
        """
        Returns the phenotype values
        """
        #self.h5file = self._open(mode="r")
        d = self._get_dataset_transformations_(phen_name, dataset)
        #self._close()
        return d

    def get_dataset_accession_ids(self, phen_name, dataset=None):

        accession_ids = {}
        if dataset is not None:
            accession_ids[dataset] = self.h5file.get_node("/phenotypes/%s/%s" % (phen_name, dataset))._v_attrs.ecotypes
        else:
            for x in self.h5file.iter_nodes('/phenotypes/%s' % (phen_name), 'Group'):
                accession_ids[x._v_attrs.id] = x._v_attrs.ecotypes
        return accession_ids

    def _get_analysis_methods_(self, phen_name, dataset, transformation):
        import math
        dict_list = []
        stats_to_return = ['step', 'chr', 'pos', 'bic', 'ebic', 'pseudo_heritability', 'max_cof_pval', 'remain_perc_gen_var', 'remain_perc_err_var', 'perc_var_expl']
        try:
            #table = self.h5file.get_node('/phenotypes/%s/%s/result_info' % (phen_name, transformation))
            for x in self.h5file.iter_nodes('/phenotypes/%s/%s/%s' % (phen_name, dataset, transformation), 'Group'):
                for res in x._f_iter_nodes('Table'):
                    #if res.name == 'results':
                    #    d = {'name': x._v_attrs.name, 'comment': x._v_attrs.comment,'snps':[]}
                    #else:
                    stats = []
                    if 'cofactors' in res._v_attrs:
                        for cofactor in res._v_attrs.cofactors:
                            stat = {}
                            for key in stats_to_return:
                                if key in cofactor:
                                    val = cofactor[key]
                                    if key in ['bic', 'ebic', 'pseudo_heritability', 'remain_perc_gen_var', 'remain_perc_err_var', 'perc_var_expl']:
                                        val = round(val, 2)
                                    elif key == 'max_cof_pval':
                                        if val > 0.0:
                                            val = round(-math.log10(val), 2)
                                    stat[key] = val
                            stats.append(stat)
                    d = {'id':res._v_attrs.name, 'name':res._v_attrs.name, 'resultName':res.name, 'comment':res._v_attrs.comment, 'type':x._v_attrs.type, 'cofactors':stats}
                    med_pval = None
                    ks_stat = None
                    ks_pval = None
                    if 'med_pval' in res._v_attrs:
                        med_pval = round(res._v_attrs.med_pval,4)
                    if 'ks_stats' in res._v_attrs:
                        ks_stat = round(res._v_attrs.ks_stats['D'],4)
                        ks_pval = res._v_attrs.ks_stats['p_val']
                    if ks_pval is not None and ks_pval != 0:
                        ks_pval = round(-math.log10(ks_pval),4)
                    d['med_pval'] = med_pval
                    d['ks_stat'] = ks_stat
                    d['ks_pval'] = ks_pval
                    d['phenotype'] = phen_name
                    d['dataset'] = str(dataset)
                    d['transformation'] = transformation
                    dict_list.append(d)
        except Exception, err_str:
            print "No results found:", err_str
        return dict_list



    def get_analysis_methods(self, phen_name, transformation):
        """
        Returns the phenotype values
        """
        #self.h5file = self._open(mode="r")
        d = self._get_analysis_methods_(phen_name, transformation)
        #self._close()
        return d



    def add_results(self, phen_name, dataset, analysis_method, name, chromosomes, positions, scores, mafs, macs,
            result_name='results', cofactors=[], analysis_comment='', transformation='raw', quantiles_dict=None, count=None,
            bh_thres=None,med_pval=None,ks_stats = None, **kwargs):
        """
        Add a result to the hdf5 file.
        """
        #h5file = self._open(mode="r+")
        try:

            trans_group = self.h5file.get_node('/phenotypes/%s/%s/%s' % (phen_name, dataset, transformation))
            if analysis_method not in trans_group:
                analysis_group = self.h5file.create_group(trans_group, analysis_method, 'Analysis method: ' + analysis_method)
                analysis_group._v_attrs.type = analysis_method
            else:
                analysis_group = trans_group._f_get_child(analysis_method)

            if analysis_method in ['amm', 'lm']:
                table = self.h5file.create_table(analysis_group, result_name, ResultRecordLM, "Regression result")
            elif analysis_method == 'kw':
                table = self.h5file.create_table(analysis_group, result_name, ResultRecordKW, "Regression result")
            else:
                raise Exception('Not implemented for analysis method %s' % analysis_method)

            #print table
            table._v_attrs.data_version = 1.0
            table._v_attrs.max_score = max(scores)
            table._v_attrs.comment = analysis_comment
            table._v_attrs.name = name
            table._v_attrs.cofactors = cofactors
            result = table.row
            
            
            #Setting the p-value threshold (not negative log transformed)
            if not bh_thres:
                table._v_attrs.pval_threshold = 0.05 / len(scores)
            else:
                table._v_attrs.pval_threshold = bh_thres


            if med_pval:
                table._v_attrs.med_pval = med_pval
            
            if ks_stats is not None:
                table._v_attrs.ks_stats = ks_stats
            
            
            quantiles_table_name = 'quantiles'
            if result_name != 'results':
                quantiles_table_name = 'quantiles_%s' % result_name
            log_quantiles_table_name = "log_%s" % quantiles_table_name

            if quantiles_dict:
                #Create a quantile array
                if quantiles_table_name not in analysis_group:
                    quantiles_array = self.h5file.create_carray(analysis_group, quantiles_table_name, tables.Float32Atom(), shape=(len(quantiles_dict['quantiles']), 2))
                else:
                    quantiles_array = self.h5file.get_node(analysis_group, quantiles_table_name)
                quantiles_array[:, 0] = quantiles_dict['exp_quantiles']
                quantiles_array[:, 1] = quantiles_dict['quantiles']

                if log_quantiles_table_name not in analysis_group:
                    log_quantiles_array = self.h5file.create_carray(analysis_group, log_quantiles_table_name, tables.Float32Atom(),
                                                               shape=(len(quantiles_dict['log_quantiles']), 2))
                else:
                    log_quantiles_array = self.h5file.get_node(analysis_group, log_quantiles_table_name)
                log_quantiles_array[:, 0] = quantiles_dict['exp_log_quantiles']
                log_quantiles_array[:, 1] = quantiles_dict['log_quantiles']

            t_list = zip(chromosomes, positions, scores, mafs, macs)
            t_list.sort(key=lambda k:(k[0], k[2] * -1))
            for i, cpsmm in enumerate(t_list):
                (result['chromosome'], result['position'], result['score'], result['maf'], result['mac']) = cpsmm
                if analysis_method == 'kw':
                    result['statistic'] = kwargs['statistics'][i]
                else: #EMMAX or LM
                    if kwargs['beta0'] is not None:
                        result['beta0'] = kwargs['beta0'][i]
                    if kwargs['beta1'] is not None:
                        result['beta1'] = kwargs['beta1'][i]
                    if 'correlation' in kwargs:

                        result['correlation'] = kwargs['correlation'][i]
                    result['genotype_var_perc'] = kwargs['genotype_var_perc'][i]
                result.append()
            table.flush()
            table.cols.chromosome.create_index()
            table.cols.score.create_index()
            table.cols.mac.create_index()
            table.flush()
            if count is not None:
                analysis_group._v_attrs.count = count
            self.h5file.flush()
        except Exception, err:
            raise err
        finally:
            test = ''
#        table.cols.chromosome.create_index()
#        table.cols.score.create_index()
#        table.cols.mac.create_index()
            #self._close(h5file)




    def get_results(self, phen_name, dataset, transformation, analysis_method, result_name, min_mac=0, max_pval=1.0):
        """
        Return results..
        """
        d = {'chromosome': [], 'position': [], 'score': [], 'maf': [], 'mac': []}
        if analysis_method == 'kw':
            d['statistic'] = []
        else:
            d['beta0'] = []
            d['beta1'] = []
            d['correlation'] = []
            d['genotype_var_perc'] = []

        #h5file = self.open(mode="r")
        try:
            info_group = self.h5file.get_node('/phenotypes/%s/%s/%s/%s' % (phen_name, dataset, transformation, analysis_method))
            table = info_group._f_get_child(result_name)
            #for x in table.where('(score<=%f) & (mac>=%d)' % (max_pval, min_mac)):
            for x in table.iterrows():
                if x['score'] <= max_pval and x['mac'] >= min_mac:
                    for k in d:
                        d[k].append(x[k])
        except Exception, err:
            raise(err)
        finally:
            test = ''
            #self._close(h5file)
        return d
    
    def get_genomestats_list(self, phen_name, dataset):
        data = []
        try:
            dataset_group = self.h5file.get_node('/phenotypes/%s/%s' % (phen_name, dataset))
            if 'stats' in dataset_group:
                stats_group = self.h5file.get_node(dataset_group, 'stats')
                for stats_table in self.h5file.iter_nodes(stats_group, classname='Table'):
                    data.append({'name':stats_table.name, 'isStackable':False, 'isCustom':True})
            retval = {'status':'OK', 'stats': data}   
        except Exception, err:
            retval = {'status':'ERROR', 'statustext':str(err)}
        return retval
    
    def get_genomestats_data(self, phen_name, dataset, stat, chr):
        data = []
        try:
            chr_numbers = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5']
            chr_num = chr_numbers.index(chr.lower())
            stats_table = self.h5file.get_node('/phenotypes/%s/%s/stats/%s/' % (phen_name, dataset, stat))
            chr_region = stats_table._v_attrs['chr_regions'][chr_num]
            data = stats_table[chr_region[0]:chr_region[1]]
        except Exception, err:
            raise(err)
        return data

    def get_results_for_csv(self, phen_name, dataset, transformation, analysis_method, result_name, min_mac=0):
        """
        Return results..
        """
        header = ['chromosome', 'position', 'score', 'maf', 'mac']
        if analysis_method == 'kw':
            header.append('statistic')
        else:
            header.extend(['beta0', 'beta1', 'correlation', 'genotype_var_perc'])
        result = [header[:]]
        #h5file = self.open(mode="r")
        try:
            info_group = self.h5file.get_node('/phenotypes/%s/%s/%s/%s' % (phen_name, dataset, transformation, analysis_method))
            table = info_group._f_get_child(result_name)
            #for x in table.where('(score<=%f) & (mac>=%d)' % (max_pval, min_mac)):
            for x in table.iterrows():
                row = []
                if x['mac'] >= min_mac:
                    for col in header:
                        row.append(x[col])
                    result.append(row)
        except Exception, err:
            raise(err)
        finally:
            test = ''
            #self._close(h5file)
        return result

    def get_results_by_chromosome(self, phen_name, dataset, analysis_method, result_name, transformation='raw', min_mac=15, min_score=0.0, \
                top_fraction=0.05, chromosomes=[1, 2, 3, 4, 5], log_transform=True):
        """
        Return results..
        """
        cd = {}
        #h5file = self._open(mode="r")
        try:
            info_group = self.h5file.get_node('/phenotypes/%s/%s/%s/%s' % (phen_name, dataset, transformation, analysis_method))
            table = info_group._f_get_child(result_name)

            max_score = table._v_attrs.max_score

            d_keys = ['score', 'position', 'maf', 'mac']
            if analysis_method == 'kw':
                d_keys.append('statistic')
            else:
                d_keys.extend(['beta0', 'beta1', 'correlation', 'genotype_var_perc'])

            for chromosome in chromosomes:
                sort_list = []
                sort_list = table.read_where('(chromosome==%d) &(score>=%f) & (mac>=%d)' % (chromosome, min_score, min_mac))
                if 'data_version' not in table._v_attrs:
                    sort_list.sort(order=['score'])
                    sort_list = sort_list[::-1]
                sort_list = sort_list[:int(top_fraction * len(sort_list))]
                sort_list['score'] = sp.around(sort_list['score'], decimals=2)
                d = {}
                for k in d_keys:
                    d[k] = sort_list[k].tolist()
                cd[chromosome] = d
            cd['chromosome_ends'] = chromosome_ends
            cd['max_score'] = max_score
            pvals = None
            
            if 'pval_threshold' not in table._v_attrs:
                pvals = 10 ** -table[:]['score']
                table._v_attrs.pval_threshold = mtcorr.get_bhy_thres(pvals, fdr_thres=0.05)['thes_pval']
            cd['pval_threshold'] = -math.log10(table._v_attrs.pval_threshold)
            if 'med_pval' not in table._v_attrs:
                if pvals is None:
                    pvals = 10 ** -table[:]['score']
                table._v_attrs.med_pval = agr.calc_median(pvals)
            cd['med_pval'] = round(table._v_attrs.med_pval,4)
            
            if 'ks_stats' not in table._v_attrs:
                if pvals is None:
                    pvals = 10 ** -table[:]['score']
                table._v_attrs.ks_stats = agr.calc_ks_stats(pvals) 
            cd['ks_stat'] =  round(table._v_attrs.ks_stats['D'],4)
            ks_pval = table._v_attrs.ks_stats['p_val']
            if  ks_pval != 0:
                ks_pval = -math.log10(ks_pval)
            cd['ks_pval'] =  round(ks_pval,4)
        

        #Calculate the Kolmogorov-Smirnov statistic
        
            
        except Exception, err:
            raise(err)
        finally:
            test = ''
            #self._close(h5file)
        return cd

    def calculate_ld(self, phen_name, dataset, transformation, analysis_method, result_name, call_method_id=82, progress_file_writer=None):
        try:
            analysis_group = self.h5file.get_node('/phenotypes/%s/%s/%s/%s' % (phen_name, dataset, transformation, analysis_method))
            if '%s_LD' % result_name in analysis_group:
                return {'status':'OK'}
            if progress_file_writer is not None:
                progress_file_writer.update_progress_bar(progress=0.95, task_status='Calculating LD')
            accession_ids = map(str, self.get_dataset_accession_ids(phen_name, dataset)[dataset])
            sd = dp.load_snps_call_method(call_method_id=call_method_id, data_format='binary', min_mac=5) #Load SNPs data
            sd.filter_accessions(accession_ids)
            res = self.get_results_by_chromosome(phen_name, dataset, analysis_method, result_name, transformation)
            snps_to_keep = []
            snps_to_store = []
            snps_to_store_regions = []
            for chr in  [1, 2, 3, 4, 5]:
                snps_to_keep_per_chr = res[chr]['position']
                snps_to_keep_per_chr.sort()
                snps_to_store_regions.append((len(snps_to_store), len(snps_to_store) + len(snps_to_keep_per_chr))) 
                snps_to_store.extend(snps_to_keep_per_chr[:])
                snps_to_keep.append(snps_to_keep_per_chr)
            snps_to_store = sp.array(snps_to_store)
            sd.filter_snps(snps_to_keep)
            snps = sp.array(sd.get_snps())
            if progress_file_writer is not None:
                progress_file_writer.update_progress_bar(progress=0.97, task_status='Calculating pairwise r2 values')
            r2_values = self._calculate_ld(snps)
            ldSNPsArray = self.h5file.create_carray(analysis_group, '%s_LD_SNPS' % result_name, tables.UInt32Atom(), snps_to_store.shape)
            ldSNPsArray[:] = snps_to_store
            ldSNPsArray._v_attrs.chr_regions = snps_to_store_regions
            ldArray = self.h5file.create_vlarray(analysis_group, '%s_LD' % result_name, tables.Float32Atom(shape=()), "ld values")
            for i in range(0, len(r2_values)):
                ldArray.append(r2_values[i][:i + 1])
            self.h5file.flush()
            if progress_file_writer is not None:
                progress_file_writer.update_progress_bar(progress=0.99, task_status='Finished calculating LD')
        except Exception, err:
            raise(err)
        finally:
            test = ''
            #self._close(h5file)
        return {'status':'OK'}
    
    def get_exact_ld(self, phen_name, dataset, transformation, analysis_method, result_name, chromosome, position, call_method_id=82):
        num_snps = 250
        accession_ids = map(str, self.get_dataset_accession_ids(phen_name, dataset)[dataset])
        sd = dp.load_snps_call_method(call_method_id=call_method_id, data_format='binary', min_mac=5)
        sd.filter_accessions(accession_ids)
        #snps_in_chr = sd.snpsDataList[int(chromosome) - 1]
        
        snps_to_keep = []
        for chr in  [1, 2, 3, 4, 5]:
            snps_to_keep_per_chr = []
            snps_in_chr = sd.snpsDataList[chr - 1]
            if chr == int(chromosome):
                pos_ix = snps_in_chr.positions.index(int(position))
                start_ix = 0 if pos_ix - num_snps < 0 else pos_ix - num_snps 
                end_ix = len(snps_in_chr.positions) if pos_ix + num_snps > len(snps_in_chr.positions) else pos_ix + num_snps
                snps_to_keep_per_chr = [ix for ix in range(start_ix, end_ix)]
            snps_in_chr.filter_snp_indices(snps_to_keep_per_chr)
        snps = sp.array(sd.get_snps())
        r2_values = self._calculate_ld(snps)
        ldArray = []
        snps = sd.get_positions()
        for i in range(0, len(r2_values)):
            ldArray.append(r2_values[i][:i + 1].tolist())
            #self._close(h5file)
        retval = {'snps':snps, 'r2':ldArray, 'start':snps[0], 'end':snps[-1]}
        return retval

    def _calculate_ld(self, genotype):
        genotype_t = sp.transpose(genotype)
        genotype_stand = sp.transpose((genotype_t - sp.mean(genotype, 1)) / sp.std(genotype, 1))
        r2_values = sp.dot(genotype_stand, sp.transpose(genotype_stand))
        r2_values *= (1.0 / genotype.shape[1])
        r2_values **= 2
        return r2_values
    
    
        
    def get_ld_for_snp(self, phen_name, dataset, transformation, analysis_method, result_name, chromosome, position):
        grp = self.h5file.get_node('/phenotypes/%s/%s/%s/%s' % (phen_name, dataset, transformation, analysis_method))
        if '%s_LD' % result_name not in grp:
            self.calculate_ld(phen_name, dataset, transformation, analysis_method, result_name)
        snp_array = self.h5file.get_node(grp, '%s_LD_SNPS' % result_name)
        snp_array_regions = snp_array._v_attrs.chr_regions
        snp_array_region = snp_array_regions[int(chromosome) - 1]
        snps = snp_array[:]
        snp_ix = sp.where(snps == int(position))[0][0]
        ld_array = self.h5file.get_node('/phenotypes/%s/%s/%s/%s/%s_LD' % (phen_name, dataset, transformation, analysis_method, result_name))
        data = []
        for chr in range(5):
            chr_data = {}
            chr_region = snp_array_regions[chr]
            snps_in_chr = snps[chr_region[0]:chr_region[1]]
            if chr_region[0] <= snp_ix:
                stop_ix = chr_region[1] if chr_region[1] <= snp_ix else snp_ix 
                r2_values_in_chr = ld_array[snp_ix][chr_region[0]:stop_ix].tolist()
                if chr_region[1] > snp_ix:
                    for r2_val_row in ld_array[snp_ix:chr_region[1]]:
                        r2_values_in_chr.append(r2_val_row[snp_ix].tolist())
            else:
                r2_values_in_chr = []
                for r2_val_row in ld_array[chr_region[0]:chr_region[1]]:
                        r2_values_in_chr.append(r2_val_row[snp_ix].tolist())
            chr_data['snps'] = snps_in_chr.tolist()
            chr_data['r2'] = r2_values_in_chr
            data.append(chr_data)
        retval = {'data':data}
        return retval  
        
    def get_ld_for_region(self, phen_name, dataset, transformation, analysis_method, result_name, chromosome, start, end):
        import decimal
        grp = self.h5file.get_node('/phenotypes/%s/%s/%s/%s' % (phen_name, dataset, transformation, analysis_method))
        if '%s_LD' % result_name not in grp:
            self.calculate_ld(phen_name, dataset, transformation, analysis_method, result_name)
        snp_array = self.h5file.get_node(grp, '%s_LD_SNPS' % result_name)
        snp_array_regions = snp_array._v_attrs.chr_regions
        snp_array_region = snp_array_regions[int(chromosome) - 1]
        snps = snp_array[snp_array_region[0]:snp_array_region[1]]
        start_ix = int(sp.where(snps == int(start))[0])
        end_ix = int(sp.where(snps == int(end))[0])
        snps_to_return = snps[start_ix:end_ix]
        ld_array = self.h5file.get_node('/phenotypes/%s/%s/%s/%s/%s_LD' % (phen_name, dataset, transformation, analysis_method, result_name))
        r2_values = ld_array[snp_array_region[0]:snp_array_region[1]]
        r2_values = r2_values[start_ix:end_ix]
        r2_values_to_return = [r2val[snp_array_region[0] + start_ix:].tolist() for r2val in r2_values]
        
        """saves a lot of space (1/3 because only 4 digits but takes 3 times longer to process"""
        #r2_values_to_return = [[decimal.Decimal('%.4f' % r2valConv) for r2valConv in r2val[snp_array_region[0] + start_ix:].tolist()] for r2val in r2_values]
        retval = {'snps':snps_to_return.tolist(), 'r2':r2_values_to_return}
        return retval
        

    def get_phenotype_bins(self, phen_name, transformation='raw', phen_vals=None, bin_number=20):
        if phen_vals is None:
            phen_data = self.get_phenotype_values(phen_name, transformation)
            phen_vals = phen_data['mean_value']
        return self._get_phenotype_bins(phen_vals, bin_number)

    def _get_phenotype_bins(self, phen_vals, bin_number=20):
        """
        Returns the 
        """
        #Get the phenotype

        min_phen_val = min(phen_vals)
        max_phen_val = max(phen_vals)
        chunk_size = ((max_phen_val - min_phen_val) / bin_number) * (1 + 1e-5 / bin_number)
        pa = (sp.array(phen_vals) - min_phen_val) / chunk_size
        bin_counts = sp.bincount(sp.array(pa , dtype='int32'))
        keys = ['x_axis', 'frequency']
        bin_list = []
        for i, bin_count in enumerate(bin_counts.tolist()):
            x1 = min_phen_val + i * chunk_size
            x2 = x1 + chunk_size
            bin_list.append({'x_axis':'%0.2f-%0.2f' % (x1, x2), 'frequency':bin_count})

        return bin_list


    def preview_transform_phenotype(self, phen_name, dataset, transformation, original_transformation='raw'):
        new_phen_vals, sp_pval = self.transform_phenotype(phen_name, dataset, transformation, original_transformation)
        return self._get_phenotype_bins(new_phen_vals), sp_pval


    def transform_phenotype(self, phen_name, dataset, transformation, original_transformation='raw', store=False):
        """
        Apply a transformation to a phenotype.
        """
        from scipy import stats
        """if store:
            if dataset != 'Fullset':
                phen_group = self.h5file.get_node('/phenotypes/%s/%s' % (phen_name, dataset))
                trans_group = self.h5file.create_group(phen_group, transformation, 'Transformation: ' + transformation)
                trans_group._v_attrs.name = transformation
                trans_group._v_attrs.description = ''
                self.h5file.flush()
            phen_data = self.get_phenotype_values(phen_name, "Fullset", original_transformation)
        else:"""
        phen_data = self.get_phenotype_values(phen_name, dataset, original_transformation)
        phen_vals = sp.array(phen_data['mean_value'])
        if transformation == 'raw':
            return
        elif transformation == 'log':
            new_phen_vals = sp.log(phen_vals - min(phen_vals) + sp.var(phen_vals) * 0.1)
        elif transformation == 'sqrt':
            new_phen_vals = sp.sqrt(phen_vals - min(phen_vals) + sp.var(phen_vals) * 0.1)
        elif transformation == 'box_cox':
            vals = phen_vals - min(phen_vals) + sp.var(phen_vals) * 0.1
            sw_pvals = []
            lambdas = sp.arange(-2.0, 2.1, 0.1)
            for l in lambdas:
                if l == 0:
                    vs = sp.log(vals)
                else:
                    vs = ((vals ** l) - 1) / l
                r = stats.shapiro(vs)
                if sp.isfinite(r[0]):
                    pval = r[1]
                else:
                    pval = 0.0
                sw_pvals.append(pval)
            print sw_pvals
            i = sp.argmax(sw_pvals)
            l = lambdas[i]
            if l == 0:
                new_phen_vals = sp.log(vals)
            else:
                new_phen_vals = ((vals ** l) - 1) / l

        #Calculating Shapiro-Wilks p-value
        sp_pval = self.calculate_sp_pval(new_phen_vals)
        if store:
            self.add_phenotype_values(phen_name, phen_data['ecotype'], new_phen_vals, dataset, transformation, sp_pval=sp_pval)
        return new_phen_vals, sp_pval

    def calculate_sp_pval(self, phen_vals):
        r = stats.shapiro(phen_vals)
        if sp.isfinite(r[0]):
            sp_pval = r[1]
        else:
            sp_pval = 0.0
        return sp_pval

    def dec_run_count(self):
        if 'run_count' not in self.h5file.root.phenotypes._v_attrs:
            self.h5file.root.phenotypes._v_attrs.run_count = 0
            return
        current_run_count = self.h5file.root.phenotypes._v_attrs.run_count
        if current_run_count > 0:
            self.h5file.root.phenotypes._v_attrs.run_count = self.h5file.root.phenotypes._v_attrs.run_count - 1

    def inc_run_count(self):
        if 'run_count' not in self.h5file.root.phenotypes._v_attrs:
            self.h5file.root.phenotypes._v_attrs.run_count = 0 
        self.h5file.root.phenotypes._v_attrs.run_count = self.h5file.root.phenotypes._v_attrs.run_count + 1
        self.h5file.root.phenotypes._v_attrs.run_date = time.time() + 10 * 60  

    def check_run_count(self):
        if 'run_count'  in self.h5file.root.phenotypes._v_attrs:
            if self.h5file.root.phenotypes._v_attrs.run_count > 0:
                if 'run_date' in self.h5file.root.phenotypes._v_attrs and self.h5file.root.phenotypes._v_attrs.run_date > time.time():
                    return 'GWAS analysis is already running'
                else:
                    self.h5file.root.phenotypes._v_attrs.run_count = 0

    def save_dataset(self, dataset):
        try:
            retval = {'status':'OK', 'statustext':''}

            if dataset['id'] is None or dataset['id'] == '' or dataset['id'] == '_NEW_SUBSET_':
                dataset_group = self._add_dataset_(dataset['phenotype'])
            else:
                dataset_group = self.h5file.get_node("/phenotypes/%s/%s" % (dataset['phenotype'], dataset['id']))
                self.h5file.remove_node(dataset_group, 'raw', recursive=True)
            dataset_group._v_attrs.name = dataset['name']
            dataset_group._v_attrs.description = dataset['description']
            ecotypes = dataset['accession_ids']
            ecotypes = sorted(ecotypes)
            dataset_group._v_attrs.ecotypes = ecotypes
            dataset_group.title = dataset['description']
            phen_vals = self._get_subset_phenotype_values(dataset['phenotype'], ecotypes)
            self._add_phenotype_values_(dataset['phenotype'], ecotypes, phen_vals, dataset_group._v_attrs.id, 'raw')
            self.h5file.flush()
            retval['statustext'] = str(dataset_group._v_attrs.id)
        except Exception, err:
            retval['status'] = 'ERROR'
            retval['statustext'] = str(err)
        return retval


    def exists_transformation(self, phen_name, dataset, transformation):
        transformations = self.get_dataset_transformations(phen_name, dataset)
        for trans in transformations:
            if trans['name'] == transformation:
                return True
        return False




    def perform_gwas(self, phen_name, dataset, transformation='raw', analysis_method='kw', call_method_id=82,
                     kinship_method='ibs', progress_file_writer=None):

        """
        Performs GWAS and updates the datastructure.
        """

        import bisect
        import gwa
        import kinship
        self.inc_run_count()
        step_wise = False
        if analysis_method not in ['lm', 'amm', 'kw']:
            raise Exception('analysis method %s not supported' % analysis_method)

        progress_file_writer.update_progress_bar(progress=0.0, task_status='Loading phenotype data')
        phen_dict = self.get_phenotype_values(phen_name, dataset, transformation) #Load phenotype
        phend = pd.phenotype_data({1:{'values':phen_dict['mean_value'], 'ecotypes':map(str, phen_dict['ecotype']), 'name':phen_name}})
        phend.convert_to_averages()
        progress_file_writer.update_progress_bar(task_status='Loading genotype data')
        sd = dp.load_snps_call_method(call_method_id=call_method_id, data_format='binary', min_mac=5) #Load SNPs data
        #Preprocessing data.
        progress_file_writer.update_progress_bar(step=0.05, task_status='Coordinating genotype and phenotype data')
        d = sd.coordinate_w_phenotype_data(phend, 1)
        phen_vals = phend.get_values(1)
        snps = sd.getSnps()
        positions = sd.getPositions()
        chromosomes = []
        for i, (s, c) in enumerate(itertools.izip(sd.snpsDataList, sd.chromosomes)):
            chromosomes.extend([c] * len(s.snps))
        maf_dict = sd.get_mafs()


#        kwargs = {}
        if analysis_method == 'amm':
            progress_file_writer.update_progress_bar(progress=0.15, task_status='Retrieving the kinship matrix')
            k = kinship.get_kinship(call_method_id=call_method_id, n_removed_snps=d['n_filtered_snps'], remain_accessions=sd.accessions)
            d = lm.emmax_step(phen_vals, sd, k, [], progress_file_writer=progress_file_writer, emma_num=200)
            progress_file_writer.update_progress_bar(progress=0.95, task_status='Processing and saving results')
            res = d['res']
            stats_dict = d['stats']
        elif analysis_method == 'lm':
            progress_file_writer.update_progress_bar(progress=0.3, task_status='Performing linear regression')
            d = lm.lin_reg_step(phen_vals, sd, [], progress_file_writer=progress_file_writer)
            progress_file_writer.update_progress_bar(progress=0.95, task_status='Processing and saving results')
            res = d['res']
            stats_dict = d['stats']
        elif analysis_method == 'kw':
            progress_file_writer.update_progress_bar(progress=0.7, task_status='Performing KW')
            res = util.kruskal_wallis(snps, phen_vals)
            progress_file_writer.update_progress_bar(progress=0.95, task_status='Processing and saving results')
            #scores = map(lambda x:-math.log10(x), kw_res['ps'])
        else:
            raise Exception('analysis method %s not supported' % analysis_method)

        #Calculate Benjamini-Hochberg threshold
        bh_thres_d = mtcorr.get_bhy_thres(res['ps'], fdr_thres=0.05)

        #Calculate Median p-value
        med_pval = agr.calc_median(res['ps'])

        #Calculate the Kolmogorov-Smirnov statistic
        ks_res = agr.calc_ks_stats(res['ps'])

        #Calculate quantiles
        num_dots = 1000
        max_val = 6
        quantiles = agr.get_quantiles(res['ps'], num_dots=num_dots)
        exp_quantiles = agr._getExpectedPvalueQuantiles_(num_dots)
        log_quantiles = agr.get_log_quantiles(res['ps'], num_dots=num_dots, max_val=max_val)
        exp_log_quantiles = sp.arange(1, num_dots + 1, dtype='single') / (num_dots + 1) * max_val

        quantiles_dict = {'quantiles':quantiles, 'exp_quantiles':exp_quantiles,
			'log_quantiles':log_quantiles, 'exp_log_quantiles':exp_log_quantiles}

       

        scores = map(lambda x:-math.log10(x), res['ps'])
        if analysis_method in ['lm', 'amm']:
            if 'betas' in res:
                betas = map(list, zip(*res['betas']))
            else:
                betas = [None, None]

            stats_dict['step'] = 0
            cofactors = [stats_dict]
            
            self.add_results(phen_name, dataset, analysis_method, analysis_method, chromosomes, positions, scores, maf_dict['marfs'],
                             maf_dict['mafs'], transformation=transformation,
                             genotype_var_perc=res['var_perc'], beta0=betas[0], beta1=betas[1],
                             cofactors=cofactors, quantiles_dict=quantiles_dict, bh_thres=bh_thres_d['thes_pval'],med_pval=med_pval,ks_stats=ks_res)
        else:
            self.add_results(phen_name, dataset, analysis_method, analysis_method, chromosomes, positions, scores, maf_dict['marfs'],
                    maf_dict['mafs'], transformation=transformation, statistics=res['ds'], quantiles_dict=quantiles_dict, bh_thres=bh_thres_d['thes_pval'],med_pval=med_pval,ks_stats=ks_res)
        self.calculate_ld(phen_name, dataset, transformation, analysis_method, 'results', progress_file_writer=progress_file_writer)
        progress_file_writer.update_progress_bar(progress=1.0, task_status='Done')
        print 'Done!'
        return analysis_method

    def perform_stepwise_gwas(self, phen_name, dataset, transformation, analysis_method, result_name, chromosome, position,
                              call_method_id=82, kinship_method='ibs', progress_file_writer=None):

        """
        Performs GWAS and updates the datastructure.
        """
        self.inc_run_count()
        #if analysis_method not in ['emmax','lm']:
        #    raise Exception("Step-Wise GWAS only possible with emmax or LM")
        snp = ((int(chromosome), int(position)))
        result_group = self.h5file.get_node('/phenotypes/%s/%s/%s/%s' % (phen_name, dataset, transformation, analysis_method))
        result = result_group._f_get_child(result_name)
        cofactors = result._v_attrs.cofactors[:]
        co_var_snps = [(int(factors['chr']), int(factors['pos'])) for factors in cofactors if 'chr' in factors and 'pos' in factors]
        if snp in co_var_snps:
            raise Exception('The SNP %s,%s is already in the result' % (chromosome, position))
        co_var_snps.append(snp)
        co_var_snps = set(co_var_snps)
        #for avail_result in result_group._f_iter_nodes(classname='Table'):
        #   if set(avail_result._v_attrs.cofactors) == co_var_snps:
        #      raise Exception("There is already a result with the selected snps") 


        if 'count' not in result_group._v_attrs:
            count = len(self.h5file.list_nodes(result_group, classname='Table'))
        else:
            count = result_group._v_attrs.count + 1

        new_result_name = "%s" % count
        name = "%s_%s" % (analysis_method, new_result_name)

        import bisect
        import gwa
        import kinship
        if analysis_method not in ['lm', 'amm', 'kw']:
            raise Exception('analysis method %s not supported' % analysis_method)
        if analysis_method == 'kw':
            analysis_method = 'amm'
        progress_file_writer.update_progress_bar(progress=0.0, task_status='Loading phenotype data')
        phen_dict = self.get_phenotype_values(phen_name, dataset, transformation) #Load phenotype
        phend = pd.phenotype_data({1:{'values':phen_dict['mean_value'], 'ecotypes':map(str, phen_dict['ecotype']), 'name':phen_name}})
        phend.convert_to_averages()
        progress_file_writer.update_progress_bar(task_status='Loading genotype data')
        sd = dp.load_snps_call_method(call_method_id=call_method_id, data_format='binary', min_mac=5) #Load SNPs data
        progress_file_writer.update_progress_bar(step=0.05, task_status='Coordinating genotype and phenotype data')
        d = sd.coordinate_w_phenotype_data(phend, 1)
        phen_vals = phend.get_values(1)
        snps = sd.get_snps()
        positions = sd.get_positions()
        chromosomes = []
        for i, (s, c) in enumerate(itertools.izip(sd.snpsDataList, sd.chromosomes)):
            chromosomes.extend([c] * len(s.snps))
        maf_dict = sd.get_mafs()


        kwargs = {}
        if analysis_method == 'amm':
            progress_file_writer.update_progress_bar(progress=0.15, task_status='Retrieving the kinship matrix')
            k = kinship.get_kinship(call_method_id=call_method_id, n_removed_snps=d['n_filtered_snps'], remain_accessions=sd.accessions)
            progress_file_writer.update_progress_bar(progress=0.20, task_status='Performing conditional AMM')
            d = lm.emmax_step(phen_vals, sd, k, list(co_var_snps), progress_file_writer=progress_file_writer, emma_num=200)
            progress_file_writer.update_progress_bar(0.95, 'Processing and saving results')
            res = d['res']
            stats_dict = d['stats']
        elif analysis_method == 'lm':
            progress_file_writer.update_progress_bar(progress=0.20, task_status='Performing Linear regression')
            d = lm.lin_reg_step(phen_vals, sd, list(co_var_snps), progress_file_writer=progress_file_writer)
            progress_file_writer.update_progress_bar(0.95, 'Processing and saving results')
            res = d['res']
            stats_dict = d['stats']
        else:
            raise Exception('analysis method %s not supported' % analysis_method)

        if analysis_method in ['lm', 'amm']:
            if 'betas' in res:
                betas = map(list, zip(*res['betas']))
            else:
                betas = [None, None]
            scores = map(lambda x:-math.log10(x), res['ps'])

            stats_dict['chr'] = snp[0]
            stats_dict['pos'] = snp[1]
            stats_dict['step'] = len(cofactors)
            cofactors.append(stats_dict)

            #Calculate Benjamini-Hochberg threshold
            bh_thres_d = mtcorr.get_bhy_thres(res['ps'], fdr_thres=0.05)
            
            #Calculate Median p-value
            med_pval = agr.calc_median(res['ps'])
            
            #Calculate the Kolmogorov-Smirnov statistic
            ks_res = agr.calc_ks_stats(res['ps'])
            
 
            num_dots = 1000
            max_val = 6
            quantiles = agr.get_quantiles(res['ps'], num_dots=num_dots)
            exp_quantiles = agr._getExpectedPvalueQuantiles_(num_dots)
            log_quantiles = agr.get_log_quantiles(res['ps'], num_dots=num_dots, max_val=max_val)
            exp_log_quantiles = sp.arange(1, num_dots + 1, dtype='single') / (num_dots + 1) * max_val
            quantiles_dict = {'quantiles':quantiles, 'exp_quantiles':exp_quantiles,
            'log_quantiles':log_quantiles, 'exp_log_quantiles':exp_log_quantiles}

            self.add_results(phen_name, dataset, analysis_method, name, chromosomes, positions, scores, maf_dict['marfs'],
                    maf_dict['mafs'], transformation=transformation,
                    genotype_var_perc=res['var_perc'], beta0=betas[0], beta1=betas[1],
                    cofactors=cofactors, result_name=new_result_name, quantiles_dict=quantiles_dict, count=count,
                    bh_thres=bh_thres_d['thes_pval'],med_pval=med_pval,ks_stats=ks_res)
        print 'Done!'
        self.calculate_ld(phen_name, dataset, transformation, analysis_method, new_result_name, progress_file_writer=progress_file_writer)
        progress_file_writer.update_progress_bar(1.0, 'Done')
        return name



#    def add_new_phenotype_file(hdf5_file_name, phenotype_file, phen_name, growth_conditions='', phenotype_scoring='',
#                method_description='', measurement_scale='', is_binary=False):
#        """
#        Initializes the phenotype group for this phenotype and inserts it into the file object.
#        """
#        #Now parsing the phenotype file
#        h5file = tables.open_file(self.filename, mode="r+")
#        print h5file
#        phend = pd.readPhenotypeFile(phenotype_file)
#        _init_phenotype_(h5file, phen_name, growth_conditions=growth_conditions, phenotype_scoring=phenotype_scoring,
#                method_description=method_description, measurement_scale=measurement_scale, is_binary=is_binary)
#        add_phenotype_values(h5file, phen_name, phend.accessions, phend.getPhenVals(1), transformation='raw',
#                accessions=phend.accessionNames, std_dev_values=None, value_comments=None)
#        h5file.flush()
#        h5file.close()
#


def _test_():
    #Load phenotype data..
    import phenotypeData as pd
    import gwaResults as gr
    phed = pd.readPhenotypeFile('/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/phen_raw_092910.tsv')
    pid1 = 1
    phed.filter_accessions_w_missing_data(pid1)
    phen_name = phed.getPhenotypeName(pid1)
    phen_vals = phed.getPhenVals(pid1)
    ecotypes = phed.accessions
    is_binary = phed.isBinary(pid1)

    #Creating the first hdf5 file
    hdf5_file_name_1 = '/Users/bjarni.vilhjalmsson/tmp/test1.hdf5'
    gwa_record = GWASRecord(hdf5_file_name_1)
    gwa_record.init_file()
    gwa_record.add_new_phenotype(phen_name, phen_vals, ecotypes, is_binary=is_binary)
    print "First file is constructed"

    print "Now testing it"
    r = gwa_record.get_phenotype_values(phen_name, 'raw')
    #print r

    phed = pd.readPhenotypeFile('/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/phen_raw_092910.tsv')
    pid2 = 5
    phed.filter_accessions_w_missing_data(pid2)
    phen_name = phed.getPhenotypeName(pid2)
    phen_vals = phed.getPhenVals(pid2)
    ecotypes = phed.accessions
    is_binary = phed.isBinary(pid2)
    gwa_record.add_new_phenotype(phen_name, phen_vals, ecotypes, is_binary=is_binary)

    print "Now testing it"
    r = gwa_record.get_phenotype_values(phen_name, 'raw')
    #print r
    r = gwa_record.get_phenotype_info(phen_name)
    print r

    gwa_record.transform_phenotype('FT10', transformation='sqrt')
    print "Now testing it"
    r = gwa_record.get_phenotype_values(phen_name, 'raw')
    #print r
    r = gwa_record.get_phenotype_info(phen_name)
    print r

    result_file = '/Users/bjarnivilhjalmsson/tmp/pi1_pid5_FT10_emmax_none.pvals'
    res = gr.Result(result_file=result_file, name='FT10')
    res.neg_log_trans()

#    for c in ['chromosomes', 'positions', 'scores', 'marfs', 'mafs', 'genotype_var_perc', 'beta0', \
#        'beta1', 'correlations']:
#        print c, res.snp_results[c][:10]


    gwa_record.add_results(phen_name, 'emmax', res.snp_results['chromosomes'], res.snp_results['positions'],
            res.scores, res.snp_results['marfs'], res.snp_results['mafs'],
            transformation='raw', genotype_var_perc=res.snp_results['genotype_var_perc'],
            beta0=res.snp_results['beta0'], beta1=res.snp_results['beta1'],
            correlation=res.snp_results['correlations'])


    print "Result added."

    print "Now fetching a result."
    res = gwa_record.get_results(phen_name, 'emmax')#, min_mac=15, max_pval=0.01)
    print "Result loaded"
#    for c in ['chromosome', 'position', 'score', 'maf', 'mac', 'genotype_var_perc', 'beta0', \
#        'beta1', 'correlation']:
#        print c, res[c][:10]
    r = gwa_record.get_phenotype_info()
    print r
    s1 = time.time()
    res = gwa_record.get_results_by_chromosome(phen_name, 'emmax')
    print "Result re-loaded"
    secs = time.time() - s1
    if secs > 60:
        mins = int(secs) / 60
        secs = secs - mins * 60
        print 'Took %d mins and %f seconds.' % (mins, secs)
    else:
        print 'Took %f seconds.' % (secs)
    for chromosome in [1, 2, 3, 4, 5]:
        for c in ['position', 'score', 'maf', 'mac', 'genotype_var_perc', 'beta0', \
            'beta1', 'correlation']:
            print c, res[chromosome][c][:10]
    print res['chromosome_ends']
    print res['max_score']
    print gwa_record.get_phenotype_bins(phen_name)
    s1 = time.time()
    gwa_record.perform_gwas('LD', analysis_method='kw')
    secs = time.time() - s1
    if secs > 60:
        mins = int(secs) / 60
        secs = secs - mins * 60
        print 'Took %d mins and %f seconds.' % (mins, secs)
    else:
        print 'Took %f seconds.' % (secs)

    gwa_record.transform_phenotype('LD', transformation='log')
    gwa_record.perform_gwas('LD', analysis_method='emmax', transformation='log')



if __name__ == '__main__':
    _test_()
