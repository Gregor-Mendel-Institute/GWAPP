"""
Contains helper functions for parsing various phenotype files.
"""
import csv
import phenotypeData as pd
from env import *

def load_phentoype_file(filename):
	"""
	Load a FLC type phenotype data file.
	"""
	print "Loading phenotype file:", filename
	f = open(filename, "r")
	reader = csv.reader(f)
	phenotype_names = reader.next()[3:]
	#print phenotype_names
	accession_names = []
	accession_ID = []
	phenotypes = [[] for i in range(len(phenotype_names))]	#[phenotype_name][acc_id]
	for row in reader:
		accession_names.append(row[2].lower())
		accession_ID.append(row[0])
		for i, phen_val in enumerate(row[3:]):
			try:
				p_val = float(phen_val)
			except Exception:
				p_val = "NA"
			#print p_val
			phenotypes[i].append(p_val)
	f.close()
	#print accession_names
	acc_dict = pd._getAccessionToEcotypeIdDict_(accession_names)
	acc_dict["cibc-5"] = 6908
	acc_dict["pla-0"] = 8357
	#print acc_dict
	#accession_names.sort()
	new_phenotypes = [[] for i in range(len(phenotype_names))]
	ecotypes = []
	for acc in acc_dict:
		acc_i = accession_names.index(acc)
		ecotypes.append(acc_dict[acc])
		for i in range(len(phenotype_names)):
			new_phenotypes[i].append(phenotypes[i][acc_i])
	#print new_phenotypes
	#print len(ecotypes)
	#return {"phenotypes":new_phenotypes,"phenotype_names":phenotype_names, "ecotypes":ecotypes}
	phenotypes = map(list, zip(*new_phenotypes))
	phend = pd.PhenotypeData(ecotypes, phenotype_names, phenotypes)
	phend.writeToFile("/tmp/FLC_phenotypes_102809.tsv", delimiter="\t")
	return phend


def load_phentoype_file_Pecinka():
	accession_file = "/Users/bjarnivilhjalmsson/Projects/Ales_Pecinka/NatVar-AP-2010-Feb.csv"
	f = open(accession_file, "r")
	reader = csv.reader(f)
	accession_names = []
	accession_ID = []
	for row in reader:
		accession_names.append(row[1].split()[0].lower())
		accession_ID.append("CS" + row[0][1:])
	f.close()
	print accession_names
	acc_dict = pd._getAccessionToEcotypeIdDict_(accession_names + ["n13", "kno-10", "kno-10", "shahdara", "nd-1"])
	acc_dict["cibc-5"] = 6908
	acc_dict["wa-1"] = 6978
	acc_dict["gu-0"] = 6922
	acc_dict["cs22491"] = acc_dict["n13"]
	acc_dict["knox-10"] = acc_dict["kno-10"]
	acc_dict["knox-18"] = acc_dict["kno-10"]
	acc_dict["shakdara"] = acc_dict["shahdara"]
	acc_dict["wd-1"] = acc_dict["nd-1"]
	print acc_dict


	filename = "/Users/bjarnivilhjalmsson/Projects/Ales_Pecinka/NatVar-AP-2010-Feb_phen.csv"
	#phenotype_names = reader.next()[2:]
	phenotype_names = ["Absolute_root_growth", "Absolute_root_growth_sd", "Percentage_of_root_elongation", "Percentage_of_bent roots",
			   "Percentage_of_dead_plants", "Percentage_of_unaffected_plants", "Percentage_of_average_survival"]
	phenotype_indices = [1, 2, 5, 8, 11, 14, 17]
	phenotype_ecotypes = [0, 0, 4, 7, 10, 13, 16]
	print phenotype_names
	ecotype_ids = [[] for i in range(len(phenotype_names))]
	phenotypes = [[] for i in range(len(phenotype_names))]	#[phenotype_name][acc_id]
	f = open(filename, "r")
	reader = csv.reader(f)
	new_ecotype_ids = set()
	for row in reader:
		print row
		for i, (pi, ei) in enumerate(zip(phenotype_indices, phenotype_ecotypes)):
			if row[ei] != "":
				acc_name = (row[ei].split()[0]).lower()
				if acc_name in acc_dict:
					eid = acc_dict[(row[ei].split()[0]).lower()]
					new_ecotype_ids.add(eid)
					pv = float(row[pi])
					ecotype_ids[i].append(eid)
					phenotypes[i].append(pv)
				else:
					print "Wrong accession name?", acc_name


	new_phenotypes = []
	new_ecotype_ids = list(new_ecotype_ids)
	for i, phen_vals in enumerate(phenotypes):
		new_phen_vals = []
		for ei in new_ecotype_ids:
			if ei in ecotype_ids[i]:
				j = ecotype_ids[i].index(ei)
				new_phen_vals.append(phen_vals[j])
			else:
				new_phen_vals.append('NA')
		new_phenotypes.append(new_phen_vals)
	phenotypes = map(list, zip(*new_phenotypes))
	ecotypes = map(str, new_ecotype_ids)
	phed = pd.PhenotypeData(ecotypes, phenotype_names, phenotypes)
	phed.writeToFile("/Users/bjarnivilhjalmsson/Projects/Ales_Pecinka/phen_pecinka_170310.tsv", delimiter='\t')


def load_phentoype_file_wilczek():
	filename = "/Users/bjarnivilhjalmsson/Projects/Amity_Wilczek/PhenotypeDataWilczek.csv"
	f = open(filename, "r")
	reader = csv.reader(f)
	phenotype_names = reader.next()[2:]
	for i in range(len(phenotype_names)):
		phenotype_names[i] = phenotype_names[i].replace(" ", "_")
	print phenotype_names
	accession_names = []
	accession_ID = []
	for row in reader:
		accession_names.append(row[1].split()[0].lower())
		accession_ID.append(row[0])
	f.close()
	print accession_names
	acc_dict = pd._getAccessionToEcotypeIdDict_(accession_names)#+["n13","kno-10","kno-10","shahdara","nd-1"])
	acc_dict["cibc-5"] = 6908
	acc_dict["wa-1"] = 6978
	acc_dict["gu-0"] = 7149
	acc_dict['Rubezhnoe-1'] = 7323
	print len(acc_dict), acc_dict
	import env
	d250k_file = env.home_dir + "Projects/Data/250k/250K_t54.csv"
	import dataParsers
	d250k_sd = dataParsers.parse_snp_data(d250k_file)
	ecotypes = []
	key_file = "/Users/bjarnivilhjalmsson/Projects/Amity_Wilczek/unique_id_to_ecotype_id.csv"
	f = open(key_file, "w")
	f.write("unique_id, accession_name, ecotype_id, in_250k_data\n")
	for acc, acc_id in zip(accession_names, accession_ID):
		if not acc in acc_dict or acc_id == 'karl27' or acc_id == 'karl05':
			print "(%s, %s) is missing" % (acc, acc_id)
		else:
			ecotype = acc_dict[acc]
			ecotypes.append(ecotype)
			f.write("%s,%s,%s,%s\n" % (acc_id, acc, str(ecotype), str(str(ecotype) in d250k_sd.accessions)))
	f.close()

	#phenotype_names = reader.next()[2:]
	phenotype_indices = range(2, len(phenotype_names) + 2)
	phenotypes = []	#[acc_id][phenotype_name]
	f = open(filename, "r")
	reader = csv.reader(f)
	reader.next()

	for row in reader:
		#print row
		if row[1].split()[0].lower() in acc_dict:
			phen_vals = []
			for pv in row[2:]:
				if pv == "":
					pv = 'NA'
				else:
					pv = float(pv)
				phen_vals.append(pv)
			phenotypes.append(phen_vals)
		else:
			print "Missing:", row[1]

	phed = pd.PhenotypeData(ecotypes, phenotype_names, phenotypes)
	phed.writeToFile("/Users/bjarnivilhjalmsson/Projects/Amity_Wilczek/phen_wilzcek_050710.tsv", delimiter='\t')
	phed.writeToFile("/Users/bjarnivilhjalmsson/Projects/Amity_Wilczek/phen_wilzcek_050710.csv", delimiter=',')



def load_phentoype_file_dilkes():
	filename = env['phen_dir'] + 'dilkes_metabolites.csv'
	f = open(filename, "r")
	print f.next()
	header = f.next()
	phen_names = map(str.strip, header.split(',')[2:])
	print phen_names
	pids = range(len(phen_names))
	accessions = []
	phen_dict = {}
	for pid, name in zip(pids, phen_names):
		phen_dict[pid] = {'name':name, 'ecotypes':[], 'values':[]}
	for line in f:
		l = map(str.strip, line.split(','))
		accessions.append(l[1].lower())
		for i in pids:
			phen_dict[i]['values'].append(float(l[2 + i]))

	f.close()
	acc_dict = pd.get_250K_accession_to_ecotype_dict()
#	acc_dict['buckhorn'] = [0, 0, 0, 0, 7033]
	acc_dict['sakhdara'] = [0, 0, 0, 0, 6962]
	ecotypes = []
	for acc in accessions:
		if not acc in acc_dict:
			print "%s is missing in dictionary" % acc
		else:
			ecotype = acc_dict[acc][4]
			ecotypes.append(ecotype)

	for pid, name in zip(pids, phen_names):
		phen_dict[pid]['ecotypes'] = ecotypes

	phed = pd.phenotype_data(phen_dict)
	phed.write_to_file(env['phen_dir'] + 'b_dilkes_metabolites.csv')




def load_phentoype_file_nc_resistance():
	filename = "/Users/bjarnivilhjalmsson/Summary_results_330Arabidopsis_accessions.csv"
	f = open(filename, "r")
	line = map(str.strip, f.next().split(','))
	phenotype_names = line[-2:]
	print phenotype_names
	phenotypes = []
	accession_names = []
	full_accession_names = []
	for l in f:
		line = map(str.strip, l.split(','))
		accession_names.append(line[0].lower())
		full_accession_names.append(line[2].lower())
	f.close()
	print accession_names
	acc_dict = pd.get_accession_to_ecotype_id_dict(accession_names)#+["n13","kno-10","kno-10","shahdara","nd-1"])
#	acc_dict["cibc-5"] = 6908
#	acc_dict["wa-1"] = 6978
#	acc_dict["gu-0"] = 7149
#	acc_dict['Rubezhnoe-1'] = 7323
	print len(acc_dict), acc_dict
	import env
	d250k_file = env.home_dir + "Projects/Data/250k/250K_t54.csv.binary"
	import dataParsers
	d250k_sd = dataParsers.parse_binary_snp_data(d250k_file)
	ecotypes = []
	key_file = "/Users/bjarnivilhjalmsson/Projects/Amity_Wilczek/unique_id_to_ecotype_id.csv"
	f = open(key_file, "w")
	f.write("unique_id, accession_name, ecotype_id, in_250k_data\n")
	for acc, acc_id in zip(accession_names, full_accession_names):
		if not acc in acc_dict or acc_id == 'karl27' or acc_id == 'karl05':
			print "(%s, %s) is missing" % (acc, acc_id)
		else:
			ecotype = acc_dict[acc]
			ecotypes.append(ecotype)
			f.write("%s,%s,%s,%s\n" % (acc_id, acc, str(ecotype), str(str(ecotype) in d250k_sd.accessions)))

	#phenotype_names = reader.next()[2:]
	phenotypes = []	#[acc_id][phenotype_name]
	f = open(filename, "r")

	for l in f:
		line = map(str.strip, l.split(','))
		if line[0].lower() in acc_dict:
			phen_vals = []
			for pv in line[-2:]:
				if pv == "NA":
					pv = 'NA'
				else:
					pv = float(pv)
				phen_vals.append(pv)
			phenotypes.append(phen_vals)
		else:
			print "Missing:", line[0]

	phed = pd.PhenotypeData(ecotypes, phenotype_names, phenotypes)
	phed.insert_into_DB(growth_condition='Field', biology_category_id='2')
	phed.writeToFile("/Users/bjarnivilhjalmsson/Projects/Amity_Wilczek/nc14_resistance_091610.tsv", delimiter='\t')



def load_phentoype_file_nc_resistance_2():
	filename = "/Users/bjarnivilhjalmsson/Projects/Albugo_laibachii_nc14.csv"
	f = open(filename, "r")
	line = map(str.strip, f.next().split(','))
	phenotype_names = line[-1:]
	print phenotype_names
	phenotypes = []
	accession_names = []
	for l in f:
		line = map(str.strip, l.split(','))
		accession_names.append(line[0].lower())
	f.close()
	print accession_names
	acc_dict = pd.get_accession_to_ecotype_id_dict(accession_names)#+["n13","kno-10","kno-10","shahdara","nd-1"])
#	acc_dict["cibc-5"] = 6908
#	acc_dict["wa-1"] = 6978
#	acc_dict["gu-0"] = 7149
#	acc_dict['Rubezhnoe-1'] = 7323
	print len(acc_dict), acc_dict
	import env
	d250k_file = env.home_dir + "Projects/Data/250k/250K_t54.csv.binary"
	import dataParsers
	d250k_sd = dataParsers.parse_binary_snp_data(d250k_file)
	ecotypes = []
	for acc in accession_names:
		if not acc in acc_dict:
			print "%s is missing" % (acc)
		else:
			ecotype = acc_dict[acc]
			ecotypes.append(ecotype)

	phenotypes = []	#[acc_id][phenotype_name]
	f = open(filename, "r")

	for l in f:
		line = map(str.strip, l.split(','))
		if line[0].lower() in acc_dict:
			phen_vals = []
			for pv in line[-1:]:
				if pv == "NA":
					pv = 'NA'
				else:
					pv = float(pv)
				phen_vals.append(pv)
			phenotypes.append(phen_vals)
		else:
			print "Missing:", line[0]

	phed = pd.PhenotypeData(ecotypes, phenotype_names, phenotypes)
	phed.insert_into_DB(growth_condition='Field', biology_category_id='2')
	phed.writeToFile("/Users/bjarnivilhjalmsson/Projects/Amity_Wilczek/nc14_resistance_96accessions_092810.tsv", \
			delimiter='\t')



def load_phentoype_file_nc_resistance_3():
	filename = "/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/20dd5_330.csv"
	with open(filename) as f:
		line = map(str.strip, f.next().split(','))
		phenotype_names = line[-1:]
		print phenotype_names
		phenotypes = []
		accession_names = []
		ecotypes = []
		full_accession_names = []
		for l in f:
			line = map(str.strip, l.split(','))
			accession_names.append(line[1].lower())
			ecotypes.append(line[0])
			full_accession_names.append(line[5].lower())
			phenotypes.append(line[4])

	print accession_names
	acc_dict = pd.get_accession_to_ecotype_id_dict(accession_names)#+["n13","kno-10","kno-10","shahdara","nd-1"])
#	acc_dict["cibc-5"] = 6908
#	acc_dict["wa-1"] = 6978
#	acc_dict["gu-0"] = 7149
#	acc_dict['Rubezhnoe-1'] = 7323
	print len(acc_dict), acc_dict
	ets = []
	phen_vals = []
	for acc1, acc2, et, pt in zip(accession_names, full_accession_names, ecotypes, phenotypes):
		if acc1 in acc_dict:
			ecotype = acc_dict[acc1]

			if str(ecotype) != et and et != 'NA':
				print "Ecotype mismatch.. %s, %s, %s, %s" % (unicode(acc1, "latin-1"),
									unicode(acc2, "latin-1"), et, ecotype)
			else:
				et = ecotype
			if et != 'NA' and et != '':
				ets.append(et)
				if not pt in ['R', 'S']: print pt
				phen_vals.append(0 if pt == 'R' else 1)
	print len(phen_vals)

	phen_dict = {1:{'name':'resistance_20dd5', 'ecotypes':ets, 'values':phen_vals}}
	phed = pd.phenotype_data(phen_dict)
	phed.write_to_file('resistance_20dd5.csv', ',')



def load_phentoype_file_bergelsson():
	import env
	filename = "/Users/bjarnivilhjalmsson/Projects/Joy_Bergelsson/bergelsson_rosette_glucs.csv"
	f = open(filename, "r")
	reader = csv.reader(f)
	phenotype_names = reader.next()[2:]
	for i in range(len(phenotype_names)):
		phenotype_names[i] = phenotype_names[i].replace(" ", "_")
		phenotype_names[i] = 'jb_' + phenotype_names[i]
	print phenotype_names
	accession_names = []
	accession_ID = []
	for row in reader:
		accession_names.append(row[0].split()[0].lower())
		accession_ID.append(row[1])
	f.close()
	print accession_names
	#acc_dict = pd._getAccessionToEcotypeIdDict_(accession_names)#+["n13","kno-10","kno-10","shahdara","nd-1"])
	e_info_dict = pd._getEcotypeIdInfoDict_()
	ei_2_tgei = pd._getEcotype2TgEcotypeDict_()
	#print len(acc_dict),acc_dict
	ecotypes = []
        uncertain_list = []
	for acc, acc_id in zip(accession_names, accession_ID):
		#if not acc in acc_dict:
		if not int(acc_id) in ei_2_tgei:
			print "(%s, %s) is missing in dictionary" % (acc, acc_id)
			a_id = int(acc_id)
			if a_id in e_info_dict:
				e_info = e_info_dict[a_id]
				print "Guessing that it's:", e_info
			else:
				print "No good guess for it.  Look it up!!\n"
			#acc_dict[acc] = acc_id
			ecotypes.append(acc_id)
		else:
			#ecotype = acc_dict[acc]
			ecotype = ei_2_tgei[int(acc_id)]
			ecotypes.append(ecotype)
	phenotype_indices = range(2, len(phenotype_names) + 2)
	phenotypes = []	#[acc_id][phenotype_name]
	f = open(filename, "r")
	reader = csv.reader(f)
	reader.next()

	print len(set(accession_ID)), len(set(ecotypes))

	for row in reader:
		#print row
		#if row[0].split()[0].lower() in acc_dict:
			phen_vals = []
			for pv in row[2:]:
				if pv == "" or pv == 'NA':
					pv = 'NA'
				else:
					pv = float(pv)
				phen_vals.append(pv)
			if len(phen_vals) != len(phenotype_names):
				import pdb;
				pdb.set_trace()
			phenotypes.append(phen_vals)
		#else:
		#	print "Missing:",row[0]


	phed = pd.PhenotypeData(ecotypes, phenotype_names, phenotypes)
	phed.writeToFile("/Users/bjarnivilhjalmsson/Projects/Joy_Bergelsson/phen_bergelsson_051710.tsv", delimiter='\t')
	phed.writeToFile("/Users/bjarnivilhjalmsson/Projects/Joy_Bergelsson/phen_bergelsson_051710.csv", delimiter=',')



def load_duszynska_file4():
	"""
	Loads the heterosis data.
	"""
	fn = env['home_dir'] + 'Projects/duszynska_data/seed_size_heterosis.csv'
	acc_dict = pd.get_250K_accession_to_ecotype_dict()
	phen_dict = {}

	name_dict = {'knox-18':'kno-18', 'knox-10':'kno-10', 'kas-1':'kas-2', 'pu-2-7':'pu2-7', 'cs22491':'n13',
			'shahdara':'sha'}

	with open(fn) as f:
		ets = []
		header = f.next()
		phen_names = map(str.strip, header.split(','))
		et_indices = [0, 3, 6, 9]
		phen_names = [phen_names[1], phen_names[4], phen_names[7], phen_names[10]]
		for i, pn in zip([1, 2, 3, 4], phen_names):
			phen_dict[i] = {'name': pn, 'values':[], 'ecotypes':[]}
		for line in f:
			l = map(str.strip, line.split(','))
			for e_i, pid in zip(et_indices, [1, 2, 3, 4]):
				acc = l[e_i].lower()
				if acc in name_dict:
					acc = name_dict[acc]
				if not acc in acc_dict:
					print "(%s) is missing in dictionary" % (acc)
				else:
					phen_dict[pid]['ecotypes'].append(acc_dict[acc][4])
					phen_dict[pid]['values'].append(float(l[e_i + 1]))

	phed = pd.phenotype_data(phen_dict)
	phed.write_to_file(env['phen_dir'] + 'duszynska_heterosis_data.csv')


def load_duszynska_file3():
	fn1 = env['home_dir'] + 'Projects/duszynska_data/ANU_proportions_4x2_male.csv'
	fn2 = env['home_dir'] + 'Projects/duszynska_data/ANU_proportions_2x4_female.csv'
	acc_dict = pd.get_250K_accession_to_ecotype_dict()
	phen_names = []
	phen_dict = {}

	name_dict = {'knox-18':'kno-18', 'knox-10':'kno-10', 'kas-1':'kas-2', 'pu-2-7':'pu2-7', 'cs22491':'n13',
			'shahdara':'sha'}

	with open(fn1) as f:
		ets = []
		header = f.next()
		phen_names = map(str.strip, header.split(',')[1:])
		for i, pn in zip([1, 2, 3, 4, 5, 6, 7], phen_names):
			phen_dict[i] = {'name':'male_4x2_' + pn, 'values':[]}
		for line in f:
			l = map(str.strip, line.split(','))
			acc = l[0].lower()
			if acc in name_dict:
				acc = name_dict[acc]
			if not acc in acc_dict:
				print "(%s) is missing in dictionary" % (acc)
			else:
				ets.append(acc_dict[acc][4])
				for pid in [1, 2, 3, 4, 5, 6, 7]:
					phen_dict[pid]['values'].append(float(l[pid]))
		for pid in [1, 2, 3, 4, 5, 6, 7]:
			phen_dict[pid]['ecotypes'] = ets[:]


	with open(fn2) as f:
		ets = []
		header = f.next()
		phen_names = map(str.strip, header.split(',')[1:])
		for i, pn in zip([8, 9, 10, 11, 12, 13, 14], phen_names):
			phen_dict[i] = {'name':'female_2x4_' + pn, 'values':[]}
		for line in f:
			l = map(str.strip, line.split(','))
			acc = l[0].lower()
			if acc in name_dict:
				acc = name_dict[acc]
			if not acc in acc_dict:
				print "(%s) is missing in dictionary" % (acc)
			else:
				ets.append(acc_dict[acc][4])
				for pid in [8, 9, 10, 11, 12, 13, 14]:
					phen_dict[pid]['values'].append(float(l[pid - 7]))
		for pid in [8, 9, 10, 11, 12, 13, 14]:
			phen_dict[pid]['ecotypes'] = ets[:]

	phed = pd.phenotype_data(phen_dict)
	phed.write_to_file(env['phen_dir'] + 'duszynska_data_new.csv')



def load_duszynska_file2():
	fn1 = env['home_dir'] + 'Projects/duszynska_data/male_data_proportion_AN.csv'
	fn2 = env['home_dir'] + 'Projects/duszynska_data/female_data_proportion_AN.csv'
	acc_dict = pd.get_250K_accession_to_ecotype_dict()
	phen_names = []
	phen_dict = {}

	name_dict = {'knox-18':'kno-18', 'knox-10':'kno-10', 'kas-1':'kas-2', 'pu-2-7':'pu2-7', 'cs22491':'n13',
			'shahdara':'sha'}

	with open(fn1) as f:
		ets = []
		header = f.next()
		phen_names = map(str.strip, header.split(',')[1:])
		for i, pn in zip([1, 2, 3], phen_names):
			phen_dict[i] = {'name':'male_' + pn, 'values':[]}
		for line in f:
			l = map(str.strip, line.split(','))
			acc = l[0].lower()
			if acc in name_dict:
				acc = name_dict[acc]
			if not acc in acc_dict:
				print "(%s) is missing in dictionary" % (acc)
			else:
				ets.append(acc_dict[acc][4])
				for pid in [1, 2, 3]:
					phen_dict[pid]['values'].append(float(l[pid]))
		for pid in [1, 2, 3]:
			phen_dict[pid]['ecotypes'] = ets[:]


	with open(fn2) as f:
		ets = []
		header = f.next()
		phen_names = map(str.strip, header.split(',')[1:])
		for i, pn in zip([4, 5, 6], phen_names):
			phen_dict[i] = {'name':'female_' + pn, 'values':[]}
		for line in f:
			l = map(str.strip, line.split(','))
			acc = l[0].lower()
			if acc in name_dict:
				acc = name_dict[acc]
			if not acc in acc_dict:
				print "(%s) is missing in dictionary" % (acc)
			else:
				ets.append(acc_dict[acc][4])
				for pid in [4, 5, 6]:
					phen_dict[pid]['values'].append(float(l[pid - 3]))
		for pid in [4, 5, 6]:
			phen_dict[pid]['ecotypes'] = ets[:]

	phed = pd.phenotype_data(phen_dict)
	phed.write_to_file(env['phen_dir'] + 'duszynska_data.csv')




def load_phentoype_file_duszynska():
	fn1 = env['phen_dir'] + 'seed_size_2n.csv'
	fn2 = env['phen_dir'] + 'seed_size_3n_2x4.csv'
	fn3 = env['phen_dir'] + 'seed_size_3n_4x2.csv'
	fn4 = env['phen_dir'] + 'seed_size_spss.csv'
	fns = [fn1, fn2, fn3, fn4]
	phen_names = ['seed_size_2n', 'seed_size_3n_2x4', 'seed_size_3n_4x2', 'seed_size_spss']
	phen_dict = {}
	for i, pn in enumerate(phen_names):
		phen_dict[i + 1] = {'name':pn }
	accs_list = []
	for i, fn in enumerate(fns):
		f = open(fn, "r")
		acc_dict = pd.get_250K_accession_to_ecotype_dict()
		ecotypes = []
		values = []
		print f.next()
		for line in f:
			l = map(str.strip, line.split(','))
			acc = l[0].lower()
			if not acc in acc_dict:
				print "(%s) is missing in dictionary" % (acc)
			else:
				ecotypes.append(acc_dict[acc][4])
				values.append(float(l[1]))
		f.close()
		phen_dict[i + 1]['ecotypes'] = ecotypes
		phen_dict[i + 1]['values'] = values


	phed = pd.phenotype_data(phen_dict)
	phed.write_to_file(env['phen_dir'] + 'seed_size.csv')





def load_phentoype_file_riha():
	filename = env['phen_dir'] + 'telomere_lengths_192_raw.csv'
	f = open(filename, "r")
	phen_name = 'telomere_length'
	accession_names = []
	accession_ids = []
	parent_ids = []
	phen_vals = []
	print f.next()
	for line in f:
		l = map(str.strip, line.split(','))
		parent_ids.append(l[0])
		acc_l = l[1].split()
		acc_name = acc_l[0]
		if len(acc_l) > 1:
			acc_id = acc_l[1]
		else:
			acc_id = ''
		accession_names.append(acc_name.lower())
		accession_ids.append(acc_id)
		phen_vals.append(float(l[2]))

	f.close()
	print accession_names
	acc_dict = pd.get_250K_accession_to_ecotype_dict()
	acc_dict['buckhorn'] = [0, 0, 0, 0, 7033]
	acc_dict['shahdara'] = [0, 0, 0, 0, 6962]
	ecotypes = []
        uncertain_list = []
        new_phen_vals = []
	for acc, par_id, acc_id, phen_val in zip(accession_names, parent_ids, accession_ids, phen_vals):
		if not acc in acc_dict:
			print "(%s, %s, %s) is missing in dictionary" % (acc, par_id, acc_id)
		else:
			ecotype = acc_dict[acc][4]
			ecotypes.append(ecotype)
			new_phen_vals.append(phen_val)

	print len(set(accession_names)), len(set(ecotypes))
	phen_dict = {1:{'name':phen_name, 'ecotypes':ecotypes, 'values':new_phen_vals}}

	phed = pd.phenotype_data(phen_dict)
	phed.write_to_file(env['phen_dir'] + 'telomere_lengths_192.csv')


def load_gene_expression_traits():
	filename = '/Users/bjarni.vilhjalmsson/Projects/Data/rna_seq/gene_expression_table_20110208.tsv'
	ecotypes = []
	et_ids = []
	import scipy as sp
	with open(filename, "r") as f:
		i = 0
		for line in f:
			if line[0] != '#': break
			l = line.split('\t')
			if l[7] == '16C':
				ecotypes.append(l[5])
				et_ids.append(i)
			i += 1
		gene_ids = line.split('\t')
		print len(ecotypes), len(set(ecotypes))
		print ecotypes
		et_dict = {'Col-0':'6909', 'Col':'6909', 'Ler':'6932', 'Ws-0':'6980'}
		ets = []
		for et in ecotypes:
			if et in et_dict:
				ets.append(et_dict[et])
			else:
				ets.append(et)
		print ets
		phen_dict = {}
		num_const_phen = 0
		phen_i = 1
		for i, line in enumerate(f): #For each gene
			l = line.split()
			phen_name = l[0]
			phen_vals = map(float, l[1:])
			phen_vals = [phen_vals[i] for i in et_ids]
			if len(phen_vals) != len(ets):
				raise Exception('Arrg')
			if len(sp.unique(phen_vals)) > 1:
				phen_dict[phen_i] = {'name':phen_name, 'ecotypes':ets, 'values':phen_vals}
				phen_i += 1
			else:
				num_const_phen += 1
	print 'Total number of gene expressions was %d, of which %d were constant (removed), leaving %d gene expressions.' \
		% ((phen_i - 1) + num_const_phen, num_const_phen, phen_i - 1)

	phed = pd.phenotype_data(phen_dict)
	phed.write_to_file(env['phen_dir'] + 'rna_seq_020811_16C.csv')



def load_gene_expression_traits_2():
	import scipy as sp
	filename = env['home_dir'] + \
			'Projects/Data/rna_seq/expression_matrices_upload_8_01_2011/' + \
			'expression_matrix_wSNPmap_7_29_2011-bioreps_combined_cov_filter-normalized.txt'
	print 'Loading file:', filename
	ets = {'10C':[], '16C':[]}
	i_map = {}
	expressions_dict = {}
	with open(filename, "r") as f:
		i = 0
		line = f.next()
		l = map(str.strip, line.split())
		for i in range(2, len(l)):
			e_t_l = l[i].split('_')
			et = e_t_l[0][1:]
			t = e_t_l[1]
			i_map[i] = t
			if int(et) == 2:
				et = '6932'
			elif int(et) == 3:
				et = '6980'
			if int(et) == 4:
				et = 'ALyr'
			elif int(et) == 5:
				et = 'ACap'
			ets[t].append(et)

		for line in f:
			l = map(str.strip, line.split())
			gene_name = l[0]
			gene_type = l[1]
			d = {'gene_type':gene_type, '10C':[], '16C':[]}
			for i in range(2, len(l)):
				t = i_map[i]
				val = float(l[i])
				d[t].append(val)
			expressions_dict[gene_name] = d
	print 'File was parsed, now constructing phenotype object'
	phen_dict_10C = {}
	phen_dict_16C = {}
	phen_i = 1
	for gene_name in expressions_dict:
		values = expressions_dict[gene_name]['10C']
		if len(sp.unique(values)) > 1:
			phen_dict_10C[phen_i] = {'name':gene_name, 'ecotypes':ets['10C'],
						'values':expressions_dict[gene_name]['10C']}
		phen_i += 1

	phen_i = 1
	for gene_name in expressions_dict:
		values = expressions_dict[gene_name]['16C']
		if len(sp.unique(values)) > 1:
			phen_dict_16C[phen_i] = {'name':gene_name, 'ecotypes':ets['16C'],
						'values':expressions_dict[gene_name]['16C']}
		phen_i += 1


	phed_10C = pd.phenotype_data(phen_dict_10C)
	print 'Phenotype object constructed with %d phenotypes, now writing to phenotype file' % len(phen_dict_10C)
	phed_10C.write_to_file(env['phen_dir'] + 'rna_seq_081411_10C.csv')
	phed_16C = pd.phenotype_data(phen_dict_16C)
	print 'Phenotype object constructed with %d phenotypes, now writing to phenotype file' % len(phen_dict_16C)
	phed_16C.write_to_file(env['phen_dir'] + 'rna_seq_081411_16C.csv')


def load_gene_expression_traits_3(temperature='10C'):
	import scipy as sp
	filename = env['home_dir'] + \
			'/Projects/Data/rna_seq/expression_variance_stabilized_11_08_11/' + \
			'expr_%s_merged.csv' % temperature
	phen_dict = {}
	phen_i = 1
	with open(filename) as f:
		header = (f.next().strip()).split(',')
		ets = map(lambda x: x[1:], header[1:])
		for i in range(len(ets)):
			if int(ets[i]) == 2:
				ets[i] = '6932'
			elif int(ets[i]) == 3:
				ets[i] = '6980'
			if int(ets[i]) == 4:
				ets[i] = 'ALyr'
			elif int(ets[i]) == 5:
				ets[i] = 'ACap'
		print ets
		for l in f:
			line = (l.strip()).split(',')
			gene_name = line[0]
			vals = map(float, line[1:])
			phen_dict[phen_i] = {'name':gene_name, 'ecotypes':ets, 'values':vals}
			phen_i += 1

	phed = pd.phenotype_data(phen_dict)
	print 'Phenotype object constructed with %d phenotypes, now writing to phenotype file' % len(phen_dict)
	phed.write_to_file(env['phen_dir'] + 'rna_seq_vs_081411_%s.csv' % temperature)




def load_total_expressions():
	import scipy as sp
	filename = env['home_dir'] + \
			'Projects/Data/rna_seq/expression_matrix_upload_data_5_10_2011/mapping_files/' + \
			'all_expression_matrix_5_09_2011_flagged_removed_libID_bioreps_combined_cov_filter.txt'
	print 'Loading file:', filename
	ets = {'10C':[], '16C':[]}
	i_map = {}
	expressions_dict = {}
	with open(filename, "r") as f:
		i = 0
		line = f.next()
		l = map(str.strip, line.split())
		for i in range(2, len(l)):
			e_t_l = l[i].split('_')
			et = e_t_l[0][1:]
			t = e_t_l[1]
			i_map[i] = t
			ets[t].append(et)

		line = f.next()
		l = map(str.strip, line.split())
		gene_name = l[0]
		gene_type = l[1]
		d = {'gene_type':gene_type, '10C':[], '16C':[]}
		for i in range(2, len(l)):
			t = i_map[i]
			val = float(l[i])
			d[t].append(val)
		expressions_dict[gene_name] = d
	print 'File was parsed, now constructing phenotype object'
	phen_dict = {}
	phen_i = 1
	for t in ['10C', '16C']:
		values = expressions_dict['GENE'][t]
		phen_dict[phen_i] = {'name':'total_expression_%s' % t, 'ecotypes':ets[t], 'values':values}
		phen_i += 1


	phed = pd.phenotype_data(phen_dict)
	print 'Phenotype object constructed with %d phenotypes, now writing to phenotype file' % len(phen_dict)
	phed.write_to_file(env['phen_dir'] + 'rna_seq_061611_total.csv')





def parse_NFBC_traits():
	phen_dict = {}
	height_file = env['data_dir'] + 'NFBC_20091001/pheno.Height'
	with open(height_file) as f:
		f.next()
		ets = []
		values = []
		for line in f:
			l = line.split()
			ets.append(int(l[0]))
			values.append(float(l[2]))
	phen_dict[1] = {'name':'height', 'ecotypes':ets, 'values':values}

	metabolite_file = env['data_dir'] + 'NFBC_20091001/MetaboPheno.txt'
	with open(metabolite_file) as f:
		line = f.next()
		l = map(str.strip, line.split())
		phen_names = l[2:]
		for i, pname in enumerate(phen_names):
			phen_dict[i + 2] = {'name':pname, 'ecotypes':[], 'values':[]}
		ets = []
		values = []
		for line in f:
			l = line.split()
			et = int(l[0])
			for i, v in enumerate(l[2:]):
				try:
					val = int(v)
				except Exception:
					val = float(v)
				if val != -9 and type(val) != int:
					phen_dict[i + 2]['ecotypes'].append(et)
					phen_dict[i + 2]['values'].append(val)
	phed = pd.phenotype_data(phen_dict)
	phed.write_to_file(env['data_dir'] + 'NFBC_20091001/phenotype.scv')
	return phed



def load_skin_color_traits():
	dir_prefix = env['home_dir'] + 'Projects/data/skin_eye_color/'
	#dir_prefix = env.env['home_dir'] + 'Projects/Data/Skin_color/'
	filename = dir_prefix + 'CV685-skin_eye_color.txt'
	d = {1:{'name':'skin_color', 'ecotypes':[], 'values':[]}, 2:{'name':'eye_color', 'ecotypes':[], 'values':[]}}
	sc_vals = []
	ec_vals = []
	sc_iids = [] #individual IDs (ecotypes)
	with open(filename) as f:
		print f.next()
		for line in f:
			l = line.split()
			if int(float(l[2])) != -9:
				d[1]['values'].append(float(l[2]))
				d[1]['ecotypes'].append(l[1])
			if int(float(l[3])) != -9:
				d[2]['values'].append(float(l[3]))
				d[2]['ecotypes'].append(l[1])
	phed = pd.phenotype_data(d)
	phed.write_to_file(dir_prefix + 'phenotypes.csv')
	return phed


def load_genome_size_factors():
	phen_file = env['phen_dir'] + 'measures_gene_model.tsv'
	phen_dict = {}
	ets = []
	with open(phen_file) as f:
		line = f.next()
		names = line.split()[1:]
		for i, name in enumerate(names):
			phen_dict[i + 1] = {'values':[], 'name':name, 'transformation':None}
		for line in f:
			l = line.split()
			ets.append(l[0])
			for i, v in enumerate(l[1:]):
				phen_dict[i + 1]['values'].append(float(v))
	for i in phen_dict:
		phen_dict[i]['ecotypes'] = ets[:]
		phen_dict[i]['raw_values'] = phen_dict[i]['values']
	phed = pd.phenotype_data(phen_dict)
	phed.write_to_file(env['phen_dir'] + 'measures_gene_model.csv')



def _run_():
	pd = load_phentoype_file("/Users/bjarnivilhjalmsson/Projects/FLC_analysis/data_102509/FLC_soil_data_102509.csv")
	phenotypes = pd["phenotypes"]
	phenotype_names = pd["phenotype_names"]
	new_phenotype_names = []
	for pn in phenotype_names:
		new_phenotype_names.append("FLC_v2_" + pn)
	phenotype_names = new_phenotype_names
	ecotypes = pd["ecotypes"]
	method_descriptions = ["without vernalization", "4 weeks after vernalization in cold", "10 days after 4 weeks of vernalization", "30 days after 4 weeks of vernalization", "Ratio between 4wT0 and NV", "Ratio between 4wT10 and NV", "Ratio between 4wT30 and 4wT0", "Ratio between 4wT30 and 4wT10"]
	data_descriptions = ["Stable soil growth condition and new Roch PCR machine." for i in range(len(phenotype_names))]
	add_phenotypes_to_db(phenotypes, phenotype_names, ecotypes, method_ids=range(352, 352 + len(phenotype_names)), method_descriptions=method_descriptions, data_descriptions=data_descriptions)


if __name__ == "__main__":
	#_run_()
	#load_phentoype_file("/Users/bjarnivilhjalmsson/Projects/FLC_analysis/data_102509/FLC_soil_data_102509.csv")
	#load_phentoype_file_Pecinka()
	#load_phentoype_file_wilczek()
	load_genome_size_factors()
	#load_gene_expression_traits_3('16C')
	#load_phentoype_file_nc_resistance_3()
	print "Done!"


