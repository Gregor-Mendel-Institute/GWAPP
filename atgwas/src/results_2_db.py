#!/usr/bin/env python2.5
"""
Usage: AddResults2DB.py [OPTIONS] RESULTS_FILE

Option:

	-z ..., --hostname=...	    the hostname, (papaya.usc.edu is default).
	-u ..., --user=...          the username, (otherwise it will ask for it).
	-p ..., --passwd=...	    the password, (otherwise it will ask for it).
	--callMethodID=...          * Use the call_method_id in the DB. 
        --resultsMethodID=...       The results_method_type in the DB. (Default is 1, i.e. association).
        --phenotypeMethodID=...     * The phenotype_method_id from the DB.
        --analysisMethodID=...      * The analysis_method_id from the DB.
	--createdBy=...             * Who created the data (created_by in the DB.)
	--shortName=...             * Short name for the result.
	--methodDescription=...     Description of the method (method_description in the DB).
	--dataDescription=...       Description of the data (data_description in the DB).
	--comment=...               Comments (comment in the DB.) 
	-h, --help	            Show this help

	* these options are required!
Examples:

	
Description:
        Add a result file to the DB.


Important method IDs: 
callMethodID:
        6 : The full dataset.  (This dataset should always be used.)
        7 : First 96 accessions.
	10 : phenotyped (data for which we have phenotypes). (Dataset is possibly outdated)
	12 : phenotyped and merged with 2010 data and MAF <5% removed. (Dataset is possibly outdated)
	13 : phenotyped and merged with 2010 data.  (Dataset is possibly outdated)

phenotypeMethodID:
        The number that the phenotype had e.g. LD is number 1, since it's called 1_LD in the phenotype file.

analysisMethodID:
        1: Kruskal-Wallis.
	2: Fisher's test (binary phenotypes)
	3: Chi-square test.
	4: Emma
	5: Margarita 
	6: Random Forest
"""

import sys, getopt, traceback, env

def _run_():
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)

	long_options_list = ["callMethodID=", "resultsMethodID=", "phenotypeMethodID=", "analysisMethodID=", "createdBy=", "shortName=", "methodDescription=", "dataDescription=", "comment=", "help", "hostname=", "user=", "passwd="]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:u:p:h", long_options_list)

	except:
		traceback.print_exc()
		print sys.exc_info()
		print __doc__
		sys.exit(2)


	hostname = 'papaya.usc.edu'
	user = 'bvilhjal'
	passwd = 'bamboo123'
	resultsFile = args[0]
	help = 0
	callMethodID = None
	resultsMethodID = 1
	phenotypeMethodID = None
	analysisMethodID = None
	createdBy = None
	shortName = None
	methodDescription = ""
	dataDescription = ""
	comment = ""

	for opt, arg in opts:
		if opt in ("-h", "--help"):
			help = 1
			print __doc__
		elif opt in ("-z", "--hostname"):
			hostname = arg
		elif opt in ("-u", "--user"):
			user = arg
		elif opt in ("-p", "--passwd"):
			passwd = arg
		elif opt in ("-z", "--hostname"):
			hostname = arg
		elif opt in ("-u", "--user"):
			user = arg
		elif opt in ("-p", "--passwd"):
			passwd = arg
		elif opt in ("--callMethodID"):
			callMethodID = int(arg)
		elif opt in ("--resultsMethodID"):
			resultsMethodID = int(arg)
		elif opt in ("--phenotypeMethodID"):
			phenotypeMethodID = int(arg)
		elif opt in ("--analysisMethodID"):
			analysisMethodID = int(arg)
		elif opt in ("--createdBy"):
			createdBy = arg
		elif opt in ("--shortName"):
			shortName = arg
		elif opt in ("--methodDescription"):
			methodDescription = arg
		elif opt in ("--dataDescription"):
			dataDescription = arg
		elif opt in ("--comment"):
			comment = arg
		else:
			if help == 0:
				print "Unkown option!!\n"
				print __doc__
			sys.exit(2)

	if not resultsFile:
		if help == 0:
			print "Result file missing!!\n"
			print __doc__
		sys.exit(2)
	raise NotImplementedError
	addResultsToDB(resultsFile, callMethodID, phenotypeMethodID, analysisMethodID,
		       createdBy, shortName, resultsMethodID, methodDescription, dataDescription, comment, commit=commit)




def add_results_to_db(results_file, short_name, call_method_id, phenotype_method_id, analysis_method_id,
		      transformation_method_id, remove_outliers=0, transformation_parameters='null',
		      results_method_type_id=1, method_description="", data_description="", comment="",
		      pseudo_heritability=None):

        #Connect to DB
        import dbutils
	conn = dbutils.connect_to_default_insert()
	cursor = conn.cursor()
	db_result_dir = "/Network/Data/250k/db/results/type_1/"
	if transformation_parameters == 'null':
		transformation_parameters_str = "(transformation_parameters = 'null' OR transformation_parameters IS NULL)"
	else:
		transformation_parameters_str = "transformation_parameters = '%s'" % transformation_parameters


	print "Checking whether result is in DB."
	sql_statement = "SELECT id FROM stock_250k.results_method \
			WHERE phenotype_method_id=%d AND call_method_id=%d AND analysis_method_id=%d \
			AND results_method_type_id=%d AND transformation_method_id=%d AND remove_outliers=%d \
			AND %s"\
			% (phenotype_method_id, call_method_id, analysis_method_id, results_method_type_id, \
				transformation_method_id, remove_outliers, str(transformation_parameters_str))
	#print sql_statement
	cursor.execute(sql_statement)
	row = cursor.fetchone()
	if row:
		results_id = int(row[0])
		print "Found result in DB, with id:", results_id
		db_file = "/Network/Data/250k/db/results/type_1/" + str(results_id) + "_results.tsv"
	else:
		print "Inserting results into DB."
		#Insert info
		sql_statement = "INSERT into stock_250k.results_method \
				(short_name, original_filename, method_description, data_description, \
				phenotype_method_id, call_method_id, results_method_type_id, comment, \
				analysis_method_id, transformation_method_id, remove_outliers, transformation_parameters) \
				values ('%s', '%s', '%s', '%s', %d, %d, %d, '%s', %d, %d, %d, '%s');"\
			       % (short_name, results_file, method_description, data_description, phenotype_method_id, \
				call_method_id, results_method_type_id, comment, analysis_method_id, transformation_method_id, \
				remove_outliers, str(transformation_parameters))
		#print sql_statement
		cursor.execute(sql_statement)

		sql_statement = "SELECT id FROM stock_250k.results_method WHERE short_name like '" + short_name + "'"
		#print sql_statement
		cursor.execute(sql_statement)
		row = cursor.fetchone()
		results_id = int(row[0])

		#Updating filename
		db_file = db_result_dir + str(results_id) + "_results.tsv"
		sql_statement = "UPDATE stock_250k.results_method SET filename='%s' WHERE id=%d" % (db_file, results_id)
		#print sql_statement
		cursor.execute(sql_statement)

		print "Committing transaction (making changes permanent)."
		conn.commit()

	print "Closing connection.\n"
        cursor.close()
	conn.close()


	print "Reading results file:", results_file
	#Convert resultsfile to a tsv file.
	f = open(results_file, "r")
	lines = f.readlines()
	f.close()

        #Write to a designated place.
	results_dir = env.env['db_results_dir'] #'/home/cmbpanfs-01/bvilhjal/data/result_files/'
	file_name = results_dir + str(results_id) + "_results.tsv"
	print "Writing tsv file:", file_name
	try:
		f = open(file_name, "w")
		for line in lines:
			line = line.replace(",", "\t")
			f.write(line)
		f.close()

	except Exception, err_str:
		print 'Failed at writing the result files to the designated results directory:', err_str
		print "Make sure the 'db_results_dir' path is correct in ~/.gwa_config\n"
	if db_result_dir != results_dir:
		print "Remember to copy the result file: scp %s user_name@arabidopsis.gmi.oeaw.ac.at:%s\n"\
			% (file_name, db_result_dir)



if __name__ == '__main__':
	#_test1_()
	_run_()


