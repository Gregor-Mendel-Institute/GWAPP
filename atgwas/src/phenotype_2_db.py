"""
Usage: phenotype_2_db.py [OPTIONS]

Option:

	-h		            		Show this help
	-f ...,					Phenotype file!
	-z ...,					The hostname, (default is specified in your .gwa_config file).
	-u ...,					The username, (default is specified in your .gwa_config file).
	-p ...,					The password, (default is specified in your .gwa_config file).
	-i ...,					Phenotype ID(s) to insert (default is all)
	-m ...,                      		DB phenotype ID(s) (overwrites the ID in DB if it exists and deletes all existing values), 
						appends by default.
	-d 					Dry run, doesn't commit changes to the db.
						
        --phenotype_scoring=...			Phenotype scoring text. 
        --method_description=...		Method description text. 
        --growth_condition=...			Growth_condition text. 
	--citations=...				Citation text. 
	--data_description=...			Data description text. 
	--transformation_description=...	Transformation description text. 
	--comment=...				Comments (comment in the DB.). 
	--biology_category_id=...		Default is '', look up in DB to see what number fits your trait.



NOTE THAT PARSER CAN'T HANDLE SPACES IN OPTIONS!!  (If you need spaces, then you can edit the DB.)

ONLY INSERTS AVERAGES... FOR NOW
"""

import env
import sys
import traceback
import phenotypeData as pd
import getopt

def _parse_ids_(id_arg_str):
	t_ids = id_arg_str.split(',')
	ids = []
	for s in t_ids:
		if '-' in s:
			id_range = map(int, s.split('-'))
		        for id in range(id_range[0], id_range[1] + 1):
		        	ids.append(id)
		else:
			ids.append(int(s))
	return ids

def parse_parameters():
	'Parse the parameters into a dict, etc.'
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(0)

	long_options_list = ["comment=", 'phenotype_scoring', 'method_description=', 'growth_condition=', 'citations=',
			'data_description', 'transformation_description=', 'data_type=', 'biology_category_id=']
	try:
		opts, args = getopt.getopt(sys.argv[1:], "f:z:u:p:i:m:hd", long_options_list)

	except:
		traceback.print_exc()
		print __doc__
		sys.exit(2)


	p_dict = {'phen_filename':None, 'pids':None, 'method_ids':None, 'db_user':env.env['db_user'],
		'db_passwd':env.env['db_passwd'], 'host':env.env['default_insert_db'], 'comment':'',
		'phenotype_scoring':'', 'method_description':'', 'growth_condition':'', 'citations':'',
		'data_description':'', 'transformation_description':'None',
		'biology_category_id':4, 'dry_run':False}


	for opt, arg in opts:
		if opt in ("-h"):
			print __doc__
			return
		elif opt in ('-f'): p_dict['phen_filename'] = arg
		elif opt in ('-i'): p_dict['pids'] = _parse_ids_(arg)
		elif opt in ('-z'): p_dict['host'] = arg
		elif opt in ('-m'): p_dict['method_ids'] = _parse_ids_(arg)
		elif opt in ('-u'): p_dict['db_user'] = arg
		elif opt in ('-p'): p_dict['db_passwd'] = arg
		elif opt in ('-d'): p_dict['dry_run'] = True
		elif opt in ("--comment"): p_dict['comment'] = arg
		elif opt in ("--phenotype_scoring"): p_dict['phenotype_scoring'] = arg
		elif opt in ("--method_description"): p_dict['method_description'] = arg
		elif opt in ("--growth_condition"): p_dict['growth_condition'] = arg
		elif opt in ("--citations"): p_dict['citations'] = arg
		elif opt in ("--data_description"): p_dict['data_description'] = arg
		elif opt in ("--transformation_description"): p_dict['transformation_description'] = arg
		elif opt in ("--biology_category_id"): p_dict['biology_category_id'] = arg
		else:
			print "Unkown option:", opt
			print __doc__
			sys.exit(2)
	return p_dict, args

if __name__ == '__main__':
	p_dict, args = parse_parameters()
	print p_dict, args
	phed = pd.parse_phenotype_file(p_dict['phen_filename'], with_db_ids=False)
	phed.convert_to_averages()
	pids = p_dict['pids']
	if not pids:
		pids = phed.phen_ids
	mids = []
	if p_dict['method_ids']:
		for pid, mid in zip(pids, p_dict['method_ids']):
			data_type = 'binary' if phed.is_binary(pid) else 'quantitative'
			new_mid = phed.insert_into_db([pid], p_dict['phenotype_scoring'], p_dict['method_description'],
						p_dict['growth_condition'], p_dict['biology_category_id'],
						p_dict['citations'], p_dict['data_description'],
						p_dict['transformation_description'], mid,
						data_type, p_dict['comment'], dry_run=p_dict['dry_run'])
			mids += new_mid
	else:
		for pid in pids:
			data_type = 'binary' if phed.is_binary(pid) else 'quantitative'
			new_mid = phed.insert_into_db([pid], p_dict['phenotype_scoring'], p_dict['method_description'],
						p_dict['growth_condition'], p_dict['biology_category_id'],
						p_dict['citations'], p_dict['data_description'],
						p_dict['transformation_description'], None,
						data_type, p_dict['comment'], dry_run=p_dict['dry_run'])
			mids += new_mid
	if len(mids):
		print 'Inserted phenotypes with DB method IDs:', mids


