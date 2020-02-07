
"""
1000 crosses GWAS
"""
import sys
sys.path.append('./..')
import env
import dataParsers as dp
import scipy as sp
import bisect
import snpsdata
import phenotypeData as pd



def create_diploid_dataset(call_method_id=76, file_name='/tmp/test.csv', coding_type='normal'):

    #Load parent list
    parents = []
    with open(env.env['data_dir'] + 'heterozygous_genotypes.csv') as f:
        f.next()
        for l in f:
            parents.append(map(str.strip, l.split(',')))

    snpsd = dp.load_snps_call_method(call_method_id)
    l = zip(snpsd.accessions, range(len(snpsd.accessions)))
    l.sort()
    l = map(list, zip(*l))
    acc_list = l[0]
    orders = l[1]
    sds = []
    for i, sd in enumerate(snpsd.snpsDataList):
        snps = sp.array(sd.snps, dtype='int8')
        snps_list = []
        p_list = []
        for ps in parents:
            f_id = bisect.bisect(acc_list, ps[0]) - 1
            m_id = bisect.bisect(acc_list, ps[1]) - 1
            if acc_list[f_id] == ps[0] and acc_list[m_id] == ps[1]:
                f_gt = snps[:, orders[f_id]].flatten()
                m_gt = snps[:, orders[m_id]].flatten()
                if coding_type == 'normal':
                    o_gt = f_gt + m_gt
                elif coding_type == 'dominant':
                    o_gt = sp.bitwise_xor(f_gt, m_gt)
                snps_list.append(o_gt)
                p_list.append('%s_%s' % (ps[0], ps[1]))
        snps_list = sp.transpose(sp.array(snps_list, dtype='int8'))
        snps = []
        for s in snps_list:
            snps.append(s)
        sds.append(snpsdata.SNPsData(snps, sd.positions, accessions=p_list, chromosome=i + 1))
    sd = snpsdata.SNPsDataSet(sds, [1, 2, 3, 4, 5])
    sd.writeToFile(file_name)




def parse_phenotype_data(file_name='/tmp/test_phen.csv'):
    phen_dict = {}
    for pid in range(1, 7):
        phen_dict[pid] = {'name':'criteria_%d' % pid, 'ecotypes':[], 'values':[],
                          'transformation':None}
    et_set = set()
    with open(env.env['data_dir'] + 'phenot_HybridIncompatibilities.csv') as f:
        f.next()
        for line in f:
            l = map(str.strip, line.split(','))
            et = '%s_%s' % (l[0], l[1])
            if not et in et_set:
                et_set.add(et)
                for pid in range(1, 7):
                    val = l[pid + 1]
                    if val != 'NA':
                        phen_dict[pid]['values'].append(int(val))
                        phen_dict[pid]['ecotypes'].append(et)
    for pid in range(1, 7):
        phen_dict[pid]['raw_values'] = phen_dict[pid]['values']
    phend = pd.phenotype_data(phen_dict=phen_dict)
    phend.write_to_file(file_name)


if __name__ == "__main__":
    #parse_phenotype_data()
    #create_diploid_dataset(call_method_id=76, file_name='/tmp/cm_76_crosses_additive.csv')
    create_diploid_dataset(call_method_id=76, file_name='/tmp/cm_76_crosses_dominant.csv', coding_type='dominant')




