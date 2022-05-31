import pandas, requests, itertools

mapping = pandas.read_csv('./mapping.csv')

def infer_properties(row):
    cols = row.local_sample_name.split('-')

    assert cols[2] == 'cBA.1'

    n_snps = int(cols[4][0])


    if len(cols) == 9:
        dropped_amplicon = True
    else:
        dropped_amplicon = False

    return pandas.Series([n_snps, dropped_amplicon])


mapping[['N_SNPS','DROPPED_AMPLICON']] = mapping.apply(infer_properties, axis=1)

samples = list(mapping.gpas_sample_name)

table=[]
for idx,guid0 in enumerate(samples):
    a = mapping.loc[mapping.gpas_sample_name == guid0]
    for guid1 in samples[idx+1:]:
        b = mapping.loc[mapping.gpas_sample_name == guid1]
        n_snps = abs(int(a.N_SNPS)+int(b.N_SNPS))
        if n_snps > 3:
            present = False
        else:
            present = True
        table.append([guid0, guid1, present, n_snps])
api_calls = pandas.DataFrame(table, columns=['guid0','guid1','present','snp_distance'])        

api_start = 'http://fn4-gpasdev-app.privsn001.zymkc.oraclevcn.com:5092/api/v2/'
api_end = '/neighbours_within/3/in_format/'

def make_api_calls(row):

    passes=False

    api_url = api_start + row.guid0 + api_end + '3'
    a = requests.get(api_url)
    if a.ok:
        if row.guid1 not in a.json():
            if not row.present:
                passes=True
        else:
            api_url = api_start + row.guid0 + api_end + '1'
            a = requests.get(api_url)
            if a.ok:
                for (guid1, distance) in a.json():
                    if row.present:
                        if row.guid1 == guid1 and row.snp_distance == distance:
                            passes = True
                
                
    else:
        print("failed")
        
    return(passes)         

api_calls['passes'] = api_calls.apply(make_api_calls, axis=1)

print(api_calls)

api_calls.to_csv('test-fn4-results.csv')

print((~api_calls.passes).sum()==0)
