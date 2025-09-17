from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import cdist
from sklearn.decomposition import PCA
import numpy as np
import pandas as pd
import os
from pathlib import Path
import shutil
import gffutils

ta_orientation = 'internal_3prime'
assert ta_orientation in ['internal_3prime', 'adjacent_5prime']

s3_top_dir = 's3://lp-2022-tn-seq/lp-diet/results/combined/'
s3_wigs = s3_top_dir + "wigs/"

s3_samplesheet = 's3://lp-2022-tn-seq/lp-diet/assets/diet_samplesheet.csv'
s3_tests_to_run = 's3://lp-2022-tn-seq/lp-diet/assets/diet_tests_to_run.csv'
s3_fastas_and_gffs = 's3://lp-2022-tn-seq/lp-diet/assets/diet_fastas_and_gffs.csv'
s3_design = 's3://lp-2022-tn-seq/lp-diet/assets/diet_designs.csv'

results_dir = Path('/home/bradf/projects/nf-core-tnseq/aggregate_reports/diet/')
working_dir = results_dir / "work"

for dir in [results_dir, working_dir]:
    if os.path.isdir(dir):
        shutil.rmtree(dir)
    os.mkdir(dir)


fastas_and_gffs_df = pd.read_csv(s3_fastas_and_gffs)
tests_to_run = pd.read_csv(s3_tests_to_run)
design = pd.read_csv(s3_design)
design.to_csv(os.path.join(results_dir, "design.csv"), index=False)

contigs = fastas_and_gffs_df.name
tests = tests_to_run.design_id

print("Collecting All Zinb Files")
zinb_files = []
failed_files = []

for contig in contigs:
    for test in tests:
        s3_path = f"{s3_top_dir}transit/{contig}_{test}-zinb_out.tsv"
        cmd = f"""aws s3 cp {s3_path} {working_dir} """
        result = os.system(cmd)
        if result != 0:
            print("failed to cllect zinb files")
            failed_files.append(s3_path)
        else: 
            zinb_files.append((contig, test, os.path.basename(s3_path)))

print("Aggregating Fold Change Files")
fc_columns=['contig', 'numerator_sample_name', 'denominator_sample_name', 
            'locus_tag', 'gene_name', 'number_of_tas', 'log2_fold_change', 
            'pvalue', 'adjusted_p_value', 'mean_numerator', 
            'mean_denominator', 'nz_mean_numerator', 'nz_percent_numerator',
            'nz_mean_denominator', 'nz_percent_denominator', 'status']

all_fc_files = []
empty_files = []

for (contig, test, file) in zinb_files:
    file_loc = os.path.join(working_dir, file)
    print(f"Adding in {file}")
    numerator = test.split('_to_')[0]
    denominator = test.split('_to_')[1]
    try:
        fcs = pd.read_csv(file_loc, skiprows=3, sep='\t', header=None)
    except pd.errors.EmptyDataError:
        empty_files.append((contig, test))
        continue
    fcs.columns = ['locus_tag', 'gene_name', 'number_of_tas',	'mean_numerator', 
                    'mean_denominator',	'log2_fold_change', 'nz_mean_numerator', 
                    'nz_mean_denominator', 'nz_percent_numerator', 'nz_percent_denominator',	
                    'pvalue', 'adjusted_p_value', 'status']
    fcs['contig'] = contig
    fcs['numerator_sample_name'] = numerator
    fcs['denominator_sample_name'] = denominator
    all_fc_files.append(fcs)

stacked_fc = pd.concat(all_fc_files)

stacked_fc = stacked_fc[fc_columns]

print(f"Total records in stacked file is {stacked_fc.shape[0]}")
stacked_fc.to_csv(os.path.join(results_dir, 'stacked_fc.csv'), index=False)

print("Getting Genes and Attributes")
def get_attribute(id, db, name):
    try:
        out = db[id][name][0]
    except KeyError:
        out = None
    return out

attributes_columns = ['contig', 'id', 'locus_tag',
                      'start', 'end', 'strand', 'product']

combined_attributes = []

for index, row in fastas_and_gffs_df.iterrows():
    gff = row['gff']
    contig = row['name']
    attributes = {}
    result = os.system(f""" aws s3 cp {gff} {working_dir}""")
    if result != 0:
        print("failed to collect gff file")
    local_gff = os.path.join(working_dir, os.path.basename(gff))
    print(f"Getting attributes from {local_gff}")
    db = gffutils.create_db(local_gff, 
                            dbfn=f"{contig}.db", 
                            force=True, 
                            keep_order=True,
                            merge_strategy='merge', 
                            sort_attribute_values=True)
    for feature in db.features_of_type('CDS'):
        attributes[feature.id] = {'locus_tag': get_attribute(feature.id, db, 'locus_tag'),
                                                        'start': feature.start,
                                                        'end': feature.end,
                                                        'strand': feature.strand,
                                                        'product': get_attribute(feature.id, db, 'product')}
    tmp_df = pd.DataFrame.from_dict(attributes, orient='index').reset_index()
    tmp_df.columns = ['id', 'locus_tag', 'start', 'end', 'strand', 'product']
    tmp_df['contig'] = contig
    combined_attributes.append(tmp_df)

attributes_df = pd.concat(combined_attributes)
attributes_df.to_csv(os.path.join(results_dir, 'feature_attributes.csv'), 
                     index=False)





samplesheet = pd.read_csv(s3_samplesheet)
samples = samplesheet['sample']

print("Aggregating Wig Files")
stacked_wig = pandas.DataFrame(columns=['experiment', 'sample_name', 'contig', 'ta_position', 'count'])

get_all_wigs_cmd = f"""aws s3 cp --recursive {s3_wigs} {working_dir} --exclude="*" --include="*.wig" """
result = os.system(get_all_wigs_cmd)

wig_list = []

for sample in samples:
    for contig in contigs:
        file_loc = os.path.join(working_dir)
        file = os.path.join(file_loc, f"{contig}_{sample}.wig")
        print(f"Adding in {file}")
        sample_name = sample
        tmp_df = pd.read_csv(file, skiprows=2, header=None, names=['ta_position', 'count'], sep=' ')
        print(tmp_df.shape)
        tmp_df['sample_name'] = sample_name
        tmp_df['contig'] = contig
        wig_list.append(tmp_df)


stacked_wig = pd.concat(wig_list)

print(f"Total records in wig file is {stacked_wig.shape[0]}")
stacked_wig.to_csv(os.path.join(results_dir, 'stacked_wig.csv'), index=False)
# samplesheet_df.to_csv(os.path.join(results_dir, 'samplesheet.csv')

gene_coord_list = []
for i, row in attributes_df.iterrows():
    coord = range(row.start, row.end+1)
    contig = row.contig
    tmp_df = pd.DataFrame({'contig': contig,
                               'id': row.id,
                               'ta_position': coord
                               })
    gene_coord_list.append(tmp_df)

gene_coord = pd.concat(gene_coord_list)

gene_coord = gene_coord.merge(stacked_wig, how='inner', on=['contig', 'ta_position'])
gene_coord.to_csv(os.path.join(results_dir, 'gene_ta_positions.csv'), index=False)

pca_data = gene_coord.groupby(["contig", "sample_name", "id"])['count'].sum().reset_index()
pca_data.columns = ['contig', 'sample_name', 'feature_id', 'count']
pca_data_pivot = pca_data.pivot(index=['contig', 'feature_id'], 
                              values=['count'], 
                              columns='sample_name').reset_index()
pca_data_pivot.columns = ['_'.join(col).strip().replace(' ', '_').strip("_").replace("count_", "") for col in pca_data_pivot.columns.values]

all_pca_list = []

# for contig in contigs:
# pca_data_sub = pca_data_pivot.loc[pca_data_pivot.contig==contig,:]
pca_data_sub = pca_data_sub.set_index(["feature_id"]).drop(columns=['contig'])
row_index = pca_data_sub.index
col_index = pca_data_sub.columns
pca_array = np.array(pca_data_sub, dtype='f')
pca_array = np.log(pca_array+1)
pca_array = np.nan_to_num(pca_array, nan=0, neginf=-10, posinf=10)
# features in columns, samples in rows
pca_array = pca_array.transpose()
scaler = StandardScaler().fit(pca_array)
pca_input = scaler.transform(pca_array)
pca = PCA()
pca_output = pca.fit(X=pca_input)
pca_transformed_data = pca_output.transform(pca_input)
pca_with_index = pd.DataFrame(data=pca_transformed_data, index=col_index)
pca_with_index.columns = [f"pc_{c+1}" for c in pca_with_index.columns]
# pca_with_index['contig'] = contig
pca_with_index['sample_name'] = pca_with_index.index
all_pca_list.append(pca_with_index)

all_pca_df = pd.concat(all_pca_list)
all_pca_df.to_csv(os.path.join(results_dir, 'pca.csv'), index=False)