import pandas as pd

def get_linear_regression_matrix():
    treatment_fimo = pd.read_csv('./TFBS Regression Modeling/trmt_fimo/fimo.tsv', header=0, sep='\t', comment='#')
    narrowpeak = pd.read_csv('./TFBS Regression Modeling/dmel_s2_bampe_q01_peaks.narrowPeak', header=None, sep='\t', comment='#')

    count_dataframe = get_wide_matrix(treatment_fimo)

    sequence_name_generated = [f'{chr.strip()}:{start}-{end}' for chr, start, end in zip(narrowpeak.iloc[:,0], narrowpeak.iloc[:,1], narrowpeak.iloc[:,2])]

    narrowpeak_only_height_name = pd.DataFrame()
    narrowpeak_only_height_name['sequence_name'] = sequence_name_generated
    # Columns follow this format: https://genome.ucsc.edu/FAQ/FAQformat.html#format12
    # We will use the 4th column (0 indexed), or score.
    narrowpeak_only_height_name['peak_height'] = narrowpeak.iloc[:,4]
    narrowpeak_only_height_name.head(3)

    joined_dataframe = count_dataframe.join(narrowpeak_only_height_name.set_index('sequence_name'), on='sequence_name')
    return joined_dataframe


def get_wide_matrix(data: pd.DataFrame):
    sequences = list(data['sequence_name'].unique())
    unique_motifs = list(data['motif_alt_id'].unique())
    unique_motifs.sort(key=lambda x: int(x.split('-')[-1]))

    motif_aggregation = dict.fromkeys(sequences, None)

    for sequence in sequences:
        motif_aggregation[sequence] = dict.fromkeys(unique_motifs, 0)

    for row_index, row in data.iterrows():
        count_dictionary = motif_aggregation[row['sequence_name']]
        count_dictionary[row['motif_alt_id']] += 1

    first_key = list(motif_aggregation.keys())[0]
    column_names = list(motif_aggregation[first_key].keys())

    all_rows = []
    for key, row_dict in motif_aggregation.items():
        row_values = [key]
        for motif_name, appearances in row_dict.items():
            row_values.append(appearances)
        
        all_rows.append(row_values)

    count_dataframe = pd.DataFrame(data=all_rows, columns=['sequence_name',*column_names])
    return count_dataframe