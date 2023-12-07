import pandas as pd

def merge_phys_databases(source_path, input_path, output_path):

    SourceDataFrame = pd.read_csv(source_path, sep="\s{2,}", header=0, index_col=0, engine='python')
    InputDataFrame = pd.read_csv(input_path)
    MergedDataFrame = pd.merge(InputDataFrame, SourceDataFrame, how='inner', on='CATAID_1')
    print(MergedDataFrame.shape)
    MergedDataFrame.to_csv(output_path, index=False)

