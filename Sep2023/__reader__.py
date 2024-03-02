import pandas as pd

def merge_phys_databases(source_path, input_path, output_path):

    SourceDataFrame = pd.read_csv(source_path, sep="\s{2,}", header=0, index_col=0, engine='python')
    InputDataFrame = pd.read_csv(input_path)
    MergedDataFrame = pd.merge(InputDataFrame, SourceDataFrame, how='inner', on='CATAID_1')
    print(MergedDataFrame.shape)
    MergedDataFrame.to_csv(output_path, index=False)
    
def merge_phys_databases_outflow(source_path, input_path, output_path):

    SourceDataFrame = pd.read_csv(source_path, sep="\s{2,}", header=0, index_col=0, engine='python', usecols=['SPECID', 'CATAID_1', 'Z_1', 'SFR_0_1Gyr_percentile16', 'SFR_0_1Gyr_percentile50', 'SFR_0_1Gyr_percentile84', 'mass_stellar_percentile16', 'mass_stellar_percentile50', 'mass_stellar_percentile84', 'ager_percentile50'])
    InputDataFrame = pd.read_csv(input_path)
    MergedDataFrame = pd.merge(InputDataFrame, SourceDataFrame, how='inner', on='SPECID')
    print(MergedDataFrame.shape)
    MergedDataFrame.to_csv(output_path, index=False)

def merge_phys_databases_all(source_path, input_path, output_path):

    SourceDataFrame = pd.read_csv(source_path, sep="\s{2,}", header=0, index_col=0, engine='python')
    InputDataFrame = pd.read_csv(input_path)
    MergedDataFrame = pd.merge(InputDataFrame, SourceDataFrame, how='inner', on='SPECID')
    print(MergedDataFrame.shape)
    MergedDataFrame.to_csv(output_path, index=False)
    
def calculating_mdms(output_path):
    DataFrame = pd.read_csv(output_path)
    MD_down = DataFrame['mass_dust_percentile16']
    MD = DataFrame['mass_dust_percentile50'] 
    MD_up = DataFrame['mass_dust_percentile84']
    
    MS_down = DataFrame['mass_stellar_percentile16']
    MS = DataFrame['mass_stellar_percentile50'] 
    MS_up = DataFrame['mass_stellar_percentile84']
    
    MDMS = MD - MS
    MDMS_up = MDMS + ((MD_up - MD)**2 + (MS - MS_down)**2)**0.5
    MDMS_down = MDMS - ((MD - MD_down)**2 + (MS_up - MS)**2)**0.5
    
    DataFrame['mdms_percentile50'] = MDMS
    DataFrame['mdms_percentile16'] = MDMS_down
    DataFrame['mdms_percentile84'] = MDMS_up
    
    DataFrame.to_csv(output_path, index=False)
    

