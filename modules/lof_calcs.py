import pandas as pd
import numpy as np

def make_dataframe(file):
    # Takes a .tsv, makes it binary, outputs a pandas dataframe
    # Input: sv_caller.py output to be converted to binary
    # Output: A binary dataframe with data from a population.
    
    #Read the input file into a pandas dataframe
    df_in = pd.read_csv(file, sep = '\t', index_col=0)
    
    #for each key and values in the dictionary of dataframes
    #make it binary, create new dataframe in clean dataframe dictionary
    #with the resulting numpy array and the old indices and columns
    
    df_array = np.where(df_in < 0.5,0,1)
    final_df = (pd.DataFrame(df_array, index = df_in.index, columns=df_in.columns))

    #final_df_dict.to_csv(file[:-4]+"_binary.tsv")
    #return merged dataframe, population sizes
    return final_df

def intra_pop_stat(df):
    #Calculates the percentage of samples for which a gene is mutated
    #by summing the binary row and dividing by the number of original columns
    #Input: binary pandas dataframe from make_dataframe.py
    #Output: pandas dataframe with percentage of samples each gene was mutated in 
    
    import pandas as pd
    import numpy as np
    
    #Create an empty pandas dataframe keeping the index, i.e. genes
    intra_pop_stats = (pd.DataFrame(index = df.index))

    #since it's a binary df, sum of rows over number of columns i.e. samples will give percentage of samples genes were mutated in
    intra_pop_stats["pc_mutated_spls"] = df.sum(axis = 1) / len(df.columns)
    
    return intra_pop_stats

def inter_pop_stat(df1, df2):
    #Calculates the absolute difference between percentages of mutated samples for each gene for two populations
    #Input: two pandas dataframes with percentage of samples each gene was mutated in 
    #Output: pandas dataframe with percentage differences for the two populations
    
    import pandas as pd
    import numpy as np

    inter_difs_df = (pd.DataFrame(index = df1.index))
    
    inter_difs_df["abs_diffs_pc_mut_spls"] = (df1["pc_mutated_spls"] - df2["pc_mutated_spls"]).abs()
    
    return inter_difs_df

def product_calc(a, b, c):
    #Calculates the absolute difference between percentages of mutated samples for each gene for two populations
    #Input: two pandas dataframes with percentage of samples each gene was mutated in 
    #Output: pandas dataframe with percentage differences for the two populations
    
    import pandas as pd
    import numpy as np

    inter_difs_df = (pd.DataFrame(index = a.index))

    inter_difs_df["abs_diffs_pc_mut_spls"] = (a["abs_diffs_pc_mut_spls"] - b["abs_diffs_pc_mut_spls"] + c["abs_diffs_pc_mut_spls"])
    
    return inter_difs_df

def add_gff_metadata(df, gff):
    #(adapted from lof_utility.py) add annotations from gff to LOF occurrence calculations per gene
    #Input: one pandas dataframe, handled by inter_pop_stat_2(), one GFF whose genes correspond to indices in df
    #Output: a pandas dataframe annotated with each gene's start and end, as well as the corresponding chromosome
    
    with open(gff, "r") as gff_file:
        for line in gff_file:
            if line.startswith(("#", "A", "C", "T", "G", ">")):
                continue
            ann = line.split("\t")[8]
            if len(ann.split(";")) > 1:
                id = ann.split(";")[1]
                gene = id.replace('geneID=','')
                chrom = line.split("\t")[0]
                start = line.split("\t")[3]
                end = line.split("\t")[4]
            
                df.loc[gene, 'chromosome'] = chrom
                df.loc[gene, 'start'] = start
                df.loc[gene, 'end'] = end
            
    return df

def melt_merge_pops(df1, df2):
    #Melt and merge two dataframes vertically.
    #Input: two pandas dataframes out of add_gff_metadata()
    #Output: a pandas dataframe, product of the input two, ready for all kinds of plotting
    
    df1['zigosity'] = "heterozyguous"
    df2['zigosity'] = "homozyguous"

    temp_df1 = pd.melt(df1, id_vars=['chromosome', 'start', 'end', 'zigosity'], 
                        value_vars=['abs_diffs_pc_mut_spls'], value_name = 'abs_diffs_pc_mut_spls',
                        ignore_index=False).drop(['variable'], axis=1)

    temp_df2 = pd.melt(df2, id_vars=['chromosome', 'start', 'end', 'zigosity'], 
                        value_vars=['abs_diffs_pc_mut_spls'], value_name = 'abs_diffs_pc_mut_spls',
                        ignore_index=False).drop(['variable'], axis=1)

    temp_df1['gene'] = temp_df1.index
    temp_df2['gene'] = temp_df2.index
    temp_df1.reset_index(drop=True, inplace=True)
    temp_df2.reset_index(drop=True, inplace=True)

    merged_temp = temp_df1.append(temp_df2, ignore_index=True)

    return merged_temp

def count_lofs(df):
    #counts numbers of times a gene has had LOF mutations detected per sample
    #Input: dataframe of LOF mutations per sample (columns) and per gene (rows, index)
    #Output: dataframe with total numbers of LOF counts per sample (rows) per bin (column)
    
    #create dummy dataframe from columns in input df, skipping 'genes' column
    out_df = pd.DataFrame()
    out_df['samples'] = df.columns
    out_df = out_df[~out_df.samples.str.contains("genes")]

    #create columns for output
    out_df['LOF_count_0'] = 0
    out_df['LOF_count_1'] = 0
    out_df['LOF_count_2'] = 0
    out_df['LOF_count_3'] = 0
    out_df['LOF_count_4'] = 0
    out_df['LOF_count_5'] = 0
    out_df['LOF_count_6_10'] = 0
    out_df['LOF_count_10_plus'] = 0

    #loop through all samples
    for sample in df.columns:
        if sample not in 'genes':
    #loop through all genes, after skipping 'genes' column
            for idx in df.genes:
                #if the value in the input df for that sample and gene is 0
                if df[sample][idx] == 0:
                #add one to output df, in column "LOF_count_0", same row as sample being analysed
                    out_df.loc[out_df.samples == sample, "LOF_count_0"] += 1
                elif df[sample][idx] == 1:
                    out_df.loc[out_df.samples == sample, "LOF_count_1"] += 1
                elif df[sample][idx] == 2:
                    out_df.loc[out_df.samples == sample, "LOF_count_2"] += 1
                elif df[sample][idx] == 3:
                    out_df.loc[out_df.samples == sample, "LOF_count_3"] += 1
                elif df[sample][idx] == 4:
                    out_df.loc[out_df.samples == sample, "LOF_count_4"] += 1
                elif df[sample][idx] == 5:
                    out_df.loc[out_df.samples == sample, "LOF_count_5"] += 1
                elif (df[sample][idx] > 5 & df[sample][idx] <= 10):
                    out_df.loc[out_df.samples == sample, "LOF_count_6_10"] += 1
                elif df[sample][idx] > 10:
                    print(df[sample][idx])
                    #out_df.loc[out_df.samples == sample, "LOF_count_10_plus"] += 1
                else:
                    break
        else:
            continue
    return out_df
