import pandas as pd
import numpy as np

def make_dataframe(file1, file2):
    # Input: Two filenames, each containing LOF data
    #        from populations that should be compared
    # Output: A dataframe where data from the two populations 
    #         have been joinde. First set of columns pertains 
    #         to population 1 and the rest to population 2. 

    # Data from each population is stored in separete input files
    # Note: The column with gene name is used an index (e.g. row name) 
    # in the resulting data frame
    df1 = pd.read_csv(file1, sep = '\t', index_col=0)
    
    # Make mutation data binary (0 = not mutated, 1 = mutated)
    df1[:] = np.where(df1 < 0.5,0,1)
    
    df2 = pd.read_csv(file2, sep = '\t', index_col=0)
    df2[:] = np.where(df2 < 0.5,0,1)
    
    # Count the number of samples in each population
    pop_size_1 = len(df1.columns)
    pop_size_2 = len(df2.columns)
    
    # Join the population data in one data frame
    df = pd.merge(df1, df2, right_index=True, left_index=True)
    
    return df, pop_size_1, pop_size_2



def prep_data(df, pop1, pop2):
    # Input: Data frame with joined data from the two populations
    #        Integer indicating the number of columns in the data 
    #            frame pertaining to population 1
    #        Integer indicating the number of columns in the data 
    #            frame pertaining to population 2
    # Output: Dictionary with population statistics in the form of boolean numpy arrays.
    #         Dictionary with lists of gene names
    
    # The resulting statistics is appended to the input data frame. 
    # However, selections of the data frame, given specif criteria 
    # (e.g "which genes in population 1 is not mutated") are stored 
    # as boolean numpy arrays. These arrays are in turned stored in 
    # a list which is returned from this function. 
    statistics = {}
    gene_names = {}

    # Summarise the number of identified LOF mutations in each gene per pupulation
    df["mutated_pop1"] = df.iloc[:,:pop1].sum(axis = 1)
    df["mutated_pop2"] = df.iloc[:,pop1:pop1 + pop2].sum(axis = 1)
    df["mutated_in_both"] = df["mutated_pop1"] + df["mutated_pop2"]

    # Ration of genes (mutated genes / total number of samples) in each population
    df["ratio_mutated_pop1"] = df["mutated_pop1"] / pop1
    df["ratio_mutated_pop2"] = df["mutated_pop2"] / pop2
    df["mutation_ratio_difference"] = (df["ratio_mutated_pop1"] - df["ratio_mutated_pop2"]).abs()
    
    # No LOF mutations in population 1
    not_mutated_pop1 = df["mutated_pop1"] == 0
    statistics["not_mutated_pop1"] = not_mutated_pop1

    # Muteted genes in population 1
    mutated_pop1 = df["mutated_pop1"] != 0
    statistics["mutated_pop1"] = mutated_pop1

    # No LOF mutations in population 2
    not_mutated_pop2 = df["mutated_pop2"] == 0
    statistics["not_mutated_pop2"] = not_mutated_pop2
    
    # Muteted genes in population 2
    mutated_pop2 = df["mutated_pop2"] != 0
    statistics["mutated_pop2"] = mutated_pop2

    # At least one sample in one of the populations has a LOF mutation in these genes
    mutated = mutated_pop1 | mutated_pop2
    statistics["mutated"] = mutated
    
    # No sample in either of the populations has a LOF mutation in these genes
    never_mutated = df["mutated_in_both"] == 0
    statistics["never_mutated"] = never_mutated
    
    # Extract the gene names in the different categories
    genes_never_mutated = list(df[never_mutated].index.values)
    gene_names["genes_never_mutated"] = genes_never_mutated

    # Genes sometimes or always mutated
    genes_mutated = list(df[mutated].index.values)
    gene_names["genes_mutated"] = genes_mutated

    # Genes mutated in population 2 but not in population 1
    genes_not_mutated_pop1 = list(df[not_mutated_pop1].index.values)
    genes_mutated_pop2 = list(df[mutated_pop2].index.values)
    unique_not_mutated_genes_pop1 = []
    for gene in genes_not_mutated_pop1:
        if gene in genes_mutated_pop2:
            unique_not_mutated_genes_pop1.append(gene)
    gene_names["unique_not_mutated_genes_pop1"] = unique_not_mutated_genes_pop1

    # Genes mutated in population 1 but not in population 2
    genes_not_mutated_pop2 = list(df[not_mutated_pop2].index.values)
    genes_mutated_pop1 = list(df[mutated_pop1].index.values)
    unique_not_mutated_genes_pop2 = []
    for gene in genes_not_mutated_pop2:
        if gene in genes_mutated_pop1:
            unique_not_mutated_genes_pop2.append(gene)
    gene_names["unique_not_mutated_genes_pop2"] = unique_not_mutated_genes_pop2
        
    # Genes with a difference in gene mutation ration >= 0.25 between populations
    ratio_025 = df["mutation_ratio_difference"] >= 0.25
    genes_ration_025 = list(df[ratio_025].index.values)
    gene_names["genes_ration_025"] = genes_ration_025
        
    # Genes with a difference in gene mutation ration >= 0.5 between populations
    ratio_05 = df["mutation_ratio_difference"] >= 0.5
    genes_ration_05 = list(df[ratio_05].index.values)
    gene_names["genes_ration_05"] = genes_ration_05

    # Genes with a difference in gene mutation ration >= 0.75 between populations
    ratio_075 = df["mutation_ratio_difference"] >= 0.75
    genes_ration_075 = list(df[ratio_075].index.values)
    gene_names["genes_ration_075"] = genes_ration_075

    # Genes with a difference in gene mutation ration >= 0.9 between populations
    ratio_09 = df["mutation_ratio_difference"] >= 0.9
    genes_ration_09 = list(df[ratio_09].index.values)
    gene_names["genes_ration_09"] = genes_ration_09
    
    return statistics, gene_names



def general_stats(df, statistics, gene_names):
    print(f"Number of genemodels: {len(df)}")
    print(f"Number of genes with mutation: {len(df[statistics['mutated']])}")
    print(f"Number of genes never mutated: {len(df[statistics['never_mutated']])}")
    print(f"Number of unique LOF genes population 1: {len(gene_names['unique_not_mutated_genes_pop1'])}")
    print(f"Number of unique LOF genes population 2: {len(gene_names['unique_not_mutated_genes_pop2'])}")
    print(f"Number of genes with mutation ration >= 0.25: {len(gene_names['genes_ration_025'])}")
    print(f"Number of genes with mutation ration >= 0.5: {len(gene_names['genes_ration_05'])}")
    print(f"Number of genes with mutation ration >= 0.75: {len(gene_names['genes_ration_075'])}")
    print(f"Number of genes with mutation ration >= 0.9: {len(gene_names['genes_ration_09'])}")
