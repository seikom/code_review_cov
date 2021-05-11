#!/usr/bin/python3
#
#
#
#
# AK (December), contact: AK@addenbrookes.nhs.uk

import pandas as pd
import argparse
from datetime import date

def parse_arguments():
    
    parser = argparse.ArgumentParser(description = 'Generate a report containing genes with suboptimal coverage')
    parser.add_argument('sambamba_output',type=str,help='Provide path to the sambamba output to be interrogated')
    parser.add_argument('--outfile',type=str,help='Optional filename for generated report')

    args = parser.parse_args()

    return args

def average(lst):
    """ Calculates the mean of a list of values """
    return sum(lst)/len(lst)


def per_gene_coverage(genes,df):
    """Given a list of gene, it aggregates the exon converage and returns the genes with sub-optimal coverage """

    sub_genes =[]

    #For every gene in the list, check the average coverage, if less than 100 add it to the final list.
    for gene in genes:
        coverage = average(df[df['GeneSymbol;Accession'] == gene]['percentage30'])

        if coverage < 100:
          sub_genes.append([gene.split(';')[0],round(coverage,2)])
    
    return sub_genes
        

def main(args):
    
    input = args.sambamba_output
    today = date.today()
    
    # Header had a mix of tabs and spaces, so assigned it manually assuming the sambamba_output would be always the same.
    correct_columns = ['chromosome','StartPosition','EndPosition','FullPosition', 'NotUsed','NotUsed2',
        'GeneSymbol;Accession','Size','readCount','meanCoverage','percentage30','sampleName']
    # Read the data in a dataframe, with the above column names and omitting the original ones.
    df = pd.read_csv(input, '\t', names = correct_columns, skiprows=1)

    # Check if user has provided a prefered output name.
    if args.outfile:
        output = args.outfile
    else:
        output = input.rstrip('sambamba_output.txt') + "_suboptimal_genes.txt"
    
    print ("Output file generated ", output)

    # Get the list of genes/transcripts 
    gene_names = list(set(df['GeneSymbol;Accession']))
    result = per_gene_coverage(gene_names,df)

    # Write the output in a .txt file
    out= open(output,'w')
    out.write("Input file:  {}\n".format(input))
    out.write("Date analysed:{}\n\n".format(today))
    out.write("Genes with less than 100% coverage at 30x: \n")
    
    for gene in result:
        out.write("{}\t {}%\n".format(gene[0],gene[1]))
    
    out.close()


if __name__ == '__main__':
    arguments = parse_arguments()
    main(arguments)

