#!/usr/bin/python3
#
#
#
#
# AK (December), contact: AK@addenbrookes.nhs.uk

'''
## SM comment:
# overall concise script with good comments throughout
# good to include README.md with a required module and how to use the script with usage
# object and variable names are meaningful for what they intend to do/be

# There are several points which can improve this script:
# mentioning which python3 version can be used for development and to use
# the general description and the purpose of the script in itself with a version or anything to keep tracking edit history can be helpful for someone who maintain it in future
# specification of an input file format would be helpful for troubleshooting. This can be included in README.md document
#
# error handling to consider - average() function may give error if a length (len() may be zero
# incorrect indentation and blank lines without indentation give an error so use editor functions to visualise whitespace and indentation may be helpful to identify this type of errors (parse_arguments() and per_gene_coverage() have blank lines which cause errors)
# comments on arguments what type they are are helpful for maintenance and troubleshooting
'''


import pandas as pd
import argparse
from datetime import date

'''
SM comment on the function parse_arguments():
blank lines (I know it is hard to tell whether indentation exists or not) would give errors when running. Use of functions in code editors to show whitespace and indentation may be helpful


'''
def parse_arguments():

    parser = argparse.ArgumentParser(description = 'Generate a report containing genes with suboptimal coverage')
    parser.add_argument('sambamba_output',type=str,help='Provide path to the sambamba output to be interrogated')
    parser.add_argument('--outfile',type=str,help='Optional filename for generated report')

    args = parser.parse_args()

    return args

'''
SM comment: 
good example of reusable function with comment which can be used within other functions

suggest to add an error handling code (try...except) for dividedbyzero error
'''
def average(lst):
    """ Calculates the mean of a list of values """
    return sum(lst)/len(lst)


def per_gene_coverage(genes,df):
    """Given a list of gene, it aggregates the exon converage and returns the genes with sub-optimal coverage """

    sub_genes =[] # SM: empty line without indentation gives an error

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

