MAIC_inputcheck.py

This script collates a set of MAIC input lists (within a single folder), checks the formatting, harmonizes gene names/symbols,
and converts to a single MAIC input files. It can also be used to convert gene names in .gmt files
(inputs for gene set enrichment analysis).

Requirements: Python 3, package pickle. Gene name alias file.

Usage:
MAIC_inputcheck.py <-d inputdir> <-a aliasfile> <-s suffix> <o outputfile> <-gmt>

or MAIC_inputcheck_multispecies.py <-d inputdir> <-a aliasfile> <-s suffix> <o outputfile> <-gmt> <-x> <-sp>


Optional arguments:
 -d input directory - directory containing gene lists. All lists within a single directory will be processed.
                      Default: current working directory
 -a alias file - .pkl file containing json object with gene alias definitions. Default: gene_name_aliases.pkl
 -s suffix - suffix identifying gene lists to process. Default: ".txt"
 -o output filename. Default: MAIC_input_{date}.txt"
 -gmt  Convert .gmt files rather than MAIC input
 -x  Remove unconverted genes from output list
 -sp  (for multispecies version only) - this option indicates that species will be specified on the 5th line of each list (after the 'NAMED_GENES' line)
 

NOTES ON INPUTS:
----------------
Inputs should be simple gene lists in the standard format for MAIC, i.e. .txt file with one gene per line, 
with the first 4 lines containing the headers: list type, name, RANKED/UNRANKED, "NAMED_GENES". They should
all be identified by a single common suffix which distinguishes them from other files in the target directory.
The script will fix some common errors/typos in the header, including variations in category descriptions, e.g.
"transcriptomic" vs "transcriptomics".

Caution if you are re-running the script in a folder after an initial attempt - make sure the suffixes of your
input files are distinct from previous output files (or delete the old ones!), otherwise the script will try and
include these as new gene lists.

If specifying species in line 5 of the input, this must (currently) be one of the following, exactly as given here:
human
mouse
rat
ferret
macaque
African green monkey
bovine
pig
chicken
horse
dog
cat
sheep
goat


NOTES ON GENE NAME CONVERSION:
------------------------------
The script attempts to convert gene IDs from a variety of sources to the associated HGNC symbol where possible.
If there is no HGNC equivalent, GenCode symbols will be used as a second choice, or failing that Refseq / other 
original annotation (e.g. Ensembl gene ID).

This relies on a dictionary of gene aliases (gene_name_aliases.pkl, last updated Jan 2021).
The following annotation systems should be supported:
    HGNC symbol
    HGNC approved name (if spelled exactly as source!!)
    HGNC number as 'HGNC:1234'
    Old HGNC symbols
    Homologues from Homologene(mouse; ferret; macaque; African green monkey. Bovine, pig and chicken
	   added in later version of file gene_name_aliases_v2species.pkl)
    GenCode gene names (including lncRNAs etc not included in HGNC)
    Ensembl gene/transcript ID (ENSG / ENST) - currently supported for human, mouse, rat,
      cat and horse, but can add others as needed in the future.
    ucsc ID (e.g. uc123abc.1)
    uniprot ID (e.g Q9UNA3)
    gi numbers formatted as gi|123456|
    Vega / havana  e.g. OTTHUMG / OTTHUMT...
    refseq ID beginning with NM_, XR_ LOC... etc.
    Miscellaneous symbols including old symbols, Cosmic symbols etc 

In case the dictionary is to be used in other scripts, it is in the following format:
  {"Gene alias":
       {"source of annotation" : "Gene symbol to return (HGNC if possible)"},
   }

Before conversion, the script will fix common errors introduced by Excel date formatting,
and remove unecessary suffixes which can make matching difficult (e.g. ".1" on the end of a transcript ID).

Where no match can be found for a gene ID, it will be left in its original form in the output lists.

Where more than one entry in the lists maps to the same gene (expected behaviour for transcriptomic data,
as there may be more than one transcript per gene), the first entry will be kept and subsequent entries discarded
to avoid duplication.

Where a gene ID in the list is ambiguous (i.e. could map to more than one gene depending on which annotation system
the input is using), this is handled as follows:
   1. If species is specified as human, remove all non-human homologue options
   2. If a non-human species is specified, select the appropriate species-specific homologue
   3. If this does not resolve ambiguity, or species not specified, exclude any options already in gene list
   2. If still >1 option, choose 'hgnc' if this is an option
   3. If not, attempt to infer the annotation system of the input, by choosing the most popular of the available
      annotation schemes among non-ambiguous genes 
   4. If still can't resolve, arbitrarily assign to first in list
   
   
NOTES ON OUTPUTS:
------------------
The main output is a single MAIC input file as a tab-delimited text file with one gene list per line.

The script will also print a list of the experiment categories detected, so this can be checked against expectations
and any undetected typos removed to avoid unintended assignment as separate categories.

Additional outputs:
  [].txt.duplicates file lists the genes which had more than one entry for each gene list.
     (Expected to have multiple for transcriptomics etc, but wouldn't be expected for some other data types)
  [].txt.unconverted file lists IDs for each list for which no match could be found.
     **** IMPORTANT TO CHECK THIS MANUALLY BEFORE RUNNING MAIC **** - It will pick up typos / other errors
	either in the original data or during gene list extraction - these can be corrected in the original lists and then 
	the script can be re-run. It will also flag any gene lists which use an annotation system not supported by the alias
	definition file. Note that it is quite common for deprecated IDs to be picked up here, sometimes in large numbers
	(esp. Uniprot IDs and putative lncRNAs) - default behaviour is for these to be left in the list in their original form,
	but the researcher may also choose to remove them.



 