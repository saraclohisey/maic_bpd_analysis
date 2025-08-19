# PROSPERO Protocol

1. Review title

    Prioritising host genes implicated in Bronchopulmonary Dysplasia: a systematic review and meta-analysis by information content of omic studies.

2. Original language title

    NA

3. Anticipated or actual start date

    01/02/2024

4. Anticipated completion date

    31/01/2025

[comment]: # always give yourself the most time, easier to finish early than amend if you're late

5. Stage of review at time of submission

[comment]: # all ok

6. Named contact

    Sara Clohisey Hendry

7. Named contact email

    sara.clohisey@ed.ac.uk

8. Named contact address

    Baillie Lab, Roslin Institute, Easter Bush Campus, Midlothian, United Kingdom, EH25 9RG

[comment]: # I wouldnt supply any personal information as it'll be visible on the internet, same goes for your phone number.

9. Named contact phone number

    0131 651 9100

10. Organisational affiliation of the review

    The University of Edinburgh

11. Review team members and their organisational affiliations

    Dr Sara Clohisey Hendry, The University of Edinburgh <sara.clohisey@roslin.ed.ac.uk>

    Dr Jonathan Millar, The University of Edinburgh <jonathan.millar@ed.ac.uk>
    
    Prof J Kenneth Baillie, The University of Edinburgh, <j.k.baillie@ed.ac.uk> 

    TBD

12. Funding sources/sponsors

    None

13. Conflicts of interest

    None

14. Collaborators

    None

15. Review question

    "In premature babies with Bronchopulmonary Dysplasia, which genes are prioritised as being implicated in the host-response?"

16. Searches

    Primary search

    OVID MEDLINE and Epub Ahead of Print, In-Process, In-Data Review & Other Non-Indexed Citations, Daily and Versions 1946 to January 31, 2024.

    AND

    Embase Classic + Embase 1947 to January 31, 2024.

    Secondary search

    - Search of [bioRxiv](https://www.biorxiv.org) and [medRiv](https://www.medrxiv.org) using relevant keywords.
    - Single-level forward and backward citation searches of the reference lists of included articles.
    - Extraction of references contained in the [ARDS DB](https://doi.org/10.3389/fgene.2021.750568).
    - Correspondence with corresponding authors of included studies.

    No language restrictions will be imposed.

    Inclusion will be restricted to studies published after 1967.

    The full search strategy is described at the link below.

17. URL to search strategy

    [Full search strategy]https://www.dropbox.com/s/jdffhs3lokw724t/SearchStrategy_20220120.txt?dl=0

[comment]: # I'll add an updated file to the dropbox for this purpose

[comment]: # filepath

18. Condition or domain being studied

    Bronchopulmonary Dysplasia.

19. Participants/population

    Inclusion:

    - Human studies: *in-vivo* or *in-vitro*
    - Premature humans (<29 weeks PMA)
    - Accepted methodologies:
      - CRISPR screen
      - RNAi screen
      - Protein-protein interaction study
      - Host proteins incorporated into virion or virus-like particle
      - Genome wide association study
      - Transcriptomic study
      - Proteomic study

    Exclusion:

    - Term babies, Children, Adults
    - Animal studies
    - Meta-analyses, *in-silico* analyses, or re-analysis of previously published data
    - Excluded methodologies:
      - Candidate *in-vivo* or *in-vitro* transcriptomic or proteomic studies (defined as those investigating < 50 genes)
      - Candidate gene association studies
      - Studies including fewer than 5 individuals in either the control or BPD arm

20. Intervention(s), exposure(s)

    NA

21. Comparator(s)/control

    NA

22. Types of study to be included

    We are seeking to include whole-genome studies which report an association between genes, transcripts, or proteins and susceptibility to BPD or with severity or outcome. 

    We will include studies employing the following methodologies: CRISPR screen, RNAi screen, protein-protein interaction, host proteins incorporated into virion or virus-like particle, genome-wide association, transcriptomic study, or proteomic study.

    We will exclude studies in which only results from non-human cells or animals are reported.

    We will exclude studies employing the following methodologies: candidate transcriptomic or proteomic studies (< 50 genes investigated) and candiate gene association studies.

    We will exclude studies including fewer than 5 patients in any arm.

23. Context

    Where applicable we will report our review in concordance with the [Human Genome Epidemiology Network (HuGE Net) Handbook of Systematic Reviews](https://www.cdc.gov/genomics/hugenet/participate.htm). 

24. Main outcomes

    Ranked list of genes associated with BPD susceptibility, severity, and/or outcome.

26. Data extraction

    Selection of studies

    Article titles and abstracts obtained using the search strategy will be stored using reference management software (endNote X9, Clarivate Analytics, United States). Intital screening of titles wil be conducted by single authors against eligibility criteria, using the Screenatron tool (Systematic Review Accelerator, Bond University, Australia). Thereafter, screening of abstracts against eligibility criteria will be conducted by two authors independently. Inconsistencies wll be resolved in discussion with a thrd author. Full text articles will be retieved for studies matching the eligibility criteria. 

    Data extraction

    Data will be extracted by two independent reviewers using a pre-piloted proforma. 

    Gene lists will be extracted and their ranking preserved if possible. Ranking may be based on magnitude of effect or signficance. Where multiple measures are available we will preference magnitude of effect. Similarly, adjusted P values will be preferred over raw P values. If studies report multiple time points we will rank genes based on their minimum P value. We will exlcude genes for which the magnitude of effect or significance fall outside the authors threshold, or when this information is not available, for which P > 0.05, or z score < 1.96, or log fold change < 1.5. Gene, transcript, or protein identifiers will be mapped to its HUGO Gene Nomenclature Committee (HGNC) symbol. If one is not available we will use an equivalent Ensembl or Refseq symbol.

    In addition, we will extract information relating to study design, methodology, tissue/cell type, demographics, ARDS aetiology, risk factors, severity, and outcomes. 

27. Risk of bias (quality) assessment

    All genome-wide association studies will be assessed for risk of bias using domain-based evaluation as described in the [Q-Genie tool](https://doi.org/10.1186/s12863-015-0211-2). Studies will be classified as low, moderate, or high quality from their overall score. 

    For each gene ranked in the "top 50" by our meta-analysis we will rate the cumulatative evidence for genetic association using the [Venice interim guidelines](https://doi.org/10.1093/ije/dym159).


28. Strategy for data synthesis

    MAIC

    We will conduct a meta-analysis by information content (MAIC) of extracted gene lists. We have previously described our MAIC methodolgy in detail (<https://doi.org/10.1038/s41467-019-13965-x>, <https://doi.org/10.1038/s41598-020-79033-3>, <https://doi.org/10.1038/s41586-020-03065-y>). All components of our core algorithm can be found at <https://github.com/baillielab/maic>. 

    Functional enrichment analysis

    We will conduct gene set enrichment analysis based on rankings by MAIC score. We will use two methods: 1. FGSEA in R (<http://bioconductor.org/packages/release/bioc/html/fgsea.html>) using the full ranked list, and 2. over-enrichment analysis using Enrichr (<https://maayanlab.cloud/Enrichr/>) on the "top 100" genes. 

    We will control for false discovery using the Benjamini-Hochberg procedure (FDR < 0.05). 

29. Analysis of subgroups or subsets

    If sufficient data are available we will conduct sub-group analyses based on BPD aetiology e.g., PMA at birth, viral infection. 

    If sufficient data are available we will conduct sub-group analyses based on ancestory. 

30. Type and method of review

    Type of review

[comment]: # should be "yes" to Epidemiologic/Meta-analysis/Systematic review

    Health area of the review

[comment]: # should be "yes" to Genetics/Respiratory disorders

31. Language

    English

32. Country

    United Kingdom

33. Other registration details

    NA

34. Reference and/or URL for published protocol

    NA

35. Dissemination plans

    Do you intend to publish the review on completion?

    Yes

    Give brief details of plans for communicating review findings?

    We will make our manuscript available on medRxiv simultaneously with submission for publication.

    We will make the curated list of references and genes, as well as the results of our meta-analysis, available at <https://baillielab.net> following the publication of our manuscript.

36. Keywords

    Acute Respiratory Distress Syndrome; Genetics; Genomics; Susceptibility. 

37. Details of any existing review of the same topic by the same authors

    Our group has previously employed similar methodology to study the host genetics of [Influenza A](https://doi.org/10.1038/s41467-019-13965-x) and [SARS-CoV-2](https://doi.org/10.1038/s41598-020-79033-3) infection. 

38. Current review status

    Review ongoing.

39. Any additional information

    None.

40. Details of final report/publication(s) or preprints if available

    NA

[comment]: # overall, we should include enough detail to have a credible protocol but not too much that we tie our hands if we need to change something at a later date.



###
"Bronchopulmonary Dysplasia"[Mesh] AND ((gene*.[Title/Abstract]) OR (genome*[Title/Abstract]) OR (transcript*[Title/Abstract]) OR (protein*[Title/Abstract]) OR ("Susceptibility"[Title/Abstract]) OR (siRNA[All Fields]))
###



