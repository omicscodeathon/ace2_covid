# Expression Level Analysis of ACE2 Receptor to Determine the Effect of Malaria Endemicity on CoVID-19 Outcomes

## Global Malaria belt

![image](https://user-images.githubusercontent.com/45264074/161214522-0ad8b18c-7530-485d-b877-875bac477a2a.png)


source: http://www.malariacampaign.gov.lk/en/travelers-guide

## Background

The global COVID-19 pandemic caused by SARS-CoV-2 has spread rapidly across the continents. While the incidence and mortality rate of COVID-19 in Africa has been reported to be lower than other continents, the Malaria rate in Africa is higher than in the non-African population. ACE2 plays a role in both Malaria and COVID-19 infections. SARS-CoV-2 uses the ACE2 enzyme to enter host cells, while downregulation of ACE2 leads to accumulation of angiotensin II substrate, which impairs Plasmodium development and thus has a protective effect in Malaria. Although the difference in COVID-19 incidence can be explained by many factors such as low testing capacity in Africa, the variable distribution of the ACE2 gene has been linked to this observed phenomenon. Little is known about ACE2 expression in African COVID-19 patients compared to non-African COVID-19 patients.

## General Objective
To analyze expression level of ACE2 gene in African COVID-19 patients compared to non- African COVID-19 patients.

## Specific objectives
1. To determine the expression level of ACE2 gene in African and non-African COVID-19 patients
2. To determine how the expression level of ACE2 gene in African COVID-19 patients compares to that of non- African COVID-19 patients

## Methods
In this study, transcriptomes from African and non-African COVID-19 patients were downloaded from the NCBI sequence read archive (SRA) and analyzed for ACE2 gene expression. The reads were aligned to the human reference genome using HISAT2 and the raw counts generated using HTseq-count. EdgeR was used to conduct differential gene expression analysis with statistical analysis being done in R.

## Analysis Workflow

![image](https://user-images.githubusercontent.com/45264074/160587704-756a14cf-982f-43ab-9508-646d8b3e8f50.png)

## Prerequisites for Analysis

The data used for the analysis was retrieved from NCBI. The accession numbers for the data used can be found under [accession](accessions/accessions.txt).

#### Tools and packages.
Prior to installation of the analysis packages you need to have a unix environment i.e linux or macOS. Moreover you need to install [R](https://cran.r-project.org/), [RStudio IDE](https://www.rstudio.com/products/rstudio/download/) and [miniconda](https://docs.conda.io/en/latest/miniconda.html).

###### Installation of unix packages.
 ```bash
 conda install -c bioconda fastqc
```

``` bash
conda install -c bioconda cutadapt
```

```bash
conda install -c bioconda hisat2
```

```bash
conda install -c bioconda htseq
```

```bash
conda install -c bioconda samtools
```


###### Installation of packages in RStudio

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install("edgeR")
```

## Team members.

1. [Marion Nyaboke](https://github.com/marionnyaboke) Student, Msc. Bioinformatics, Pwani University, Kenya | Project/Tech Lead.
2. [Kauthar M. Omar](https://github.com/Kauthar-Omar) Msc. Bioinformatics student, Pwani University, Kenya | Writers Lead.
3. [Ayorinde F. Fayehun](https://github.com/Ayor1) Biotechnology unit, Institute of Child Health, University of Ibadan
4. [Oumaima Dachi]() National School of Arts and Crafts, University of Hassan II Casablanca, Morocco.
5. [Billiah Bwana]() University of Embu
