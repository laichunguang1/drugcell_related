---
title: "Collecting mutation gene list"
author: "Chunguang Lai"
date: "2022-02-18"
output:
  html_document:
    toc: true
    toc_float: true
    keep_md: true
editor_options: 
  chunk_output_type: console
---



# 1.Preparation 
## 1.1 Enviroment setup

```r
requiredPackages = c('here','compare',"dplyr")
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}

requiredPackages2=c()
for (k in requiredPackages2){
  if (!requireNamespace(k)) BiocManager::install(k)
  library(k,character.only = TRUE)
  }
```

## 1.2 Loading data

```r
mutation_18Q2=read.csv2(here::here("data","cell_mutation","CCLE_DepMap_18Q2_maf_20180502.txt"),header=TRUE,sep="\t")

mutation_22Q1=read.csv2(here::here("data","cell_mutation","CCLE_mutations_22Q1.csv"),header=TRUE,sep=",")

cellmutation_DC=read.csv2(here::here("data","cell2mutation.txt"),header=FALSE,sep=",")
cellname_DC=read.csv2(here::here("data","cell2ind.txt"),header=FALSE,sep="\t")

# total gene numbers and total cell line number respectively

# 18704 genes, 19284 Hugo_Symbol, 1549 cell lines
dim(table(mutation_18Q2[,2]))
```

```
## [1] 18704
```

```r
dim(table(mutation_18Q2[,16]))
```

```
## [1] 1549
```

```r
# 18784 genes, 19537 Hugo_Symbol, 1759 cell lines
dim(table(mutation_22Q1[,2]))
```

```
## [1] 18784
```

```r
dim(table(mutation_22Q1[,16]))
```

```
## [1] 1759
```


## 1.2 Extract top 15%

```r
table_18Q2=as.data.frame(table(mutation_18Q2[,1]))
table_18Q2[,1]=as.character(table_18Q2[,1])

table_18Q2=dplyr::arrange(table_18Q2,desc(table_18Q2[,2]))
table_18Q2_15P <- table_18Q2[1:(dim(table_18Q2)[1]*0.15),]
table_18Q2_15P[grepl("CACNA1B|BRAF",table_18Q2_15P[,1],ignore.case = TRUE),]
```

```
##        Var1 Freq
## 138 CACNA1B  357
## 282    BRAF  273
```

```r
table_22Q1=as.data.frame(table(mutation_22Q1[,1]))
table_22Q1[,1]=as.character(table_22Q1[,1])

table_22Q1=dplyr::arrange(table_22Q1,desc(table_22Q1[,2]))
table_22Q1_15P <- table_22Q1[1:(dim(table_22Q1)[1]*0.15),]
table_22Q1_15P[grepl("CACNA1B|BRAF",table_22Q1_15P[,1],ignore.case = TRUE),]
```

```
##        Var1 Freq
## 224 CACNA1B  307
## 337    BRAF  264
```

```r
# somehow, the frequncy of these 2 genes are decreased
```

## 1.3 check common genes to determine threshold(top 15% are not really top 15%)

```r
genename_DC=read.csv2(here::here("data","gene2ind.txt"),header=FALSE,sep="\t")
genename_DC=as.data.frame(genename_DC[,2])
genename_DC_sorted=dplyr::arrange(genename_DC,genename_DC)

table_18Q2_genename=as.data.frame(table_18Q2)


# matched1=table_18Q2_genename[grepl(paste(genename_DC_sorted[1:2008,1],collapse = "|"),table_18Q2_genename[,1],fixed=TRUE),]

matched1=table_18Q2_genename[table_18Q2_genename[,1] %in% genename_DC_sorted[,1],]


table_22Q1_genename=as.data.frame(table_22Q1)

matched2=table_22Q1_genename[table_22Q1_genename[,1] %in% genename_DC_sorted[,1],]
```



## 1.4 check paper mutation gene with cell name


```r
cellmutation_name_DC=cbind(cellname_DC[,2],cellmutation_DC)
cellmutation_number_DC=data.frame(cellname_DC[,2],as.numeric(rowSums(cellmutation_name_DC[,2:3009])))

cellmutation_number_DC=dplyr::arrange(cellmutation_number_DC,cellmutation_number_DC[,2])



# CORL51_LUNG_DC,ACH-001047
DC_CORL51_LUNG=cellmutation_name_DC[which(cellmutation_name_DC[,1]=="CORL51_LUNG"),]

DC_CORL51_LUNG_genelist=data.frame(genename_DC[which(DC_CORL51_LUNG[,2:3009]==1),])

R18Q2_CORL51_LUNG_genelist=unique(data.frame(mutation_18Q2[grepl("CORL51_LUNG",mutation_18Q2[,16],ignore.case = TRUE),][,1]))
rownames(R18Q2_CORL51_LUNG_genelist)=c()

comparison=compare(R18Q2_CORL51_LUNG_genelist[,1],DC_CORL51_LUNG_genelist[,1],allowAll=TRUE)
comparison
```

```
## FALSE
##   shortened model
##   sorted
##   ignored case
```

```r
setdiff(R18Q2_CORL51_LUNG_genelist[,1],DC_CORL51_LUNG_genelist[,1])
```

```
##  [1] "RBM15"  "ADCK3"  "TRIM67" "FMN2"   "TTN"    "SPEG"   "GATA2"  "TRIM36"
##  [9] "BAG6"   "PKHD1"  "GLI3"   "TRRAP"  "LETM2"  "COPS5"  "PSKH2"  "CSMD3" 
## [17] "EIF3H"  "EXT1"   "ERCC6"  "MPPED2" "NIN"    "KTN1"   "CSPG4"  "TEX14" 
## [25] "LAMA1"  "ALPK2"  "NRK"
```

```r
all(R18Q2_CORL51_LUNG_genelist[,1] %in% DC_CORL51_LUNG_genelist[,1])
```

```
## [1] FALSE
```

```r
all(DC_CORL51_LUNG_genelist[,1] %in% R18Q2_CORL51_LUNG_genelist[,1])
```

```
## [1] TRUE
```

```r
library(sets)
```

```
## 
## Attaching package: 'sets'
```

```
## The following object is masked from 'package:dplyr':
## 
##     %>%
```



