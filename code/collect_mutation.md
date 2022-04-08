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

requiredPackages2=c("org.Hs.eg.db")
for (k in requiredPackages2){
  if (!requireNamespace(k)) BiocManager::install(k,update=FALSE)
  library(k,character.only = TRUE)
  }
```

## 1.2.1 Loading data

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


## 1.2.2 Get unique gene frequncy

```r
# there's duplication in mutated gene for cell lines
cell_list_18Q2=data.frame(as.character(data.frame(table(mutation_18Q2[,16]))[,1]))
cell_list_22Q1=data.frame(as.character(data.frame(table(mutation_22Q1[,16]))[,1]))



#####################################
# custom function for lapply
# grepl not good
# get_unique_gene_for_each_cell =function(cell_list,mutation){
#   gene_list_18Q2=mutation[grepl(cell_list,mutation[,16],ignore.case = TRUE),c(1,2,16)]
#   gene_list_18Q2=dplyr::distinct(gene_list_18Q2, gene_list_18Q2[,1],.keep_all=TRUE)
#   return(gene_list_18Q2)
# }
# # use which instead
# get_unique_gene_for_each_cell =function(cell_list,mutation){
#   gene_list_single=mutation[which(mutation[,16]==cell_list),c(1,2,16)]
#   gene_list_single=dplyr::distinct(gene_list_18Q2, gene_list_18Q2[,1],.keep_all=TRUE)
#   return(gene_list_single)
# }
# 
# 
#   
# # lapply
# gene_list_18Q2=lapply(X=cell_list_18Q2[,1],FUN=get_unique_gene_for_each_cell2,mutation=mutation_18Q2)
# gene_list_22Q1=lapply(X=cell_list_22Q1[,1],FUN=get_unique_gene_for_each_cell2,mutation=mutation_22Q1)
# # convert list to a whole dataframe
# gene_list_18Q2=as.data.frame(do.call(rbind,test))
# gene_list_22Q1=as.data.frame(do.call(rbind,test2))
################################




#######################
# get difference one, and check, found grepl made mistakes
# dim1_test=as.data.frame(do.call(rbind,lapply(test,dim))[,1])
# dim1_test2=as.data.frame(do.call(rbind,lapply(test2,dim))[,1])
# difference_test12=dim1_test-dim1_test2
# 
# cell_list_diff_18Q2_name=as.data.frame(cell_list_18Q2[which(difference_test12!=0),])
# cell_list_diff_18Q2=which(difference_test12!=0)
######################





## if alternative
# 18Q2 first
gene_list_18Q2=list()
gene_list_18Q2_truegene=list()

for (i in 1:length(cell_list_18Q2[,1])){
  separate_cellline=mutation_18Q2[which(mutation_18Q2[,16]==cell_list_18Q2[,1][i]),c(1,2,16)]
  # separate_cellline=mutation_18Q2[grepl(cell_list_18Q2[,1][i],mutation_18Q2[,16],ignore.case = TRUE),c(1,2,16)]
  separate_cellline=dplyr::distinct(separate_cellline, separate_cellline[,1],.keep_all=TRUE)
  gene_list_18Q2[[i]]=separate_cellline
  separate_cellline_truegene=separate_cellline[separate_cellline[,2]!=0,]
  gene_list_18Q2_truegene[[i]]=separate_cellline_truegene
}



# View(gene_list_18Q2[[2]])
# View(gene_list_18Q2_truegene[[2]])

gene_list_18Q2_all=as.data.frame(do.call(rbind,gene_list_18Q2))

gene_list_18Q2_truegene_all=as.data.frame(do.call(rbind,gene_list_18Q2_truegene))

freq_18Q2_unique=as.data.frame(table(gene_list_18Q2_all[,1]))
freq_18Q2_unique_truegene=as.data.frame(table(gene_list_18Q2_truegene_all[,1]))


## 22Q1 next
gene_list_22Q1=list()
gene_list_22Q1_truegene=list()
for (i in 1:length(cell_list_22Q1[,1])){
  separate_cellline=mutation_22Q1[which(mutation_22Q1[,16]==cell_list_22Q1[,1][i]),c(1,2,16)]
  # separate_cellline=mutation_22Q1[grepl(cell_list_22Q1[,1][i],mutation_22Q1[,16],ignore.case = TRUE),c(1,2,16)]
  separate_cellline=dplyr::distinct(separate_cellline, separate_cellline[,1],.keep_all=TRUE)
  gene_list_22Q1[[i]]=separate_cellline
  separate_cellline_truegene=separate_cellline[separate_cellline[,2]!=0,]
  gene_list_22Q1_truegene[[i]]=separate_cellline_truegene
}

# View(gene_list_22Q1[[2]])
# View(gene_list_22Q1_truegene[[2]])

gene_list_22Q1_all=as.data.frame(do.call(rbind,gene_list_22Q1))

gene_list_22Q1_truegene_all=as.data.frame(do.call(rbind,gene_list_22Q1_truegene))

freq_22Q1_unique=as.data.frame(table(gene_list_22Q1_all[,1]))
freq_22Q1_unique_truegene=as.data.frame(table(gene_list_22Q1_truegene_all[,1]))


# convert gene symbol to Entrez ID
hs <- org.Hs.eg.db

freq_22Q1_fakegene=freq_22Q1[!freq_22Q1[,1] %in% freq_22Q1_unique_truegene[,1],]
```

```
## Error in eval(expr, envir, enclos): object 'freq_22Q1' not found
```

```r
my.symbols <- freq_22Q1_fakegene[,1]
```

```
## Error in eval(expr, envir, enclos): object 'freq_22Q1_fakegene' not found
```

```r
freq_22Q1_fakegene_noentrezid=select(hs, 
       keys = my.symbols,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")
```

```
## Error in .testForValidKeys(x, keys, keytype, fks): object 'my.symbols' not found
```

```r
write.table(file=here::here("data","22Q1_fakegene_noentrezid.csv"),freq_22Q1_fakegene_noentrezid[,1] ,col.names = F,row.names = F,sep =",",quote=F)
```

```
## Error in is.data.frame(x): object 'freq_22Q1_fakegene_noentrezid' not found
```

```r
write.table(file=here::here("data","22Q1_gene_list_all.txt"),freq_22Q1[,1] ,col.names = F,row.names = F,sep ="\t",quote=F)
```

```
## Error in is.data.frame(x): object 'freq_22Q1' not found
```

## 1.3 extract top 15%

```r
# extract top 15%

freq_18Q2=as.data.frame(table(mutation_18Q2[,1]))
freq_18Q2[,1]=as.character(freq_18Q2[,1])

freq_18Q2=dplyr::arrange(freq_18Q2,desc(freq_18Q2[,2]))
freq_18Q2_15P <- freq_18Q2[1:(dim(freq_18Q2)[1]*0.15),]
freq_18Q2_15P[grepl("CACNA1B|BRAF",freq_18Q2_15P[,1],ignore.case = TRUE),]
```

```
##        Var1 Freq
## 138 CACNA1B  357
## 282    BRAF  273
```

```r
freq_22Q1=as.data.frame(table(mutation_22Q1[,1]))
freq_22Q1[,1]=as.character(freq_22Q1[,1])

freq_22Q1=dplyr::arrange(freq_22Q1,desc(freq_22Q1[,2]))
freq_22Q1_15P <- freq_22Q1[1:(dim(freq_22Q1)[1]*0.15),]
freq_22Q1_15P[grepl("CACNA1B|BRAF",freq_22Q1_15P[,1],ignore.case = TRUE),]
```

```
##        Var1 Freq
## 224 CACNA1B  307
## 337    BRAF  264
```

```r
# somehow, the frequncy of these 2 genes are decreased
```

## 1.4 check common genes to determine threshold(top 15% are not really top 15%)

```r
genename_DC=read.csv2(here::here("data","gene2ind.txt"),header=FALSE,sep="\t")
genename_DC=as.data.frame(genename_DC[,2])
genename_DC_sorted=dplyr::arrange(genename_DC,genename_DC)

freq_18Q2=as.data.frame(freq_18Q2)
matched_18Q2=freq_18Q2[freq_18Q2[,1] %in% genename_DC_sorted[,1],]
matched_18Q2_unique_truegene=freq_18Q2_unique_truegene[freq_18Q2_unique_truegene[,1] %in% genename_DC_sorted[,1],]


freq_22Q1=as.data.frame(freq_22Q1)
matched_22Q1=freq_22Q1[freq_22Q1[,1] %in% genename_DC_sorted[,1],]
matched_22Q1_unique_truegene=freq_22Q1_unique_truegene[freq_22Q1_unique_truegene[,1] %in% genename_DC_sorted[,1],]

# only not matched is GSTT1, because of in both release of CCLE mutation, the entrez ID is missing.
notmatched_18Q2_unique_truegene=genename_DC_sorted[!genename_DC_sorted[,1] %in% matched_18Q2_unique_truegene[,1],]
notmatched_22Q1_unique_truegene=genename_DC_sorted[!genename_DC_sorted[,1] %in% matched_22Q1_unique_truegene[,1],]
```



## 1.5 check paper mutation gene with cell name, CORL51_LUNG_DC,ACH-001047


```r
cellmutation_name_DC=cbind(cellname_DC[,2],cellmutation_DC)
cellmutation_number_DC=data.frame(cellname_DC[,2],as.numeric(rowSums(cellmutation_name_DC[,2:3009])))

cellmutation_number_DC=dplyr::arrange(cellmutation_number_DC,cellmutation_number_DC[,2])



# CORL51_LUNG_DC,ACH-001047
DC_CORL51_LUNG=cellmutation_name_DC[which(cellmutation_name_DC[,1]=="CORL51_LUNG"),]

DC_CORL51_LUNG_genelist=data.frame(genename_DC[which(DC_CORL51_LUNG[,2:3009]==1),])

R18Q2_CORL51_LUNG_genelist_unique_notunique=mutation_18Q2[grepl("CORL51_LUNG",mutation_18Q2[,16],ignore.case = TRUE),]

R18Q2_CORL51_LUNG_genelist_unique=unique(data.frame(mutation_18Q2[grepl("CORL51_LUNG",mutation_18Q2[,16],ignore.case = TRUE),][,1]))
rownames(R18Q2_CORL51_LUNG_genelist_unique)=c()

comparison=compare(R18Q2_CORL51_LUNG_genelist_unique[,1],DC_CORL51_LUNG_genelist[,1],allowAll=TRUE)


setdiff(R18Q2_CORL51_LUNG_genelist_unique[,1],DC_CORL51_LUNG_genelist[,1])
```

```
##  [1] "RBM15"  "ADCK3"  "TRIM67" "FMN2"   "TTN"    "SPEG"   "GATA2"  "TRIM36"
##  [9] "BAG6"   "PKHD1"  "GLI3"   "TRRAP"  "LETM2"  "COPS5"  "PSKH2"  "CSMD3" 
## [17] "EIF3H"  "EXT1"   "ERCC6"  "MPPED2" "NIN"    "KTN1"   "CSPG4"  "TEX14" 
## [25] "LAMA1"  "ALPK2"  "NRK"
```

```r
all(R18Q2_CORL51_LUNG_genelist_unique[,1] %in% DC_CORL51_LUNG_genelist[,1])
```

```
## [1] FALSE
```

```r
all(DC_CORL51_LUNG_genelist[,1] %in% R18Q2_CORL51_LUNG_genelist_unique[,1])
```

```
## [1] TRUE
```



