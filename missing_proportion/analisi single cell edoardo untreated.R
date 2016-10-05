#test the distribution of failed genes
setwd("C:/Users/edoardo/Desktop/seattle experiments/single cells analysis/160720 single cell vincristine 10 nM")
dir()
df<-read.csv("1361334033 untreated.csv",header = T,skip = 11)
library(purrr)
#filter the df
df<-df[,c("ID","Name","rConc","Name.1","Value")]
#add the missing labeò
missing<-vector(mode = "numeric",length = nrow(df))
missing[df$Value!=999]<-1
df$missing<-missing
#separate the column if the ID
library(tidyr)
df2<-separate(data = df,col = ID,into = c("sample","assay"),sep = "-")
df2$sample<-as.character(df2$sample)

d<-split(x = df2,f = match(df2$sample,unique(df2$sample)))

#define the value of expression per gene. if one of the repicate is missing take the other 
#if bot are missing or present meke the mean

#make the mean of the value
d2<-map(d,function(p){
  dxx<-data.frame()
  genes<-as.character(unique(p$Name.1))
  for(i in genes){
    dx<-data.frame()
    dg<-p[p$Name.1==i,]
    if(sum(dg$missing)==0){
      dx<-dg[1,]
    }else if(sum(dg$missing)==1){
      dx<-dg[dg$missing==1,]
    }else if(sum(dg$missing)==2){
      avg<-mean(dg$Value)
      dx<-dg[1,]
      dx$Value<-avg
    }
    dxx<-rbind.data.frame(dxx,dx)
  }
  dxx
  })
#find the all human rna
nam<-map_chr(d2,function(x){
  as.character(unique(x$Name))
})
#identify the empty and the all human RNA

empty<-which(nam=="empty")
RNA<-which(nam=="all human rna")
rem<-c(empty,RNA)

d_cells<-d2[-rem]

d_control<-d2[rem]

#identify the percentage of missing gene per sample
per_cell_missing<-map_dbl(d_cells,function(x){
  1-mean(x$missing)
})
library(ggplot2)
gg<-data.frame("value"=per_cell_missing,"group"=rep(1,length(per_cell_missing)))
ggplot(gg,aes(y=value,x=as.factor(group)))+
  geom_boxplot()+
  geom_jitter(width = 0.2)

#make the list back to df
df_cells<-data.frame()
for(i in 1:length(d_cells)){
  df_cells<-rbind.data.frame(df_cells,d_cells[[i]])
}

df_control<-data.frame()
for(i in 1:length(d_control)){
  df_control<-rbind.data.frame(df_control,d_control[[i]])
}

#split the df per genes and identify the % working per genes

d_genes<-split(df_cells,f = match(df_cells$Name.1,unique(df_cells$Name.1)))
names(d_genes)<-unique(df_cells$Name.1)
data.frame(sort(map_dbl(d_genes,function(x){
  mean(x$missing)
})))

#check the same thing in the control df
d_genes_control<-split(df_control,f = match(df_control$Name.1,unique(df_control$Name.1)))
names(d_genes_control)<-unique(df_control$Name.1)
data.frame(sort(map_dbl(d_genes_control,function(x){
  mean(x$missing)
})))
