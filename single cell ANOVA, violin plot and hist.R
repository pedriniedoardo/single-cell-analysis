setwd("C:/Users/edoardo/Desktop/seattle experiments/single cells analysis/160420 single cell HL60 experiment/single cell analysis vincristine")
df<-read.csv("data mock vs 48h.csv",skip = 11)
library(ggplot2)
library(reshape)
#is better to be sure to remove the factor and transform the variable in character
df$ID<-as.character(df$ID)
df$Name<-as.character(df$Name)
df$Name.1<-as.character(df$Name.1)
##################################
##################################
##################################
a<-{}
b<-{}
c<-0
v<-{}
#how is organized the basic df?
for(i in 1:96){
  print(i)
  a<-i+c
  b<-96+c
  v<-rbind(v,unique(t(as.data.frame(strsplit(df[a:b,1],split = "-")))[,1]))
  c<-96*i
}

a<-{}
b<-{}
c<-0
x<-{}
#how is organized the basic df?
for(i in 1:96){
  print(i)
  a<-i+c
  b<-96+c
  x<-rbind(x,unique(t(as.data.frame(strsplit(df[a:b,2],split = "-")))[,1]))
  c<-96*i
}

#the data frame is ordered for cells
#prepare a fake label on the bases of the 
l<-rep(101:196,each=96)
l2<-rep(1:96,each=96)
#part from the df
t(as.data.frame(strsplit(df$Name,split = "_")))
first<-t(as.data.frame(strsplit(df$Name,split = "_")))[,1]
second<-t(as.data.frame(strsplit(df$Name,split = "_")))[,2]

label<-paste0(first,second,"_",l,"-",l2)

df_2<-df
df_2$treatment<-df_2$Name
df_2$Name<-as.character(label)


#sort out all the 100 cell and the all human rna control
df_3<-df_2[df_2$rConc==1,]
#make all the miss and double peak as fail
#table(df_3$Comments)
#df_3$Call[df_3$Comments=="double" | df_3$Comments=="miss"]<-"Fail"

#collapse the tecnical replicate
px<-{}
#ciclo1
for(i in unique(df_3$Name)){
  print(i)
  p<-df_3[df_3$Name==i,]
  #ciclo2
  for(l in unique(df_3$Name.1)){
    print(l)
    p2<-p[p$Name.1==l,]
    if(p2$Call[1]=="Fail" & p2$Call[2]=="Fail"){
      px<-rbind.data.frame(px,p2[1,])
      }
    else if(p2$Call[1]=="Fail" & p2$Call[2]!="Fail"){
      px<-rbind.data.frame(px,p2[2,])
      }
    else if(p2$Call[1]!="Fail" & p2$Call[2]=="Fail"){
      px<-rbind.data.frame(px,p2[1,])
      }
    else if(p2$Call[1]!="Fail" & p2$Call[2]!="Fail"){
      p3<-p2
      p3[1,7]<-mean(p2$Value)
      px<-rbind.data.frame(px,p3[1,])
    }
  }
}
write.csv(px,"C:/Users/edoardo/Desktop/seattle experiments/single cells analysis/160420 single cell HL60 experiment/single cell analysis vincristine/dataframe fluidigm.csv")

#I cannot make it work I downloaded the manual to procede with the regular analysis
#I downloaded the guidelines for the quantification
#definisci il LOD
LOD<-24
dat<-px
#how many value are over the lod
dat$Value[dat$Value>LOD & dat$Comments=="ok"]
#tarnsform the ok data above the LOD in the LOD
dat$Value[dat$Value>LOD & dat$Comments=="ok"]<-LOD

#what is the shape of the 999
table(dat[dat$Value==999,7],dat[dat$Value==999,15])
#what is the shape of comments
table(dat$Comments)

#new variable
dat$label[dat$Comments=="ok"]<-"ok"
dat$label[dat$Comments=="no peak"]<-"no peak"
dat$label[(dat$Comments=="miss" | dat$Comments=="double") & dat$Value>LOD]<-"almost no peak"
dat$label[dat$Comments=="double" & dat$Value<LOD]<-"double"
dat$label[dat$Comments=="miss" & dat$Value<LOD]<-"miss"

#check NAs
sum(is.na(dat$label))

#is better to wait to transform the miss 999 in LOD
#there are miss non 999
#transform the missing whit 999 not miss in LOD 
#dat$Value[dat$Value==999]<-LOD
#deal with the double and the missing now different from LOD
t<-table(dat[,16],dat[,5])
#there genes with a lot of double peak is it safe to merge the value?

#remove the genes where there is not even an ok
df_f<-{}
for(i in unique(dat$Name.1)){
  print(i)
  d<-dat[dat$Name.1==i,]
  if(sum(d$label=="ok")!=0){
    df_f<-rbind.data.frame(df_f,d)
  }
}
#control they have been removed
table(df_f$label,df_f$Name.1)
#convert all the no peak and almost no peak in LOD
df_f$Value[df_f$label=="no peak" | df_f$label=="almost no peak"]<-24
head(table(df_f$label,df_f$Value))
#convert the double and the miss in the mean of the ok
#note that I am transforming in the mean Ct of the ok only whithin the same sample group
for(g in unique(df_f$treatment)){
  print(g)
for(i in unique(df_f$Name.1)){
  print(i)
  df_f$Value[df_f$Name.1==i &
               (df_f$label=="double" | df_f$label=="miss") &
               df_f$treatment==g]<-
    mean(df_f$Value[df_f$Name.1==i &
                      df_f$label=="ok"&
                      df_f$treatment==g])
}
}
##
#now i am ready to transform all the value in log2 of expression
dfGeneExp<-df_f
dfGeneExp$Value<-24-df_f$Value
#verify there is no negative number in the df
sum(dfGeneExp$Value<0)
#######################
#######################
#look at the distribution of a random genes
s<-sample(dfGeneExp$Name.1,1)
boxplot(dfGeneExp$Value[dfGeneExp$Name.1==s]~dfGeneExp$treatment[dfGeneExp$Name.1==s])
#for the available genes as expected there is the same amount of sample
table(dfGeneExp$Name.1)

###########################
#remove the sample with an overall low expression
###########################
#look at the distribution per sample
#reorder the data in decreasing median
dfGeneExp$Name<-as.character(dfGeneExp$Name)
med<-{}
for(i in 1:length(unique(dfGeneExp$Name))){
  print(i)
  l<-as.character(unique(dfGeneExp$Name)[i])
  med[i]<-median(dfGeneExp$Value[dfGeneExp$Name==l])
}
#define the indexing for decreasing median
index<-order(med,decreasing = T)
#redifine the list of name on the base of the index
nameIndex<-unique(dfGeneExp$Name)[index]
#set up a df to produce the boxplot
df_box<-{}
for(i in nameIndex){
  print(i)
  df_box<-rbind.data.frame(df_box,dfGeneExp[dfGeneExp$Name==i,])
}
df_box$Name<-factor(df_box$Name,levels = unique(df_box$Name))
boxplot(df_box$Value~df_box$Name)
#######################
#######################
#define the threshold
#######################
#######################
#make a chart of the overall expression of the genes
e<-{}
for(i in 1:length(unique(dfGeneExp$Name.1))){
  l<-unique(dfGeneExp$Name.1)[i]
  print(l)
  e[i]<-sum(dfGeneExp$Value[dfGeneExp$Name.1==l])
}
df_e<-cbind.data.frame("gene"=unique(dfGeneExp$Name.1),"exp"=e)
df_e<-df_e[order(df_e$exp,decreasing = F),]
chart<-as.character(df_e$gene)
#set-up a df with the expression defined by the order for every cell
df_out<-{}
for(i in unique(dfGeneExp$Name)){
  print(i)
  z<-dfGeneExp$Value[dfGeneExp$Name==i]
  names(z)<-dfGeneExp$Name.1[dfGeneExp$Name==i]
  df_out<-rbind.data.frame(df_out,z[chart])
}
colnames(df_out)<-chart
rownames(df_out)<-unique(dfGeneExp$Name)

abline(h=5,col="red")
#############################
#start the trimming
for(i in 1:ncol(df_out)){
  print(i)
  df_out2<-df_out[,i:ncol(df_out)]
  #the 95% of gene are raunded to 46 gene different from 0
  su<-rowSums(df_out2==0)
  if (sum(su>=2)<length(su)/2){
    break
  }
}
#identify the sample used as control
nameCtrlSamples<-names(su[su==0])
##############################
#the threshold is defined as the 15 percentile of the all the data of the control samples combined
out<-{}
for(i in nameCtrlSamples){
  print(i)
  out<-c(out,dfGeneExp$Value[dfGeneExp$Name==i])
}
hist(out)
#finally define the quantile
q<-quantile(out,0.15)
#######################
#it's definitely better to remove some cell that have clear problem
#so remove from the df the sample where the median expression is below the treshold
#
df_filter<-cbind.data.frame("median"=med,"sample"=unique(dfGeneExp$Name))
SampleKeep<-as.character(df_filter$sample[df_filter$median>q])
boxplot(df_box$Value~df_box$Name)
abline(h=q,col="red")
#######################
#######################
#remove the sample that are identified as outlier
dfGeneExp2<-{}
for(i in SampleKeep){
  print(i)
  dfGeneExp2<-rbind.data.frame(dfGeneExp2,dfGeneExp[dfGeneExp$Name==i,])
}
table(dfGeneExp2$treatment)
#perform a one way anova for every gene
intGene<-{}
nonGene<-{}
noExp<-{}
tukey<-{}
val<-{}
for(i in unique(dfGeneExp2$Name.1)){
  print(i)
  d<-dfGeneExp2[dfGeneExp2$Name.1==i,c("Value","treatment")]
  if(sum(d$Value)==0){
    noExp<-c(noExp,i)
  }else{
    a<-aov(d$Value~d$treatment)
    s<-summary(a)
    p_val<-s[[1]][[5]][1]
    if(p_val<0.05){
      #save the name of the gene of interest
      intGene<-c(intGene,i)
      val<-c(val,p_val)
      t<-TukeyHSD(a)
      tukey[[i]]<-t$`d$treatment`[t$`d$treatment`[,4]<0.05,]
    }else{
      nonGene<-c(nonGene,i)
    } 
  }
}
#create a df with ordered gene name and p value
intGene<-cbind.data.frame(intGene,val)
#order the df in term of increasing p_val
intGene<-intGene[order(intGene$val),]

#build a violin plot only for the intGene
his<-{}
viol<-{}
for(i in 1:length(intGene$intGene)){
  print(i)
  l<-intGene$intGene[i]
  d<-dfGeneExp2[dfGeneExp2$Name.1==l,c("Value","treatment")]
  his[[i]]<-ggplot(d,aes(x=Value,fill=as.factor(treatment)))+
    geom_density(alpha=0.5)+
    labs(x="log2(Exp)",title=l)+
    scale_fill_discrete(name="treatment")
  viol[[i]]<-ggplot(d,aes(factor(treatment),Value))+
    geom_violin(aes(fill=treatment))+
    labs(y="log2(Exp)",x="",title=l)+
    geom_jitter(shape=16, position=position_jitter(0.2))
}
library(gridExtra)
do.call("grid.arrange",c(viol,ncol=3))
do.call("grid.arrange",c(his,ncol=3))
#use the value of expression of the 100 cell to understand if the proportion of missing is resonable


