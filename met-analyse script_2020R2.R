rm(list=ls())

#Packages
  library(car)
  library(caper)
  library(cluster)
  library(dplyr)
  library(ecodist)
  library(effects)
  library(ggplot2)
  library(glmulti)
  library(igraph)
  library(leaps)
  library(multcomp)
  library(MuMIn)
  library(relaimpo)
  library(scales)
  library(vegan)
  library(ape)
  library(nlme)
  library(car)
  library(rcompanion)
  library(multcomp)
  library(Matrix)
  library(MuMIn)
  library(extrafont)
  library(quantreg)
  library(FSA)
  library(ecodist)


#Reading in dataset
datasetmeta<-read.table("D:/OneDrive/OneDrive - UGent/Thomasdoc/Projects/Ongoing/Meta-analyse/2020/metaanalysis2020R2.txt",header=TRUE)  
rownames(datasetmeta)<-datasetmeta[,1]
edit(datasetmeta)
tdat<-t(datasetmeta[,10:219]) #host-symbiont matrix without variables
table(datasetmeta$FUNCTIONAL_GROUP2)
edit(tdat)
edit(datasetmeta)
#######################################################################
#######################################################################
########## 1.  CIRCULAR PLOT PER GENUS ################################
#######################################################################
#######################################################################

  
#group species per genus
#########################


Aph<-tdat[grep("Aphaenogaster", rownames(tdat)), ]
Aphaenogaster<-apply(Aph, 2, max)

Bothriomyrmex<-tdat[grep("Bothr", rownames(tdat)), ]

Cam<-tdat[grep("Camponotus", rownames(tdat)), ]
Camponotus<-apply(Cam, 2, max)

Cat<-tdat[grep("Catag", rownames(tdat)), ]
Cataglyphis<-apply(Cat, 2, max)

Car<-tdat[grep("Cardioc", rownames(tdat)), ]
Cardiocondyla<-apply(Car, 2, max)

Cha<-tdat[grep("Chale", rownames(tdat)), ]
Chalepoxenus<-apply(Cha, 2, max)

Cre<-tdat[grep("Crematogaster", rownames(tdat)), ]
Crematogaster<-apply(Cre, 2, max)

For<-tdat[grep("ormica", rownames(tdat)), ]
Formica<-apply(For, 2, max)

Gon<-tdat[grep("Gonio", rownames(tdat)), ]
Goniomma<-apply(Gon, 2, max)

Har<-tdat[grep("Harpag", rownames(tdat)), ]
Harpagoxenus<-apply(Har, 2, max)

Iberoformica<-tdat[grep("Ibero", rownames(tdat)), ]

Las<-tdat[grep("Lasius", rownames(tdat)), ]
Lasius<-apply(Las, 2, max)

Lin<-tdat[grep("Linepithema", rownames(tdat)), ]
Linepithema<-apply(Lin, 2, max)

Lep<-tdat[grep("Lepto", rownames(tdat)), ]
Leptothorax<-apply(Lep, 2, max)

Lio<-tdat[grep("Liometopum", rownames(tdat)), ]
Liometopum<-apply(Lio, 2, max)

Man<-tdat[grep("Manica", rownames(tdat)), ]
Manica<-apply(Man, 2, max)

Mes<-tdat[grep("Messor", rownames(tdat)), ]
Messor<-apply(Mes, 2, max)

Mon<-tdat[grep("Monomorium", rownames(tdat)), ]
Monomorium<-apply(Mon, 2, max)

Myrmecina<-tdat[grep("Myrmecina", rownames(tdat)), ]

Myr<-tdat[grep("Myrmica", rownames(tdat)), ]
Myrmica<-apply(Myr, 2, max)

Phe<-tdat[grep("Phei", rownames(tdat)), ]
Pheidole<-apply(Phe, 2, max)

Pla<-tdat[grep("Plagio", rownames(tdat)), ]
Plagiolepis<-apply(Pla, 2, max)

Pon<-tdat[grep("Pone", rownames(tdat)), ]
Ponera<-apply(Pon, 2, max)

Prof<-tdat[grep("Prof", rownames(tdat)), ]
Proformica<-apply(Prof, 2, max)

Prenolepis<-tdat[grep("Prenolepis", rownames(tdat)), ]
Prenolepis

Strumigen<-tdat[grep("Strumigen", rownames(tdat)), ]
Strumigenys<-apply(Strumigen, 2, max)

Sol<-tdat[grep("Solenopsis", rownames(tdat)), ]
Solenopsis<-apply(Sol, 2, max)

Stenamma<-tdat[grep("Stenamma", rownames(tdat)), ]

Strongylognathus<-tdat[grep("Strongyl", rownames(tdat)), ]

Tap<-tdat[grep("Tapinoma", rownames(tdat)), ]
Tapinoma<-apply(Tap, 2, max)

Tem<-tdat[grep("Temnothorax", rownames(tdat)), ]
Temnothorax<-apply(Tem, 2, max)

Tet<-tdat[grep("Tetra", rownames(tdat)), ]
Tetramorium<-apply(Tet, 2, max)

#summarize totals per genus
############################
compile<-t(cbind(Camponotus,Plagiolepis,Cataglyphis,Lasius,Prenolepis,Iberoformica,Proformica,Formica,Liometopum,Linepithema, Bothriomyrmex,Tapinoma,Leptothorax,Harpagoxenus,Monomorium,Goniomma, Aphaenogaster,Cardiocondyla,Chalepoxenus,Myrmecina,Crematogaster,Messor,Tetramorium,Strongylognathus,Manica,Myrmica,Pheidole,Strumigenys,Solenopsis,Temnothorax,Stenamma,Ponera))
#totalen per functional groep
edit(datasetmeta)
compilefun<-compile[,709:722]
compilemyrm<-compile[,c(1:535)]
compilehelm<-compile[,687:708]
compilesocial<-compile[,536:606]
compiletroph<-compile[,c(607:686)]

speciestot<-c(rowSums(compile))
speciesfun<-c(rowSums(compilefun))
speciesmyrm<-c(rowSums(compilemyrm))
speciessocial<-c(rowSums(compilesocial))
speciestroph<-c(rowSums(compiletroph))
specieshelm<-c(rowSums(compilehelm))

tablespecies<-t(rbind(speciesmyrm,speciestroph,speciessocial,specieshelm,speciesfun))
# edit(tablespecies)
colSums(tablespecies)
colMeans(tablespecies/rowSums(tablespecies))
q<-(tablespecies/rowSums(tablespecies))
colMeans(q)
pie.values <- split(tablespecies, seq(nrow(tablespecies)))
pie.values <- setNames(split(tablespecies, seq(nrow(tablespecies))), rownames(tablespecies))
#edit(pie.values)

#pivot table unique species
tcompile<-t(compile)
dftcompile<-data.frame(tcompile)
dftcompile$FUNCTIONAL_GROUP<-datasetmeta$FUNCTIONAL_GROUP
dftcompilemin = subset(dftcompile, select = -c(dftcompile$FUNCTIONAL_GROUP) )

dftcompilemin<-data.frame(tcompile)
# edit(dftcompilemin)

uniqueall<-dftcompile[rowSums(dftcompilemin) == 1,]
# edit(uniqueall)

pivottabledf<-uniqueall %>% 
  group_by(FUNCTIONAL_GROUP) %>% 
  summarise_all(funs(sum))
# edit(pivottabledf)

minpivottabledf<-pivottabledf[,-1]

speciesfun2<-minpivottabledf[1,]
specieshelm2<-minpivottabledf[2,]
speciesmyrm2<-minpivottabledf[3,]
speciessocial2<-minpivottabledf[4,]
speciestroph2<-minpivottabledf[5,]

tablespecies2<-t(rbind(speciesmyrm2,speciestroph2,speciessocial2,specieshelm2,speciesfun2))
tablespecies2[tablespecies2 == 0] <- 0.001
pie.values2 <- split(tablespecies2, seq(nrow(tablespecies2)))
uniquetot<-c(rowSums(tablespecies2))
1-mean(round(uniquetot)/rowSums(tablespecies))


#graph
  #all
co_occurrence <- compile %*% t(compile)
# edit(co_occurrence)
graphall <- graph.adjacency(co_occurrence,
                         weighted=TRUE,
                         mode="undirected",
                         diag=FALSE)
ball<-layout_in_circle(graphall)

#egdge shared myrmecophiles
co_occurrencemyrm <- compilemyrm %*% t(compilemyrm)
# edit(co_occurrencemyrm)
graphmyrm<- graph.adjacency(co_occurrencemyrm,
                            weighted=TRUE,
                            mode="undirected",
                            diag=FALSE)


  #egdge shared trophobionts
co_occurrencetroph <- compiletroph %*% t(compiletroph)
# edit(co_occurrencetroph)
graphtroph<- graph.adjacency(co_occurrencetroph,
                            weighted=TRUE,
                            mode="undirected",
                            diag=FALSE) 

#egdge shared social parasites
co_occurrencesocial <- compilesocial %*% t(compilesocial)
# edit(co_occurrencesocial)
graphsocial<- graph.adjacency(co_occurrencesocial,
                              weighted=TRUE,
                              mode="undirected",
                              diag=FALSE)

#egdge shared fungi
co_occurrencefun <- compilefun %*% t(compilefun)
# edit(co_occurrencefun)
graphfun <- graph.adjacency(co_occurrencefun,
                            weighted=TRUE,
                            mode="undirected",
                            diag=FALSE)


  #egdge shared helminths
co_occurrencehelm<- compilehelm %*% t(compilehelm)
# edit(co_occurrencenemat)
graphhelm<- graph.adjacency(co_occurrencehelm,
                             weighted=TRUE,
                             mode="undirected",
                             diag=FALSE)

#expanding labels in circle layout
radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}
lab.locs <- radian.rescale(x=1:35, direction=-1, start=0)
colpie<-c("azure3", "goldenrod1", "hotpink","aquamarine1","dodgerblue")
valuecolpie<-list(colpie)
# edit(ball)


#FORMICINAE
        #Camponotus
      ball[1,1]<-cos(2*pi/39)
      ball[1,2]<-sin(2*pi/39)
        #Plagio
      ball[2,1]<-cos(2*2*pi/39)
      ball[2,2]<-sin(2*2*pi/39)
        #Cata
      ball[3,1]<-cos(2.65*2*pi/39)
      ball[3,2]<-sin(2.65*2*pi/39)
        #Lasius
      ball[4,1]<-cos(4*2*pi/39)
      ball[4,2]<-sin(4*2*pi/39)
        #Prenolepis
      ball[5,1]<-cos(5.2*2*pi/39)
      ball[5,2]<-sin(5.2*2*pi/39)
        #Iberoformica
      ball[6,1]<-cos(6*2*pi/39)
      ball[6,2]<-sin(6*2*pi/39)
        #Proformica
      ball[7,1]<-cos(6.7*2*pi/39)
      ball[7,2]<-sin(6.7*2*pi/39)
        #Formica
      ball[8,1]<-cos(8*2*pi/39)
      ball[8,2]<-sin(8*2*pi/39)
        #Liometopum
      ball[9,1]<-cos(11*2*pi/39)
      ball[9,2]<-sin(11*2*pi/39)
        #Linepithema
      ball[10,1]<-cos(12*2*pi/39)
      ball[10,2]<-sin(12*2*pi/39)
        #Bothr
      ball[11,1]<-cos(13*2*pi/39)
      ball[11,2]<-sin(13*2*pi/39)
        #Tapinoma
      ball[12,1]<-cos(14*2*pi/39)
      ball[12,2]<-sin(14*2*pi/39)
        #MYRMICINAE"
      ball[13,1]<-cos(29*2*pi/39)
      ball[13,2]<-sin(29*2*pi/39)
      ball[14,1]<-cos(18*2*pi/39)
      ball[14,2]<-sin(18*2*pi/39)
      ball[15,1]<-cos(19*2*pi/39)
      ball[15,2]<-sin(19*2*pi/39)
      ball[16,1]<-cos(20*2*pi/39)
      ball[16,2]<-sin(20*2*pi/39)
      ball[17,1]<-cos(21*2*pi/39)
      ball[17,2]<-sin(21*2*pi/39)
      ball[18,1]<-cos(22*2*pi/39)
      ball[18,2]<-sin(22*2*pi/39)
      ball[19,1]<-cos(23*2*pi/39)
      ball[19,2]<-sin(23*2*pi/39)
      ball[20,1]<-cos(24*2*pi/39)
      ball[20,2]<-sin(24*2*pi/39)
      ball[21,1]<-cos(25*2*pi/39)
      ball[21,2]<-sin(25*2*pi/39)
      ball[22,1]<-cos(25.8*2*pi/39)
      ball[22,2]<-sin(25.8*2*pi/39)
      ball[23,1]<-cos(27*2*pi/39)
      ball[23,2]<-sin(27*2*pi/39)
      ball[24,1]<-cos(31*2*pi/39)
      ball[24,2]<-sin(31*2*pi/39)
      ball[25,1]<-cos(28*2*pi/39)
      ball[25,2]<-sin(28*2*pi/39)
      ball[26,1]<-cos(17*2*pi/39)
      ball[26,2]<-sin(17*2*pi/39)
      ball[27,1]<-cos(30*2*pi/39)
      ball[27,2]<-sin(30*2*pi/39)
      ball[28,1]<-cos(32*2*pi/39)
      ball[28,2]<-sin(32*2*pi/39)
      ball[29,1]<-cos(33*2*pi/39)
      ball[29,2]<-sin(33*2*pi/39)
      ball[30,1]<-cos(34*2*pi/39)
      ball[30,2]<-sin(34*2*pi/39)
      ball[31,1]<-cos(35*2*pi/39)
      ball[31,2]<-sin(35*2*pi/39)

    #Ponera
      ball[32,1]<-cos(-2*2*pi/39)
      ball[32,2]<-sin(-2*2*pi/39)

diversity_per_genus<-read.table("D:/OneDrive/OneDrive - UGent/Thomasdoc/Projects/Ongoing/Meta-analyse/2020/ants_per_genus.txt",header=TRUE)  
totspecpergenus<-diversity_per_genus$species_dataset+diversity_per_genus$species_not_dataset

svg("D:/te.svg")#provide additional arguments for height, width, etc
cor.test(log(totspecpergenus),log(speciestot),method= "pearson")
q<-lm(log(totspecpergenus)~log(speciestot))
plot(log(totspecpergenus),log(speciestot),pch=16,xlab="ln(total number of species per ant genus)",ylab="ln(total number of associated symbionts)")
abline(q)
dev.off()

allprop<-round(uniquetot)/speciestot
socprop<-speciessocial2/speciessocial
myrmprop<-speciesmyrm2/speciesmyrm
trophprop<-speciestroph2/speciestroph
proptable<-t(rbind(allprop,socprop,myrmprop,trophprop))
proptable

svg("D:/examplefile6.svg")#provide additional arguments for height, width, etc
par(family = "sans") 
plot(graphall,vertex.frame.color="black",layout=ball,xlim=c(-1.2,1.2),ylim=c(-1.2,1.2), vertex.pie=pie.values,,vertex.size=sqrt(3*speciestot),vertex.pie.color=valuecolpie,vertex.shape="pie",edge.curved=0.1,rescale=F,edge.width=sqrt(E(graphall)$weight),edge.color="dimgrey",vertex.label=NA)
par(new=TRUE)
plot(graphall,vertex.frame.color="black",layout=ball,xlim=c(-1.2,1.2),ylim=c(-1.2,1.2),vertex.label.cex=(1.2/81)*(totspecpergenus-82) + 2,vertex.label.color="black",vertex.label.dist=1,vertex.label.degree=lab.locs,vertex.pie=pie.values2,vertex.pie.color=valuecolpie, vertex.size=sqrt(3*uniquetot),edge.color=NA,vertex.shape="pie",rescale=F,vertex.label=row.names(compile),vertex.label.font=4,vertex.label=F)
points(0,0)
dev.off()


  #######################################################################  
  #######################################################################
  ########## 2. HOST SPECIFICITY ########################################
  #######################################################################
  #######################################################################

datasetmetab<-datasetmeta[, -grep("spx", colnames(datasetmeta))]
datasetmetared<-datasetmetab[rowSums(abs(datasetmetab[,-1:-9]))!= 0,]

#histogram met host specificity gevonden in literatuur
p<-hist(rowSums(datasetmetared[,-1:-9]),breaks=seq(0,35,by=1),xlab="number of host species",col="grey",main="")
mean(rowSums(datasetmetared[,-1:-9])) #3.60 vs ca. 1 voor obligates in Glasier et al. 2018 zie fig 4
length(which(rowSums(datasetmetared[,-1:-9]) ==1)) #255/647 soorten bij 1 soort waargenomen
length(which(rowSums(datasetmetared[,-1:-9]) !=1)) #392 symbionten bij meer dan 1 soort waargenomen
255/(392+255) #39.4% bij 1 soort waargenomen <--> 80% bij Glasier et al. 2018
length(which(rowSums(datasetmetared[,-1:-9]) > 5)) #122/ 620 --> 19.7%
#edit(helminth)
fungired <- subset(datasetmetared, FUNCTIONAL_GROUP=="fungus")
trophobiontred<-subset(datasetmetared, FUNCTIONAL_GROUP=="trophobiont")
socialred<-subset(datasetmetared, FUNCTIONAL_GROUP=="social_parasite")
helminthred<-subset(datasetmetared, FUNCTIONAL_GROUP=="helminth")
myrmecophilered<-subset(datasetmetared, FUNCTIONAL_GROUP=="myrmecophile")
coleored<-subset(datasetmetared, TAXONOMIC_GROUP=="Coleoptera")
zygentomared<-subset(datasetmetared, TAXONOMIC_GROUP=="Zygentoma")
datfun<-fungired[,-1:-9]
dattro<-trophobiontred[,-1:-9]
datsoc<-socialred[,-1:-9]
datmyrm<-myrmecophilered[,-1:-9]
dathelm<-helminthred[,-1:-9]
datcol<-coleored[,-1:-9]
datzyg<-zygentomared[,-1:-9]

tothosts<-rowSums(datasetmetared[,-1:-9])
totgenera<-rowSums(t(compile))



con <- datasetmetared
conxx<-datasetmetared[,-1:-9]
edit(conxx)
tconxx<-t(conxx)
sums<-colSums((tconxx))
conxx$total<- sums
conxx$group<-con$FUNCTIONAL_GROUP2
levels(conxx$group)
edit(conxx)
conxx$hits<-con[,9]
edit(conxx)
conxxred<-conxx[-634,] #min mut fungus
generalist0<-glm((total)~group,data=conxxred)

generalist1<-glm(log(total)~log(hits+1),data=conxxred)
generalist2<-kruskal.test((residuals(generalist1))~group,data=conxxred)
hist(residuals(generalist1))
shapiro.test(residuals((generalist1)))
Anova(generalist2)
summary(glht(generalist2, mcp(group = "Tukey")))
hist(residuals(generalist2))
DT<-dunnTest((residuals(generalist1))~group, data=conxxred,method="bh")

DT
cldList(P.unadj ~ Comparison, data=DT,threshold  = 0.05)
row.names(phylallorderpurg2)
conxxred$total
#Plot met error bars zonder controle voor sampling effort
library(ggplot2)
levels(conxxred$group)
conxxred$group2<-factor(conxxred$group,c("social_parasite","parasitic_fungi","helminth","trophobiont","parasitoid","unspecialized_myrmecophile","specialized_myrmecophile"))
svg("D:/S2.svg")
par(family = "sans")
g0 <- ggplot(conxxred, aes(x=conxxred$group2, y=(total),fill=group2,color=group2)) +stat_summary(fun.y=mean,geom="point")+ylab("number of host species")+xlab("")
g1<-stat_summary(geom = "errorbar",width=0.2, fun.data = mean_se, position = "dodge", color=c("hotpink",rep("dodgerblue",1),"aquamarine1","goldenrod1",rep("azure3",3)))
g0+g1+ scale_fill_manual(values=c("hotpink",rep("dodgerblue",1),"aquamarine1","goldenrod1",rep("azure3",3)))+scale_color_manual(values=c("hotpink",rep("dodgerblue",1),"aquamarine1","goldenrod1",rep("azure3",3))) +theme_classic()+ theme(legend.position="none")
dev.off()



#Plot met error bars met controle voor sampling effort
library(ggplot2)
conxxred$group2<-factor(conxxred$group,c("social_parasite","parasitic_fungi","helminth","trophobiont","parasitoid","unspecialized_myrmecophile","specialized_myrmecophile"))
svg("D:/fig3.svg")
par(family = "sans")
g0 <- ggplot(conxxred, aes(x=conxxred$group2, y=residuals(generalist1),fill=group2,color=group2)) +stat_summary(fun.y=mean,geom="point")+ylab("residuals(number of host species ~ sample effort)")+xlab("")
g1<-stat_summary(geom = "errorbar",width=0.2, fun.data = mean_se, position = "dodge", color=c("hotpink",rep("dodgerblue",1),"aquamarine1","goldenrod1",rep("azure3",3)))
g0+g1+ scale_fill_manual(values=c("hotpink",rep("dodgerblue",1),"aquamarine1","goldenrod1",rep("azure3",3)))+scale_color_manual(values=c("hotpink",rep("dodgerblue",1),"aquamarine1","goldenrod1",rep("azure3",3))) +theme_classic()+ theme(legend.position="none")
--dev.off()


  ###Hier de myrmecophiles opdelen en vergelijken
    #social parasites group (xenobiosis maar 1 soort en dus weggelaten)
    myrgroup <- subset(datasetmetared, FUNCTIONAL_GROUP=="myrmecophile")
                       edit(myrgroup)
    myrgroupxx<-myrgroup[,-1:-9]
    tmyrgroupxx<-t(myrgroupxx)
    sumsmyrgroup<-colSums((tmyrgroupxx))
    myrgroupxx$totalmyr<- sumsmyrgroup
    myrgroupxx$hits<-myrgroup[,9]
    generalistsoc1<-glm(log10(myrgroupxx$totalmyr)~log10(hits+1),data=myrgroupxx)
    hist(residuals(generalistsoc1))
    shapiro.test(residuals(generalistsoc1))
    generalistsoc2<-glm(residuals(generalistsoc1)~GROUP2,data=myrgroup)
    Anova(generalistsoc2)
    
    #Box plot met controle voor sampling effort
    g0 <- ggplot(generalistsoc2,aes(x=GROUP2,y=residuals(generalistsoc1)))
    g_box<-g0 + geom_boxplot(fill="grey",colour="black")+theme_bw()
    g_box+ theme_classic()
    #Plot met error bars met controle voor sampling effort
    g_mean<-g0+stat_summary(fun.y=mean,geom="point")
    g_mean+  stat_summary(geom = "errorbar",width=0.2, fun.data = mean_se, position = "dodge")+theme_bw()+ theme_classic()

    summary(glht(generalistsoc2, mcp(groupmultisoc = "Tukey")))
    
  #######################################################################    
  #######################################################################
  ########## 3.  TAXONOMIC DISTANCE HOSTS ###############################
  #######################################################################
  #######################################################################    

datatax<-read.table("D:/OneDrive/OneDrive - UGent/Thomasdoc/Projects/Ongoing/Meta-analyse/2020/taxonomy2020.txt",header=TRUE)  #181soorten
edit(datxx)
#i. Host x symbiont matrix enkel met mierensoorten
edit(dat2nozerobis)
dat2<-datasetmeta[, -grep("spx", colnames(datasetmeta))]
dat2nozero<-dat2[rowSums(abs(dat2[,-1:-9]))!= 0,]
edit(dat2nozero)
dat2nozerobis<-dat2nozero[-634,] #min mut fungus
# edit(datxx)
datxx<-dat2nozerobis[,-1:-9][, which(colSums(dat2nozerobis[,-1:-9]) != 0)]
tdatxx<-t(datxx)

rowSums(tdatxx==0)

tdatxx<-t(datxx)
AxM<-t(tdatxx)
# edit(AxM)

AxA<-datatax[,-1]+t(datatax[,-1])
AxA<-data.frame(id=c(1:181),AxA)


denumerator <- 0                
for (i in 1:646) #646 komt van dat2nozerobis)
{ 
  denumerator[i] <- nnzero(AxA[(AxA$id %in% (which(AxM[i,] != 0))),][,-1][,(AxA$id %in% (which(AxM[i,] != 0)))])
} 
# edit(denumerator)

numerator <- 0                
for (i in 1:646) #
{ 
  numerator[i] <- sum(AxA[(AxA$id %in% (which(AxM[i,] != 0))),][,-1][,(AxA$id %in% (which(AxM[i,] != 0)))])
} 

# edit(numerator)

taxdist <- 0
for (i in 1:646) #
{ 
  taxdist[i] <- sum(AxA[(AxA$id %in% (which(AxM[i,] != 0))),][,-1][,(AxA$id %in% (which(AxM[i,] != 0)))])/nnzero(AxA[(AxA$id %in% (which(AxM[i,] != 0))),][,-1][,(AxA$id %in% (which(AxM[i,] != 0)))])
} 
#
# edit(dattaxdist)

dattaxdist<-dat2nozerobis[,1:9]
# edit(dattaxdist)
dattaxdist$taxdist<-taxdist
# edit(dattaxdist)
dattaxdistNA<-subset(dattaxdist,dattaxdist$taxdist=="NaN")
table(dattaxdistNA$FUNCTIONAL_GROUP2)
dattaxdistnoNA<-subset(dattaxdist,dattaxdist$taxdist!="NaN")
# edit(dattaxdistnoNA)
spdat<-c(rep("parasitoid",15),rep("helminth",7),rep("parasitic_fungi",5),rep("social_parasite",30),rep("specialized_myrmecophile",30),rep("unspecialized_myrmecophile",155),rep("trophobiont",12))
edit(spdat)
integ<-rep.int(1,254)
pop<-cbind(spdat,integ)
colnames(pop)<-colnames((dattaxdistnoNA[,c(4,10)]))
dattaxdistnoNAext<-rbind((dattaxdistnoNA[,c(4,10)]),pop)
# edit(dattaxdistnoNAext)
modeltaxdist<-lm((dattaxdistnoNA$taxdist) ~FUNCTIONAL_GROUP2, data = dattaxdistnoNA)
modeltaxdistext<-lm((dattaxdistnoNAext$taxdist) ~FUNCTIONAL_GROUP2, data = dattaxdistnoNAext)

shapiro.test(residuals(modeltaxdistext))
v<-kruskal.test(dattaxdistnoNA$taxdist ~FUNCTIONAL_GROUP2, data = dattaxdistnoNA)
v
levels(dattaxdistnoNA$FUNCTIONAL_GROUP2)
#BH voor letter codes
DT = dunnTest(dattaxdistnoNA$taxdist ~ FUNCTIONAL_GROUP2,data=dattaxdistnoNA,method="bh")
DT = DT$res
DT
cldList(P.adj ~ Comparison, data      = DT,threshold = 0.05)
       
        summary(modeltaxdist)


# edit(modeltaxdist)
tabletaxdist<-summary(modeltaxdist)[[4]]

edit(tabletaxdist[,2])
tabletaxdist
summary(modeltaxdist)
#Plot met error bars 
levels(dattaxdistnoNA$group)
edit(dattaxdistnoNA)
dattaxdistnoNA$group = factor(dattaxdistnoNA$FUNCTIONAL_GROUP2,levels=c("social_parasite","parasitic_fungi","helminth","trophobiont","parasitoid","unspecialized_myrmecophile","specialized_myrmecophile"))

dattaxdistnoNAext$groupext = factor(dattaxdistnoNAext$FUNCTIONAL_GROUP2,levels=c("social_parasite","parasitic_fungi","helminth","trophobiont","parasitoid","unspecialized_myrmecophile","specialized_myrmecophile"))
as.integer(dattaxdistnoNAext$taxdist)


table(dattaxdistnoNA$group)/table(dattaxdistnoNAext$groupext)

svg("D:/figS3.svg")
par(family = "sans")
ga <- ggplot(modeltaxdist,aes(x=dattaxdistnoNA$group,y=as.integer(dattaxdistnoNA$taxdist)))
ga_mean<-ga+stat_summary(fun.y=mean,geom="point",color=c("hotpink",rep("dodgerblue",1),"aquamarine1","goldenrod1",rep("azure3",3)))
ga_mean+stat_summary(geom = "errorbar",width=0.2, fun.data = mean_se, position = "dodge",color=c("hotpink",rep("dodgerblue",1),"aquamarine1","goldenrod1",rep("azure3",3))) + ylab("taxonomic distance")+xlab("")+theme_bw()+ theme_classic()
dev.off()

svg("D:/figS4.svg")
par(family = "sans")
ga <- ggplot(modeltaxdistext,aes(x=dattaxdistnoNAext$group,y=as.integer(dattaxdistnoNAext$taxdist)))
ga_mean<-ga+stat_summary(fun.y=mean,geom="point",color=c("hotpink",rep("dodgerblue",1),"aquamarine1","goldenrod1",rep("azure3",3)))
ga_mean+stat_summary(geom = "errorbar",width=0.2, fun.data = mean_se, position = "dodge",color=c("hotpink",rep("dodgerblue",1),"aquamarine1","goldenrod1",rep("azure3",3))) + ylab("taxonomic distance")+xlab("")+theme_bw()+ theme_classic()
dev.off()


    ###Hier de myrmecophiles opdelen en vergelijken
    myrgrouptax <- subset(dattaxdistnoNA, FUNCTIONAL_GROUP=="myrmecophile")
    edit(myrgrouptax)

    
     taxmyrmodel <-glm(taxdist~GROUP2,data=myrgrouptax)
    Anova(taxmyrmodel)
    
    #Box plot met controle voor sampling effort
    g0 <- ggplot(taxmyrmodel,aes(x=GROUP2,y=taxdist))
    g_box<-g0 + geom_boxplot(fill="grey",colour="black")+theme_bw()
    g_box+ theme_classic()
    #Plot met error bars met controle voor sampling effort
    g_mean<-g0+stat_summary(fun.y=mean,geom="point")
    g_mean+  stat_summary(geom = "errorbar",width=0.2, fun.data = mean_se, position = "dodge")+theme_bw()+ theme_classic()
    
    summary(glht(generalistsoc2, mcp(groupmultisoc = "Tukey")))
    

############################################################################
#phylogenetic matrix in stead of taxonomic distance matrix
edit(tdatxx)
    diversity<-read.table("D:/OneDrive/OneDrive - UGent/Thomasdoc/Projects/Ongoing/Meta-analyse/2020/anttraits2020.txt",header=TRUE) 
    row.names(diversity)<-diversity$species
    row.names(tdatxx)<-diversity$species #row.names zitten wat foutjes in, gelijk maken
    diversity$all<-rowSums(tdatxx)
    options(na.action = "na.fail") 
    diversitycompl<-diversity[complete.cases(diversity), ] #mieren waar we alle gegevens voor hebben
    datasetmetaredy<-subset(tdatxx,rownames(tdatxx) %in% rownames(diversitycompl)) #95*640
    
    
    
    #iii. Phylogenetic distance matrix based on Arman et al. 2016
    distAman<-read.table("D:/OneDrive/OneDrive - UGent/Thomasdoc/Projects/Ongoing/Meta-analyse/2020/distAman.txt",header=TRUE) 
    # edit(distAman) 155 soorten mieren
    phylall<-subset(distAman,distAman$species %in% rownames(tdatxx))
    #edit(phylall) 108 mieren met fylogenie
    tphylall<-t(phylall[,-1])
    #edit(phylall)
    phylall2<-subset(tphylall,rownames(tphylall) %in% rownames(tdatxx))
    phylallorder<-phylall2[ order(row.names(phylall2)), ] #108x108
    # edit(phylall)    
    
edit(tdatxx)
phyldatasetmeta<-subset(tdatxx,rownames(tdatxx) %in% rownames(phylallorder))
sum(datasetmetaredx)/sum(tdatxx)

sum(phyldatasetmeta)/sum(datasetmetared[,-1:-9])

edit(phyldatasetmeta)

phylnozero<-phyldatasetmeta[, which(colSums(phyldatasetmeta) != 0)] #108*616
# edit(phylnozero)
sum(rowSums(t(phylnozero))==1)

rowSums(phylnozero)
phylxx<- phylnozero[ order(row.names(phylnozero)), ]
s<-phylxx[,colSums(phylxx)!=1]

ts<-t(s)
edit(s)
library(xlsx)
write.xlsx(diversitycomplball,file="D:/y.xlsx")
AxMp<-t(phylxx)
# edit(AxMp)

AxAp<-data.frame(id=c(1:108),phylallorder)
# edit(AxAp)



phyldist <- 0
for (i in 1:623) #623 komt van phylnozero)
{ 
  phyldist[i] <- sum(AxAp[(AxAp$id %in% (which(AxMp[i,] != 0))),][,-1][,(AxAp$id %in% (which(AxMp[i,] != 0)))])/nnzero(AxAp[(AxAp$id %in% (which(AxMp[i,] != 0))),][,-1][,(AxAp$id %in% (which(AxMp[i,] != 0)))])
} 
#
edit(phyldist)
q<-subset(phyldist,phyldist=="NaN")

# edit(dat2nozero)
phylnozerotraits<-subset(dat2nozero,rownames(dat2nozero)  %in% rownames(AxMp))
# edit(phylnozerotraits)

# edit(dattaxdist)
phylnozerotraits$phyldist<-phyldist
# edit(phylnozerotraits)
phylnozerotraitsnoNA<-subset(phylnozerotraits,phylnozerotraits$phyldist!="NaN")
edit(phylnozerotraitsnoNA[335:338,])
modelphyldist<-lm(phylnozerotraitsnoNA$phyldist ~FUNCTIONAL_GROUP2, data = phylnozerotraitsnoNA)
shapiro.test(residuals(modelphyldist))
Anova(modelphyldist,test.statistic=c("LR"))
levels(phylnozerotraitsnoNA$FUNCTIONAL_GROUP2)

#Plot met error bars 
phylnozerotraitsnoNA$group = factor(phylnozerotraitsnoNA$FUNCTIONAL_GROUP2,levels=c("social_parasite","parasitic_fungi","helminth","trophobiont","parasitoid","unspecialized_myrmecophile","specialized_myrmecophile"))
phylnozerotraitsnoNA$group
svg("D:/figS4.svg")
par(family = "sans")
ga_phyl <- ggplot(modelphyldist,aes(x=phylnozerotraitsnoNA$group,y=phylnozerotraitsnoNA$phyldist))
ga_mean_phyl<-ga_phyl+stat_summary(fun.y=mean,geom="point",color=c("hotpink",rep("dodgerblue",1),"aquamarine1","goldenrod1",rep("azure3",3)))

ga_mean_phyl+stat_summary(geom = "errorbar",width=0.2, fun.data = mean_se, position = "dodge",color=c("hotpink",rep("dodgerblue",1),"aquamarine1","goldenrod1",rep("azure3",3))) +xlab("")+ylab("phylogenetic distance (node counts)")+theme_bw()+ theme_classic()
dev.off()
sumpothoc_phyl<-summary(glht(modelphyldist, mcp(FUNCTIONAL_GROUP2 = "Tukey"),test = adjusted('BH')))
cld(sumpothoc_phyl)

#???###Hier de myrmecophiles opdelen en vergelijken
#myrgroupphyl <- subset(phylnozerotraitsnoNA, FUNCTIONAL_GROUP=="myrmecophile")


#myrgroupphyl$GROUP2
#phylmyrmodel <-glm(phyldist~GROUP2,data=myrgroupphyl)
#summary(phylmyrmodel)
#Anova(phylmyrmodel)
#shapiro.test(residuals(phylmyrmodel))
#Box plot met controle voor sampling effort
#g0 <- ggplot(phylmyrmodel,aes(x=GROUP2,y=phyldist))
#g_box<-g0 + geom_boxplot(fill="grey",colour="black")+theme_bw()
#g_box+ theme_classic()
#Plot met error bars met controle voor sampling effort
#myrgroupphyl$GROUP2 = factor(myrgroupphyl$GROUP2,levels=c("Lepidoptera","Coleoptera","Coleoptera_SP", "Diptera","Hymenoptera","Araneae","Acari","Zygentoma","Hemiptera","Isopoda","Orthoptera","Collembola"))
#g_mean<-ggplot(myrgroupphyl,aes(x=GROUP2,y=phyldist))+stat_summary(fun.y=mean,geom="point")

#g_mean+  stat_summary(geom = "errorbar",width=0.2, fun.data = mean_se, position = "dodge")+theme_bw()+ theme_classic()

#summary(glht(phylmyrmodel, mcp(GROUP2 = "Tukey")))
#cld(glht(phylmyrmodel, mcp(GROUP2 = "Tukey")))


  #####################################################################
  #####################################################################
  ########## 4.  Spatial autocorrelation  #############################
  #####################################################################
  #####################################################################

  
#A. MIDPOINT DISTRIBUTION
##########################
  diversity<-read.table("D:/Postdocdata/meta-analyse/2019/anttraits2019.txt",header=TRUE)  
  coordinates<-cbind(diversity$mean_y,diversity$mean_x)
    
  eucldistancematrix<-(rdist.earth(coordinates,miles = F, R = NULL))
  hist(sqrt(eucldistancematrix))

#B. SYMPATRY
#############
  landen<-read.table("D:/OneDrive/OneDrive - UGent/Thomasdoc/Projects/Ongoing/Meta-analyse/2020/landen2020.txt",header=TRUE) 
  # edit(landen)
  Mlanden<-data.matrix(landen[,-1], rownames.force = NA)
  dim(Mlanden)
  sympatry<-matrix( rep( 0),ncol=181,nrow = 181)
  for (i in 1:181) {
    for(j in 1:181) {
      sympatry[i,j]<- sum(Mlanden[i,]*Mlanden[j,])
      print(sympatry[i,j])
    } 
  }
  
  #zelfde met functie
  # sympfunc=function(mat){
  #   for(i in 1:dim(mat)[1]){
  #     for (j in 1:dim(mat)[1]){
  #       sympatry[i,j]=sum(mat[i,]*mat[j,])
  #     }
  #   }
  #   print(sympatry)
  # }
  # 
  # smat=sympfunc(Mlanden)
  # edit(smat)
  
  
 
  
  # edit(sympatry) #sympatrie = aantal landen gemeenschappelijke distributie
  
  sympatrybin<-ifelse(sympatry > 0,1,0) #overlap in distributie 1, geen overlap 0
  sympatric_ants<-(rowSums(sympatrybin))
edit(sympatric_ants)
  ###################################################################
  ###################################################################
  ########## 5.  SIMILARITY SYMBIONT COMMUNITY ######################
  ###################################################################
  ###################################################################

#A. ALL SPECIES
#***************
diversitycomplball<-subset(diversitycompl,rownames(diversitycompl) %in% rownames(phylallorder))
x<-nrow(diversitycomplball) #96 soorten
 edit(diversitycomplball)

# edit(diversitycomplball) 
  a<-diversitycomplball$lncolonysize
  b<-diversitycomplball$habitat
  c<-diversitycomplball$nest
  d<-diversitycomplball$workersize
  e<-diversitycomplball$biogeographic_region
  f<-diversitycomplball$hitsgooglescholar
  g<-diversitycomplball$humidity
  h<-diversitycomplball$countries
    
  lncolonysize.dist<-matrix( rep( 0),ncol=x,nrow = x)
  for (i in 1:x) {
    for(j in 1:x) {
      lncolonysize.dist[i,j]<-log(abs(exp(a[i])-exp(a[j]))+1)
      print(lncolonysize.dist[i,j])
    } 
  }
  
  habitat.dist<-matrix( rep( 0),ncol=x,nrow = x)
  for (i in 1:x) {
    for(j in 1:x) {
      habitat.dist[i,j]<- ifelse(b[i]==b[j],0,1)
      print(habitat.dist[i,j])
    } 
  }
  
  nest.dist<-matrix( rep( 0),ncol=x,nrow = x)
  for (i in 1:x) {
    for(j in 1:x) {
      nest.dist[i,j]<- ifelse(c[i]==c[j],0,1)
      print(nest.dist[i,j])
    } 
  }
  
  workersize.dist<-matrix( rep( 0),ncol=x,nrow = x)
  for (i in 1:x) {
    for(j in 1:x) {
      workersize.dist[i,j]<- abs(d[i]-d[j])
      print(workersize.dist[i,j])
    } 
  }
  
  
  biogeographic_region.dist<-matrix( rep( 0),ncol=x,nrow = x)
  for (i in 1:x) {
    for(j in 1:x) {
      biogeographic_region.dist[i,j]<- ifelse(e[i]==e[j],0,1)
      print(biogeographic_region.dist[i,j])
    } 
  }
  
  lnhits.dist<-matrix( rep( 0),ncol=x,nrow = x)
  for (i in 1:x) {
    for(j in 1:x) {
      lnhits.dist[i,j]<-log(f[i]+1)*log(f[j]+1)
      print(lnhits.dist[i,j])
    } 
  }
  
  humidity.dist<-matrix( rep( 0),ncol=x,nrow = x)
  for (i in 1:x) {
    for(j in 1:x) {
      humidity.dist[i,j]<-ifelse(g[i]==g[j],0,1)
      print(humidity.dist[i,j])
    } 
  }
  
  lncountries.dist<-matrix( rep( 0),ncol=x,nrow = x)
  for (i in 1:x) {
    for(j in 1:x) {
      lncountries.dist[i,j]<-log(abs(h[i]-h[j])+1)
      print(lncountries.dist[i,j])
    } 
  }
  
  #hier met mutualistic fungi !!!!
  dat2<-datasetmeta[, -grep("spx", colnames(datasetmeta))]
  dat2nozero<-dat2[rowSums(abs(dat2[,-1:-9]))!= 0,]
  # edit(dat2nozero)
  
  datxxmut<-dat2nozero[,-1:-9][, which(colSums(dat2nozero[,-1:-9]) != 0)]
  library(xlsx)
  write.xlsx((diversitycomplballx),file="D:/p.xlsx")
  sum(datxxmut)
  tdatxxmut<-t(datxxmut)
  sum(dat2nozero[,-1:-9])
colSums(phylxx)
tdistall<-subset(tdatxxmut,rownames(tdatxxmut) %in% rownames(diversitycomplball))
datasetmetaredymut<-subset(tdatxxmut,rownames(tdatxxmut) %in% rownames(diversitycompl)) #95*640
row.names(datasetmetaredymut)
jacdistanceall<-vegdist(datasetmetaredymut,method="jaccard") #in principe foutje want ook nog symbionten met 0 interacties, maar geen invloed op distance matrix!
edit(datasetmetaredymut)
piep<-datasetmetaredymut[, which(colSums(datasetmetaredymut) != 0)]
piep<-subset(t(datasetmetaredymut),colSums(datasetmetaredymut) == 0)
edit(piep)
(match(row.names(phylallorder),row.names(datasetmetaredymut)))=="NA"

distmatrixall<-as.matrix(jacdistanceall)
jacdistanceallcor<-vegdist(phylnozero,method="jaccard") #hier alle symbionten met 0 interacties weggelaten, gelijk aan jacdistanceall

row.names(sympatry)<-diversity$species
sympall<-subset(sympatry,rownames(sympatry) %in% rownames(diversitycomplball))
tsympall<-t(sympall)
row.names(tsympall)<-row.names(tdatxxmut)
sympall2<-subset(tsympall,rownames(tsympall) %in% rownames(diversitycomplball))

phylallorderpurg<-subset(phylallorder,rownames(phylallorder) %in% rownames(diversitycomplball))
tphylallorderpurg<-t(phylallorderpurg)
row.names(tphylallorderpurg)<-row.names(phylallorder)
phylallorderpurg2<-subset(tphylallorderpurg,rownames(tphylallorderpurg) %in% rownames(diversitycomplball))
row.names(phylallorderpurg)==row.names(phylallorder)
match(row.names(phylallorderpurg),row.names(phylallorder))
# edit(sympall2)
# hist(log(sympall2))

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

d.distmatrix<-as.dist(distmatrixall)
edit(d.distmatrix)
d.countries<-as.dist(range01(lncountries.dist)) #niet opnemen, maakt niet echt biologisch zin zie ook sympatrie
d.workersize<-as.dist(range01(workersize.dist))
d.hits<-as.dist(range01(lnhits.dist))
d.biogeography<-as.dist(range01(biogeographic_region.dist))
d.phylogeny<-as.dist(range01(phylallorderpurg2))
d.colonysize<-as.dist(range01(lncolonysize.dist))
d.habitat<-as.dist(range01(habitat.dist))
d.nest<-as.dist(range01(nest.dist))
d.sympatry<-as.dist(1-range01(log(sympall2+1)))
hist(d.sympatry)
d.humidity<-as.dist(range01(humidity.dist))
hist((1-range01((sympall2))))
#MRM test
  full<-MRM(d.distmatrix~d.workersize+d.hits+d.biogeography+d.phylogeny+d.colonysize+d.habitat+d.nest+d.sympatry,nperm=9999)
  full<-MRM(d.distmatrix~d.workersize+d.hits+d.biogeography+d.phylogeny+d.colonysize+d.habitat,nperm=9999)
  
  
  #stepwise removing non-significant predictors cf Martiny et al. 2011
  full2def<-MRM(d.distmatrix~d.workersize+d.hits+d.biogeography+d.phylogeny+d.colonysize+d.habitat,nperm=9999)
  full2def
 
  
  
  full2<-lm(d.distmatrix~d.workersize+d.hits+d.biogeography+d.phylogeny+d.colonysize+d.habitat)
  vif(full2)
  plot(full2)
  library(relaimpo)
  calc.relimp(full2,type = c("lmg"),rela=TRUE)
  bootall <- boot.relimp(full2, b = 1000, type = c("lmg"), rank = TRUE, diff = TRUE, rela = TRUE)
  booteval.relimp(bootall) # print result
  plot(booteval.relimp(bootall,sort=TRUE)) 
  


  #B. MYRMECOPHILES
#############################
myrmecophile <- subset(datasetmeta, FUNCTIONAL_GROUP=="myrmecophile")
  
#edit(datasetmeta)  
 
datmyr2<-myrmecophile[, -grep("spx", colnames(myrmecophile))]
dat2myrnozero<-datmyr2[rowSums(abs(datmyr2[,-1:-9]))!= 0,]
datmyrxx<-dat2myrnozero[,-1:-9][, which(colSums(dat2myrnozero[,-1:-9]) != 0)]
tdatmyrxx<-t(datmyrxx)
row.names(tdatmyrxx)
tdatmyrxx2<-cbind(tdatmyrxx, total = rowSums(tdatmyrxx))
# edit(tdatmyrxx2)

diversitycomplbmyr<-subset(diversitycompl,rownames(diversitycompl) %in% row.names(tdatmyrxx2)) #selecteren in diversity trait set welke mieren met myrmecofielen compleet zijn
diversitycomplbmyr2<-subset(diversitycomplbmyr,rownames(diversitycomplbmyr) %in% row.names(phylallorder)) #selecteren in diversity trait set welke mieren met myrmecofielen compleet zijn
# edit(diversitycomplbmyr2)
#totaal aantal myrmecofielen bijvoegen
summyrmecophiles<-subset(tdatmyrxx2,rownames(tdatmyrxx2) %in% row.names (diversitycomplbmyr2)) 
# edit(summyrmecophiles)
diversitycomplbmyr3<-cbind(diversitycomplbmyr2, totalmyrm=summyrmecophiles[,468])   
# edit(diversitycomplbmyr3) #78 mierensoorten
#  
  a<-diversitycomplbmyr2$lncolonysize
  b<-diversitycomplbmyr2$habitat
  c<-diversitycomplbmyr2$nest
  d<-diversitycomplbmyr2$workersize
  e<-diversitycomplbmyr2$biogeographic_region
  f<-diversitycomplbmyr2$hitsgooglescholar
  g<-diversitycomplbmyr2$humidity

  myrxx<-subset(tdatmyrxx,rownames(tdatmyrxx) %in% rownames(diversitycomplbmyr2))
  colSums(myrxx)
  jacdistancemyr<-vegdist(myrxx,method="jaccard") #erg belangrijk om dit juist te zetten
  distmatrixmyr<-as.matrix(jacdistancemyr)
  # edit(distmatrixmyr)
  
  sympmyr<-subset(sympatry,rownames(sympatry) %in% rownames(diversitycomplbmyr2))
  tsympmyr<-t(sympmyr)
  row.names(tsympmyr)<-row.names(tdatxx)
  sympmyr2<-subset(tsympmyr,rownames(tsympmyr) %in% rownames(diversitycomplbmyr2))
  ncol(sympmyr2)
  
  phylmyr<-subset(distAman,distAman$species %in% rownames(diversitycomplbmyr2))
  # edit(phylmyr)
  tphylmyr<-t(phylmyr[,-1])
  phylmyr2<-subset(tphylmyr,rownames(tphylmyr) %in% rownames(diversitycomplbmyr2))
  row.names(phylmyr2)<-row.names(diversitycomplbmyr2)
  phylmyrorder<-phylmyr2[ order(row.names(phylmyr2)), ]

  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  range01(nest.dist)
  lncolonysize.dist<-matrix( rep( 0),ncol=78,nrow = 78)
  for (i in 1:78) {
    for(j in 1:78) {
      lncolonysize.dist[i,j]<-log(abs(exp(a[i])-exp(a[j]))+1)
      print(lncolonysize.dist[i,j])
    } 
  }
  
  habitat.dist<-matrix( rep( 0),ncol=78,nrow = 78)
  for (i in 1:78) {
    for(j in 1:78) {
      habitat.dist[i,j]<- ifelse(b[i]==b[j],0,1)
      print(habitat.dist[i,j])
    } 
  }
  
  nest.dist<-matrix( rep( 0),ncol=78,nrow = 78)
  for (i in 1:78) {
    for(j in 1:78) {
      nest.dist[i,j]<- ifelse(c[i]==c[j],0,1)
      print(nest.dist[i,j])
    } 
  }
  
  workersize.dist<-matrix( rep( 0),ncol=78,nrow = 78)
  for (i in 1:78) {
    for(j in 1:78) {
      workersize.dist[i,j]<- abs(d[i]-d[j])
      print(workersize.dist[i,j])
    } 
  }
  
  
  biogeographic_region.dist<-matrix( rep( 0),ncol=78,nrow = 78)
  for (i in 1:78) {
    for(j in 1:78) {
      biogeographic_region.dist[i,j]<- ifelse(e[i]==e[j],0,1)
      print(biogeographic_region.dist[i,j])
    } 
  }
  
  lnhits.dist<-matrix( rep( 0),ncol=78,nrow = 78)
  for (i in 1:78) {
    for(j in 1:78) {
      lnhits.dist[i,j]<-log(f[i]+1)*log(f[j]+1)
      print(lnhits.dist[i,j])
    } 
  }
  
  humidity.dist<-matrix( rep( 0),ncol=78,nrow = 78)
  for (i in 1:78) {
    for(j in 1:78) {
      humidity.dist[i,j]<-ifelse(g[i]==g[j],0,1)
      print(humidity.dist[i,j])
    } 
  }
  
  
  
  
  d.distmatrix<-as.dist(distmatrixmyr)
  
  d.workersize<-as.dist(range01(workersize.dist))
  d.hits<-as.dist(range01(lnhits.dist))
  d.biogeography<-as.dist(range01(biogeographic_region.dist))
  d.phylogeny<-as.dist(range01(phylmyrorder))
  d.colonysize<-as.dist(range01(lncolonysize.dist))
  d.habitat<-as.dist(range01(habitat.dist))
  d.nest<-as.dist(range01(nest.dist))
  d.sympatry<-as.dist(1-range01(log(sympmyr2+1)))
  
  d.humidity<-as.dist(range01(humidity.dist))
  
  
  # edit(d.phylogenymyr)

#MRM  
  fullmyr<-MRM(d.distmatrix~d.workersize+d.hits+d.biogeography+d.phylogeny+d.colonysize+d.habitat+d.nest+d.sympatry,nperm=9999)
  fullmyr
  
  # Int             0.947889369 0.00080008
  # d.workersize    0.031557114 0.00180018
  # d.hits         -0.073406142 0.00020002
  # d.biogeography  0.041313951 0.00010001
  # d.phylogeny    -0.010023902 0.24932493
  # d.colonysize    0.009935095 0.34703470
  # d.habitat       0.003097875 0.34363436
  # d.nest          0.009195807 0.01010101
  # d.sympatry      0.021201616 0.05220522
  # 
  fullmyrdef<-MRM(d.distmatrix~d.workersize+d.hits+d.biogeography+d.nest+d.sympatry,nperm=9999)
  

  # $`coef`
  # d.distmatrix       pval
  # Int              0.94761709 0.00020002
  # d.workersize     0.03108153 0.00190019
  # d.hits          -0.07133719 0.00050005
  # d.biogeography   0.04166273 0.00010001
  # d.nest           0.01040816 0.00410041
  # d.sympatry       0.02060916 0.05540554
  
  fullmyr3<-lm(d.distmatrix~d.workersize+d.hits+d.biogeography+d.nest+d.sympatry,nperm=9999)
  calc.relimp(fullmyr3,type=c("lmg"),rela=TRUE)
  bootmyr <- boot.relimp(fullmyr3, b = 1000, type = c("lmg"), rank = TRUE, diff = TRUE, rela = TRUE)
  booteval.relimp(bootmyr)# print result
  plot(booteval.relimp(bootmyr,sort=TRUE)) 

#C. TROPHOBIONT #including caterpillars Lepidoptera
################
  trophobiont <- subset(datasetmeta, FUNCTIONAL_GROUP2=="trophobiont")
  edit(trophobiont)
dataphid2<-trophobiont[, -grep("spx", colnames(trophobiont))]
dat2aphidnozero<-dataphid2[rowSums(abs(dataphid2[,-1:-9]))!= 0,]
dataphidxx<-dat2aphidnozero[,-1:-9][, which(colSums(dat2aphidnozero[,-1:-9]) != 0)]
tdataphidxx<-t(dataphidxx)
tdataphidxxorder<- tdataphidxx[ order(row.names(tdataphidxx)), ]  
match(row.names(tdatmyrxx),row.names(tdatmyrxxorder))
totaphid<-rowSums(tdataphidxxorder)
tdataphidxxorder2<-cbind(tdataphidxxorder, total = rowSums(tdataphidxxorder))
# edit(tdataphidxxorder2) # 68 soorten mieren met 79 trofobionts
  
diversitycomplbaphid<-subset(diversitycompl,rownames(diversitycompl) %in% row.names(tdataphidxxorder2)) #selecteren in diversity trait set welke mieren met aphids compleet zijn
diversitycomplbaphid2<-subset(diversitycomplbaphid,rownames(diversitycomplbaphid) %in% rownames(phylallorder)) #selecteren in diversity trait set welke mieren met myrmecofielen compleet zijn
# edit(diversitycomplbaphid2) 
#totaal aantal aphids bijvoegen
sumaphids<-subset(tdataphidxxorder2,rownames(tdataphidxxorder2) %in% row.names (diversitycomplbaphid2)) 
diversitycomplbaphid3<-cbind(diversitycomplbaphid2, totalaphids=sumaphids[,79]) #49 mierensoorten
  
# edit(diversitycomplbaphid3) #49 mierensoorten
a<-diversitycomplbaphid3$lncolonysize
b<-diversitycomplbaphid3$habitat
c<-diversitycomplbaphid3$nest
d<-diversitycomplbaphid3$workersize
e<-diversitycomplbaphid3$biogeographic_region
f<-diversitycomplbaphid3$hitsgooglescholar
g<-diversitycomplbaphid$humidity

phylaphid<-subset(distAman,distAman$species %in% rownames(diversitycomplbaphid3))
tphylaphid<-t(phylaphid[,-1]) #heeft al row.names
phylaphid2<-subset(tphylaphid,rownames(tphylaphid) %in% rownames(diversitycomplbaphid3))
phylaphidorder<-phylaphid2[ order(row.names(phylaphid2)), ]
edit(tdatxxorder)
aphidxxorder<-subset(tdataphidxxorder,rownames(tdataphidxxorder) %in% rownames(diversitycomplbaphid3))
jacdistanceaphid<-vegdist(aphidxxorder,method="jaccard") #erg belangrijk om dit juist te zetten
distmatrixaphid<-as.matrix(jacdistanceaphid)

sympaphid<-subset(sympatry,rownames(sympatry) %in% rownames(diversitycomplbaphid3))
sympaphid
tsympaphid<-t(sympaphid)
row.names(tsympaphid)<-row.names(tdatxx)
sympaphid2<-subset(tsympaphid,rownames(tsympaphid) %in% rownames(diversitycomplbaphid3))



lncolonysize.dist<-matrix( rep( 0),ncol=49,nrow = 49)
for (i in 1:49) {
  for(j in 1:49) {
    lncolonysize.dist[i,j]<-log(abs(exp(a[i])-exp(a[j]))+1)
    print(lncolonysize.dist[i,j])
  } 
}

habitat.dist<-matrix( rep( 0),ncol=49,nrow = 49)
for (i in 1:49) {
  for(j in 1:49) {
    habitat.dist[i,j]<- ifelse(b[i]==b[j],0,1)
    print(habitat.dist[i,j])
  } 
}

nest.dist<-matrix( rep( 0),ncol=49,nrow = 49)
for (i in 1:49) {
  for(j in 1:49) {
    nest.dist[i,j]<- ifelse(c[i]==c[j],0,1)
    print(nest.dist[i,j])
  } 
}

workersize.dist<-matrix( rep( 0),ncol=49,nrow = 49)
for (i in 1:49) {
  for(j in 1:49) {
    workersize.dist[i,j]<- abs(d[i]-d[j])
    print(workersize.dist[i,j])
  } 
}


biogeographic_region.dist<-matrix( rep( 0),ncol=49,nrow = 49)
for (i in 1:49) {
  for(j in 1:49) {
    biogeographic_region.dist[i,j]<- ifelse(e[i]==e[j],0,1)
    print(biogeographic_region.dist[i,j])
  } 
}

lnhits.dist<-matrix( rep( 0),ncol=49,nrow = 49)
for (i in 1:49) {
  for(j in 1:49) {
    lnhits.dist[i,j]<-log(f[i]+1)*log(f[j]+1)
    print(lnhits.dist[i,j])
  } 
}

humidity.dist<-matrix( rep( 0),ncol=49,nrow = 49)
for (i in 1:49) {
  for(j in 1:49) {
    humidity.dist[i,j]<-ifelse(g[i]==g[j],0,1)
    print(humidity.dist[i,j])
  } 
}
d.distmatrixaphid<-as.dist(distmatrixaphid)
d.workersizeaphid<-as.dist(range01(workersize.dist))
d.hitsaphid<-as.dist(range01(lnhits.dist))
d.humidityaphid<-as.dist(range01(humidity.dist))
d.biogeographyaphid<-as.dist(range01(biogeographic_region.dist))
d.phylogenyaphid<-as.dist((phylaphidorder))
d.colonysizeaphid<-as.dist(range01(lncolonysize.dist))
d.habitataphid<-as.dist(range01(habitat.dist))
d.nestaphid<-as.dist(range01(nest.dist))
d.sympatryaphid<-as.dist(1-range01(log(sympaphid2+1)))
hist(1-range01((sympaphid2)))
#MRM
  fullaphid<-MRM(d.distmatrixaphid~d.workersizeaphid+d.hitsaphid+d.biogeographyaphid+d.phylogenyaphid+d.colonysizeaphid+d.habitataphid+d.nestaphid+d.sympatryaphid,nperm=1000)
  fullaphid<-MRM(d.distmatrixaphid~d.hitsaphid+d.biogeographyaphid+d.phylogenyaphid,nperm=1000)
  
  fullaphiddef<-MRM(d.distmatrixaphid~d.hitsaphid+d.biogeographyaphid+d.phylogenyaphid,nperm=9999)
  


  # $`coef`
  # d.distmatrixaphid  pval
  # Int                       0.916018265 0.001
  # d.hitsaphid              -0.147222879 0.001
  # d.biogeographyaphid       0.032526434 0.001
  # d.phylogenyaphid          0.004198789 0.001
  # d.sympatryaphid           0.042966488 0.022

  fullaphid4<-lm(d.distmatrixaphid~d.hitsaphid+d.biogeographyaphid+d.phylogenyaphid)
    calc.relimp(fullaphid4,type=c("lmg"),rela=TRUE)
  bootaphid <- boot.relimp(fullaphid4, b = 1000, type = c("lmg"), rank = TRUE, diff = TRUE, rela = TRUE)
  booteval.relimp(bootaphid) # print result
  plot(booteval.relimp(bootaphid,sort=TRUE)) 


#D. SOCIAL PARASITE
####################
social <- subset(datasetmeta, FUNCTIONAL_GROUP=="social_parasite")
datsoc2<-social[, -grep("spx", colnames(social))]
dat2socnozero<-datsoc2[rowSums(abs(datsoc2[,-1:-9]))!= 0,]
datsocxx<-dat2socnozero[,-1:-9][, which(colSums(dat2socnozero[,-1:-9]) != 0)]
tdatsocxx<-t(datsocxx)
tdatsocxxorder<- tdatsocxx[ order(row.names(tdatsocxx)), ]


diversitycomplbsoc<-subset(diversitycompl,rownames(diversitycompl) %in% row.names(tdatsocxxorder)) #selecteren in diversity trait set welke mieren met myrmecofielen compleet zijn
diversitycomplbsoc2<-subset(diversitycomplbsoc,rownames(diversitycomplbsoc) %in% rownames(phylallorder)) 
edit(diversitycomplbsoc2)
##47 mierensoorten

a<-diversitycomplbsoc2$lncolonysize
b<-diversitycomplbsoc2$habitat
c<-diversitycomplbsoc2$nest
d<-diversitycomplbsoc2$workersize
e<-diversitycomplbsoc2$biogeographic_region
f<-diversitycomplbsoc2$hitsgooglescholar
g<-diversitycomplbsoc2$humidity

phylsoc<-subset(distAman,distAman$species %in% rownames(diversitycomplbsoc2))
tphylsoc<-t(phylsoc[,-1])
phylsoc2<-subset(tphylsoc,rownames(tphylsoc) %in% rownames(diversitycomplbsoc2))
phylsocorder<-phylsoc2[ order(row.names(phylsoc2)), ]


socxxorder<-subset(tdatsocxxorder,rownames(tdatsocxxorder) %in% rownames(diversitycomplbsoc2))
jacdistancesoc<-vegdist(socxxorder,method="jaccard") #erg belangrijk om dit juist te zetten
distmatrixsoc<-as.matrix(jacdistancesoc)

sympsoc<-subset(sympatry,rownames(sympatry) %in% rownames(diversitycomplbsoc2))
tsympsoc<-t(sympsoc)
row.names(tsympsoc)<-row.names(tdatxx)
sympsoc2<-subset(tsympsoc,rownames(tsympsoc) %in% rownames(diversitycomplbsoc2))
ncol(sympsoc2)


lncolonysize.dist<-matrix( rep( 0),ncol=47,nrow = 47)
for (i in 1:47) {
  for(j in 1:47) {
    lncolonysize.dist[i,j]<-log(abs(exp(a[i])-exp(a[j]))+1)
    print(lncolonysize.dist[i,j])
  } 
}

habitat.dist<-matrix( rep( 0),ncol=47,nrow = 47)
for (i in 1:47) {
  for(j in 1:47) {
    habitat.dist[i,j]<- ifelse(b[i]==b[j],0,1)
    print(habitat.dist[i,j])
  } 
}

nest.dist<-matrix( rep( 0),ncol=47,nrow = 47)
for (i in 1:47) {
  for(j in 1:47) {
    nest.dist[i,j]<- ifelse(c[i]==c[j],0,1)
    print(nest.dist[i,j])
  } 
}

workersize.dist<-matrix( rep( 0),ncol=47,nrow = 47)
for (i in 1:47) {
  for(j in 1:47) {
    workersize.dist[i,j]<- abs(d[i]-d[j])
    print(workersize.dist[i,j])
  } 
}


biogeographic_region.dist<-matrix( rep( 0),ncol=47,nrow = 47)
for (i in 1:47) {
  for(j in 1:47) {
    biogeographic_region.dist[i,j]<- ifelse(e[i]==e[j],0,1)
    print(biogeographic_region.dist[i,j])
  } 
}

lnhits.dist<-matrix( rep( 0),ncol=47,nrow = 47)
for (i in 1:47) {
  for(j in 1:47) {
    lnhits.dist[i,j]<-log(f[i]+1)*log(f[j]+1)
    print(lnhits.dist[i,j])
  } 
}

humidity.dist<-matrix( rep( 0),ncol=47,nrow = 47)
for (i in 1:47) {
  for(j in 1:47) {
    humidity.dist[i,j]<-ifelse(g[i]==g[j],0,1)
    print(humidity.dist[i,j])
  } 
}



d.distmatrixsoc<-as.dist(distmatrixsoc)
d.workersizesoc<-as.dist(range01(workersize.dist))
d.hitssoc<-as.dist(range01(lnhits.dist))
d.humiditysoc<-as.dist(range01(humidity.dist))
d.biogeographysoc<-as.dist(range01(biogeographic_region.dist))
d.phylogenysoc<-as.dist((phylsocorder))
d.colonysizesoc<-as.dist(range01(lncolonysize.dist))
d.habitatsoc<-as.dist(range01(habitat.dist))
d.nestsoc<-as.dist(range01(nest.dist))
d.sympatrysoc<-as.dist(1-range01(log(sympsoc2+1)))
# edit(d.phylogenysoc)

#MRM  
fullsoc<-MRM(d.distmatrixsoc~d.workersizesoc+d.hitssoc+d.biogeographysoc+d.phylogenysoc+d.colonysizesoc+d.habitatsoc+d.nestsoc+d.sympatrysoc,nperm=9999)

# $r.squared
fullsocdef<-MRM(d.distmatrixsoc~d.workersizesoc+d.biogeographysoc+d.phylogenysoc,nperm=9999)


fullsoc4<-lm(d.distmatrixsoc~d.workersizesoc+d.biogeographysoc+d.phylogenysoc,nperm=1000)
calc.relimp(fullsoc4,type = c("lmg"),rela=TRUE)
bootsoc <- boot.relimp(fullsoc4, b = 1000, type = c("lmg"), rank = TRUE, diff = TRUE, rela = TRUE)
plot(booteval.relimp(bootsoc,sort=TRUE)) 

#D. HELMINTH
####################
helminth <- subset(datasetmeta, FUNCTIONAL_GROUP2=="helminth")

edit(helminth)
dathel2<-helminth[, -grep("spx", colnames(helminth))]
dat2helnozero<-dathel2[rowSums(abs(dathel2[,-1:-9]))!= 0,]
dathelxx<-dat2helnozero[,-1:-9][, which(colSums(dat2helnozero[,-1:-9]) != 0)]
tdathelxx<-t(dathelxx)


diversitycomplbhel<-subset(diversitycompl,rownames(diversitycompl) %in% row.names(tdathelxx)) #selecteren in diversity trait set welke mieren met myrmecofielen compleet zijn
diversitycomplbhel2<-subset(diversitycomplbhel,rownames(diversitycomplbhel) %in% rownames(phylallorder)) 
##34

a<-diversitycomplbhel2$lncolonysize
b<-diversitycomplbhel2$habitat
c<-diversitycomplbhel2$nest
d<-diversitycomplbhel2$workersize
e<-diversitycomplbhel2$biogeographic_region
f<-diversitycomplbhel2$hitsgooglescholar
g<-diversitycomplbhel2$humidity

phylhel<-subset(distAman,distAman$species %in% rownames(diversitycomplbhel2))
tphylhel<-t(phylhel[,-1])
phylhel2<-subset(tphylhel,rownames(tphylhel) %in% rownames(diversitycomplbhel2))
phylhelorder<-phylhel2[ order(row.names(phylhel2)), ]


helxxorder<-subset(tdathelxx,rownames(tdathelxx) %in% rownames(diversitycomplbhel2))
jacdistancehel<-vegdist(helxxorder,method="jaccard") #erg belangrijk om dit juist te zetten
distmatrixhel<-as.matrix(jacdistancehel)

symphel<-subset(sympatry,rownames(sympatry) %in% rownames(diversitycomplbhel2))
tsymphel<-t(symphel)
row.names(tsymphel)<-row.names(tdatxx)
symphel2<-subset(tsymphel,rownames(tsymphel) %in% rownames(diversitycomplbhel2))
ncol(sympsoc2)


lncolonysize.dist<-matrix( rep( 0),ncol=34,nrow = 34)
for (i in 1:34) {
  for(j in 1:34) {
    lncolonysize.dist[i,j]<-log(abs(exp(a[i])-exp(a[j]))+1)
    print(lncolonysize.dist[i,j])
  } 
}

habitat.dist<-matrix( rep( 0),ncol=34,nrow = 34)
for (i in 1:34) {
  for(j in 1:34) {
    habitat.dist[i,j]<- ifelse(b[i]==b[j],0,1)
    print(habitat.dist[i,j])
  } 
}

nest.dist<-matrix( rep( 0),ncol=34,nrow = 34)
for (i in 1:34) {
  for(j in 1:34) {
    nest.dist[i,j]<- ifelse(c[i]==c[j],0,1)
    print(nest.dist[i,j])
  } 
}

workersize.dist<-matrix( rep( 0),ncol=34,nrow = 34)
for (i in 1:34) {
  for(j in 1:34) {
    workersize.dist[i,j]<- abs(d[i]-d[j])
    print(workersize.dist[i,j])
  } 
}


biogeographic_region.dist<-matrix( rep( 0),ncol=34,nrow = 34)
for (i in 1:34) {
  for(j in 1:34) {
    biogeographic_region.dist[i,j]<- ifelse(e[i]==e[j],0,1)
    print(biogeographic_region.dist[i,j])
  } 
}

lnhits.dist<-matrix( rep( 0),ncol=34,nrow = 34)
for (i in 1:34) {
  for(j in 1:34) {
    lnhits.dist[i,j]<-log(f[i]+1)*log(f[j]+1)
    print(lnhits.dist[i,j])
  } 
}

humidity.dist<-matrix( rep( 0),ncol=34,nrow = 34)
for (i in 1:34) {
  for(j in 1:34) {
    humidity.dist[i,j]<-ifelse(g[i]==g[j],0,1)
    print(humidity.dist[i,j])
  } 
}



d.distmatrixhel<-as.dist(distmatrixhel)
d.workersizehel<-as.dist(range01(workersize.dist))
d.hitshel<-as.dist(range01(lnhits.dist))
d.humidityhel<-as.dist(range01(humidity.dist))
d.biogeographyhel<-as.dist(range01(biogeographic_region.dist))
d.phylogenyhel<-as.dist((phylhelorder))
d.colonysizehel<-as.dist(range01(lncolonysize.dist))
d.habitathel<-as.dist(range01(habitat.dist))
d.nesthel<-as.dist(range01(nest.dist))
d.sympatryhel<-as.dist(1-range01(log(symphel2+1)))
# edit(d.phylogenysoc)

#MRM  
fullhel<-MRM(d.distmatrixhel~d.workersizehel+d.hitshel+d.biogeographyhel+d.phylogenyhel+d.colonysizehel+d.habitathel+d.nesthel+d.sympatryhel,nperm=9999)

fullheldef<-MRM(d.distmatrixhel~d.workersizehel+d.phylogenyhel+d.habitathel,nperm=9999)




fullhel4<-lm(d.distmatrixhel~d.workersizehel+d.habitathel+d.phylogenyhel,nperm=9999)
calc.relimp(fullhel4,type = c("lmg"),rela=TRUE)
boothel <- boot.relimp(fullhel4, b = 1000, type = c("lmg"), rank = TRUE, diff = TRUE, rela = TRUE)
plot(booteval.relimp(boothel,sort=TRUE)) 

#F. parasitic fungi
####################
funendomecophile <- subset(datasetmeta, FUNCTIONAL_GROUP2=="parasitic_fungi")
edit(funendomecophile)
datfunendo2<-funendomecophile[, -grep("spx", colnames(funendomecophile))]
dat2funendonozero<-datfunendo2[rowSums(abs(datfunendo2[,-1:-9]))!= 0,]
datfunendoxx<-dat2funendonozero[,-1:-9][, which(colSums(dat2funendonozero[,-1:-9]) != 0)]
tdatfunendoxx<-t(datfunendoxx)
match(row.names(tdatfunendoxx),row.names(tdatfunendoxxorder))
tdatfunendoxx2<-cbind(tdatfunendoxx, total = rowSums(tdatfunendoxx))
# edit(tdatfunendoxx2)

diversitycomplbfunendo<-subset(diversitycompl,rownames(diversitycompl) %in% row.names(tdatfunendoxx2)) #selecteren in diversity trait set welke mieren met funendomecofielen compleet zijn
diversitycomplbfunendo2<-subset(diversitycomplbfunendo,rownames(diversitycomplbfunendo) %in% row.names(phylallorder)) #selecteren in diversity trait set welke mieren met funendomecofielen compleet zijn
# edit(diversitycomplbfunendo2)
#totaal aantal funendomecofielen bijvoegen
sumfunendomecophiles<-subset(tdatfunendoxx2,rownames(tdatfunendoxx2) %in% row.names (diversitycomplbfunendo2)) 

#  
a<-diversitycomplbfunendo2$lncolonysize
b<-diversitycomplbfunendo2$habitat
c<-diversitycomplbfunendo2$nest
d<-diversitycomplbfunendo2$workersize
e<-diversitycomplbfunendo2$biogeographic_region
f<-diversitycomplbfunendo2$hitsgooglescholar
g<-diversitycomplbfunendo2$humidity

funendoxx<-subset(tdatfunendoxx,rownames(tdatfunendoxx) %in% rownames(diversitycomplbfunendo2))
jacdistancefunendo<-vegdist(funendoxx,method="jaccard") #erg belangrijk om dit juist te zetten
distmatrixfunendo<-as.matrix(jacdistancefunendo)
# edit(distmatrixfunendo)

sympfunendo<-subset(sympatry,rownames(sympatry) %in% rownames(diversitycomplbfunendo2))
tsympfunendo<-t(sympfunendo)
row.names(tsympfunendo)<-row.names(tdatxx)
sympfunendo2<-subset(tsympfunendo,rownames(tsympfunendo) %in% rownames(diversitycomplbfunendo2))
ncol(sympfunendo2)

phylfunendo<-subset(distAman,distAman$species %in% rownames(diversitycomplbfunendo2))
# edit(phylfunendo)
tphylfunendo<-t(phylfunendo[,-1])
phylfunendo2<-subset(tphylfunendo,rownames(tphylfunendo) %in% rownames(diversitycomplbfunendo2))
row.names(phylfunendo2)<-row.names(diversitycomplbfunendo2)
phylfunendoorder<-phylfunendo2[ order(row.names(phylfunendo2)), ]

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
range01(nest.dist)
lncolonysize.dist<-matrix( rep( 0),ncol=43,nrow = 43)
for (i in 1:43) {
  for(j in 1:43) {
    lncolonysize.dist[i,j]<-log(abs(exp(a[i])-exp(a[j]))+1)
    print(lncolonysize.dist[i,j])
  } 
}

habitat.dist<-matrix( rep( 0),ncol=43,nrow = 43)
for (i in 1:43) {
  for(j in 1:43) {
    habitat.dist[i,j]<- ifelse(b[i]==b[j],0,1)
    print(habitat.dist[i,j])
  } 
}

nest.dist<-matrix( rep( 0),ncol=43,nrow = 43)
for (i in 1:43) {
  for(j in 1:43) {
    nest.dist[i,j]<- ifelse(c[i]==c[j],0,1)
    print(nest.dist[i,j])
  } 
}

workersize.dist<-matrix( rep( 0),ncol=43,nrow = 43)
for (i in 1:43) {
  for(j in 1:43) {
    workersize.dist[i,j]<- abs(d[i]-d[j])
    print(workersize.dist[i,j])
  } 
}


biogeographic_region.dist<-matrix( rep( 0),ncol=43,nrow = 43)
for (i in 1:43) {
  for(j in 1:43) {
    biogeographic_region.dist[i,j]<- ifelse(e[i]==e[j],0,1)
    print(biogeographic_region.dist[i,j])
  } 
}

lnhits.dist<-matrix( rep( 0),ncol=43,nrow = 43)
for (i in 1:43) {
  for(j in 1:43) {
    lnhits.dist[i,j]<-log(f[i]+1)*log(f[j]+1)
    print(lnhits.dist[i,j])
  } 
}

humidity.dist<-matrix( rep( 0),ncol=43,nrow = 43)
for (i in 1:43) {
  for(j in 1:43) {
    humidity.dist[i,j]<-ifelse(g[i]==g[j],0,1)
    print(humidity.dist[i,j])
  } 
}




d.distmatrix<-as.dist(distmatrixfunendo)

d.workersize<-as.dist(range01(workersize.dist))
d.hits<-as.dist(range01(lnhits.dist))
d.biogeography<-as.dist(range01(biogeographic_region.dist))
d.phylogeny<-as.dist(range01(phylfunendoorder))
d.colonysize<-as.dist(range01(lncolonysize.dist))
d.habitat<-as.dist(range01(habitat.dist))
d.nest<-as.dist(range01(nest.dist))
d.sympatry<-as.dist(1-range01(log(sympfunendo2+1)))

d.humidity<-as.dist(range01(humidity.dist))


# edit(d.phylogenyfunendo)

#MRM  
fullfunendo<-MRM(d.distmatrix~d.workersize+d.hits+d.biogeography+d.phylogeny+d.colonysize+d.habitat+d.nest+d.sympatry,nperm=9999)

fullfunendodef<-MRM(d.distmatrix~d.workersize+d.hits+d.biogeography+d.colonysize,nperm=9999)


fullfunendo3<-lm(d.distmatrix~d.workersize+d.hits+d.biogeography+d.colonysize,nperm=9999)
calc.relimp(fullfunendo3,type=c("lmg"),rela=TRUE)
bootfunendo <- boot.relimp(fullfunendo3, b = 1000, type = c("lmg"), rank = TRUE, diff = TRUE, rela = TRUE)
booteval.relimp(bootfunendo)# print result
plot(booteval.relimp(bootfunendo,sort=TRUE)) 

#D. MYRMECOPHILES unspecialized
####################
myrnidmecophile <- subset(datasetmeta,FUNCTIONAL_GROUP2=="unspecialized_myrmecophile" )
edit(myrnidmecophile)
datmyrnid2<-myrnidmecophile[, -grep("spx", colnames(myrnidmecophile))]
dat2myrnidnozero<-datmyrnid2[rowSums(abs(datmyrnid2[,-1:-9]))!= 0,]
datmyrnidxx<-dat2myrnidnozero[,-1:-9][, which(colSums(dat2myrnidnozero[,-1:-9]) != 0)]
tdatmyrnidxx<-t(datmyrnidxx)
row.names(tdatmyrnidxx)
tdatmyrnidxxorder<- tdatmyrnidxx[ order(row.names(tdatmyrnidxx)), ]
tdatmyrnidxxorder2<-cbind(tdatmyrnidxxorder, total = rowSums(tdatmyrnidxxorder))
# edit(tdatmyrnidxxorder2)

diversitycomplbmyrnid<-subset(diversitycompl,rownames(diversitycompl) %in% row.names(tdatmyrnidxxorder2)) #selecteren in diversity trait set welke mieren met myrnidmecofielen compleet zijn
diversitycomplbmyrnid2<-subset(diversitycomplbmyrnid,rownames(diversitycomplbmyrnid) %in% row.names(phylallorder)) #selecteren in diversity trait set welke mieren met myrnidmecofielen compleet zijn
# edit(diversitycomplbmyrnid2)
#totaal aantal myrnidmecofielen bijvoegen
summyrnidmecophiles<-subset(tdatmyrnidxxorder2,rownames(tdatmyrnidxxorder2) %in% row.names (diversitycomplbmyrnid2)) 

#  
a<-diversitycomplbmyrnid2$lncolonysize
b<-diversitycomplbmyrnid2$habitat
c<-diversitycomplbmyrnid2$nest
d<-diversitycomplbmyrnid2$workersize
e<-diversitycomplbmyrnid2$biogeographic_region
f<-diversitycomplbmyrnid2$hitsgooglescholar
g<-diversitycomplbmyrnid2$humidity

myrnidxxorder<-subset(tdatmyrnidxxorder,rownames(tdatmyrnidxxorder) %in% rownames(diversitycomplbmyrnid2))
colSums(myrnidxxorder)
jacdistancemyrnid<-vegdist(myrnidxxorder,method="jaccard") #erg belangrijk om dit juist te zetten
distmatrixmyrnid<-as.matrix(jacdistancemyrnid)
# edit(distmatrixmyrnid)

sympmyrnid<-subset(sympatry,rownames(sympatry) %in% rownames(diversitycomplbmyrnid2))
tsympmyrnid<-t(sympmyrnid)
row.names(tsympmyrnid)<-row.names(tdatxx)
sympmyrnid2<-subset(tsympmyrnid,rownames(tsympmyrnid) %in% rownames(diversitycomplbmyrnid2))
ncol(sympmyrnid2)

phylmyrnid<-subset(distAman,distAman$species %in% rownames(diversitycomplbmyrnid2))
# edit(phylmyrnid)
tphylmyrnid<-t(phylmyrnid[,-1])
phylmyrnid2<-subset(tphylmyrnid,rownames(tphylmyrnid) %in% rownames(diversitycomplbmyrnid2))
row.names(phylmyrnid2)<-row.names(diversitycomplbmyrnid2)
phylmyrnidorder<-phylmyrnid2[ order(row.names(phylmyrnid2)), ]

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
range01(nest.dist)
lncolonysize.dist<-matrix( rep( 0),ncol=69,nrow = 69)
for (i in 1:69) {
  for(j in 1:69) {
    lncolonysize.dist[i,j]<-log(abs(exp(a[i])-exp(a[j]))+1)
    print(lncolonysize.dist[i,j])
  } 
}

habitat.dist<-matrix( rep( 0),ncol=69,nrow = 69)
for (i in 1:69) {
  for(j in 1:69) {
    habitat.dist[i,j]<- ifelse(b[i]==b[j],0,1)
    print(habitat.dist[i,j])
  } 
}

nest.dist<-matrix( rep( 0),ncol=69,nrow = 69)
for (i in 1:69) {
  for(j in 1:69) {
    nest.dist[i,j]<- ifelse(c[i]==c[j],0,1)
    print(nest.dist[i,j])
  } 
}

workersize.dist<-matrix( rep( 0),ncol=69,nrow = 69)
for (i in 1:69) {
  for(j in 1:69) {
    workersize.dist[i,j]<- abs(d[i]-d[j])
    print(workersize.dist[i,j])
  } 
}


biogeographic_region.dist<-matrix( rep( 0),ncol=69,nrow = 69)
for (i in 1:69) {
  for(j in 1:69) {
    biogeographic_region.dist[i,j]<- ifelse(e[i]==e[j],0,1)
    print(biogeographic_region.dist[i,j])
  } 
}

lnhits.dist<-matrix( rep( 0),ncol=69,nrow = 69)
for (i in 1:69) {
  for(j in 1:69) {
    lnhits.dist[i,j]<-log(f[i]+1)*log(f[j]+1)
    print(lnhits.dist[i,j])
  } 
}

humidity.dist<-matrix( rep( 0),ncol=69,nrow = 69)
for (i in 1:69) {
  for(j in 1:69) {
    humidity.dist[i,j]<-ifelse(g[i]==g[j],0,1)
    print(humidity.dist[i,j])
  } 
}




d.distmatrix<-as.dist(distmatrixmyrnid)

d.workersize<-as.dist(range01(workersize.dist))
d.hits<-as.dist(range01(lnhits.dist))
d.biogeography<-as.dist(range01(biogeographic_region.dist))
d.phylogeny<-as.dist(range01(phylmyrnidorder))
d.colonysize<-as.dist(range01(lncolonysize.dist))
d.habitat<-as.dist(range01(habitat.dist))
d.nest<-as.dist(range01(nest.dist))
d.sympatry<-as.dist(1-range01(log(sympmyrnid2+1)))

d.humidity<-as.dist(range01(humidity.dist))


# edit(d.phylogenymyrnid)

#MRM  
fullmyrnid<-MRM(d.distmatrix~d.workersize+d.hits+d.biogeography+d.phylogeny+d.colonysize+d.habitat+d.nest+d.sympatry,nperm=9999)
fullmyrnid

# Int             0.947889369 0.00080008
# d.workersize    0.031557114 0.00180018
# d.hits         -0.073406142 0.00020002
# d.biogeography  0.041313951 0.00010001
# d.phylogeny    -0.010023902 0.24932493
# d.colonysize    0.009935095 0.34703470
# d.habitat       0.003097875 0.34363436
# d.nest          0.009195807 0.01010101
# d.sympatry      0.021201616 0.05220522
# 

fullmyrniddef<-MRM(d.distmatrix~d.biogeography+d.nest+d.sympatry,nperm=9999)

# $`coef`
# d.distmatrix       pval
# Int              0.94761709 0.00020002
# d.workersize     0.03108153 0.00190019
# d.hits          -0.07133719 0.00050005
# d.biogeography   0.04166273 0.00010001
# d.nest           0.01040816 0.00410041
# d.sympatry       0.02060916 0.05690569

fullmyrnid3<-lm(d.distmatrix~d.nest+d.biogeography+d.sympatry,nperm=999)
calc.relimp(fullmyrnid3,type=c("lmg"),rela=TRUE)
bootmyrnid <- boot.relimp(fullmyrnid3, b = 1000, type = c("lmg"), rank = TRUE, diff = TRUE, rela = TRUE)
booteval.relimp(bootmyrnid)# print result
plot(booteval.relimp(bootmyrnid,sort=TRUE)) 
graph2ppt(file="D:/Postdocdata/meta-analyse/figuren.pptx", aspectr=1,append=T)

#D. MYRMECOPHILES unspecialized
####################
myrnidmecophile <- subset(datasetmeta,FUNCTIONAL_GROUP2=="specialized_myrmecophile" )
edit(myrnidmecophile)
datmyrnid2<-myrnidmecophile[, -grep("spx", colnames(myrnidmecophile))]
dat2myrnidnozero<-datmyrnid2[rowSums(abs(datmyrnid2[,-1:-9]))!= 0,]
datmyrnidxx<-dat2myrnidnozero[,-1:-9][, which(colSums(dat2myrnidnozero[,-1:-9]) != 0)]
tdatmyrnidxx<-t(datmyrnidxx)
row.names(tdatmyrnidxx)
tdatmyrnidxxorder<- tdatmyrnidxx[ order(row.names(tdatmyrnidxx)), ]
tdatmyrnidxxorder2<-cbind(tdatmyrnidxxorder, total = rowSums(tdatmyrnidxxorder))
# edit(tdatmyrnidxxorder2)

diversitycomplbmyrnid<-subset(diversitycompl,rownames(diversitycompl) %in% row.names(tdatmyrnidxxorder2)) #selecteren in diversity trait set welke mieren met myrnidmecofielen compleet zijn
diversitycomplbmyrnid2<-subset(diversitycomplbmyrnid,rownames(diversitycomplbmyrnid) %in% row.names(phylallorder)) #selecteren in diversity trait set welke mieren met myrnidmecofielen compleet zijn
# edit(diversitycomplbmyrnid2)
#totaal aantal myrnidmecofielen bijvoegen
summyrnidmecophiles<-subset(tdatmyrnidxxorder2,rownames(tdatmyrnidxxorder2) %in% row.names (diversitycomplbmyrnid2)) 

#  
a<-diversitycomplbmyrnid2$lncolonysize
b<-diversitycomplbmyrnid2$habitat
c<-diversitycomplbmyrnid2$nest
d<-diversitycomplbmyrnid2$workersize
e<-diversitycomplbmyrnid2$biogeographic_region
f<-diversitycomplbmyrnid2$hitsgooglescholar
g<-diversitycomplbmyrnid2$humidity

myrnidxxorder<-subset(tdatmyrnidxxorder,rownames(tdatmyrnidxxorder) %in% rownames(diversitycomplbmyrnid2))
colSums(myrnidxxorder)
jacdistancemyrnid<-vegdist(myrnidxxorder,method="jaccard") #erg belangrijk om dit juist te zetten
distmatrixmyrnid<-as.matrix(jacdistancemyrnid)
# edit(distmatrixmyrnid)

sympmyrnid<-subset(sympatry,rownames(sympatry) %in% rownames(diversitycomplbmyrnid2))
tsympmyrnid<-t(sympmyrnid)
row.names(tsympmyrnid)<-row.names(tdatxx)
sympmyrnid2<-subset(tsympmyrnid,rownames(tsympmyrnid) %in% rownames(diversitycomplbmyrnid2))
ncol(sympmyrnid2)

phylmyrnid<-subset(distAman,distAman$species %in% rownames(diversitycomplbmyrnid2))
# edit(phylmyrnid)
tphylmyrnid<-t(phylmyrnid[,-1])
phylmyrnid2<-subset(tphylmyrnid,rownames(tphylmyrnid) %in% rownames(diversitycomplbmyrnid2))
row.names(phylmyrnid2)<-row.names(diversitycomplbmyrnid2)
phylmyrnidorder<-phylmyrnid2[ order(row.names(phylmyrnid2)), ]

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
range01(nest.dist)
lncolonysize.dist<-matrix( rep( 0),ncol=54,nrow = 54)
for (i in 1:54) {
  for(j in 1:54) {
    lncolonysize.dist[i,j]<-log(abs(exp(a[i])-exp(a[j]))+1)
    print(lncolonysize.dist[i,j])
  } 
}

habitat.dist<-matrix( rep( 0),ncol=54,nrow = 54)
for (i in 1:54) {
  for(j in 1:54) {
    habitat.dist[i,j]<- ifelse(b[i]==b[j],0,1)
    print(habitat.dist[i,j])
  } 
}

nest.dist<-matrix( rep( 0),ncol=54,nrow = 54)
for (i in 1:54) {
  for(j in 1:54) {
    nest.dist[i,j]<- ifelse(c[i]==c[j],0,1)
    print(nest.dist[i,j])
  } 
}

workersize.dist<-matrix( rep( 0),ncol=54,nrow = 54)
for (i in 1:54) {
  for(j in 1:54) {
    workersize.dist[i,j]<- abs(d[i]-d[j])
    print(workersize.dist[i,j])
  } 
}


biogeographic_region.dist<-matrix( rep( 0),ncol=54,nrow = 54)
for (i in 1:54) {
  for(j in 1:54) {
    biogeographic_region.dist[i,j]<- ifelse(e[i]==e[j],0,1)
    print(biogeographic_region.dist[i,j])
  } 
}

lnhits.dist<-matrix( rep( 0),ncol=54,nrow = 54)
for (i in 1:54) {
  for(j in 1:54) {
    lnhits.dist[i,j]<-log(f[i]+1)*log(f[j]+1)
    print(lnhits.dist[i,j])
  } 
}

humidity.dist<-matrix( rep( 0),ncol=54,nrow = 54)
for (i in 1:54) {
  for(j in 1:54) {
    humidity.dist[i,j]<-ifelse(g[i]==g[j],0,1)
    print(humidity.dist[i,j])
  } 
}




d.distmatrix<-as.dist(distmatrixmyrnid)

d.workersize<-as.dist(range01(workersize.dist))
d.hits<-as.dist(range01(lnhits.dist))
d.biogeography<-as.dist(range01(biogeographic_region.dist))
d.phylogeny<-as.dist(range01(phylmyrnidorder))
d.colonysize<-as.dist(range01(lncolonysize.dist))
d.habitat<-as.dist(range01(habitat.dist))
d.nest<-as.dist(range01(nest.dist))
d.sympatry<-as.dist(1-range01(log(sympmyrnid2+1)))

d.humidity<-as.dist(range01(humidity.dist))


# edit(d.phylogenymyrnid)

#MRM  
fullmyrnid<-MRM(d.distmatrix~d.workersize+d.hits+d.biogeography+d.phylogeny+d.colonysize+d.habitat+d.nest+d.sympatry,nperm=9999)
fullmyrnid


fullmyrniddef<-MRM(d.distmatrix~d.workersize+d.hits+d.biogeography+d.colonysize,nperm=9999)

# $`coef`
# d.distmatrix       pval
# Int              0.94761709 0.00020002
# d.workersize     0.03108153 0.00190019
# d.hits          -0.07133719 0.00050005
# d.biogeography   0.04166273 0.00010001
# d.nest           0.01040816 0.00410041
# d.sympatry       0.02060916 0.05690569

fullmyrnid3<-lm(d.distmatrix~d.workersize+d.hits+d.colonysize+d.biogeography,nperm=999)
calc.relimp(fullmyrnid3,type=c("lmg"),rela=TRUE)
bootmyrnid <- boot.relimp(fullmyrnid3, b = 1000, type = c("lmg"), rank = TRUE, diff = TRUE, rela = TRUE)
booteval.relimp(bootmyrnid)# print result
plot(booteval.relimp(bootmyrnid,sort=TRUE)) 
graph2ppt(file="D:/Postdocdata/meta-analyse/figuren.pptx", aspectr=1,append=T)


#D. MYRMECOPHILES  endoparasitoids
####################
myrendomecophile <- subset(datasetmeta, FUNCTIONAL_GROUP2=="endoparasitic_myrmecophile")
edit(myrendomecophile)
datmyrendo2<-myrendomecophile[, -grep("spx", colnames(myrendomecophile))]
dat2myrendonozero<-datmyrendo2[rowSums(abs(datmyrendo2[,-1:-9]))!= 0,]
datmyrendoxx<-dat2myrendonozero[,-1:-9][, which(colSums(dat2myrendonozero[,-1:-9]) != 0)]
tdatmyrendoxx<-t(datmyrendoxx)
row.names(tdatmyrendoxx)
tdatmyrendoxxorder<- tdatmyrendoxx[ order(row.names(tdatmyrendoxx)), ]
tdatmyrendoxxorder2<-cbind(tdatmyrendoxxorder, total = rowSums(tdatmyrendoxxorder))
# edit(tdatmyrendoxxorder2)

diversitycomplbmyrendo<-subset(diversitycompl,rownames(diversitycompl) %in% row.names(tdatmyrendoxxorder2)) #selecteren in diversity trait set welke mieren met myrendomecofielen compleet zijn
diversitycomplbmyrendo2<-subset(diversitycomplbmyrendo,rownames(diversitycomplbmyrendo) %in% row.names(phylallorder)) #selecteren in diversity trait set welke mieren met myrendomecofielen compleet zijn
# edit(diversitycomplbmyrendo2)
#totaal aantal myrendomecofielen bijvoegen
summyrendomecophiles<-subset(tdatmyrendoxxorder2,rownames(tdatmyrendoxxorder2) %in% row.names (diversitycomplbmyrendo2)) 

#  
a<-diversitycomplbmyrendo2$lncolonysize
b<-diversitycomplbmyrendo2$habitat
c<-diversitycomplbmyrendo2$nest
d<-diversitycomplbmyrendo2$workersize
e<-diversitycomplbmyrendo2$biogeographic_region
f<-diversitycomplbmyrendo2$hitsgooglescholar
g<-diversitycomplbmyrendo2$humidity

myrendoxxorder<-subset(tdatmyrendoxxorder,rownames(tdatmyrendoxxorder) %in% rownames(diversitycomplbmyrendo2))
colSums(myrendoxxorder)
jacdistancemyrendo<-vegdist(myrendoxxorder,method="jaccard") #erg belangrijk om dit juist te zetten
distmatrixmyrendo<-as.matrix(jacdistancemyrendo)
# edit(distmatrixmyrendo)

sympmyrendo<-subset(sympatry,rownames(sympatry) %in% rownames(diversitycomplbmyrendo2))
tsympmyrendo<-t(sympmyrendo)
row.names(tsympmyrendo)<-row.names(tdatxx)
sympmyrendo2<-subset(tsympmyrendo,rownames(tsympmyrendo) %in% rownames(diversitycomplbmyrendo2))
ncol(sympmyrendo2)

phylmyrendo<-subset(distAman,distAman$species %in% rownames(diversitycomplbmyrendo2))
# edit(phylmyrendo)
tphylmyrendo<-t(phylmyrendo[,-1])
phylmyrendo2<-subset(tphylmyrendo,rownames(tphylmyrendo) %in% rownames(diversitycomplbmyrendo2))
row.names(phylmyrendo2)<-row.names(diversitycomplbmyrendo2)
phylmyrendoorder<-phylmyrendo2[ order(row.names(phylmyrendo2)), ]

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
range01(nest.dist)
lncolonysize.dist<-matrix( rep( 0),ncol=26,nrow = 26)
for (i in 1:26) {
  for(j in 1:26) {
    lncolonysize.dist[i,j]<-log(abs(exp(a[i])-exp(a[j]))+1)
    print(lncolonysize.dist[i,j])
  } 
}

habitat.dist<-matrix( rep( 0),ncol=26,nrow = 26)
for (i in 1:26) {
  for(j in 1:26) {
    habitat.dist[i,j]<- ifelse(b[i]==b[j],0,1)
    print(habitat.dist[i,j])
  } 
}

nest.dist<-matrix( rep( 0),ncol=26,nrow = 26)
for (i in 1:26) {
  for(j in 1:26) {
    nest.dist[i,j]<- ifelse(c[i]==c[j],0,1)
    print(nest.dist[i,j])
  } 
}

workersize.dist<-matrix( rep( 0),ncol=26,nrow = 26)
for (i in 1:26) {
  for(j in 1:26) {
    workersize.dist[i,j]<- abs(d[i]-d[j])
    print(workersize.dist[i,j])
  } 
}


biogeographic_region.dist<-matrix( rep( 0),ncol=26,nrow = 26)
for (i in 1:26) {
  for(j in 1:26) {
    biogeographic_region.dist[i,j]<- ifelse(e[i]==e[j],0,1)
    print(biogeographic_region.dist[i,j])
  } 
}

lnhits.dist<-matrix( rep( 0),ncol=26,nrow = 26)
for (i in 1:26) {
  for(j in 1:26) {
    lnhits.dist[i,j]<-log(f[i]+1)*log(f[j]+1)
    print(lnhits.dist[i,j])
  } 
}

humidity.dist<-matrix( rep( 0),ncol=26,nrow = 26)
for (i in 1:26) {
  for(j in 1:26) {
    humidity.dist[i,j]<-ifelse(g[i]==g[j],0,1)
    print(humidity.dist[i,j])
  } 
}




d.distmatrix<-as.dist(distmatrixmyrendo)

d.workersize<-as.dist(range01(workersize.dist))
d.hits<-as.dist(range01(lnhits.dist))
d.biogeography<-as.dist(range01(biogeographic_region.dist))
d.phylogeny<-as.dist(range01(phylmyrendoorder))
d.colonysize<-as.dist(range01(lncolonysize.dist))
d.habitat<-as.dist(range01(habitat.dist))
d.nest<-as.dist(range01(nest.dist))
d.sympatry<-as.dist(1-range01(log(sympmyrendo2+1)))

d.humidity<-as.dist(range01(humidity.dist))


# edit(d.phylogenymyrendo)

#MRM  
fullmyrendo<-MRM(d.distmatrix~d.workersize+d.hits+d.biogeography+d.phylogeny+d.colonysize+d.habitat+d.nest+d.sympatry,nperm=9999)
fullmyrendo

# Int             0.942689369 0.00080008
# d.workersize    0.031557114 0.00180018
# d.hits         -0.073406142 0.00020002
# d.biogeography  0.041313951 0.00010001
# d.phylogeny    -0.010023902 0.24932493
# d.colonysize    0.009935095 0.34703470
# d.habitat       0.003097875 0.34363436
# d.nest          0.009195807 0.01010101
# d.sympatry      0.021201616 0.05220522
# 

fullmyrendodef<-MRM(d.distmatrix~d.workersize+d.colonysize,nperm=9999)

# $`coef`
# d.distmatrix       pval
# Int              0.94761709 0.00020002
# d.workersize     0.03108153 0.00190019
# d.hits          -0.07133719 0.00050005
# d.biogeography   0.04166273 0.00010001
# d.nest           0.01040816 0.00410041
# d.sympatry       0.02060916 0.05690554

fullmyrendo3<-lm(d.distmatrix~d.workersize+d.colonysize,nperm=999)
calc.relimp(fullmyrendo3,type=c("lmg"),rela=TRUE)
bootmyrendo <- boot.relimp(fullmyrendo3, b = 1000, type = c("lmg"), rank = TRUE, diff = TRUE, rela = TRUE)
booteval.relimp(bootmyrendo)# print result
plot(booteval.relimp(bootmyrendo,sort=TRUE)) 

#####
#PLOTS
print(booteval.relimp(bootall,sort=TRUE), width=25)
print(booteval.relimp(bootaphid,sort=TRUE))
print(booteval.relimp(bootsoc,sort=TRUE)) 
print(booteval.relimp(bootmyr,sort=TRUE)) 
print(booteval.relimp(boothel,sort=TRUE)) 
print(booteval.relimp(bootfunendo,sort=TRUE)) 


xf<-c((calc.relimp(full2,type=c("lmg"),rela=TRUE))$namen[-1])
xf<-c("d worker size","d sample effort","d biogeography","d phylogeny","d colony size","d habitat")
yf<-c((booteval.relimp(bootall))$lmg)
lf<-c((booteval.relimp(bootall))$lmg.lower)
uf<-c((booteval.relimp(bootall))$lmg.upper)
dff<-data.frame(xf,yf,lf,uf)

xa<-c((calc.relimp(fullaphid4,type=c("lmg"),rela=TRUE))$namen[-1])
xa<-c("d sample effort","d biogeography","d phylogeny")
ya<-c((booteval.relimp(bootaphid))$lmg)
la<-c((booteval.relimp(bootaphid))$lmg.lower)
ua<-c((booteval.relimp(bootaphid))$lmg.upper)
dfa<-data.frame(xa,ya,la,ua)
dum1<-data.frame("",0,0,0)
dum2<-data.frame(".",0,0,0)
dum3<-data.frame("..",0,0,0)
names(dum1)<-names(dfa)
names(dum2)<-names(dfa)
names(dum3)<-names(dfa)
dfa<-rbind(dfa,dum1,dum2,dum3)

xs<-c((calc.relimp(fullsoc4,type=c("lmg"),rela=TRUE))$namen[-1])
xs<-c("d worker size","d biogeography", "d phylogeny")
ys<-c((booteval.relimp(bootsoc))$lmg)
ls<-c((booteval.relimp(bootsoc))$lmg.lower)
us<-c((booteval.relimp(bootsoc))$lmg.upper)
dfs<-data.frame(xs,ys,ls,us)
dum1<-data.frame("",0,0,0)
dum2<-data.frame(".",0,0,0)
names(dum1)<-names(dfs)
names(dum2)<-names(dfs)
dfs<-rbind(dfs,dum1,dum2)

xm<-c((calc.relimp(fullmyr3,type=c("lmg"),rela=TRUE))$namen[-1])
xm<-c("d worker size", "d sample effort","d biogeography","d nest","1-sympatry")
ym<-c((booteval.relimp(bootmyr))$lmg)
lm<-c((booteval.relimp(bootmyr))$lmg.lower)
um<-c((booteval.relimp(bootmyr))$lmg.upper)
dfm<-data.frame(xm,ym,lm,um)
dum1<-data.frame("",0,0,0)
names(dum1)<-names(dfm)
dfm<-rbind(dfm,dum1)


xh<-c((calc.relimp(fullhel4,type=c("lmg"),rela=TRUE))$namen[-1])
xh<-c("d worker size","d habitat","d phylogeny")
yh<-c((booteval.relimp(boothel))$lmg)
lh<-c((booteval.relimp(boothel))$lmg.lower)
uh<-c((booteval.relimp(boothel))$lmg.upper)
dfh<-data.frame(xh,yh,lh,uh)
dum1<-data.frame("",0,0,0)
dum2<-data.frame(".",0,0,0)
dum3<-data.frame("..",0,0,0)
names(dum1)<-names(dfh)
names(dum2)<-names(dfh)
names(dum3)<-names(dfh)
dfhel<-rbind(dfh,dum1,dum2,dum3)


xfun<-c((calc.relimp(fullfunendo3,type=c("lmg"),rela=TRUE))$namen[-1])
xfun<-c("d worker size","d sample effort","d biogeography","d colony size")
yfun<-c((booteval.relimp(bootfunendo))$lmg)
lfun<-c((booteval.relimp(bootfunendo))$lmg.lower)
ufun<-c((booteval.relimp(bootfunendo))$lmg.upper)
dfun<-data.frame(xfun,yfun,lfun,ufun)
dum1<-data.frame("",0,0,0)
dum2<-data.frame(".",0,0,0)
names(dum1)<-names(dfun)
names(dum2)<-names(dfun)
dffun<-rbind(dfun,dum1,dum2)

dff$xf <- factor(dff$xf, levels = dff$xf[order(-dff$yf)])
f<-ggplot(data=dff, aes(dff$xf, y=dff$yf,rela=TRUE)) +geom_bar(stat="identity",fill=c("dodgerblue","azure3","aquamarine1","hotpink", "brown","goldenrod1"))+geom_errorbar(aes(ymin=lf, ymax=uf), width=.2,position=position_dodge(.9))
ff<-f+xlab("")+ylab("R²")+theme_bw()+ theme_classic()
ff

dfa$xa<-factor(dfa$xa, levels = dfa$xa[order(-dfa$ya)])
a<-ggplot(data=dfa, aes(dfa$xa,, y=dfa$ya,rela=TRUE)) +geom_bar(stat="identity",fill=c("aquamarine1","dodgerblue","azure3",NA,NA,NA))+geom_errorbar(aes(ymin=dfa$la, ymax=dfa$ua), width=.2,position=position_dodge(.9))
aa<-a+xlab("")+ylab("R²")+theme_bw()+ theme_classic()
aa

dfs$xs<-factor(dfs$xs, levels = dfs$xs[order(-dfs$ys)])
s<-ggplot(data=dfs, aes(xs, y=ys,rela=TRUE)) +geom_bar(stat="identity",fill=c("dodgerblue","hotpink","azure3",NA,NA))+geom_errorbar(aes(ymin=ls, ymax=us), width=.2,position=position_dodge(.9))
ss<-s+xlab("")+ylab("R²")+theme_bw()+ theme_classic()
ss

dfm$xm<-factor(dfm$xm, levels = dfm$xm[order(-dfm$ym)])
m<-ggplot(data=dfm, aes(xm, y=ym,rela=TRUE)) +geom_bar(stat="identity",fill=c("azure3","purple","azure4","aquamarine1","hot pink",NA))+geom_errorbar(aes(ymin=lm, ymax=um), width=.2,position=position_dodge(.9))+scale_colour_gradient2()
mm<-m+xlab("")+ylab("R²")+theme_bw()+ theme_classic()
mm

dfhel$xh<-factor(dfhel$xh, levels = dfhel$xh[order(-dfhel$yh)])
h<-ggplot(data=dfhel, aes(xh, y=yh,rela=TRUE)) +geom_bar(stat="identity",fill=c("dodgerblue","hotpink","goldenrod1",NA,NA,NA))+geom_errorbar(aes(ymin=lh, ymax=uh), width=.2,position=position_dodge(.9))
hh<-h+xlab("")+ylab("R²")+theme_bw()+ theme_classic()
hh

dffun$xfun<-factor(dffun$xfun, levels = dffun$xfun[order(-dffun$yfun)])
fun<-ggplot(data=dffun,aes( xfun, y=yfun,rela=TRUE)) +geom_bar(stat="identity",fill=c("azure3","aquamarine1","hotpink","brown",NA,NA))+geom_errorbar(aes(ymin=lfun, ymax=ufun), width=.2,position=position_dodge(.9))
funfun<-fun+xlab("")+ylab("R²")+theme_bw()+ theme_classic()
funfun


q<-c(signif.num(c(summary(fullsoc4))$coefficients[-1,4]))
q

signif.num <- function(x) {
  symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
         cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
         symbols = c("***", "**", "*", ".", " "))
}


library(ggpubr)
library(grid)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

plot1 <-ff+xlab("")+ylab("Proportion of R²")+theme_bw()+ theme_classic()+ggtitle("all symbionts")+theme(axis.line.x=element_blank(),axis.ticks.x=element_blank())+scale_y_continuous(limits=c(0.0,1))+theme(axis.text.x = element_text(angle = 60, hjust = 1))+annotation_custom(grid.text("R² = 11.9 %", x=0.7,  y=0.8))
plot2 <- aa+xlab("")+ylab("Proportion of R²")+theme_bw()+ theme_classic()+ggtitle("trophobionts")+theme(axis.line.x=element_blank(),axis.ticks.x=element_blank())+scale_y_continuous(limits=c(0.0,1))+theme(axis.text.x = element_text(angle = 60, hjust = 1))+annotation_custom(grid.text("R² = 14.9 %", x=0.7,  y=0.8)) 
plot3 <-ss+xlab("")+ylab("Proportion of R²")+theme_bw()+ theme_classic()+ggtitle("social parasites")+theme(axis.line.x=element_blank(),axis.ticks.x=element_blank())+scale_y_continuous(limits=c(0.0,1))+theme(axis.text.x = element_text(angle = 60, hjust = 1))+annotation_custom(grid.text("R² = 22.6 %", x=0.7,  y=0.8))
plot4 <- mm+xlab("")+ylab("Proportion of R²")+theme_bw()+ theme_classic()+ggtitle("myrmecophiles")+theme(axis.line.x=element_blank(),axis.ticks.x=element_blank())+scale_y_continuous(limits=c(0.0,1))+theme(axis.text.x = element_text(angle = 60, hjust = 1))+annotation_custom(grid.text("R² = 13.0 %", x=0.7,  y=0.8))
plot5 <-hh+xlab("")+ylab("Proportion of R²")+theme_bw()+ theme_classic()+ggtitle("helminths")+theme(axis.line.x=element_blank(),axis.ticks.x=element_blank())+scale_y_continuous(limits=c(0.0,1))+theme(axis.text.x = element_text(angle = 60, hjust = 1))+annotation_custom(grid.text("R² = 19.1 %", x=0.7,  y=0.8))
plot6 <- funfun+xlab("")+ylab("Proportion of R²")+theme_bw()+ theme_classic()+ggtitle("parasitic fungi")+theme(axis.line.x=element_blank(),axis.ticks.x=element_blank())+scale_y_continuous(limits=c(0.0,1))+theme(axis.text.x = element_text(angle = 60, hjust = 1))+annotation_custom(grid.text("R² = 13.3 %", x=0.7,  y=0.8))
grid.arrange(plot1, plot4,plot2,plot5,plot3,plot6, nrow=3)+theme(axis.line.x=element_blank(),axis.ticks.x=element_blank())+scale_y_continuous(limits=c(0.0,0.7))

svg("D:/examplefilfig5.svg")
par(family = "sans")
grid.arrange(plot1, plot4,plot2,plot5,plot3,plot6, nrow=3)
dev.off()


grid.arrange(plot1,plot2,plot4,plot3, plot5,plot6)

  #######################################################################
  ########## 6.  DIVERSITY ~ ANT TRAITS #################################
  #######################################################################
  #######################################################################

diversity<-read.table("D:/OneDrive/OneDrive - UGent/Thomasdoc/Projects/Ongoing/Meta-analyse/2020/anttraits2020.txt",header=TRUE) 
row.names(diversity)<-diversity$species
diversity$all<-rowSums(tdatxxmut)
edit(datxxmut)
#diversity$totalmyrm<-colSums(myrmecophile[,-1:-9]) 
#diversity$totalhelm<-colSums(helminth[,-1:-9]) 
#diversity$totalfun<-colSums(fun[,-1:-9]) 
#diversity$totaltrophobionts<-colSums(trophobiont[,-1:-9]) 
#diversity$totalsoc<-colSums(datsoc2[,-1:-9]) 




diversity$sympatric_ants<-sympatric_ants

# edit(diversity)
diversitycompl<-diversity[complete.cases(diversity), ] #96 mieren
#Min exotic species
#???diversitycomplb <- subset(diversitycompl, habitat!="exotic")

datasetmetaredx<-subset(tdatxxmut,rownames(tdatxxmut) %in% rownames(diversitycompl))
# edit(dataset)
sum(datasetmetaredx)/sum(tdatxx) # 85.5 % interacties gecovered
sum(phyldatasetmeta)/sum(tdatxxmut)
sum(colSums(datasetmetaredx)==0)
                                            ########
                                            ##PGLS##
                                            ########

#A. PGLS total diversity
########################
library(ape)
ar<-read.nexus("D:/OneDrive/OneDrive - UGent/Thomasdoc/Projects/Ongoing/Meta-analyse/2018/Arnan.nex")
tiplabels<-data.frame(ar$tip.label)

ar$tip.label
# edit(diversitycomplb)
#selecteren in diversity trait set welke mieren fylogenie hebben
#95 mieren --> tree aanpassen in mesquite
diversitycomplball<-subset(diversitycompl,rownames(diversitycompl) %in% tiplabels[,1]) #
# edit(diversitycomplball)
treeall<-read.nexus("D:/OneDrive/OneDrive - UGent/Thomasdoc/Projects/Ongoing/Meta-analyse/2020/treeall2020")
# edit(treeall)
diversitycomplballx<-diversitycomplball

#diversitycomplballx<-diversitycomplball[-56,] #min Linepithema
sympatric_ants
attach(diversitycomplballx)
    countries_stand<-scale(sqrt(countries))
    hist((countries_stand))
    lncolonysize_stand<-scale(lncolonysize)
    logsympatric_ants_stand<-scale(log(diversitycomplballx$sympatric_ants))
    loghits_stand<-scale(log(hitsgooglescholar.230119.))
    workersize_stand<-scale((workersize))
    hist((workersize))
    habitat<-factor(habitat)
    nest<-factor(nest)
    biogeographic_region<-factor(biogeographic_region)
library(caper)
      shorebirdall <- comparative.data(treeall, diversitycomplballx, species, vcv=TRUE)
  mod <- pgls((all)^0.5 ~ countries_stand+lncolonysize_stand+logsympatric_ants_stand+loghits_stand+workersize_stand+habitat+nest+biogeographic_region,  lambda='ML', shorebirdall)
  shapiro.test(residuals(mod)) #ok !!!
  plot(mod)
  mod.l <- pgls.profile(modred, 'lambda')
  plot(mod.l)
  modred<- pgls((all)^0.5 ~countries_stand+lncolonysize_stand+logsympatric_ants_stand+loghits_stand+workersize_stand,  lambda='ML', shorebirdall)
  anova.pgls.fixed(mod)
  mod_rank<-dredge(mod,rank = "AICc")
  mod_sel<-subset(mod_rank,delta < 2)
  library(xlsx)
  write.xlsx(mod_sel,file="D:/test.xlsx")
  anova(mod)
  anova.pgls.fixed(get.models(mod_sel, subset =4)[[1]])
  hist(residuals(get.models(mod_sel, subset = 1)[[1]]))
  shapiro.test(residuals(get.models(mod_sel, subset = 1)[[1]]))
  edit(diversitycompl)
  pgls.profile(get.models(mod_sel, subset = 3)[[1]], 'lambda')
write.xlsx(mod_sel,"D:/Postdocdata/meta-analyse/2019/table.xlsx")

modbest <- pgls((all)^0.5 ~countries_stand+habitat+lncolonysize_stand+loghits_stand,  lambda='ML', shorebirdall)
summary(modbest)
modbestred <- pgls((all)^0.5 ~countries_stand+habitat+lncolonysize_stand,  lambda='ML', shorebirdall)

AICc(modbest)-AICc(modbestred)
anova(modbest)
summary(modbest)
(residuals(modbest))


#B. PGLS myrmecophile diversity
################################

tdatmyrxx2<-cbind(tdatmyrxx, total = rowSums(tdatmyrxx))
# edit(tdatmyrxx2) # 125 soorten mieren met myrmecofielen
  
diversitycomplbmyr<-subset(diversitycompl,rownames(diversitycompl) %in% row.names(tdatmyrxx2)) #selecteren in diversity trait set welke mieren met myrmecofielen compleet zijn
diversitycomplbmyr2<-subset(diversitycomplbmyr,rownames(diversitycomplbmyr) %in% rownames(phylallorder)) #selecteren in diversity trait set welke mieren met myrmecofielen compleet zijn
# edit(diversitycomplbmyr2)
# edit(summyrmecophiles)
dropmyr<-c(which(!(treeall$tip.label%in%rownames(diversitycomplbmyr2))))
treemyr<-drop.tip(treeall,treeall$tip.label[dropmyr])

#totaal aantal myrmecofielen bijvoegen
summyrmecophiles<-subset(tdatmyrxx2,rownames(tdatmyrxx2) %in% row.names (diversitycomplbmyr2)) 
# edit(summyrmecophiles)
diversitycomplbmyr3<-cbind(diversitycomplbmyr2, totalmyrm=summyrmecophiles[,468]) #71 soorten

attach(diversitycomplbmyr3)
# edit(diversitycomplbmyr3)
  
  countries_stand<-scale(sqrt(diversitycomplbmyr3$countries))
  lncolonysize_stand<-scale(diversitycomplbmyr3$lncolonysize)
  logsympatric_ants_stand<-scale(log(diversitycomplbmyr3$sympatric_ants))
  loghits_stand<-scale(log(diversitycomplbmyr3$hitsgooglescholar))
  workersize_stand<-scale((diversitycomplbmyr3$workersize))
  hist(sympatric_ants)
  habitat<-factor(diversitycomplbmyr3$habitat)
  nest<-factor(diversitycomplbmyr3$nest)
  biogeographic_region<-factor(diversitycomplbmyr3$biogeographic_region)

shorebird <- comparative.data(treemyr, diversitycomplbmyr3, species, vcv=TRUE)
modmyr <- pgls((totalmyrm)^0.5 ~ countries_stand+lncolonysize_stand+logsympatric_ants_stand+loghits_stand+workersize_stand+habitat+nest+biogeographic_region,  lambda='ML', shorebird)
summary(modmyr)

modmyr_rank<-dredge(modmyr,rank = "AICc")
modmyr_sel<-subset(modmyr_rank,delta < 2)
write.xlsx(modmyr_sel,"D:/table.xlsx")
summary(get.models(modmyr_sel, subset = 4)[[1]])
(pgls.profile(get.models(modmyr_sel, subset = 3)[[1]], 'lambda'))
anova.pgls.fixed(get.models(modmyr_sel, subset =5)[[1]])


modmyrbest <- pgls((totalmyrm)^0.5 ~ workersize_stand+countries_stand+lncolonysize_stand+loghits_stand+habitat,  lambda='ML', shorebird)
modmyrbestred <- pgls((totalmyrm)^0.5~ workersize_stand+countries_stand+lncolonysize_stand+loghits_stand,  lambda='ML', shorebird)
AICc(modmyrbest)-AICc(modmyrbestred)
hist(residuals(modmyrbest))
shapiro.test(residuals(modmyrbest)) #plus minus correct 0.049

#B.1 PGLS Helminth diversity
################################
tdathelxx
totcol<-rowSums(tdathelxx)
tdatcolxx2<-cbind(tdathelxx, total = rowSums(tdathelxx))
# edit(tdatcolxxorder2) 

diversitycomplbcol<-subset(diversitycompl,rownames(diversitycompl) %in% row.names(tdatcolxx2)) #selecteren in diversity trait set welke mieren met colmecofielen compleet zijn
diversitycomplbcol2<-subset(diversitycomplbcol,rownames(diversitycomplbcol) %in% rownames(phylallorder)) #selecteren in diversity trait set welke mieren met colmecofielen compleet zijn
# edit(diversitycomplbcol2)
# edit(sumcolmecophiles)
dropcol<-c(which(!(treeall$tip.label%in%rownames(diversitycomplbcol2))))
treecol<-drop.tip(treeall,treeall$tip.label[dropcol])

#totaal aantal colmecofielen bijvoegen
sumcolmecophiles<-subset(tdatcolxx2,rownames(tdatcolxx2) %in% row.names (diversitycomplbcol2)) 
diversitycomplbcol3<-cbind(diversitycomplbcol2, totalcol=sumcolmecophiles[,22]) #58 soorten
# edit(diversitycomplbcol3)
attach(diversitycomplbcol3)


countries_stand<-scale(sqrt(diversitycomplbcol3$countries))
lncolonysize_stand<-scale(diversitycomplbcol3$lncolonysize)
logsympatric_ants_stand<-scale(log(diversitycomplbcol3$sympatric_ants))
loghits_stand<-scale(log(diversitycomplbcol3$hitsgooglescholar))
workersize_stand<-scale(diversitycomplbcol3$workersize)

habitat<-factor(diversitycomplbcol3$habitat)
nest<-factor(diversitycomplbcol3$nest)
biogeographic_region<-factor(diversitycomplbcol3$biogeographic_region)

  shorebird <- comparative.data(treecol, diversitycomplbcol3, species, vcv=TRUE)
  modcol <- pgls((totalcol)^0.5 ~ countries_stand+lncolonysize_stand+logsympatric_ants_stand+loghits_stand+workersize_stand+habitat+nest+biogeographic_region,  lambda='ML', shorebird)
summary(modcol)
  modcol_rank<-dredge(modcol,rank = "AICc")
  modcol_sel<-subset(modcol_rank,delta < 2)
  write.xlsx(modcol_sel,"D:/table.xlsx")
  anova.pgls.fixed(get.models(modcol_sel, subset = 5)[[1]])
  
  
  modcolbest <- pgls((totalcol)^0.5 ~ habitat+loghits_stand,  lambda='ML', shorebird)
  AICc(modcolbest)
  modcolred <- pgls((totalcol)^0.5 ~  habitat+loghits_stand+lncolonysize_stand,  lambda='ML', shorebird)
  AICc(modcolbest)-AICc(modcolred)
  hist(residuals(modcol))
  shapiro.test(residuals(modcol))
  
  modbest <- pgls((totalcol)^0.5 ~ habitat+loghits_stand,  lambda='ML', shorebird)
  summary(modbest)
  anova(modbest)
  mod.lbest <- pgls.profile(modbest, 'lambda')
  plot(mod.lbest)
  shapiro.test(residuals(modcol))
  
  anova.pgls.fixed(get.models(modcol_sel, subset =8)[[1]])
  
  pgls.profile(get.models(modcol_sel, subset = 8)[[1]], 'lambda')
  
plot(diversitycomplbcol3$habitat,diversitycomplbcol3$totalcol)
  
  #B.1 PGLS Fungi diversity
  ################################  

totfungi<-rowSums(tdatfunendoxx)
tdatfungixx2<-cbind(tdatfunendoxx, total = totfungi)

diversitycomplbfungi<-subset(diversitycomplb,rownames(diversitycomplb) %in% row.names(tdatfungixx2)) #selecteren in diversity trait set welke mieren met colmecofielen compleet zijn
diversitycomplbfungi2<-subset(diversitycomplbfungi,rownames(diversitycomplbfungi) %in% rownames(phylallorder)) #selecteren in diversity trait set welke mieren met colmecofielen compleet zijn
# edit(diversitycomplbfungi2)
# edit(sumcolmecophiles)
dropfungi<-c(which(!(treeall$tip.label%in%rownames(diversitycomplbfungi2))))
treefungi<-drop.tip(treeall,treeall$tip.label[dropfungi])

#totaal aantal  bijvoegen
sumfungi<-subset(tdatfungixx2,rownames(tdatfungixx2) %in% row.names (diversitycomplbfungi2)) 
diversitycomplbfungi3<-cbind(diversitycomplbfungi2, totalfungi=sumfungi[,171]) #36 soorten
# edit(diversitycomplbcol3)
attach(diversitycomplbfungi3)

countries_stand<-scale(sqrt(diversitycomplbfungi3$countries))
lncolonysize_stand<-scale(diversitycomplbfungi3$lncolonysize)
logsympatric_ants_stand<-scale(log(diversitycomplbfungi3$sympatric_ants))
loghits_stand<-scale(log(diversitycomplbfungi3$hitsgooglescholar))
workersize_stand<-scale(diversitycomplbfungi3$workersize)

habitat<-factor(diversitycomplbfungi3$habitat)
nest<-factor(diversitycomplbfungi3$nest)
biogeographic_region<-factor(diversitycomplbfungi3$biogeographic_region)

shorebird <- comparative.data(treefungi, diversitycomplbfungi3, species, vcv=TRUE)
modfungi <- pgls((totalfungi)^0.25 ~ countries_stand+lncolonysize_stand+logsympatric_ants_stand+loghits_stand+workersize_stand+habitat+nest+biogeographic_region,  lambda='ML', shorebird)
summary(mod)
anova(mod)
plot(mod)

modfungi_rank<-dredge(modfungi,rank = "AICc")
modfungi_sel<-subset(modfungi_rank,delta < 2)
write.xlsx(modfungi_sel,"D:/Postdocdata/meta-analyse/2019/table.xlsx")
anova(get.models(modfungi_sel, subset = 1)[[1]])


modfungibest <- pgls((totalfungi)^0.25 ~ lncolonysize_stand+loghits_stand+habitat+countries_stand,  lambda='ML', shorebird)
AICc(modfungibest)
modfungibestred <- pgls((totalfungi)^0.25 ~ lncolonysize_stand+loghits_stand+habitat+countries_stand,  lambda='ML', shorebird)
AICc(modfungibest)-AICc(modfungibestred)
hist(residuals(modfungi))
shapiro.test(residuals(modfungibest))


  
#C. PGLS aphid diversity
#########################
totaphid<-rowSums(tdataphidxx)
tdataphidxx2<-cbind(tdataphidxx, total = rowSums(tdataphidxx))
# edit(tdataphidxx2) # 68 soorten mieren met trofobionts

diversitycomplbaphid<-subset(diversitycompl,rownames(diversitycompl) %in% row.names(tdataphidxx2)) #selecteren in diversity trait set welke mieren met aphids compleet zijn
diversitycomplbaphid2<-subset(diversitycomplbaphid,rownames(diversitycomplbaphid) %in% rownames(phylallorder)) #selecteren in diversity trait set welke mieren met myrmecofielen compleet zijn
# edit(diversitycomplbaphid2) 
dropaphid<-c(which(!(treeall$tip.label%in%rownames(diversitycomplbaphid2))))
treeaphid<-drop.tip(treeall,treeall$tip.label[dropaphid])

#totaal aantal aphids bijvoegen
sumaphids<-subset(tdataphidxx2,rownames(tdataphidxx2) %in% row.names (diversitycomplbaphid2)) 
diversitycomplbaphid3<-cbind(diversitycomplbaphid2, totalaphids=sumaphids[,79]) #40 mierensoorten
edit(sumaphids)
attach(diversitycomplbaphid3)
  
  countries_stand<-scale(sqrt(diversitycomplbaphid3$countries)) #hier square root ipv llncolonysize_stand<-scale(diversitycomplbaphid3$lncolonysize)
  logsympatric_ants_stand<-scale(log(diversitycomplbaphid3$sympatric_ants))
  loghits_stand<-scale(log(diversitycomplbaphid3$hitsgooglescholar))
  workersize_stand<-scale(diversitycomplbaphid3$workersize)
  lncolonysize_stand<-scale(diversitycomplbaphid3$lncolonysize)
  
  habitat<-factor(diversitycomplbaphid3$habitat)
  nest<-factor(diversitycomplbaphid3$nest)
  biogeographic_region<-factor(diversitycomplbaphid3$biogeographic_region)

  shorebirdA <- comparative.data(treeaphid, diversitycomplbaphid3, species, vcv=TRUE)
  modaphids <- pgls((totalaphids)^0.5 ~ countries_stand+lncolonysize_stand+logsympatric_ants_stand+loghits_stand+workersize_stand+habitat+biogeographic_region,  lambda='ML', shorebirdA)
  
  modaphids_rank<-dredge(modaphids,rank = "AICc")
  modaphids_sel<-subset(modaphids_rank,delta < 2)
  write.xlsx(modaphids_sel,"D:/table.xlsx")
  anova.pgls.fixed(get.models(modaphids_sel, subset =4)[[1]])
  pgls.profile(get.models(modaphids_sel, subset = 4)[[1]], 'lambda')
  
  
  modaphidsbest <- pgls((totalaphids)^0.5 ~habitat+loghits_stand,  lambda='ML', shorebirdA)
  AICc(modaphidsbest)
  hist(residuals(modaphidsbest))
  modaphidsbestred <- pgls((totalaphids)^0.5 ~biogeographic_region+habitat,  lambda='ML', shorebirdA)
  AICc(modaphidsbest)-AICc(modaphidsbestred)
  shapiro.test(residuals(modaphids))
  
  
#D. PGLS social parasite diversity
###################################
  
  totsoc<-rowSums(tdatsocxx)
  tdatsocxx2<-cbind(tdatsocxx, total = rowSums(tdatsocxx))
  # edit(tdatsocxxorder2) # 79 soorten mieren met socs
  # edit(tdatxxorder)
  diversitycomplbsoc<-subset(diversitycompl,rownames(diversitycompl) %in% row.names(tdatsocxx2)) #44 mierenselecteren in diversity trait set welke mieren met socs compleet zijn
  diversitycomplbsoc2<-subset(diversitycomplbsoc,rownames(diversitycomplbsoc) %in% row.names(phylallorder)) #selecteren in diversity trait set welke mieren met myrmecofielen compleet zijn
  # edit(diversitycomplbsoc2) 
  
  dropsoc<-c(which(!(treeall$tip.label%in%rownames(diversitycomplbsoc2))))
  treesoc<-drop.tip(treeall,treeall$tip.label[dropsoc])
  
  
  #totaal aantal socs bijvoegen
  sumsocs<-subset(tdatsocxx2,rownames(tdatsocxx2) %in% row.names (diversitycomplbsoc2)) 
  # edit(sumsocs)
  diversitycomplbsoc3<-cbind(diversitycomplbsoc2, totalsocs=sumsocs[,68]) #47 mierensoorten
  # edit(diversitycomplbsoc3)
  attach(diversitycomplbsoc3)
  
  countries_stand<-scale(sqrt(diversitycomplbsoc3$countries))
  lncolonysize_stand<-scale(diversitycomplbsoc3$lncolonysize)
  logsympatric_ants_stand<-scale(log(diversitycomplbsoc3$sympatric_ants))
  loghits_stand<-scale(log(diversitycomplbsoc3$hitsgooglescholar))
  workersize_stand<-scale(diversitycomplbsoc3$workersize)
  
  habitat<-factor(diversitycomplbsoc3$habitat)
  nest<-factor(diversitycomplbsoc3$nest)
  biogeographic_region<-factor(diversitycomplbsoc3$biogeographic_region)

  shorebird <- comparative.data(treesoc, diversitycomplbsoc3, species, vcv=TRUE)
  modsoc <- pgls((totalsocs)^0.5 ~ countries_stand+lncolonysize_stand+logsympatric_ants_stand+loghits_stand+workersize_stand+habitat+nest+biogeographic_region,  lambda='ML', shorebird)
  summary(modsoc)
 sqrt(totalsocs)
  modsoc_rank<-dredge(modsoc,rank = "AICc")
  modsoc_sel<-subset(modsoc_rank,delta < 2)
  write.xlsx(modsoc_sel,"D:/table.xlsx")
  (get.models(modsoc_sel, subset = 1)[[1]])
  
  hist((totalsocs)^0.25)
  modsocbest <- pgls((totalsocs)^0.5~ countries_stand,  lambda='ML', shorebird)
summary(modsocbest)
  modsocbestred <- pgls((totalsocs)^0.5 ~ countries_stand+nest,  lambda='ML', shorebird)
  AICc(modsocbest)-AICc(modsocbestred)
  shapiro.test(residuals(modsoc))
  
  anova.pgls.fixed(get.models(modsoc_sel, subset =7)[[1]])
  
  pgls.profile(get.models(modsoc_sel, subset = 7)[[1]], 'lambda')
  
  #####"AIC plot
AICtable<-read.table("D:/OneDrive/OneDrive - UGent/Thomasdoc/Projects/Ongoing/Meta-analyse/2020/AICtable2020.txt",header=TRUE) 


library(plyr)
AICtable$trait<-revalue(AICtable$trait, c("loghits"="sample effort", "lncolonysize"="colony size","logsympatric_ants"="sympatric ants","countries"="distribution range","workersize"="worker size","biogeographic_region"="biogeographic region"))
edit(AICtable)
AICall<-subset(AICtable,AICtable$type=="all_symbionts")
AICmyr<-subset(AICtable,AICtable$type=="myrmecophiles")
AIChelminth<-subset(AICtable,AICtable$type=="helminth")
AICaphid<-subset(AICtable,AICtable$type=="trophobionts")
AICsoc<-subset(AICtable,AICtable$type=="social_parasites")
dev.off()

AICall$trait <- factor(AICall$trait, levels = AICall$trait[order(AICall$diffAICc)])
AICmyr$trait <- factor(AICmyr$trait, levels = AICmyr$trait[order(AICmyr$diffAICc)])
AIChelminth$trait <- factor(AIChelminth$trait, levels = AIChelminth$trait[order(AIChelminth$diffAICc)])
AICaphid$trait <- factor(AICaphid$trait, levels = AICaphid$trait[order(AICaphid$diffAICc)])
AICsoc$trait <- factor(AICsoc$trait, levels = AICsoc$trait[order(AICsoc$diffAICc)])






P1<-ggplot()  +geom_abline(intercept = 0, slope = 0, color="black", size=0.5)+ geom_point(data = AICall, aes(trait, y=diffAICc),stat = "identity") +coord_flip()+ylab("??AICc")+xlab("")+theme_bw()+ggtitle("all symbionts")
P2<-ggplot()+geom_abline(intercept = 0, slope = 0, color="black", size=0.5)  +  geom_point(data = AICmyr, aes(x=trait, y=diffAICc),stat = "identity") +coord_flip()+ylab("??AICc")+xlab("")+theme_bw()+ggtitle("myrmecophiles")
P3<-ggplot() +geom_abline(intercept = 0, slope = 0, color="black", size=0.5) +  geom_point(data = AIChelminth, aes(x=trait, y=diffAICc),stat = "identity") +coord_flip()+ylab("??AICc")+xlab("")+theme_bw()+ggtitle("helminths")
P5<-ggplot() +geom_abline(intercept = 0, slope = 0, color="black", size=0.5) +  geom_point(data = AICaphid, aes(x=trait, y=diffAICc),stat = "identity") +coord_flip()+ylab("??AICc")+xlab("")+theme_bw()+ggtitle("trophobionts")
P6<-ggplot() +geom_abline(intercept = 0, slope = 0, color="black", size=0.5) +  geom_point(data = AICsoc, aes(x=trait, y=diffAICc),stat = "identity") +coord_flip()+ylab("??AICc")+xlab("")+theme_bw()+ggtitle("social parasites")
library(gridExtra)

svg("D:/peo.svg")
par(family = "sans")
grid.arrange(P1, P2,P5,P3,P6, nrow=3)
dev.off()

grid.arrange(P1,P2,P5,P3,P6, nrow=3)+theme(axis.line.x=element_blank(),axis.ticks.x=element_blank())+scale_y_continuous(limits=c(0.0,0.7))


#--------------------------
  
vfull<-dat2nozero[, which(colSums(dat2nozero[,-1:-9]) != 0)]
edit(vfull)
edit(phylallorder)
 