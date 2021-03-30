#pop gen of hermaphrodite "moms" only

library("polysat")

test<-read.table("C:/Users/andre/Documents/Documents/Documents/Icerya/Icerya_hem_full_geno_for_polysat.txt", header=T)

STRdata <- read.STRand("C:/Users/andre/Documents/Documents/Documents/Icerya/Icerya_hem_full_geno_for_polysat.txt")


###using full genos
library("data.table")

#what's our most common genotype?
test<-fread("C:/Users/andre/Documents/Documents/Documents/Icerya/Icerya_12gen_countrypop.txt",header=T)
df<-test[,3:14]
df<-test[,c(1,3:14)]
test[duplicated(df),]

count.duplicates <- function(DF){
x <- do.call('paste', c(DF, sep = '\r'))
  ox <- order(x)
  rl <- rle(x[ox])
  cbind(DF[ox[cumsum(rl$lengths)],,drop=FALSE],count = rl$lengths)
}

count.duplicates(df)
# 183/183 171/171 173/173  99/99 131/131 176/176 159/159 107/107 126/126 215/215 110/110 121/121   count = 93
# 183/183 171/171 173/173  99/99 134/134 176/176 159/159 107/107 126/126 215/215 110/110 121/121   count = 96 


#looking for evidence of linkage between loci
#based on one pedigree 4 loci (plcp37, plp39, plcp62 plpc93) are completely linked, and one (plcp58) is seriously linked
#how true is this overall?
table(test$PIcp37)
#206/215 215/215 
#      2     207 
table(test$PIcp39)
#165/171  171/171 
#      5      204
table(test$PIcp62)
#173/173 173/176 176/176 
#      1       1     207
table(test$PIcp93)
#115/115 115/121 121/121 124/124 
#      1       1     206       1 
#stepping out this individual
#        Pop    Ind  PIcp21  PIcp39  PIcp45 PIcp71  PIcp58  PIcp62  PIcp75 PIcp101  PIcp29  PIcp37   Picp68  PIcp93
# California CCS1-3 183/183 165/171 173/179  99/99 131/131 173/173 159/159  98/107 126/126 206/215  107/107 115/115

table(test$PIcp37)
#206/215 215/215 
#      1     207 
table(test$PIcp39)
#165/171  171/171 
#      4      204
table(test$PIcp62)
# 173/176 176/176 
#       1     207
table(test$PIcp93)
#115/121 121/121 124/124 
#      1     206       1 

#stepping out this one:
# Pop Ind     PIcp21  PIcp39  PIcp45 PIcp71  PIcp58  PIcp62  PIcp75 PIcp101  PIcp29  PIcp37  Picp68  PIcp93
# France FR1 183/186 171/165 173/179 99/108 131/137 173/176 159/159  98/107 126/126 206/215 107/107  115/121

table(test$PIcp37)
#215/215 
#    207 
table(test$PIcp39)
#165/171 171/171 
#      3     204
table(test$PIcp62)
#176/176 
#    207
table(test$PIcp93)
#121/121 124/124 
#    206       1 

#let's back up and ask about heterozygosity
#vectors follow PIcp21  PIcp39  PIcp45 PIcp71  PIcp58  PIcp62  PIcp75 PIcp101  PIcp29  PIcp37  Picp68 PIcp93
hetf<-c(1, 5, 2, 2, 1, 1, 1, 2, 0, 2, 0, 1)/209
nall<-c(3, 2, 4, 4, 3, 2, 2, 3, 2, 2, 2, 3)


banana <- read.STRand("C:/Users/andre/Documents/Documents/Documents/Icerya/Icerya_12gen_countrypop.txt")
93
sumfull <- fread("C:/Users/andre/Documents/Documents/Documents/Icerya/Icerya_12gen_countrypop.txt")

Ploidies(banana)<-2
#now we need to know repeat lengths
Usatnts(banana) <- c(3, 3, 3, 3, 3, 3, 3, 3, 3,3,3,3)

#simfreq <- deSilvaFreq(banana, self = 0.9, initNull = 0.01,samples = Samples(banana))
simfreq<-simpleFreq(banana, samples = Samples(banana))


testmat4 <- meandistance.matrix2(banana, Samples(banana), freq=simfreq, self=0.9)
clones <- assignClones(testmat4,  threshold=0.01)
#samples=paste("A", 1:100, sep=""),
clones



write.csv(simfreq, "C:/Users/andre/Documents/Documents/Documents/Icerya/icerya_simple_freqs.csv",row.names=T)

#the UK throws errors, perhaps because there's only one sample?


simFst <- calcPopDiff(simfreq, metric = "Fst")
simFst

#              California     France      Turkey        Korea     Mexico South Africa  Spain
#California   0.000000000 
#France       0.032799452 0.00000000 
#Turkey       0.033961091 0.08652313 0.000000000 
#Korea        0.067408654 0.09203047 0.004856395 0.0000000000 
#Mexico       0.022489905 0.03734174 0.067064063 0.0627305683 0.00000000 
#South Africa 0.059530896 0.11184852 0.004731410 0.0001156464 0.08464355 0.0000000000  
#Spain        0.001981244 0.04690972 0.057712795 0.0935279933 0.03546773 0.0908161192 0.0000000000

#let's look at isolation by distance
fstvec<-c(0.032799452,0.033961091, 0.08652313,0.067408654, 0.09203047, 0.004856395, 0.022489905, 0.03734174, 0.067064063, 0.0627305683,
0.059530896, 0.11184852, 0.004731410, 0.0001156464, 0.08464355,0.001981244, 0.04690972, 0.057712795, 0.0935279933, 0.03546773, 0.0908161192)
#

#this time with 2/3rds
#            California   UK         France     Turkey     Korea      Mexico    SouthAfrica Spain      Australia
#California  0.00000000 0.09471236 0.29228521 0.11676164 0.34950011 0.08056277  0.26881482  0.0247270  0.2879996    
#UK          0.00000000 0.00000000 0.04785312 0.26680948 0.47321633 0.31728945  0.43869990  0.3267566  0.2987322  
#France      0.00000000 0.00000000 0.00000000 0.34724067 0.48075689 0.13877442  0.43332900  0.4183241  0.1887788    
#Turkey      0.00000000 0.00000000 0.00000000 0.00000000 0.12824765 0.31998031  0.09443106  0.4354201  0.5016451      
#Korea       0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.86725664  0.01928375  0.8974403  0.6747477      
#Mexico      0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000  0.77468938  0.4840802  0.4718577      
#SouthAfrica 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000  0.00000000  0.8330380  0.6451677    
#Spain       0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000  0.00000000  0.0000000  0.6066772 
#Australia



fstvec23<-c(0.09471236, 0.29228521, 0.11676164, 0.34950011, 0.08056277,  0.26881482,  0.0247270,  0.2879996, 
            0.04785312, 0.26680948, 0.47321633, 0.31728945,  0.43869990,  0.3267566,  0.2987322,
            0.34724067, 0.48075689, 0.13877442,  0.43332900,  0.4183241,  0.1887788, 0.12824765,
            0.31998031,  0.09443106,  0.4354201,  0.5016451, 0.86725664,  0.01928375,  0.8974403,  0.6747477,
            0.77468938,  0.4840802,  0.4718577,0.8330380,  0.6451677,0.6066772)
#there's a better way to do this for sure but I'm hardcoding to push over the finish line
fstmat23<-matrix(nrow=9,ncol=9)
rownames(fstmat23)<-c("Cali","UK","FR","TU","KR","MX","SA","SP","AU")
colnames(fstmat23)<-c("Cali","UK","FR","TU","KR","MX","SA","SP","AU")
fstmat23[1,]<-c(0,0.09471236, 0.29228521, 0.11676164, 0.34950011, 0.08056277,  0.26881482,  0.0247270,  0.2879996)
fstmat23[2,]<-c(0.09471236, 0,0.04785312, 0.26680948, 0.47321633, 0.31728945,  0.43869990,  0.3267566,  0.2987322)
fstmat23[3,]<-c(0.29228521,0.04785312,0,0.34724067, 0.48075689, 0.13877442,  0.43332900,  0.4183241,  0.1887788)
fstmat23[,1]<-c(0,0.09471236, 0.29228521, 0.11676164, 0.34950011, 0.08056277,  0.26881482,  0.0247270,  0.2879996)
fstmat23[,2]<-c(0.09471236, 0,0.04785312, 0.26680948, 0.47321633, 0.31728945,  0.43869990,  0.3267566,  0.2987322)
fstmat23[,3]<-c(0.29228521,0.04785312,0,0.34724067, 0.48075689, 0.13877442,  0.43332900,  0.4183241,  0.1887788)
fstmat23[4,4:9]<-c(0,0.12824765, 0.31998031,  0.09443106,  0.4354201,  0.5016451)
fstmat23[4:9,4]<-c(0,0.12824765, 0.31998031,  0.09443106,  0.4354201,  0.5016451)
fstmat23[5,5:9]<-c(0, 0.86725664,  0.01928375,  0.8974403,  0.6747477)
fstmat23[5:9,5]<-c(0, 0.86725664,  0.01928375,  0.8974403,  0.6747477)
fstmat23[6,6:9]<-c(0,0.77468938,  0.4840802,  0.4718577)
fstmat23[6:9,6]<-c(0,0.77468938,  0.4840802,  0.4718577)
fstmat23[7,7:9]<-c(0,0.8330380,  0.6451677)
fstmat23[7:9,7]<-c(0,0.8330380,  0.6451677)
fstmat23[8,8:9]<-c(0, 0.6066772)
fstmat23[9,8:9]<-c( 0.6066772,0)
#approximate coords: lat, lon 
#Cali: 33.67, 117.19 and 37.77, 122.42 and 33.98, 117.38 mean = 35.14, 119
#UK: 51.43, -0.10
#France: 43.17, -5.61 and 43.61, -3.88    mean = 43.39, -4.75
#Turkey:: 38.89, -29.4 and 40.86, -38.2   mean = 39.88 , -33.8
#Korea: 33.50, -126.53
#Mexico: 20.03, 100.72
#South Africa: -33.23, -21.9
#Spain: 39.47, 0.37
#Australia: -35.30, -149.11 and -34.84, -149.96 mean = -35.07, -149.54

#via google searching webtools
#specifically http://www.meridianoutpost.com/resources/etools/calculators/calculator-latitude-longitude-distance.php?
distvec<-c(9591.09, 9591.09, 2432.79, 9782.61, 9624.50, 7908.97, 2455.14, 9668.09, 11854.24, 12235.86,
16434.55, 8695.78, 8220.72, 13186.71, 14197.39,9602.85, 609.81, 2906.61, 10234.25, 9449.40, 8402.52)

#            California   UK       France     Turkey     Korea      Mexico      SouthAfrica Spain      Australia
#California  0.00000000 8711       9591       11219      9783       2458        16435       9603       12266    
#UK          0.00000000 0.00000000 960        2883       9223       9000        9652        1330       16973  
#France      0.00000000 0.00000000 0.00000000 2433       9625       9670        8696        610        16866    
#Turkey      0.00000000 0.00000000 0.00000000 0.00000000 7909       11857       8221        2907       14438      
#Korea       0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 12238       13187       10234      7985      
#Mexico      0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000  14197       9449       13026      
#SouthAfrica 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000  0.00000000  8403       10672    
#Spain       0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000  0.00000000  0.0000000  17317  
#AUstralia

distvec2<-c(8711,9591,11219,9783,2458,16435,9603,12266,960,2883,9223,9000,9652,1330,16973,2433,9625,9670,8696,610,16866,
           7909,11857,8221, 2907, 14438, 12238, 13187,10234, 7985,14197,9449, 13026,8403,10672, 17317)

#likewise for the distance matrix
distmat23<-matrix(nrow=9,ncol=9)
rownames(distmat23)<-c("Cali","UK","FR","TU","KR","MX","SA","SP","AU")
colnames(distmat23)<-c("Cali","UK","FR","TU","KR","MX","SA","SP","AU")
distmat23[1,]<-c(0,8711,9591,11219,9783,2458,16435,9603,12266)
distmat23[2,]<-c(8711,0,960,2883,9223,9000,9652,1330,16973)
distmat23[3,]<-c(9591,960,0,2433,9625,9670,8696,610,16866)
distmat23[4,]<-c(11219,2883,2433,0,7909,11857,8221, 2907, 14438)
distmat23[5,]<-c(9783,9223, 9625,7909,0,12238, 13187,10234, 7985)
distmat23[6,]<-c(2458,9000,9670,11857,12238,0,14197,9449, 13026)
distmat23[7,]<-c(16435,9652,8696,8221,13187,14197,0,8403,10672)
distmat23[8,9]<-17317
distmat23[8,8]<-0
distmat23[9,9]<-0
distmat23[8,]<-distmat23[,8]
distmat23[9,]<-distmat23[,9]
distmat23[9,8]<-17317

plot(distvec2,fstvec23,xlab="Distance (km)",ylab="Fst", pch=16)
#is there a relationship?
ibd<-lm(fstvec23~distvec2)
summary(ibd)
#Multiple R-squared:  0.003151,  Adjusted R-squared:  -0.04932 
#F-statistic: 0.06005 on 1 and 19 DF,  p-value: 0.809

#biggest diff is france vs south africa

fstd<-dist(fstmat23)
disd<-dist(distmat23)

mantel.rtest(fstd,disd,nrepet=10000)
#Based on 9999 replicates
#Simulated p-value: 0.3069


library("adegenet")
library("data.table")
hem23<-read.csv("C:/Users/andre/Documents/Documents/Documents/Icerya/hem23_data.csv",header=T)
hem_names<-read.csv("C:/Users/andre/Documents/Documents/Documents/Icerya/hem23_names.csv",header=T)

hemfullest<-read.csv("C:/Users/andre/Documents/Documents/Documents/Icerya/all_hems_retry_review_countrypop_data.csv",header=T)
hemfnames<-read.csv("C:/Users/andre/Documents/Documents/Documents/Icerya/all_hems_retry_review_countrypop_names.csv",header=T)
hem_gen<-df2genind(hemfullest,sep="/",pop=hemfnames[,1])
popNames(hem_gen)

cal <- seppop(hem_gen)$California
kor <- seppop(hem_gen)$Korea
aus <- seppop(hem_gen)$Australia
fra <- seppop(hem_gen)$France
mex <- seppop(hem_gen)$Mexico
sou <- seppop(hem_gen)$SouthAfrica
spa <- seppop(hem_gen)$Spain
tur <- seppop(hem_gen)$Turkey
uni <- seppop(hem_gen)$UK



korf <- inbreeding(kor,N=100)
korF <- sapply(korf, mean)
calf <- inbreeding(cal,N=1000)
ausf <- inbreeding(aus,N=1000)
fraf <- inbreeding(fra,N=1000)
mexf <- inbreeding(mex,N=1000)
souf <- inbreeding(sou,N=1000)
spaf <- inbreeding(spa,N=1000)
turf <- inbreeding(tur,N=1000)
unif <- inbreeding(uni,N=1000)
calF <- sapply(calf, mean)

ausF <- sapply(ausf, mean)
fraF <- sapply(fraf, mean)
mexF <- sapply(mexf, mean)
souF <- sapply(souf, mean)
spaF <- sapply(spaf, mean)
turF <- sapply(turf, mean)
uniF <- sapply(unif, mean)

#Korea 56 samples...all with identical genotypes
korF<-c(rep(1,56),rep(-2,20))

#FIGURE 4

#pdf("C:/Users/andre/Documents/Documents/Documents/Icerya/figure_4.pdf",width=8.2, height=6.5)
par(mfrow=c(3,3))
#B,L,T,R
par(mai=c(0.5,0.6,0.5,0.1))
hist(calF,xlim=c(0,1),breaks=seq(from=0.0,to=1,by=.05),main="",xlab="",las=1,col="black",cex.lab=1.2,cex.axis=1.2)
text(0.25,65,"California",cex=1.8)
par(mai=c(0.5,0.1,0.5,0.1))
hist(korF,xlim=c(0,1.07),breaks=100,main="",xlab="",las=1,ylab="",col="black",cex.lab=1.2,cex.axis=1.2)
text(0.2,40,"Korea",cex=1.8)
mtext(side=3, adj=0.5, "Global rates of selfing", cex=1.9, line=1.4)
hist(ausF,xlim=c(0,1),breaks=seq(from=0.0,to=1,by=.05),main="",xlab="",las=1,ylab="",col="black",cex.lab=1.2,cex.axis=1.2)
text(0.2,2.2,"Australia",cex=1.8)
par(mai=c(0.5,0.6,0.05,0.1))
hist(fraF,xlim=c(0,1),breaks=seq(from=0.0,to=1,by=.05),main="",xlab="",las=1,col="black",cex.lab=1.2,cex.axis=1.2)
text(0.2,8,"France",cex=1.8)
par(mai=c(0.5,0.1,0.05,0.1))
hist(mexF,xlim=c(0,1),breaks=seq(from=0.0,to=1,by=.05),main="",xlab="",las=1,ylab="",col="black",cex.lab=1.2,cex.axis=1.2)
text(0.2,2.4,"Mexico",cex=1.8)
hist(souF,xlim=c(0,1),breaks=seq(from=0.0,to=1,by=.05),main="",xlab="",las=1,ylab="",col="black",cex.lab=1.2,cex.axis=1.2)
text(0.27,28.5,"South Africa",cex=1.8)
par(mai=c(0.75,0.6,0.05,0.1))
hist(spaF,xlim=c(0,1),breaks=seq(from=0.0,to=1,by=.05),main="",xlab=expression("F"[IS]),las=1,col="black",cex.lab=1.2,cex.axis=1.2)
text(0.2,24,"Spain",cex=1.8)
par(mai=c(0.75,0.1,0.05,0.1))
hist(turF,xlim=c(0,1),breaks=seq(from=0,to=1,by=.05),main="",xlab=expression("F"[IS]),las=1,ylab="",col="black",cex.lab=1.2,cex.axis=1.2)
text(0.2,22,"Turkey",cex=1.8)
hist(uniF,xlim=c(0,1),breaks=seq(from=0.0,to=1,by=.05),main="",xlab=expression("F"[IS]),las=1,ylab="",col="black",cex.lab=1.2,cex.axis=1.2)
text(0.35,5,"United Kingdom",cex=1.8)
#dev.off()

###testing structure plotting
stru<-fread("C:/Users/andre/Documents/Documents/Documents/Icerya/struct_k3_output.txt",header=T)

hem_names<-read.csv("C:/Users/andre/Documents/Documents/Documents/Icerya/hemfull_names_for_struct.csv",header=T)
hem_names$Label <- paste(hem_names$POP,hem_names$Label,sep="")
head(hem_names)

stru2<-merge(stru,hem_names,by="Label")

#123,132,213,312,321
stru2<-stru2[order(stru2$Q2,stru2$Q3,stru2$Q1),]

#supplementary structure plot
stru3<-as.matrix(stru2[,5:7])
barplot(t(stru3),col=c("violetred3","dodgerblue1","goldenrod1"),names=stru2$Population,las=2)

















###################fresh parent-offspring analysis
#what's our most common genotype?
pareg<-fread("C:/Users/andre/Documents/Documents/Documents/Icerya/icerya_par_egg_fresh.csv",header=T)
#we have 238 families...

paregd<-pareg[which(pareg$Diffs>0),]

length(unique(paregd$Family))
#44 families have at least one parent offspring diff

famX <- split(pareg, pareg$Family)
str(famX)









barplot(c(stru2$Q1,stru2$Q2,stru2$Q3))


#
###troubleshooting the reading of missing data: 

banana <- read.STRand("C:/Users/andre/Documents/Documents/Documents/Icerya/hems_2_3rds_for_fst.txt")
#bbug<-fread("C:/Users/andre/Documents/Documents/Documents/Icerya/all_hems_retry_review_countrypop.txt")

Ploidies(banana)<-2
#now we need to know repeat lengths
Usatnts(banana) <- c(3, 3, 3, 3, 3, 3, 3, 3, 3,3,3,3)

#simfreq <- deSilvaFreq(banana, self = 0.9, initNull = 0.01,samples = Samples(banana))
simfreq<-simpleFreq(banana, samples = Samples(banana))

simFst <- calcPopDiff(simfreq, metric = "Fst")
simFst


simgen <- read.STRand("C:/Users/andre/Documents/Documents/Documents/Icerya/hem23rds_poppop.txt")
Ploidies(simgen)<-2

write.Structure(simgen, ploidy = 2, file="C:/Users/andre/Documents/Documents/Documents/Icerya/iceryaStruct23pop.txt")
#this has 296 samples


#structure runs for reviewer comments
#k = 2 use run 18...looks to be two common haplotypes vs everyone else

stru2<-fread("C:/Users/andre/Documents/Documents/Documents/Icerya/struct_k2_output.txt",header=T)

hem_names<-read.csv("C:/Users/andre/Documents/Documents/Documents/Icerya/hem23rds_poppop.csv",header=T)
hem_names$Label <- paste(hem_names$Pop,hem_names$Label,sep="")

stru2<-merge(stru2,hem_names,by="Label")

#123,132,213,312,321
stru2<-stru2[order(stru2$Q2,stru2$Q1),]
#how did 1 Pakistan sample get through the filter?
stru2<-stru2[which(stru2$Population!="Pakistan"),]
#supplementary structure plot
stru21<-as.matrix(stru2[,5:6])
par(mai=c(1,0.7,0.7,0.1))
#tiff("fig3topright.tiff",height = )
barplot(t(stru21),col=c("violetred3","dodgerblue1"),names=stru2$Population,las=2, main="",cex.names=0.7, cex.main=3.5, cex.axis=1.5)
mtext("K = 2", side = 3, line = 0.5, adj=0.5,cex=3.7)

#K = 3
stru3<-fread("C:/Users/andre/Documents/Documents/Documents/Icerya/struct_k3_output2.txt",header=T)

#hem_names<-read.csv("C:/Users/andre/Documents/Documents/Documents/Icerya/hem23rds_poppop.csv",header=T)
#hem_names$Label <- paste(hem_names$Pop,hem_names$Label,sep="")

stru3<-merge(stru3,hem_names,by="Label")

#123,132,213,312,321
stru3<-stru3[order(stru3$Q2,stru3$Q3,stru3$Q1),]
#how did 1 Pakistan sample get through the filter?
stru3<-stru3[which(stru3$Population!="Pakistan"),]

#supplementary structure plot
stru31<-as.matrix(stru3[,5:7])
barplot(t(stru31),col=c("violetred3","dodgerblue1","goldenrod1"),names=stru3$Population,las=2, main="",cex.names=0.7, cex.main=3.5, cex.axis=1.5)
mtext("K = 3", side = 3, line = 0.5, adj=0.5,cex=3.7)


#k = 4
#use "C:/Users/andre/Documents/Documents/Documents/Icerya/stuct_k41.txt"
stru4<-fread("C:/Users/andre/Documents/Documents/Documents/Icerya/stuct_k41.txt",header=T)
colnames(stru4)[2]<-"Ind"
stru4<-cbind.data.frame(stru4,hem_names)

#123,132,213,312,321
stru4<-stru4[order(stru4$Q2,stru4$Q3,stru4$Q1),]
#how did 1 Pakistan sample get through the filter?
stru4<-stru4[which(stru4$Population!="Pakistan"),]

#supplementary structure plot
stru41<-as.matrix(stru4[,5:8])
#par(mai=c(2,0.7,0.7,0))
barplot(t(stru41),col=c("goldenrod1","violetred3","dodgerblue1","red2"),names=stru4$Population,las=2, main="",cex.names=0.7, cex.main=3.5, cex.axis=1.5)
mtext("K = 4", side = 3, line = 0.5, adj=0.5,cex=3.7)
#barplot(t(stru41),col=c("goldenrod1","violetred3","dodgerblue1","red2"),names=stru4$Pop,las=2, main="K = 4",cex.names=0.7, xlim=c(120,195))

#admixed indivuals here:

och<-read.table("C:/Users/andre/Documents/Documents/Documents/Icerya/iceryaStruct23popn.txt",header=T)

och[276*2,]
#potentially admixed individuals:
#SFSF01
#CE62CE62-2 
#SFSF36  
#CE65CE65-1
#CHCh2 - this is the UK outbred!
#AUAU1-2


###assessing the probability of seeing low-F individuals by genotype error
hprob<-63/3615
hetput<-read.csv("C:/Users/andre/Documents/Documents/Documents/Icerya/icerya_het_counts.csv",header=T)
hetprob<-pbinom(hetput$Het, size=hetput$Loci, prob=hprob,lower=F)
hist(hetprob, breaks = 50, xlab = "Probability of observed heterozygosity", las=1,
     main="Binomial probability of sample heterozygosities",col=c("goldenrod1","white","white","white","white","white","white","white","white","white","white","white","white","white","white","white","white","white","white","white","white"))

#and doing calcs for RMES file generation
rmesdat<-read.csv("C:/Users/andre/Documents/Documents/Documents/Icerya/icerya_genotypes_for_RMES_count_loc.csv",header=T)

table(rmesdat$POP)

length(unique(rmesdat$POP))

library('ape')
