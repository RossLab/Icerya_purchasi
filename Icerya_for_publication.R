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



#approximate coords: lat, lon 
#Cali: 33.67, 117.19 and 37.77, 122.42 and 33.98, 117.38 mean = 35.14, 119
#France: 43.17, -5.61 and 43.61, -3.88    mean = 43.39, -4.75
#Turkey:: 38.89, -29.4 and 40.86, -38.2   mean = 39.88 , -33.8
#Korea: 33.50, -126.53
#Mexico: 20.03, 100.72
#South Africa: -33.23, -21.9
#Spain: 39.47, 0.37

#via google searching webtools
#specifically http://www.meridianoutpost.com/resources/etools/calculators/calculator-latitude-longitude-distance.php?
distvec<-c(9591.09, 9591.09, 2432.79, 9782.61, 9624.50, 7908.97, 2455.14, 9668.09, 11854.24, 12235.86,
16434.55, 8695.78, 8220.72, 13186.71, 14197.39,9602.85, 609.81, 2906.61, 10234.25, 9449.40, 8402.52)

plot(distvec,fstvec,xlab="Distance (km)",ylab="Fst", pch=16)
#is there a relationship?
ibd<-lm(fstvec~distvec)
summary(ibd)
#Multiple R-squared:  0.003151,  Adjusted R-squared:  -0.04932 
#F-statistic: 0.06005 on 1 and 19 DF,  p-value: 0.809

#biggest diff is france vs south africa



library("adegenet")
library("data.table")
hem23<-read.csv("C:/Users/andre/Documents/Documents/Documents/Icerya/hem23_data.csv",header=T)
hem_names<-read.csv("C:/Users/andre/Documents/Documents/Documents/Icerya/hem23_names.csv",header=T)

hem_gen<-df2genind(hem23,sep="/",pop=hem_names[,1])
popNames(hem_gen)

cal <- seppop(hem_gen)$California
kor <- seppop(hem_gen)$Korea
aus <- seppop(hem_gen)$Australia
fra <- seppop(hem_gen)$France
mex <- seppop(hem_gen)$Mexico
pak <- seppop(hem_gen)$Pakistan
sou <- seppop(hem_gen)$SouthAfrica
spa <- seppop(hem_gen)$Spain
tur <- seppop(hem_gen)$Turkey
uni <- seppop(hem_gen)$UK


calf <- inbreeding(cal,N=100)
korf <- inbreeding(kor,N=100)
ausf <- inbreeding(aus,N=100)
fraf <- inbreeding(fra,N=100)
mexf <- inbreeding(mex,N=100)
pakf <- inbreeding(pak,N=100)
souf <- inbreeding(sou,N=100)
spaf <- inbreeding(spa,N=100)
turf <- inbreeding(tur,N=100)
unif <- inbreeding(uni,N=100)

calF <- sapply(calf, mean)
korF <- sapply(korf, mean)
ausF <- sapply(ausf, mean)
fraF <- sapply(fraf, mean)
mexF <- sapply(mexf, mean)
pakF <- sapply(pakf, mean)
souF <- sapply(souf, mean)
spaF <- sapply(spaf, mean)
turF <- sapply(turf, mean)
uniF <- sapply(unif, mean)

#two of these are broken, Pakistan which has n=1 and Korea, which has 56 samples...all with identical genotypes
korF<-c(rep(1,56),rep(-2,20))

#FIGURE 3
par(mfrow=c(3,3))
#B,L,T,R
par(mai=c(0.5,0.6,0.5,0.1))
hist(calF,xlim=c(0,1),breaks=10,main="",xlab="",las=1,col="black",cex.lab=1.2,cex.axis=1.2)
text(0.2,45,"California",cex=1.8)
par(mai=c(0.5,0.1,0.5,0.1))
hist(korF,xlim=c(0,1.07),breaks=100,main="",xlab="",las=1,ylab="",col="black",cex.lab=1.2,cex.axis=1.2)
text(0.2,47,"Korea",cex=1.8)
mtext(side=3, adj=0.5, "Global rates of selfing", cex=1.9, line=1.4)
hist(ausF,xlim=c(0,1),breaks=10,main="",xlab="",las=1,ylab="",col="black",cex.lab=1.2,cex.axis=1.2)
text(0.2,1.75,"Australia",cex=1.8)
par(mai=c(0.5,0.6,0.05,0.1))
hist(fraF,xlim=c(0,1),breaks=10,main="",xlab="",las=1,col="black",cex.lab=1.2,cex.axis=1.2)
text(0.2,9.3,"France",cex=1.8)
par(mai=c(0.5,0.1,0.05,0.1))
hist(mexF,xlim=c(0,1),breaks=10,main="",xlab="",las=1,ylab="",col="black",cex.lab=1.2,cex.axis=1.2)
text(0.2,1.75,"Mexico",cex=1.8)
hist(souF,xlim=c(0,1),breaks=10,main="",xlab="",las=1,ylab="",col="black",cex.lab=1.2,cex.axis=1.2)
text(0.22,7,"South Africa",cex=1.8)
par(mai=c(0.5,0.6,0.05,0.1))
hist(spaF,xlim=c(0,1),breaks=10,main="",xlab="F",las=1,col="black",cex.lab=1.2,cex.axis=1.2)
text(0.2,28,"Spain",cex=1.8)
par(mai=c(0.5,0.1,0.05,0.1))
hist(turF,xlim=c(0,1),breaks=10,main="",xlab="F",las=1,ylab="",col="black",cex.lab=1.2,cex.axis=1.2)
text(0.2,2.75,"Turkey",cex=1.8)
hist(uniF,xlim=c(0,1),breaks=10,main="",xlab="F",las=1,ylab="",col="black",cex.lab=1.2,cex.axis=1.2)
text(0.2,2.75,"UK",cex=1.8)


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



barplot(c(stru2$Q1,stru2$Q2,stru2$Q3))

