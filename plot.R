## plot the equivalent of figure 1 in Charlesworth and Charlesworth (1983)
## 1. run `modelCC83`
## 2. rename the output file into `data.csv`

inFile <- "data.csv"
d <- read.table( inFile, header=T, sep="\t" )

table(d$simu)
range(d$gen)

simu <- 4
plot( d$gen[d$simu==simu], d$meanC[d$simu==simu], type="l" )

png( paste(inFile,".png",sep=""), width=900, height=600 )
par( mar=c(5,5,3,2), font=2, font.axis=2, font.lab=2, cex=1.5, lwd=2 )
plot( unique(d$gen), seq(0,max(d$meanC),length.out=length(unique(d$gen))),
     type="n",
     xlab="Generations", ylab="Mean TE copy number",
     main="Model from Charlesworth and Charlesworth (1983)" )
for( s in unique(d$simu) ){
  points( unique(d$gen[seq(1,max(d$gen),10)]),
         d$meanC[d$simu==s][seq(1,max(d$gen),10)],
         type="l", lwd=0.5 )
}
abline( h=mean(d$meanC[d$gen==max(d$gen)]), lty=2 )
dev.off()



## distribution of the nb of TEs per individual
par( mfrow=c(2,2) )
hist( d$n[d$gen==25], xlab="x individuals", ylab="nb of TEs in exactly x individuals", xlim=c(min(d$n),max(d$n)) )
hist( d$n[d$gen==50], xlab="x individuals", ylab="nb of TEs in exactly x individuals", xlim=c(min(d$n),max(d$n)) )
hist( d$n[d$gen==75], xlab="x individuals", ylab="nb of TEs in exactly x individuals", xlim=c(min(d$n),max(d$n)) )
hist( d$n[d$gen==100], xlab="x individuals", ylab="nb of TEs in exactly x individuals", xlim=c(min(d$n),max(d$n)) )
