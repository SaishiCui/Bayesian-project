rm(list = ls())
setwd("C:/Users/cuisa/Desktop/Biostat_PhD_First_year/Biostat_PhD_19_Spring_Quarter/BST-604_Applied_Bayesian_Analysis")

drug<-read.table("drug.txt",header = T,sep='\t',stringsAsFactors=FALSE)[,-1]


drug[drug=="Suppressed"]=NA

drug$Population[is.na(drug$Population)]=10

drug$Deaths=as.numeric(drug$Deaths)
drug$Population=as.numeric(drug$Population)

drug$Race.Code[drug$Race=="Black or African American"]<-1
drug$Race.Code[drug$Race=="White"]<-2

drug$Gender.Code[drug$Gender=="Female"]<-1
drug$Gender.Code[drug$Gender=="Male"]<-2

drug$Year.Code[drug$Year==2013]<-1
drug$Year.Code[drug$Year==2014]<-2
drug$Year.Code[drug$Year==2015]<-3
drug$Year.Code[drug$Year==2016]<-4
drug$Year.Code[drug$Year==2017]<-5
drug$Year.Code[drug$Year==2018]<-6
drug$County.Code=rep(1:58,each=24)
drug$Crude.Rate<-as.numeric(drug$Crude.Rate)
drug$County[which(drug$Crude.Rate==max(drug$Crude.Rate[!is.na(drug$Crude.Rate)]))]


drug$Y_cen<-drug$Deaths
drug$Y_cen[is.na(drug$Deaths)]<-10



results<-read.table("WinBUGS_results.txt",header = T,sep='\t',stringsAsFactors=FALSE)[,-c(1,2)]

drug$adj_lambda_median<-results$median


drug$County[which(drug$adj_lambda_median==max(drug$adj_lambda_median))]
 


BF2013=drug$adj_lambda_median[drug$Race=="Black or African American"& drug$Gender=="Female"&drug$Year==2013 ]
BF2014=drug$adj_lambda_median[drug$Race=="Black or African American"& drug$Gender=="Female"&drug$Year==2014 ]
BF2015=drug$adj_lambda_median[drug$Race=="Black or African American"& drug$Gender=="Female"&drug$Year==2015 ]
BF2016=drug$adj_lambda_median[drug$Race=="Black or African American"& drug$Gender=="Female"&drug$Year==2016 ]
BF2017=drug$adj_lambda_median[drug$Race=="Black or African American"& drug$Gender=="Female"&drug$Year==2017 ]
BF2018=drug$adj_lambda_median[drug$Race=="Black or African American"& drug$Gender=="Female"&drug$Year==2018 ]
BF2018_2013=(BF2018-BF2013)/BF2013*100


BM2013=drug$adj_lambda_median[drug$Race=="Black or African American"& drug$Gender=="Male"&drug$Year==2013 ]
BM2014=drug$adj_lambda_median[drug$Race=="Black or African American"& drug$Gender=="Male"&drug$Year==2014 ]
BM2015=drug$adj_lambda_median[drug$Race=="Black or African American"& drug$Gender=="Male"&drug$Year==2015 ]
BM2016=drug$adj_lambda_median[drug$Race=="Black or African American"& drug$Gender=="Male"&drug$Year==2016 ]
BM2017=drug$adj_lambda_median[drug$Race=="Black or African American"& drug$Gender=="Male"&drug$Year==2017 ]
BM2018=drug$adj_lambda_median[drug$Race=="Black or African American"& drug$Gender=="Male"&drug$Year==2018 ]
BM2018_2013=(BM2018-BM2013)/BM2013*100



WF2013=drug$adj_lambda_median[drug$Race=="White"& drug$Gender=="Female"&drug$Year==2013 ]
WF2014=drug$adj_lambda_median[drug$Race=="White"& drug$Gender=="Female"&drug$Year==2014 ]
WF2015=drug$adj_lambda_median[drug$Race=="White"& drug$Gender=="Female"&drug$Year==2015 ]
WF2016=drug$adj_lambda_median[drug$Race=="White"& drug$Gender=="Female"&drug$Year==2016 ]
WF2017=drug$adj_lambda_median[drug$Race=="White"& drug$Gender=="Female"&drug$Year==2017 ]
WF2018=drug$adj_lambda_median[drug$Race=="White"& drug$Gender=="Female"&drug$Year==2018 ]
WF2018_2013=(WF2018-WF2013)/WF2013*100



WM2013=drug$adj_lambda_median[drug$Race=="White"& drug$Gender=="Male"&drug$Year==2013 ]
WM2014=drug$adj_lambda_median[drug$Race=="White"& drug$Gender=="Male"&drug$Year==2014 ]
WM2015=drug$adj_lambda_median[drug$Race=="White"& drug$Gender=="Male"&drug$Year==2015 ]
WM2016=drug$adj_lambda_median[drug$Race=="White"& drug$Gender=="Male"&drug$Year==2016 ]
WM2017=drug$adj_lambda_median[drug$Race=="White"& drug$Gender=="Male"&drug$Year==2017 ]
WM2018=drug$adj_lambda_median[drug$Race=="White"& drug$Gender=="Male"&drug$Year==2018 ]
WM2018_2013=(WM2018-WM2013)/WM2013*100


library(spdep)      # spatial dependencies functions
library(rgdal)      # manages spatial projections
library(maptools)   # general functions for map/spatial manipulation
library(sp)         # foundational definition of spatial classes in R
library(arm)        # Loads R2WinBugs, Lattice, Matrix, others
library(raster)
library(rgeos)
library(RColorBrewer)

load(file='counties.rdata')
#load(file='states.rdata')

Ns=length(us)
keep=(1:Ns)[us$STATE_FIPS=='06']
CA=us[keep,]
Ns=length(CA)



ncols=10
cols=brewer.pal(ncols,'RdYlBu')[ncols:1]
tcuts=quantile(BF2018_2013,1:(ncols-1)/ncols)
tcolb=array(rep(BF2018_2013,each=ncols-1) > tcuts,
            dim=c(ncols-1,Ns))
tcol =apply(tcolb,2,sum)+1

png('CA_BF_change.png',height=520,width=1000)
par(mar=c(0,0,0,10),cex=1)
plot(CA,col=cols[tcol],border='lightgray',lwd=.5)
legend('right',inset=c(-.15,0),xpd=TRUE,
       legend=c(paste(
         c('Below',round(tcuts[-(ncols-1)],0),'Over'),
         c(' ',rep( ' - ',ncols-2),' '),
         c(round(tcuts,0),round(tcuts[ncols-1],0)),sep='')),
       fill=cols,title='Mortality rate change in percent',bty='n',cex=1.5,
       border='lightgray')
dev.off()



ncols=10
cols=brewer.pal(ncols,'RdYlBu')[ncols:1]
tcuts=quantile(BM2018_2013,1:(ncols-1)/ncols)
tcolb=array(rep(BM2018_2013,each=ncols-1) > tcuts,
            dim=c(ncols-1,Ns))
tcol =apply(tcolb,2,sum)+1

png('CA_BM_change.png',height=520,width=1000)
par(mar=c(0,0,0,10),cex=1)
plot(CA,col=cols[tcol],border='lightgray',lwd=.5)
legend('right',inset=c(-.15,0),xpd=TRUE,
       legend=c(paste(
         c('Below',round(tcuts[-(ncols-1)],0),'Over'),
         c(' ',rep( ' - ',ncols-2),' '),
         c(round(tcuts,0),round(tcuts[ncols-1],0)),sep='')),
       fill=cols,title='Mortality rate change in percent',bty='n',cex=1.5,
       border='lightgray')
dev.off()





ncols=10
cols=brewer.pal(ncols,'RdYlBu')[ncols:1]
tcuts=quantile(WF2018_2013,1:(ncols-1)/ncols)
tcolb=array(rep(WF2018_2013,each=ncols-1) > tcuts,
            dim=c(ncols-1,Ns))
tcol =apply(tcolb,2,sum)+1

png('CA_WF_change.png',height=520,width=1000)
par(mar=c(0,0,0,10),cex=1)
plot(CA,col=cols[tcol],border='lightgray',lwd=.5)
legend('right',inset=c(-.15,0),xpd=TRUE,
       legend=c(paste(
         c('Below',round(tcuts[-(ncols-1)],0),'Over'),
         c(' ',rep( ' - ',ncols-2),' '),
         c(round(tcuts,0),round(tcuts[ncols-1],0)),sep='')),
       fill=cols,title='Mortality rate change in percent',bty='n',cex=1.5,
       border='lightgray')
dev.off()





ncols=10
cols=brewer.pal(ncols,'RdYlBu')[ncols:1]
tcuts=quantile(WM2018_2013,1:(ncols-1)/ncols)
tcolb=array(rep(WM2018_2013,each=ncols-1) > tcuts,
            dim=c(ncols-1,Ns))
tcol =apply(tcolb,2,sum)+1

png('CA_WM_change.png',height=520,width=1000)
par(mar=c(0,0,0,10),cex=1)
plot(CA,col=cols[tcol],border='lightgray',lwd=.5)
legend('right',inset=c(-.15,0),xpd=TRUE,
       legend=c(paste(
         c('Below',round(tcuts[-(ncols-1)],0),'Over'),
         c(' ',rep( ' - ',ncols-2),' '),
         c(round(tcuts,0),round(tcuts[ncols-1],0)),sep='')),
       fill=cols,title='Mortality rate change in percent',bty='n',cex=1.5,
       border='lightgray')
dev.off()
















ncols=7
cols=brewer.pal(ncols,'RdYlBu')[ncols:1]
tcuts=quantile(BF2013,1:(ncols-1)/ncols)
tcolb=array(rep(BF2013,each=ncols-1) > tcuts,
            dim=c(ncols-1,Ns))
tcol =apply(tcolb,2,sum)+1

png('CA_BF_2013.png',height=520,width=1000)
par(mar=c(0,0,0,10),cex=1)
plot(CA,col=cols[tcol],border='lightgray',lwd=.5)
legend('right',inset=c(-.15,0),xpd=TRUE,
       legend=c(paste(
         c('Below',round(tcuts[-(ncols-1)],0),'Over'),
         c(' ',rep( ' - ',ncols-2),' '),
         c(round(tcuts,0),round(tcuts[ncols-1],0)),sep='')),
       fill=cols,title='Deaths per 100,000',bty='n',cex=1.5,
       border='lightgray')
dev.off()




ncols=7
cols=brewer.pal(ncols,'RdYlBu')[ncols:1]
tcuts=quantile(BF2014,1:(ncols-1)/ncols)
tcolb=array(rep(BF2014,each=ncols-1) > tcuts,
            dim=c(ncols-1,Ns))
tcol =apply(tcolb,2,sum)+1

png('CA_BF_2014.png',height=520,width=1000)
par(mar=c(0,0,0,10),cex=1)
plot(CA,col=cols[tcol],border='lightgray',lwd=.5)
legend('right',inset=c(-.15,0),xpd=TRUE,
       legend=c(paste(
         c('Below',round(tcuts[-(ncols-1)],0),'Over'),
         c(' ',rep( ' - ',ncols-2),' '),
         c(round(tcuts,0),round(tcuts[ncols-1],0)),sep='')),
       fill=cols,title='Deaths per 100,000',bty='n',cex=1.5,
       border='lightgray')
dev.off()



ncols=7
cols=brewer.pal(ncols,'RdYlBu')[ncols:1]
tcuts=quantile(BF2015,1:(ncols-1)/ncols)
tcolb=array(rep(BF2015,each=ncols-1) > tcuts,
            dim=c(ncols-1,Ns))
tcol =apply(tcolb,2,sum)+1

png('CA_BF_2015.png',height=520,width=1000)
par(mar=c(0,0,0,10),cex=1)
plot(CA,col=cols[tcol],border='lightgray',lwd=.5)
legend('right',inset=c(-.15,0),xpd=TRUE,
       legend=c(paste(
         c('Below',round(tcuts[-(ncols-1)],0),'Over'),
         c(' ',rep( ' - ',ncols-2),' '),
         c(round(tcuts,0),round(tcuts[ncols-1],0)),sep='')),
       fill=cols,title='Deaths per 100,000',bty='n',cex=1.5,
       border='lightgray')
dev.off()


ncols=7
cols=brewer.pal(ncols,'RdYlBu')[ncols:1]
tcuts=quantile(BF2016,1:(ncols-1)/ncols)
tcolb=array(rep(BF2016,each=ncols-1) > tcuts,
            dim=c(ncols-1,Ns))
tcol =apply(tcolb,2,sum)+1

png('CA_BF_2016.png',height=520,width=1000)
par(mar=c(0,0,0,10),cex=1)
plot(CA,col=cols[tcol],border='lightgray',lwd=.5)
legend('right',inset=c(-.15,0),xpd=TRUE,
       legend=c(paste(
         c('Below',round(tcuts[-(ncols-1)],0),'Over'),
         c(' ',rep( ' - ',ncols-2),' '),
         c(round(tcuts,0),round(tcuts[ncols-1],0)),sep='')),
       fill=cols,title='Deaths per 100,000',bty='n',cex=1.5,
       border='lightgray')
dev.off()



ncols=7
cols=brewer.pal(ncols,'RdYlBu')[ncols:1]
tcuts=quantile(BF2017,1:(ncols-1)/ncols)
tcolb=array(rep(BF2017,each=ncols-1) > tcuts,
            dim=c(ncols-1,Ns))
tcol =apply(tcolb,2,sum)+1

png('CA_BF_2017.png',height=520,width=1000)
par(mar=c(0,0,0,10),cex=1)
plot(CA,col=cols[tcol],border='lightgray',lwd=.5)
legend('right',inset=c(-.15,0),xpd=TRUE,
       legend=c(paste(
         c('Below',round(tcuts[-(ncols-1)],0),'Over'),
         c(' ',rep( ' - ',ncols-2),' '),
         c(round(tcuts,0),round(tcuts[ncols-1],0)),sep='')),
       fill=cols,title='Deaths per 100,000',bty='n',cex=1.5,
       border='lightgray')
dev.off()



ncols=7
cols=brewer.pal(ncols,'RdYlBu')[ncols:1]
tcuts=quantile(BF2018,1:(ncols-1)/ncols)
tcolb=array(rep(BF2018,each=ncols-1) > tcuts,
            dim=c(ncols-1,Ns))
tcol =apply(tcolb,2,sum)+1

png('CA_BF_2018.png',height=520,width=1000)
par(mar=c(0,0,0,10),cex=1)
plot(CA,col=cols[tcol],border='lightgray',lwd=.5)
legend('right',inset=c(-.15,0),xpd=TRUE,
       legend=c(paste(
         c('Below',round(tcuts[-(ncols-1)],0),'Over'),
         c(' ',rep( ' - ',ncols-2),' '),
         c(round(tcuts,0),round(tcuts[ncols-1],0)),sep='')),
       fill=cols,title='Deaths per 100,000',bty='n',cex=1.5,
       border='lightgray')
dev.off()




ncols=7
cols=brewer.pal(ncols,'RdYlBu')[ncols:1]
tcuts=quantile(BM2013,1:(ncols-1)/ncols)
tcolb=array(rep(BM2013,each=ncols-1) > tcuts,
            dim=c(ncols-1,Ns))
tcol =apply(tcolb,2,sum)+1

png('CA_BM_2013.png',height=520,width=1000)
par(mar=c(0,0,0,10),cex=1)
plot(CA,col=cols[tcol],border='lightgray',lwd=.5)
legend('right',inset=c(-.15,0),xpd=TRUE,
       legend=c(paste(
         c('Below',round(tcuts[-(ncols-1)],0),'Over'),
         c(' ',rep( ' - ',ncols-2),' '),
         c(round(tcuts,0),round(tcuts[ncols-1],0)),sep='')),
       fill=cols,title='Deaths per 100,000',bty='n',cex=1.5,
       border='lightgray')
dev.off()



ncols=7
cols=brewer.pal(ncols,'RdYlBu')[ncols:1]
tcuts=quantile(BM2014,1:(ncols-1)/ncols)
tcolb=array(rep(BM2014,each=ncols-1) > tcuts,
            dim=c(ncols-1,Ns))
tcol =apply(tcolb,2,sum)+1

png('CA_BM_2014.png',height=520,width=1000)
par(mar=c(0,0,0,10),cex=1)
plot(CA,col=cols[tcol],border='lightgray',lwd=.5)
legend('right',inset=c(-.15,0),xpd=TRUE,
       legend=c(paste(
         c('Below',round(tcuts[-(ncols-1)],0),'Over'),
         c(' ',rep( ' - ',ncols-2),' '),
         c(round(tcuts,0),round(tcuts[ncols-1],0)),sep='')),
       fill=cols,title='Deaths per 100,000',bty='n',cex=1.5,
       border='lightgray')
dev.off()




ncols=7
cols=brewer.pal(ncols,'RdYlBu')[ncols:1]
tcuts=quantile(BM2015,1:(ncols-1)/ncols)
tcolb=array(rep(BM2015,each=ncols-1) > tcuts,
            dim=c(ncols-1,Ns))
tcol =apply(tcolb,2,sum)+1

png('CA_BM_2015.png',height=520,width=1000)
par(mar=c(0,0,0,10),cex=1)
plot(CA,col=cols[tcol],border='lightgray',lwd=.5)
legend('right',inset=c(-.15,0),xpd=TRUE,
       legend=c(paste(
         c('Below',round(tcuts[-(ncols-1)],0),'Over'),
         c(' ',rep( ' - ',ncols-2),' '),
         c(round(tcuts,0),round(tcuts[ncols-1],0)),sep='')),
       fill=cols,title='Deaths per 100,000',bty='n',cex=1.5,
       border='lightgray')
dev.off()




ncols=7
cols=brewer.pal(ncols,'RdYlBu')[ncols:1]
tcuts=quantile(BM2016,1:(ncols-1)/ncols)
tcolb=array(rep(BM2016,each=ncols-1) > tcuts,
            dim=c(ncols-1,Ns))
tcol =apply(tcolb,2,sum)+1

png('CA_BM_2016.png',height=520,width=1000)
par(mar=c(0,0,0,10),cex=1)
plot(CA,col=cols[tcol],border='lightgray',lwd=.5)
legend('right',inset=c(-.15,0),xpd=TRUE,
       legend=c(paste(
         c('Below',round(tcuts[-(ncols-1)],0),'Over'),
         c(' ',rep( ' - ',ncols-2),' '),
         c(round(tcuts,0),round(tcuts[ncols-1],0)),sep='')),
       fill=cols,title='Deaths per 100,000',bty='n',cex=1.5,
       border='lightgray')
dev.off()


ncols=7
cols=brewer.pal(ncols,'RdYlBu')[ncols:1]
tcuts=quantile(BM2017,1:(ncols-1)/ncols)
tcolb=array(rep(BM2017,each=ncols-1) > tcuts,
            dim=c(ncols-1,Ns))
tcol =apply(tcolb,2,sum)+1

png('CA_BM_2017.png',height=520,width=1000)
par(mar=c(0,0,0,10),cex=1)
plot(CA,col=cols[tcol],border='lightgray',lwd=.5)
legend('right',inset=c(-.15,0),xpd=TRUE,
       legend=c(paste(
         c('Below',round(tcuts[-(ncols-1)],0),'Over'),
         c(' ',rep( ' - ',ncols-2),' '),
         c(round(tcuts,0),round(tcuts[ncols-1],0)),sep='')),
       fill=cols,title='Deaths per 100,000',bty='n',cex=1.5,
       border='lightgray')
dev.off()


ncols=7
cols=brewer.pal(ncols,'RdYlBu')[ncols:1]
tcuts=quantile(BM2018,1:(ncols-1)/ncols)
tcolb=array(rep(BM2018,each=ncols-1) > tcuts,
            dim=c(ncols-1,Ns))
tcol =apply(tcolb,2,sum)+1

png('CA_BM_2018.png',height=520,width=1000)
par(mar=c(0,0,0,10),cex=1)
plot(CA,col=cols[tcol],border='lightgray',lwd=.5)
legend('right',inset=c(-.15,0),xpd=TRUE,
       legend=c(paste(
         c('Below',round(tcuts[-(ncols-1)],0),'Over'),
         c(' ',rep( ' - ',ncols-2),' '),
         c(round(tcuts,0),round(tcuts[ncols-1],0)),sep='')),
       fill=cols,title='Deaths per 100,000',bty='n',cex=1.5,
       border='lightgray')
dev.off()



ncols=7
cols=brewer.pal(ncols,'RdYlBu')[ncols:1]
tcuts=quantile(WF2013,1:(ncols-1)/ncols)
tcolb=array(rep(WF2013,each=ncols-1) > tcuts,
            dim=c(ncols-1,Ns))
tcol =apply(tcolb,2,sum)+1

png('CA_WF_2013.png',height=520,width=1000)
par(mar=c(0,0,0,10),cex=1)
plot(CA,col=cols[tcol],border='lightgray',lwd=.5)
legend('right',inset=c(-.15,0),xpd=TRUE,
       legend=c(paste(
         c('Below',round(tcuts[-(ncols-1)],0),'Over'),
         c(' ',rep( ' - ',ncols-2),' '),
         c(round(tcuts,0),round(tcuts[ncols-1],0)),sep='')),
       fill=cols,title='Deaths per 100,000',bty='n',cex=1.5,
       border='lightgray')
dev.off()



ncols=7
cols=brewer.pal(ncols,'RdYlBu')[ncols:1]
tcuts=quantile(WF2014,1:(ncols-1)/ncols)
tcolb=array(rep(WF2014,each=ncols-1) > tcuts,
            dim=c(ncols-1,Ns))
tcol =apply(tcolb,2,sum)+1

png('CA_WF_2014.png',height=520,width=1000)
par(mar=c(0,0,0,10),cex=1)
plot(CA,col=cols[tcol],border='lightgray',lwd=.5)
legend('right',inset=c(-.15,0),xpd=TRUE,
       legend=c(paste(
         c('Below',round(tcuts[-(ncols-1)],0),'Over'),
         c(' ',rep( ' - ',ncols-2),' '),
         c(round(tcuts,0),round(tcuts[ncols-1],0)),sep='')),
       fill=cols,title='Deaths per 100,000',bty='n',cex=1.5,
       border='lightgray')
dev.off()



ncols=7
cols=brewer.pal(ncols,'RdYlBu')[ncols:1]
tcuts=quantile(WF2015,1:(ncols-1)/ncols)
tcolb=array(rep(WF2015,each=ncols-1) > tcuts,
            dim=c(ncols-1,Ns))
tcol =apply(tcolb,2,sum)+1

png('CA_WF_2015.png',height=520,width=1000)
par(mar=c(0,0,0,10),cex=1)
plot(CA,col=cols[tcol],border='lightgray',lwd=.5)
legend('right',inset=c(-.15,0),xpd=TRUE,
       legend=c(paste(
         c('Below',round(tcuts[-(ncols-1)],0),'Over'),
         c(' ',rep( ' - ',ncols-2),' '),
         c(round(tcuts,0),round(tcuts[ncols-1],0)),sep='')),
       fill=cols,title='Deaths per 100,000',bty='n',cex=1.5,
       border='lightgray')
dev.off()





ncols=7
cols=brewer.pal(ncols,'RdYlBu')[ncols:1]
tcuts=quantile(WF2016,1:(ncols-1)/ncols)
tcolb=array(rep(WF2016,each=ncols-1) > tcuts,
            dim=c(ncols-1,Ns))
tcol =apply(tcolb,2,sum)+1

png('CA_WF_2016.png',height=520,width=1000)
par(mar=c(0,0,0,10),cex=1)
plot(CA,col=cols[tcol],border='lightgray',lwd=.5)
legend('right',inset=c(-.15,0),xpd=TRUE,
       legend=c(paste(
         c('Below',round(tcuts[-(ncols-1)],0),'Over'),
         c(' ',rep( ' - ',ncols-2),' '),
         c(round(tcuts,0),round(tcuts[ncols-1],0)),sep='')),
       fill=cols,title='Deaths per 100,000',bty='n',cex=1.5,
       border='lightgray')
dev.off()



ncols=7
cols=brewer.pal(ncols,'RdYlBu')[ncols:1]
tcuts=quantile(WF2017,1:(ncols-1)/ncols)
tcolb=array(rep(WF2017,each=ncols-1) > tcuts,
            dim=c(ncols-1,Ns))
tcol =apply(tcolb,2,sum)+1

png('CA_WF_2017.png',height=520,width=1000)
par(mar=c(0,0,0,10),cex=1)
plot(CA,col=cols[tcol],border='lightgray',lwd=.5)
legend('right',inset=c(-.15,0),xpd=TRUE,
       legend=c(paste(
         c('Below',round(tcuts[-(ncols-1)],0),'Over'),
         c(' ',rep( ' - ',ncols-2),' '),
         c(round(tcuts,0),round(tcuts[ncols-1],0)),sep='')),
       fill=cols,title='Deaths per 100,000',bty='n',cex=1.5,
       border='lightgray')
dev.off()



ncols=7
cols=brewer.pal(ncols,'RdYlBu')[ncols:1]
tcuts=quantile(WF2018,1:(ncols-1)/ncols)
tcolb=array(rep(WF2018,each=ncols-1) > tcuts,
            dim=c(ncols-1,Ns))
tcol =apply(tcolb,2,sum)+1

png('CA_WF_2018.png',height=520,width=1000)
par(mar=c(0,0,0,10),cex=1)
plot(CA,col=cols[tcol],border='lightgray',lwd=.5)
legend('right',inset=c(-.15,0),xpd=TRUE,
       legend=c(paste(
         c('Below',round(tcuts[-(ncols-1)],0),'Over'),
         c(' ',rep( ' - ',ncols-2),' '),
         c(round(tcuts,0),round(tcuts[ncols-1],0)),sep='')),
       fill=cols,title='Deaths per 100,000',bty='n',cex=1.5,
       border='lightgray')
dev.off()


ncols=7
cols=brewer.pal(ncols,'RdYlBu')[ncols:1]
tcuts=quantile(WM2013,1:(ncols-1)/ncols)
tcolb=array(rep(WM2013,each=ncols-1) > tcuts,
            dim=c(ncols-1,Ns))
tcol =apply(tcolb,2,sum)+1

png('CA_WM_2013.png',height=520,width=1000)
par(mar=c(0,0,0,10),cex=1)
plot(CA,col=cols[tcol],border='lightgray',lwd=.5)
legend('right',inset=c(-.15,0),xpd=TRUE,
       legend=c(paste(
         c('Below',round(tcuts[-(ncols-1)],0),'Over'),
         c(' ',rep( ' - ',ncols-2),' '),
         c(round(tcuts,0),round(tcuts[ncols-1],0)),sep='')),
       fill=cols,title='Deaths per 100,000',bty='n',cex=1.5,
       border='lightgray')
dev.off()




ncols=7
cols=brewer.pal(ncols,'RdYlBu')[ncols:1]
tcuts=quantile(WM2014,1:(ncols-1)/ncols)
tcolb=array(rep(WM2014,each=ncols-1) > tcuts,
            dim=c(ncols-1,Ns))
tcol =apply(tcolb,2,sum)+1

png('CA_WM_2014.png',height=520,width=1000)
par(mar=c(0,0,0,10),cex=1)
plot(CA,col=cols[tcol],border='lightgray',lwd=.5)
legend('right',inset=c(-.15,0),xpd=TRUE,
       legend=c(paste(
         c('Below',round(tcuts[-(ncols-1)],0),'Over'),
         c(' ',rep( ' - ',ncols-2),' '),
         c(round(tcuts,0),round(tcuts[ncols-1],0)),sep='')),
       fill=cols,title='Deaths per 100,000',bty='n',cex=1.5,
       border='lightgray')
dev.off()



ncols=7
cols=brewer.pal(ncols,'RdYlBu')[ncols:1]
tcuts=quantile(WM2015,1:(ncols-1)/ncols)
tcolb=array(rep(WM2015,each=ncols-1) > tcuts,
            dim=c(ncols-1,Ns))
tcol =apply(tcolb,2,sum)+1

png('CA_WM_2015.png',height=520,width=1000)
par(mar=c(0,0,0,10),cex=1)
plot(CA,col=cols[tcol],border='lightgray',lwd=.5)
legend('right',inset=c(-.15,0),xpd=TRUE,
       legend=c(paste(
         c('Below',round(tcuts[-(ncols-1)],0),'Over'),
         c(' ',rep( ' - ',ncols-2),' '),
         c(round(tcuts,0),round(tcuts[ncols-1],0)),sep='')),
       fill=cols,title='Deaths per 100,000',bty='n',cex=1.5,
       border='lightgray')
dev.off()




ncols=7
cols=brewer.pal(ncols,'RdYlBu')[ncols:1]
tcuts=quantile(WM2016,1:(ncols-1)/ncols)
tcolb=array(rep(WM2016,each=ncols-1) > tcuts,
            dim=c(ncols-1,Ns))
tcol =apply(tcolb,2,sum)+1

png('CA_WM_2016.png',height=520,width=1000)
par(mar=c(0,0,0,10),cex=1)
plot(CA,col=cols[tcol],border='lightgray',lwd=.5)
legend('right',inset=c(-.15,0),xpd=TRUE,
       legend=c(paste(
         c('Below',round(tcuts[-(ncols-1)],0),'Over'),
         c(' ',rep( ' - ',ncols-2),' '),
         c(round(tcuts,0),round(tcuts[ncols-1],0)),sep='')),
       fill=cols,title='Deaths per 100,000',bty='n',cex=1.5,
       border='lightgray')
dev.off()


ncols=7
cols=brewer.pal(ncols,'RdYlBu')[ncols:1]
tcuts=quantile(WM2017,1:(ncols-1)/ncols)
tcolb=array(rep(WM2017,each=ncols-1) > tcuts,
            dim=c(ncols-1,Ns))
tcol =apply(tcolb,2,sum)+1

png('CA_WM_2017.png',height=520,width=1000)
par(mar=c(0,0,0,10),cex=1)
plot(CA,col=cols[tcol],border='lightgray',lwd=.5)
legend('right',inset=c(-.15,0),xpd=TRUE,
       legend=c(paste(
         c('Below',round(tcuts[-(ncols-1)],0),'Over'),
         c(' ',rep( ' - ',ncols-2),' '),
         c(round(tcuts,0),round(tcuts[ncols-1],0)),sep='')),
       fill=cols,title='Deaths per 100,000',bty='n',cex=1.5,
       border='lightgray')
dev.off()


ncols=7
cols=brewer.pal(ncols,'RdYlBu')[ncols:1]
tcuts=quantile(WM2018,1:(ncols-1)/ncols)
tcolb=array(rep(WM2018,each=ncols-1) > tcuts,
            dim=c(ncols-1,Ns))
tcol =apply(tcolb,2,sum)+1

png('CA_WM_2018.png',height=520,width=1000)
par(mar=c(0,0,0,10),cex=1)
plot(CA,col=cols[tcol],border='lightgray',lwd=.5)
legend('right',inset=c(-.15,0),xpd=TRUE,
       legend=c(paste(
         c('Below',round(tcuts[-(ncols-1)],0),'Over'),
         c(' ',rep( ' - ',ncols-2),' '),
         c(round(tcuts,0),round(tcuts[ncols-1],0)),sep='')),
       fill=cols,title='Deaths per 100,000',bty='n',cex=1.5,
       border='lightgray')
dev.off()
















Nr=2        #two race groups
Ng=2        #two gender groups
Ns=length(unique(drug$County.Code))
Ns=58     #58 counties in CA
clabs=unique(drug$County)
Nt=6  # 2013-2018, 6 years 


nsupp<-array(0,dim=c(Nt,Ng,Nr))






#Next we organize things a bit...
Y=array(drug$Deaths,dim=c(Nt,Ng,Nr,Ns),dimnames = list(c("2013","2014","2015","2016","2017","2018"),
                                                       c("F","M"),c("B","W") ))
n=array(drug$Population,dim=c(Nt,Ng,Nr,Ns),dimnames = list(c("2013","2014","2015","2016","2017","2018"),
                                                           c("F","M"),c("B","W") ))

dY<-!is.na(Y)

###################


tau<-array(1000,dim = c(Ng,Nr)) #big prior variances for beta
qbeta=0.9 # Candidate distribution variance 
nits=10000 # Number of interations 
beta0<-array(0,dim = c(Ng,Nr,nits))
beta1<-array(0,dim = c(Ng,Nr,nits))




for (it in 2:nits) {
  
  
    
    
    for (g in 1:Ng) {
      
        for (r in 1:Nr) {
          
             
          for (t in 1:Nt) {
          ## update missing Y ##
          nsupp[t,g,r]=sum(is.na(Y[t,g,r,]))
    
          
          if(nsupp[t,g,r]>0){
            for (i in 1:nsupp[t,g,r]) {
              
              u=runif(1,0,ppois(9.5,n[t,g,r,][!dY[t,g,r,]][i] * exp(beta0[g,r,it-1]+beta1[g,r,it-1]*(t-1) )  ) )
              Y[t,g,r,][!dY[t,g,r,]][i]=qpois(u,n[t,g,r,][!dY[t,g,r,]][i] * exp(beta0[g,r,it-1]+beta1[g,r,it-1]*(t-1) ) )
            }
          }
         } 
      
          # update beta0
          beta0s=rnorm(1,beta0[g,r,it-1],qbeta)
          ra=sum(Y[,g,r,])*(beta0s-beta0[g,r,it-1])
          rb=0
          
          for (t in 1:Nt) {
            for (s in 1:Ns) {
              
            
          
          rb=rb+(-n[t,g,r,s]*exp(beta0s+beta1[g,r,it-1]*(t-1))+n[t,g,r,s]*exp(beta0[g,r,it-1]+beta1[g,r,it-1]*(t-1)))
            }
          
          }
          
         
          rc= -1/(2*tau[g,r])*( beta0s^2-beta0[g,r,it-1]^2)
          
          R=exp(ra+rb+rc)
          beta0[g,r,it]=ifelse(R>runif(1),beta0s,beta0[g,r,it-1])
          
          
          
          
        
          # update beta1
          beta1s=rnorm(1,beta1[g,r,it-1],qbeta)
          ra=0
          rb=0
          
          for (t in 1:Nt) {
            for (s in 1:Ns) {
              ra=ra+Y[t,g,r,s]*(t-1)
              
              rb=rb+(  -n[t,g,r,s]*exp(beta0[g,r,it]+beta1s*(t-1)) + n[t,g,r,s]*exp(beta0[g,r,it]+beta1[g,r,it-1]*(t-1))   )
              
              
            }
          } 
          ra=ra*(beta1s-beta1[g,r,it-1])
          rc= -1/(2*tau[g,r])*(beta1s^2-beta1[g,r,it-1]^2)
          R=exp(ra+rb+rc)
          beta1[g,r,it]=ifelse(R>runif(1),beta1s,beta1[g,r,it-1])
          
          
         
          
        }
      
    }
    
  }


par(mfrow=c(2,4))
plot(beta0[1,1,][-c(1:200)],type = "l",ylab = "beta 0",main = "Black Female")
plot(beta0[1,2,][-c(1:400)],type = "l",ylab = "beta 0",main = "White Female")
plot(beta0[2,1,][-c(1:200)],type = "l",ylab = "beta 0",main = "Black Male")
plot(beta0[2,2,][-c(1:200)],type = "l",ylab = "beta 0",main = "White Male")
plot(beta1[1,1,][-c(1:200)],type = "l",ylab = "beta 1",main = "Black Female")
plot(beta1[1,2,][-c(1:400)],type = "l",ylab = "beta 1",main = "White Female")
plot(beta1[2,1,][-c(1:200)],type = "l",ylab = "beta 1",main = "Black Male")
plot(beta1[2,2,][-c(1:200)],type = "l",ylab = "beta 1",main = "White Male")
dev.off()



par(mfrow=c(2,4))
hist(beta0[1,1,],breaks=100,main = "Black Female",xlab = "beta0")
hist(beta0[1,2,],breaks=100,main = "White Female",xlab = "beta0")
hist(beta0[2,1,],breaks=100,main = "Black Male",xlab = "beta0")
hist(beta0[2,2,],breaks=100,main = "White Male",xlab = "beta0")
hist(beta1[1,1,],breaks=100,main = "Black Female",xlab = "beta1")
hist(beta1[1,2,],breaks=100,main = "White Female",xlab = "beta1")
hist(beta1[2,1,],breaks=100,main = "Black Male",xlab = "beta1")
hist(beta1[2,2,],breaks=100,main = "White Male",xlab = "beta1")
dev.off()

par(mfrow=c(2,4))
acf(beta0[1,1,],main = "beta_0 Black Female")
acf(beta0[1,2,],main = "beta_0 White Female")
acf(beta0[2,1,],main = "beta_0 Black Male")
acf(beta0[2,2,],main = "beta_0 White Male")
acf(beta1[1,1,],main = "beta_1 Black Female")
acf(beta1[1,2,],main = "beta_1 White Female")
acf(beta1[2,1,],main = "beta_1 Black Male")
acf(beta1[2,2,],main = "beta_1 White Male")






quantile(beta0[1,1,],c(0.025,0.5,0.975))
quantile(beta0[1,2,],c(0.025,0.5,0.975))
quantile(beta0[2,1,],c(0.025,0.5,0.975))
quantile(beta0[2,2,],c(0.025,0.5,0.975))
quantile(beta1[1,1,],c(0.025,0.5,0.975))
quantile(beta1[1,2,],c(0.025,0.5,0.975))
quantile(beta1[2,1,],c(0.025,0.5,0.975))
quantile(beta1[2,2,],c(0.025,0.5,0.975))


drug$Y_cen[which(drug$Y_cen==10)]<-9.5


source("bugs_data_function.R")

bugs.data(vari = list(Y=drug$Y_cen),file = 'Bugs_cen.txt')

bugs.data(vari = list(Y=drug$Deaths,
                      n=drug$Population,
                      year=drug$Year.Code,
                      race=drug$Race.Code,
                      sex=drug$Gender.Code,
                      county=drug$County.Code,
                      Y.cen=drug$Y_cen),file = 'BUGS_data.txt')





length(drug$Population[is.na(drug$Deaths)])
