singlevariant<-function(num,gdsfile,varfile=NULL,varidfile,phenfile,nullfile,stat="Score",outfile,mcount=10){

######
###### read gds file
gds <- seqOpen(gdsfile, allow.duplicate=T)


##### samples
phen1<-fread(phenfile,header=T,data.table=F,sep="\t")
names(phen1)[1]<-"sample.id"
samid0<-phen1$sample.id


######
###### QCed variants
if(!is.null(varfile)){
vardata<-fread(varfile,header=T,sep="\t",data.table=F,select=c(1:3))
vardata$minorcount<-apply(vardata[,c("total.ref","total.alt")],1,min)
vardata2<-subset(vardata,minorcount>=mcount)
varid0<-vardata2$variant.id
}else{
varid0 <- seqGetData(gds, "variant.id")
}

######
###### matching samples
samples <- seqGetData(gds, "sample.id")
missamples<-samples[!samples %in% samid0]
misphen<-data.frame(matrix(NA,nrow=length(missamples),ncol=ncol(phen1)))
colnames(misphen)<-names(phen1)
misphen$sample.id<-missamples
combphen<-rbind(phen1,misphen)
rownames(combphen)<-as.character(combphen$sample.id)
combphen2<-combphen[as.character(samples),]

######
###### # construct a SeqVarData object
seqData <- SeqVarData(gds, sampleData=AnnotatedDataFrame(combphen2))

######
###### filter the gdsfile with minor count >= 10
seqSetFilter(seqData, sample.id=samid0, variant.id=varid0)

if(!is.null(varidfile)){
###### high moderate variants
varids<-fread(cmd = paste0("grep -v '^#' ",varidfile),header=F,data.table=F)
print("varids loaded")

library(tidyr)
varids<-separate(varids,col=V1,into=c("chr","pos","ref","alt"),remove=F)
#######
seqSetFilterPos(seqData,chr=varids$chr,pos=varids$pos,intersect=T,multi.pos=T)
print("filter_applied")
######
###### variantInfo
varinfo1<-variantInfo(seqData)
varinfo2<-merge(varinfo1,varids,by=c("chr","pos","ref","alt"))

###### filters matching thier ref and alternatives
seqSetFilter(seqData, sample.id=samid0,variant.id=varinfo2$variant.id)
}

####### frequency < 2%
afreq0<-SeqVarTools::alleleFrequency(seqData,use.names=T)
afreq0<-ifelse(afreq0>0.5,1-afreq0,afreq0)
lowfreq0<-names(which(afreq0<0.02))

######
###### freq < 2%
print("freq < 2%")
seqSetFilter(seqData, sample.id=samid0,variant.id=lowfreq0)

###### not necessary in TOPMed
# mcount0<-SeqVarTools::minorAlleleCount(seqData,use.names=T)
# minmcount0<-names(which(mcount0>=mcount))

###### mcount >= 10
# seqSetFilter(seqData, sample.id=samid0,variant.id=minmcount0)

######
iterator <- SeqVarBlockIterator(seqData,variantBlock=100)

######
###### load null model
nullmod<-get(load(nullfile))

###### perfrom assocation test
assoc <- assocTestSingle(iterator, nullmod, test=stat,verbose=TRUE)

#####
save(assoc,file=outfile)
#####
seqClose(gds)
}
