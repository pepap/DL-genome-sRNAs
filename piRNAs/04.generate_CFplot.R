library(GenomicRanges)
library(GenomicAlignments)
library(data.table)

load( "pst01-Robjs.rda",verbose=T )
load( "PST.mrg05kb.frac.dt.rda",verbose=T )

BAMFILE="/storage/brno12-cerit/home/pepap/brno1/TobiasBer/08.DL.new.20241117/UNCOLLAPSED_BAMS/DLovo_mrg3.se.sortCrd.bam"
RLRANGE=c(27,30)
NORMRAN=c(18,32)
ANNOTGR=PST.mrg05kb.gr
OUTFRAC="PST.mrg05kb.frac.dt"

#> PLOT FUN
plotCF <- function( INDT=get(OUTFRAC),PLOT.VAR="RPM",PERC=0.90,XLIM=c(0,100),BY.X=50,BY.YL=0.5,BY.YR=1000,TIT=OUTFRAC,OUT.PDF=NULL,NMDT=norm.DLovo_mrg3.dt ) {

 if ( is.null(OUT.PDF) ) {
 OUT.PDF <- paste0( "outCFplot-",PLOT.VAR,"-",round(PERC*100,0),"thPerc" )
 }

 VS <- INDT[ order( get(PLOT.VAR),decreasing=T ) ][[PLOT.VAR]]
 CS <- cumsum( VS )/sum( VS )
 PC <- (seq(nrow(INDT))[ CS>=PERC ])[1]

 if ( is.null(XLIM) ) {
 XX <- seq( from=0,to=(((nrow(INDT)-(nrow(INDT)%%BY.X))/BY.X)+1)*BY.X,by=BY.X )
 } else {
 XX <- seq( from=0,to=max(XLIM),by=BY.X )
 }
 YL <- seq( from=0,to=1,by=BY.YL )
 YR <- seq( from=0,to=max(VS),by=BY.YR )

 pdf( file=paste0(OUT.PDF,"-",pepDate(),".pdf"),width=10,height=10,pointsize=20 )

 par( bg="white",mar=c(5,5,3,5),pty="s" )
 plot(
  c(0,CS),xlim=XLIM,ylim=c(0,max(CS)),
#!#  c(0,CS),xlim=range(XX),ylim=c(0,max(CS)),
  type="b",pch=16,col="darkorange",cex=1.00,
  main="", xlab="",ylab="",xaxt="n",yaxt="n"
 )
 axis( side=2,at=seq(0,1,0.1),labels=rep.int("",length(seq(0,1,0.1))) )
 axis( side=2,at=YL,labels=YL,col.axis="darkorange" )
 mtext( text="piRNAs (CF)",side=2,col="darkorange",line=3,font=2 )
 arrows( x0=PC,y0=0,x1=PC,y1=CS[PC],code=0,col="red",lty=5,lwd=1 )
 mtext( text=PC,side=1,at=PC,col="red",line=1 )
 text( x=PC,y=(CS[PC])/3,labels=paste0(round(PERC*100,0),"th percentile"),col="red",pos=4,offset=+1,srt=90,font=2 )
 par( new=T )
 plot(
  VS,xlim=XLIM,ylim=c(0,max(VS)),
#!#  VS,xlim=range(XX),ylim=c(0,max(VS)),
  type="p",pch=16,col="gray50",    cex=0.75,
  main=TIT,xlab="",ylab="",xaxt="n",yaxt="n"
 )
 axis( side=4,at=YR,labels=YR,col.axis="gray50" )
 mtext( text=paste0("piRNAs ",PLOT.VAR),side=4,col="gray50",line=3,font=2 )
 axis( side=1,at=XX,labels=XX,col.axis="black" )
 mtext( text=paste0("Ranked piCs (1-",nrow(INDT),")"),side=1,col="black",line=3,font=2 )

 dev.off()

 return( data.table( XX=seq(nrow(INDT)),YL=CS,YR=VS ) )

}
#<

if (T) {

dir.create( path="PDFs-CumSumFcs/",showWarnings=F )

outCFtab.list <- list()
qqq <-
 plotCF(
  INDT=PST.mrg05kb.frac.dt,PLOT.VAR="RPM", PERC=0.90,XLIM=c(0,100),BY.YR=1000,
  TIT=paste0("DL : piRNA clusters (",RLRANGE[1],"-",RLRANGE[2],"nt)"),
  OUT.PDF="PDFs-CumSumFcs/CFplot-piC-mrg5kB-DLovo3reps-RPM-90thPerc"
 )
outCFtab.list <- c( outCFtab.list,list( DLovo3reps_RPM_90perc= qqq ) )
qqq <-
 plotCF(
  INDT=PST.mrg05kb.frac.dt,PLOT.VAR="RPM", PERC=0.95,XLIM=c(0,100),BY.YR=1000,
  TIT=paste0("DL : piRNA clusters (",RLRANGE[1],"-",RLRANGE[2],"nt)"),
  OUT.PDF="PDFs-CumSumFcs/CFplot-piC-mrg5kB-DLovo3reps-RPM-95thPerc"
 )
outCFtab.list <- c( outCFtab.list,list( DLovo3reps_RPM_95perc= qqq ) )
qqq <-
 plotCF(
  INDT=PST.mrg05kb.frac.dt,PLOT.VAR="RPKM",PERC=0.90,XLIM=c(0,180),BY.YR=0500,
  TIT=paste0("DL : piRNA clusters (",RLRANGE[1],"-",RLRANGE[2],"nt)"),
  OUT.PDF="PDFs-CumSumFcs/CFplot-piC-mrg5kB-DLovo3reps-RPKM-90thPerc"
 )
outCFtab.list <- c( outCFtab.list,list( DLovo3reps_RPKM_90perc=qqq ) )
qqq <-
 plotCF(
  INDT=PST.mrg05kb.frac.dt,PLOT.VAR="RPKM",PERC=0.95,XLIM=c(0,180),BY.YR=0500,
  TIT=paste0("DL : piRNA clusters (",RLRANGE[1],"-",RLRANGE[2],"nt)"),
  OUT.PDF="PDFs-CumSumFcs/CFplot-piC-mrg5kB-DLovo3reps-RPKM-95thPerc"
 )
outCFtab.list <- c( outCFtab.list,list( DLovo3reps_RPKM_95perc=qqq ) )

}

if (T) {

library(openxlsx)

#> create a new open work-book
 xwb <- createWorkbook()
#> add new sheet to work-book "xwb" : max 31 characters !!!
 addWorksheet(   wb=xwb,sheetName="DLovo3reps_RPM"  )
 addWorksheet(   wb=xwb,sheetName="DLovo3reps_RPKM" )
#> close work-book "xwb"
 writeDataTable( wb=xwb,sheet=    "DLovo3reps_RPM" , x=outCFtab.list[["DLovo3reps_RPM_90perc"]], colNames=T,rowNames=F,headerStyle=createStyle(textDecoration="bold") )
 writeDataTable( wb=xwb,sheet=    "DLovo3reps_RPKM", x=outCFtab.list[["DLovo3reps_RPKM_90perc"]],colNames=T,rowNames=F,headerStyle=createStyle(textDecoration="bold") )

#> save work-book "xwb" into file
 saveWorkbook(   wb=xwb,file=paste0("PDFs-CumSumFcs/CFdata-piRNAclusters-27to30nt-",pepDate(),".xlsx"),overwrite=T )

}

