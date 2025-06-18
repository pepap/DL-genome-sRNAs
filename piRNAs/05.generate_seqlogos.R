library(GenomicRanges)
library(GenomicAlignments)
library(data.table)
library(pbapply)
library(seqLogo)

# @pepap : load piC ranges in the format of GRanges (manually curated & merged)
load( "PST.mrg05kb.frac.gr.rda",verbose=T )

# @pepap : R-object name of the piC ranges
ANNOTGR="PST.mrg05kb.frac.gr"
# @pepap : minimal overlap for "findOverlaps" command
MINOVRL=20
# @pepap : merged UNCOLLAPSED ovotestis replicates
BAMFILE="/storage/brno12-cerit/home/pepap/brno1/TobiasBer/08.DL.new.20241117/UNCOLLAPSED_BAMS/DLovo_mrg3.se.sortCrd.bam"
# @pepap : read length range - used for seqLogos
RLRANGE=c(27,30)
# @pepap : read length range - used for normalisation
NORMRAN=c(18,32)
# @pepap : which part of the read is going to be visuelized in the output seqLogo figure. Default : nucleotides 1 .. 15
RCRD=c(1,15)
# @pepap : the units of the y-axes : Informational Content (IC) = <0;2> in Bits
IC=T

cat( "\n ** Reading BAM file ...\n",sep="" )
if ( !exists("tmp.gr") ) {
 tmp.tga  <- readGAlignments( file=BAMFILE,param=ScanBamParam( tag=c("NH","HI","nM"),what=c("qname","seq") ),use.names=T )
 cat( "   > Count all 18-32nt primary reads ...\n",sep="" )
 TMP.NORM <- length( tmp.tga[ ( width(tmp.tga)>=NORMRAN[1] ) & ( width(tmp.tga)<=NORMRAN[2] ) & ( mcols(tmp.tga)[["HI"]]==1 ) ] )
 cat( "   > Removing spliced reads ... \n",sep="" )
 tmp.tga  <- tmp.tga[ njunc(tmp.tga)==0 ]
 cat( "   > Converting GA => GR object ...\n",sep="" )
 tmp.gr   <- as( object=tmp.tga,Class="GRanges" )
 mcols(tmp.gr)[["cigar"]]  <- cigar(tmp.tga)
 mcols(tmp.gr)[["rcount"]] <- as.numeric( sub( "[.][0-9]*$","",sub( "^[0-9]*[-]","",mcols(tmp.gr)[["qname"]] ) ) )
 norm.DLovo_mrg3.dt <- data.table( all_18to32nt=TMP.NORM,all_27to30nt=length( tmp.gr[ ( width(tmp.gr)>=RLRANGE[1] ) & ( width(tmp.gr)<=RLRANGE[2] ) & ( mcols(tmp.gr)[["HI"]]==1 ) ] ) )
}

cat( " ** Selecting reads by width ...\n",sep="" )
sele.all.gr      <- tmp.gr[ ( width(tmp.gr)>=RLRANGE[1] ) & ( width(tmp.gr)<=RLRANGE[2] ) ]
# @pepap : extract only mapped nucleotides, remove soft-clipped nucleotides
seq.irl          <- cigarRangesAlongQuerySpace( cigar=mcols(sele.all.gr)[["cigar"]],after.soft.clipping=F,ops=c("M","I","D","H","P","=","X"),with.ops=T,reduce.ranges=T )
sele.all.dt      <-
 data.table(
  seq    = as.character( unlist(extractAt(          x=mcols(sele.all.gr)[["seq"]],at=seq.irl )) ),
  strand = as.character(strand(sele.all.gr)),
  cigar  = mcols(sele.all.gr)[["cigar"]],
  rname  = sub( "[.][0-9]*$","",mcols(sele.all.gr)[["qname"]] ),
  rcount = mcols(sele.all.gr)[["rcount"]]
 )
sele.all.dt[ strand=="-" ][["seq"]] <- as.character(reverseComplement(DNAStringSet(sele.all.dt[ strand=="-" ][["seq"]])))

cat( " ** Selecting reads by overlaps ...\n",sep="" )
tmp.ol           <- findOverlaps( query=sele.all.gr,subject=get(ANNOTGR),minoverlap=MINOVRL,ignore.strand=F )
sele.piC.gr      <- sele.all.gr[ queryHits(tmp.ol) ]
sele.piC.dt      <- sele.all.dt[ queryHits(tmp.ol) ]

cat( " ** Removing duplicates ...\n",sep="" )
sele.all.dt      <- sele.all.dt[ duplicated( sele.all.dt[,paste0(seq,":",rname)] ) ]
sele.piC.dt      <- sele.piC.dt[ duplicated( sele.piC.dt[,paste0(seq,":",rname)] ) ]

cat( " ** Calc logos ...\n",sep="" )
all.nfm <-
 pbsapply(
  X   = seq( RCRD[1],RCRD[2] ),
  FUN = function(X){
   nfm.mat <- oligonucleotideFrequency( x=DNAStringSet( x=sele.all.dt[["seq"]],start=X,end=X ),width=1 )*sele.all.dt[["rcount"]]
   return( colSums(nfm.mat)/sum(colSums(nfm.mat)) )
  }
 )
cat( "\n",sep="" )
all.pwm               <- makePWM( pwm=all.nfm,alphabet="RNA" )
colnames(all.pwm@pwm) <- as.character( seq(RCRD[1],RCRD[2]) )
piC.nfm <-
 pbsapply(
  X   = seq( RCRD[1],RCRD[2] ),
  FUN = function(X){
   nfm.mat <- oligonucleotideFrequency( x=DNAStringSet( x=sele.piC.dt[["seq"]],start=X,end=X ),width=1 )*sele.piC.dt[["rcount"]]
   return( colSums(nfm.mat)/sum(colSums(nfm.mat)) )
  }
 )
cat( "\n",sep="" )
piC.pwm               <- makePWM( pwm=piC.nfm,alphabet="RNA" )
colnames(piC.pwm@pwm) <- as.character( seq(RCRD[1],RCRD[2]) )


dir.create( path="PDFs-seqLogos-refine/",showWarnings=F )

#if ( IC ) {
#1#pdf("PDFs-seqLogos/all.27to30nt.IC.allRead.pdf",height=5,width=10)
pdf("PDFs-seqLogos-refine/all.27to30nt.IC.rmvDups.pdf",height=5,width=10)
 barplot(
  t(t( all.pwm@pwm )*all.pwm@ic),ylim=c(0,2),ylab="Information content",xlab="Position in read",col=c(A='#61D04F',C='#2297E6',G='#F5C710',U='#DF536B'),
  legend.text=c("A","C","G","U"),args.legend=list( bty="n" )
 )
dev.off()
#} else    {
#1#pdf("PDFs-seqLogos/all.27to30nt.PP.allRead.pdf",height=5,width=10)
pdf("PDFs-seqLogos-refine/all.27to30nt.PP.rmvDups.pdf",height=5,width=10)
 barplot(
       all.pwm@pwm              ,ylim=c(0,1),ylab="Probability",        xlab="Position in read",col=c(A='#61D04F',C='#2297E6',G='#F5C710',U='#DF536B'),
  legend.text=c("A","C","G","U"),args.legend=list( bty="n" )
 )
dev.off()
#}

#if ( IC ) {
#1#pdf("PDFs-seqLogos/piC.27to30nt.IC.allRead.pdf",height=5,width=10)
pdf("PDFs-seqLogos-refine/piC.27to30nt.IC.rmvDups.pdf",height=5,width=10)
 barplot(
  t(t( piC.pwm@pwm )*piC.pwm@ic),ylim=c(0,2),ylab="Information content",xlab="Position in read",col=c(A='#61D04F',C='#2297E6',G='#F5C710',U='#DF536B'),
  legend.text=c("A","C","G","U"),args.legend=list( bty="n" )
 )
dev.off()
#} else    {
#1#pdf("PDFs-seqLogos/piC.27to30nt.PP.allRead.pdf",height=5,width=10)
pdf("PDFs-seqLogos-refine/piC.27to30nt.PP.rmvDups.pdf",height=5,width=10)
 barplot(
       piC.pwm@pwm              ,ylim=c(0,1),ylab="Probability",        xlab="Position in read",col=c(A='#61D04F',C='#2297E6',G='#F5C710',U='#DF536B'),
  legend.text=c("A","C","G","U"),args.legend=list( bty="n" )
 )
dev.off()
#}

pdf( file="PDFs-seqLogos-refine/Nucleotides--all.27to30nt.IC.rmvDups.pdf" )
 qqq <- seqLogo::seqLogo( pwm=all.pwm,ic.scale=T )
 print(qqq)
dev.off()
pdf( file="PDFs-seqLogos-refine/Nucleotides--piC.27to30nt.IC.rmvDups.pdf" )
 qqq <- seqLogo::seqLogo( pwm=piC.pwm,ic.scale=T )
 print(qqq)
dev.off()

