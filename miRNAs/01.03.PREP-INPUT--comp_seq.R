library(data.table)
library(Biostrings)

STOR="/storage/brno12-cerit/home/pepap/brno1/TobiasBer/02.smallRNAseq-DL-newKit-20240504/FASTQ/CUTADAPT/"
SUFF="_L1-all.re-collapsed.fa.gz"
CFAS=system( command=paste0("/bin/ls ",STOR,"DL*",SUFF),intern=T )

RLENRAN=c(19,32)
MINFREQ=1000

if (T) {
if ( !exists("sSEQ") ) {

sSEQ    <- c()
cat( " > GET SEQ\n",sep="" )
for ( icfa in seq_along(CFAS) ) {

cat( " >> ",sub("^.*[/]","",sub(SUFF,"",CFAS[icfa])),"\n",sep="" )

all.dss <- readDNAStringSet( CFAS[icfa] )
all.cnt <- as.numeric(sub("^[0-9]*[-]","",names(all.dss)))
sel.dss <- all.dss[ width(all.dss)>=RLENRAN[1] & width(all.dss)<=RLENRAN[2] ]
sel.cnt <- as.numeric(sub("^[0-9]*[-]","",names(sel.dss)))

sSEQ    <- unique(c( sSEQ,as.vector(as.character(sel.dss[ sel.cnt>=MINFREQ ])) ))

}
sSEQ    <- data.table( SEQ=sSEQ,WIDTH=nchar(sSEQ) )
}

cat( " > GET DT\n",sep="" )
for ( icfa in seq_along(CFAS) ) {

cat( " >> ",sub("^.*[/]","",sub(SUFF,"",CFAS[icfa])),"\n",sep="" )

all.dss                                     <- readDNAStringSet( CFAS[icfa] )
all.cnt                                     <- as.numeric(sub("^[0-9]*[-]","",names(all.dss)))
all.dt                                      <- data.table( SEQ=as.character(all.dss),TMP=(all.cnt/sum(all.cnt))*1e06 )
colnames(all.dt)[ colnames(all.dt)=="TMP" ] <- sub("^.*[/]","",sub(SUFF,"",CFAS[icfa]))
sSEQ                                        <- merge( sSEQ,all.dt,by="SEQ",sort=F,all.x=T )

}

cat( "\n",sep="" )

sSEQ <- sSEQ[ , lapply( X=.SD,FUN=function(C){ C[is.na(C)] <- as.numeric("0"); return(C) } ) , .SDcols=colnames(sSEQ) ]
save( sSEQ,file=paste0("sSEQ-",pepDate(),".rda") )

}

library(GenomicAlignments)

tmp.tga <- readGAlignments( "03.all-DL/all-DL.se.Aligned.sortedByCoord.out.bam",param=ScanBamParam( what=c("seq"),tag=c("NH","HI","nM") ) )
tmp.gr  <- as( tmp.tga,"GRanges" )
tmp.gr  <- tmp.gr[ mcols(tmp.gr)[["HI"]]==1 ]
tmp.seq <- mcols(tmp.gr)[["seq"]]
tmp.seq[ as.character(strand(tmp.gr))=="-" ] <- reverseComplement(tmp.seq[ as.character(strand(tmp.gr))=="-" ])
crd.dt  <- data.table( SEQ=as.character(tmp.seq),CRD=paste0( as.character(seqnames(tmp.gr)),":",start(tmp.gr),"-",end(tmp.gr) ),NH=mcols(tmp.gr)[["NH"]] )

load( "sSEQ-20250508.rda",verbose=T )

aSEQ <- sSEQ[ , {
 list(
  WIDTH        = WIDTH,
  DL_FOOT      = mean(c(DLfoot1,DLfoot2,DLfoot3)),
  DL_HEAD      = mean(c(DLhead1,DLhead2,DLhead3)),
  DL_JUVENILE  = mean(c(DLjuv1, DLjuv2, DLjuv3 )),
  DL_OVOTESTES = mean(c(DLovo1, DLovo2, DLovo3 ))
 )
                } , by="SEQ" ]
library(pbapply)

aSEQ[["MAX"]] <-
 pbapply::pbapply( X=as.matrix(aSEQ[,c("SEQ","DL_FOOT","DL_HEAD","DL_JUVENILE","DL_OVOTESTES"),with=F],rownames="SEQ"),MARGIN=1,FUN=max )
aSEQ <- merge( aSEQ,crd.dt,by="SEQ",all.x=T,sort=F )
aSEQ <- aSEQ[ !is.na(CRD) ]

sSEQ <- sSEQ[ SEQ %in% aSEQ[["SEQ"]] ]

library(openxlsx)

#> create a new open work-book
xwb <- createWorkbook()
#> add new sheet to work-book "xwb" : max 31 characters !!!
addWorksheet(   wb=xwb,sheetName="avr-smallRNAs-seqs" )
addWorksheet(   wb=xwb,sheetName="frequent-smallRNAs-seqs" )
#> close work-book "xwb"
writeDataTable( wb=xwb,sheet=    "avr-smallRNAs-seqs",      x=aSEQ,colNames=T,rowNames=F,headerStyle=createStyle(textDecoration="bold") )
writeDataTable( wb=xwb,sheet=    "frequent-smallRNAs-seqs", x=sSEQ,colNames=T,rowNames=F,headerStyle=createStyle(textDecoration="bold") )

#> save work-book "xwb" into file
saveWorkbook(   wb=xwb,file=paste0("LI-smallRNAseq-seqs-",pepDate(),".xlsx"),overwrite=T )

