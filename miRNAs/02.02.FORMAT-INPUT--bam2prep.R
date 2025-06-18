library(data.table)
library(Biostrings)
library(GenomicAlignments)

STOR="/storage/brno12-cerit/home/pepap/brno1/TobiasBer/10.DLgenome-CZ-1.0/BAM.smallRNAseq/00.stats/04.mireap/01.merge_input_reads/"
BAMFILES="03.all-DL/all-DL.se.Aligned.sortedByCoord.out.bam"

RLRAN=c(21,23)
MINRFREQ=50
#TABDELIMINATOR="\t"
#TABDELIMINATOR=" "

#>>> FUN
gr2smfa <- function( igr ) {

 readID <- mcols(igr)[["qname"]]
 readNN <- as.numeric(sub( "^.*[-]","",readID ))
 readID <- as.numeric(sub( "[-].*$","",readID ))
 NSIZE  <- nchar(max(readID))
 readID <- sprintf( fmt=paste0("t%0",NSIZE,"d"),readID )
 readHH <- paste( readID,readNN,sep=" " )
 readSS <- mcols(igr)[["seq"]]
 readSS[ as.character(strand(igr))=="-" ] <- reverseComplement(readSS[ as.character(strand(igr))=="-" ]) 
 names( readSS ) <- readHH

 out.dss <- readSS[ order(names(readSS)) ]

 return( unique(out.dss) )

}
gr2map  <- function( igr ) {

 readID <- mcols(igr)[["qname"]]
 readID <- as.numeric(sub( "[-].*$","",readID ))
 NSIZE  <- nchar(max(readID))
 readID <- sprintf( fmt=paste0("t%0",NSIZE,"d"),readID )
 readCH <- as.character(seqnames(igr))
 readBB <-                 start(igr)
 readEE <-                   end(igr)
 readST <- as.character(  strand(igr))

 out.dt <- data.table( read_ID=readID,chr_ID=readCH,start=readBB,end=readEE,strand=readST )
# out.dt <- out.dt[ order(read_ID) ]
 out.dt <- out.dt[ order(chr_ID,start,end) ]

 return( out.dt )

}
#<<<

cat( "\n >> reading BAM file ...\n",sep="" )
tmp.ga <- readGAlignments( file=paste0(STOR,BAMFILES),param=ScanBamParam( what=c("qname","seq"),tag=c("HI","NH","nM") ) )
tmp.ga <- tmp.ga[ njunc(tmp.ga)==0 ]

cat( " >> converting GAlignment to GRanges ...\n",sep="" )
tmp.gr <- as( tmp.ga,"GRanges" )
tmp.gr <- tmp.gr[ mcols(tmp.gr)[["qname"]] %in% unique(mcols(tmp.gr[ ( width(tmp.gr)>=RLRAN[1] ) & ( width(tmp.gr)<=RLRAN[2] ) ])[["qname"]]) ]

tmp.gr <- tmp.gr[ mcols(tmp.gr)[["qname"]] %in% unique(mcols(tmp.gr[ as.numeric( gsub( "^.*[-]","",mcols(tmp.gr)[["qname"]] ) ) >= MINRFREQ ])[["qname"]]) ]

cat( " >> formatting smfa file ...\n",sep="" )
tmp.ofas <- gr2smfa( igr=tmp.gr )
writeXStringSet( x=tmp.ofas,filepath="all.smrna.fa",compress=F,format="FASTA",append=F )
cat( " >> formatting map file ...\n",sep="" )
tmp.omap <- gr2map(  igr=tmp.gr )
fwrite( x=tmp.omap,file="all.map.txt",sep="\t",quote=F,row.names=F,col.names=F,append=F )

