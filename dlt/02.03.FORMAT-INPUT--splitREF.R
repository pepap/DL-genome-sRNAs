library(Biostrings)

qqq        <- readDNAStringSet("/PATH/TO/INPUT/GENOMIC/FASTA")
names(qqq) <- gsub( pattern="_",replacement=".",x=names(qqq) )

if ( sum(duplicated(names(qqq)))!=0 ) { stop( "\n\n !!! duplicated names in the reference fasta !!!\n\n" ) }

dir.create( path="SINGLEREF/",showWarnings=F )

j=1
for( iNAME in names(qqq) ) {
 cat( sprintf( "%6s / %6s\n",j,length(qqq) ),sep="" )
 writeXStringSet( x=qqq[iNAME],filepath=paste0("SINGLEREF/",iNAME,".fa"),append=F,compress=F,format="FASTA" )
 j <- (j+1)
}
cat("\n",sep="")

