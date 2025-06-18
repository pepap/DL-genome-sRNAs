library(data.table)
library(seqinr)
library(openxlsx)
library(dplyr)
library(IRanges)
library(GenomicRanges)
library(GenomicAlignments)
library(Rsamtools)
library(Biostrings)
library(GenomeInfoDb)
library(BSgenome)
library(rtracklayer)

# @pepap : installation & documentation https://github.com/HaaseLab/PICB ( "remotes::install_github("HaaseLab/PICB")" worked for me )
library(PICB)

# @pepap : generate the actual date in the text format "YYYYMMDD"
pepDate <- function() { return( format(Sys.time(), "%Y%m%d") ) }

INPUT_GENOME_FASTA="/storage/brno12-cerit/home/pepap/brno1/TobiasBer/06.Mexico-genoma/26-Genoma/derLae1_hic.FINAL.fasta"
STOR="/storage/brno12-cerit/home/pepap/brno1/TobiasBer/08.DL.new.20241117/UNCOLLAPSED_BAMS/"
SUFF=".se.sortCrd.bam"
#!#SUFF=".se.Aligned.sortedByCoord.out.bam"
# @pepap : BAM root name. The expected full path of a individual BAM file : STOR,BAMS[i],SUFF
BAMS <- c(
 "DLfoot_mrg3","DLhead_mrg3","DLjuv_mrg3","DLovo_mrg3"
)
# @pepap : the range of the read widths used in the clustering ( I expect that Petr will ask to adjust it )
READ_WIDTH_RANGE=c(27,30)

# @pepap : load the genome ( ~ 6 minutes )
if ( !exists("myGenome") ) {
cat( "\n ++ Loading Input Genome Fasta ++\n",sep="" )
myGenome     <- PICBloadfasta( FASTA.NAME=INPUT_GENOME_FASTA )
}

# @pepap : load & analyze the individual BAM files ( ~ 46 minutes for the 12 BAM files )
cat( " ++ Loading Input Alignments && running PICB ++\n",sep="" )
OUT.CLUSTERS <- list()
for ( ib in seq_along(BAMS) ) {

cat( "\n  ** ",BAMS[ib]," **\n\n",sep="" )
myAlignments <-
 PICBload(
  BAMFILE                  = paste0(STOR,BAMS[ib],SUFF),
  REFERENCE.GENOME         = myGenome,
# @pepap : TRUE keeps only reads with "M" in cigar-record. Splicing in not allowed in the STAR alignment and I wanted to keep soft-clipped reads ( == cigar contains "S" ) => I turned it off
  SIMPLE.CIGAR             = F,
# @pepap : for standard annotation, turned off just to be sure
  STANDARD.CONTIGS.ONLY    = F,
# @pepap : allow read-width filtering + the range on the next line 
  USE.SIZE.FILTER          = T,
  READ.SIZE.RANGE          = READ_WIDTH_RANGE,
# @pepap : necessary for the "COMPUTE.1U.10A.FRACTIONS=TRUE" option, otherwise ends up with error
  GET.ORIGINAL.SEQUENCE    = T
 )
myClusters   <-
 PICBbuild(
  IN.ALIGNMENTS                               = myAlignments,
  REFERENCE.GENOME                            = myGenome,
# @pepap : default parameters which supposed to be optimize
  UNIQUEMAPPERS.SLIDING.WINDOW.WIDTH          = 350,
  PRIMARY.MULTIMAPPERS.SLIDING.WINDOW.WIDTH   = 350,
  SECONDARY.MULTIMAPPERS.SLIDING.WINDOW.WIDTH = 1000,
  MIN.OVERLAP                                 = 5,
# @pepap : include possible identifier of the ping-pong mechanism
  COMPUTE.1U.10A.FRACTIONS                    = T
 )
tmp.list        <- list( XXX=myClusters )
# @pepap : the R-object "myClusters" is a list containing "seeds","cores" and "clusters" slots, the separate slots are GRanges
# @pepap : all the resulted lists are save as a one comprehensive list - names are the "BAMS" names
names(tmp.list) <- BAMS[ib]
OUT.CLUSTERS    <- c( OUT.CLUSTERS,tmp.list )

}

save( OUT.CLUSTERS,file=paste0("PICB-OUT.CLUSTERS-",READ_WIDTH_RANGE[1],"-",READ_WIDTH_RANGE[2],"nt--",pepDate(),".rda") )

