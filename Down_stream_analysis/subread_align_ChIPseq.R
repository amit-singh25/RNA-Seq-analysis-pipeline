rm(list=ls()) 

binquire <- function(packagename) {
  pckgname <- toString(substitute(packagename))
  eval(substitute(
    if (!require(pckgname)) {
      source("https://bioconductor.org/biocLite.R")
      biocLite(pckgname)
      require(pckgname)
    }
  ))
}

inquire <- function(packagename) {
  pckgname <- toString(substitute(packagename))
  eval(substitute(
    if (!require(pckgname)) {
      install.packages(pckgname)
      require(pckgname)
    }
  ))
}

binquire(rbamtools)
binquire(Rsamtools)
binquire(GenomicAlignments)
binquire(rtracklayer)

require(Rsubread)

refpath <-  "~/Chipseq_raw_data/mm9.fasta"
gtfpath <-  "~/Chipseq_raw_data/mm9.gtf"
indexdir <- "~/Chipseq_raw_data/genome/subread-index-mm9"
indir <- "~/Chipseq_raw_data" # where are the fastq files?
align.outdir <- "~/Chipseq_raw_data/alignment/Chipseq/subread" # where is output directory for mapping?

# ensure existence of output directories
dir.create(indexdir, recursive = TRUE)
dir.create(align.outdir, recursive = TRUE)

##########################################################
# build subread index (once done, not necessary)
##########################################################
# jump to target index folder specified above
current_dir <- getwd()
setwd(dir=indexdir)
buildindex(basename="reference_index", reference=refpath)
setwd(current_dir) # jump to original directory

##########################################################
# detect fastq files for paired-end read
##########################################################
reads_R1 <- list.files(path = indir, pattern = "*-R1.fastq")
reads_R2 <- list.files(path = indir, pattern = "*-R2.fastq")

##########################################################
# mapping/alignment using index
##########################################################
align(index=file.path(indexdir, "reference_index"), 
      readfile1 = file.path(indir, reads_R1),
      readfile2 = file.path(indir, reads_R2),
      type = "dna", # or 0 for RNA, 1 for DNA or "dna"
      output_format = "BAM",
      output_file = paste(align.outdir, "/", gsub("-R1.fastq", "", reads_R1), "-align.bam", sep=''),
      phredOffset = 33,
      nthreads = 6)

##########################################################
# sort bam
##########################################################
binquire(rbamtools)
bamfiles <- dir(align.outdir, "bam$", full.names = T)
bamnames <- dir(align.outdir, "bam$")

for (i in 1:length(bamfiles)) {
  bampath <- bamfiles[i]
  bamname <- bamnames[i]
  bam.obj <- bamReader(bampath,idx=FALSE,verbose=0)
  bam_prefix <- gsub(".bam", "-sorted", bamname)
  bamSort(bam.obj, prefix=bam_prefix, byName = FALSE, maxmem = 1e+9)
}

##########################################################
# bam 2 bigwig
##########################################################
binquire(GenomicAlignments)
binquire(rtracklayer)
bamsortedfiles <- dir(align.outdir, "-sorted.bam$", full.names = T)
bamsortednames <- dir(align.outdir, "-sorted.bam$")

for (i in 1:length(bamsortedfiles)) {
  bampath <- bamsortedfiles[i]
  bamname <- bamsortednames[i]
  alignment <- readGAlignments(bampath)
  reads_coverage <- coverage(alignment)
  bigwigpath <- gsub(".bam", "-bigWig.bw", bampath)
  export.bw(reads_coverage, con = bigwigpath)
}




