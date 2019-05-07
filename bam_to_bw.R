library(GenomicAlignments)
library(rtracklayer)
library(chipseq)
library(stringr)

args <- commandArgs(trailingOnly = T)
message(args)
if(length(args) > 2){stop("Give a maximum of 1 argument")}
if(!(args[1] %in% c("stranded", "notstranded"))){
        stop("Invalid strandedness parameter, allowed are stranded, nonstranded")
} else{sp <- args[1]}

#if (sp == "both") {sp <- c("stranded", "notstranded")}

chroms <- paste0("chr", 1:22)
chroms <- c(chroms, "chrX", "chrY")
strandtype <- c("+", "-")

bamfiles <- list.files(getwd(), pattern = "bam$")
message(paste0("Files are: ", bamfiles))

if (length(bamfiles) == 0) {stop("No bam files detected")}

for (i in seq_along(bamfiles)){
        message(paste0("Loading file ", bamfiles[i]))
        current_bm <- readGAlignments(bamfiles[i])
        message("Processing...")
        if ("stranded" %in% sp){
                current_list <- list(pos = current_bm[strand(current_bm) == "+"], neg = current_bm[strand(current_bm) == "-"])
                name_list <- c("_pos", "_neg")
        } else if ("notstranded" %in% sp){
                current_list <- list(current_bm)
                name_list <- "_full"
        }

        #current_list <- lapply(current_list, function(x) {seqlevelsStyle(x) <- "UCSC"})
        current_list <- lapply(current_list, coverage)
        for (j in seq_along(current_list)){
                namebase <- str_remove(bamfiles[i], ".bam") 
                export.bw(current_list[[j]], paste0(namebase, name_list[j], ".bw")) 

        }}

~                                                                                            
