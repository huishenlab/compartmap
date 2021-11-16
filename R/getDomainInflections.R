#' A wrapper function to generate a GRanges object of chromatin domain inflection points
#'
#' @name getDomainInflections
#'
#' @param gr Input GRanges object with mcols column corresponding to chromatin domains
#' @param what The name of the column containing the chromatin domain information
#' @param res What resolution the domains were called
#' @param chrs Which chromosomes to work on
#' @param genome Which genome does the input data come from
#'
#' @return A GRanges object of compartment inflection points
#' 
#' @import GenomicRanges
#' @export
#'
#' @examples
#' data("k562_scrna_chr14", package = "compartmap")
#' chr14_domains <- scCompartments(k562_scrna_chr14,
#'                                 res = 1e6, genome = "hg19",
#'                                 group = TRUE, bootstrap = FALSE)
#' chr14_domain_inflections <- getDomainInflections(chr14_domains, what = "pc")

getDomainInflections <- function(gr, what = "score", res = 1e6,
                                      chrs = c(paste0("chr", 1:22), "chrX"),
                                      genome = c("hg19", "hg38", "mm9", "mm10")) {
  #find the compartment inflection points
  stopifnot(is(gr, "GenomicRanges"))
  stopifnot(what %in% names(mcols(gr)))
  #determine which genome we are working with
  seqlengths <- getSeqLengths(genome = genome, chr = chrs)
  
  #we may not be able to assume continuous compartment structure here
  #so somehow, we need to find continuous runs
  message("Tiling genome.")
  tiles <- tileGenome(seqlengths = seqlengths,
                      tilewidth = res,
                      cut.last.tile.in.chrom = TRUE)
  #reset to 0-based
  start(tiles) <- suppressWarnings(start(tiles) - 1)
  end(tiles) <- suppressWarnings(end(tiles) - 1)
  
  #add a column for continuous runs!
  mcols(tiles)$run <- seq(1:length(tiles))
  mcols(tiles)$score <- NA
  
  #overlap
  message("Finding overlaps.")
  ol <- findOverlaps(tiles, gr)
  
  #tack on the eigenvalues
  message("Injecting eigenvalues.")
  mcols(tiles)$score[queryHits(ol)] <- mcols(gr)[[what]][subjectHits(ol)]
  tiles.sub <- tiles[!is.na(mcols(tiles)$score),]
  
  #look at continuous sequence only
  #this is the tricky part...
  d <- diff(mcols(tiles.sub)$run)
  #theoretically, continuous sequence will always be 1
  contig <- which(d == 1)
  
  #need to do this in chunks
  #technically, we really only need the places that are non-contiguous
  #to bookend our vectors to calculate
  non_contig <- which(d != 1)
  
  #convert to pos and negative signs
  .getInflections <- function(gr, what = "score") {
    gr.signs <- sign(mcols(gr)[[what]])
    #find the inflected compartment
    gr.signs.bool <- gr.signs < 0
    #big list of things that are less than zero
    #where the flip will be FALSE
    gr.inflect <- gr[!gr.signs.bool,]
    if (length(gr.inflect) == 0) {
      return(GRanges())
    }
    #merge them
    #shift
    end(gr.inflect) <- end(gr.inflect) + 1
    #reduce
    gr.inflect <- reduce(gr.inflect)
    #shift back
    end(gr.inflect) <- end(gr.inflect) - 1
    #gr.signs.diff <- diff(gr.signs) != 0
    #subset the original granges
    #gr.inflect <- gr[gr.signs.diff,]
    if (length(gr.signs) == 2 & any(gr.signs.bool)) {
      #special case
      gr.inflect <- gr[2,]
      gr.inflect.new <- GRanges(seqnames = seqnames(gr.inflect),
                                ranges = IRanges(start = start(gr.inflect),
                                                 end = start(gr.inflect)),
                                strand = "*")
      return(gr.inflect.new)
    }
    gr.inflect.new.start <- GRanges(seqnames = seqnames(gr.inflect),
                                    ranges = IRanges(start = start(gr.inflect),
                                                     end = start(gr.inflect)),
                                    strand = "*")
    gr.inflect.new.end <- GRanges(seqnames = seqnames(gr.inflect),
                                  ranges = IRanges(start = end(gr.inflect),
                                                   end = end(gr.inflect)),
                                  strand = "*")
    gr.inflect.new <- sort(c(gr.inflect.new.start,
                             gr.inflect.new.end))
    #NOTE: the start of the inflected compartment is what we want
    return(gr.inflect.new)
  }
  
  #if we have fully continuous data, short circuit
  if (!length(non_contig)) {
    message("Contiguous runs. Finding inflections.")
    grl.new <- .getInflections(gr, what = what)
    return(grl.new)
  }
  
  #loop through the non-contiguous space
  grl <- lapply(1:length(non_contig), function(i) {
    message("Block processing non-contiguous space for block ", i)
    #get block of contiguous sequence
    if (non_contig[i] == non_contig[1]) {
      block <- tiles.sub[contig[1]:non_contig[i],]
    } else {
      block <- tiles.sub[(non_contig[i-1] + 1):non_contig[i],]
    }
    return(block)
  })
  
  #we now need to add the final block to the grl
  #but what if we have a non-contig block at the end all by itself?
  fin.block <- tiles.sub[(tail(non_contig, n=1L) + 1):length(tiles.sub),]
  grl <- c(grl, fin.block)
  
  #run it
  grl.new <- lapply(grl, function(x) {
    message("Finding inflections.")
    return(.getInflections(x, what = what))
  })
  message("Returning inflections.")
  return(unlist(as(grl.new, "GRangesList")))
}
