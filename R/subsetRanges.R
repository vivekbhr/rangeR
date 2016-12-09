

#' Filter a GRanges by the direction of genes
#'
#' @param ranges A GRanges object to subset
#' @param direction Direction of Gene-pairs to subset. Options are 'convergent','divergent', or 'tandem'
#'
#' @return A list of two GRanges objects
#' @export
#'
#' @examples
#'
#' divergent <- subsetRange_byGeneDirection(testRanges, direction = "divergent")
#'


subsetRange_byGeneDirection <- function(ranges, direction = "divergent"){

        if(direction == "divergent") {
                # Fix everything to start/end 1bp and check the strand of upstream genes
                ranges.plus <- resize(ranges[GenomicRanges::strand(ranges) == "+"], width = 1, fix = "start")
                upstream.toplus <- follow(ranges.plus, ranges, ignore.strand = TRUE)
                # remove na Entries from query and return entries from subject that have a hit
                ranges.plus <- ranges.plus[-which(is.na(upstream.toplus))]
                upstream.toplus <- ranges[na.omit(upstream.toplus)]
                # out of these, only keep -ve strand entries
                minus <- which(GenomicRanges::strand(upstream.toplus) == "-")
                ranges.minus <- resize(upstream.toplus[minus], width = 1, fix = "start")
                # only keep queries which have -ve strand in subject
                ranges.plus <- ranges.plus[minus]
        } else
                if(direction == "convergent") {
                        # same as above
                        ranges.minus <- resize(ranges[GenomicRanges::strand(ranges) == "-"],width = 1, fix = "end")
                        upstream.tominus <- follow(ranges.minus, ranges, ignore.strand = TRUE)
                        ranges.minus <- ranges.minus[-which(is.na(upstream.tominus))]
                        upstream.tominus <- ranges[na.omit(upstream.tominus)]

                        plus <- which(GenomicRanges::strand(upstream.tominus) == "+")
                        ranges.plus <- resize(upstream.tominus[plus], width = 1, fix = "end")
                        ranges.minus <- ranges.minus[plus]
                } else
                        if(direction == "tandem") {

                                tandem <- lapply(c("+","-"), function(strand){

                                        if(strand == "+") {
                                                fix1 <- "start"
                                                fix2 <- "end"
                                        } else {
                                                fix1 <- "end"
                                                fix2 <- "start"
                                        }

                                        # for + , get + genes upstream and for - get -genes
                                        ranges.plus <- resize(ranges[GenomicRanges::strand(ranges) == strand], width = 1, fix = fix1)
                                        upstream.toplus <- follow(ranges.plus, ranges, ignore.strand = TRUE)
                                        # remove na Entries from query and return entries from subject that have a hit
                                        ranges.plus <- ranges.plus[-which(is.na(upstream.toplus))]
                                        upstream.toplus <- ranges[na.omit(upstream.toplus)]
                                        # out of these, only keep +ve strand entries
                                        plusUP <- which(GenomicRanges::strand(upstream.toplus) == strand)
                                        ranges.plusUP <- resize(upstream.toplus[plusUP], width = 1, fix = fix2)
                                        # only keep queries which have -ve strand in subject
                                        ranges.plus <- ranges.plus[plusUP]
                                        return(list(plus = ranges.plus, UPtoplus = ranges.plusUP))
                                })

                                ranges.plus <- c(tandem[[1]]$plus, tandem[[2]]$plus)
                                ranges.minus <- c(tandem[[1]]$UPtoplus, tandem[[2]]$UPtoplus)
                        }

        return(list(GR.plus = ranges.plus, GR.minus = ranges.minus))

}


#' Merge Ranges from the given list of GRanges
#'
#' @param rlist A list of GRanges
#' @param filterRanges GRanges, If provided , the regions overlapping with these regions are removed
#' @param lengthCutoff Only return regions upto this length (bp)
#'
#' @return GRanges object
#' @export
#'
#' @examples
#'
#' divergent <- subsetRange_byGeneDirection(testRanges, direction = "divergent")
#' divergent_merged <- mergeRanges(divergent)
#'

mergeRanges <- function(rlist, filterGenes = NA, lengthCutoff = 2000) {
        chroms <- seqnames(rlist$GR.plus)
        starts <- pmin(min(ranges(rlist$GR.plus)), min(ranges(rlist$GR.minus)) )
        ends <- pmax(max(ranges(rlist$GR.plus)), max(ranges(rlist$GR.minus)) )
        out <- GRanges(seqnames = chroms, IRanges(starts, ends))
        # make intergenic distance as scores and sort by it
        out$score <- distance(rlist$GR.plus, rlist$GR.minus, ignore.strand=TRUE)
        out <- out[order(out$score)]
        out <- out[out$score <= lengthCutoff]
        # remove regions overlapping with any gene (>5bp)
        if(!is.na(filterGenes)) {
                remove <- GenomicRanges::findOverlaps(out, filterGenes, minoverlap = 5) %>% queryHits()
                out <- out[-remove]
        }
        return(out)
}
