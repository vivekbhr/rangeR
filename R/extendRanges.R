

#' Make bins extending away from a GRanges (which will be anchored at the start)
#'
#' @param anchor_gr GRanges object to anchor on (this becomes the first bin)
#' @param window_size Size of the bins
#' @param num_windows Number of bins
#' @param ntop subset input GRanges to take top N regions (based on "score" column)
#'
#' @return New GRanges object with <anchor>---------> extended windows
#' @export
#'
#' @examples
#' 
#' gr100k <- make_anchored_bins(test_ranges, 10000, 10)
#' 
make_anchored_bins <- function(anchor_gr, window_size, num_windows, ntop=NA) {
  
  ## select top 1000 peaks and intersect with genes
  if(!is.na(ntop)) {
    prom_peaks <- prom_peaks[order(prom_peaks$score, decreasing = T)[1:ntop]]
  }
  ## make 50x 10kb bins starting from promoter peaks
  prom_peaks <- resize(prom_peaks, window_size, fix="center")
  prom_peaks <- prom_peaks[countOverlaps(prom_peaks, prom_peaks) == 1] # remove overlaps
  
  # 250kb window starting from each element
  windows_ext <- num_windows + 6 # to account for overlapping windows that would be removed
  nkeep <- num_windows + 1
  window <- GRanges(seqnames(prom_peaks), IRanges(start(prom_peaks), start(prom_peaks) + (windows_ext*window_size)))
  grl <- tile(window, width = window_size)
  
  # remove windows that overlap with another prom_peak and take first 51
  lapply(seq_along(grl), function(n) {
    subset <- prom_peaks[-n]
    gr <- grl[[n]]
    gr <- gr[!overlapsAny(gr, subset)]
    if(length(gr)>nkeep) {
      return(gr[1:nkeep])
    } else {
      return(gr)
    }
  }) -> grl_extended
  
  ## make bed 
  names(grl_extended) <- as.character(prom_peaks)
  grl_extended <- grl_extended[sapply(grl_extended, length)==nkeep]
  gr_extended <- unlist(GRangesList(grl_extended))
  
  return(gr_extended)
}