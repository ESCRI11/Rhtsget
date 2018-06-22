#' @importFrom GenomicRanges GRanges
#' @importMethodsFrom GenomicRanges resize seqnames

.to_0_start_halfopen <-
    function(gr)
{
    resize(gr, width(gr) + 1L, fix = "end")
}

.from_0_start_halfopen <-
    function(gr)
{
    resize(gr, width(gr) - 1L, fix = "end")
}

#' @importFrom base64enc base64decode

.as_bam <-
    function(urls, output)
{
    bam <- file(output, "wb")
    on.exit(close(bam))
    for (url in urls) {
        if (startsWith(url$url, "data")) {
            data <- sub("data:application/octet-stream;base64,", "", url$url)
            base64decode(data, bam)
        } else {
            .htsstream(bam, url$url, url$headers)
        }
    }
    output
}

#' Retrieve reads (BAM files)
#'
#' @param granges `GRanges()` instance with genomic coordinates for
#'     retrieval. Coordinates follow Bioconductor standards (1-based,
#'     closed intervals) and are translated to GA4GH coordinates
#'     (0-based, half-open intervals).
#' @param sample_id `character(1)` sample identifier
#' @param endpoint `character(1)` data resource url, excluding
#'     endpoint.
#' @param destination `character(1)` file for output; must not exist.
#'
#' @examples
#' gr <- GRanges("chr12:111766922-111817529")
#' sample_id <- "platinum/NA12878"
#' endpoint <- "https://htsnexus.rnd.dnanex.us/v1"
#' fl <- file.path(tempdir(), paste0(basename(url), ".bam"))
#' unlink(fl)
#' bam <- htsget_reads(gr, sample_id, endpoint, fl)
#' Rsamtools::countBam(bam)
#'
#' @export
htsget_reads <-
    function(granges, sample_id, endpoint, destination)
{
    queries <- sprintf(
        "%s/reads/%s?format=BAM&referenceName=%s&start=%s&end=%s",
        endpoint, sample_id, seqnames(gr), start(gr), end(gr)
    )
    content <- .htsget(queries[[1]])
    .as_bam(content$urls, output)
}
