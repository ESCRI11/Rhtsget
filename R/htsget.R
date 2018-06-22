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

.as_file <-
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

#' Retrieve reads (BAM files) or variants (VCF files)
#' @rdname htsget
#'
#' @param granges `GRanges()` instance with genomic coordinates for
#'     retrieval. Coordinates follow Bioconductor standards (1-based,
#'     closed intervals) and are translated to GA4GH coordinates
#'     (0-based, half-open intervals).
#' @param sample_id `character(1)` sample identifier
#' @param url `character(1)` data resource url, excluding endpoint and
#'     sample information.
#' @param destination `character(1)` file for output; must not exist.
#'
#' @examples
#' gr <- GenomicRanges::GRanges("chr12:111766922-111817529")
#' sample_id <- "platinum/NA12878"
#' url <- "https://htsnexus.rnd.dnanex.us/v1"
#' fl <- tempfile(fileext=".bam")
#' bam <- htsget_reads(gr, sample_id, url, fl)
#' Rsamtools::countBam(bam)
#'
#' @importMethodsFrom GenomicRanges seqnames start end
#' @export
htsget_reads <-
    function(granges, sample_id, url, destination)
{
    stopifnot(
        .is_single_granges(granges), .is_single_string(sample_id),
        .is_single_string(url), .is_nonexistent_file(destination)
    )

    queries <- sprintf(
        "%s/reads/%s?format=BAM&referenceName=%s&start=%s&end=%s",
        url, sample_id,
        seqnames(granges), start(granges), end(granges)
    )
    content <- .htsget(queries[[1]])
    .as_file(content$urls, destination)
}

#' @rdname htsget
#'
#' @examples
#' gr <- GenomicRanges::GRanges("12:112204691-112247789")
#' sample_id <- "1000genomes/20130502_autosomes"
#' url <- "https://htsnexus.rnd.dnanex.us/v1"
#' fl <- tempfile(fileext=".vcf")
#' vcf <- htsget_variants(gr, sample_id, url, fl)
#' VariantAnnotation::readVcf(vcf)
#' @export
htsget_variants <-
    function(granges, sample_id, url, destination)
{
    stopifnot(
        .is_single_granges(granges), .is_single_string(sample_id),
        .is_single_string(url), .is_nonexistent_file(destination)
    )

    queries <- sprintf(
        "%s/variants/%s?format=VCF&referenceName=%s&start=%s&end=%s",
        url, sample_id,
        seqnames(granges), start(granges), end(granges)
    )
    content <- .htsget(queries[[1]])
    .as_file(content$urls, destination)
}
