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
#' @param fields `character()` of BAM fields to be included; default
#'     (`NULL`) includes all fields.
#' @param destination `character(1)` file for output; must not exist.
#'
#' @return `htsget_reads()` return the path to a BAM file containing
#'     the requested information.
#'
#' @examples
#' gr <- GenomicRanges::GRanges("chr12:111766922-111817529")
#' sample_id <- "platinum/NA12878"
#' url <- "https://htsnexus.rnd.dnanex.us/v1"
#' fl <- tempfile(fileext=".bam")
#' bam <- htsget_reads(gr, sample_id, url, destination = fl)
#' Rsamtools::countBam(bam)
#'
#' @importMethodsFrom GenomicRanges seqnames start end
#' @export
htsget_reads <-
    function(granges, sample_id, url, token, fields = NULL, destination)
{
    stopifnot(
        .is_single_granges(granges), .is_single_string(sample_id),
        .is_single_string(url), .is_nonexistent_file(destination),
        all(fields %in% names(.BAM_fields))
    )

    fields <-
        if (is.null(fields)) {
            ""
        } else {
            sprintf("&fields=%s", paste(fields, collapse=","))
        }
    queries <- sprintf(
        "%s/files/%s?format=BAM&referenceName=%s&start=%s&end=%s%s",
        url, sample_id,
        seqnames(granges), start(granges), end(granges),
        fields
    )
    content <- .htsget(queries[[1]], token)
    .as_file(content$urls, destination)
    }


#' @rdname htsget_get_token
#'
#' @return
#' @export
#'
#' @examples
htsget_get_token <- function(url, user, pass){

    res <- httr::POST(url, add_headers("Content-Type" = "application/x-www-form-urlencoded"),
                      body = list(grant_type = "password",
                                  client_id = "f20cd2d3-682a-4568-a53e-4262ef54c8f4",
                                  scope = "openid",
                                  client_secret = "AMenuDLjVdVo4BSwi0QD54LL6NeVDEZRzEQUJ7hJOM3g4imDZBHHX0hNfKHPeQIGkskhtCmqAJtt_jm7EKq-rWw",
                                  username = user,
                                  password = pass), encode = "form")
    
    return(paste("Bearer", content(res)$access_token))
}


#' @rdname htsget
#'
#' @return `BAMfields()` returns a data.frame describing available
#'     fields.
#'
#' @examples
#' BAMfields()
#' fl <- tempfile(fileext=".bam")
#' fields <- c("RNAME", "POS", "CIGAR")
#' bam <- htsget_reads(gr, sample_id, url, fields, destination = fl)
#' names(Rsamtools::scanBam(bam)[[1]])
#'
#' @export
BAMfields <-
    function()
{
    as.data.frame(.BAM_fields)
}

.BAM_fields <- c(
    QNAME = "Read names",
    FLAG = "Read bit flags",
    RNAME = "Reference sequence name",
    POS = "Alignment position",
    MAPQ = "Mapping quality score",
    CIGAR = "CIGAR string",
    RNEXT = "Reference sequence name of the next fragment template",
    PNEXT = "Alignment position of the next fragment in the template",
    TLEN = "Inferred template size",
    SEQ = "Read bases",
    QUAL = "Base quality scores"
)

#' @rdname htsget
#'
#' @return `htsget_variants()` returns a `character(1)` path to a VCF file
#'     containing requested information.
#'
#' @examples
#' gr <- GenomicRanges::GRanges("12:112204691-112247789")
#' sample_id <- "1000genomes/20130502_autosomes"
#' url <- "https://htsnexus.rnd.dnanex.us/v1"
#' fl <- tempfile(fileext=".vcf")
#' vcf <- htsget_variants(gr, sample_id, url, destination = fl)
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
