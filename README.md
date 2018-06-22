- specification: https://samtools.github.io/hts-specs/htsget.html
- python client: https://htsget.readthedocs.io/
- ga4gh news: https://www.ga4gh.org/news/-9msBlhISDK_ltjA7Vt6aA.article

Example

```
if (!"Rhtsget" %in% rownames(installed.packages()))
    BiocManager::install("Bioconductor/Rhtsget")
library(Rhtsget)
gr <- GenomicRanges::GRanges("chr12:111766922-111817529")
sample_id <- "platinum/NA12878"
endpoint <- "https://htsnexus.rnd.dnanex.us/v1"
fl <- tempfile(fileext=".bam")
bam <- htsget_reads(gr, sample_id, endpoint, fl)
Rsamtools::countBam(bam)
```

Two endpoints

- GET /reads/<id>
- GET /variants/<id>

Note from the specification: 

We use the following pan-GA4GH standards:
- 0 start, half open coordinates
- The structuring of POST inputs, redirects and other non-reads data will be protobuf3 compatible JSON
