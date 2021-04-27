#' @importFrom httr GET stop_for_status content

.htsget <-
    function(query, token)
{
    if(!is.null(token)){
        response <- GET(query, add_headers("Authorization" = token))
    }
    else{
        response <- GET(query)
    }
    stop_for_status(response)
    content(response, type="application/json")[["htsget"]]
}

#' @importFrom httr write_stream add_headers progress

.htsstream_writer <- function(bam) {
    write_stream(function(x) {
        writeBin(x, bam)
    })
}

.htsstream <-
    function(bam, query, headers)
{
    headers <- add_headers(unlist(headers))
    response <- GET(query, headers, .htsstream_writer(bam), progress())
    stop_for_status(response)
    }

.fetch_file <-
    function(query, token)
    {
        if(!is.null(token)){
            response <- GET(query, add_headers("Authorization" = token))
        }
        else{
            response <- GET(query)
        }
        stop_for_status(response)
        content(response, type="text/csv")#[["htsget"]]
    }
