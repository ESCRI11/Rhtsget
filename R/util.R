.is_single_string <-
    function(x)
{
    is.character(x) && length(x) == 1L && !is.na(x)
}

.is_single_granges <-
    function(x)
{
    is(x, "GRanges") && length(x) == 1L
}

.is_nonexistent_file <-
    function(x)
{
    .is_single_string(x) && !file.exists(x)
}
