isWholeNumber <- function(x) {
    return(x %% 1 == 0)
}


# Shifts to the right all strings of a list with a number of blanks
shift <- \(strL, nBlanks = 0) {
    pre <- paste(rep(" ", nBlanks), collapse = "")
    lapply(strL, \(x) sprintf("%s%s", pre, x)) |> unlist()
}

#' Check and keep valid index only
#'
#' @param display Numeric or character. Display index
#' @param elecNames Character. All electrode names
checkIndex <- function(indices, names) {
    if (length(names) == 0) {
        return()
    }
    if (length(indices) == 0) {
        return()
    }
    if (is(indices, "numeric")) {
        allIndices <- seq_along(names)
        diffIndices <- setdiff(indices, allIndices)
        indicesFiltered <- indices[!indices %in% diffIndices]
        result <- indicesFiltered
    } else {
        diffIndices <- setdiff(indices, names)
        indicesFiltered <- indices[!indices %in% diffIndices]
        result <- which(names %in% indicesFiltered)
    }
    if (length(diffIndices)) {
        indicesMissing <- paste(diffIndices, collapse = ", ")
        indicesExist <- paste(indicesFiltered, collapse = ", ")
        warning(
            glue("Indices {indicesMissing} are out of range. I will keep the valid values {indicesExist}.")
        )
    }
    result
}

# Try to convert a vector to numeric
# If not possible, this does not raise a warning or error
# but returns NULL
tryToNum <- function(x) {
    x <- tryCatch(as.numeric(x), error = function(e) x, warning = function(w) x)
    if (is.numeric(x)) {
        return(x)
    } else {
        return(NULL)
    }
}
