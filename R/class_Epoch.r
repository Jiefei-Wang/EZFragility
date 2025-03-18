#' @title Epoch Class
#' @description S4 class to handle epoch data with electrodes and time points
#' @slot data a tibble containing epoch data (columns=time points, rows=electrodes)
#' @slot times Numeric vector containing time range
.Epoch <- setClass("Epoch",
    slots = list(
        data = "matrix",
        times = "numericOrNULL"
    )
)



#' Constructor for Epoch class
#' @param data Matrix containing epoch data (rows=electrodes, columns=time points)
#' @param electrodes Optional character vector for electrode names, if not provided, column names of data are used. If both are NULL, electrodes are named E1, E2, ...
#' @param timeRanges Optional numeric vector of 2 containing start 
#' and end time points. Only one of times or timeRanges can be non-null
#' @param times Optional numeric vector of time points. Only one of times or
#' timeRanges can be non-null
#' @export
#' @return An Epoch object
Epoch <- function(data, electrodes = NULL, timeRanges = NULL, times = NULL) {
    if (!is.null(times) && !is.null(timeRanges)) {
        stop("Only one of times or timeRanges can be non-null")
    }
    if (!is.null(timeRanges) && length(timeRanges) != 2) {
        stop("timeRanges must be a numeric vector of length 2")
    }
    if (!is.null(times) && length(times) != ncol(data)) {
        stop("Length of times must be equal to number of columns in data")
    }
    if (!is.null(electrodes) && nrow(data) != length(electrodes)) {
        stop("Length of electrodes must be equal to number of rows in data")
    }

    # set default time points if not provided
    if (is.null(times)) {
        if (is.null(timeRanges)) {
            ## check if data has colnames as time points
            if (!is.null(colnames(data))) {
                times <- tryToNum(colnames(data))
            }
        } else {
            times <- seq(timeRanges[1], timeRanges[2], length.out = nrow(data))
        }
    } else {
        times <- as.numeric(times)
    }


    # Set default electrode names if not provided
    if (is.null(electrodes)) {
        electrodes <- if (!is.null(rownames(data))) {
            rownames(data)
        } else {
            paste0("E", seq_len(nrow(data)))
        }
    }

    rownames(data) <- electrodes
    colnames(data) <- NULL

    # Create new Epoch object
    .Epoch(
        data = data,
        times = times
    )
}

###############################
## getter and setter
###############################

#' Epoch Methods
#'
#'
#' @description
#' `$electrodes`: Get or set electrode names
#' `$times`: Get or set time points
#' `$timeRange`: Get time range if time points are defined
#' `$data`: Get or set data matrix
#'
#' @param x Epoch object
#' @param name a value name, must be one of 'electrodes', 'times',
#' 'timeRange', 'data'
#' @param value Value to set
#' @rdname Epoch-method
#' @export
setMethod("$", "Epoch", function(x, name) {
    if (!name %in% names(x)) {
        stop(glue("Invalid field name: {name}, must be one of {paste0(names(x), collapse = ', ')}"))
    }
    switch(name,
        electrodes = rownames(x@data),
        times = x@times,
        timeRange = if (!is.null(x@times)) range(x@times) else NULL,
        data = x@data,
        stop(glue("Unexpected field name {name}"))
    )
})


#' @rdname Epoch-method
setMethod("$<-", "Epoch", function(x, name, value) {
    if (!name %in% names(x)) {
        stop(glue("Invalid field name: {name}, must be one of {paste0(names(x), collapse = ', ')}"))
    }
    if (name == "electrodes") {
        rownames(x@data) <- value
    }

    if (name == "times") {
        x@times <- value
    }

    if (name == "timeRange") {
        if (length(value) != 2) {
            stop("timeRange must be a numeric vector of length 2")
        }
        x@times <- seq(value[1], value[2], length.out = nrow(x@data))
    }

    if (name == "data") {
        colNms <- rownames(value)
        rowNms <- colnames(value)

        rownames(value) <- x$electrodes
        colnames(value) <- x$times

        x <- Epoch(value, electrodes = colNms, times = rowNms)
    }

    invisible(x)
})



#' @description `[`: Subset an Epoch object using matrix indexing syntax
#'
#' @param i Row (electrode) indices
#' @param j Column (time) indices
#' @rdname Epoch-method
#' @export
setMethod("[", "Epoch", function(x, i, j) {
    if (!missing(i)){
        i <- checkIndex(i, x$electrodes)
    }
    
    new_data <- x@data[i, j, drop = FALSE]

    if (missing(j)) {
        newTimes <- x@times
    } else {
        newTimes <- x@times[j]
    }

    Epoch(
        data = new_data,
        times = newTimes
    )
})



#' @description `nrow`, `ncol`, `colnames`, `rownames`, `names`: Getting the data properties,
#'  similar to base R functions.
#' @rdname Epoch-method
#' @return nrow: Number of rows in the data
#' @export
setMethod("nrow", "Epoch", function(x) {
    nrow(x@data)
})

#' @rdname Epoch-method
#' @return ncol: Number of columns in the data
#' @export
setMethod("ncol", "Epoch", function(x) {
    ncol(x@data)
})

#' @rdname Epoch-method
#' @return colnames: electrode names of the data
#' @export
setMethod("colnames", "Epoch", function(x) {
    x$times
})

#' @rdname Epoch-method
#' @export
setMethod("colnames<-", "Epoch", function(x, value) {
    x$times <- value
    x
})

#' @rdname Epoch-method
#' @return rownames: time points of the data
#' @export
setMethod("rownames", "Epoch", function(x) {
    x$electrodes
})

#' @rdname Epoch-method
#' @export
setMethod("rownames<-", "Epoch", function(x, value) {
    x$electrodes <- value
    x
})


#' @rdname Epoch-method
#' @return names: Return all available properties for an Epoch object
#' @export
setMethod("names", "Epoch", function(x) {
    c("electrodes", "times", "data", "timeRange")
})

#' @rdname Epoch-method
#' @export
setMethod("names<-", "Epoch", function(x, value) {
    stop("Cannot set names for Epoch object")
    invisible(x)
})


###############################
## Data Conversion Methods
###############################
#' @rdname Epoch-method
#' @export
setMethod("as.matrix", "Epoch", function(x) {
    dt <- x@data
    colnames(dt) <- x$times
    dt
})

#' @rdname Epoch-method
#' @param ... additional arguments
#' @export
setMethod("as.data.frame", "Epoch", function(x, ...) {
    as.data.frame(as.matrix(x), ...)
})




###############################
## other Methods
###############################
#' @description `truncateTime`: Truncating time range
#'
#' @param from Numeric value specifying start of new time range
#' @param to Numeric value specifying end of new time range
#' @return truncateTime: Truncated object
#' @rdname Epoch-method
#' @export
setGeneric("truncateTime", function(x, from, to) standardGeneric("truncateTime"))

#' @rdname Epoch-method
#' @export
setMethod("truncateTime", "Epoch", function(x, from, to) {
    if (is.null(x$times)) {
        if (!isWholeNumber(from) || !isWholeNumber(to)) {
            stop("Time points is not defined for this Epoch object, from and to must be whole numbers")
        }
        indices <- seq(from, to)
        newTimes <- NULL
    } else {
        # current time points
        times <- x$times
        # Find indices within new time range
        indices <- which(times >= from & times <= to)
        newTimes <- times[indices]
    }

    newData <- x$data[, indices, drop = FALSE]

    # Create new Epoch object with truncated data
    Epoch(
        data = newData,
        times = newTimes
    )
})


#' @param object Epoch object
#' @rdname Epoch-method
#' @export
setMethod("show", "Epoch", function(object) {
    dt <- object$data
    pprint(dt, rowdots = 4, coldots = 4, digits = 3)
    cat("\n")
    if (!is.null(object$times)) {
        timeRange <- range(object$times)
        cat(glue("Time range: {timeRange[1]} to {timeRange[2]}"))
        cat("\n")
    }
    cat("Use $ to access its methods. see \"?`Epoch-method`\"\n")
    invisible(object)
})
