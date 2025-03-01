FragStat <- setClass(
    "FragStat",
    slots = list(
        qmatrix = "matrixOrNULL",
        cmeansoz = "numericOrNULL",
        cmeansozc = "numericOrNULL",
        csdsoz = "numericOrNULL",
        csdsozc = "numericOrNULL"
    )
)



#' Getters and Setters for S4 object
#' 
#' @param x S4 object
#' @param name Slot name
#' @param value Value to set
#' @return S4 object itself or slot value
#' @export
setMethod("$", "FragStat", function(x, name) {
    slot(x, name)
})

#' @rdname cash-FragStat-method
setMethod("$<-", "FragStat", function(x, name, value) {
    slot(x, name) <- value
    invisible(x)
})

setMethod("show", "FragStat", function(object) {
    cat("\nFragStat object (Summary Statistics by Step)\n")
    printSlots(object)
    cat("Use '$attr' to access the data\n")
    invisible(object)
})
