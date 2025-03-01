.Fragility <- setClass(
    "Fragility",
    slots = list(
        ieegts = "matrixOrNULL",
        adj = "arrayOrNULL",
        frag = "matrixOrNULL",
        frag_ranked = "matrixOrNULL",
        R2 = "matrixOrNULL",
        lambdas = "numericOrNULL"
    )
)

Fragility <- function(ieegts, adj, frag, frag_ranked, R2, lambdas) {
    if (!pkgData$debug){
        ieegts <- NULL
        adj <- NULL
    }
    .Fragility(
        ieegts = ieegts,
        adj = adj,
        frag = frag,
        frag_ranked = frag_ranked,
        R2 = R2,
        lambdas = lambdas
    )
}

#' @rdname cash-FragStat-method
setMethod("$", "Fragility", function(x, name) {
    slot(x, name)
})

#' @rdname cash-FragStat-method
setMethod("$<-", "Fragility", function(x, name, value) {
    slot(x, name) <- value
    invisible(x)
})


## Define the print method
setMethod("show", "Fragility", function(object) {
    cat("\nFragility object\n")
    if(pkgData$debug){
        slots <- c("ieegts", "adj", "frag", "frag_ranked", "R2", "lambdas")
    }else{
        slots <- c("frag", "frag_ranked", "R2", "lambdas")
    }
    printSlots(object, slots = slots)
    cat("Use '$attr' to access the data\n")
    invisible(object)
})
