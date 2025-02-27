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



## Allow users to get and set the slots via $ operator
setMethod("$", "FragStat", function(x, name) {
    slot(x, name)
})

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
