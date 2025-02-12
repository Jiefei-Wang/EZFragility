Fragility <- setClass(
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

## Allow users to get and set the slots via $ operator
setMethod("$", "Fragility", function(x, name) {
    slot(x, name)
})

setMethod("$<-", "Fragility", function(x, name, value) {
    slot(x, name) <- value
    invisible(x)
})


## Define the print method
setMethod("show", "Fragility", function(object) {
    cat("Fragility object:\n")
    printSlotValue(object, "ieegts")
    printSlotValue(object, "adj")
    printSlotValue(object, "frag")
    printSlotValue(object, "frag_ranked")
    printSlotValue(object, "R2")
    printSlotValue(object, "lambdas")
    cat("Use '$attr' to access the data\n")
    
    invisible(object)
})