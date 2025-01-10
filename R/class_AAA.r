## This file will take precedence over all the other class files
setClassUnion("matrixOrNULL", c("matrix", "NULL"))
setClassUnion("arrayOrNULL", c("array", "NULL"))
setClassUnion("numericOrNULL", c("numeric", "NULL"))

printSlotValue <- function(object, slotName, k = 3) {
    val <- methods::slot(object, slotName)
    if (is.null(val)) {
        msg <- glue::glue("{slotName}: NULL")
    } else if (is.matrix(val) || is.array(val)) {
        truncated <- if (length(as.vector(val)) > k) "..." else ""
        msg <- glue::glue(
                "{slotName} ({paste(dim(val), collapse=' x ')}): ",
                "{paste(head(as.vector(val), k), collapse=', ')}{truncated}"
            )
    } else if (is.numeric(val)) {
        if(length(val)>k){
            msg <- glue::glue(
                    "{slotName} ({length(val)}): ",
                    "{paste(head(val, k), collapse=', ')}..."
                )
        }else{
            msg <- 
                glue::glue(
                    "{slotName}: ",
                    "{paste(val, collapse=', ')}"
                )
        }
    } else {
        msg <- glue::glue("{slotName}: {val}")
    }
    cat(msg)
    cat("\n")
}
