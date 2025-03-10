isWholeNumber <- function(x) {
  return(x %% 1 == 0)
}


# Shifts to the right all strings of a list with a number of blanks
shift <- \(strL, nBlanks = 0) {
  pre <- paste(rep(" ", nBlanks), collapse = "")
  lapply(strL, \(x) sprintf("%s%s", pre, x)) |> unlist()
}

#' Check and filter display index
#' 
#' @param display Numeric or character. Display index
#' @param elecNames Character. All electrode names
checkDisplayIndex <- function(display, elecNames){
  if(is(display, "numeric")){  
    displayTot <- seq_along(elecNames)
    diffDisplayTot<- setdiff(display, displayTot)
    displayFiltered <- display[!display%in%diffDisplayTot]
    displayid <- displayFiltered
  }else{
    diffDisplayTot <- setdiff(display, elecNames)
    displayFiltered <- display[!display%in%diffDisplayTot]
    displayid <- which(elecNames%in%displayFiltered)
  }
    if(length(diffDisplayTot)){
      listDisplayMissing <- paste(diffDisplayTot,collapse=", ")
      displayExist <- paste(displayFiltered,collapse=", ")
      warning(
        glue("Display values {listDisplayMissing} are out of electrode range. I will keep the valid values {displayExist}.")
        )
    } 
    displayid
}