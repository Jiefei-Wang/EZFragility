#' Pt01 seizure 1 (-10:10s) around seizure onset
#' 
#' This data corresponds to the first seizure of patient PT01 from 
#' the Fragility Data Set.
#' The data contains only the good channels. 
#' It has been notch filtered and common average referenced in RAVE. 
#' The full data (pt01Epoch) has been epoched -10:10s around the seizure onset.
#' The acquisition frequency is 1000 Hz
#' EcoG recording gathered in collaboration with the National Institue of Health
#'
#' @docType data
#'
#' @usage 
#' ## EEG data
#' data(pt01Epoch)
#' 
#' @format 
#' pt01Epoch: A Matrix with 20001 rows (time points) and 84 columns (electrodes)
#' 
#' pt01Frag: A fragility object result of example 3 for 
#' calc_adj_frag function help with and 84 columns (electrodes)
#'  \[-1:2\]s around the seizure onset
#'
#' @keywords datasets
#'
#' @source Fragility Multi-Center Retrospective Study
#' (\href{https://openneuro.org/datasets/ds003029/versions/1.0.0}{OpenNeuro})
#' 
#' @aliases pt01Epoch pt01Frag
"pt01Epoch"

