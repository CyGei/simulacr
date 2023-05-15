#' Compute serial interval values
#'
#' Compute the serial interval for the transmission pairs.
#'
#' @param df The input data frame.
#' @return The modified data frame with the added si column.
#' @keywords internal
#' @export
compute_si <- function(df) {
  source_df <- df[, c("id", "date_onset")]
  names(source_df) <- c("source", "source_onset")
  
  df$si <-
    df$date_onset - source_df$source_onset[match(df$source, source_df$source)]
  return(df)
  
}



#' Compute generation time values
#'
#' Compute the generation time values for the transmission pairs.
#'
#' @param df The input data frame.
#' @return The modified data frame with the added gt column.
#' @keywords internal
#' @export
compute_gt <- function(df) {
  source_df <- df[, c("id", "date_infection")]
  names(source_df) <- c("source", "source_infection")
  df$gt <-
    df$date_infection - source_df$source_infection[match(df$source, source_df$source)]
  return(df)
  
}