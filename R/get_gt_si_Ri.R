#' Compute generation time values
#'
#' Compute the generation time values for the transmission pairs.
#'
#' @param data The data element of `simulate_groups`.
#' @return The modified data frame with the added gt column.
#' @keywords internal
#' @export
get_gt <- function(data) {
  data$gt <- data$date_infection - data$date_infection[match(data$source, data$id)]

}


#' Compute serial interval
#'
#' Compute the serial interval values for the transmission pairs.
#'
#' @param data The data element of `simulate_groups`. Please, name `date_onset` for the date of symptom onset column.
#' @return The modified data frame with the added si column.
#' @keywords internal
#' @export
get_si <- function(data) {
  data$si <- data$date_onset - data$date_onset[match(data$source, data$id)]

}


#' Returns the number of cases each ID generated in each group
#'
#' This function takes in a data frame with columns for source ID and group, and returns the data frame with the number of cases each ID generated in each group.
#'
#' @param data A data frame with columns for source ID and group.
#'
#' @return The data frame with the number of cases each ID generated in each group.
#' @export
#' @examples
#' data <- data.frame(id = c(1, 2, 3, 4, 5),
#'                   source = c(NA, 1, 2, 2, 3))
#' get_Ri(data)

get_Ri <- function(data){
  
  #count the number of times someone appears as a source 
  Ri_mat <- as.matrix(table(data$source))
  
  Ri_df <- data.frame(id = row.names(Ri_mat),
                      Ri = Ri_mat,
                      row.names = NULL)
  
  #Merge back to data
  #data <- merge(data, Ri_df, by = "id", all.x = TRUE)
  data$Ri <- Ri_df$Ri[match(data$id, Ri_df$id)]  
  # Where Ri is NA replace with 0
  data["Ri"][is.na(data["Ri"])] <- 0
  
  return(data)
  
}

