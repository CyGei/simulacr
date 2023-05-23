#' Simulate outbreaks with transmission trees
#'
#' This function implements an individual-based outbreak simulator which
#' generates transmission trees. A Poisson branching process is used to generate
#' new cases in time, using the distribution of the reproduction number (R) and
#' of the duration of the infectious period (or generation time) to determine 
#' rates of infection. A density-dependence term is used so that individual 
#' infectiousness decreaseswith the proportion of susceptible individuals in the 
#' population (see details). For each simulated case, the simulator generates an 
#' individual reproduction number, date of infection, symptom onset, and reporting 
#' using user-specified distributions. Distributions can be provided either as
#' `distcrete` objects, as functions computing the probability mass functions
#' (PMF), or as vectors of numbers taken as the PMF for values on 0, 1, ...,
#' `length(input) - 1`
#'
#' @param duration the number of days to run the simulation for
#'
#' @param population_size the number of susceptible hosts to use in the
#'   simulation; defaults to 100
#'
#' @param R_values a vector of values to be used as reproduction number; values
#'   will be drawn at random from this vector to determine the expected numbero
#'   of secondary cases for each new case
#'
#' @param dist_incubation the distribution of the incubation period, i.e. the
#'   time interval between infection and onset of symtpoms
#'
#' @param dist_infectious_period the distribution of the infectious period,
#'   i.e. the time interval between the moment a case starts showing symptoms
#'   (onset) and the moment they infect new secondary cases. If specified, leave
#'   `dist_generation_time` to its default (`NULL`).
#'   
#' @param dist_generation_time the distribution of the generation time,
#'   i.e. the time interval between the moment a primary case is infected
#'    and the moment they infect new secondary cases. If specified, leave
#'   `dist_infectious_period` to its default (`NULL`).
#'   
#' @param dist_reporting (optional) the distribution of the reporting delay,
#'   i.e. the time interval between symptom onset and the date at which the case
#'   is notified; if `NULL` (default) reporting dates will not be simulated
#'
#' @details
#' Individual infectiousness is determined as \eqn{R x w(t - t_{onset})} where
#' \eqn{R} is the individual reproduction number, \eqn{w} is the PMF of the
#' duration of the infectious period, \eqn{t} is current time, and
#' \eqn{t_{onset}} is the date of symptom onset for the considered individual.
#' The number of new cases at time *t* is then taken from a Poisson distribution
#' with a rate of infection \eqn{\lambda_t n_s / n} where \eqn{lambda_t} is the
#' sum of all individual infectiousness at time *t*, \eqn{n_s} is the number of
#' susceptible individuals in the population, and *n* is the total population
#' size.
#' 
#' @export
#' 
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#'
#' @examples
#' 
#' ## make toy distributions (showing the 3 possible input types)
#' incubation <- c(0, 1, 1, 1, 1) # numbers = unscaled PMF
#' infectious_period <- make_disc_gamma(5, 1) # distcrete object
#' reporting <- function(x) dpois(x, 5) # PMF function
#' set.seed(1)
#' x <- simulate_outbreak(R = runif(100, 1, 3), # random values on [1;3]
#'                        dist_incubation= incubation,
#'                        dist_infectious_period = infectious_period,
#'                        dist_reporting = reporting)$data
#' dim(x)
#' head(x)
#' tail(x)
#' if (require(epicontacts)) {
#'   plot(x)
#' }

simulate_outbreak <- function(duration = 100, # duration of the simulation
                              population_size = 100,
                              R_values, # average secondary cases
                              dist_incubation, # infection -> onset
                              dist_infectious_period = NULL, # onset -> new infections
                              dist_generation_time = NULL, # infection -> new infections
                              dist_reporting = NULL  # onset -> reporting
                              ) {

  ## General strategy to handle the distribution inputs

  ## We want to be a bit flexible as to how inputs are provided for the
  ## distributions. Inputs should always correspond to a probability mass
  ## function (pmf), from which we derive random number generators if
  ## needed. The function `make_pmf` will accept `distcrete` objects, functions
  ## computing the pmf, or a vector of numbers taken as the relative
  ## probabilities on 0:(n-1) (rescaled to sum to 1 if needed). To create random
  ## number generators, we evaluate the pmf on 0:1000, and pass this on to
  ## make_number_generator.
  
  ## determine if we need to simulate reporting dates
  simulate_reporting <- !is.null(dist_reporting)
  
  ## determine if we need to simulate with dist_infectious_period or dist_generation_time
  simulate_infectious_period <- !is.null(dist_infectious_period)
  simulate_generation_time <- !is.null(dist_generation_time)
  
  # check that either dist_infectious_period or dist_generation_time has been specified
  if ((!simulate_infectious_period && !simulate_generation_time) ||
      (simulate_infectious_period && simulate_generation_time)) {
    stop("Either dist_infectious_period or dist_generation_time must be specified")
  }
  
  ## make random distribution from inputs, from pmf to random numbers
  x_values <- 0:1000
  pmf_incubation <- make_pmf(dist_incubation)
  r_incubation <- make_number_generator(pmf_incubation(x_values))
  
  # dist_infectious_period or dist_generation_time PMF
  if(simulate_infectious_period){
    pmf_infectious_period <- make_pmf(dist_infectious_period)
  } else {
    pmf_generation_time <- make_pmf(dist_generation_time)
  }
  
  assert_R(R_values)
  r_R <- function(n = 1) sample_(R_values, n, replace = TRUE)


  ## make index case
  index_case <- make_index_case(id = draw_labels(1),
                                date_infection = 0,
                                date_onset = r_incubation(1),
                                R = r_R(1))

  if (simulate_reporting) {
    pmf_reporting <- make_pmf(dist_reporting)
    r_reporting <- make_number_generator(pmf_reporting(x_values))
    index_case$date_report <- index_case$date_onset + r_reporting(1)
  }

  n_susceptibles <- population_size - 1
  
  ## make vectors for the outputs
  out <- index_case
  stats <- data.frame()
#browser()
  for (t in seq_len(duration)) {
    
    ## Individual infectiousness at time 't' is defined as the infectious period 
    ## or generation time pmf multiplied by the individual R; the global force 
    ## of infection is a rate defined as the sum of individual infectiousnesses. 
    ## This rate is weighted by the proportion of susceptible individuals in the 
    ## population to introduce density-dependence, the result being used to draw 
    ## the number of secondary cases from a Poisson distribution
    
    if (n_susceptibles > 0) {
      if(simulate_infectious_period){
        individual_infectiousness <- pmf_infectious_period(t - out$date_onset) * out$R
      } else{
        individual_infectiousness <- pmf_generation_time(t - out$date_infection) * out$R
      }
      #rate_infection <- sum(individual_infectiousness) * (n_susceptibles / population_size)
      # n_new_cases <- stats::rpois(1, lambda = rate_infection)
      # n_new_cases <- min(n_susceptibles, n_new_cases)
      rate_infection <- sum(individual_infectiousness) / population_size
      proba_infection <- 1 - exp(-rate_infection)
      n_new_cases <- round(n_susceptibles * proba_infection) #deterministic
    } else {
      n_new_cases <- 0
    }

    ## handle new cases
    if(n_new_cases > 0){

      ## remove susceptibles
      n_susceptibles <- n_susceptibles - n_new_cases

      
      ## build new cases: id, different dates, etc.

      ## infectors are picked from infectious cases, weighted by their
      ## respective infectiousness
      new_infectors <- sample(out$id,
                              size = n_new_cases,
                              replace = TRUE,
                              prob = individual_infectiousness)
      out$source <- c(out$source,
                      new_infectors)

      ## id
      out$id <- c(out$id,
                  draw_labels(n_new_cases))
      
      ## new dates of new infections are the current day
      new_date_infection <- rep(t, n_new_cases)
      out$date_infection <- c(out$date_infection,
                              new_date_infection)

      ## new dates of onset picked using the incubation period
      new_date_onset <- t + r_incubation(n_new_cases)
      out$date_onset <- c(out$date_onset,
                          new_date_onset)

      if (simulate_reporting) {
        ## new dates of report picked using the delay to reporting
        new_date_report <- new_date_onset + r_reporting(n_new_cases)
        out$date_report <- c(out$date_report,
                             new_date_report)
      }
      
      ## new R values are drawn at random
      new_R <- r_R(n_new_cases)
      out$R <- c(out$R,
                 new_R)
      
    }
    stats_t <- data.frame(
      time = t,
      size = population_size, 
      n_susceptible = n_susceptibles,
      prop_susceptible = n_susceptibles/population_size,
      new_cases = n_new_cases,
      infection_rate = rate_infection,
      infection_prob = proba_infection
    )
    stats <- rbind(stats, stats_t)
    
  }


  ## output will be a data.frame with a special type for e.g. printing and
  ## conversion
  if (!simulate_reporting) {
    out$date_report <- NULL
  }
  out <- data.frame(out, stringsAsFactors = FALSE)
  class(out) <- c("outbreak", class(out))
  attr(out, "has_contacts") <- FALSE
  
  return(list(data = out,
              stats = stats))
}

