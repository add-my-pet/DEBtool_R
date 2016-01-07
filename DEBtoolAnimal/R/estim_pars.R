#' Estimates parameters
#'
#' @description Runs the entire estimation procedures: gets the parameters, gets the data, initiates the estimation procedure and sends the results for handling
#' @family add-my-pet functions
#' @examples estim_pars()
#' @export
estim_pars <- function(){

  list[data, auxData, metaData, txtData, weights] <- mydata_pets();

  petsNumber <- length(pets)

  if(petsNumber == 1) {
    pars.initnm <- paste("pars_init_", pets, sep = "")
    resultsnm <- paste("results_", pets, ".mat", sep = "")
  } else {
    pars.initnm <- "pars_init_group"
    resultsnm <- "results_group.mat"
  }

  # set comments
  switch(toString(pars_init_method),
         "0" = {
           if(petsNumber == 1)
             list[par, metaPar, txtPar, flag] <- get_pars(data.pet1, auxData.pet1, metaData.pet1)
           else
             stop("    For multispecies estimation get_pars cannot be used (pars_init_method cannot be 0)")
         },
         "1" = {
           temp <- readMat(resultsnm)
           par <- temp$par; rm(temp)
           list[par2, metaPar, txtPar] <- do.call(pars.initnm, eval(parse(text = paste("metaData$", pets[1], sep = ""))))
           if(length(names(par$free)) != length(names(par2$free)))
             stop("The number of parameters in pars.free in the pars_init and in the .mat file are not the same.")
           par$free <- par2$free
         },
         "2" = {
           list[par, metaPar, txtPar] <- do.call(pars.initnm, eval(parse(text = paste("metaData$", pets[1], sep = ""))))
         }
  )

  if(petsNumber > 1)
    cov.rulesnm <- paste("cov_rules_", metaPar$covRules, sep = "")
  else
    cov.rulesnm <- "cov_rules_1species"

  # check parameter set if you are using a filter
  if(filter) {
    filternm <- paste("filter_", metaPar.model, sep = "")
    pass <- TRUE
    for(i in petsNumber){
      list[passSpec, flag] <- do.call(pars.initnm, eval(parse(text = paste("metaData$", pets[1], sep = ""))))
      if(!passSec) {
        cat("The seed parameter set is not realistic. \n")
        print_filterflag(flag)
      }
      pass <- pass * passSpec
    }
    if(!pass)
      stop("The seed parameter set is not realistic.")
  }


  return(data)

}


