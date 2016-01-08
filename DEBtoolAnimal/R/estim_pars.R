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
           list[par2, metaPar, txtPar] <- eval(strVec2Exp(c(pars.initnm, "(metaData$", pets[1], ")")))
           if(length(names(par$free)) != length(names(par2$free)))
             stop("The number of parameters in pars.free in the pars_init and in the .mat file are not the same.")
           par$free <- par2$free
         },
         "2" = {
           list[par, metaPar, txtPar] <- eval(strVec2Exp(c(pars.initnm, "(metaData$", pets[1], ")")))
         }
  )

  if(petsNumber > 1)
    cov.rulesnm <- paste("cov_rules_", metaPar$covRules, sep = "")
  else
    cov.rulesnm <- "cov_rules_1species"

  # check parameter set if you are using a filter
  if(filter) {
    filternm <- paste("filter_", metaPar$model, sep = "")
    pass <- TRUE
    for(i in petsNumber){
      list[passSpec, flag] <- do.call(filternm, list(do.call(cov.rulesnm, list(par, i))))
      if(!passSpec) {
        cat("The seed parameter set is not realistic. \n")
        print_filterflag(flag)
      }
      pass <- pass * passSpec
    }
    if(!pass)
      stop("The seed parameter set is not realistic.")
  }

  if(method != "no")
    if(method == "nm")
      if(petsNumber == 1)
        par <- petregr_f("predict_pets", par, data, auxData, weights, filternm)     # WLS estimate parameters using overwrite
      else
        par <- groupregr_f("predict_pets", par, data, auxData, weights, filternm, covRulesnm);  # WLS estimate parameters using overwrite

  # Results
  #results_pets(par, metaPar, txtPar, data, auxData, metaData, txtData, weights);
  predict_my_pet(par, data$my_pet, auxData$my_pet)

  if(filter) {
    warningnm <- paste("warning_", metaPar$model, sep = "")
    if(petsNumber == 1)
      do.call(warningnm, list(par))
    else
      for(i in petsNumber){
        do.call(warningnm, list(par, i))
      }
  }

  return(data)

}


cov_rules_1species <- function(par, i) par


