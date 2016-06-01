#' Sets options for estim_pars
#'
#' @description Sets options for estimation one by one
#' @family regression functions
#' @param key string with option to set
#' @param val value of the option
#' @details no input: print values to screen
#'
#' one input:
#'
#' * "default": sets options at default values
#' * any other key (see below): print value to screen
#'
#' two inputs:
#'
#' * "filter": 1 - use filter (default); 0 - do not;
#' * "pars_init_method":
#'    0 - get initial estimates from automatized computation (default)
#'    1 - read initial estimates from .mat file (for continuation)
#'    2 - read initial estimates from pars_init file
#' * "pseudodata_pets":
#'    0 - put pseudodata together with data (default)
#'    1 - put it apart (only for multispecies estimation)
#' * "results_output":
#'    0 - prints results to screen (default)
#'    1 - prints results to screen, saves to .mat file
#'    2 - saves data to .mat file and graphs to .png files
#'    (prints results to screen using a customized results file when it exists)
#'  * "method": "nm" - use Nelder-Mead method; "no" - do not estimate;
#'
#'  for other options see corresponding options file of the method (e.g. nmregr_options)
#' @examples estim_options("default")
#' @export
estim_options <- function(key = "inexistent", val = ""){

  switch(key,
         default = {
           assign("filter", 1, envir = .GlobalEnv)
           assign("cov_rules", "1species", envir = .GlobalEnv)
           assign("pars_init_method", 0, envir = .GlobalEnv)
           assign("pseudodata_pets", 0, envir = .GlobalEnv)
           assign("results_output", 0, envir = .GlobalEnv)
           assign("method", "nm", envir = .GlobalEnv)
           nmregr_options("default");
         },
         filter = {
           if(is.numeric(val)){
             assign("filter", val, envir = .GlobalEnv)
           } else {
             if(is.numeric(filter))
               cat("filter = ", filter, "\n", sep = "")
             else
               cat("filter = unknown \n")
             cat("0 - do not use filter \n");
             cat("1 - use filter \n");
           }
         },
         cov_rules = {
           if(nchar(val) != 0){
             assign("method", val, envir = .GlobalEnv)
           } else {
             if(exists("method"))
               cat("cov_rules = ", cov_rules, "\n", sep = "")
             else
               cat("cov_rules = unknown \n")
             cat("'1species' - one species only \n");
             cat("'basic' - basic body-size scaling relationships \n");
           }
         },
         pars_init_method = {
           if(is.numeric(val)){
             assign("pars_init_method", val, envir = .GlobalEnv)
           } else {
             if(exists("pars_init_method"))
               cat("pars_init_method = ", pars_init_method, "\n", sep = "")
             else
               cat("pars_init_method = unknown \n")
             cat("0 - get initial estimates from automatized computation \n");
             cat("1 - read initial estimates from .mat file \n");
             cat("2 - read initial estimates from pars_init file \n");
             cat("3 - prints results using a customized results file \n");
           }
         },
         pseudodata_pets = {
           if(is.numeric(val)){
             assign("pseudodata_pets", val, envir = .GlobalEnv)
           } else {
             if(exists("pseudodata_pets"))
               cat("pseudodata_pets = ", pseudodata_pets, "\n", sep = "")
             else
               cat("pseudodata_pets = unknown \n")
             cat("0 - put pseudodata together with data \n");
             cat("1 - put it apart (for multispecies estimation) \n");
           }
         },
         results_output = {
           if(is.numeric(val)){
             assign("results_output", val, envir = .GlobalEnv)
           } else {
             if(exists("results_output"))
               cat("results_output = ", results_output, "\n", sep = "")
             else
               cat("results_output = unknown \n")
             cat("0 - prints results to screen \n");
             cat("1 - prints results, saves to matlab file and produces html \n");
             cat("2 - read initial estimates from pars_init file \n");
             cat("2 - saves to matlab file and produces html \n");
           }
         },
         method = {
           if(nchar(val) != 0){
             assign("method", val, envir = .GlobalEnv)
           } else {
             if(exists("method"))
               cat("method = ", method, "\n", sep = "")
             else
               cat("method = unknown \n")
             cat("'no' - do not estimate \n");
             cat("'nm' - use Nelder-Mead method \n");
           }
         },
         inexistent = {
           if(is.numeric(filter))
             cat("filter = ", filter, "\n", sep = "")
           else
             cat("filter = unknown \n")
           if(exists("cov_rules"))
             cat("cov_rules = ", cov_rules, sep = "")
           else
             cat("cov_rules = unknown \n")
           if(exists("pars_init_method"))
             cat("pars_init_method = ", pars_init_method, "\n", sep = "")
           else
             cat("pars_init_method = unknown \n")
           if(exists("pseudodata_pets"))
             cat("pseudodata_pets = ", pseudodata_pets, "\n", sep = "")
           else
             cat("pseudodata_pets = unknown \n")
           if(exists("results_output"))
             cat("results_output = ", results_output, "\n", sep = "")
           else
             cat("results_output = unknown \n")
           if(exists("method")) {
             cat("method = ", method, sep = "")
             if(method != "no")
               eval(parse(text = paste(method, "regr_options", sep = "")))
           }
           else
             cat("method = unknown")
         },
         {
           if(exists("method")) {
             if(method != "no")
               if(is.numeric(val))
                 eval(parse(text = paste("info <- ", method, "regr_options(key = '", key, "', val = ", val, ")", sep = "")))
               else
                 eval(parse(text = paste("info <- ", method, "regr_options(key = '", key, "', val = '", val, "')", sep = "")))
           }
           else
             cat("method = unknown \n")
           if(info == 0) {
             cat(paste("key ", key, " is unkown \n", sep = ""))
             cat("method = ", method, "\n", sep = "")
             if(is.numeric(filter))
               cat("filter = ", filter, "\n", sep = "")
             else
               cat("filter = unknown \n")
             if(exists("cov_rules"))
               cat("cov_rules = ", cov_rules, sep = "")
             else
               cat("cov_rules = unknown \n")
             if(exists("pars_init_method"))
               cat("pars_init_method = ", pars_init_method, "\n", sep = "")
             else
               cat("pars_init_method = unknown \n")
             if(exists("pseudodata_pets"))
               cat("pseudodata_pets = ", pseudodata_pets, "\n", sep = "")
             else
               cat("pseudodata_pets = unknown \n")
             if(exists("results_output"))
               cat("results_output = ", results_output, "\n", sep = "")
             else
               cat("results_output = unknown \n")
           }
          }
  )

  }
