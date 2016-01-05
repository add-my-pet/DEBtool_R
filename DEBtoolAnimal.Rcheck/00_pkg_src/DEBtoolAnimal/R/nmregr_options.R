#' Sets options for function nmregr
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
#' * "report": 1 - to report steps to screen; 0 - not to;
#' * "max_step_number": maximum number of steps
#' * "max_fun_evals": maximum number of function evaluations
#' * "tol_simplex": tolerance for how close the simplex points must be together to call them the same
#' * "tol_tun": tolerance for how close the loss-function values must be together to call them the same
#' * "simplex_size": fraction added (subtracted if negative) to the free parameters when building the simplex
#' @return 1 if input is valid key, 0 if input is unknown key
#' @examples nmregr_options("default")
#' @export
nmregr_options <- function(key = "inexistent", val = ""){

  switch(key,
         default = {
           assign("report", 1, envir = .GlobalEnv)
           assign("max_step_number", 500, envir = .GlobalEnv)
           assign("max_fun_evals", 2000, envir = .GlobalEnv)
           assign("tol_simplex", 1e-4, envir = .GlobalEnv)
           assign("tol_fun", 1e-4, envir = .GlobalEnv)
           assign("simplex_size", 0.05, envir = .GlobalEnv)
           return(1)
         },
         report = {
           if(is.numeric(val)){
             assign("report", val, envir = .GlobalEnv)
           } else {
             if(exists("report"))
               cat("report = ", report, "\n", sep = "")
             else
               cat("report = unknown \n")
             cat("0 - do not report the regression steps \n");
             cat("1 - report the regression steps \n");
           }
           return(1)
         },
         max_step_number = {
           if(is.numeric(val)){
             assign("max_step_number", val, envir = .GlobalEnv)
           } else {
             if(exists("max_step_number"))
               cat("max_step_number = ", max_step_number, "\n", sep = "")
             else
               cat("max_step_number = unknown \n")
           }
           return(1)
         },
         max_fun_evals = {
           if(is.numeric(val)){
             assign("max_fun_evals", val, envir = .GlobalEnv)
           } else {
             if(exists("max_fun_evals"))
               cat("max_fun_evals = ", max_fun_evals, "\n", sep = "")
             else
               cat("max_fun_evals = unknown \n")
           }
           return(1)
         },
         tol_simplex = {
           if(is.numeric(val)){
             assign("tol_simplex", val, envir = .GlobalEnv)
           } else {
             if(exists("tol_simplex"))
               cat("tol_simplex = ", tol_simplex, "\n", sep = "")
             else
               cat("tol_simplex = unknown \n")
           }
           return(1)
         },
         tol_fun = {
           if(is.numeric(val)){
             assign("tol_fun", val, envir = .GlobalEnv)
           } else {
             if(exists("tol_fun"))
               cat("tol_fun = ", tol_fun, "\n", sep = "")
             else
               cat("tol_fun = unknown \n")
           }
           return(1)
         },
         simplex_size = {
           if(is.numeric(val)){
             assign("simplex_size", val, envir = .GlobalEnv)
           } else {
             if(exists("simplex_size"))
               cat("simplex_size = ", simplex_size, "\n", sep = "")
             else
               cat("simplex_size = unknown \n")
           }
           return(1)
         },
         inexistent = {
           if(exists("report"))
             cat("report = ", pars_init_method, "\n", sep = "")
           else
             cat("report = unknown \n")
           if(exists("max_step_number"))
             cat("max_step_number = ", pars_init_method, "\n", sep = "")
           else
             cat("max_step_number = unknown \n")
           if(exists("max_fun_evals"))
             cat("max_fun_evals = ", pars_init_method, "\n", sep = "")
           else
             cat("max_fun_evals = unknown \n")
           if(exists("tol_simplex"))
             cat("tol_simplex = ", pseudodata_pets, "\n", sep = "")
           else
             cat("tol_simplex = unknown \n")
           if(exists("tol_fun"))
             cat("tol_fun = ", results_output, "\n", sep = "")
           else
             cat("simplex_size = unknown \n")
           if(exists("simplex_size"))
             cat("simplex_size = ", results_output, "\n", sep = "")
           else
             cat("simplex_size = unknown \n")
           return(1)
         },
         {
         cat(paste("key ", key, " is unkown \n", sep = ""))
         if(exists("report"))
           cat("report = ", pars_init_method, "\n", sep = "")
         else
           cat("report = unknown \n")
         if(exists("max_step_number"))
           cat("max_step_number = ", pars_init_method, "\n", sep = "")
         else
           cat("max_step_number = unknown \n")
         if(exists("max_fun_evals"))
           cat("max_fun_evals = ", pars_init_method, "\n", sep = "")
         else
           cat("max_fun_evals = unknown \n")
         if(exists("tol_simplex"))
           cat("tol_simplex = ", pseudodata_pets, "\n", sep = "")
         else
           cat("tol_simplex = unknown \n")
         if(exists("tol_fun"))
           cat("tol_fun = ", results_output, "\n", sep = "")
         else
           cat("simplex_size = unknown \n")
         if(exists("simplex_size"))
           cat("simplex_size = ", results_output, "\n", sep = "")
         else
           cat("simplex_size = unknown \n")
         return(0)
         }
  )

  }
