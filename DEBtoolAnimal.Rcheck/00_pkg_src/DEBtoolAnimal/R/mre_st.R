#' Computes mean relative error
#'
#' @description Computes relative errors and mean relative error for using data and predictions
#' @family mre functions
#' @param func string with predict file name
#' @param par data frame with parameter values
#' @param data data frame with data values
#' @param auxData data frame with auxiliary data values
#' @param weights data frame with values of weights
#' @examples results_pets(par, metaPar, txtPar, data, auxData, metaData, txtData, weights)
#' @export
mre_st <- function(func, par, data, auxData, weights){

  nm <- fieldnmnst(data);
  list[prdData, prdInfo] <- do.call(func, list(par, data, auxData)); # call predicted values for all of the data
  if(prdInfo == 0)
    return(list(merr = NA, rerr = NA))

  rerr <- matrix(0, length(nm), 2);  # prepare output

  for(currentData in nm) {   # first we remove independent variables from uni-variate data sets
    var <- as.matrix(data[[currentData]]);   # scalar, vector or matrix with data in field nm{i}
    k <- dim(var)[2];
    if (k > 1)
      data[[currentData]] <- var[,2];
  }

  i <- 1;
  for(currentData in nm) {  # next we compute the weighted relative error of each data set
    var    <- data[[currentData]];
    prdVar <- prdData[[currentData]];
    w      <- weights[[currentData]];
    meanval <- abs(mean(var));
    diff <- abs(prdVar - var);

    if(sum(diff) > 0 && meanval > 0)
      if(sum(w) != 0)
        rerr[i,1] <- sum(w * abs(prdVar - var)/ meanval)/ sum(w)
      else
        rerr[i,1] <- sum(abs(prdVar - var)/ meanval)
    else
      rerr[i,1] = 0;

    rerr[i,2] = (sum(w)!=0); # weight 0 if all of the data points in a data set were given wieght zero, meaning that that data set was effectively excluded from the estimation procedure
    i <- i + 1;
    }
  end

  merr <- sum(apply(rerr, 1, prod))/ sum(rerr[,2]);


  return(list(merr, rerr, info = prdInfo))
}
