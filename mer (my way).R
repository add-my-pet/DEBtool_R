## mre
# Calculates the mean absolute relative error
#mre('predict_Oncorhynchus_mykiss', pars, Data, tW)
##
mre=function(FUN, p, data, tW) {
  #  created: 2010/05/08 by Bas Kooijman, modified 2011/05/02
  
  ## Syntax
  # [merr rerr] = mre(func, p, varargin)
  
  ## Description
  # Calculates the mean absolute relative error, used in add_my_pet
  #
  # Input
  #
  # * func: character string with name of user-defined function;
  #    see <nrregr.html *nrregr*>
    # * p: (np,nc) matrix with p(:,1) parameter values
  # * xywi: (ni,3) matrix with
  #
  #     xywi(:,1) independent variable
  #     xywi(:,2) dependent variable
  #     xywi(:,3) weight coefficients (optional)
  #     The number of data matrices xyw1, xyw2, ... is optional
  #     The first data matrix is assumed to be zero-variate, 
  #       the others uni-variate, which are first reduced to zero-variate data
  #       if all weight coefficients in a uni-variate data-set are zero,
  #       that relative error gets weight zero
  #
  # Output
  #
  # * merr: scalar with mean absolute relative error
  # * rerr: vector with absolute relative errors
  
  FUN=match.fun(FUN)
  
  # get function values
  estimate=FUN(p[,1], data, tW)
  
  # abs relative error for zero-variate data
  rerr=cbind(abs(estimate[[2]]-tW[,2])/tW[,2], tW[,3])
  merr = sum(rerr[,1]*rerr[,2])/ sum(rerr[,2])
  return(merr)
}