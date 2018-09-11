# nmregr
# Calculates least squares estimates using Nelder Mead's simplex method

##
#function [q, info] = 

nmregr=function(func, p, varargin){
  
  # created 2001/09/07 by Bas Kooijman; modified 2013/10/04
  
  ## Syntax
  # [q, info] = <../nmregr.m *nmregr*> (func, p, varargin)
  
  ## Description
  # Calculates least squares estimates using Nelder Mead's simplex method
  #
  # Input
  #
  # * func: string with name of user-defined function
  #
  #     f = func (p, xyw) with
  #       p: k-vector with parameters; xyw: (n,c)-matrix or (n,c)-dataframe; f: n-vector
  #     [f1, f2, ...] = func (p, xyw1, xyw2, ...) with  p: k-vector  and
  #      xywi: (ni,k)-matrix; fi: ni-vector with model predictions
  #     The dependent variable in the output f; For xyw see below.
  #
  # * p: (k,2)-matrix with
  #
  #     p(:,1) initial guesses for parameter values
  #     p(:,2) binaries with yes or no iteration (optional)
  #
  # * xyzi (read as xyw1, xyw2, .. ): (ni,3) matrix with
  #
  #     xywi(:,1) independent variable i
  #     xywi(:,2) dependent variable i
  #     xywi(:,3) weight coefficients i (optional)
  #     xywi(:,>3) data-pont specific information data (optional)
  #     The number of data matrices xyw1, xyw2, ... is optional but >0
  #     Default for xywi(:,3) is (number of data points in the set i)^-1
  #
  # Output
  #
  # * q: matrix like p, but with least squares estimates
  # * info: 1 if convergence has been successful; 0 otherwise
  
  ## Remarks
  # Calls user-defined function 'func'
  # Set options with <nmregr_options.html *nmregr_options*>
    # See
  # <nrregr.html  *nrregr*> for Newton-Raphson method, 
  # <garegr.html  *garegr*> for genetic algorithm,
  # <nrregr2.html *nrregr2*> for 2 independent variables, and 
  # <nmvcregr.html *nmvcregr*> for standard deviations proportional to the mean.
  # It is usually a good idea to run <nrregr.html *nrregr*> on the result of nmregr.
  
  ## Example of use
  # See <../mydata_regr.m *mydata_regr*>
  
  func=match.fun(func)
  
  i = 1 # initiate data set counter
  info = 1 # initiate info setting
  ci = as.character(i) # character string with value of i
  if (!inherits(varargin, "list")){
    varargin<-list(as.matrix(varargin)) #Change the matrix or dataframe to a 1-list
  }
  
  nxyw = length(varargin)
  
  while (i <= nxyw){ # loop across data sets
    if (i == 1){
      listxyw = c('xyw1') # initiate list xyw
#       listf = c("f1") # initiate list f
    } else{ listxyw = paste(listxyw, paste('xyw', ci,sep=""), sep=',') # append list xyw
#             listf = c(listf, paste('f', ci,sep=""))
#               paste(listf, paste('f', ci,sep=""), sep=',') # append list f
          }
  i = i + 1
  ci = as.character(i) # character string with value of i
  }
  
#   nl = length(listxyw)
#   listxyw = listxyw[1:(nl - 1)] # remove last ','
#   nl = length(listf)
#   listf = listf[1:(nl - 1)] # remove last ','
#   
#   global_txt = strrep(['global ', listxyw], ',', ' ');
#   eval(global_txt); # make data sets global
#   
  library(pracma)
  N = zeros(nxyw, 1) # initiate data counter
  Y = vector(length=0)
  W = vector(length=0)  # initiate observations and weight coefficients

  for (i in c(1:nxyw)) { # loop across data sets
    # assing unnamed arguments to xyw1, xyw2, etc
    eval(parse(text=paste('xyw', i, ' = varargin[[', i,']]', sep="")))
    eval(parse(text=paste('N[', i, ']  = nrow(xyw', i, ')', sep=""))) # number of data points
    eval(parse(text=paste('k= ncol(xyw', i, ')', sep="")))
    eval(parse(text=paste('Y = c(Y,xyw', i, '[,2])', sep=""))) # append dependent variables
    if (k > 2) {
      eval(parse(text=paste('W = c(W, xyw', i, '[,3])', sep=""))) # append weight coefficients
    } else {W = c(W, (ones(N[i],1)/ N[i]))} # append weight coefficients

  }
  
  p <- as.matrix(p)
  q = p # copy input parameter matrix to output
  info = 1 # convergence has been successful
  
  np = nrow(p) # k: number of parameters
  k = ncol(p)
  index = c(1:np)
  if (k > 1) {
    index = index[p[,2]==1] # indices of iterated parameters
  }
  n_par = max(length(index))  # number of parameters that must be iterated
  
  # set options if necessary
  if (!exists('max_step_number')){
    nmregr_options(key='max_step_number', 200 * n_par)
  }
  if (!exists('max_fun_evals')){
    nmregr_options(key='max_fun_evals', 200 * n_par);
  }
  if (!exists('tol_simplex')){
    nmregr_options(key='tol_simplex', 1e-4);
  }
  if (!exists('tol_fun')){
    nmregr_options(key='tol_fun', 1e-4)
  }
  if (!exists('report')){
    nmregr_options(key='report', 1)
  }
  if (is.na(max_step_number)){
    nmregr_options(key='max_step_number', 200 * n_par)
  }
  if (is.na(max_fun_evals)){
    nmregr_options(key='max_fun_evals', 200 * n_par)
  }
  if (is.na(tol_simplex)){
    nmregr_options(key='tol_simplex', 1e-4)
  }
  if (is.na(tol_fun)){
    nmregr_options(key='tol_fun', 1e-4)
  }
  if (is.na(report)){
    nmregr_options(key='report', 1)
  }
  
  # Initialize parameters
  rho = 1
  chi = 2
  psi = 0.5
  sigma = 0.5
  onesn = ones(1, n_par)
  two2np1 = 2:(n_par + 1)
  one2n = 1:n_par
  np1 = n_par + 1
  
  # Set up a simplex near the initial guess.
  v = zeros(n_par, np1)
  fv = zeros(1,np1)
  xin = q[index, 1]    # Place input guess in the simplex
  v[,1] = xin
  eval(parse(text=paste('listf = func(q[,1],', listxyw, ')')))
  if (nxyw == 1) {
    fv[,1] = t(W) %*% (listf[[1]] - Y)^2
  } else  {fv[,1] = t(W) %*% ((unlist(listf)-Y)^2)}

  # Following improvement suggested by L.Pfeffer at Stanford
  usual_delta = 0.05             # 5 percent deltas for non-zero terms
  zero_term_delta = 0.00025      # Even smaller delta for zero elements of q
  # Why choose these values (1.05 times each non-zero parameter) to compute the initial simplex?
  for (j in  c(1:n_par)){
    y = xin
    if (y[j] != 0){
      y[j] = (1 + usual_delta) * y[j]
    }  else {y[j] = zero_term_delta}
    v[,j+1] = y
    q[index,1] = y
    eval(parse(text=paste('listf = func(q[,1],', listxyw, ')')))
    if (nxyw == 1) {
      fv[1, j + 1] = t(W) %*% ((listf[[1]] - Y) ^ 2)
    } else {fv[1, j + 1] = t(W) %*% ((unlist(listf) - Y) ^ 2)}
  }   
  
  # sort so v(1,:) has the lowest function value 
  j=order(fv)
  fv=matrix(fv[1,j], nrow=1, byrow=T)
  v = v[,j]
  
  how = 'initial'
  itercount = 1
  func_evals = n_par + 1
  if (report == 1){
    print(paste('step ', itercount, ' ssq ', min(fv), '-',
           max(fv), ' ', how))
  }
  info = 1; # I propose to remove this line, as info as already been assigned the value 1 before (l.109) and as not been modified since. Else remove it before.
  
  # Main algorithm ---------------------------------------------------------------------------------------------------------------------------------------
  # Iterate until the diameter of the simplex is less than tol_simplex
  #   AND the function values differ from the min by less than tol_fun,
  #   or the max function evaluations are exceeded. (Cannot use OR instead of AND.)

  while (func_evals < max_fun_evals && itercount < max_step_number){
    if ( max(max(abs(v[,two2np1] - v[,onesn]))) <= tol_simplex && max(abs(fv[1] - fv[two2np1])) <= tol_fun ){
      break
    }
  
    how = ''
  
  # Compute the reflection point ----
  # xbar = average of the n (NOT n+1) best points
    xbar = rowSums(v[,one2n])/ n_par
    xr = (1 + rho) * xbar - rho * v[,np1]
    q[index,1] = xr
    eval(parse(text=paste('listf = func(q[,1],', listxyw, ')')))
    if (nxyw == 1){
      if(!inherits(varargin, "list")){listf <- list(listf)}
      fxr = t(W) %*% ((listf[[1]] - Y) ^ 2)
    } else { fxr = t(W) %*% ((unlist(listf) - Y) ^ 2) }

    func_evals = func_evals + 1
  
    if (fxr < fv[,1]){
    # Calculate the expansion point ----
      xe = (1 + rho * chi) * xbar - rho * chi * v[, np1]
      q[index,1] = xe
      eval(parse(text=paste('listf = func(q[,1],', listxyw, ')')))
      if (nxyw == 1){
        if(!inherits(varargin, "list")){listf <- list(listf)}
        fxe = t(W) %*% ((listf[[1]] - Y) ^ 2)
        } else {fxe = t(W) %*% ((unlist(listf) - Y) ^2)}
    
      func_evals = func_evals + 1
      if (fxe < fxr){
        v[,np1] = xe
        fv[,np1] = fxe
        how = 'expand'
      } else {
        v[,np1] = xr
        fv[,np1] = fxr
        how = 'reflect'
      }
    } else {# fv[,1] <= fxr 
        if (fxr < fv[,n_par]){
          v[,np1] = xr
          fv[,np1] = fxr
          how = 'reflect'
        } else{
            # fxr >= fv(:,n_par) 
            # Perform contraction
            if (fxr < fv[,np1]){
              # Perform an outside contraction
              xc = (1 + psi * rho) * xbar - psi * rho * v[,np1]
              q[index,1] = xc
              eval(parse(text=paste('listf = func(q[,1],', listxyw, ')')))
              if (nxyw == 1){
                if(!inherits(varargin, "list")){listf <- list(listf)}
                fxc = t(W) %*% ((listf[[1]] - Y) ^ 2)
                } else{fxc = t(W) %*% ((unlist(listf) - Y) ^ 2)}
              func_evals = func_evals + 1
    
              if (fxc <= fxr){
                v[,np1] = xc
                fv[,np1] = fxc
                how = 'contract outside'
              } else { how = 'shrink' } # perform a shrink
            } else {
                # Perform an inside contraction
                xcc = (1 - psi) * xbar + psi * v[,np1]
                q[index,1] = xcc
                eval(parse(text=paste('listf = func(q[,1],', listxyw, ')')))
                
                if (nxyw == 1){
                  if(!inherits(varargin, "list")){listf <- list(listf)}
                  fxcc = t(W) %*% ((listf[[1]] - Y) ^ 2)
                }  else{fxcc = t(W) %*% ((unlist(listf) - Y) ^ 2)}
  
                func_evals = func_evals + 1
    
                if (fxcc < fv[,np1]){
                  v[,np1] = xcc
                  fv[,np1] = fxcc
                  how = 'contract inside'
                } else {
                  # perform a shrink
                how = 'shrink'
                }
            }
  
            if ( how == 'shrink' ){
              for (j in two2np1){
                v[,j] = v[,1] + sigma * (v[,j] - v[,1])
                q[index,1] = v[,j]
                eval(parse(text=paste('listf = func(q[,1],', listxyw, ')')))
                  if (nxyw == 1){
                    if(!inherits(varargin, "list")){listf <- list(listf)}
                    fv[,j] = t(W) %*% ((listf[[1]] - Y) ^ 2)
                    } else {fv[,j] = t(W) %*% ((unlist(listf) - Y)^2)}
              }
              
              func_evals = func_evals + n_par
            }
        }
    }
  
    j=order(fv)
    fv=matrix(fv[1,j], nrow=1, byrow=T)
    v = v[,j]
    itercount = itercount + 1
    if (report == 1){
      print(paste('step ', itercount, ' ssq ', min(fv), '-', max(fv), ' ', how, '\n'))
    }

  }  # while of main loop
  
  
  q[index,1] = v[,1]
  
  fval = min(fv)
  if (func_evals >= max_fun_evals){
    if (report > 0){
      print(paste('No convergences with ', max_fun_evals, ' function evaluations'))
    }
  info = 0
  } else if (itercount >= max_step_number ) {
      if (report > 0){
      print(paste('No convergences with ', max_step_number, ' steps', sep=''))
      }
    info = 0
  } else {
      if (report > 0){
        print('Successful convergence \n')       
      }
  }
  info = 1

return(list(q, info))
}