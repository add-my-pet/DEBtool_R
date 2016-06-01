#' Calculates least squares estimates using Nelder Mead's simplex method using a filter
#'
#' @description Calculates least squares estimates using Nelder Mead's simplex method using a filter
#' @family add-my-pet functions
#' @param func character string with name of user-defined function
#' @param par list with parameters
#' @param data list with data
#' @param auxData list with auxiliary data
#' @param weights list with weights
#' @param filternm character string with name of user-defined filter function
#' @return list with list with parameters resulting from estimation procedure (par)
#' and indicator 1 if convergence has been successful or 0 otherwise (info)
#' @export
petregr_f <- function(func, par, data, auxData, weights, filternm){

  # option settings
  info <- 1; # initiate info setting

  # prepare variable
  #   st: structure with dependent data values only
  st <- data;
  nm <- fieldnmnst(st); # nst: number of data sets

  for(currentData in nm) {   # makes st only with dependent variables
    auxVar <- st[[currentData]];   # data in field currentData
    k <- dim(as.matrix(auxVar))[2];
    if(k > 1)
      st[[currentData]] <- auxVar[,2];
  }

  # Y: vector with all dependent data
  # W: vector with all weights
  Y <- struct2vector(st, nm);
  W <- struct2vector(weights, nm);

  parnm <- names(par$free);
  np <- length(parnm);
  n_par <- sum(unlist(par$free));
  if(n_par == 0)
    stop(); # no parameters to iterate
  index <- unlist(par$free) == 1;  # indices of free parameters

  free <- par$free; # free is here removed, and after iteration added again
  q <- par;  q$free <- NULL; # copy input parameter matrix into output
  qvec <- unlist(q);
  info <- 1; # convergence has been successful

  # set options if necessary
  if(!exists("max_step_number") || length(max_step_number) == 0)
    nmregr_options("max_step_number", 200 * n_par);
  if(!exists("max_fun_evals") || length(max_fun_evals) == 0)
    nmregr_options("max_fun_evals", 200 * n_par);
  if(!exists("tol_simplex") || length(tol_simplex) == 0)
    nmregr_options("tol_simplex", 1e-4);
  if(!exists("tol_fun") || length(tol_fun) == 0)
    nmregr_options("tol_fun", 1e-4);
  if(!exists("simplex_size") || length(simplex_size) == 0)
    nmregr_options("simplex_size", 0.05);
  if(!exists("report") || length(report) == 0)
    nmregr_options("report", 1);

  # Initialize parameters
  rho <- 1; chi <- 2; psi <- 0.5; sigma <- 0.5;
  onesn <- rep(1, n_par);
  two2np1 <- 2:(n_par + 1);
  one2n <- 1:n_par;
  np1 <- n_par + 1;

  # Set up a simplex near the initial guess.
  xin <- qvec[index];    # Place input guess in the simplex
  v <- as.matrix(xin);
  list[f] <- do.call(func, list(q, data, auxData));
  fv <- W %*% (struct2vector(f, nm) - Y)^2;
  # Following improvement suggested by L.Pfeffer at Stanford
  usual_delta <- simplex_size;     # 5 percent deltas is the default for non-zero terms
  zero_term_delta <- 0.00025;      # Even smaller delta for zero elements of q

  for(j in 1:n_par) {
    y <- xin;
    f_test <- FALSE;
    step_reducer <- 1;             # step_reducer will serve to reduce usual_delta if the parameter set does not pass the filter
    y_test <- y;
    while(!f_test) {
      if(y[j] != 0)
        y_test[j] <- (1 + usual_delta / step_reducer) * y[j]
      else
        y_test[j] <- zero_term_delta / step_reducer;
      qvec[index] <- y_test; q <- relist(qvec, q);
      list[f_test] <- do.call(filternm, list(q));
      if(!f_test) {
        cat("The parameter set for the simplex construction is not realistic. \n");
        step_reducer <- 2 * step_reducer;
      } else {
        list[f, f_test] <- do.call(func, list(q, data, auxData));
        if(!f_test) {
          cat("The parameter set for the simplex construction is not realistic. \n");
          step_reducer <- 2 * step_reducer;
        }
      }
    }
    v <- cbind(v, y_test);   colnames(v) <- NULL;
    fv <- cbind(fv, W %*% (struct2vector(f, nm) - Y)^2);
  }

  # sort so v[,1] has the lowest function value
  v <- v[,order(fv)];
  fv <- t(as.matrix(sort(fv)));

  how <- "initial";
  itercount <- 1;
  func_evals <- n_par + 1;
  if(report == 1)
    cat("step ", itercount, " ssq ", min(fv), "-", max(fv), " ", how, "\n", sep = "");
  end
  info <- 1;

  ## Main algorithm
  ## Iterate until the diameter of the simplex is less than tol_simplex
  ##   AND the function values differ from the min by less than tol_fun,
  ##   or the max function evaluations are exceeded. (Cannot use OR instead of AND.)
  while(func_evals < max_fun_evals && itercount < max_step_number) {
    if(max(abs(v[,two2np1]-v[,onesn])) <= tol_simplex && max(abs(fv[1]-fv[two2np1])) <= tol_fun)
      break
    how = "";

  # Compute the reflection point

    # xbar = average of the n (NOT n+1) best points
    xbar <- rowMeans(v[,one2n]);
    xr <- (1 + rho) * xbar - rho * v[,np1];
    qvec[index] <- xr; q <- relist(qvec, q);
    list[f_test] <- do.call(filternm, list(q));
    if(!f_test)
      fxr <- fv[,np1] + 1
    else {
      list[f, f_test] <- do.call(func, list(q, data, auxData));
      if(!f_test)
        fxr <- fv[,np1] + 1
      else
        fxr <- W %*% (struct2vector(f, nm) - Y)^2;
    }
    func_evals <- func_evals + 1;

    if(fxr < fv[,1]) {
    # Calculate the expansion point
      xe <- (1 + rho * chi) * xbar - rho * chi * v[,np1];
      qvec[index] <- xe; q <- relist(qvec, q);
      list[f_test] <- do.call(filternm, list(q));
      if(!f_test)
        fxe <- fxr + 1
      else {
        list[f, f_test] <- do.call(func, list(q, data, auxData));
        if(!f_test)
          fxe <- fv[,np1] + 1
        else
          fxe = W %*% (struct2vector(f, nm) - Y)^2;
      }
      func_evals <- func_evals + 1;
      if(fxe < fxr) {
        v[,np1] <- xe;
        fv[,np1] <- fxe;
        how = "expand";
      } else {
        v[,np1] <- xr;
        fv[,np1] <- fxr;
        how <- "reflect";
      }
    } else { # fv(:,1) <= fxr
      if(fxr < fv[,n_par]) {
        v[,np1] <- xr;
        fv[,np1] <- fxr;
        how <- "reflect";
      } else { # fxr >= fv(:,n_par)
        # Perform contraction
        if(fxr < fv[,np1]) {
          # Perform an outside contraction
          xc <- (1 + psi * rho) * xbar - psi * rho * v[,np1];
          qvec[index] <- xc; q <- relist(qvec, q);
          list[f_test] <- do.call(filternm, list(q));
          if(!f_test)
            fxc <- fxr + 1
          else {
            list[f, f_test] <- do.call(func, list(q, data, auxData));
            if(!f_test)
              fxc <- fv[,np1] + 1
            else
              fxc <- W %*% (struct2vector(f, nm) - Y)^2;
          }
          func_evals <- func_evals + 1;

          if(fxc <= fxr) {
            v[,np1] <- xc;
            fv[,np1] <- fxc;
            how <- "contract outside";
          } else
            # perform a shrink
            how = "shrink";
        } else {
          # Perform an inside contraction
          xcc <- (1 - psi) * xbar + psi * v[,np1];
          qvec[index] <- xcc; q <- relist(qvec, q);
          list[f_test] <- do.call(filternm, list(q));
          if(!f_test)
            fxcc <- fv[,np1] + 1
          else {
            list[f, f_test] <- do.call(func, list(q, data, auxData));
            if(!f_test)
              fxcc <- fv[,np1] + 1
            else
              fxcc = W %*% (struct2vector(f, nm) - Y)^2;
          }
          func_evals <- func_evals + 1;

          if(fxcc < fv[,np1]) {
            v[,np1] <- xcc;
            fv[,np1] <- fxcc;
            how = "contract inside";
          } else
            # perform a shrink
            how = "shrink";
        }
        if(how == "shrink") {
          for(j in two2np1) {
            f_test <- 0;
            step_reducer <- 1;             # step_reducer will serve to reduce usual_delta if the parameter set does not pass the filter
            while(!f_test) {
              v_test <- v[,1] + sigma / step_reducer * (v[,j] - v[,1]);
              qvec[index] <- v_test; q <- relist(qvec, q);
              list[f_test] <- do.call(filternm, list(q));
              if(!f_test) {
                cat("The parameter set for the simplex shrinking is not realistic. \n");
                step_reducer <- 2 * step_reducer;
              } else {
                list[f, f_test] <- do.call(func, list(q, data, auxData));
                if(!f_test) {
                  cat("The parameter set for the simplex shrinking is not realistic. \n");
                  step_reducer <- 2 * step_reducer;
                }
              }
            }
            v[,j] <- v_test;
            fv[,j] <- W %*% (struct2vector(f, nm) - Y)^2;
          }
          func_evals <- func_evals + n_par;
        }
      }
    }
    v <- v[,order(fv)];
    fv <- t(as.matrix(sort(fv)));
    itercount <- itercount + 1;
    if(report == 1 && itercount %% 10 == 0)
      cat("step ", itercount, " ssq ", min(fv), "-", max(fv), " ", how, "\n", sep = "");
  }   # while

  qvec[index] <- v[,1]; q <- relist(qvec, q);
  q$free <- free; # add substructure free to q,

  fval <- min(fv);
  if(func_evals >= max_fun_evals) {
    cat("No convergences with ", max_fun_evals, " function evaluations\n", sep = "");
    info <- 0;
  } else
    if(itercount >= max_step_number) {
      cat("No convergences with ", max_step_number, " steps\n");
      info <- 0;
    } else {
      cat("Successful convergence \n");
      info = 1;
    }

  return(list(par = q, info))
}


struct2vector <- function(struct, fieldNames) {
  vec <- vector('numeric');
  for(currentField in fieldNames)
    vec <- append(vec, struct[[currentField]]);
  return(vec);
}

