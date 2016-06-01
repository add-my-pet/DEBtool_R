#' Prints results of estimation
#'
#' @description Prints the results of the esimation procedure in the screen, .mat file and makes figures of graphs
#' @family add-my-pet functions
#' @param par data frame with parameter values
#' @param metaPar data frame with metainformation on models
#' @param txtPar data frame with information on parameters
#' @param data data frame with data values
#' @param auxData data frame with auxiliary data values
#' @param metaData data frame with metainformation on the entry
#' @param txtData data frame with infromation on data
#' @param weights data frame with values of weights
#' @examples results_pets(par, metaPar, txtPar, data, auxData, metaData, txtData, weights)
#' @export
results_pets <- function(par, metaPar, txtPar, data, auxData, metaData, txtData, weights){

  petsNumber <- length(pets)

  cov.rulesnm <- paste("cov_rules_", metaPar$covRules, sep = "")

  weightsMRE <- weights;  # define a weights structure with weight 1 for every data point and 0 for the pseudodata
  for(currentPet in pets) {   # makes st only with dependent variables
    if("psd" %in% names(weights[[currentPet]])) {
      psdSets = names(weights[[currentPet]]$psd);
      for(currentPsd in psdSets)
        weightsMRE[[currentPet]]$psd[[currentPsd]] <- 0 * weights[[currentPet]]$psd[[currentPsd]];
    }
  }

  if(petsNumber == 1) {
    list[MRE, RE, info] <- mre_st("predict_pets", par, data, auxData, weightsMRE); # WLS-method
    if(info == 0)
      stop(  "One parameter set did not pass the customized filters in the predict file")
    metaPar[[pets]]$MRE <- MRE; metaPar[[pets]]$RE = RE;
  }

  data2plot <- data;

  if(results_output < 3) {

#    if(results_output == 2) # to avoid saving figures generated prior the current run
#      close all


    univarX <- list();
    for(currentPet in pets) {
      st <- data2plot[[currentPet]];
      nm <- fieldnmnst(st);
      for(currentField in nm) {  # replace univariate data by plot data
        varData <- st[[currentField]];   # scalar, vector or matrix with data in field nm{i}
        k <- dim(as.matrix(varData))[2];
        if(k == 2) {
          auxDataFields <- names(auxData[[currentPet]]);
          dataCode <- tail(currentField);
          univarAuxData <- list();
          for(currentAuxData in auxDataFields) # add to univarAuxData all auxData for the data set that has length > 1
            if(dataCode %in% auxData[[currentPet]][[currentAuxData]] && length(auxData[[currentPet]][[currentAuxData]][[dataCode]]) > 1)
              univarAuxData[end + 1] <- currentAuxData;
          dataVec <- st[[currentField]][,1];
          if(length(univarAuxData) == 0) { # if there is no univariate auxiliary data the axis can have 100 points otherwise it will have the same points as in data
            xAxis <- as.matrix(seq(min(dataVec), max(dataVec), len = 100)); #';
            univarX[[dataCode]] <- "dft"; # "dft": default number of plotting points for the x-axis
          } else {
            xAxis <- dataVec;
            univarX[[dataCode]] <- "usr"; # "usr": user defined number of plotting points for the x-axis
          }
        st[[currentField]] <- xAxis;
        }
      }
    }

    data2plot[[currentPet]] <- st;
  }

  list[prdData] = predict_pets(par, data2plot, auxData);

  petNm <- 1;
  for(currentPet in pets) {
    if(file.exists(paste("custom_results_", currentPet, ".R", sep = "")))
      do.call(paste(C("custom_results_", currentPet), sep = ""), list(par, metaPar, data[[currentPet]], txtData[[currentPet]], auxData[[currentPet]]))
    else {
      st <- data[[currentPet]];
      nm <- fieldnmnst(st);
      counter <- 0;
      for(currentData in nm) {
        varData <- st[[currentData]];   # scalar, vector or matrix with data in field nm{i}
        k <- dim(as.matrix(varData))[2];
        if(k == 2){
          xData = st[[currentData]][,1];
          yData = st[[currentData]][,2];
          plot(xData, yData,
               xlab = paste(txtData[[currentPet]]$label[[currentData]][1], ", ", txtData[[currentPet]]$units[[currentData]][1], sep = ""),
               ylab = paste(txtData[[currentPet]]$label[[currentData]][2], ", ", txtData[[currentPet]]$units[[currentData]][2], sep = ""),
               col = "red")
          xPred = data2plot[[currentPet]][[currentData]][,1];
          yPred = prdData[[currentPet]][[currentData]];
          lines(xPred, yPred, col = "blue");
          title(txtData[[currentPet]]$bibkey[[currentData]]);
        }
      }
    }
    if(results_output < 2) {
      cat(currentPet, "\n"); # print the species name
      cat("COMPLETE = ", metaData[[currentPet]]$COMPLETE ,"\n")
      cat("MRE = ", metaPar[[currentPet]]$MRE, "\n\n")

      cat("\n");
      printprd(data[[currentPet]], txtData[[currentPet]], prdData[[currentPet]], metaPar[[currentPet]]$RE);

      free <- par$free;
      corePar <- par;  corePar$free <- NULL;
      coreTxtPar <- list();  coreTxtPar$units <- txtPar$units; coreTxtPar$label <- txtPar$label;
      parFields <- names(corePar);
      for(currentPar in parFields)
        if(grepl("n.", currentPar) == 1 || grepl("mu.", currentPar) == 1 || grepl("d.", currentPar) == 1) {
          corePar[[currentPar]] <- NULL;
          coreTxtPar$units[[currentPar]] <- NULL;
          coreTxtPar$label[[currentPar]] <- NULL;
          free[[currentPar]] <- NULL;
        }
      corePar <- do.call(paste("cov_rules_", cov_rules, sep = ""), list(corePar, petNm));
      parFreenm <- names(free);
      for(currentPar in parFreenm)
        if(length(free[[currentPar]]) != 1)
          free[[currentPar]] = free[[currentPar]][petNm];
      corePar.free <- free;
      printpar(corePar, coreTxtPar);
      cat("\n")



    }

  }

}

cov_rules_1species <- function(par, i) par
