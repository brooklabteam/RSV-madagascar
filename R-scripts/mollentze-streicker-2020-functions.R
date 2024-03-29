library(rphylopic)
library(ggtree)
library(cowplot)
library(ggnewscale)
library(ggimage)

get_partial_resids <- function(gamFit, terms, seWithMean) {
  predType <- ifelse(seWithMean, 'iterms', 'terms')  # Doesn't have much meaning here, but included for consistency with get_partial_preds
  
  linearTerm <- predict(gamFit, type = predType, terms = terms) %>% 
    rowSums()
  
  partialResids <- residuals(gamFit) + linearTerm     # TODO: unclear if using deviance residuals is the best idea (differs from plot.gam)
  
  partialResids
}
get_partial_preds <- function(gamFit, newdata, terms, seWithMean) {
  predType <- ifelse(seWithMean, 'iterms', 'terms')
  
  predict(gamFit, newdata = newdata, se.fit = T,
          type = predType, terms = terms) %>% 
    lapply(rowSums) %>% 
    as.data.frame() %>% 
    rename(y = fit, se = se.fit) %>% 
    mutate(ylower = y - 1.96*se,
           yupper = y + 1.96*se)
}
get_partial_effects_interaction <- function(gamFit, var1, var2, seWithMean = TRUE, fixedEffect = FALSE) {
  ## Term names: wrap in s():
  if (!is.null(var2)) {
    if (fixedEffect) stop('Non-smooth interactions not implemented')
    
    termnames <- c(paste0('s(', var1, ',', var2, ')'))
    
  } else {
    if (!fixedEffect) {
      termnames <- paste0('s(', var1, ')')
    } else {
      termnames <- var1
    }
  }
  
  
  ## Variables not part of the effect / interaction are kept constant:
  modelData <- gamFit$model
  responseIndex <- attr(modelData, 'terms') %>% attr('response')
  responseName <- colnames(modelData)[responseIndex]
  
  otherData <- modelData %>% 
    select(-one_of(responseName, var1, var2))
  
  numericData <- otherData %>% 
    summarise_if(is.numeric, ~ median(.))
  
  factorData <- otherData %>% 
    summarise_if(is.factor, ~ names(which.max(table(.))))		
  
  stopifnot(all(colnames(otherData) %in% c(colnames(numericData), colnames(factorData))))  # Would indicate unhandled column types
  
  
  ## Calculate partial residuals
  partialDat <- modelData %>% 
    data.frame() %>% 
    select(one_of(var1, var2))
  
  if (length(numericData) > 0) partialDat <- cbind(partialDat, numericData)
  if (length(factorData) > 0) partialDat <- cbind(partialDat, factorData)
  
  
  partialResids <- get_partial_resids(gamFit, termnames, seWithMean)
  partialResids <- cbind(partialDat,
                         Residual = partialResids)
  
  
  ## Predictions
  # - Make a prediction for each level of the interaction var (so for all interactions that occur)
  newData <- modelData %>% 
    data.frame() %>% 
    select(one_of(var1, var2)) %>% 
    unique()
  
  # - All other data get set to their median (or the most common factor level)
  if (length(numericData) > 0) newData <- cbind(newData, numericData)
  if (length(factorData) > 0) newData <- cbind(newData, factorData)
  
  
  # - Make predictions
  newPredictions <- get_partial_preds(gamFit, newData, termnames, seWithMean) %>% 
    mutate(IsSignificant = if_else(ylower <= 0 & yupper >= 0, 'No', 'Yes')) %>%    # Check if CI crosses zero
    cbind(newData)
  
  # Add significance to residuals (for plotting):
  partialResids <- newPredictions %>% 
    select(one_of(var1, var2), IsSignificant) %>% 
    right_join(partialResids)
  
  # Return:
  list(effects = newPredictions, partialResiduals = partialResids)
}
get_partial_effects <- function(fit, var, seWithMean = TRUE) {
  get_partial_effects_interaction(fit, var, NULL, seWithMean)
}
get_partial_effects_binary_single <- function(fit, var, seWithMean = TRUE, fixedEffect = TRUE, removeNegatives = TRUE) {
  plotData <- get_partial_effects_interaction(fit, var1 = var, NULL, seWithMean, fixedEffect)
  
  # Remove negatives
  if (removeNegatives) {
    plotData$effects <- plotData$effects[plotData$effects[[var]] == 1, ]
    plotData$partialResiduals <- plotData$partialResiduals[plotData$partialResiduals[[var]] == 1, ]
  }
  
  # Add a column containing var as a label
  plotData$effects$variable <- var
  plotData$partialResiduals$variable <- var
  
  # Return
  plotData
}
get_partial_effects_binary <- function(fit, vars, seWithMean = TRUE, fixedEffect = TRUE, removeNegatives = TRUE) {
  allData <- lapply(vars, get_partial_effects_binary_single, fit = fit, 
                    seWithMean = seWithMean, 
                    fixedEffect = fixedEffect, 
                    removeNegatives = removeNegatives)
  
  extract_by_name <- function(x, name) x[[name]]
  effects <- lapply(allData, extract_by_name, 'effects')
  partialResiduals <- lapply(allData, extract_by_name, 'partialResiduals')
  
  effects <- do.call(rbind, effects)
  partialResiduals <- do.call(rbind, partialResiduals)
  
  list(effects = effects, partialResiduals = partialResiduals)
}
get_partial_effects_continuous <- function(gamFit, var, resolution = 1, seWithMean = TRUE) {
  ## Term names: wrap in s():
  termnames <- paste0('s(', var, ')')
  
  
  ## Data not part of effect kept constant:
  modelData <- gamFit$model
  responseIndex <- attr(modelData, 'terms') %>% attr('response')
  responseName <- colnames(modelData)[responseIndex]
  
  otherData <- modelData %>% 
    select(-one_of(responseName, var))
  
  numericData <- otherData %>% 
    summarise_if(is.numeric, ~ median(.))
  
  factorData <- otherData %>% 
    summarise_if(is.factor, ~ names(which.max(table(.))))		
  
  stopifnot(all(colnames(otherData) %in% c(colnames(numericData), colnames(factorData))))  # Would indicate unhandled column types
  
  
  ## Calculate partial residuals
  partialDat <- modelData %>% 
    data.frame() %>% 
    select(one_of(var))
  
  if (length(numericData) > 0) partialDat <- cbind(partialDat, numericData)
  if (length(factorData) > 0) partialDat <- cbind(partialDat, factorData)
  
  
  partialResids <- get_partial_resids(gamFit, termnames, seWithMean)
  partialResids <- cbind(partialDat,
                         Residual = partialResids)
  
  
  ## Predictions
  # - Predictions over a smooth range of values spanning the range of var:
  newData <- seq(min(modelData[, var]), max(modelData[, var]), by = resolution) %>% 
    data.frame()
  
  colnames(newData) <- var
  
  # - All other data get set to their median (or the most common factor level)
  if (length(numericData) > 0) newData <- cbind(newData, numericData)
  if (length(factorData) > 0) newData <- cbind(newData, factorData)
  
  # - Make predictions
  newPredictions <- get_partial_preds(gamFit, newData, termnames, seWithMean) %>% 
    mutate(NotSignificant = ylower <= 0 & yupper >= 0,
           IsSignificant = if_else(all(NotSignificant), 'No', 'Yes')) %>%    # Check if CI crosses 0 over entire range
    cbind(newData)
  
  partialResids$IsSignificant <- unique(newPredictions$IsSignificant)
  
  # Return:
  list(effects = newPredictions, partialResiduals = partialResids)
}

