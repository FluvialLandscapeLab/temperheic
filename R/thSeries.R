# myObservations = zoo(
#   data.frame(
#     x0 = cos(seq(0, 6*pi, 0.1))*5 + 12,
#     x10 = cos(seq(-1, 6*pi-1, 0.1))*2.5 + 12,
#     x30 = cos(seq(-3, 6*pi-3, 0.1))*1.25 + 12
#     ),
#   order.by = seq(0, 6*pi, 0.1)*365*86400/(2*pi)
# )

#'@export
thObservations = function(empiricalData, period, nmin, specificUnits = specificUnits()) {
  if(!inherits(empiricalData, "zoo")) stop("empiricalSignals must be a zoo object.")
#  if(length(xVals) != ncol(empiricalData)) stop("Please specify one distance for each column on the empiricalData.")
#  if(!identical(names(empiricalData), paste0("x", xVals))) stop("The names of the columns in empirical data must be: ", paste0(paste0("x", xVals), collapse = ", "))
  nSeries = ncol(empiricalData)

  # some vectors are calculated later on.  These vectors will be converted to
  # matricies that are ncol = nSeries, nrow = nSeries.  The following line of
  # code claculates the vector locations of the matrix diagnal -- where well
  # data is compared to itself
  diagnalLocs = (0:(nSeries-1))*nSeries + 1:nSeries

  combos = expand.grid(from = names(empiricalData), to = names(empiricalData), stringsAsFactors = F)
  #combos = combos[combos$from < combos$to,]
  #combos = data.frame(apply(combos, 2, function(.x) paste0("x", .x)), stringsAsFactors = F)

  # calculate phases using phase shifting approach
  dphase = t(
    sapply(
      1:nrow(combos),
      function(row) {
        if(row %in% diagnalLocs) {
          result = list(minimum = 0, objective = 0)
        } else {
          results = optimize(f = laggedMSR, c(-period/8, period*7/8), thSeriesPair = empiricalData[,as.character(combos[row,])], nmin)
        }
      }
    )
    )

  dphase[dphase[,"objective"] < 0,"minimum"] = NA
  dphase = unlist(dphase[,"minimum"])
  deltaPhaseRadians = dphase*2*pi/period

  ampRatio = sapply(
    1:nrow(combos),
    function(rowIndex) {
      if(is.na(dphase[rowIndex])) {
        result = NA
      } else {
        result = laggedData(dphase[rowIndex], empiricalData[,as.character(combos[rowIndex,])])
        result = laggedModel(result)
        result = coefficients(result)[2]
      }
      return(result)
    }
  )

  eta = -log(ampRatio)/deltaPhaseRadians
  eta[diagnalLocs] = NaN

  derivedVals = derivedArray(ampRatio, deltaPhaseRadians, eta, names(empiricalData))

  myZoo = structure(
    empiricalData,
#    phaseRadians = phaseRadians,
#    amplitude = amplitude,
    derivedVals = derivedVals
    # empiricalVals = empiricalVals,
  )
  class(myZoo) = c("thObservations", class(myZoo))

  return(myZoo)
}

#' @export
thSeries = function(sg, xVals, tVals, specificUnits, ayeScale = 1, beeScale = 1) {
  if(!identical(attr(sg, "specificUnits"), specificUnits)) stop ("The units of the signal are not the same as those specified in specificUnits.")

  # make variables from required attributes.
  assign("aq", attr(sg, "aquifer"))
  assign("bd", attr(sg, "boundary"))
  assign("hy", attr(sg, "hydro"))

  # variable in eqn 8, Luce et al 2013
  a = (1 / (2 * hy$diffusivity_effective)) * (sqrt(((sqrt((hy$advectiveThermVel^4) + (4 * 2 * pi * bd$frequency * hy$diffusivity_effective)^2) + hy$advectiveThermVel^2)/2)) - hy$advectiveThermVel)
  b = (1 / (2 * hy$diffusivity_effective)) * sqrt(((sqrt((hy$advectiveThermVel^4) + (4 * 2 * pi * bd$frequency * hy$diffusivity_effective)^2) - hy$advectiveThermVel^2)/2))

  a = a * ayeScale
  b = b * beeScale

  # create the individual time series from the vector of xVals ("well locations") specified by the user.
  seriesList =
    lapply(
      xVals,
      function(x) bd$mean + bd$amplitude * exp(-a * x) * cos((2 * pi / bd$period) * tVals - (b * x) - (2 * pi * bd$phase / bd$period))
    )
  names(seriesList) = paste0("x", xVals)

  # create a vector of amplitude values -- one for each xVal ("well location")
  amplitude =
    sapply(
      xVals,
      function(x)  bd$amplitude * exp(-a *x)
    )
  names(amplitude) = names(seriesList)

  # create a vector of phases, in radians, for each xVal ("well location")
  phaseRadians =
    sapply(
      xVals,
      function(x) (b * x + (2 * pi * bd$phase / bd$period))
    )
  names(phaseRadians) = names(seriesList)

  ##### CREATE GRIDS OF VALUES FOR PERMUTATIVE PAIRS OF xVals ("WELL LOCATIONS")

  #create factorial combinations of all wells as a data.frame
  factorialDists = expand.grid(from = names(seriesList), to = names(seriesList))

  ampRatio = amplitude[factorialDists$to] / amplitude[factorialDists$from]
  deltaPhaseRadians = phaseRadians[factorialDists$to] - phaseRadians[factorialDists$from]
  eta = -log(ampRatio) / deltaPhaseRadians


  names(xVals) = names(seriesList)
  deltaXvals = (xVals[factorialDists$to] - xVals[factorialDists$from])
#   diffusivity_effectiveEmp = (eta * (2*pi / bd$period) * deltaXvals^2) / (log(ampRatio)^2 + deltaPhaseRadians^2)
#   advectiveThermVelEmp = (((2*pi / bd$period) * deltaXvals) / sqrt((log(ampRatio)^2 + deltaPhaseRadians^2))) * ((1 - eta^2) / sqrt(1 + eta^2))

  derivedVals = derivedArray(ampRatio, deltaPhaseRadians, eta, names(seriesList))
#   empiricalVals = array(
#     data = c(diffusivity_effectiveEmp, advectiveThermVelEmp),
#     dim = c(length(xVals), length(xVals), 2),
#     dimnames = list(from = names(seriesList), to = names(seriesList), value = c("diffusivity_effectiveEmp", "advectiveThermVelEmp"))
#   )

  myZoo = structure(
    zoo(do.call(data.frame, args = c(seriesList)), order.by = tVals),
    signal = sg,
    phaseRadians = phaseRadians,
    amplitude = amplitude,
    derivedVals = derivedVals
    # empiricalVals = empiricalVals,
  )
  class(myZoo) = c("thSeries", "thObservations", class(myZoo))

  return(myZoo)
}

# #' @export
# as.xts.thSeries = function(x, POSIXctRange, xRange = c(1, nrow(x)), trim = F, ...) {
#   if(trim) {
#     x = x[xRange[1]:xRange[2],]
#   }
#   numericRange = as.numeric(POSIXctRange)
#   time = x$time[xRange]
#   time2Date = lm(POSIXctRange ~ time)
#   POSIXctIndex = as.POSIXct(round(predict(time2Date, newdata = x)), origin = "1970-01-01 00:00.00 UTC")
#   class(x) = "data.frame"
#   return(xts::as.xts(x, order.by = POSIXctIndex))
# }
