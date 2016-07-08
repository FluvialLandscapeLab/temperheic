# myObservations = zoo(
#   data.frame(
#     x0 = cos(seq(0, 6*pi, 0.1))*5 + 12,
#     x10 = cos(seq(-1, 6*pi-1, 0.1))*2.5 + 12,
#     x30 = cos(seq(-3, 6*pi-3, 0.1))*1.25 + 12
#     ),
#   order.by = seq(0, 6*pi, 0.1)*365*86400/(2*pi)
# )

#'@rdname thSeries
#'@export
thObservedSeries = function(empiricalData, xVals, period, aquifer, nmin, optimizeRange = c(-period/8, period*7/8), specificUnits = thUnits()) {
  if(!identical(specificUnits, attr(aquifer, "specificUnits"))) stop("specificUnits of aquifer is not equal to specificUnits of empirical data (passed as 'specifiUnits' argument).")

  if(!inherits(empiricalData, "zoo")) stop("empiricalData must be a zoo object.")
#  if(length(xVals) != ncol(empiricalData)) stop("Please specify one distance for each column on the empiricalData.")
#  if(!identical(names(empiricalData), paste0("x", xVals))) stop("The names of the columns in empirical data must be: ", paste0(paste0("x", xVals), collapse = ", "))
  nSeries = ncol(empiricalData)

  # some vectors are calculated later on.  These vectors will be converted to
  # matricies that are ncol = nSeries, nrow = nSeries.  The following line of
  # code calculates the vector locations of the matrix diagnal -- where well
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
          results = optimize(f = laggedMSR, optimizeRange, thSeriesPair = empiricalData[,as.character(combos[row,])], nmin)
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

  ampRatio = derived2DArray(ampRatio, names(empiricalData))
  deltaPhaseRadians = derived2DArray(deltaPhaseRadians, names(empiricalData))
  eta = derived2DArray(eta, names(empiricalData))
  # derivedVals = derivedArray(ampRatio, deltaPhaseRadians, eta, names(empiricalData))

  # remove anything you don't want as an element of the thObservedSeries
  rm(nSeries, dphase, combos, diagnalLocs)

  newObservedSeries = .temperheic(
    thEnvir = environment(),
    thClass = "thObservedSeries",
    generalUnits = c(
      "T"
    ),
    specificUnits = specificUnits
  )
  class(newObservedSeries) = c("thObservedSeries", "thSeries", "temperheic")
  return(newObservedSeries)
}


#   myZoo = structure(
#     empiricalData,
# #    phaseRadians = phaseRadians,
# #    amplitude = amplitude,
#     derivedVals = derivedVals
#     # empiricalVals = empiricalVals,
#   )
#   class(myZoo) = c("thObservations", class(myZoo))
#
#   return(myZoo)
# }

#' Create thSeries and thObservation Objects.
#'
#' thObservations are \code{\link{zoo}} objects.  Each column contains a
#' temperature record for a different well in the aquifer.  thSeries
#' @export
thSeries = function(signal, xVals, tVals, specificUnits) { #, ayeScale = 1, beeScale = 1) {
  if(!identical(attr(signal, "specificUnits"), specificUnits)) stop ("The units of the signal are not the same as those specified in specificUnits.")

  # variable in eqn 8, Luce et al 2013
  a = (1 / (2 * signal$hydro$diffusivity_effective)) * (sqrt(((sqrt((signal$hydro$advectiveThermVel^4) + (4 * 2 * pi * signal$boundary$frequency * signal$hydro$diffusivity_effective)^2) + signal$hydro$advectiveThermVel^2)/2)) - signal$hydro$advectiveThermVel)
  b = (1 / (2 * signal$hydro$diffusivity_effective)) * sqrt(((sqrt((signal$hydro$advectiveThermVel^4) + (4 * 2 * pi * signal$boundary$frequency * signal$hydro$diffusivity_effective)^2) - signal$hydro$advectiveThermVel^2)/2))

  #a = a * ayeScale
  #b = b * beeScale

  # create the individual time series from the vector of xVals ("well locations") specified by the user.

  timeSeries =
    zoo(
      do.call(
        data.frame,
        args =
          c(
            lapply(
              xVals,
              function(x) signal$boundary$mean + signal$boundary$amplitude * exp(-a * x) * cos((2 * pi / signal$boundary$period) * tVals - (b * x) - (2 * pi * signal$boundary$phase / signal$boundary$period))
            )
          )
        ),
      order.by = tVals
      )
  names(timeSeries) = paste0("x", xVals)

  # create a vector of amplitude values -- one for each xVal ("well location")
  amplitude =
    sapply(
      xVals,
      function(x)  signal$boundary$amplitude * exp(-a *x)
    )
  names(amplitude) = names(timeSeries)

  # create a vector of phases, in radians, for each xVal ("well location")
  phaseRadians =
    sapply(
      xVals,
      function(x) (b * x + (2 * pi * signal$boundary$phase / signal$boundary$period))
    )
  names(phaseRadians) = names(timeSeries)

  ##### CREATE GRIDS OF VALUES FOR PERMUTATIVE PAIRS OF xVals ("WELL LOCATIONS")

  #create factorial combinations of all wells as a data.frame
  factorialDists = expand.grid(from = names(timeSeries), to = names(timeSeries))

  ampRatio = derived2DArray(amplitude[factorialDists$to] / amplitude[factorialDists$from], names(timeSeries))
  deltaPhaseRadians = derived2DArray(phaseRadians[factorialDists$to] - phaseRadians[factorialDists$from], names(timeSeries))
  eta = -log(ampRatio) / deltaPhaseRadians


  names(xVals) = names(timeSeries)
  deltaXvals = derived2DArray(xVals[factorialDists$to] - xVals[factorialDists$from], names(timeSeries))
#   diffusivity_effectiveEmp = (eta * (2*pi / signal$boundary$period) * deltaXvals^2) / (log(ampRatio)^2 + deltaPhaseRadians^2)
#   advectiveThermVelEmp = (((2*pi / signal$boundary$period) * deltaXvals) / sqrt((log(ampRatio)^2 + deltaPhaseRadians^2))) * ((1 - eta^2) / sqrt(1 + eta^2))

#  derivedVals = derivedArray(ampRatio, deltaPhaseRadians, eta, names(seriesList))
#   empiricalVals = array(
#     data = c(diffusivity_effectiveEmp, advectiveThermVelEmp),
#     dim = c(length(xVals), length(xVals), 2),
#     dimnames = list(from = names(seriesList), to = names(seriesList), value = c("diffusivity_effectiveEmp", "advectiveThermVelEmp"))
#   )

#   myZoo = structure(
#     zoo(do.call(data.frame, args = c(seriesList)), order.by = tVals),
#     signal = signal,
#     phaseRadians = phaseRadians,
#     amplitude = amplitude,
#     derivedVals = derivedVals
#     # empiricalVals = empiricalVals,
#   )
#   class(myZoo) = c("thSeries", "thObservations", "temporheic", class(myZoo))

  #get rid of a few variables we don't want in the thSeries object.
  rm(a,b, tVals, factorialDists)

  newSeries = .temperheic(
    thEnvir = environment(),
    thClass = "thSeries",
    generalUnits = c(
      "T"
    ),
    specificUnits = specificUnits
  )
  class(newSeries) = c("thSpecifiedSeries", class(newSeries))
  return(newSeries)
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
