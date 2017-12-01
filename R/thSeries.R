
#'@rdname thSeries
#'@export
thObservedSeries = function(empiricalData,
                            xVals,
                            aquifer,
                            hydro,
                            period,
                            headGrad,
                            nmin,
                            freq = (2*pi)/period,
                            optimizeRange = c(-1/8, 7/8),
                            specificUnits = thUnits(),
                            laggedLinearFit = T) {

  if((optimizeRange[2] - optimizeRange[1]) != 1) stop("max optimize range - min optimize range must = 1.0")
  if(!inherits(empiricalData, "zoo")) stop("EmpiricalData must be a zoo object")
  # Check to be sure that empiricalData column names are the same as the names of xVals.
  if(!identical(sort(names(empiricalData)), sort(names(xVals)))) stop("The names of the columns in empirical data must match the names of xvals.")
  # Ensure that xVals are in ascending order.
  xVals = xVals[order(xVals)]
  if(!xVals[1] == 0) stop("xVals must include zero to designate input signal.")
  # Reorder empiricalData columns to match xVals.
  empiricalData = empiricalData[,names(xVals)]

  if(!identical(specificUnits, attr(aquifer, "specificUnits"))) stop("specificUnits of aquifer is not equal to specificUnits of empirical data (passed as 'specifiUnits' argument).")

  if(!inherits(empiricalData, "zoo")) stop("empiricalData must be a zoo object.")
  if(length(xVals) != ncol(empiricalData)) stop("Please specify one distance for each column on the empiricalData.")

  # some vectors are calculated later on.  These vectors will be converted to
  # matricies that are ncol = nSeries, nrow = nSeries.  The following line of
  # code calculates the vector locations of the matrix diagnal -- where well
  # data is compared to itself
  nSeries = ncol(empiricalData)
  diagnalLocs = (0:(nSeries-1))*nSeries + 1:nSeries

  if(laggedLinearFit) {
    results <- lagLinFit(empiricalData, period, optimizeRange, nmin)
    relativePhase = NA
    amplitude = NA
  } else {
    results <- fitCosine(empiricalData, period, optimizeRange, nmin)
    relativePhase = attr(results, "phases")
    amplitude = attr(results, "amplitudes")
  }
  ## ampRatio, deltaPhaseRadians, diagnalLocs, empiricalData, xVals, freq
  eta = -log(results$ampRatio)/results$deltaPhaseRadians  # See Luce et al. 2013
  eta[diagnalLocs] = NaN

  factorialDists = expand.grid(from = names(empiricalData), to = names(empiricalData))
  names(xVals) = names(empiricalData)
  deltaXvals = derived2DArray(xVals[factorialDists$to] - xVals[factorialDists$from], names(empiricalData))
  ampRatio = derived2DArray(results$ampRatio, names(empiricalData))
  deltaPhaseRadians = derived2DArray(results$deltaPhaseRadians, names(empiricalData))
  eta = derived2DArray(eta, names(empiricalData))

  # See Luce et al. 2013 equation #60
  advectiveThermVelEmpirical = ((freq * deltaXvals) / (sqrt((log(ampRatio)^2) + deltaPhaseRadians^2))) * ((1 - eta^2) / sqrt(1 + eta^2))
  # See Luce et al. 2013 equation #59
  diffusivity_effective_empirical = (eta * freq * deltaXvals^2) / ((log(ampRatio)^2) + deltaPhaseRadians^2)

  darcyFlux = advectiveThermVelEmpirical * ((aquifer$volHeatCap_bulk) / (aquifer$density_h2o * aquifer$spHeat_h2o) )

  velocity_h2o = darcyFlux / aquifer$porosity

  exponent = log((diffusivity_effective_empirical - hydro$diffusivity_cond)/hydro$dispersivity)/log(advectiveThermVelEmpirical)

  dispersivity = ((diffusivity_effective_empirical * aquifer$volHeatCap_bulk - aquifer$thermCond_bulk)  / (advectiveThermVelEmpirical * aquifer$volHeatCap_bulk))

  hydraulicCond = darcyFlux / headGrad

  rad1A = (advectiveThermVelEmpirical^4 + (4 * 2 * pi * freq * diffusivity_effective_empirical)^2)^0.5
  rad22A = sqrt(2)/(((rad1A + advectiveThermVelEmpirical^2))^0.5 - sqrt(2)*advectiveThermVelEmpirical)
  thermDecayDist = 2*diffusivity_effective_empirical * rad22A

  rad1B = (advectiveThermVelEmpirical^4 + (4 * 2 * pi * freq * diffusivity_effective_empirical)^2)^0.5
  rad22B = sqrt(2)/abs((rad1B - advectiveThermVelEmpirical^2))^0.5
  phaseVel = 2*diffusivity_effective_empirical * rad22B  * 2 * pi * freq

  pecletNumber = (advectiveThermVelEmpirical * thermDecayDist)/diffusivity_effective_empirical

  rm(results, factorialDists, deltaXvals, envir = environment())

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



#' Create thSeries and thObservation Objects.
#'
#' thObservations are \code{\link{zoo}} objects.  Each column contains a
#' temperature record for a different well in the aquifer.  thSeries
#' @export
thSeries = function(signal, xVals, tVals, specificUnits) { #, ayeScale = 1, beeScale = 1) {
  if(!identical(attr(signal, "specificUnits"), specificUnits)) stop ("The units of the signal are not the same as those specified in specificUnits.")

  # variable in eqn 8, Luce et al 2013
  #a = (1 / (2 * signal$hydro$diffusivity_effective)) * (sqrt(((sqrt((signal$hydro$advectiveThermVel^4) + (4 * 2 * pi * signal$boundary$frequency * signal$hydro$diffusivity_effective)^2) + signal$hydro$advectiveThermVel^2)/2)) - signal$hydro$advectiveThermVel)
  #b = (1 / (2 * signal$hydro$diffusivity_effective)) * sqrt(((sqrt((signal$hydro$advectiveThermVel^4) + (4 * 2 * pi * signal$boundary$frequency * signal$hydro$diffusivity_effective)^2) - signal$hydro$advectiveThermVel^2)/2))

  # create the individual time series from the vector of xVals ("well locations") specified by the user.

  timeSeries =
    zoo::zoo(
      do.call(
        data.frame,
        args =
          c(
            lapply(
              xVals,
              function(x) signal$boundary$mean + signal$boundary$amplitude * exp(-x/signal$thermDecayDist) * cos((2 * pi / signal$boundary$period) * ((tVals-tVals[1]) - signal$boundary$phase - (x/signal$phaseVel)))
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
      function(x)  signal$boundary$amplitude * exp(-x/signal$thermDecayDist)
    )
  names(amplitude) = names(timeSeries)

  # create a vector of phases, in radians, for each xVal ("well location")
  phaseRadians =
    sapply(
      xVals,
      function(x) ((x/signal$phaseVel + signal$boundary$phase) * (2 * pi / signal$boundary$period))
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

  #get rid of a few variables we don't want in the thSeries object.
  rm(tVals, factorialDists)

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
