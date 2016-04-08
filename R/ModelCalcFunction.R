modelCalc <- function (k, dh, n, kt.bulk, c.vol, dispersivity, p.w, cw, annual = T) {
  if(annual) {
    freq = 86400 * 365
  } else {
    freq = 86400
  }
  q = k*dh / 86400 #unit conversion from m/d to m/s
  # v = q * n
  ## diffusive + dispersive ke standard linear model
  ## ke = (kt.bulk / c.vol)  + (dispersivity) * ((p.w * cw ) / c.vol) * q
  vt = q * ((n * p.w * cw) / c.vol) # areally averaged rate of heat movement (eqn 5, Luce et al 2013)
  w = (2 * pi) / freq # annual frequency
  
  ## Thermal Pelet number eqn 23 Luce et al. 2013
  cond.diffusivity = kt.bulk / c.vol # diffusive spreading
  disp.diffusivity = dispersivity * ((p.w * cw ) / c.vol)*q # dispersive standard linear model
  Dt = cond.diffusivity + disp.diffusivity 
  
  zdisp <- sqrt((2*disp.diffusivity)/w)  ## variable zd diffusive decay depth
  zcond <- sqrt((2*cond.diffusivity)/w)  ## fixed zc condictive decay depth
  zDt = sqrt((2*Dt)/w)
  zRatio = zdisp/rep(zcond, 6)
  
  vdisp = sqrt(2*w*disp.diffusivity)  ## variable phase velocity
  vcond = sqrt(2*w*cond.diffusivity)  ## variable phase velocity
  vDt = sqrt(2*w*Dt)
  vRatio = vdisp/vcond
  vStar = vt/vDt
  
  PeCond = (vt*zcond)/cond.diffusivity # advective/diffusive Peclet
  PeDisp = (vt*zDt)/cond.diffusivity
  r = disp.diffusivity/cond.diffusivity # disp/cond ratio
  return(list(PeCond = PeCond, PeDisp = PeDisp, r = r, zdisp = zdisp, 
              zcond = zcond, zDt = zDt, vdisp = vdisp, vcond = vcond, 
              vDt = vDt, vStar = vStar, zRatio = zRatio, vt = vt))
}

