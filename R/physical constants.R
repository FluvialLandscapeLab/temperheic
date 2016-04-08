k =  200
dh = 10/1000 # m/m
q = k*dh / 86400 # Darcy velocity m/s based on m/day(s) formulation. Q/n = q; q/86400 = Darcy velocity in m/s
n = 0.34 #0.3 #porosity
v.x = (q / n); # m / sec
kt.s = 8.37 #1.5 #2.8 # w / m * C thermal conductivity (2.2.trattori sand table 1B, Stonestrom & Constantz, 2007)
kt.w =  0.598 #0.58
kt.bulk = (n * kt.w) + ((1 - n) * kt.s) # Langevin et al, 2008
cs = 1078 #850; # J /  m^3 * C heat capacity (2.6e6 trattori sand table 1B, Stonestrom & Constantz)
cw = 4186 # J /  m^3 * C heat capacity of water
p.s = 2650 #3000
p.w = 998 #1000
c.vol = ((n * p.w * cw) + ((1-n) * p.s * cs)) # volumetric heat capacity of the bulk
vt = q * ((n * p.w * cw) / c.vol) # areally averaged rate of heat movement (eqn 5, Luce et al 2013)
