dget_l_V1= function(vH, l, pars){
  # vH: -, scaled maturity
  # l: -, scaled length
  # lref: -,scaled length when acceleration starts (can be lb or ls)
  # dl: d/dvH l during exponential growth
  # called from get_lp, get_lj, get_ls
  k=pars[1]
  lT=pars[2]
  g=pars[3]
  f=pars[4]
  lref=pars[5]
  r   = g * (f - lT - lref)/ lref/ (g + f)    # specific growth rate
  dl  = r * l/ 3                              # d/dt l
  dvH = f * l^3 * (1/ lref - r/ g) - k * vH   # d/dt vH
  dl  = dl/ dvH                               # dl/ dvH
  return(list(dl))
}