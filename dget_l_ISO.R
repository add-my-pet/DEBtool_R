dget_l_ISO= function(vH, l, pars){
# vH: scalar with scaled maturity
# l: scalar with scaled length
# dl: scalar with d/dvH l
# called from get_lp, get_lj, get_ls
  
  k=pars[1]
  lT=pars[2]
  g=pars[3]
  f=pars[4]
  sM=pars[5]
  
  
  r = g * (f * sM - lT * sM - l)/ l/ (f + g) # specific growth rate
  dl = l * r/ 3                              # d/dt l
  dvH = f * l^2 * (sM - l * r/ g) - k * vH   # d/dt vH
  dl = dl/ dvH                               # dl/ dvH
  return(list(dl))
}
