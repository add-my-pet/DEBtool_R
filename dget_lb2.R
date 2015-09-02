
dget_lb2=function(x, y, pars){
  lb=Re(pars[1])
  xb=pars[2]
  xb3=pars[3]
  g=pars[4]
  k=pars[5]
  # y = x e_H; x = g/(g + e)
  # (x,y): (0,0) -> (xb, xb eHb) 
  
  x3 = x^(1/ 3)
  l = x3/ (xb3/ lb - beta0(x, xb)/ 3/ g)
  dy = l + g - y * (k - x)/ (1 - x) * l/ g/ x
  return(list(dy))
}