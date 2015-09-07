
fnget_lb2=function(lb, xbxb3gvHbk){
  
  xb=xbxb3gvHbk[1]
  xb3=xbxb3gvHbk[2]
  g=xbxb3gvHbk[3]
  vHb=xbxb3gvHbk[4]
  k=xbxb3gvHbk[5]
  
  # f = y(x_b) - y_b = 0; x = g/ (e + g); x_b = g/ (e_b + g)
  # y(x) = x e_H = x g u_H/ l^3 and y_b = x_b g u_H^b/ l_b^3
  
  tspan=seq(1e-10, xb, length=100)
  yini=c(y=0)
  xy = as.data.frame(ode(yini, tspan, dget_lb2, parms=c(lb, xb, xb3, g, k), method="ode23"))
  f = xy$y[length(xy$y)] - xb * g * vHb/ lb^3
  return(f)
}