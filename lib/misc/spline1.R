spline1 = function(x, knots, Dy1, Dyk) {
  #  created at 2007/03/29 by Bas Kooijman; modified 2009/09/29
  #
  ## Description
  #  First order splines connect knots by straight lines, and is linear outside the knots.
  #  Calculates the ordinates and the first derivatives of a first order spline, given abcissa and knot coordinates. 
  #  The spline is interpolating. 
  #
  ## Input
  #  x: n-vector with abcissa values
  #  knots: (nk,2)-matrix with coordinates of knots; we must have nk > 3
  #         knots(:,1) must be ascending
  #  Dy1: scalar with first derivative at first knot (optional)
  #       empty means: zero
  #  Dyk: scalar with first derivative at last knot (optional)
  #       empty means: zero
  #
  ## Ouput
  #  y: n-vector with spline values (ordinates)
  #  dy: n-vector with derivatives
  #  index: n-vector with indices of first knot-abcissa smaller than x
  #
  ## Remarks
  #  See ispline1 for integration, rspline1 for roots, espline1 for local extremes. 
  #
  ## Example of use
  #  x = (1:10)'; y = 3*(x+.1*rand(10,1)).^2; [Y, dY] = spline1([x,y],k); iY = ispline1(x,k); rspline1(k,5) 
  #  See mydata_smooth for further illustration. 

  x = as.matrix(x); nx = length(x); nk = nrow(knots);
  y = matrix(0,nx,1); dy = y; index = matrix(0,nx,1); # initiate output
  
  if(exists('Dy1') == FALSE){ # make sure that left clamp is specified
    Dy1 = 0;
  }
  if(exists('Dyk') == FALSE){ # make sure that right clamp is specified
    Dyk = 0; 
  }

  ## derivatives right of knot-abcissa
  Dy =c((knots[2:nk,2] - knots[1:nk-1,2]) / 
       (knots[2:nk,1] - knots[1:nk-1,1]), Dyk);
  for(i in 1:nx){ # loop across abcissa values
    j = 1;
    while(x[i] > knots[min(nk,j),1] & j <= nk){
      j = j + 1;
    }
    j = j - 1;
    if(j == 0){      
      y[i] = knots[1,2] - Dy1 * (knots[1,1] - x[i]);
      dy[i] = Dy1;
    }else{
      y[i] = knots[j,2] - Dy[j] * (knots[j,1] - x[i]);
      dy[i] = Dy[j];
    }
    index[i] = j;
  }

return(list(y, dy, index))
}