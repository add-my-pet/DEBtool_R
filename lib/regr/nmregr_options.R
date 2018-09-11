## nmregr_options
# Sets options for function <nmregr.html *nmregr*>
  
  ##
nmregr_options = function (key='inexistent', val=NA) {
  #  created at 2002/02/10 by Bas Kooijman; modified 2015/01/16, 2015/02/27 Goncalo Marques
  
  ## Syntax
  # <../nmregr_options .m *nmregr_options*> (key, val)
  
  ## Description
  # Sets options for function 'nmregr' one by one
  #
  # Input
  #
  # * no input: print values to screen
  # * one input: 
    #
  #    'default' sets options at default values
  #    other keys (see below) to print value to screen
  #
  # * two inputs
  #
  #    'report': 1 - to report steps to screen; 0 - not to;
  #    'max_step_number': maximum number of steps 
  #    'max_fun_evals': maximum number of function evaluations
  #    'tol_simplex': tolerance for how close the simplex points must be
  #       together to call them the same
  #    'tol_tun': tolerance for how close the loss-function values must be
  #       together to call them the same
  #
  # Output
  #
  # * no output, but globals are set to values or values printed to screen
  
  ## Example of use
  # nmregr_options('default'); nmregr_options('report', 0)
  

  if(key=='default'){
      report <<- 1
      max_step_number <<- 500
      max_fun_evals <<- 2000
      tol_simplex <<- 1e-4
      tol_fun <<- 1e-4
    
  }else if (key=='report'){
      if (!exists('val')){
        if (numel(report) != 0){
          cat('report = ', report)
        } else {print('report = unknown \n')}
    	} else {report <<- val}
    
  } else if (key=='max_step_number'){
      if (!exists('val')){
        if (numel(max_step_number) != 0){
          cat('max_step_number = ', max_step_number)
        } else {print('max_step_number = unknown \n')}
      } else{max_step_number <<- val}
    
  } else if (key== 'max_fun_evals'){
      if (!exists('val')){
        if (numel(max_fun_evals) != 0){
          cat('max_fun_evals = ', max_fun_evals)
        } else {print('max_fun_evals = unkown \n')}
      } else {max_fun_evals <<- val}
    
  } else if (key== 'tol_simplex'){
      if (!exists('val')){
        if (numel(tol_simplex) != 0){
          cat('tol_simplex = ', tol_simplex)
        } else {print('tol_simplex = unknown \n')}
      } else {tol_simplex <<- val}
    
  } else if (key == 'tol_fun'){
      if (!exists('val')){
        if (numel(tol_fun) != 0){
          cat('tol_fun = ', tol_fun)}
        else {print('tol_fun = unknown \n')}
      }else {tol_fun <<- val}
   
  } else {
      if (key!='inexistent'){   
        cat('key ', key, ' is unkown')
      }
      if (numel(report) != 0){
        cat('report = ', report, '\n')
      } else {print('report = unknown \n')}
      
      if (numel(max_step_number) != 0){
        cat('max_step_number = ', max_step_number, '\n')
      } else {print('max_step_number = unkown \n')}

      if (numel(max_fun_evals) != 0){
        cat('max_fun_evals = ', max_fun_evals, '\n')
      } else {print('max_fun_evals = unkown')}

      if (numel(tol_simplex) != 0){
        cat('tol_simplex = ', tol_simplex, '\n')
      } else {print('tol_simplex = unknown \n')}

      if (numel(tol_fun) != 0){
        cat('tol_fun = ', tol_fun)
      } else {print('tol_fun = unknown \n')}

    }

}