#' Gets scaled age at puberty
#'
#' @description Obtains scaled age at puberty.
#' @family scaled get functions
#' @param pars 5-vector with parameters: g, k, l_T, v_H^b, v_H^p
#' @param f optional scalar with functional response (default f = 1)
#' @param lb0 optional scalar with scaled length at birth (default lb0: lb for k = 1)
#' @details Food density is assumed to be constant. Multiply the result with the
#' somatic maintenance rate coefficient to arrive at age at puberty.
#' @return list with scaled age at puberty tau_p = a_p k_M (tp),
#' scaled age at birth tau_b = a_b k_M (tb),
#' scaled length at puberty (lp), scaled length at birth (lb)
#' and indicator equals 1 if successful convergence, 0 otherwise (info)
#' @examples
#' get_tp(c(g = 10, k = 1, lT = 0, vHb = 0.5, vHp = 10), 1)
#' @export
get_tp = function(pars, f = 1, lb0 = as.numeric(pars[4]^(1/3))){

  with(as.list(pars), {
    pars.lb <- c(g = g, k = k, vHb = vHb);

    if(k == 1 && f * (f - lT)^2 > vHp * k) {
      lb <- vHb^(1/3);
      tb <- get_tb(pars.lb, f, lb);
      lp <- vHp^(1/3);
      li <- f - lT;
      rB <- 1/ 3/ (1 + f/g);
      tp <- tb + log((li - lb)/ (li - lp))/ rB;
      info <- 1;
    } else if(f * (f - lT)^2 <= vHp * k) {  # reproduction is not possible
      list[tb, lb, info] <- get_tb (pars.lb, f, lb0);
      tp <- 1e20; # tau_p is never reached
      lp <- 1;    # lp is nerver reached
    } else {  # reproduction is possible
      li <- f - lT;
      irB <- 3 * (1 + f/ g); # k_M/ r_B
      list[lp, lb, info] <- get_lp(pars, f, lb0);
      if(length(lb0) != 2) {  # lb0 = l_b
        tb <- get_tb(pars.lb, f, lb);
        tp <- tb + irB * log((li - lb)/ (li - lp));
      } else {                # lb0 = l and t for a juvenile
        tb = NA;
        l = lb0(1);
        tp = irB * log((li - l)/ (li - lp));
      }
    }

    if(!is.numeric(tp))
      info <- 0
    else if(tp < 0)
      info <- 0

    return(data.frame(tp, tb, lp, lb, info))
  })
}
