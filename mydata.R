#-------------------
#------------------- Translated from matlab code of Kooijman, see http://www.bio.vu.nl/thb/deb/deblab/add_my_pet_old/html/Oncorhynchus_mykiss.html
#-------------------

Dir='C:\\Users\\starrlight\\Documents\\GitHub\\Kooijman Rainbow trout'

## mydata_Oncorhynchus_mykiss: rainbow trout, steelhead
# Bas Kooijman 2014/09/26

## References
# * <http://www.bio.vu.nl/thb/deb/deblab/add_my_pet/add_my_pet.pdf *add_my_pet manual*> 
  # * <http://www.bio.vu.nl/thb/research/bib/LikaKear2011.html *LikaKear2011*>
  #   Lika, K., Kearney, M. R., Freitas, V., Veer, H. W. van der, Meer, J. van der, Wijsman, J. W. M., Pecquerie, L. and Kooijman, S. A. L. M.
#   The `covariation method' for estimating the parameters of the standard Dynamic Energy Budget model {I}: philosophy and approach.   
#   J. Sea Res. 66: 270--277.
# * <http://www.bio.vu.nl/thb/research/bib/Kooy2010.html *Kooy2010*>
#   Kooijman, S.A.L.M. 2010 
#   Dynamic Energy Budget theory for metabolic organisation. Cambridge Univ. Press
#   Table 8.1, page 300
# * *YaniHisa2002*
#   T. Yanik, S. A. Hisar and C. Bölükbas (2002)
#   EARLY DEVELOPMENT AND GROWTH OF ARCTIC CHARR (SALVELINUS ALPINUS) AND RAINBOW TROUT (ONCORHYNCHUS MYKISS) AT A LOW WATER TEMPERATURE
#   The Israeli Journal of Aquaculture - Bamidgeh 54(2), 2002, 73-
# * <http://www.fishbase.org/summary/239 *fishbase*>
# * <en.wikipedia.org/wiki/Oncorhynchus_mykiss *Wiki*>

## Facts
# * wiki: many subspecies, e.g. O. m. irideus  (coastal rainbow trout), O. m. gairdneri (Columbia River redband trout)
# * YaniHisa2002: best culturing temp 15-16 C
# * wiki: able to spawn several times, each time separated by months


# Call all the functions

setwd(Dir)
names.files=list.files(Dir)
for (i in c(1:length(names.files))){
  if (!(names.files[i] %in% c("mydata.R", "mer.R","nmrerg_options.R", "mydata 2.R" , "verifs.R", "verifs 2.R"))){
    source(names.files[i])  
  }
}

COMPLETE = 2.5    # judge the level using LikaFrei2011  adjust if you have less or more data
FIT = 9.4         # compute after having obtained the estimates


# set data

# zero-variate data
# real data
ah = 33        #  1 d, age at hatch at f (YaniHisa2002: 30-36 d)
T_ah = 273 + 8.5  # K, temperature for ab
ab = ah + 20   #  2 d, age at birth at f (YaniHisa2002: ah + 19-21 d)
T_ab = 273 + 8.5  # K, temperature for ab
ap = 2.5*365   #  3 d, age at puberty at f (fishbase)
T_ap = 273 + 5  # K, temperature for ap
am = 11*365    #  4 d, life span at f (fishbase)  
T_am = 273 + 5  # K, temperature for am
L0 = 0.45     #       cm, diam egg (wiki)
Lb = 2        #  5 cm, physical length at birth at f (from Salmo trutta)
Lp = 15       #  6 cm, physical length at puberty at f (fishbase)
Li = 120      #  7 cm, ultimate total length at f (fishbase)
Wi = 25400    #  8 g, ultimate wet weight at f (fishbase)
Ri = Wi * 2.5/ 365 #  9 #/d, reproduction rate(wiki: 2000 till 3000 eggs per kg
T_Ri = 273 + 5  # K, temperature for Ri

# pseudo-data from pars_my_pet at T_ref  don't change these data
v = 0.02      # 10 cm/d, energy conductance
kap = 0.80    # 11 -, allocation fraction to soma = growth + somatic maintenance
kap_R = 0.95  # 12 -, reproduction efficiency
p_M = 18      # 13 J/d.cm^3, [p_M] vol-specific somatic maintenance
p_T =  0      # 14 J/d.cm^2, {p_T} surface-specific som maintenance
k_J = 0.002   # 15 1/d, maturity maintenance rate coefficient
kap_G = .8    # 16 -, growth efficiency

# pack data
data = c(ah,  ab,  ap,  am,  Lb,  Lp,  Li,  Wi,  Ri,    # 01:09 real data
        v,  kap,  kap_R,  p_M,  p_T,  k_J,  kap_G)  # 10:16 pseudo data

# weight coefficients for WLS and ML criteria differ
Data = data.frame(data, data, pmin(100,1 / pmax(rep(1e-6, length(data)), data) ^ 2))  # nmregr, nrregr (WLS criterion)
Data[1:9,3] = 10 * Data[1:9,3]  # Data(1:10,3) = 10 * Data(1:10,3)  # give real data more weight
Data[16,3] = 200 * Data[16,3]     # more weight to kap_G

# insert temperature data for rates and times in the first column
Data[c(1 ,2 ,3 ,4 ,9), 1] = c(T_ah , T_ab , T_ap , T_am  ,T_Ri)

txt_data = c(# for presentation of predictions
            '1 ah, d, age at hatch ' ,
            '2 ab, d, age at birth ' ,
            '3 ap, d, age at puberty ' ,
            '4 am, d, life span ' ,
            '5 Lb, cm, physical length at birth ' ,
            '6 Lp, cm, physical length at puberty ' ,
            '7 Li, cm, ultimate physical length ', 
            '8 Wi, g, ultimate wet weight ' ,
            '9 Ri, #/d, maximum reprod rate ' ,
            '10 v, cm/d, energy conductance ' ,
            '11 kap, -, allocation fraction to soma ' ,
            '12 kap_R, -, reproduction efficiency ' ,
            '13 [p_M], J/d.cm^3, vol-spec som maint ', 
            '14 {p_T}, J/d.cm^2, sur-spec som maint ' ,
            '15 k_J, 1/d, maturity maint rate coefficient ' ,
            '16 kap_G, -, growth efficiency')

# uni-variate data 
# t-Ww data from YaniHisa2002 at T = 273 + 8.5
# initial weight 1.54 g

tW = matrix(c(          # time since birth (d), wet weight (g)
      0.202 , 1.471,
      15.609,	2.434,
      31.397,	3.448,
      46.426,	4.716,
      62.217,	6.136,
      77.253,	8.215,
      92.859,	10.294,
      108.271,	11.917,
      123.688	,14.148,
      139.104,	16.227,
      154.725,	20.335),
      ncol=2, byrow=T)
tW = cbind(tW, tW[,2])
tW[,3] = 5e-1/mean(tW[,2])^2  # append weight coefficients for WLS criterion

# conversion coefficients (selected copy-paste from pars_my_pet)

# chemical indices
#       X     V     E     P
n_O = matrix(c(1.00, 1.00, 1.00, 1.00,   # C/C, equals 1 by definition
       1.80, 1.80, 1.80, 1.80,   # H/C  these values show that we consider dry-mass
       0.50, 0.50, 0.50, 0.50,   # O/C
       0.15, 0.15, 0.15, 0.15), ncol=4, byrow=T)  # N/C
#       C     H     O     N
n_M = matrix(c( 1    , 0  ,   0 ,    0  ,   # C/C, equals 0 or 1
        0  ,   2    , 0   ,  3   ,  # H/C
        2  ,   1  ,   2   ,  0  ,   # O/C
        0  ,   0   ,  0   ,  1), ncol=4, byrow=T)    # N/C

# specific densities
#       X     V     E     P
d_O = c(0.2,   0.2,   0.2,   0.2)     # g/cm^3, specific densities for organics

# chemical potentials
#        X     V     E     P
mu_O = c(525,  500 ,  550  , 480) * 1000  # J/mol, chemical potentials for organics

# molecular weights
w_O = t(n_O) %*% c(12,  1 , 16 , 14)   # g/mol, mol-weights for organics

# pack coefficients
dwm = cbind(d_O, w_O, mu_O)  # g/cm^3, g/mol, kJ/mol spec density, mol weight, chem pot

# parameters: initial values at T_ref
T_ref  = 293       # 1 K, temp for which rate pars are given  don't change this value
T_A  = 8000        # 2 K, Arrhenius temp 
f = 1              # 3 -, scaled functional response
z = 9.882          # 4 -, zoom factor  for z = 1: L_m = 1 cm
del_M = 0.1461     # 5 -, shape coefficient
F_m = 6.5          # 6 l/d.cm^2, {F_m} max spec searching rate
kap_X = 0.8        # 7 -, digestion efficiency of food to reserve
v = 0.0486         # 8 cm/d, energy conductance
kap = 0.31         # 9 -, allocation fraction to soma = growth + somatic maintenance
kap_R = 0.95       #10 -, reproduction efficiency
p_M = 18.28        #11 J/d.cm^3, [p_M] vol-specific somatic maintenance
p_T = 0            #12 J/d.cm^2, {p_T} surface-specific som maintenance
k_J = 0.002        #13 1/d, maturity maint rate coefficient
E_G = 5265         #14 J/cm^3, [E_G], spec cost for structure
E_Hh = 3.351e1     #15 J, E_H^h maturity at hatch
E_Hb = 1.510e2     #15 J, E_H^b maturity at birth
E_Hj = 1.107e3     #16 J, E_H^j maturity at metamorphosis
E_Hp = 1.554e5     #17 J, E_H^p maturity at puberty
h_a = 1.115e-6     #18 1/d^2, Weibull aging acceleration
s_G = 1e-4         #19 -, Gompertz stress coefficient
f_tL = 1.97        #20 -, scaled functional response for tL-data
W_0 = 1.5          #21 g, wet weight at zero

# pack parameters and fix T_ref and f and possibly other as well
#   in second column: 0 = fix  1 = release
pars = matrix(c(T_ref, 0,  T_A, 0 , f,    0,  z ,    1 , del_M ,1  ,F_m , 0  , 
        kap_X ,0 , v  , 1 , kap  ,1  ,kap_R ,0 , p_M  , 1 , p_T,  0   ,
        k_J ,  0 , E_G ,0 , E_Hh ,1 , E_Hb ,1 , E_Hj,  1,  E_Hp,  1 , h_a ,  1  ,s_G , 0 , f_tL ,1  ,W_0, 0), ncol=2, byrow=T)

txt_pars = c(    # for presentation of parameter estimates
            'T_ref, K' , 'T_A, K'       ,    'f, -'  ,
            'z, -'  ,    'del_M, -'   ,      '{F_m}, l/d.cm^2'  ,
            'kap_X, -' , 'v, cm/d'   ,       'kap, -'  ,
            'kap_R, -'  ,'[p_M], J/d.cm^3',  '{p_T}, J/d.cm^2'  ,
            'k_J, 1/d'  ,'[E_G], J/cm^3'  ,  'E_Hh, J' , 'E_Hb, J' , 'E_Hj, J' ,
            'E_Hp, J'  , 'h_a, 1/d^2'    ,   's_G, -',  'f_tL, -',  'W_0, g')

# estimate parameters
# 
# nmregr_options('default')  # set options for parameter estimation                                #!!!!!!!!!!!!!!!!!!!!!!!!!!! nmregr_options ???
# nmregr_options('max_step_number',1000)  # set options for parameter estimation
# nmregr_options('max_fun_evals',2e4)    # set options for parameter estimation
# 
# pars = nmregr('predict_Oncorhynchus_mykiss', pars, list(Data,tW))    # WLS estimate parameters using overwrite
# sd = 0 * pars[,1]                                  # initiate standard deviations

#  [cov cor sd] = pregr('predict_Oncorhynchus_mykiss', pars, Data, tW)  # get standard deviation for WLS

# get FIT

# Data[,3] = 0 
# Data[1:9,3] = 1  # give unit weight to real data, zero to pseudo-data
# merr = mre('predict_Oncorhynchus_mykiss', pars, Data, tW)  # WLS-method
# FIT = 10 * (1 - merr)  # get mark for goodness of fit 

# get predictions

tspan = data.frame(0:380)  # times for plotting length data
predict=predict_Oncorhynchus_mykiss(pars[,1], Data, tspan)  # notice use of first column of pars only

# PLOT


plot(tW[,2]~tW[,1], col="red", xlim=c(0,max(tspan[,1])), ylim=c(0, max(predict[[2]])))
points(predict[[2]]~tspan[,1], type="l", col="green")

## Results
# 
# fprintf('FIT = #3.1f \n', FIT)
# fprintf('COMPLETE = #3.1f \n\n', COMPLETE)
# 
# printpar(txt_data, data, Edata, 'data and expectations')  # for zero-variate data
# fprintf('\n') # insert a blank line in screen output
# printpar(txt_pars, pars, sd) 
# 
# close all # remove existing figures, else you get more and more if you retry
# 
# figure # one figure to show results of uni-variate data
# plot(tW(:,1), tW(:,2), '.r', t, EW, 'g', 'linewidth', 2, 'markersize', 20)
# set(gca, 'Fontsize', 15, 'Box', 'on')
# xlabel('time, d')
# ylabel('wet weight, g')
# 
# # options.showCode = false  publish('mydata_Oncorhynchus_mykiss', options) 