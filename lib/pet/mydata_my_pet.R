## mydata_my_pet
# Sets referenced data


mydata_my_pet <- function() {
  
  metaData <- list(); data <- list(); units <- list(); label <- list(); bibkey <- list(); comment <- list();
  
  ## set metaData
  
  metaData$phylum     <- "my_pet_phylum"
  metaData$class      <- "my_pet_class"
  metaData$order      <- "my_pet_order"
  metaData$family     <- "my_pet_family"
  metaData$species    <- "my_pet"
  metaData$species_en <- "my_pet_english_name"
  #metaData$T_typical  <- C2K(20)  # K, body temp
  metaData$data_0     <- c("ab", "ap", "am", "Lb", "Lp", "Li", "Wdb", "Wdp", "Wdi", "Ri")   # tags for different types of zero-variate data
  metaData$data_1     <- c("t-L", "L-W")  # tags for different types of uni-variate data
  
  metaData$COMPLETE = 2.5   # using criteria of LikaKear2011
  
  metaData$author   = "FirstName1 LastName1"              # put names as authors as separate strings:  {'FirstName1 LastName2','FirstName2 LastName2'} , with corresponding author in first place 
  metaData$date_subm = c(2015, 04, 20)                    # [year month day], date at which the entry is submitted
  metaData$email    = "myname@myuniv.univ"                # e-mail of corresponding author
  metaData$address  = "affiliation, zipcode, country"     # affiliation, postcode, country of the corresponding author
  
  # uncomment and fill in the following fields when the entry is updated:
  # metaData$author_mod_1  = "FirstName3 LastName3"            # put names as authors as separate strings:  {'author1','author2'} , with corresponding author in first place 
  # metaData$date_mod_1    = c(2015, 04, 20)                   # [year month day], date modified entry is accepted into the collection
  # metaData$email_mod_1   = "myname@myuniv.univ"              # e-mail of corresponding author
  # metaData$address_mod_1 = "affiliation, zipcode, country"   # affiliation, postcode, country of the corresponding author
  
  ## set data
  # zero-variate data;
  # typically depend on scaled functional response f.
  # here assumed to be equal for all real data; the value of f is specified in pars_init_my_pet.
  # add an optional comment structure to give any additional explanations on
  # how the value was chosen, see the last column of the ab data set for an
  # example
  
  # age 0 is at onset of embryo development
  data$ab = 15;      units$ab = 'd';    label$ab = 'age at birth';  bibkey$ab = 'MollCano2010';   comment$ab  = 'mean value taken from several measurements'; 
  
  return(list(data, metaData))
}

# [data, auxData, metaData, txtData, weights] 