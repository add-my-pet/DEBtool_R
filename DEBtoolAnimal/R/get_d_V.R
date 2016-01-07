#' Sets chemical parameters and text for units and labels
#'
#' @description Sets chemical parameters and text for units and labels
#' @family add-my-pet auxiliary functions
#' @param phylum string with species phylum
#' @param class string with species class
#' @details Calls get_d_V to set specific density of structure. For a specific density of wet mass of 1 g/cm^3,
#' a specific density of d_E = d_V = 0.1 g/cm^3 means a dry-over-wet weight ratio of 0.1
#' @return list with d_V and info
#' @examples get_d_V("Chordata", "Mammalia")
#' @export
get_d_V <- function(phylum, class){

  info = 1;
  d.V <- switch(phylum,
         "Porifera" = 0.11,
         "Ctenophora" =, "Cnidaria" = 0.01,
         "Gastrotricha" = 0.05,              # Platyzoa
         "Rotifera" = 0.06,
         "Platyhelminthes" =, "Acanthocephala" =, "Chaetognatha" = 0.07,
         "Bryozoa" =, "Entoprocta" =, "Phoronida" =, "Brachiopoda" = 0.07,  # Spiralia
         "Annelida" = 0.16,
         "Sipuncula" = 0.11,
         "Mollusca" =
           switch(class,
                  "Cephalopoda" = 0.21,
                  "Gastropoda" = 0.15,
                  "Bivalvia" = 0.09,
                  0.1
           ),
         "Tardigrada" =, "Priapulida" = 0.07,    # Ecdysozoa
         "Arthropoda" = 0.17,
         "Echinodermata" = 0.09,                 # deuterostomata
         "Hemichordata" = 0.07,
         "Chordata" =
           switch(class,
                  "Mammalia" =, "Reptilia" = 0.3,
                  "Aves" =, "Amphibia" = 0.28,
                  "Chondrichthyes" =, "Actinopterygii" =, "Sarcopterygii" = 0.2,
                  "Myxini" = 0.17,
                  "Cephalaspidomorphi" = 0.125,
                  "Appendicularia" = 0.045,
                  "Thaliacea" = 0.08,
                  0.06, # Ascidiacea
           ),
         "my_pet_phylum" = 0.1,
         {
           cat("warning from get_d_V: taxon could not be identified: d_V = 0.1 g/cm^3\n")
           info = 0
           0.1
         }
  )

  return(list(d.V, info))
}
