#' Check for available distribution data
#'
#' This function will compare the species or genera in your species list with the
#' entries of the FLORKART-database and returns the matching entries.
#' @param x   A \code{data.frame} with one column containing the species names
#'  in the format 'genus_epithet'.
#' @param level A character string describing wether the check should be conducted
#'  on the "species" or "genus" level.
#' @details If input data consists only of genus names the function will call
#' \code{level = "genus"}.
#' @examples
#' \dontrun{
#' library("distdrawr")
#' specieslist <- data.frame(c("Bellis_perennis","Taraxacum_officinale",
#'                             "Capsella_bursa-pastoris"))
#' check_species(specieslist, level="species")
#' }
#' @return  Returns a data frame containing the matching taxon number and species
#' names of x in the FloraWeb database. Mismatches are saved in the attributes of
#' the data frame. The output can be modified and used in \code{\link{get.dist}}
#' with \code{input = TRUE}.
#' @references Datenbank FLORKART der Floristischen Kartierung Deutschlands,
#' Stand 2013, Bundesamt fuer Naturschutz (BfN) und Netzwerk Phytodiversitaet
#' Deutschland (NetPhyD): \url{http://www.floraweb.de}


check_species <- function(x,level="species"){
    spec <- speciesidref
    if (is.data.frame(x)==FALSE){
      stop("input data needs to be a data.frame")
    }
    if (level!="species"&level!="genus"){
      stop("improper value for level")
    }
    if (is.logical(x[1,])==T|is.numeric(x[1,]==T)){
      stop("input data needs to be of type strings or factors")
    }
    if(length(grep(pattern="_",as.character(x[,1]),ignore.case=T))==0){
      level="genus"
    }
    if ((length(unlist(strsplit(as.character(x[1:nrow(x),1]),split="_",fixed=T)))==
         length(x[,1])*2)==FALSE & length(grep(pattern="_",as.character(x[,1]),ignore.case=T))>0){
      stop("All taxon names must include genus and epithet separated by only one \"_\".
           Hint: Sometimes species like Capsella_bursa-pastoris are written as
           Capsella_bursa_pastoris. This pattern can cause this error.")
    }

    if(length(grep(pattern="_",as.character(x[,1]),ignore.case=T))==0){
      split <- t(matrix(unlist(strsplit(as.character(x[1:nrow(x),]),
                                        split="_",fixed=T))))
    }else{
      split <- matrix(unlist(strsplit(as.character(x[1:nrow(x),]),
                                      split="_",fixed=T)),nrow=2)
    }

    matches <- cbind(NAMNR=NA,species=NA)
    attr(matches,"mismatches") <- NA
    if (level=="species"){
      for (i in 1:ncol(split)){
        matches <- rbind(matches,spec[grep(split[1,i],spec[,2])[which(
                                      grep(split[1,i],spec[,2])%in%
                                      grep(split[2,i],spec[,2]))],])
       if(i==1){ matches <- matches[-1,]}

        if (length(grep(split[1,i],spec[,2])[which(
          grep(split[1,i],spec[,2])%in%
          grep(split[2,i],spec[,2]))]) == 0){
          attr(matches,"mismatches") <- c(attr(matches,"mismatches"),
                                              paste(split[1,i],split[2,i],sep="_"))
          }
      }

    }


    if (level=="genus"){

      for (i in 1:ncol(split)){
        matches <- rbind(matches,spec[grep(split[1,i],spec[,2]),])
        if(i==1){ matches <- matches[-1,]}
        if (length(grep(split[1,i],spec[,2])) == 0){
          attr(matches,"mismatches") <- c(attr(matches,"mismatches"),
                                          paste(split[1,i],split[2,i],sep="_"))
        }
      }
    }

    row.names(matches) <- NULL
    matches$NAMNR <- as.integer(matches$NAMNR)


    if (length(attributes(matches)$mismatches)>0){
    message(length(attributes(matches)$mismatches)," ", "input entries without matching data. For more information use attributes()$mismatches")
    }
    return(matches)
  }


