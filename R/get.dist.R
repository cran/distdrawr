#' Get distribution data for your species list.
#'
#' This function downloads distribution data from the FloraWeb-database.
#'
#'
#' @param x   A \code{data.frame} with one column containing the species names
#'  in the format 'genus_epithet'.
#' @param input  logical. If \code{TRUE x} is treated as the output of \code{check_species},
#' if \code{FALSE x} is treated as the input for \code{check_species} with
#' \code{level = "species"}
#' @param output A character string describing whether the output format should
#' be a \code{"list"} or a \code{"matrix"}.
#' @param matrix.res A character string describing the resolution of the produced
#' \code{matrix}, either \code{"TK25"} or \code{"quarterTK"}.
#' @details This function applies functions from the get.dist-family on your species
#'  list, depending on the value of \code{output}. The default value is \code{list}.
#' @examples
#' \dontrun{
#' library("distdrawr")
#' specieslist <- data.frame(c("Bellis_perennis","Abies_alba",
#'                             "Capsella_bursa-pastoris"))
#' get.dist(specieslist, output = "list")
#' }
#' @return Output depends on the \code{output} value. \code{output = "list"}
#' returns a \code{list} of one \code{data.frame} per taxon, including all
#' information found in the FloraWeb for the particular taxa. \code{output= "matrix"}
#' returns a presence/absence matrix for the TK25 plane survey sheets, when
#' \code{matrix.res = "TK25"} or the quarter TK25 plane survey sheets, when
#' \code{matrix.res = "quarterTK"} including every matching taxon in \code{x}.
#' @references Datenbank FLORKART der Floristischen Kartierung Deutschlands,
#' Stand 2013, Bundesamt fuer Naturschutz (BfN) und Netzwerk Phytodiversitaet
#' Deutschland (NetPhyD): \url{http://www.floraweb.de}
#' @export

get.dist <- function(x, input = FALSE, output = "list", matrix.res = "TK25"){


##checkpoint internet connection
  testCon <- function(url = "http://www.google.com") {

    # test the http capabilities of the current R build
    http <- as.logical(capabilities(what = "http/ftp"))
    if (!http) return(FALSE)

    # test connection by trying to read first line at url
    test <- try(suppressWarnings(readLines(url, n = 1)), silent = TRUE)  # silent errors

    # return FALSE if test is class 'try-error'
    ifelse(inherits(test, "try-error"), FALSE, TRUE)
  }
  con <- testCon()
  if (con==FALSE){
    stop("You have no working internet connection.")
  }

  conFloraWeb <- testCon("http://floraweb.de/pflanzenarten/download_tkq.xsql?suchnr=1")
  if (conFloraWeb==FALSE){
    stop("No respond from the Server of FloraWeb, try again later.")
  }
  #rm(c("con","conFloraWeb","testCon")
##checkpoint data input
if (is.logical(input)==F){
  stop ("Value for input must be logical.")
}
if (input==F){
  if(length(grep(pattern="_",as.character(x[,1]),ignore.case=T))<length(x[,1])){
    stop("All taxon names must be seperated into genus and epithet by \"_\".")
  }
  if ((length(unlist(strsplit(as.character(x[1:nrow(x),1]),split="_",fixed=T)))==
      length(x[,1])*2)==FALSE){
    stop("All taxon names must include genus and epithet seperatet by \"_\".
         Hint: Sometimes species like Capsella_bursa-pastoris are written as
           Capsella_bursa_pastoris. This pattern can cause this error.")
  }
}

##checkpoint "output" value
  if (output!="list"&output!="matrix"){
    stop("improper value for output")
  }
##checkpoint "matrix.res" value
  if (matrix.res!="TK25" & matrix.res!= "quarterTK"){
    stop("Value for matrix.res must be either \"TK25\" or \"quarterTK\" ")
  }
##check_species
  if (input == F){
  matches <- suppressMessages(check_species(x,level="species"))
  }else{
    matches <- x
  }
  if (length(attributes(matches)$mismatches)!=0){
    message(length(attributes(matches)$mismatches)," ", "input entries without matching data:")
    message(attributes(matches)$mismatches)
  }
##Liste
  maps <- as.list(rep(NA,length(matches$NAMNR)))
  if (output=="list"){
    for (i in 1:nrow(matches)){
      url <- url(paste0("http://floraweb.de/pflanzenarten/download_tkq.xsql?suchnr=",
                        matches$NAMNR[i]))
      a <- suppressWarnings(utils::read.table(url, skip=43,sep=","))
      names(a) <- c("ID","NAMNR","TAXONNAME","TK","LAT","LON","FLAECHE",
                    "UNSCHAERFERADIUS","STATUSCODE","STATUS","ZEITRAUM","ZEITRAUM_TEXT")
      maps[[i]] <- data.frame(a)
      names(maps)[i] <- as.character(unique(a$TAXONNAME))
      closeAllConnections()
      message(paste0("  Downloaded  ", paste(i), "/",nrow(matches)))
    }
    return(maps)
  }

##matrix TK25-Plots
  if (output == "matrix" & matrix.res == "TK25"){
    rowsmatrixTK25 <- data.frame(unique(round(rowsmatrix$TK,-1)))
    names(rowsmatrixTK25) <- "TK"
    distmatrix <- matrix(NA,nrow=length(rowsmatrixTK25$TK),ncol=length(matches$NAMNR),
                         dimnames=list(rowsmatrixTK25$TK,matches$species))
    for (i in 1:nrow(matches)){
      url <- url(paste0("http://floraweb.de/pflanzenarten/download_tkq.xsql?suchnr=",
                        matches$NAMNR[i]))
      a <- suppressWarnings(utils::read.table(url, skip=43,sep=",")[,c(3,4)])
      names(a) <-c("TAXONNAME", "TK")
      #if(length(which((a$TK-round(a$TK,-1))==0)) != 0){
     #   a <- a[-which((a$TK-round(a$TK,-1))==0),]
      #}
      distmatrix[,i] <- ifelse(row.names(distmatrix)%in%as.character(round(a$TK,-1)),1,0)
      closeAllConnections()
      message(paste0("  Downloaded  ", paste(i), "/",nrow(matches)))
    }
    return(distmatrix)
  }


  if (output =="matrix" & matrix.res == "quarterTK"){

    distmatrix <- matrix(NA,nrow=length(rowsmatrix$TK),ncol=length(matches$NAMNR),
                         dimnames=list(rowsmatrix$TK,matches$species))
    for (i in 1:nrow(matches)){
      url <- url(paste0("http://floraweb.de/pflanzenarten/download_tkq.xsql?suchnr=",matches$NAMNR[i]))
      a <- suppressWarnings(utils::read.table(url, skip=43,sep=",")[,c(3,4)])
      names(a) <-c("TAXONNAME", "TK")
      distmatrix[,i] <- ifelse(row.names(distmatrix)%in%as.character(a$TK),1,0)
      closeAllConnections()
      message(paste0("  Downloaded  ", paste(i), "/",nrow(matches)))
    }
    return(distmatrix)
  }
} #End of function
