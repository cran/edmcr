#' Create Molecule Metadata
#' 
#' \code{makeA} Produces a list containing all of the meta data necessary
#' to complete a protein molecule.
#' 
#' @details 
#' Creates a list containing the following information about the amino acid structure: \cr
#' anchors \cr
#' atom \cr
#' cliqs \cr
#' cliqs_woh \cr
#' dihedral \cr
#' first_meta \cr
#' info \cr
#' last_meta \cr
#' meta_bonds \cr
#' omission_map \cr
#' X \cr
#' 
#' All of the raw data is available in the package source
#' 
#' @param Go Set to TRUE to create the list
#' 
#' @return
#' \item{A}{A list containing the required meta data to complete a protein molecule}
#' 
#' @examples 
#' A <- makeA()
#' 
#' @export
makeA <- function(Go=TRUE){
  
  filepath <- paste(system.file('extdata', package = 'edmcr'),"/AMatrix/",sep="")
  
  if(Go){
    InsertA <- list(anchors=list(),
                    atom = list(name = c()),
                    bonds = list(),
                    cliqs = list(),
                    cliqs_woh = list(),
                    code = list(),
                    description = list(),
                    dihedral = list(name = c(),
                                    atoms = c(),
                                    last_atom = c()),
                    enlarged_atoms = list(),
                    first_meta = list(),
                    info = list(),
                    last_meta=list(),
                    meta_bonds=list(),
                    name=list(),
                    num_atoms = list(),
                    num_diangles = list(),
                    omission_map=list(),
                    pseudo_names = list(),
                    pseudos = list(),
                    X = list())
    
    A <- list(InsertA,
              InsertA,
              InsertA,
              InsertA,
              InsertA,
              InsertA,
              InsertA,
              InsertA,
              InsertA,
              InsertA,
              InsertA,
              InsertA,
              InsertA,
              InsertA,
              InsertA,
              InsertA,
              InsertA,
              InsertA,
              InsertA,
              InsertA,
              InsertA)
    
    #Atom Names
    for(i in 1:21){
      num <- c(14,32,19,16,15,20,23,11,21,25,26,30,23,27,20,15,18,28,28,22,14)
      for(j in 1:num[i]){
        out <- as.character(read.table(paste(filepath,"atom\\name\\A",i,"atom_name",j,".csv",sep=""),as.is=TRUE))
        if(nchar(out) > 1){
          out2 <- character()
          for(k in seq(1,nchar(out),2)){
            out2 <- paste(out2,substr(out,k,k),sep="")
          }
          out <- out2
        }
        A[[i]]$atom$name <- c(A[[i]]$atom$name,out)
      }
    }
    
    #X
    for(i in 1:21){
      out <- read.csv(paste(filepath,"X\\A",i,"X.csv",sep=""), header=FALSE)
      A[[i]]$X <- out
    }
    
    #last_meta
    for(i in 1:21){
      out <- read.csv(paste(filepath,"last_meta\\A",i,"last_meta.csv",sep=""), header=FALSE)
      A[[i]]$last_meta <- as.numeric(out)
    }
    
    #first_meta
    for(i in 1:21){
      out <- read.csv(paste(filepath,"first_meta\\A",i,"first_meta.csv",sep=""), header=FALSE)
      A[[i]]$first_meta <- as.numeric(out)
    }
    
    #omission_map
    for(i in 1:21){
      out <- read.csv(paste(filepath,"omission_map\\A",i,"omission_map.csv",sep=""), header=FALSE)
      A[[i]]$omission_map <- out
    }
    
    #dihedral
    for(i in 1:21){
      out1 <- read.csv(paste(filepath,"dihedral\\name\\A",i,"dihedral_name.csv",sep=""), header=FALSE,as.is=TRUE)
      out <- matrix("A", nrow=nrow(out1),ncol=1)
      for(j in 1:nrow(out1)){
        out[j,] <- trimws(paste(out1[j,],sep="",collapse=""))
      }
      
      out2 <- read.csv(paste(filepath,"dihedral\\atoms\\A",i,"dihedral_atoms.csv",sep=""), header=FALSE)
      out3 <- read.csv(paste(filepath,"dihedral\\last_atom\\A",i,"dihedral_last_atom.csv",sep=""), header=FALSE)
      
      A[[i]]$dihedral$name <- out
      A[[i]]$dihedral$atoms <- out2
      A[[i]]$dihedral$last_atom <- out3
    }
    
    #cliqs
    num_cliqs <- c(2,5,3,3,3,4,4,1,3,5,5,6,5,3,1,3,4,3,4,4,2)
    for(i in 1:21){
      for(j in 1:num_cliqs[i]){
        out <- read.csv(paste(filepath,"cliqs\\A",i,"cliq",j,".csv",sep=""), header=FALSE)
        A[[i]]$cliqs[[j]] <- out
      }
    }
    
    #cliqs_woh
    num_cliqs <- c(1,5,3,3,2,4,4,1,3,3,3,5,4,3,1,2,2,3,3,2,2)
    for(i in 1:21){
      for(j in 1:num_cliqs[i]){
        out <- read.csv(paste(filepath,"cliqs_woh\\A",i,"cliq",j,".csv",sep=""), header=FALSE)
        A[[i]]$cliqs_woh[[j]] <- out
      }
    }
    
    #info
    for(i in 1:21){
      out <- read.csv(paste(filepath,"info\\A",i,"info.csv",sep=""), header=FALSE)
      A[[i]]$info <- out
    }
    
    #meta_bond
    for(i in 1:21){
      out <- read.csv(paste(filepath,"meta_bonds\\A",i,"meta_bonds.csv",sep=""), header=FALSE)
      A[[i]]$meta_bonds <- out
    }
    
    #anchors
    numel <- c(1,3,1,1,2,2,2,0,1,3,4,5,3,1,0,2,2,1,2,3,1)
    for(i in 1:21){
      if(numel[i] > 0){
        A[[i]]$anchors <- matrix(list(), nrow=numel[i], ncol=2)
        for(j in 1:numel[i]){
          out1 <- as.matrix(read.csv(paste(filepath,"anchors\\A",i,j,"1.csv",sep=""), header=FALSE))
          out2 <- as.matrix(read.csv(paste(filepath,"anchors\\A",i,j,"2.csv",sep=""), header=FALSE))
          A[[i]]$anchors[[j,1]] <- out1
          A[[i]]$anchors[[j,2]] <- out2
        }
      }
    }
    return(A)
    
  }
}
