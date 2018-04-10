#check where is the directory you are working
path <- getwd()
setwd(path) # put the directory into the brackets
# In order to run pedigreesim you need to make sure that this folder has lib and PedigreeSim.java file
Check <- any(dir() %in% "lib") & any(dir() %in%"PedigreeSim.jar") # Make sure Check is TRUE and then continue
library(stringi) # The function depends on this package. you install the package if you don't have it

#In order to run Pedigree Sim
run.PedigreeSim <- function(parfile, path.to.PedigreeSim="PedigreeSim.jar") {
  ##To execute (e.g.), type: run.PedigreeSim("myparfile.par")
  ps <- system2(command = "java",
                args = c("-jar",
                         path.to.PedigreeSim,
                         parfile),
                stdout = TRUE,
                stderr = TRUE)
} #run.PedigreeSim()

#Function to use pedigreesim to run and generate the files you need   
SNPhaploblock_fixed <- function(chrnum = 2, len = 100, cen = 50, gaplength=1, ploidy=6,
                                filename, popnum = 200, mtype_list=list(c(1,0), c(2,0), c(1,1), c(3,0)),
                                pairing = 0){
  for(i in 1:chrnum) {
    chrname <- as.character(as.roman(i)) 
    markerpositions <- seq(0,len,gaplength)
    posnames <- sprintf("%03d", markerpositions)
    ContigNme <- paste0("Cm", stri_rand_strings(len/gaplength+1,5,'[0-9]'))
    Conti_pos <- paste0(ContigNme, "_", posnames)
    
    phase_randomizer <- function(x, mtype) {
      c(sample(
        c(rep(
          1, mtype[1]),
          rep(
            0, ploidy - mtype[1]
          ))),
        sample(
          c(rep(
            1, mtype[2]),
            rep(
              0, ploidy - mtype[2]
            )))
      )
    }
    
    SNPset_mat <- matrix(integer(), ncol = 2*ploidy)
    
    for(mtype in mtype_list){
      posnames <- sprintf("%03d", markerpositions)
      markernames <- paste0(Conti_pos, "_", mtype[1], "x", mtype[2])
      SNPset <- t(sapply(markerpositions, phase_randomizer, mtype))
      rownames(SNPset) <- markernames
      SNPset_mat <- rbind(SNPset_mat, SNPset)
    }
    
    
    
    positions <- rep(markerpositions, length(mtype_list))
    linkagemap <- data.frame(marker = rownames(SNPset_mat), chromosome = "I", position = positions, 
                             stringsAsFactors = FALSE)
    
    posorder <- order(linkagemap$position)
    
    sortedSNPset <- SNPset_mat[posorder,]
    sortedlinkagemap <- linkagemap[posorder,]
    
    sortedSNPset <- cbind(rownames(sortedSNPset), as.data.frame(sortedSNPset))
    
    colnames(sortedSNPset) <- c("marker", paste("P1", 1:ploidy, sep="_"), paste("P2", 1:ploidy, sep="_"))
    
    
    chromfile <- data.frame(chromosome = NULL, length = NULL, centromere = NULL, prefPairing = NULL, quadrivalents = NULL)
    
    p <- format(round(pairing, 2), nsmall = 2)
    # chromfile <- data.frame(chromosome = "I", length = len, centromere = cen, prefPairing = p, quadrivalents = "0.00")
    
    ## Updates the .chrom file 
    tempchromfile <- data.frame(chromosome = chrname, length = len, centromere = cen, prefPairing = p, quadrivalents = "0.00")
    chromfile <- rbind(chromfile, tempchromfile)
  }
  
  ## Create data-frame for the chromosome file:
  #chromfile <- data.frame(chromosome = "I", length = len, centromere = cen, prefPairing = "0.00", quadrivalents = "0.00")
  
  write.table(sortedSNPset, paste(filename,c(".gen"),sep=""), quote = FALSE, sep = " ", row.names = FALSE)
  write.table(sortedlinkagemap, paste(filename,c(".map"),sep=""), quote = FALSE, sep = " ", row.names = FALSE)
  write.table(chromfile, paste(filename,c(".chrom"),sep=""), quote = FALSE, sep = " ", row.names = FALSE)
  
  ## Creation of the .par file with usual defaults
  parfilename <- paste(filename,c(".par"),sep="")
  fileConn<-file(parfilename)
  writeLines(c(paste0("PLOIDY = ", ploidy), "MAPFUNCTION = HALDANE", "MISSING = NA", "HAPLOSTRUCT = myhaplostruct", paste("CHROMFILE = ",filename,".chrom",sep=""),
               "POPTYPE = F1", paste("POPSIZE = ", popnum, sep=""), paste("MAPFILE = ", filename, ".map",sep=""),  paste("PEDFILE = ", filename, ".ped",sep=""),
               paste("FOUNDERFILE = ", filename, ".gen", sep=""), paste("OUTPUT = ", filename, "_out", sep=""), "NATURALPAIRING = 0"), fileConn)
  close(fileConn)
  
  message("Thank you. Please check your R working directory.")
}

#Example
SNPhaploblock_fixed(chrnum = 1, #1 chromosome
                    len = 100, # 100 cM length
                    ploidy=4, # ploidy level
                    gaplength = 1, # each maker has 1 cM gap in between
                    mtype_list= list(c(1,0),c(1,1)), # specify the marker type, here we simulate SN and SS marker
                    popnum = 200, #population number
                    cen = 50, # the location of the centromere
                    pairing = 0, #preferential pairing
                    filename= "FixedPo_101SNSS") #the file which will be generated with your specified filename
run.PedigreeSim(paste0("FixedPo_101SNSS",".par")) #make sure the filename is the same
#after you run it, you will have different files within the folder generated
