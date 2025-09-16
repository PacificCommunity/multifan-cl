# Code processing to replace calls to save_identifier_string()
#
drv <- c("C:/Nick")  # PC
#basedir <- paste(drv,"/MFCL/MFCL_VS2019/upgrd_vs2019",sep="")
basedir <- paste(drv,"/MFCL/MFCL_VS2019/mfcl_devvsn16",sep="")
#basedir <- paste(drv,"/MFCL/MFCL_VS2019/saved/mfcl_devvsn16",sep="")
setwd(basedir)

# Files requiring modification
flnms <- readLines("flnms.txt")

# Save original files to backup directory
bak <- paste(basedir,"/backup",sep="")
dir.create(bak)
for(i in 1:length(flnms)){
  file.copy(from=flnms[i],to=paste(bak,flnms[i],sep="/"))
}

for(ii in 1:length(flnms)){
  flnm <- flnms[ii]
  #
  # All lines containing calls to save_identifier_string()
  a <- readLines(flnm)
  #str <- c("save_identifier_string(\"")
  str <- c("save_identifier_string")
  ptr <- grep(str,a)
  
  # Test the generation of the line using the character string in the formal arg.
  #test <- a[ptr[10]]
  
  # Loop over all instances within file to comment out old, insert new

  # - place top of file in buffer and do first instance
  bb <- a[1:(ptr[1]-1)]
  # - check for commented out line
  tmpchk <- unlist(strsplit(a[ptr[1]],split="[[:punct:]]+"))
  if (length(tmpchk)==4){
    cmm <- paste("//  ",a[ptr[1]],sep="")
    bb <- append(bb,cmm)
    # - create the new lines to be inserted
    tmp <- unlist(strsplit(a[ptr[1]],split="[[:punct:]]+"))[4]
    tmp2 <- paste("  const char * str",1,";",sep="")     # create the const char*
    tmp3 <- paste("  str",1,"=\"",tmp,"\";",sep="")              # assign it the string
    tmp4 <- paste("  char* strx1=const_cast <char*> (str",1,");",sep="") #cast as a char*
    tmp5 <- c("  save_identifier_string(strx1);")         #supply this to the function
    c <- c(tmp2,tmp3,tmp4,tmp5)
    bb <- append(bb,c)
  } else {
    c <- a[ptr[1]]
    bb <- append(bb,c)
  }
  for(i in 2:length(ptr)){
    b <- a[(ptr[i-1]+1):(ptr[i]-1)]   # get code to line preceding next ptr
    # - check for commented out line
    tmpchk <- unlist(strsplit(a[ptr[i]],split="[[:punct:]]+"))
    if (length(tmpchk)==4){
      cmm <- paste("//  ",a[ptr[i]],sep="")
      b <- append(b,cmm)
      # - create the new lines to be inserted
      tmp <- unlist(strsplit(a[ptr[i]],split="[[:punct:]]+"))[4]
      tmp2 <- paste("  const char * str",i,";",sep="")     # create the const char*
      tmp3 <- paste("  str",i,"=\"",tmp,"\";",sep="")              # assign it the string
      tmp4 <- paste("  char* strx",i,"=const_cast <char*> (str",i,");",sep="") #cast as a char*
      tmp5 <- paste("  save_identifier_string(strx",i,");",sep="")         #supply this to the function
      c <- c(tmp2,tmp3,tmp4,tmp5)
      b <- append(b,c)
      bb <- append(bb,b)
    } else {
      c <- a[ptr[i]]
      b <- append(b,c)
      bb <- append(bb,b)
    }
  }
  # - append the bottom of the file
  b <- a[(ptr[length(ptr)]+1):length(a)]   # get code to bottom of file
  bb <- append(bb,b)
  writeLines(bb,con = flnms[ii])
}

################################################################################
