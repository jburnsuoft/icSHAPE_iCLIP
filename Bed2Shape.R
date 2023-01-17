library(seqinr)
library(Biostrings)

#Define iCLIP and icSHAPE Coordinate Files
#clipfile <- "SRSF1.peaks.bed"
#shapeFile <- "SRSF1_icSHAPE.bed"

###Prepare Functions for getWriteSequence###

#Extract iCLIP peak coordinates from .bed file
getClipSites <- function(clipFile) {
  zfSites <- read.table(clipFile,
                        sep = "\t",
                        header = FALSE)
  colnames(zfSites) <- c("chrom", "start", "end", "score1",
                         "score2", "strand")
  return(zfSites)
}

#Extract icSHAPE range coverage coordinates
getShapeSites <- function(shapeFile) {
  file1Tab <- read.table(shapeFile,
                         sep = "\t",
                         header = FALSE)
  colnames(file1Tab) <- c("chrom", "start", "end", "score")
  return(file1Tab)
}

#Identify Overlapping iCLIP and icSHAPE Coverage
getShapeRange <- function(zfPeaks, shapeSites, strandSh, range) {
  runSites <- data.frame()
  for (x in 1:nrow(zfPeaks)) {
    for (i in -range:range) {
      site <- zfPeaks[x, ]$"start" + i
      hasRun <- TRUE
      if (!(site %in% shapeSites$"start") ||
          !(zfPeaks[x, ]$"chrom" %in% shapeSites[which(site == shapeSites$"start"), ]$"chrom")
          || !(strandSh == zfPeaks[x, ]$"strand")) {
        hasRun <- FALSE
        break
      }
    }   
    if (hasRun) {
      runSites <- rbind(runSites, zfPeaks[x, ])
    }
  }
  return(runSites)
}


#Extract Overlapping Coordinates from iCLIP and icSHAPE
getDataSites <- function(clipFile, plusShape, minusShape, range) {
  plusSites <- getShapeRange(clipFile, plusShape, "+", range)
  minusSites <- getShapeRange(clipFile, minusShape, "-", range)
  allSites <- rbind(plusSites, minusSites)
  
  return(allSites)
}

#Extract Corresponding icSHAPE Scores from icSHAPE Output File Folder Containing .SHAPE files
writeShape <- function(posChrom, clusPos, shapeData, initRange, runSize, seqName, shapeFolderName) {
  seqName <- paste(shapeFolderName, seqName, ".SHAPE", sep="")
  shapeTrimmed <- shapeData[which(((clusPos - initRange) <= shapeData$"start") & (shapeData$"start"<= (clusPos + runSize)) & (posChrom == shapeData$"chrom")), ]
  scores <- shapeTrimmed$"score"
  shapeOut <- data.frame(1:length(scores), scores)
  write.table(shapeOut, file = seqName, quote = FALSE,
              sep = "\t",
              row.names = FALSE,
              col.names = FALSE)
}

### Run getWriteSequence Function - Extracting icSHAPE Scores and Binding Sequences for Protein of Interest from iCLIP Peak Coordinate File and invivo/invitro icSHAPE Coordinates ###
## Example Command ##
#getWriteSequence("ZRANB2.bam.bed", 25, 0, "HEK293_NP_VIVO_PLUS_HG19.bed", "HEK293_NP_VIVO_MINUS_HG19.bed", "ZRANB2_NEW25.bed_VIVO_sequences_fasta.fa", "ZRANB2_NEW25.bed_SEQ_VIVO_20_0/", "ZRANB2_NEW25.bed_FULL_SHAPE_VIVO_20_0/", "ZRANB2_NEW25.bed_CT_VIVO_20_0/", "ZRANB2_NEW25.bed_FULL_DOT_VIVO_20_0/")

getWriteSequence <- function(clipFileName, shapeRange, clusterDistance, plusShapeName,
                             minusShapeName, outFastaName, seqFolderName,
                             shapeFolderName, ctFolderName, dbFolderName) {
  #load full genome fasta
  chrSeq <- readDNAStringSet("GRCh37.p13.genome.fa")
  
  dir.create(substr(seqFolderName, 1, (nchar(seqFolderName) - 1)))
  dir.create(substr(shapeFolderName, 1, (nchar(shapeFolderName) - 1)))
  dir.create(substr(ctFolderName, 1, (nchar(ctFolderName) - 1)))
  dir.create(substr(dbFolderName, 1, (nchar(dbFolderName) - 1)))
  
  #load clip file and shape files into data frames
  fullClip <- getClipSites(clipFile = clipFileName)
  shapePlus <- getShapeSites(plusShapeName)
  shapeMinus <- getShapeSites(minusShapeName)
  
  #trim clip sites to sites that have shape data available
  clip <- getDataSites(clipFile = fullClip, plusShape = shapePlus,
                     minusShape = shapeMinus, range = shapeRange)
  
  if (length(clip) == 0) {
    return("No binding sites with shape data found")
  }
  
  clip <- clip[order(clip$chrom, clip$strand, clip$start), ]
  clip <- data.table::data.table(clip)
  clip <- clip[ , diff := start - data.table::shift(start, type="lead"), by = list(chrom, strand)]
  clip$"diff" <- -clip$"diff"
  data.table::setDF(clip)
  
  #initialize vectors to hold sequences for output fasta and sequence names
  seqList <- vector("list")
  commandList1 <- c()
  commandList2 <- c()
  nameVec <- c()
  n <- 0
  iter <- 0
  #loop over sites in trimmed clip file and get sequences
  for (i in 1:nrow(clip)) {
    #skip later iterations if sequence added to list in previous
    if (n > 0) {
      n <- n - 1
      next
    }
    #convert chromosome names to same format as in genome fasta
    chrom <- clip[i, ]$"chrom"
    if (chrom != "chrM") {
      chromFasta <- paste(chrom, substr(chrom, 4, 5))
    } else {
      chromFasta <- "chrM MT"
    }
    #prepare to get sequence
    strand <- clip[i, ]$"strand"
    pos <- clip[i, ]$"start"
    sequence <- c()
    nameVec <- c(nameVec, paste(chrom, strand, pos))
    n <- 0
    num <- i
    count <- shapeRange
    while ((num < nrow(clip)) && (!is.na(clip[num, ]$"diff")) && (!is.na(clip[i, ]$"diff")) && (clip[num, ]$"diff" < clusterDistance)) {
      count <- count + shapeRange - clip[num, ]$"diff"
      n <- n + 1
      num <- i + n
      nameVec[length(nameVec)] <- paste(nameVec[length(nameVec)],
                                        clip[(num), ]$"start")
    }
    iter <- iter + 1
    x <- -shapeRange
    #get sequence for current clip site and adjacent clip sites
    while ((x <= count)) {
      fileNamePart <- paste(chrom, ".", pos, ".", strand, sep = "")
      if (strand == "+") {
        sequence <- c(sequence, as.character(complement(chrSeq[[chromFasta]][pos + x + 1])))
        if (x == count) {
          writeShape(posChrom = chrom, clusPos = pos, shapeData = shapePlus, initRange = shapeRange, runSize = count, seqName = fileNamePart, shapeFolderName)
        }
      } else {
        sequence <- c(sequence, as.character((chrSeq[[chromFasta]][pos + x + 1])))
        if (x == count) {
          writeShape(posChrom = chrom, clusPos = pos, shapeData = shapeMinus, initRange = shapeRange, runSize = count, seqName = fileNamePart, shapeFolderName)
        }
      }
      sequence[which(sequence == "T")] = "U"
      x <- x + 1
    }
    if ((pos - shapeRange) > 1) {
      if (strand == "+") {
        seqList[[iter]] <- reverse(sequence)
      } else {
        seqList[[iter]] <- sequence
      }
      write.fasta(sequence, nameVec[length(nameVec)], paste(seqFolderName, fileNamePart, ".seq", sep = ""))
      #command1 <- paste("./Fold.exe", paste(seqFolderName, fileNamePart, ".seq", sep = ""), paste(ctFolderName,fileNamePart , ".ct", sep = ""), "-sh", paste(shapeFolderName, fileNamePart, ".SHAPE", sep = ""), "-usm", "-3.6", "-si", "-2.16", "-sm", "4.86")
      #command1 <- paste("./ShapeKnots.exe", paste(seqFolderName, fileNamePart, ".seq", sep = ""), paste(ctFolderName,fileNamePart , ".ct", sep = ""), "-sh", paste(shapeFolderName, fileNamePart, ".SHAPE", sep = ""), "-si", "-2.16", "-sm", "4.86")
      #commandList1 <- c(commandList1, command1)
      #command2 <- paste("./ct2dot.exe", paste(ctFolderName, fileNamePart, ".ct", sep=""), -1, paste(dbFolderName, fileNamePart, ".dot", sep=""))
      #commandList2 <- c(commandList2, command2)
    }
  }
  #write a fasta with sequences around clip sites
  write.fasta(seqList, nameVec, outFastaName, nbchar = 200)
  #run structural predictions
  for (command in commandList1) {
    print(command)
    system(command)
  }
  for (command in commandList2) {
    print(command)
    system(command)
  }
}
