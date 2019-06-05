library(tidyverse)

library(stringr)   # 0-pad
library(Hmisc)     # rcorr()
library(ggplot2)   # plots
  library(ggcorrplot)
library(DescTools) # Fisher R-to-Z
library(RNifti)    # Read NIFTI header

library(R.matlab)  # Save .mat files for import to MATLAB

# Setup
project.home <- "/mnt/praxic/pdnetworksr01"
setwd(paste0(project.home, "/bin/flowanalysis/actflowmapping"))

mefc <- "_e00213_mefc_reoriented_toMNI.nii.gz"

# The seeds are given indices 1-264
seeds.n <- 1:264
seeds.i <- str_pad(seeds.n, 3, pad = 0)

# Get subjects from text file, to skip over those with old temp IDs
all.subjects <- read.table(paste0(project.home, 
                                   "/data/IDs/currentIDs.txt"))[, 1] %>%
                  sort()

# Return NA instead of failing on missing subject
read.tablex <- function(file, n) {
  
  if (file.exists(file)) {
    x <- read.table(file)[, 1]
  } else{
    x <- rep(NA, n)
  }
  
  return(x)
  
}

# Return values in array: (components)x(seeds); seeds have to be in the columns
# for rcorr()
read.values <- function(subj, scan = "rs", N = NA) {
  
  sdir <- paste0(project.home, "/subjects/", subj, "/session1/")
  
  # Number of components for missing files
  if (is.na(N)) {
    mefc <- readNifti(paste0(sdir, "/mrest/mrest", mefc))
    n.comp <- dim(mefc)[4]
  } else {
    n.comp <- N
  }
  
  table <- sapply(seeds.i, function(x) 
                    read.tablex(paste0(sdir, "flow/", scan, "/", x, ".txt"), 
                                n = n.comp)) %>%
            as.matrix()
  
  return(table)
  
}

# Check whether a subject is ready to process
check.readiness <- function(subj, file) {
  
  file <- paste0(project.home, "/subjects/", subj, "/session1/", file)
  
  return(file.exists(file))
  
}

# Remove all-0 columns from RS ROI table
remove.0 <- function(table) {
  
  which.0s <- apply(table, 2, mean) == 0
  n <- sum(which.0s)
  
  if (n > 0) {
    
    message("Removing ", n, " all-0 columns")
    the.range <- range(table) / 100
    
    new.rows <- matrix(runif(nrow(table) * n, 
                              min = the.range[1], max = the.range[2]),
                        nrow = nrow(table), ncol = n)
    
    table[, which.0s] <- new.rows
  
  }
  
  return(table)
  
}

################################################################################
# Select subjects: rest, mcvsa, mcvsm must be exactly equal

# Select subjects who have completed mrest MEICA and registered to MNI
good.mrest <- sapply(all.subjects, check.readiness, 
                     file = "mrest/mrest_e00213_mefc_reoriented_toMNI.nii.gz")
mrest.subjects <- all.subjects[good.mrest]

# Get which subjects have data; have to be ready for both 
mcvsa.ready <- sapply(all.subjects, check.readiness, 
                      file = "flow/mcvsa00/001.txt")
mcvsm00.ready <- sapply(all.subjects, check.readiness, 
                        file = "flow/mcvsm00/001.txt")
mcvsm12.ready <- sapply(all.subjects, check.readiness, 
                        file = "flow/mcvsm12/001.txt")
mcvs_.subjects <- all.subjects[(mcvsa.ready & mcvsm00.ready) & mcvsm12.ready]

# Everyone who has all three modalities
subjects <- intersect(mrest.subjects, mcvs_.subjects)

################################################################################
# RS data

# Resting state data for 100023
rest <- lapply(subjects, read.values) 

# If any subjects have an all-0 column, that means that ROI was outside of the 
# brain. Instead of tossing the subject, replace it with fake data two degrees
# of magnitude below the rest of the data. This should be uncorrelated enough
# to allow the subject to be included
rest.clean <- lapply(rest, remove.0)

# Rcorr requries > 4 rows; drop anyone too short, keep just R matrices
rest.nrow <- sapply(rest.clean, nrow)
rest.r <- lapply(rest.clean[rest.nrow > 4], function(x) rcorr(x)$r)

# Create and save 4D RS array as MATLAB object
# Third dim is 1 because there's only one correlation modality
rest.arr <- array(unlist(rest.r), dim = c(264, 264, 1, length(rest.r)))
writeMat("RS.mat", x = rest.arr)

# Save what RS subjects we used for the task analysis
# This is only removing; so RS/task will be the same
subjects <- subjects[rest.nrow > 4]
n.subjects <- length(subjects)

message("Using ", n.subjects, " subjects")

# Who is missing?
rest.na <- apply(rest.arr, 1:4, is.na)
rest.na.subj <- apply(rest.na, 4, sum) > 0

for (i in 1:58) {
  
  x <- sum(rest.na[[i]])
  message(x)
  
}

which(rest.na.subj)
subjects[rest.na.subj]

################################################################################
# Task 

tasks <- paste0("mcvs", c("a", "m"))

# Each one has only one value, so override read.values()
mcvsa <- sapply(subjects, read.values, scan = "mcvsa00", N = 1)
mcvsm00 <- sapply(subjects, read.values, scan = "mcvsm00", N = 1)
mcvsm12 <- sapply(subjects, read.values, scan = "mcvsm12", N = 1)

# Create new array with proper dimensions: region x task x subj
task.arr <- array(NA, dim = c(264, 3, n.subjects))

# Reshape array: rois, modality, subject
for (roi in seeds.n) {
  for (subj in 1:n.subjects) {
    
    task.arr[roi, 1, subj] <- mcvsa[roi, subj]
    task.arr[roi, 2, subj] <- mcvsm00[roi, subj]
    task.arr[roi, 3, subj] <- mcvsm12[roi, subj]
    
  }
}

# Check to see if NAs ahave been included
sum(is.na(task.arr)) / prod(dim(task.arr))

writeMat("task.mat", x = task.arr)

################################################################################
# Some pretty plots

# RS

# Reduce along subjects
rest.mean <- apply(rest.arr, 1:2, mean)
rest.mean[rest.mean == 1] <- NA
rest.mean.Z <- FisherZ(rest.mean)
rest.mean.Z.melt <- melt(rest.mean.Z)

ggplot(rest.mean.Z.melt, aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  scale_fill_gradient2(low = "red4", high = "navyblue", mid = "white", 
                       midpoint = 0, na.value = "white") +
  theme_bw() +
  theme(axis.title = element_blank(), axis.ticks = element_blank(), 
        axis.text = element_blank(), panel.border = element_blank()) 

# Task

mcvsa.mean <- apply(mcvsa, 1, mean)
mcvsm00.mean <- apply(mcvsm00, 1, mean)
mcvsm12.mean <- apply(mcvsm12, 1, mean)

task.mean.values <- cbind(mcvsa.mean, mcvsm00.mean, mcvsm12.mean)
task.mean.vals.melt <- melt(task.mean.values)

ggplot(task.mean.vals.melt, aes(x = Var1, y = value, color = Var2)) +
  geom_line() +
  labs(x = "Power ROI Index", y = "t Statistic", color = "Stimulus") +
  scale_color_discrete(labels = c("0-back (att.)", "0-back (mem.)", 
                                  "n-back"))
