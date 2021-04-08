#####
# APLC getting post 2017 data
#####
rm(list=ls())
library(data.table)
# Put in your actual path where the text files are saved

files <- list.files(recursive = T)



# Read the files in, assuming comma separator
txt_files_df <- lapply(files, function(x) {read.table(file = x, header = T, sep =",")})





DT <- rbindlist(sapply(files, fread, simplify = FALSE),
                use.names = FALSE, idcol = "FileName")
summary(as.factor(DT$Instar1st))

# Combine them
combined_df <- do.call("rbind", lapply(txt_files_df, as.data.frame)) 






list.files(R.home())
## Only files starting with a-l or r
## Note that a-l is locale-dependent, but using case-insensitive
## matching makes it unambiguous in English locales
dir("../..", pattern = "^[a-lr]", full.names = TRUE, ignore.case = TRUE)

list.dirs(R.home("doc"))
list.dirs(R.home("doc"), full.names = FALSE)