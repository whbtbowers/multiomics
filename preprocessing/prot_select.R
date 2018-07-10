setwd("/project/home17/whb17/Documents/project3/project_files/preprocessing/")

inp.prot_data <- read.csv("../data/protein_data.csv")
inp.inf_id <- read.csv("../data/prot_id-inf.csv")

library("stringr")

# Remove uncategorised patients
for (i in rev(1:nrow(inp.inf_id))){
  if (inp.inf_id$inf_status[i] == "Unassigned"){
    inp.inf_id <- inp.inf_id[-i,]
  } else if (inp.inf_id$inf_status[i] == "Excluded"){
    inp.inf_id <- inp.inf_id[-i,]
  }
}

# Get indices of patients in either category
ind.inf.hiv_pos <- c()
ind.inf.hiv_neg <- c()


# Filter patients based on HIV status
for (i in 1:nrow(inp.inf_id)){
  if (str_detect(inp.inf_id$inf_status[i], "HIV-")){
    ind.inf.hiv_neg <- c(ind.inf.hiv_neg, i)
  } else if (str_detect(inp.inf_id$inf_status[i], "HIV+")){
    ind.inf.hiv_pos <- c(ind.inf.hiv_pos, i)
  }
}

inf.hiv_pos <- inp.inf_id[ind.inf.hiv_pos,]
inf.hiv_neg <- inp.inf_id[ind.inf.hiv_neg,]

#Get indices of prot info
ind.prot.all <-c()
ind.prot.hiv_pos <- c()
ind.prot.hiv_neg <- c()

#Get ordered lists of infection category
prot.inf_cat.all <-c()
prot.inf_cat.hiv_pos <- c()
prot.inf_cat.hiv_neg <- c()

#Extract all patients based on ID
for (i in 1:nrow(inp.prot_data)){
  for (j in 1:nrow(inp.inf_id)){
    if (str_detect(inp.prot_data$Pt[i], toString(inp.inf_id$ID[j]))){
      ind.prot.all <- c(ind.prot.all, i)
      prot.inf_cat.all <- c(prot.inf_cat.all, as.character(inp.inf_id$inf_status[j]))
    }
  }
}

#Extract HIV+ patients based on ID
for (i in 1:nrow(inp.prot_data)){
  for (j in 1:nrow(inf.hiv_pos)){
    if (str_detect(inp.prot_data$Pt[i], toString(inf.hiv_pos$ID[j]))){
      ind.prot.hiv_pos <- c(ind.prot.hiv_pos, i)
      prot.inf_cat.hiv_pos <- c(prot.inf_cat.hiv_pos, as.character(inf.hiv_pos$inf_status[j]))
    }
  }
}

#Extract HIV- patients based on ID
for (i in 1:nrow(inp.prot_data)){
  for (j in 1:nrow(inf.hiv_neg)){
    if (str_detect(inp.prot_data$Pt[i], toString(inf.hiv_neg$ID[j]))){
      ind.prot.hiv_neg <- c(ind.prot.hiv_neg, i)
      prot.inf_cat.hiv_neg <- c(prot.inf_cat.hiv_neg, as.character(inf.hiv_neg$inf_status[j]))
    }
  }
}

# Get infection statuses of chosen patients
inf.status.all <- data.frame(prot.inf_cat.all)
inf.status.hiv_pos <- data.frame(prot.inf_cat.hiv_pos)
inf.status.hiv_neg <- data.frame(prot.inf_cat.hiv_neg)

# Grab data using indices
prot.data.all <- inp.prot_data[ind.prot.all,]
prot_data.hiv_pos <- inp.prot_data[ind.prot.hiv_pos,]
prot_data.hiv_neg <- inp.prot_data[ind.prot.hiv_neg,]

# Add infection status column
df.all <- cbind(inf.status.all, prot.data.all)
df.hiv_pos <- cbind(inf.status.hiv_pos, prot_data.hiv_pos)
df.hiv_neg <- cbind(inf.status.hiv_neg, prot_data.hiv_neg)

#Rename all infection status columns
colnames(df.all)[1] <- c("inf.status")
colnames(df.hiv_pos)[1] <- c("inf.status")
colnames(df.hiv_neg)[1] <- c("inf.status")

##Convert to category number
#df.hiv_plus$inf.status <- as.character(df.hiv_plus$inf.status)
#df.hiv_plus$inf.status[df.hiv_plus$inf.status == "TB+/HIV+"] <- 1
#df.hiv_plus$inf.status[df.hiv_plus$inf.status == "LTBI/HIV+"] <- 2
#df.hiv_plus$inf.status[df.hiv_plus$inf.status == "HIV+/Inf Not TB"] <- 3

#Create CSV
write.csv(df.all, file="../data/extraction_retry/protein_data_hiv_all.csv")
write.csv(df.hiv_pos, file="../data/extraction_retry/protein_data_hiv_pos.csv")
write.csv(df.hiv_neg, file="../data/extraction_retry/protein_data_hiv_neg.csv")

