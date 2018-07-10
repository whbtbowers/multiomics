setwd("/project/home17/whb17/Documents/project3/project_files/preprocessing/p_ex/")

inp.prot_data <- read.csv("../../data/protein_data.csv")
inp.inf_id <- read.csv("../../data/id_info.csv")

#Remove weird superfluous columns
inp.prot_data <- inp.prot_data[,-c(25,26)]

#Separate off info columns
Pt <- inp.prot_data$Pt
inp.prot_data.body <- inp.prot_data[,-c(1,2)]

library("stringr")

# Remove uncategorised patients and those with missing labels
for (i in rev(1:nrow(inp.inf_id))){
  if (sum(is.na(inp.inf_id[i,])) > 0){
    inp.inf_id <- inp.inf_id[-i,]
  } else if (inp.inf_id$inf.status[i] == "Unassigned"){
    inp.inf_id <- inp.inf_id[-i,]
  } else if (inp.inf_id$inf.status[i] == "Excluded"){
    inp.inf_id <- inp.inf_id[-i,]
  } else if (inp.inf_id$array.id[i] == "no"){
    inp.inf_id <- inp.inf_id[-i,]
  }
}

# Get relevent indices

ind.inf_id <- c()
ind.prot_data <- c()

for (i in 1:nrow(inp.inf_id)){
  for (j in 1:length(Pt)){
    if (str_detect(toString(Pt[j]),toString(inp.inf_id$prot.id[i]))){
      ind.inf_id <- c(ind.inf_id, i)
      ind.prot_data <- c(ind.prot_data, j)
    }  
  }
}

#Get rows selected by indices
df.sel.inf_id <- inp.inf_id[ind.inf_id,]
df.sel.prot_data <-inp.prot_data.body[ind.prot_data,]

#Set up individual columns for disease classification
hiv.status <- c()
tb.status <- c()
group <- c()
nulabel <- c()

# Modify infection status to make easier to categorise
inf.status <- df.sel.inf_id$inf.status
inf.status <- as.character(inf.status)
inf.status[inf.status == "TB_HIV-"] <- "TB+/HIV-"
inf.status[inf.status == "S_TB_HIV-"] <- "TB+/HIV-"
inf.status[inf.status == "S_TB_HIV+"] <- "TB+/HIV+"
inf.status[inf.status == "TB_HIV+"] <- "TB+/HIV+"
inf.status[inf.status == "LTBI_HIV-"] <- "LTBI/HIV-"
inf.status[inf.status == "LTBI_long_term_HIV-"] <- "LTBI/HIV-"
inf.status[inf.status == "LTBI_HIV+"] <- "LTBI/HIV+"
inf.status[inf.status == "HIV+/Inf Not TB"] <- "OD/HIV+"
inf.status[inf.status == "Sick_control_HIV+"] <- "OD/HIV+"
inf.status[inf.status == "Excl_well_LTBI-_HIV+"] <- "OD/HIV+"
inf.status[inf.status == "HIV-/Inf Not TB"] <- "OD/HIV-"
inf.status[inf.status == "Excl_long_term_HIV-"] <- "OD/HIV-"
inf.status[inf.status == "Sick_control_HIV-"] <- "OD/HIV-"
inf.status[inf.status == "HIV-/Inf Not TB"] <- "OD/HIV-"
inf.status[inf.status == "Excl_long_term_HIV-"] <- "OD/HIV-"
inf.status[inf.status == "Excl_well_LTBI-_HIV-"] <- "OD/HIV-"

# Populate HIV column
for (i in 1:nrow(df.sel.inf_id)){
  if (str_detect(toString(inf.status[i]), "HIV-")){
    hiv.status <- c(hiv.status, "HIV-")
  } else {
    hiv.status <- c(hiv.status, "HIV+")
  }
}

# Populate TB column
for (i in 1:nrow(df.sel.inf_id)){
  if (str_detect(toString(inf.status[i]), "OD")){
    tb.status <- c(tb.status, "OD")
  } else if (str_detect(toString(df.sel.inf_id$inf.status[i]), "LTBI")){
    tb.status <- c(tb.status, "LTBI")
  } else {
    tb.status <- c(tb.status, "TB")
  }
}

# Populate group column
for (i in 1:length(hiv.status)){
  if (hiv.status[i] == "HIV-"){
    if (tb.status[i] == "TB"){
      group <- c(group, 1)
    } else if (tb.status[i] == "LTBI"){
      group <- c(group, 3)
    } else if (tb.status[i] == "OD"){
      group <- c(group, 6)
    }
  } else if (hiv.status[i] == "HIV+"){
    if (tb.status[i] == "TB"){
      group <- c(group, 2)
    } else if (tb.status[i] == "LTBI"){
      group <- c(group, 4)
    } else if (tb.status[i] == "OD"){
      group <- c(group, 5)
    }
  }
}

# Create new patient labels
for (i in 1:nrow(df.sel.inf_id)){
  label <- paste(group[i], "_", df.sel.inf_id$site[i], "_", df.sel.inf_id$prot.id[i] ,sep="")
  nulabel <- c(nulabel, label)
}

df.prot_data.ex <- cbind(data.frame(row.names=nulabel, df.sel.inf_id$prot.id, df.sel.inf_id$array.id, hiv.status, tb.status, group, df.sel.inf_id$site, df.sel.inf_id$sex), df.sel.prot_data)

colnames(df.prot_data.ex)[1] <- "prot.id"
colnames(df.prot_data.ex)[2] <- "array.id"
colnames(df.prot_data.ex)[6] <- "site"
colnames(df.prot_data.ex)[7] <- "sex"