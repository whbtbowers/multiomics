setwd("/project/home17/whb17/Documents/project3/project_files/preprocessing/ex_4/")

library(stringr)

inp.gene_data <- read.csv("../../data/gene_data.csv", header=TRUE)
inp.inf_id <- read.csv("../../data/id_info.csv")

init.id <- inp.gene_data$array.id
inp.gene_data.body <- inp.gene_data[,-1]

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
ind.gene_data <- c()

for (i in 1:nrow(inp.inf_id)){
  for (j in 1:length(init.id)){
    if (toString(init.id[j]) == toString(inp.inf_id$array.id[i])){
      ind.inf_id <- c(ind.inf_id, i)
      ind.gene_data <- c(ind.gene_data, j)
    }  
  }
}

#Get rows selected by indices
df.sel.inf_id <- inp.inf_id[ind.inf_id,]
df.sel.gene_data <-inp.gene_data.body[ind.gene_data,]

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
  ind_tag <- paste(group[i], "_", df.sel.inf_id$site[i], "_", df.sel.inf_id$prot.id[i] ,sep="")
  nulabel <- c(nulabel, ind_tag)
}

df.gene_data.ex <- cbind(data.frame(row.names=nulabel, df.sel.inf_id$prot.id, df.sel.inf_id$array.id, hiv.status, tb.status, group, df.sel.inf_id$site, df.sel.inf_id$sex), df.sel.gene_data)

colnames(df.gene_data.ex)[1] <- "prot.id"
colnames(df.gene_data.ex)[2] <- "array.id"
colnames(df.gene_data.ex)[6] <- "site"
colnames(df.gene_data.ex)[7] <- "sex"

#Remove Malawian patients
for (i in rev(1:nrow(df.gene_data.ex))){
  if (df.gene_data.ex$site[i] == "ML"){
    df.gene_data.ex <- df.gene_data.ex[-c(i),]
  }
}

#sum(
#length(df.gene_data.ex$group[df.gene_data.ex$group==1]) * 0.3
#,
#length(df.gene_data.ex$group[df.gene_data.ex$group==2]) * 0.3
#,
#length(df.gene_data.ex$group[df.gene_data.ex$group==3]) * 0.3
#,
#length(df.gene_data.ex$group[df.gene_data.ex$group==4]) * 0.3
#,
#length(df.gene_data.ex$group[df.gene_data.ex$group==5]) * 0.3
#,
#length(df.gene_data.ex$group[df.gene_data.ex$group==6]) * 0.3
#)

# Write to .csv file
write.csv(df.gene_data.ex[,-c(1:7)],"../../data/ex_4/gene_data_body.csv",row.names=TRUE)
write.csv(df.gene_data.ex[,c(1:7)],"../../data/ex_4/gene_data_meta.csv",row.names=TRUE)

