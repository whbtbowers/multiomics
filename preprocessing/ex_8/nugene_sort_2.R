setwd("/home/whb17/Documents/project3/project_files/preprocessing/ex_8/")

library(stringr)

#Probe file

print("Opening gene expression data. This will take a while")

df.gene.init <- read.csv("../../data/ex_8/tot_gene_data_phase2.csv", header=TRUE, row.names = 1)

#Gene info file
inp.inf_id <- read.csv("../../data/id_info.csv", header=TRUE)

# Maybe set as '.sub' a little later
#df.gene.init.sub <- df.gene_init[1:100, 1:100]

#Create list of array IDs
init.id.wX <- rownames(df.gene.init)
init.id <- c()

for (i in 1:length(init.id.wX)){
  init.id <- c(init.id, strsplit(init.id.wX[i], "X")[[1]][2])
}

# Remove uncategorised patients and those with missing labels

print("Removing uncategorised patients and those with missing labels")

for (i in rev(1:nrow(inp.inf_id))){
  if (inp.inf_id$inf.status[i] == "Unassigned"){
    inp.inf_id <- inp.inf_id[-i,]
  } else if (inp.inf_id$inf.status[i] == "Excluded"){
    inp.inf_id <- inp.inf_id[-i,]
  } else if (inp.inf_id$array.id[i] == "no"){
    inp.inf_id <- inp.inf_id[-i,]
  }
}

# Get relevent indices
print("Getting indices of relevent patients")
ind.inf_id <- c()
ind.gene_data <- c()

for (i in 1:length(init.id)){
  if (init.id[i] %in% inp.inf_id$array.id){
    sgl.ind.inf_id <- match(init.id[i], inp.inf_id$array.id)
    ind.inf_id <- c(ind.inf_id, sgl.ind.inf_id)
    ind.gene_data <- c(ind.gene_data, match(inp.inf_id$array.id[sgl.ind.inf_id], init.id))
  }  
}

#Get rows selected by indices
df.sel.inf_id <- inp.inf_id[ind.inf_id,]
df.sel.gene_data <-df.gene.init[ind.gene_data,]

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

print("Populating new HIV status column")

for (i in 1:nrow(df.sel.inf_id)){
  if (str_detect(toString(inf.status[i]), "HIV-")){
    hiv.status <- c(hiv.status, "HIV-")
  } else {
    hiv.status <- c(hiv.status, "HIV+")
  }
}

# Populate TB column

print("Populating new TB status column")

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

print("Populating new group number column")

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

print("Constructing new patient labels")

for (i in 1:nrow(df.sel.inf_id)){
  ind_tag <- paste(group[i], "_", df.sel.inf_id$site[i], "_", df.sel.inf_id$prot.id[i] ,sep="")
  nulabel <- c(nulabel, ind_tag)
}

print("Constructing new dataset and metadata")

df.gene_data.ex <- cbind(data.frame(row.names=nulabel, df.sel.inf_id$prot.id, df.sel.inf_id$array.id, hiv.status, tb.status, group, df.sel.inf_id$site, df.sel.inf_id$sex), df.sel.gene_data)

colnames(df.gene_data.ex)[1] <- "prot.id"
colnames(df.gene_data.ex)[2] <- "array.id"
colnames(df.gene_data.ex)[6] <- "site"
colnames(df.gene_data.ex)[7] <- "sex"

df.gene_ex.body <- df.gene_data.ex[,-c(1:7)]
df.gene_ex.meta <- df.gene_data.ex[,c(1:7)]

# Write to .csv file

print("Writing data and meta to CSV...")

write.csv(df.gene_ex.body,"../../data/ex_8/gene_init_data_body.csv",row.names=TRUE)
write.csv(df.gene_ex.meta,"../../data/ex_8/gene_data_meta.csv",row.names=TRUE)

print("DONE")

for (i in c("HIV+", "HIV-")){
  print(paste(i, ": ", length(df.gene_ex.meta$hiv.status[df.gene_ex.meta$hiv.status == i]), sep=""))
}

for (j in 1:6){
  print(paste(j, ": ", length(df.gene_ex.meta$group[df.gene_ex.meta$group == j])))
}