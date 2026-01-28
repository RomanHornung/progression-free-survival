# 0. Dataset and libraries----
## 0.1 Loading libraries----
source("./R_Code/libraries.R")

## 0.2 Loading dataset----
# df <- read_excel("./Data/Blutwerte_PFS.xlsx", skip = 1)
# df <- as.data.frame(df)
load(file = "./Data/pseudo_data.RData")

# 1.  Preparation----
## 1.1 Renaming----
### 1.1.1 PFS_d / progression free survival [days]----
names(df)[names(df) == "PFS_d"] <- "time"

### 1.1.2 status----
#### survival = überleben
#### progrediente Erkrankung = progressive disease
#### censored = zensiert
names(df)[names(df) == "Status_PFS (0=überleben;  1= progrediente Erkrankung; 3=unbekannt (zensiert))"] <- "status"
df$status <- factor(df$status, labels = c("survival", "progressive disease", "survival"))

### 1.1.3 Leucocytes----
names(df)[names(df) == "Leukozyten"] <- "leuc"

### 1.1.4 Lymphocytes----
names(df)[names(df) == "Lymphozyten"] <- "lymph"

### 1.1.5 neutrophil granulocytes----
names(df)[names(df) == "neutrophile Granulozyten"] <- "n_gran"

### 1.1.6 monocytes----
names(df)[names(df) == "Monozyten"] <- "mono"

### 1.1.7 eosinophil granulocytes----
names(df)[names(df) == "eosinophile Granulozyten"] <- "e_gran"

### 1.1.8 CRP----
names(df)[names(df) == "CRP"] <- "crp"

### 1.1.9 Albumin----
names(df)[names(df) == "Albumin"] <- "albumin"

### 1.1.10 Protein----
names(df)[names(df) == "Protein"] <- "protein"

### 1.1.11 LDH----
names(df)[names(df) == "LDH"] <- "ldh"

### 1.1.12 magnesium----
names(df)[names(df) == "Magnesium"] <- "mg"

# 2. Saving----
save(df, file = "./Data/prep.RData")
