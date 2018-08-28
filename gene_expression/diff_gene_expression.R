library('edgeR')

df <- read.table("sample_expression_matrix.csv", header = TRUE, sep = "," )

rownames(df) <- df$X

df <- subset(df, select = -c(X))

cpm_log <- cpm(df, log = TRUE)
median_log2_cpm <- apply(cpm_log, 1, median)
expr_cutoff <- -2.5
df_clean <- df[median_log2_cpm > expr_cutoff, ]
cpm_log <- cpm(df_clean, log = TRUE)

pca <- prcomp(t(cpm_log), scale. = TRUE)

pca_df <- data.frame(pca$x)

group_membership <- c()

group1 = read.table("group3.csv", header = T, sep=",")
group1_samples <- as.character(group1$Sample_Name)

group2 = read.table("group4.csv", header = T, sep=",")
group2_samples <- as.character(group2$Sample_Name)

all_samples = row.names(pca_df)

for (row in all_samples ){
  
  row <- strsplit(row, "_")[[1]][1]
  
  if (row %in% group1_samples) {
    
    group_membership <- c(group_membership, 1)
    
  } else if (row %in% group2_samples) {
    
    group_membership <- c(group_membership, 2)
  } else {
    
    group_membership <- c(group_membership, 0)
    
  }
  
}

pca_df$group_membership <- group_membership

pca_df <- pca_df[pca_df$group_membership != 0,]




group_membership <- c()

all_samples = colnames(df)

for (row in all_samples ){
  
  row <- strsplit(row, "_")[[1]][1]
  
  if (row %in% group1_samples) {
    
    group_membership <- c(group_membership, 1)
    
  } else if (row %in% group2_samples) {
    
    group_membership <- c(group_membership, 2)
  } else {
    
    group_membership <- c(group_membership, 0)
    
  }
  
}


df_clean <- df_clean[,group_membership!=0]

group_membership <-group_membership[group_membership!=0]
y <- DGEList(counts = df_clean, group = group_membership)
y <- calcNormFactors(y)
y <- estimateDisp(y)
sqrt(y$common.dispersion)
et <- exactTest(y)
results_edgeR <- topTags(et, n = nrow(df_clean), sort.by = "none")

sig_df <- results_edgeR$table[(results_edgeR$table$FDR <0.05) &  (results_edgeR$table$PValue <0.05), ]




