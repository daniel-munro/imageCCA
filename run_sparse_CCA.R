library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(forcats)

load_image_data <- function(file_name) {
  d <- read_tsv(file_name, col_types=cols(.default='d'), col_names=F) %>% as.matrix()
  colnames(d) <- 1:ncol(d)
  d
}

load_expr_data <- function(file_name, genes_file_name) {
  genes <- read_lines(genes_file_name)
  read_tsv(file_name, col_types=cols(.default='d'), col_names=genes) %>% as.matrix()
}

SCCA <- function(datasets, penalty_x, penalty_z, K) {
  datasets[[1]] <- scale(datasets[[1]])
  datasets[[2]] <- scale(datasets[[2]])
  x <- PMA::CCA(datasets[[1]], datasets[[2]], penaltyx=penalty_x, penaltyz=penalty_z, K=K, standardize=F)
  rownames(x$u) <- colnames(datasets[[1]])
  rownames(x$v) <- colnames(datasets[[2]])
  x$group_names <- names(datasets)
  x$CCA_var1 <- datasets[[1]] %*% x$u
  x$CCA_var2 <- datasets[[2]] %*% x$v
  x
}

SCCA_coefs <- function(cca) {
  data_frame(CCA_var = 1:cca$K) %>%
    group_by(CCA_var) %>%
    do({
      k <- .$CCA_var
      nonzero_1 <- which(cca$u[, k] != 0)
      nonzero_2 <- which(cca$v[, k] != 0)
      bind_rows(
        data_frame(
          type = cca$group_names[1],
          name = rownames(cca$u)[nonzero_1],
          coefficient = cca$u[nonzero_1, k]
        ),
        
        data_frame(
          type = cca$group_names[2],
          name = rownames(cca$v)[nonzero_2],
          coefficient = cca$v[nonzero_2, k]
        )
      )
    }) %>%
    ungroup
}

SCCA_cors <- function(cca) data_frame(k = 1:length(cca$d), d = cca$d)


##### BRCA

###  SIFT

image <- load_image_data("BRCA/SIFT/CCA_input_image.txt")
expr <- load_expr_data("BRCA/SIFT/CCA_input_expr.txt", "BRCA/SIFT/CCA_genes.txt")

# high sparsity
cca <- SCCA(list(gene = expr, image = image), penalty_x=0.01, penalty_z=0.15, K=100)
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("BRCA/SIFT/CCA_hs_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("BRCA/SIFT/CCA_hs_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("BRCA/SIFT/CCA_hs_coefs.txt")

# low sparsity
cca <- SCCA(list(gene = expr, image = image), penalty_x=0.05, penalty_z=0.15, K=100)
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("BRCA/SIFT/CCA_ls_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("BRCA/SIFT/CCA_ls_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("BRCA/SIFT/CCA_ls_coefs.txt")


### SIFT LDA

image <- load_image_data("BRCA/SIFT_LDA/CCA_input_image.txt")
expr <- load_expr_data("BRCA/SIFT_LDA/CCA_input_expr.txt", "BRCA/SIFT_LDA/CCA_genes.txt")

# high sparsity
cca <- SCCA(list(gene = expr, image = image), penalty_x=0.01, penalty_z=0, K=2)
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("BRCA/SIFT_LDA/CCA_hs_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("BRCA/SIFT_LDA/CCA_hs_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("BRCA/SIFT_LDA/CCA_hs_coefs.txt")

# low sparsity
cca <- SCCA(list(gene = expr, image = image_LDA), penalty_x=0.05, penalty_z=0, K=2)
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("BRCA/SIFT_LDA/CCA_ls_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("BRCA/SIFT_LDA/CCA_ls_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("BRCA/SIFT_LDA/CCA_ls_coefs.txt")


### CAE

image <- load_image_data("BRCA/CAE/CCA_input_image.txt")
expr <- load_expr_data("BRCA/CAE/CCA_input_expr.txt", "BRCA/CAE/CCA_genes.txt")

cca <- SCCA(list(gene = expr, image = image), penalty_x=0.05, penalty_z=0.015, K=100)
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("BRCA/CAE/CCA_ls_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("BRCA/CAE/CCA_ls_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("BRCA/CAE/CCA_ls_coefs.txt")


### CAE PCA

image <- load_image_data("BRCA/CAE_PCA/CCA_input_image.txt")
expr <- load_expr_data("BRCA/CAE_PCA/CCA_input_expr.txt", "BRCA/CAE_PCA/CCA_genes.txt")

cca <- SCCA(list(gene = expr, image = image), penalty_x=0.05, penalty_z=0.03, K=100)
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("BRCA/CAE_PCA/CCA_ls_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("BRCA/CAE_PCA/CCA_ls_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("BRCA/CAE_PCA/CCA_ls_coefs.txt")


### CAE2 PCA

image <- load_image_data("BRCA/CAE2_PCA/CCA_input_image.txt")
expr <- load_expr_data("BRCA/CAE2_PCA/CCA_input_expr.txt", "BRCA/CAE2_PCA/CCA_genes.txt")

cca <- SCCA(list(gene = expr, image = image), penalty_x=0.05, penalty_z=0.1, K=10)
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("BRCA/CAE2_PCA/CCA_ls_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("BRCA/CAE2_PCA/CCA_ls_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("BRCA/CAE2_PCA/CCA_ls_coefs.txt")


### CAE3

image <- load_image_data("BRCA/CAE3/CCA_input_image.txt")
expr <- load_expr_data("BRCA/CAE3/CCA_input_expr.txt", "BRCA/CAE3/CCA_genes.txt")

cca <- SCCA(list(gene = expr, image = image), penalty_x=0.05, penalty_z=0.15, K=100)
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("BRCA/CAE3/CCA_ls_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("BRCA/CAE3/CCA_ls_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("BRCA/CAE3/CCA_ls_coefs.txt")


### CAE3 PCA

image <- load_image_data("BRCA/CAE3_PCA/CCA_input_image.txt")
expr <- load_expr_data("BRCA/CAE3_PCA/CCA_input_expr.txt", "BRCA/CAE3_PCA/CCA_genes.txt")

cca <- SCCA(list(gene = expr, image = image), penalty_x=0.05, penalty_z=0.15, K=100)
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("BRCA/CAE3_PCA/CCA_ls_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("BRCA/CAE3_PCA/CCA_ls_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("BRCA/CAE3_PCA/CCA_ls_coefs.txt")


### CAE3 discriminative

image <- load_image_data("BRCA/CAE3_discrim/CCA_input_image.txt")
expr <- load_expr_data("BRCA/CAE3_discrim/CCA_input_expr.txt", "BRCA/CAE3_discrim/CCA_genes.txt")

cca <- SCCA(list(gene = expr, image = image), penalty_x=0.05, penalty_z=0.15, K=ncol(image))
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("BRCA/CAE3_discrim/CCA_ls_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("BRCA/CAE3_discrim/CCA_ls_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("BRCA/CAE3_discrim/CCA_ls_coefs.txt")


### CAE3 discriminative PCA

image <- load_image_data("BRCA/CAE3_discrim_PCA/CCA_input_image.txt")
expr <- load_expr_data("BRCA/CAE3_discrim_PCA/CCA_input_expr.txt", "BRCA/CAE3_discrim_PCA/CCA_genes.txt")

cca <- SCCA(list(gene = expr, image = image), penalty_x=0.05, penalty_z=0.15, K=ncol(image))
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("BRCA/CAE3_discrim_PCA/CCA_ls_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("BRCA/CAE3_discrim_PCA/CCA_ls_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("BRCA/CAE3_discrim_PCA/CCA_ls_coefs.txt")

##### trying other discrim variants
image <- load_image_data("BRCA/CAE3_discrim4_PCA/CCA_input_image.txt")
expr <- load_expr_data("BRCA/CAE3_discrim4_PCA/CCA_input_expr.txt", "BRCA/CAE3_discrim4_PCA/CCA_genes.txt")

cca <- SCCA(list(gene = expr, image = image), penalty_x=0.05, penalty_z=0.15, K=ncol(image))
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("BRCA/CAE3_discrim4_PCA/CCA_ls_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("BRCA/CAE3_discrim4_PCA/CCA_ls_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("BRCA/CAE3_discrim4_PCA/CCA_ls_coefs.txt")
#####


##### LGG

### SIFT

image <- load_image_data("LGG/SIFT/CCA_input_image.txt")
expr <- load_expr_data("LGG/SIFT/CCA_input_expr.txt", "LGG/SIFT/CCA_genes.txt")

# high sparsity
cca <- SCCA(list(gene = expr, image = image), penalty_x=0.01, penalty_z=0.15, K=100)
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("LGG/SIFT/CCA_hs_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("LGG/SIFT/CCA_hs_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("LGG/SIFT/CCA_hs_coefs.txt")

# low sparsity
cca <- SCCA(list(gene = expr, image = image), penalty_x=0.05, penalty_z=0.15, K=100)
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("LGG/SIFT/CCA_ls_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("LGG/SIFT/CCA_ls_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("LGG/SIFT/CCA_ls_coefs.txt")


### SIFT PCA

image <- load_image_data("LGG/SIFT_PCA/CCA_input_image.txt")
expr <- load_expr_data("LGG/SIFT_PCA/CCA_input_expr.txt", "LGG/SIFT_PCA/CCA_genes.txt")

# high sparsity
cca <- SCCA(list(gene = expr, image = image), penalty_x=0.01, penalty_z=0.15, K=49)
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("LGG/SIFT_PCA/CCA_hs_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("LGG/SIFT_PCA/CCA_hs_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("LGG/SIFT_PCA/CCA_hs_coefs.txt")

# low sparsity
cca <- SCCA(list(gene = expr, image = image), penalty_x=0.05, penalty_z=0.15, K=49)
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("LGG/SIFT_PCA/CCA_ls_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("LGG/SIFT_PCA/CCA_ls_canonvar_image.txt", col_names=F)
cca %>%SCCA_coefs() %>% write_tsv("LGG/SIFT_PCA/CCA_ls_coefs.txt")


### CAE3

image <- load_image_data("LGG/CAE3/CCA_input_image.txt")
expr <- load_expr_data("LGG/CAE3/CCA_input_expr.txt", "LGG/CAE3/CCA_genes.txt")

cca <- SCCA(list(gene = expr, image = image), penalty_x=0.05, penalty_z=0.15, K=100)
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("LGG/CAE3/CCA_ls_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("LGG/CAE3/CCA_ls_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("LGG/CAE3/CCA_ls_coefs.txt")


### CAE3 PCA

image <- load_image_data("LGG/CAE3_PCA/CCA_input_image.txt")
expr <- load_expr_data("LGG/CAE3_PCA/CCA_input_expr.txt", "LGG/CAE3_PCA/CCA_genes.txt")

cca <- SCCA(list(gene = expr, image = image), penalty_x=0.05, penalty_z=0.15, K=100)
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("LGG/CAE3_PCA/CCA_ls_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("LGG/CAE3_PCA/CCA_ls_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("LGG/CAE3_PCA/CCA_ls_coefs.txt")

# Trying different sparsity to improve CCA correlations:
cca <- SCCA(list(gene = expr, image = image), penalty_x=0.1, penalty_z=0.5, K=100)
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("LGG/CAE3_PCA/CCA_xs_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("LGG/CAE3_PCA/CCA_xs_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("LGG/CAE3_PCA/CCA_xs_coefs.txt")


### CAE3 discriminative

image <- load_image_data("LGG/CAE3_discrim/CCA_input_image.txt")
expr <- load_expr_data("LGG/CAE3_discrim/CCA_input_expr.txt", "LGG/CAE3_discrim/CCA_genes.txt")

cca <- SCCA(list(gene = expr, image = image), penalty_x=0.05, penalty_z=0.15, K=ncol(image))
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("LGG/CAE3_discrim/CCA_ls_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("LGG/CAE3_discrim/CCA_ls_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("LGG/CAE3_discrim/CCA_ls_coefs.txt")


### CAE3 discriminative PCA

image <- load_image_data("LGG/CAE3_discrim_PCA/CCA_input_image.txt")
expr <- load_expr_data("LGG/CAE3_discrim_PCA/CCA_input_expr.txt", "LGG/CAE3_discrim_PCA/CCA_genes.txt")

cca <- SCCA(list(gene = expr, image = image), penalty_x=0.05, penalty_z=0.15, K=ncol(image))
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("LGG/CAE3_discrim_PCA/CCA_ls_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("LGG/CAE3_discrim_PCA/CCA_ls_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("LGG/CAE3_discrim_PCA/CCA_ls_coefs.txt")


##### trying other discrim variants
image <- load_image_data("LGG/CAE3_discrim2_PCA/CCA_input_image.txt")
expr <- load_expr_data("LGG/CAE3_discrim2_PCA/CCA_input_expr.txt", "LGG/CAE3_discrim2_PCA/CCA_genes.txt")

cca <- SCCA(list(gene = expr, image = image), penalty_x=0.05, penalty_z=0.15, K=9)#ncol(image))
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("LGG/CAE3_discrim2_PCA/CCA_ls_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("LGG/CAE3_discrim2_PCA/CCA_ls_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("LGG/CAE3_discrim2_PCA/CCA_ls_coefs.txt")


image <- load_image_data("LGG/CAE3_discrim3_PCA/CCA_input_image.txt")[, 1]  # 2nd column is all NA
expr <- load_expr_data("LGG/CAE3_discrim3_PCA/CCA_input_expr.txt", "LGG/CAE3_discrim3_PCA/CCA_genes.txt")
# doesn't work: requires >1 image feature
cca <- SCCA(list(gene = expr, image = image), penalty_x=0.05, penalty_z=0.15, K=ncol(image))
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("LGG/CAE3_discrim3_PCA/CCA_ls_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("LGG/CAE3_discrim3_PCA/CCA_ls_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("LGG/CAE3_discrim3_PCA/CCA_ls_coefs.txt")
#####


##### GTEx

### CAE3

image <- load_image_data("GTEx/CAE3/CCA_input_image.txt")
expr <- load_expr_data("GTEx/CAE3/CCA_input_expr.txt", "GTEx/CAE3/CCA_genes.txt")

cca <- SCCA(list(gene = expr, image = image), penalty_x=0.05, penalty_z=0.15, K=100)
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("GTEx/CAE3/CCA_ls_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("GTEx/CAE3/CCA_ls_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("GTEx/CAE3/CCA_ls_coefs.txt")


### CAE3 PCA

image <- load_image_data("GTEx/CAE3_PCA/CCA_input_image.txt")
expr <- load_expr_data("GTEx/CAE3_PCA/CCA_input_expr.txt", "GTEx/CAE3_PCA/CCA_genes.txt")

cca <- SCCA(list(gene = expr, image = image), penalty_x=0.05, penalty_z=0.15, K=100)
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("GTEx/CAE3_PCA/CCA_ls_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("GTEx/CAE3_PCA/CCA_ls_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("GTEx/CAE3_PCA/CCA_ls_coefs.txt")


### CAE3 discriminative

image <- load_image_data("GTEx/CAE3_discrim/CCA_input_image.txt")
expr <- load_expr_data("GTEx/CAE3_discrim/CCA_input_expr.txt", "GTEx/CAE3_discrim/CCA_genes.txt")

cca <- SCCA(list(gene = expr, image = image), penalty_x=0.05, penalty_z=0.15, K=ncol(image))
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("GTEx/CAE3_discrim/CCA_ls_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("GTEx/CAE3_discrim/CCA_ls_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("GTEx/CAE3_discrim/CCA_ls_coefs.txt")


### CAE3 discriminative PCA

image <- load_image_data("GTEx/CAE3_discrim_PCA/CCA_input_image.txt")
expr <- load_expr_data("GTEx/CAE3_discrim_PCA/CCA_input_expr.txt", "GTEx/CAE3_discrim_PCA/CCA_genes.txt")

cca <- SCCA(list(gene = expr, image = image), penalty_x=0.05, penalty_z=0.15, K=ncol(image))
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("GTEx/CAE3_discrim_PCA/CCA_ls_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("GTEx/CAE3_discrim_PCA/CCA_ls_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("GTEx/CAE3_discrim_PCA/CCA_ls_coefs.txt")


##### trying other discrim variants
image <- load_image_data("GTEx/CAE3_discrim3_PCA/CCA_input_image.txt")
expr <- load_expr_data("GTEx/CAE3_discrim3_PCA/CCA_input_expr.txt", "GTEx/CAE3_discrim3_PCA/CCA_genes.txt")

cca <- SCCA(list(gene = expr, image = image), penalty_x=0.05, penalty_z=0.15, K=9)#ncol(image))
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("GTEx/CAE3_discrim3_PCA/CCA_ls_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("GTEx/CAE3_discrim3_PCA/CCA_ls_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("GTEx/CAE3_discrim3_PCA/CCA_ls_coefs.txt")


image <- load_image_data("GTEx/CAE3_discrim4_PCA/CCA_input_image.txt")
expr <- load_expr_data("GTEx/CAE3_discrim4_PCA/CCA_input_expr.txt", "GTEx/CAE3_discrim4_PCA/CCA_genes.txt")

cca <- SCCA(list(gene = expr, image = image), penalty_x=0.05, penalty_z=0.15, K=9)#ncol(image))
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("GTEx/CAE3_discrim4_PCA/CCA_ls_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("GTEx/CAE3_discrim4_PCA/CCA_ls_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("GTEx/CAE3_discrim4_PCA/CCA_ls_coefs.txt")
#####










# ##### tune sparseness parameters
# 
# params <- bind_rows(lapply(seq(from=0, to=0.05, length.out=6), function(i) {
#   bind_rows(lapply(seq(from=0, to=0.2, length.out=6), function(j) {
#     print(paste(i, "and", j))
#     cca <- CCA(expr, image, penaltyx=i, penaltyz=j, K=9)
#     data.frame(expr_sparseness = i, 
#                image_sparseness = j, 
#                k = 1:ncol(cca$u), 
#                expr_selected = colSums(cca$u != 0),
#                image_selected = colSums(cca$v != 0))
#   }))
# }))
# 
# ggplot(params, aes(x=expr_sparseness, y=image_sparseness, fill=expr_selected, label=expr_selected)) + 
#   geom_tile() + 
#   geom_text(color='white') + 
#   facet_wrap(~ k, ncol=3) + 
#   xlab("expression feature sparseness") + 
#   ylab("image feature sparseness") + 
#   ggtitle("number of selected genes")
# ggsave(paste(dataset, "/CCA/params.expr.pdf", sep=''))
# 
# ggplot(params, aes(x=expr_sparseness, y=image_sparseness, fill=image_selected, label=image_selected)) + 
#   geom_tile() + 
#   geom_text(color='white') +
#   facet_wrap(~ k, ncol=3) + 
#   xlab("expression feature sparseness") + 
#   ylab("image feature sparseness") + 
#   ggtitle("number of selected image features")
# ggsave(paste(dataset, "/results/params.image.pdf", sep=''))
# 
# # now for LDA image features, assuming image penalty will be 0 to always select 1 of the 2 LDA image features
# 
# params_LDA <- bind_rows(lapply(seq(from=0, to=0.05, length.out=20), function(i) {
#   cca <- CCA(expr, image_LDA, penaltyx=i, penaltyz=0, K=2)
#   data.frame(expr_sparseness = i, 
#              image_sparseness = j, 
#              k = 1:ncol(cca$u), 
#              expr_selected = colSums(cca$u != 0),
#              image_selected = colSums(cca$v != 0))
# }))
# 
# ggplot(params_LDA, aes(x=expr_sparseness, y=image_sparseness, fill=expr_selected, label=expr_selected)) + 
#   geom_tile() + 
#   geom_text(color='white') + 
#   facet_wrap(~ k, ncol=1) + 
#   xlab("expression feature sparseness") + 
#   ylab("image feature sparseness") + 
#   ggtitle("number of selected genes")
# ggsave(paste(dataset, "/CCA/params-LDA.expr.pdf", sep=''), width=9, height=3)
# 
# ggplot(params_LDA, aes(x=expr_sparseness, y=image_sparseness, fill=image_selected, label=image_selected)) + 
#   geom_tile() + 
#   geom_text(color='white') +
#   facet_wrap(~ k, ncol=1) + 
#   xlab("expression feature sparseness") + 
#   ylab("image feature sparseness") + 
#   ggtitle("number of selected image features")
# ggsave(paste(dataset, "/CCA/params-LDA.image.pdf", sep=''), width=9, height=3)


##### shuffle data and see how SCCA correlation is affected

data_mods <- list('Original' = function(m) m,
                  'Shuffle Samples' = function(m) m[sample(nrow(m)), ],
                  'Shuffle Within Genes' = function(m) apply(m, 2, sample),
                  'Normal Random' = function(m) matrix(rnorm(nrow(m) * ncol(m)), nrow=nrow(m), ncol=ncol(m)))

cca_mod <- crossing(mod = factor(names(data_mods), levels=names(data_mods)),
                    iter = 1:10) %>%
  group_by(mod, iter) %>%
  do({
    cca <- SCCA(list(gene = data_mods[[.$mod]](expr), image = image), penalty_x=0.05, penalty_z=0.15, K=1)
    SCCA_cors(cca) %>%
      rowwise() %>%
      mutate(cor = cor(cca$CCA_var1[, k], cca$CCA_var2[, k])) %>%
      ungroup()
  }) %>%
  ungroup

write_tsv(cca_mod, "manuscript/CCA_expr_rand.txt")

cca_mod <- read_tsv("manuscript/CCA_expr_rand.txt", col_types='ciidd') %>%
  mutate(mod = mod %>%
           fct_recode('Original' = 'original',
                      'Shuffle samples' = 'shuffle samples',
                      'Shuffle within genes' = 'shuffle within genes',
                      'Normal random' = 'normal random') %>%
           fct_relevel('Original',
                       'Shuffle samples',
                       'Shuffle within genes',
                       'Normal random'))

cca_mod %>%
  mutate(mod = fct_rev(mod)) %>%
  ggplot(aes(x=mod, y=d)) +
  geom_boxplot() +
  xlab("Input randomization") +
  ylab("Dot product of first pair of CCA variables") +
  ylim(0, NA) +
  coord_flip()
ggsave("manuscript/CCA_expr_rand.pdf", width=6, height=2)

cca_mod %>%
  group_by(mod) %>%
  summarise(d_mean = mean(d),
            d_stdev = sd(d),
            cor_mean = mean(cor),
            cor_stdev = sd(cor))
