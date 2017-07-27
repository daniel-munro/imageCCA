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

### CAE3 PCA

image <- load_image_data("BRCA/CCA_input_image.txt")
expr <- load_expr_data("BRCA/CCA_input_expr.txt", "BRCA/CCA_genes.txt")

cca <- SCCA(list(gene = expr, image = image), penalty_x=0.05, penalty_z=0.15, K=100)
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("BRCA/CCA_ls_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("BRCA/CCA_ls_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("BRCA/CCA_ls_coefs.txt")


### CAE3 discriminative PCA

image <- load_image_data("BRCA_discrim/CCA_input_image.txt")
expr <- load_expr_data("BRCA_discrim/CCA_input_expr.txt", "BRCA_discrim/CCA_genes.txt")

cca <- SCCA(list(gene = expr, image = image), penalty_x=0.05, penalty_z=0.15, K=ncol(image))
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("BRCA_discrim/CCA_ls_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("BRCA_discrim/CCA_ls_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("BRCA_discrim/CCA_ls_coefs.txt")


##### LGG

### CAE3 PCA

image <- load_image_data("LGG/CCA_input_image.txt")
expr <- load_expr_data("LGG/CCA_input_expr.txt", "LGG/CCA_genes.txt")

cca <- SCCA(list(gene = expr, image = image), penalty_x=0.05, penalty_z=0.15, K=100)
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("LGG/CCA_ls_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("LGG/CCA_ls_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("LGG/CCA_ls_coefs.txt")


### CAE3 discriminative PCA

image <- load_image_data("LGG_discrim/CCA_input_image.txt")
expr <- load_expr_data("LGG_discrim/CCA_input_expr.txt", "LGG_discrim/CCA_genes.txt")

cca <- SCCA(list(gene = expr, image = image), penalty_x=0.05, penalty_z=0.15, K=ncol(image))
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("LGG_discrim/CCA_ls_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("LGG_discrim/CCA_ls_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("LGG_discrim/CCA_ls_coefs.txt")


##### GTEx

### CAE3 PCA

image <- load_image_data("GTEx/CCA_input_image.txt")
expr <- load_expr_data("GTEx/CCA_input_expr.txt", "GTEx/CCA_genes.txt")

cca <- SCCA(list(gene = expr, image = image), penalty_x=0.05, penalty_z=0.15, K=100)
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("GTEx/CCA_ls_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("GTEx/CCA_ls_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("GTEx/CCA_ls_coefs.txt")


### CAE3 discriminative PCA

image <- load_image_data("GTEx_discrim/CCA_input_image.txt")
expr <- load_expr_data("GTEx_discrim/CCA_input_expr.txt", "GTEx_discrim/CCA_genes.txt")

cca <- SCCA(list(gene = expr, image = image), penalty_x=0.05, penalty_z=0.15, K=ncol(image))
cca$CCA_var1 %>% as_data_frame() %>% write_tsv("GTEx_discrim/CCA_ls_canonvar_expr.txt", col_names=F)
cca$CCA_var2 %>% as_data_frame() %>% write_tsv("GTEx_discrim/CCA_ls_canonvar_image.txt", col_names=F)
cca %>% SCCA_coefs() %>% write_tsv("GTEx_discrim/CCA_ls_coefs.txt")




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

write_tsv(cca_mod, "CCA_expr_rand.txt")

cca_mod <- read_tsv("CCA_expr_rand.txt", col_types='ciidd') %>%
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
ggsave("CCA_expr_rand.pdf", width=6, height=2)

cca_mod %>%
  group_by(mod) %>%
  summarise(d_mean = mean(d),
            d_stdev = sd(d),
            cor_mean = mean(cor),
            cor_stdev = sd(cor))
