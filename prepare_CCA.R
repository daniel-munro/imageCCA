library(tidyverse)
library(stringr)

prepare_CCA_files <- function(expr, images, samples, expr_outfile, image_outfile, gene_outfile, image_name_outfile) {
  ### only use samples with expression and image data ###
  samples <- samples %>%
    filter(sample %in% rownames(expr),
           image %in% rownames(images))

  images <- images[samples$image, ]
  expr <- expr[samples$sample, ]
  
  # required for CCA:
  images <- images[, sapply(1:ncol(images), function(x) sd(images[, x]) > 0)]
  expr <- expr[, sapply(1:ncol(expr), function(x) sd(expr[, x]) > 0)]
  
  # log and scale expression data, scale image data
  expr <- log(expr + 1) %>% scale()
  images <- scale(images)

  ### save gene and image file order for analyzing CCA output ###
  write_lines(colnames(expr), gene_outfile)
  write_lines(rownames(images), image_name_outfile)
  
  ### save image and expression matrices ###
  expr %>% as_data_frame %>% write_tsv(expr_outfile, col_names=F)
  images %>% as_data_frame %>% write_tsv(image_outfile, col_names=F)
}

PCA_whiten <- function(x) {
  pca <- prcomp(x)$x
  sapply(1:ncol(pca), function(c) pca[, c] / sd(pca[, c]))
}


### BRCA

samples <- read_tsv("BRCA/data/images/image_files.txt", col_types='cccccc')

expr <- read.table("BRCA/data/expression/RNA-Seq_BRCA.txt") %>% as.matrix()

images_SIFT <- read.csv("BRCA/SIFT/repkmeans.csv", header=F, row.names=1) %>% as.matrix()
prepare_CCA_files(
  expr = expr,
  images = images_SIFT,
  samples = samples,
  expr_outfile = "BRCA/SIFT/CCA_input_expr.txt",
  image_outfile = "BRCA/SIFT/CCA_input_image.txt",
  gene_outfile = "BRCA/SIFT/CCA_genes.txt",
  image_name_outfile = "BRCA/SIFT/CCA_images.txt"
)

images_CAE <- read.csv("BRCA/CAE/rep.csv", header=F, row.names=1) %>% as.matrix()
prepare_CCA_files(
  expr = expr,
  images = images_CAE,
  samples = samples,
  expr_outfile = "BRCA/CAE/CCA_input_expr.txt",
  image_outfile = "BRCA/CAE/CCA_input_image.txt",
  gene_outfile = "BRCA/CAE/CCA_genes.txt",
  image_name_outfile = "BRCA/CAE/CCA_images.txt"
)

images_CAE_PCA <- read.csv("BRCA/CAE_PCA/repPCA.csv", header=F, row.names=1) %>% as.matrix()
prepare_CCA_files(
  expr = expr,
  images = images_CAE_PCA,
  samples = samples,
  expr_outfile = "BRCA/CAE_PCA/CCA_input_expr.txt",
  image_outfile = "BRCA/CAE_PCA/CCA_input_image.txt",
  gene_outfile = "BRCA/CAE_PCA/CCA_genes.txt",
  image_name_outfile = "BRCA/CAE_PCA/CCA_images.txt"
)

images_CAE2_PCA <- read.csv("BRCA/CAE2_PCA/transformed.csv", header=F, row.names=1) %>% as.matrix()
prepare_CCA_files(
  expr = expr,
  images = images_CAE2_PCA,
  samples = samples,
  expr_outfile = "BRCA/CAE2_PCA/CCA_input_expr.txt",
  image_outfile = "BRCA/CAE2_PCA/CCA_input_image.txt",
  gene_outfile = "BRCA/CAE2_PCA/CCA_genes.txt",
  image_name_outfile = "BRCA/CAE2_PCA/CCA_images.txt"
)

images_CAE3 <- read.csv("BRCA/CAE3/cae_pca_brca_rep.csv", header=F, row.names=1) %>% as.matrix()
prepare_CCA_files(
  expr = expr,
  images = images_CAE3,
  samples = samples,
  expr_outfile = "BRCA/CAE3/CCA_input_expr.txt",
  image_outfile = "BRCA/CAE3/CCA_input_image.txt",
  gene_outfile = "BRCA/CAE3/CCA_genes.txt",
  image_name_outfile = "BRCA/CAE3/CCA_images.txt"
)

images_CAE3_PCA <- read.csv("BRCA/CAE3_PCA/cae_pca_brca_repPCA.csv", header=F, row.names=1) %>% as.matrix()# %>% .[, 1:10]
prepare_CCA_files(
  expr = expr,
  images = images_CAE3_PCA,
  samples = samples,
  expr_outfile = "BRCA/CAE3_PCA/CCA_input_expr.txt",
  image_outfile = "BRCA/CAE3_PCA/CCA_input_image.txt",
  gene_outfile = "BRCA/CAE3_PCA/CCA_genes.txt",
  image_name_outfile = "BRCA/CAE3_PCA/CCA_images.txt"
)

images_CAE3_discrim <- read.csv("BRCA/CAE3_discrim/BRCA_rep_discrim.csv", header=F, row.names=1) %>% as.matrix()
prepare_CCA_files(
  expr = expr,
  images = images_CAE3_discrim,
  samples = samples,
  expr_outfile = "BRCA/CAE3_discrim/CCA_input_expr.txt",
  image_outfile = "BRCA/CAE3_discrim/CCA_input_image.txt",
  gene_outfile = "BRCA/CAE3_discrim/CCA_genes.txt",
  image_name_outfile = "BRCA/CAE3_discrim/CCA_images.txt"
)

images_CAE3_discrim_PCA <- PCA_whiten(images_CAE3_discrim)
prepare_CCA_files(
  expr = expr,
  images = images_CAE3_discrim_PCA,
  samples = samples,
  expr_outfile = "BRCA/CAE3_discrim_PCA/CCA_input_expr.txt",
  image_outfile = "BRCA/CAE3_discrim_PCA/CCA_input_image.txt",
  gene_outfile = "BRCA/CAE3_discrim_PCA/CCA_genes.txt",
  image_name_outfile = "BRCA/CAE3_discrim_PCA/CCA_images.txt"
)

##### trying other discrim variants
images_CAE3_discrim4 <- read.csv("BRCA/CAE3_discrim4_PCA/BRCA_rep_discrim_4.csv", header=F, row.names=1) %>% as.matrix()
images_CAE3_discrim4_PCA <- PCA_whiten(images_CAE3_discrim4)
prepare_CCA_files(
  expr = expr,
  images = images_CAE3_discrim4_PCA,
  samples = samples,
  expr_outfile = "BRCA/CAE3_discrim4_PCA/CCA_input_expr.txt",
  image_outfile = "BRCA/CAE3_discrim4_PCA/CCA_input_image.txt",
  gene_outfile = "BRCA/CAE3_discrim4_PCA/CCA_genes.txt",
  image_name_outfile = "BRCA/CAE3_discrim4_PCA/CCA_images.txt"
)
#####

### LGG

samples <- read_tsv("LGG/data/images/image_files.txt", col_types='cccccc')

expr <- read.table("LGG/data/expression/RNA-Seq_LGG.txt") %>% as.matrix()

images_SIFT <- read.csv("LGG/SIFT/kmeansSet.csv", header=F, row.names=1) %>% as.matrix()
prepare_CCA_files(
  expr = expr,
  images = images_SIFT,
  samples = samples,
  expr_outfile = "LGG/SIFT/CCA_input_expr.txt",
  image_outfile = "LGG/SIFT/CCA_input_image.txt",
  gene_outfile = "LGG/SIFT/CCA_genes.txt",
  image_name_outfile = "LGG/SIFT/CCA_images.txt"
)

images_SIFT_PCA <- read.csv("LGG/SIFT_PCA/kmeansSetPCA.csv", header=F, row.names=1) %>% as.matrix()
prepare_CCA_files(
  expr = expr,
  images = images_SIFT_PCA,
  samples = samples,
  expr_outfile = "LGG/SIFT_PCA/CCA_input_expr.txt",
  image_outfile = "LGG/SIFT_PCA/CCA_input_image.txt",
  gene_outfile = "LGG/SIFT_PCA/CCA_genes.txt",
  image_name_outfile = "LGG/SIFT_PCA/CCA_images.txt"
)

images_CAE3 <- read.csv("LGG/CAE3/cae_pca_lgg_rep.csv", header=F, row.names=1) %>% as.matrix()
prepare_CCA_files(
  expr = expr,
  images = images_CAE3,
  samples = samples,
  expr_outfile = "LGG/CAE3/CCA_input_expr.txt",
  image_outfile = "LGG/CAE3/CCA_input_image.txt",
  gene_outfile = "LGG/CAE3/CCA_genes.txt",
  image_name_outfile = "LGG/CAE3/CCA_images.txt"
)

images_CAE3_PCA <- read.csv("LGG/CAE3_PCA/cae_pca_lgg_repPCA.csv", header=F, row.names=1) %>% as.matrix()# %>% .[, 1:10]
prepare_CCA_files(
  expr = expr,
  images = images_CAE3_PCA,
  samples = samples,
  expr_outfile = "LGG/CAE3_PCA/CCA_input_expr.txt",
  image_outfile = "LGG/CAE3_PCA/CCA_input_image.txt",
  gene_outfile = "LGG/CAE3_PCA/CCA_genes.txt",
  image_name_outfile = "LGG/CAE3_PCA/CCA_images.txt"
)

images_CAE3_discrim <- read.csv("LGG/CAE3_discrim/LGG_rep_discrim.csv", header=F, row.names=1) %>% as.matrix()
prepare_CCA_files(
  expr = expr,
  images = images_CAE3_discrim,
  samples = samples,
  expr_outfile = "LGG/CAE3_discrim/CCA_input_expr.txt",
  image_outfile = "LGG/CAE3_discrim/CCA_input_image.txt",
  gene_outfile = "LGG/CAE3_discrim/CCA_genes.txt",
  image_name_outfile = "LGG/CAE3_discrim/CCA_images.txt"
)

images_CAE3_discrim_PCA <- PCA_whiten(images_CAE3_discrim)
prepare_CCA_files(
  expr = expr,
  images = images_CAE3_discrim_PCA,
  samples = samples,
  expr_outfile = "LGG/CAE3_discrim_PCA/CCA_input_expr.txt",
  image_outfile = "LGG/CAE3_discrim_PCA/CCA_input_image.txt",
  gene_outfile = "LGG/CAE3_discrim_PCA/CCA_genes.txt",
  image_name_outfile = "LGG/CAE3_discrim_PCA/CCA_images.txt"
)


##### trying other discrim variants
images_CAE3_discrim2 <- read.csv("LGG/CAE3_discrim2_PCA/LGG_rep_discrim_2.csv", header=F, row.names=1) %>% as.matrix()
images_CAE3_discrim2_PCA <- PCA_whiten(images_CAE3_discrim2)
prepare_CCA_files(
  expr = expr,
  images = images_CAE3_discrim2_PCA,
  samples = samples,
  expr_outfile = "LGG/CAE3_discrim2_PCA/CCA_input_expr.txt",
  image_outfile = "LGG/CAE3_discrim2_PCA/CCA_input_image.txt",
  gene_outfile = "LGG/CAE3_discrim2_PCA/CCA_genes.txt",
  image_name_outfile = "LGG/CAE3_discrim2_PCA/CCA_images.txt"
)

images_CAE3_discrim3 <- read.csv("LGG/CAE3_discrim3_PCA/LGG_rep_discrim_3.csv", header=F, row.names=1) %>% as.matrix()
images_CAE3_discrim3_PCA <- PCA_whiten(images_CAE3_discrim3)
prepare_CCA_files(
  expr = expr,
  images = images_CAE3_discrim3_PCA,
  samples = samples,
  expr_outfile = "LGG/CAE3_discrim3_PCA/CCA_input_expr.txt",
  image_outfile = "LGG/CAE3_discrim3_PCA/CCA_input_image.txt",
  gene_outfile = "LGG/CAE3_discrim3_PCA/CCA_genes.txt",
  image_name_outfile = "LGG/CAE3_discrim3_PCA/CCA_images.txt"
)
#####

### GTEx

expr <- read.table("GTEx/data/expression/RNA-Seq_GTEx.txt") %>% as.matrix()

samples <- data_frame(sample = rownames(expr),
                      image = str_c(sample, '.jpg'))

images_CAE3 <- read.csv("GTEx/CAE3/cae_pca_gtex_rep.csv", header=F, row.names=1) %>% as.matrix()
prepare_CCA_files(
  expr = expr,
  images = images_CAE3,
  samples = samples,
  expr_outfile = "GTEx/CAE3/CCA_input_expr.txt",
  image_outfile = "GTEx/CAE3/CCA_input_image.txt",
  gene_outfile = "GTEx/CAE3/CCA_genes.txt",
  image_name_outfile = "GTEx/CAE3/CCA_images.txt"
)

images_CAE3_PCA <- read.csv("GTEx/CAE3_PCA/cae_pca_gtex_repPCA.csv", header=F, row.names=1) %>% as.matrix()# %>% .[, 1:10]
prepare_CCA_files(
  expr = expr,
  images = images_CAE3_PCA,
  samples = samples,
  expr_outfile = "GTEx/CAE3_PCA/CCA_input_expr.txt",
  image_outfile = "GTEx/CAE3_PCA/CCA_input_image.txt",
  gene_outfile = "GTEx/CAE3_PCA/CCA_genes.txt",
  image_name_outfile = "GTEx/CAE3_PCA/CCA_images.txt"
)

images_CAE3_discrim <- read.csv("GTEx/CAE3_discrim/GTEx_rep_discrim.csv", header=F, row.names=1) %>% as.matrix()
prepare_CCA_files(
  expr = expr,
  images = images_CAE3_discrim,
  samples = samples,
  expr_outfile = "GTEx/CAE3_discrim/CCA_input_expr.txt",
  image_outfile = "GTEx/CAE3_discrim/CCA_input_image.txt",
  gene_outfile = "GTEx/CAE3_discrim/CCA_genes.txt",
  image_name_outfile = "GTEx/CAE3_discrim/CCA_images.txt"
)

images_CAE3_discrim_PCA <- PCA_whiten(images_CAE3_discrim)
prepare_CCA_files(
  expr = expr,
  images = images_CAE3_discrim_PCA,
  samples = samples,
  expr_outfile = "GTEx/CAE3_discrim_PCA/CCA_input_expr.txt",
  image_outfile = "GTEx/CAE3_discrim_PCA/CCA_input_image.txt",
  gene_outfile = "GTEx/CAE3_discrim_PCA/CCA_genes.txt",
  image_name_outfile = "GTEx/CAE3_discrim_PCA/CCA_images.txt"
)

##### trying discrim variants
images_CAE3_discrim3 <- read.csv("GTEx/CAE3_discrim3_PCA/GTEx_rep_discrim_3.csv", header=F, row.names=1) %>% as.matrix()
images_CAE3_discrim3_PCA <- PCA_whiten(images_CAE3_discrim3)
prepare_CCA_files(
  expr = expr,
  images = images_CAE3_discrim3_PCA,
  samples = samples,
  expr_outfile = "GTEx/CAE3_discrim3_PCA/CCA_input_expr.txt",
  image_outfile = "GTEx/CAE3_discrim3_PCA/CCA_input_image.txt",
  gene_outfile = "GTEx/CAE3_discrim3_PCA/CCA_genes.txt",
  image_name_outfile = "GTEx/CAE3_discrim3_PCA/CCA_images.txt"
)

images_CAE3_discrim4 <- read.csv("GTEx/CAE3_discrim4_PCA/GTEx_rep_discrim_4.csv", header=F, row.names=1) %>% as.matrix()
images_CAE3_discrim4_PCA <- PCA_whiten(images_CAE3_discrim4)
prepare_CCA_files(
  expr = expr,
  images = images_CAE3_discrim4_PCA,
  samples = samples,
  expr_outfile = "GTEx/CAE3_discrim4_PCA/CCA_input_expr.txt",
  image_outfile = "GTEx/CAE3_discrim4_PCA/CCA_input_image.txt",
  gene_outfile = "GTEx/CAE3_discrim4_PCA/CCA_genes.txt",
  image_name_outfile = "GTEx/CAE3_discrim4_PCA/CCA_images.txt"
)
#####





# ### BCCA format (combined matrix, rows are features) ###
# cbind(images, expr_match) %>%
#   t() %>%
#   write.table(file="../data/BCCA_input.txt", sep='\t', quote=F, row.names=F, col.names=F)
