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

samples <- read_tsv("BRCA_data/images/image_files.txt", col_types='cccccc')

expr <- read.table("BRCA_data/expression/RNA-Seq_BRCA.txt") %>% as.matrix()

images <- read.csv("BRCA/CAE_output_BRCA_repPCA.csv", header=F, row.names=1) %>% as.matrix()
prepare_CCA_files(
  expr = expr,
  images = images,
  samples = samples,
  expr_outfile = "BRCA/CCA_input_expr.txt",
  image_outfile = "BRCA/CCA_input_image.txt",
  gene_outfile = "BRCA/CCA_genes.txt",
  image_name_outfile = "BRCA/CCA_images.txt"
)

images_discrim <- PCA_whiten(images_CAE3_discrim)
prepare_CCA_files(
  expr = expr,
  images = images_discrim,
  samples = samples,
  expr_outfile = "BRCA_discrim/CCA_input_expr.txt",
  image_outfile = "BRCA_discrim/CCA_input_image.txt",
  gene_outfile = "BRCA_discrim/CCA_genes.txt",
  image_name_outfile = "BRCA_discrim/CCA_images.txt"
)


### LGG

samples <- read_tsv("LGG_data/images/image_files.txt", col_types='cccccc')

expr <- read.table("LGG_data/expression/RNA-Seq_LGG.txt") %>% as.matrix()

images <- read.csv("LGG/CAE_output_LGG_repPCA.csv", header=F, row.names=1) %>% as.matrix()
prepare_CCA_files(
  expr = expr,
  images = images,
  samples = samples,
  expr_outfile = "LGG/CCA_input_expr.txt",
  image_outfile = "LGG/CCA_input_image.txt",
  gene_outfile = "LGG/CCA_genes.txt",
  image_name_outfile = "LGG/CCA_images.txt"
)

images_discrim <- PCA_whiten(images_CAE3_discrim)
prepare_CCA_files(
  expr = expr,
  images = images_discrim,
  samples = samples,
  expr_outfile = "LGG_discrim/CCA_input_expr.txt",
  image_outfile = "LGG_discrim/CCA_input_image.txt",
  gene_outfile = "LGG_discrim/CCA_genes.txt",
  image_name_outfile = "LGG_discrim/CCA_images.txt"
)


### GTEx

expr <- read.table("GTEx_data/expression/RNA-Seq_GTEx.txt") %>% as.matrix()

samples <- data_frame(sample = rownames(expr),
                      image = str_c(sample, '.jpg'))

images <- read.csv("GTEx/CAE_output_GTEx_repPCA.csv", header=F, row.names=1) %>% as.matrix()
prepare_CCA_files(
  expr = expr,
  images = images,
  samples = samples,
  expr_outfile = "GTEx/CCA_input_expr.txt",
  image_outfile = "GTEx/CCA_input_image.txt",
  gene_outfile = "GTEx/CCA_genes.txt",
  image_name_outfile = "GTEx/CCA_images.txt"
)

images_discrim <- PCA_whiten(images_CAE3_discrim)
prepare_CCA_files(
  expr = expr,
  images = images_discrim,
  samples = samples,
  expr_outfile = "GTEx_discrim/CCA_input_expr.txt",
  image_outfile = "GTEx_discrim/CCA_input_image.txt",
  gene_outfile = "GTEx_discrim/CCA_genes.txt",
  image_name_outfile = "GTEx_discrim/CCA_images.txt"
)
