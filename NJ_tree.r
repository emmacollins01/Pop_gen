## MAKE NJ TREE

#### NJ TREE ####
workdir <- "/mnt/storage12/emma/NJ_tree_files" # Working directory with plink files
prefix <- "PR_WGS" # Prefix for plink files
metadata <- "/mnt/storage12/emma/PR_combine/all_sample_country_metadata.csv" # File path to metadata

#### DIST#
## run PLINK to make the distance matrix
## plink --vcf <vcf-file> --distance square --double-id --allow-extra-chr --out <prefix>

dist <- read.table(file.path(workdir, paste0(prefix, ".dist")), header = FALSE)
id <- read.table(file.path(workdir, paste0(prefix, ".dist.id")))
met <- read.table(metadata, sep = ",", stringsAsFactors = FALSE, header = TRUE)

desc <- id %>% left_join(met, by = c("V1" = "ID"))

dist_m <- as.matrix(dist)

# Export dist_m to .newick to make neighbour joining tree
tree <- nj(dist_m)
write.tree(phy = tree, file = file.path(workdir, paste0(prefix, ".newick.tree")))

require("ape")

tree_dat<-nj(dist_dat)
write.tree(tree_dat, file="/mnt/storage12/emma/NJ_tree/PR_WGS.tree")

## open in iTOL