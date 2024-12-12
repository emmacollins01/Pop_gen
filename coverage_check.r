
coverage <- read.csv("/mnt/storage12/emma/PR_coverage/output_file", header = FALSE, sep = " ")

head(coverage)
colnames(coverage) <- c("sample", "coverage")
nrow(coverage)

samples_to_remove <- c("SRR11006622", "SRR11006875", "SRR9959052","SRR9986057")
coverage <- coverage[!(coverage$sample %in% samples_to_remove),]

metadata <- read.csv("/mnt/storage12/emma/PR_combine/all_sample_country_metadata.tsv", sep = "\t")

head(metadata)

coverage$country <- metadata$Country[match(coverage$sample, metadata$ID)]

coverage_sum <- coverage %>% group_by(country) %>% summarise(mean = mean(coverage))
