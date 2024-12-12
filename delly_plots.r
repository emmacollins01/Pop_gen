# Variants in genes of interest that PASS only 
library(data.table)



# Specify the path to your file
file_path <- "/mnt/storage11/emma/raw_data_WGS_PR/delly/June24/merged.sample_filt.BNDfiltered.pass.ann.of_interest.vcf"
file_path <- "/mnt/storage11/emma/raw_data_WGS_PR/delly/June24/filtered.vcf.gz"
file_path <- "/mnt/storage11/emma/raw_data_WGS_PR/delly/June24/filtered_of_interest.vcf.gz"
file_path <- "/mnt/storage12/emma/delly/June24/merged.sample_filt.BNDfiltered.pass.of_interest_expanded.vcf"

# Read the file into R, line by line
lines <- readLines(file_path)

# Filter out lines that start with '##' or '#'
filtered_lines <- lines[!grepl("^#", lines)]

# Write the filtered lines to a temporary file
temp_file <- tempfile()
writeLines(filtered_lines, temp_file)

# Read the data into a data.table
del <- fread(temp_file, sep = "\t", header = FALSE)

# Delete the temporary file
unlink(temp_file)

# Print the data table
head(del)

colnames(del) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","BAY_F1_Out","BAY_F2_Out", "BAY_F3_Out", "BAY_F4_Out", "BAY_M1_Out", "BAY_M2_Out", "BAY_M3_Out", "BAY_M4_Out", "CUL_F1_Out", "CUL_M1_Out", "DOR_F1_Out", "DOR_F2_Out", "DOR_F3_Out", "DOR_M1_Out", "DOR_M2_Out", "GUA_F1_Out", "GUA_F2_Out", "GUA_F4_Out", "GUA_M1_Out", "GUA_M4_Out", "GUA_M5_Out", "PON_M1_Out", "PON_M2_Out", "PON_M3_Out", "PON_M4_Out", "SJU_F1_Out", "SJU_F3_Out", "SJU_F4_Out", "SJU_F5_Out", "SJU_M1_Out", "SJU_M3_Out", "SJU_M4_Out", "SJU_M5_Out")

del <- del[del$CHROM %in% c("NC_035107.1", "NC_035108.1", "NC_035109.1"),]

table(del$CHROM)
table(del$ALT)

dels <- del[grep("SVTYPE=DEL", del$INFO),]
invs <- del[grep("SVTYPE=INV", del$INFO),]
dups <- del[grep("SVTYPE=DUP", del$INFO),]
inss <- del[grep("SVTYPE=INS", del$INFO),]
nrow(dups)

# make SV type variable
del$SVtype <- NA
for (i in 1:nrow(del)){
    if(grepl("SVTYPE=DEL", del$INFO[i])){
        del$SVtype[i] <- "Deletion" 
    }
    if(grepl("SVTYPE=INV", del$INFO[i])){
        del$SVtype[i] <- "Inversion" 
    }
    if(grepl("SVTYPE=INS", del$INFO[i])){
        del$SVtype[i] <- "Insertion" 
    }
    if(grepl("SVTYPE=DUP", del$INFO[i])){
        del$SVtype[i] <- "Duplication" 
    }
}

unique(del$SVtype)

## look at inversions
inv <- del[del$SVtype == "Inversion",]
nrow(inv)

## extact ANN column
library(stringr)

for (i in 1:nrow(del)){
    del$ANN[i] <- str_extract(del$INFO[i], "ANN=[^;]+")
}

for (i in 1:nrow(del)){
    del$csq_type[i] <- unlist(str_split(del$ANN[i], "\\|"))[2]
    del$csq_impact[i] <- unlist(str_split(del$ANN[i], "\\|"))[3]
    del$gene[i] <- str_extract(del$ANN[i], "LOC[0-9]+")
}

unique(del$csq_type)
unique(del$csq_impact)
length(unique(del$POS))
length(unique(del$gene))

table(del$CHROM)
table(del$csq_type)
sort(table(del$csq_type))

library(dplyr)
top_csq <- del %>% group_by(csq_type) %>% summarise(n=n()) %>% top_n(n = 10)

for (i in 1:nrow(del)){
    if(del$csq_type[i] %in% top_csq$csq_type){
        del$csq_simple[i] <- del$csq_type[i]
    }
    else{
        del$csq_simple[i] <- "Other"
    } 
}
unique(del$csq_simple)
table(del$csq_simple)

table(del$SVtype)
length(unique(paste(del$POS, del$CHROM)))

########################
## make summary plots ##
########################

library(dplyr)
library(ggplot2)

pie <- del %>% group_by(SVtype) %>% summarise(n= n()) %>% mutate(perc = n / sum(n) * 100)

library(ggplot2)

ggplot(pie, aes(x=1, y=perc, fill=factor(SVtype))) +
    geom_bar(stat="identity", width=1) +
    scale_fill_manual(values = c("#40476D", "#749cb4",  "#A6C48A", "#51A3A3","#EE964B"), 
                      breaks = unique(pie$SVtype)) +
    coord_polar("y", start=0) +
    theme_void() + # Optional: to clean up the plot
    theme(legend.position = "right") # Optional: to position the legend

ggplot() +
    geom_bar(data = pie, aes(x = SVtype, y = n, fill = factor(SVtype)), stat = "identity") +
    scale_fill_manual(values = c("#40476D", "#749cb4",  "#A6C48A", "#51A3A3","#EE964B"), breaks = unique(pie$SVtype)) +
    #coord_flip() +
    ylab("Frequency of Structural Variants") +
    xlab("") +
    #xlab("Structural Variant Type") +
    theme_minimal() +    
    theme(legend.position = "None",
            axis.text = element_text(size = 30, angle = 35),
            axis.title = element_text(size = 32),
            legend.text = element_text(size = 30)) # Optional: to position the legend





pie2 <- del %>% group_by(csq_simple) %>% summarise(n= n()) %>% mutate(perc = n / sum(n) * 100)

library(ggplot2)
library(unikn)

cols <- unikn::usecol(c("goldenrod2", "coral2", "darkmagenta", "darkgreen", "darkorchid3", "dodgerblue"), n = 11)
cols <- unikn::usecol(c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#ffffbf","#e6f598", "#abdda4","#66c2a5","#3288bd", "darkorchid3", "#5e4fa2"),n = 13)

pie2$csq_simple[pie2$csq_simple == "splice_acceptor_variant&splice_donor_variant&splice_region_variant&5_prime_UTR_variant&intron_variant"] <- "combination"

plot <- ggplot(pie2, aes(x=1, y=perc, fill=factor(csq_simple))) +
    geom_bar(stat="identity", width=1) +
    #scale_fill_manual(values = cols, breaks = unique(pie2$csq_simple), name = "SV type")+    
    #scale_fill_manual(values = c("#9E0142", "#D0384D", "#EE6445", "#FA9C58", 
                   # "#FDCD7B", "#FEF0A7", "#F3FAAD","#D0EC9C", "#98D5A4", 
                    #"#5CB7A9", "#3682BA", "darkorchid3", "#5E4FA2"), 
                    #breaks = c("Other", "downstream_gene_variant", "duplication", 
                    #"exon_loss_variant", "feature_ablation","gene_fusion","intragenic_variant",  
                    #"intron_variant",  "inversion", "combination","transcript_ablation", 
                    #"upstream_gene_variant"), name = "SV type") +
    scale_fill_manual(values = c("#9E0142", "#D0384D", "#EE6445", "#FA9C58",
                   "#FDCD7B", "#FEF0A7", "#F3FAAD","#D0EC9C", "#98D5A4",
                    "#5CB7A9", "#3682BA", "darkorchid3", "#5E4FA2"),
                    breaks = c("Other", "downstream_gene_variant", "duplication",  # nolint
                    "exon_loss_variant", "feature_ablation","gene_fusion","intragenic_variant",   # nolint
                    "intron_variant",  "inversion", "bidirectional_gene_fusion",  # nolint
                    "upstream_gene_variant", "3_prime_UTR_variant", "5_prime_UTR_variant"), name = "SV type") + # nolint
    coord_polar("y", start=0) +
    theme_void() + # Optional: to clean up the plot
    theme(legend.position = "None",legend.box="vertical", legend.margin=margin(),
            legend.text = element_text(size= 24),
            legend.title = element_text(size = 24)) # Optional: to position the legend

plot + guides(fill=guide_legend(byrow=TRUE))

plot2<- ggplot(pie2[pie2$csq_simple != "intron_variant",], aes(x=1, y=perc, fill=factor(csq_simple))) +
    geom_bar(stat="identity", width=1) +
    #scale_fill_manual(values = c("#9E0142", "#D0384D", "#EE6445", "#FA9C58", "#FDCD7B", "#FEF0A7", "#F3FAAD","#D0EC9C", "#98D5A4", "#5CB7A9", "#3682BA", "#5E4FA2"), breaks = c("Other", "downstream_gene_variant", "duplication", "exon_loss_variant", "feature_ablation","gene_fusion","intragenic_variant",  "intron_variant",  "inversion", "combination","transcript_ablation", "upstream_gene_variant"), name = "SV type") +
    scale_fill_manual(values = c("#9E0142", "#D0384D", "#EE6445", "#FA9C58", 
                   "#FDCD7B", "#FEF0A7", "#F3FAAD","#D0EC9C", "#98D5A4", 
                    "#5CB7A9", "#3682BA", "darkorchid3", "#5E4FA2"), 
                    breaks = c("Other", "downstream_gene_variant", "duplication",  # nolint
                    "exon_loss_variant", "feature_ablation","gene_fusion","intragenic_variant",   # nolint
                    "intron_variant",  "inversion", "bidirectional_gene_fusion",  # nolint
                    "upstream_gene_variant", "3_prime_UTR_variant", "5_prime_UTR_variant"), name = "SV type") + # nolint
    #scale_fill_manual(values = cols, breaks = unique(pie2$csq_simple), name = "SV type")+    
    coord_polar("y", start=0) +
    theme_void() + # Optional: to clean up the plot
    theme(legend.position = "None",legend.box="vertical", legend.margin=margin(),
            legend.text = element_text(size= 24),
            legend.title = element_text(size = 24)) 


png(filename = "/mnt/storage12/emma/delly/June24/SVs_delly_all.png",
    width = 3,
    height = 3,
    units = "in",
    res=1200,)
plot
graphics.off()

png(filename = "/mnt/storage12/emma/delly/June24/SVs_delly_no_intron.png",
    width = 6,
    height = 3,
    units = "in",
    res=1200,)
plot2
graphics.off()


plot3<- ggplot(pie2[pie2$csq_simple != "intron_variant",], aes(x=1, y=perc, fill=factor(csq_simple))) +
    geom_bar(stat="identity", width=1) +
    #scale_fill_manual(values = c("#9E0142", "#D0384D", "#EE6445", "#FA9C58", "#FDCD7B", "#FEF0A7", "#F3FAAD","#D0EC9C", "#98D5A4", "#5CB7A9", "#3682BA", "#5E4FA2"), breaks = c("Other", "downstream_gene_variant", "duplication", "exon_loss_variant", "feature_ablation","gene_fusion","intragenic_variant",  "intron_variant",  "inversion", "combination","transcript_ablation", "upstream_gene_variant"), name = "SV type") +
    scale_fill_manual(values = c("#9E0142", "#D0384D", "#EE6445", "#FA9C58", 
                   "#FDCD7B", "#FEF0A7", "#F3FAAD","#D0EC9C", "#98D5A4", 
                    "#5CB7A9", "#3682BA", "darkorchid3", "#5E4FA2"), 
                    breaks = c("Other", "downstream_gene_variant", "duplication",  # nolint
                    "exon_loss_variant", "feature_ablation","gene_fusion","intragenic_variant",   # nolint
                    "intron_variant",  "inversion", "bidirectional_gene_fusion",  # nolint
                    "upstream_gene_variant", "3_prime_UTR_variant", "5_prime_UTR_variant"), name = "SV type") + # nolint
    #scale_fill_manual(values = cols, breaks = unique(pie2$csq_simple), name = "SV type")+    
    coord_polar("y", start=0) +
    theme_void() + # Optional: to clean up the plot
    theme(legend.position = "left",legend.box="vertical", legend.margin=margin(),
            legend.text = element_text(size= 12),
            legend.title = element_text(size = 12)) 
png(filename = "/mnt/storage12/emma/delly/June24/pie_legend_only.png",
    width = 6,
    height = 4,
    units = "in",
    res=1200,)
plot3
graphics.off()

######################
## most common genes??
######################


mode_gene <- del %>% group_by(CHROM, gene) %>% summarise(n= n()) 
mode_gene_top <- del %>% group_by(CHROM, gene) %>% summarise(n= n()) %>% top_n(n= 15)
nrow(mode_gene_top)

mode_gene_top <- mode_gene_top[mode_gene_top$n > 15,]
mode_gene_top <- mode_gene_top[order(mode_gene_top$n, decreasing = TRUE),]
print(mode_gene_top, n = 15)

gene_list <- del %>% group_by(CHROM, gene) %>% summarise(n= n()) %>% top_n(n= 50)
nrow(gene_list)
gene_list2 <- gene_list$gene
write.csv(gene_list2, "/mnt/storage12/emma/delly/June24/gene_list_delly.csv", row.names = FALSE)

test <- del[del$gene == "LOC5572215",]

for (j in 1:nrow(del)){
    gene <- unique(unlist(str_extract_all(del$ANN[j], "LOC[0-9]+")))
    product <- unique(unlist(str_extract_all(del$SNN[j], "product=[^;]+")))
    del$gene[j] <- gene
    del$product[j] <- product
}
del$ANN[i] <- str_extract(del$INFO[1i], "ANN=[^;]+")

## IR genes

vgsc <- del[grep("LOC5567355", del$ANN),]
vgsc <- del[del$gene == "LOC5567355",]
nrow(vgsc)
table(vgsc$csq_type, vgsc$csq_impact)

rdl <- del[grep("LOC5570466", del$ANN),]
nrow(rdl)
rdl <- del[del$gene == "LOC5570466",]
nrow(rdl)
table(rdl$csq_type, rdl$csq_impact)

ace <- del[grep("LOC5578456", del$ANN),]
ace <- del[del$gene == "LOC5578456",]
nrow(ace)
table(ace$csq_type, ace$csq_impact)
table(ace$csq_impact)
ace_high <- ace[grepl("HIGH", ace$ANN),]
extract <- ace_high[grepl("CE", ace_high$ANN),]

gste <- del[grep("LOC110676855", del$ANN),]
table(gste$csq_type)

## look at high impact SVs

high <- del[grepl("HIGH", del$ANN),]
unique(high$csq_simple)
length(unique(paste(high$POS, high$CHROM)))
mod <- del[grepl("MODIFIER", del$ANN),]
length(unique(paste(mod$POS, mod$CHROM)))
med <- del[grepl("MODERATE", del$ANN),]
length(unique(paste(med$POS, med$CHROM)))
low <- del[grepl("LOW", del$ANN),]
length(unique(paste(low$POS, low$CHROM)))

high_top <- high %>% group_by(CHROM, POS, gene) %>% summarise(n = n())
high_top_top <- high_top %>% top_n(n, n = 20)
summary(high_top$n)

high_top2 <- high %>% group_by(CHROM, gene) %>% summarise(n = n())
summary(high_top2$n)
high_top2_top <- high_top2 %>% top_n(n, n = 10)
print(high_top2_top, n = 20)

### read in text file to look at frequency of SVs in samples