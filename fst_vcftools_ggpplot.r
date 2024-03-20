library(dplyr)
library(ggplot2)
options(scipen = 999)

rm(list = ls())
setwd("/mnt/storage12/emma/PR_fst/fst_output/")

#a<- read.delim("American_Samoa_LD.hap.ld", stringsAsFactors = FALSE)
#b<- read.delim("Brazil_LD.hap.ld", stringsAsFactors = FALSE)

files <- list.files(pattern = "*weir.fst")

## remove files comparing to self
files_remove <- c()

for (i in 1:length(files)){
    if(unlist(strsplit(unlist(strsplit(files[i], "\\."))[1], "_vs_"))[1] == 
                unlist(strsplit(unlist(strsplit(files[i], "\\."))[1], "_vs_"))[2]){
                    files_remove <- c(files_remove, i)
                }
    if(unlist(strsplit(unlist(strsplit(files[i], "\\."))[1], "_vs_"))[2] == 
                unlist(strsplit(unlist(strsplit(files[i], "\\."))[1], "_vs_"))[1]){
                    files_remove <- c(files_remove, i)
}}

files2 <- files[-files_remove]
length(files_remove)
files <- files2

countries <- c()
dataframes <- c()

# loop through the files
for (i in 1:length(files)){
  # save the clean filename as a char var because it will be called again
  #fnClean = substr(files[i], start = 1, stop = nchar(files[i])-10)
  fnClean = unlist(strsplit(files[i], "\\."))[1]
  dataframes <- c(dataframes, fnClean)
  # create a cmd as a string to be parsed and evaluated on-the-fly
  # the point here is that you can use the 'fnClean' var in the string
  # without knowing what it is - assign is expecting an inline string
  # ...not a string saved as a var, so it can't be used in this case
  loadFileCMD = paste0(fnClean,' = read.delim(files[i], stringsAsFactors = FALSE)')
  print(loadFileCMD) # check the cmd
  eval(parse(text=loadFileCMD))
  
  # create another string command to be evaluated to insert the file name
  # to the 'FileFrom' field
  countries <- c(countries, fnClean)
  if(nrow(eval(parse(text=loadFileCMD)) >0)){
    country_name <- unlist(strsplit(files[i], "\\."))[1]
    addFnCMD = paste0(fnClean,'$FileFrom = "', country_name, '"')
    print(addFnCMD) # check the cmd
    eval(parse(text=addFnCMD)) }
}

# check added column for country
head(Burkina_Faso_vs_Kenya)

## combine all dataframes together - now they have country column
dataframes <- Filter(function(x) is(x, "data.frame"), mget(ls()))
all <- do.call(rbind, dataframes)
head(all)

### write to csv so can load in from here ###

#write.csv(all, paste0("all_countries_fst_", Sys.Date(), ".csv"))
#df <- read.delim("/mnt/storage12/emma/PR_fst/fst_output/Burkina_Faso_vs_Kenya.weir.fst")
#library(dplyr)
#library(ggplot2)
#options(scipen = 999)
#rm(list = ls())
#setwd("/mnt/storage12/emma/PR_fst/fst_output/")
all <- read.csv("all_countries_fst_2024-03-18.csv", sep = ",", header = TRUE)

# Assuming dataframes is a list containing data frames
# Convert list elements to objects in the global environment
#list2env(dataframes, envir = .GlobalEnv)
# Remove data frames from the global environment
#rm(list = names(dataframes))


df_sub <- all

unique(df_sub$CHROM)
class(df_sub$CHROM)
df_sub <- df_sub[df_sub$CHROM %in% c(35107, 35108, 35109, 35159),]
df_sub$CHROM2 <- NA
df_sub$CHROM2[df_sub$CHROM == 35107] <- "Chromosome 1"
df_sub$CHROM2[df_sub$CHROM == 35108] <- "Chromosome 2"
df_sub$CHROM2[df_sub$CHROM == 35109] <- "Chromosome 3"
df_sub$CHROM2[df_sub$CHROM == 35159] <- "Mitochondria"

### fst above 0.8
df_sub$WEIR_AND_COCKERHAM_FST <- as.numeric(df_sub$WEIR_AND_COCKERHAM_FST)

## remove duplicate combos and summarise
summary_df <- df_sub %>%
  group_by(FileFrom) %>%
  summarize(
    Above_0.8 = sum(WEIR_AND_COCKERHAM_FST > 0.8, na.rm = TRUE),
    Below_0.8 = sum(WEIR_AND_COCKERHAM_FST < 0.8, na.rm = TRUE),
    total_sites = as.numeric(sum(!is.na(WEIR_AND_COCKERHAM_FST))),
    Equals_1 = sum(WEIR_AND_COCKERHAM_FST == 1, na.rm = TRUE),
    proportion_above_0.8 = (Above_0.8 / total_sites)*100,
    proportion_equals_1 = (Equals_1 /total_sites)*100
    
  )


summary_df %>% top_n(proportion_above_0.8, n = 10)
summary_df %>% slice_min(proportion_above_0.8, n = 10)



# Min-Max normalization
min_value <- min(summary_df$total_sites)
max_value <- max(summary_df$total_sites)
summary_df$Normalized_Value <- (summary_df$total_sites - min_value) / (max_value - min_value)

test <- summary_df[!duplicated(summary_df$Above_0.8),]
unique(test$FileFrom)
df_sub <- df_sub[df_sub$FileFrom %in% test$FileFrom,]
unique(df_sub$FileFrom)
summary_df <- summary_df[summary_df$FileFrom %in% test$FileFrom,]

meta <- read.csv("/mnt/storage12/emma/PR_fst/all_sample_country_metadata.csv")
head(meta)

for (i in 1:nrow(summary_df)){
    summary_df$country1[i] <- unlist(strsplit(summary_df$FileFrom[i], "_vs_"))[1]
    summary_df$country2[i] <- unlist(strsplit(summary_df$FileFrom[i], "_vs_"))[2] 
}

summary_df$Region_1 <- meta$Region[match(summary_df$country1, meta$Country)]
summary_df$Region_2 <- meta$Region[match(summary_df$country2, meta$Country)]

summary_df$region_match <- NA
summary_df$region_match[summary_df$Region_1 == summary_df$Region_2] <- "Same Region"
summary_df$region_match[summary_df$Region_1 != summary_df$Region_2] <- "Different Region"

summary_df$combination <- paste0(summary_df$country1, " vs ", summary_df$country2)





plot_both <- ggplot(data = summary_df) +
    geom_bar(aes(x = reorder(FileFrom, Above_0.8, desc), y = Above_0.8), stat = "identity", fill = "#40476D", alpha = 0.5) +
    geom_bar(aes(x = reorder(FileFrom, Above_0.8, desc), y = Equals_1), stat = "identity", fill = "#40476D", alpha = 1) +
    ylab("Number of positions with fst > 0.8") +
    xlab("") +
    theme_classic() +
    theme(axis.text.x = element_text(size = 20, angle = 90, hjust = 0.9, vjust = 0.5),
            axis.text.y = element_text(size = 20, hjust = 1, vjust = -0.5),
            axis.title = element_text(size = 20))

plot_both_prop <- ggplot(data = summary_df) +
    geom_bar(aes(x = reorder(combination, proportion_above_0.8, desc), y = proportion_above_0.8), stat = "identity", fill = "#40476D", alpha = 0.5) +
    geom_bar(aes(x = reorder(combination, proportion_above_0.8, desc), y = proportion_equals_1), stat = "identity", fill = "#40476D", alpha = 1) +
    ylab("Proportion of positions with fst > 0.8 or equaling 1") +
    xlab("") +
    theme_classic() +
    facet_grid(~region_match, scales = "free_x", space= "free") +
    theme(axis.text.x = element_text(size = 20, angle = 90, hjust = 0.9, vjust = 0.5),
            axis.text.y = element_text(size = 20, hjust = 1, vjust = -0.5),
            strip.text = element_text(size = 20),
            axis.title = element_text(size = 20))
plot_both_prop

plot_1 <- ggplot(data = summary_df) +
    geom_bar(aes(x = reorder(FileFrom, Equals_1, desc), y = Equals_1), stat = "identity", fill = "#40476D", alpha = 1) +
    ylab("Number of positions with fst = 1") +
    xlab("") +
    theme_classic() +
    theme(axis.text.x = element_text(size = 20, angle = 90, hjust = 0.9, vjust = 0.5),
            axis.text.y = element_text(size = 20, hjust = 1, vjust = -0.5),
            axis.title = element_text(size = 20))





plot_2 <- ggplot(data = summary_df) +
    geom_bar(aes(x = reorder(combination, proportion_above_0.8, desc), y = proportion_above_0.8), stat = "identity", fill = "#40476D", alpha = 1) +
    ylab("Proportion of sites with Fst > 0.8") +
    xlab("") +
    theme_classic() +
    facet_grid(~region_match, scales = "free_x", space= "free") +
    theme(axis.text.x = element_text(size = 20, angle = 90, hjust = 0.9, vjust = 0.5),
            axis.text.y = element_text(size = 20, hjust = 1, vjust = -0.5),
            strip.text = element_text(size = 20),
            axis.title = element_text(size = 20))
plot_2

ggplot(data = summary_df) +
    geom_tile(aes(x = FileFrom, fill = proportion_above_0.8), stat = "identity", alpha = 1) +
    ylab("Proportion of sites with Fst > 0.8") +
    xlab("") +
    theme_classic() +
    theme(axis.text.x = element_text(size = 20, angle = 90, hjust = 0.9, vjust = 0.5),
            axis.text.y = element_text(size = 20, hjust = 1, vjust = -0.5),
            axis.title = element_text(size = 20))


library(cowplot)
library(patchwork)
combined_plot <- plot_both + plot_1 + plot_layout(ncol = 2)

######################################
## checking difference in number of snps for each country comparison
######################################
ggplot()+
    geom_bar(data = summary_df, aes(x = reorder(FileFrom, total_sites, desc), y = total_sites), stat = "identity") +
    theme_minimal() +
    theme(axis.text = element_text(angle = 45))

ggplot()+
    geom_bar(data = summary_df, aes(x = reorder(FileFrom, Normalized_Value, desc), y = Normalized_Value), stat = "identity") +
    theme_minimal() +
    theme(axis.text = element_text(angle = 45))

## not much difference

######################################
### high positions only 
######################################
nrow(df_sub)

df_sub <- df_sub[!(is.na(df_sub$WEIR_AND_COCKERHAM_FST)),]
nrow(df_sub)

summary(df_sub$WEIR_AND_COCKERHAM_FST)

# mean per country?
a <- summary(df_sub$WEIR_AND_COCKERHAM_FST[df_sub$FileFrom == "Puerto_Rico_vs_Uganda"])

for (i in unique(df_sub$FileFrom)){

    df_i <- df_sub[df_sub$FileFrom == i,]
    print(df_i$FileFrom[1])
    print(summary(df_i$WEIR_AND_COCKERHAM_FST))
}




top_values <- df_sub %>%
  filter(WEIR_AND_COCKERHAM_FST > 0.8) %>%
  group_by(FileFrom) %>%
  top_n(10, WEIR_AND_COCKERHAM_FST) %>%
  arrange(FileFrom, desc(WEIR_AND_COCKERHAM_FST))
## nearly all 1 

## most common positions
mostcom <- df_sub %>% filter(WEIR_AND_COCKERHAM_FST > 0.8) %>%
                dplyr::group_by(CHROM, POS) %>% summarise(n = n())


## IR positions
#vgsc	3   315926360   316405639
#Rdl	2	41628484	41861946
#Ace-1	3	161486025	161486026
#Gste2	2	351633367	351634957

IR_vgsc <- df_sub[df_sub$POS >= 315926360 & df_sub$POS <= 316405639,]

IR_vgsc_top <- IR_vgsc %>% top_n(WEIR_AND_COCKERHAM_FST, n = 20)

IR_rdl <- df_sub[df_sub$CHROM == 035108 & df_sub$POS >= 41628484 & df_sub$POS <= 41861946,] 

IR_rdl_top <- IR_rdl %>% top_n(WEIR_AND_COCKERHAM_FST, n = 20)
summary(IR_rdl$WEIR_AND_COCKERHAM_FST)

#write.csv(summary_df, "high_fst_0.8.csv")

##### plots
library(ggplot2)

ggplot()+ 
    geom_point(data= df_sub[df_sub$FileFrom == "Burkina_Faso_vs_Kenya",], aes(x = POS, y = WEIR_AND_COCKERHAM_FST), alpha = 0.5) +
    geom_hline(yintercept = 0.8, linetype = "dashed") +
    ylim(c(0,1)) +
    ylab("Fst value") +
    xlab("Position") +
    ggtitle("Burkina Faso") +
    facet_wrap(~CHROM2, ncol = 4, scales = "free_x") +
    theme_minimal() +
    theme(axis.text = element_text(angle =45, size = 16),
            axis.title = element_text(size = 16),
            title = element_text(size = 18),
            strip.text = element_text(size = 16))

## find position with highest fst
k_k <- df_sub[df_sub$FileFrom == "Burkina_Faso_vs_Kenya",]
k_k <- k_k[complete.cases(k_k),]
k_k[k_k$WEIR_AND_COCKERHAM_FST == max(k_k$WEIR_AND_COCKERHAM_FST),]

vgsc_test <- df_sub[df_sub$FileFrom == "Burkina_Faso_vs_Kenya" & df_sub$CHROM == 35108 & (df_sub$POS >315926360 & df_sub$POS <316405639) ,]
vgsc_test <- df_sub[df_sub$CHROM == 35108 & (df_sub$POS >315926360 & df_sub$POS <316405639),]

top_values_vgsc <- vgsc_test %>%
  filter(WEIR_AND_COCKERHAM_FST > 0.8) %>%
  group_by(FileFrom) %>%
  top_n(10, WEIR_AND_COCKERHAM_FST) %>%
  arrange(FileFrom, desc(WEIR_AND_COCKERHAM_FST))



ggplot()+ 
    geom_point(data= vgsc_test[vgsc_test$FileFrom %in% top_values_vgsc$FileFrom,], aes(x = POS, y = WEIR_AND_COCKERHAM_FST), alpha = 0.5) +
    geom_hline(yintercept = 0.8, linetype = "dashed") +
    ylim(c(0.7,1)) +
    ylab("Fst value") +
    xlab("Position") +
    facet_wrap(~CHROM2 + FileFrom, ncol = 3, scales = "free_x") +
    theme_minimal() +
    theme(axis.text = element_text(angle =45, size = 10),
            axis.title = element_text(size = 10))


vgsc_test[vgsc_test$FileFrom %in% top_values_vgsc$FileFrom & vgsc_test$WEIR_AND_COCKERHAM_FST > 0.8,]
nrow(vgsc_test[vgsc_test$FileFrom %in% top_values_vgsc$FileFrom & vgsc_test$WEIR_AND_COCKERHAM_FST > 0.8,])
unique(vgsc_test$POS[vgsc_test$FileFrom %in% top_values_vgsc$FileFrom & vgsc_test$WEIR_AND_COCKERHAM_FST > 0.8])


## loop through each comparison and plot

pdf("output2.pdf",
    width = 15,
    height = 17)

plots <- list()
for (i in unique(df_sub$FileFrom)){

    country_df <- df_sub[df_sub$FileFrom == i,]

    plot[i] <- ggplot()+ 
    geom_point(data= country_df, aes(x = POS, y = WEIR_AND_COCKERHAM_FST), alpha = 0.5) +
    geom_hline(yintercept = 0.8, linetype = "dashed") +
    ylim(c(0,1)) +
    ylab("Fst value") +
    xlab("Position") +
    ggtitle(i) +
    facet_wrap(~CHROM2, ncol = 4, scales = "free_x") +
    theme_minimal() +
    theme(axis.text = element_text(angle =45, size = 16),
            axis.title = element_text(size = 18),
            title = element_text(size = 18),
            strip.text = element_text(size = 16))

    #ggsave(paste0("fst_plots/", i, ".tiff"), plot)
    print(plot)
    plots <- c(plots, plot[i])
}
graphics.off()



