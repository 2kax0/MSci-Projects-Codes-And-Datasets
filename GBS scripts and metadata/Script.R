###############################
# Describing the population structure of minor lineages of Group B Streptococcus and their potential risk post-vaccination
###############################

# load packages
library(tidyverse)
library(ggpubr)
library(ggalluvial)

complete <- read_csv(file = "metadata.csv")

##########
# raw number
##########

complete |>
  group_by(Lineage) |>
  tally() |>
  ggplot(aes(x = as.factor(Lineage), y = n)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), vjust = -1) +
  labs(
    title = "Number of isolates within each minor lineage",
    x = "Lineage",
    y = "Number of isolates"
  )

##########
# sub-region
##########

# read sub-region file
subregion <- read_csv(file = "subregion.csv")

Lineage_subregion <- left_join(complete, subregion, by = "Country") |>
  select(!n)

Lineage_subregion |>
  ggplot(aes(x = as.factor(Lineage), fill = Subregion)) +
  geom_bar() +
  scale_fill_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", 
                               "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", 
                               "#FF69B4")) +
  labs(
    x = "Lineage",
    y = "Number of isolates",
    title = "Subregion of origin for isolates from each minor lineage"
  )

##########
# CPS, host and AMR figure
##########

# read csv and manually arrange cps order
complete$cps.final <- factor(complete$cps.final, levels = c("Ia","Ib","II","III","IV","V","VI","VII","VIII","IX","NT"))

# barchart
bar <- complete |>
  group_by(Lineage) |>
  ggplot(aes(x = as.factor(Lineage), fill = cps.final)) +
  scale_fill_manual(values = c("#E69F00","#56B4E9","#0072B2","#F0E442","#009E73","#D55E00","#BCBDDC", "#FA9FB5","#BFFFEB","#808080")) +
  theme(
    plot.margin = unit(c(0.2, 0, -2, 0), "lines"),
    axis.ticks = element_blank()
  ) +
  geom_bar() +
  labs(
    x = "Lineage",
    y = "Number of isolates",
    fill = "Serotype",
    title = "Serotype, host and AMR gene distribution within minor lineages"
  )

### tile data prep
# vector storing genes to be analysed
gene <- c("aph(3'-III)", "aadE", "aac(6')-aph(2'')", "ant(6-Ia)", "cat(pc194)", "catQ", "ermT", "ermA", "ermB", "lnuC",
          "lsaC", "mefA", "msrD", "tetL", "tetM", "tetO", "tetS")

# empty dataframe to store loop output
forloop_output <- data.frame()


for (i in gene) {
  # table grouped by lineage and gene presence/absence, then pivot wider
  tile <- complete |>
    group_by(Lineage, !!sym((i))) |>
    tally() |>
    pivot_wider(names_from = !!sym((i)), values_from = n)
  # add a column beside previous table showing the total number of isolates in each lineage
  tile$total <- rowSums(tile[, !colnames(tile) %in% "Lineage"], na.rm = TRUE)
  # mutate columns showing the gene and percentage of isolates carrying the gene
  row_for_dataframe <- tile |>
    mutate(
      gene = i,
      proportion = pos/total
    ) |>
    select(Lineage, gene, proportion)
  # set NAs into 0
  row_for_dataframe[is.na(row_for_dataframe)] <- 0
  # rbind outputs for individual genes into a complete table for tile plot
  forloop_output <- rbind(row_for_dataframe, forloop_output)
}

# tile plot
tile_plot <- forloop_output |>
  ggplot(aes(x = as.factor(Lineage), y = gene, fill = proportion)) +
  geom_tile(colour = "black") +
  scale_fill_gradient(low = "white", high = "red") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    plot.margin = margin(t = 0)
  ) +
  labs(
    fill = "Proportion", 
    y = "AMR gene"
  )
plot(tile_plot)

#host plot
host_bar <- complete |>
  ggplot(aes(x = as.factor(Lineage), fill = Host_species)) +
  geom_bar() +
  scale_fill_manual(values = c("#8975ca", "#cb5582", "#71a659", "#ff7837", "#808080")) +
  theme(
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "lines")
  ) +
  labs(
    y = "Number of isolates",
    fill = "Host species",
    x = "Lineage"
  ) 

# ggarrange serotype, host and tile plot
library(ggpubr)

figure <- ggarrange(bar + rremove("xlab"), host_bar, tile_plot + rremove("xlab"),
                    labels = NULL,
                    ncol = 1, align = "hv") 
annotate_figure(figure, bottom = textGrob("Lineage"))
plot(figure)

##########
# ST table
##########

# group failed, NAs and uncertain alleles (?) and mismatches (*) into NTs
tidy_ST <- complete |>
  mutate(
    ST_tidy = ifelse(grepl("\\*", complete$ST), sub("^.*", "NF", complete$ST),
                     ifelse(grepl("\\?", complete$ST), sub("^.*", "NF", complete$ST),
                            ifelse(grepl("\\*\\?", complete$ST), sub("^.*", "NF", complete$ST), 
                                   ifelse(ST %in% "failed", "NA", 
                                          ifelse(is.na(ST), "NA", complete$ST
                                          )))))
  ) 

# STs in each lineage
st_count <- tidy_ST |>
  filter(!ST_tidy %in% "NA") |>
  group_by(Lineage, ST_tidy) |>
  summarise(
    n = n()
  )

# total number of isolates in each lineage
st_total <- tidy_ST |>
  filter(!ST_tidy %in% "NA") |>
  group_by(Lineage) |>
  summarise(
    total = n()
  )

# join
st_complete <- left_join(st_count, st_total, by = "Lineage")

# calculate percentage, then filter out the dominant ST in each lineage
st_percentage <- st_complete |>
  group_by(Lineage) |>
  mutate(
    percentage = round(n/total*100, 1)
  ) |>
  slice_max(n)

# add brackets showing percentage of isolates in each lineage that belongs to the dominant ST
st_percentage <- st_percentage |>
  mutate(
    Dominant_ST = paste0(ST_tidy, " (", percentage, "%)")
  ) |>
  select(Lineage, total, Dominant_ST)

# number of unique STs in each lineage
st_unique <- st_count |>
  group_by(Lineage) |>
  filter(!ST_tidy %in% c("NF", "NA")) |>
  summarise(
    Unique_ST = length(unique(ST_tidy))
  ) |>
  select(Lineage, Unique_ST)

# number of NTs in each lineage
NTs <- tidy_ST |>
  group_by(Lineage) |>
  filter(ST_tidy %in% "NF") |>
  summarise(
    NT_count = n()
  )

# joining tables and write CSV file
st_table <- left_join(st_percentage, st_unique, by = "Lineage")
st_table_nt <- left_join(st_table, NTs, by = "Lineage")
st_table_nt[is.na(st_table_nt)] <- 0
write_csv(st_table_nt, file = "st table.csv")

##########
# alpha-like proteins
##########

# vector of each ALP
alp_protein <- c("rib", "alpha", "alp2/3", "alp1")

# empty dataframe to store loop results
alp_forloop_output <- data.frame()

for (i in alp_protein) {
  # table grouped by lineage and ALP presence/absence, then pivot wider
  tile <- complete |>
    group_by(Lineage, !!sym((i))) |>
    tally() |>
    pivot_wider(names_from = !!((i)), values_from = n)
  # add a column beside previous table showing the total number of isolates in each lineage
  tile$total <- rowSums(tile[, !colnames(tile) %in% "Lineage"], na.rm = TRUE)
  # mutate columns showing the gene and percentage of isolates carrying the gene
  row_for_dataframe <- tile |>
    mutate(
      ALP = i,
      proportion = pos/total
    ) |>
    select(Lineage, ALP, proportion)
  # set NAs to 0
  row_for_dataframe[is.na(row_for_dataframe)] <- 0
  # rbind outputs for individual genes into a complete table for tile plot
  alp_forloop_output <- rbind(row_for_dataframe, alp_forloop_output)
}

# relevelling ALPs
alp_forloop_output$ALP <- factor(alp_forloop_output$ALP, levels = c("rib", "alpha", "alp2/3", "alp1"))

# tile plot
alp_forloop_output |>
  ggplot(aes(x = as.factor(Lineage), y = ALP, fill = proportion)) +
  geom_tile(colour = "black") +
  scale_fill_gradient(low = "white", high = "red") +
  theme(
    axis.ticks.x = element_blank()
  ) +
  labs(
    title = "Alpha-like protein distribution in minor lineages",
    x = "Lineage",
    y = "alp-like proteins",
    fill = "Proportion"
  )

##########
# other surface markers
##########

# marker vector
markers <- c("PI2B", "PI2A2", "PI2A1", "PI1", "srr2", "srr1", "hvgA")

# empty dataframe to store results
markers_forloop_output <- data.frame()

for (i in markers) {
  # table grouped by lineage and surface marker presence/absence, then pivot wider
  tile <- complete |>
    group_by(Lineage, !!sym((i))) |>
    tally() |>
    pivot_wider(names_from = !!((i)), values_from = n)
  # add a column beside previous table showing the total number of isolates in each lineage
  tile$total <- rowSums(tile[, !colnames(tile) %in% "Lineage"], na.rm = TRUE)
  # mutate columns showing the gene and percentage of isolates carrying the gene
  row_for_dataframe <- tile |>
    mutate(
      markers = i,
      proportion = pos/total
    ) |>
    select(Lineage, markers, proportion)
  # set NAs to 0
  row_for_dataframe[is.na(row_for_dataframe)] <- 0
  # rbind outputs for individual genes into a complete table for tile plot
  markers_forloop_output <- rbind(row_for_dataframe, markers_forloop_output)
}

# relevelling surface markers
markers_forloop_output$markers <- factor(markers_forloop_output$markers, levels = c("PI2B", "PI2A2", "PI2A1", "PI1", "srr2", "srr1", "hvgA"))

# tile plot
markers_forloop_output |>
  ggplot(aes(x = as.factor(Lineage), y = markers, fill = proportion)) +
  geom_tile(colour = "black") +
  scale_fill_gradient(low = "white", high = "red") +
  labs(
    title = "Distribution of hvgA surface adhesin, serine rich repeats and pilus islands in minor lineages",
    x = "Lineage",
    y = "Markers",
    fill = "Proportion"
  ) +
  theme(
    axis.ticks.x = element_blank()
  )

##########
# alluvial plot
##########

# data tidying and color assignment
grouped_tally <- complete |>
  group_by(Lineage, cps.final, Host_species, Age_group, Host_status) |>
  tally()

plot_tally <- grouped_tally |>
  mutate(
    host = ifelse(is.na(Age_group), Host_species, Age_group)
  ) |>
  mutate(
    manifest = ifelse(
      host %in% c("fish", "cattle", "camel"), "animal", Host_status)
  )

cps.final <- c("Ia","Ib","II","III","IV","V","VI","VII","VIII","IX","NT")
cps.colours <- c("#E69F00","#56B4E9","#0072B2","#F0E442","#009E73","#D55E00","#BCBDDC","#756BB1","#FA9FB5","#BFFFEB","#808080")
cps_color <- data.frame(cps.final, cps.colours)
plot_tally <- left_join(plot_tally, cps_color, by = "cps.final")

plot_tally$host[plot_tally$host == "human"] <- "human unknown"
plot_tally[is.na(plot_tally)] = "unknown"

# ordering factors
plot_tally$cps.final <- factor(plot_tally$cps.final, levels = c("Ia","Ib","II","III","IV","V","VI","VII","VIII","IX","NT"))
plot_tally$host <- factor(plot_tally$host, levels = c("unknown", "infant", "child", "adult", "human unknown", "camel", "cattle", "fish"))
plot_tally$manifest <- factor(plot_tally$manifest, levels = c("invasive disease", "non-invasive disease", "carriage", "animal", "unknown"))

# plot
plot_tally |>
  ggplot(aes(axis1 = cps.final, axis2 = as.factor(Lineage), axis3 = host, axis4 = manifest, y = n)) +
  geom_alluvium(aes(fill = cps.final), alpha = 0.5) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 4.5) +
  scale_x_discrete(limits = c("Serotype", "Lineage", "Host age/species", "Clinical manifestation")) +
  scale_fill_manual(values = c("#E69F00","#56B4E9","#0072B2","#F0E442","#009E73","#D55E00","#BCBDDC", "#FA9FB5","#BFFFEB","#808080")) +
  labs(
    y = "Number of isolates",
    fill = "Serotype",
    title = "Host and clinical manifestation of serotypes in minor lineages"
  ) + 
  theme(
    legend.text = element_text(size = 16),
    legend.key.size = unit(0.8, "cm"),
    title = element_text(size = 18),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 16)
  )


###############################
# Supplementary figures
###############################


##########
# cps table
##########

# table grouped by lineage and cps, pivot wider
cps_table <- complete |>
  group_by(Lineage, cps.final) |>
  tally() |>
  pivot_wider(names_from = cps.final, values_from = n)

# reorder columns, set NAs to 0, then write CSV
cps_table <- cps_table[, c("Lineage", "Ia","Ib","II","III","IV","V","VI","VIII","IX","NT")]
cps_table[is.na(cps_table)] <- 0
write_csv(cps_table, file = "cps table.csv")


##########
# host table
##########

# table grouped by lineage and host species, pivot wider
host <- complete |>
  group_by(Lineage, Host_species) |>
  tally() |>
  pivot_wider(names_from = Host_species, values_from = n)

# reorder columns, set NAs to 0s, then write CSV
host <- host[, c("Lineage", "camel", "cattle", "fish", "human", "NA")]
host[is.na(host)] <- 0
write_csv(host, file = "host.csv")


##########
# marker table
##########

# mutate columns showing isolates negative for any ALP/SRR/PI
complete_surfaceneg <- complete |>
  mutate(
    srr_neg = ifelse(srr1 %in% "neg" & srr2 %in% "neg", "pos", "neg"),
    pi_neg = ifelse(PI1 %in% "neg" & PI2A1 %in% "neg" & PI2A2 %in% "neg" & PI2B %in% "neg", "pos", "neg"),
    alp_neg = ifelse(alp1 %in% "neg" & `alp2/3` %in% "neg" & alpha %in% "neg" & rib %in% "neg", "pos", "neg")
  )

# number of isolates in each lineage that carry each marker, with percentage in brackets
marker_table <- complete_surfaceneg |>
  group_by(Lineage) |>
  summarise(
    across(!c(Sample_ID:tkt),
           ~if (sum(. %in% "pos") > 0) {
             paste0(sum(. %in% "pos"), " (", round(sum(. %in% "pos")/n() * 100, 1), "%)")
           } else{
             "0"
           }
      )
    )


##########
# Fisher's test
##########

# gene vector
gene <- c("tetS", "tetO", "tetM", "tetL", "msrD", "mefA", "lsaC", "lnuC", 
          "ermT", "ermB", "ermA", "catQ", "cat(pc194)", "aph(3'-III)", "ant(6-Ia)", 
          "aadE", "aac(6')-aph(2'')", "alp1", "alp2/3", "alpha", 
          "rib", "hvgA", "srr1", "srr2", "PI1", "PI2A1", "PI2A2", "PI2B")

# define lineages
unique_lineage <- unique(complete$Lineage)

# create dataframe with specific columns to store results
results <- data.frame(Lineage = character(), Gene = character(), PValue = numeric(), stringsAsFactors = FALSE)

for (Lineage in unique_lineage) {
  # for each loop, set the looped lineage as "lineage" and others as "others"
  complete$Lineage_binary <- ifelse(complete$Lineage == Lineage, "target", "other")
  
  for (i in gene) {
    # for each looped gene, create contingency table showing gene presence/absence in target/other lineage
    contingency_table <- table(complete[[i]], complete$Lineage_binary)
    # perform two-sided fishers test on the contingency table
    fisher_test <- fisher.test(contingency_table, alternative = "t")
    # rbind results into pre-defined dataframe
    results <- rbind(
      results, data_frame(
        Lineage = Lineage,
        Gene = i,
        PValue = fisher_test$p.value
      )
    )
  }
}

# adjust p-value using bonferroni method
results$PAdjusted <- p.adjust(results$PValue, method = "bonferroni")

# tidy up results into table format suitable for report
tidy_results <- results |>
  select(Lineage, Gene, PAdjusted) |>
  pivot_wider(names_from = Gene, values_from = PAdjusted)

write_csv(tidy_results, file = "fisher ref.csv")


##########
# Host distribution of ALP negative isolates
##########

# mutate alp negative column
complete_alpneg <- complete |>
  mutate(
    alp_neg = ifelse(alp1 %in% "neg" & `alp2/3` %in% "neg" & alpha %in% "neg" & rib %in% "neg", "pos", "neg")
  )

# filter alp negative isolates and plot bar chart
complete_alpneg |>
  filter(alp_neg %in% "pos") |>
  group_by(Lineage) |>
  ggplot(aes(x = as.factor(Lineage), fill = Host_species)) +
  geom_bar()+
  scale_fill_manual(values = c("#8975ca", "#cb5582", "#71a659", "#ff7837", "#808080")) +
  labs(
    title = "Host distribution of alpha-like protein negative isolates",
    x = "Lineage", 
    y = "Number of isolates",
    fill = "Host species"
  )
