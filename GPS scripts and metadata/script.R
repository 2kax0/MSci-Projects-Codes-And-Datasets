###############################
# Global Whole-genome sequencing analysis of antimicrobial resistance in Streptococcus pneumoniae 
# Author: Max Choong
###############################


# load packages and read file
library(tidyverse)
library(scatterpie)
library(ggpubr)
library(meta)
library(maps)
library(ggpattern)


complete <- read_csv(file = "metadata.csv")


# defining vectors for nvt and vt list for geom pattern, antibiotics examined
PCV13 <- c("1", "3", "4", "5", "6A", "6B", "7F", "9V", "14", "18C", "19A", "19F", "23F")
list <- c("NVT" = "stripe", "VT" = "none")
antibiotic <- c("TAX (meningital)", "CFT (meningital)", "MER", "ERY")



##########
# heat map
##########

# initialise map and removing antarctica
world_map <- map_data("world")
world_map <- subset(world_map, region != "Antarctica")

# country table
tally <- complete |>
  group_by(Country) |>
  tally()

# need to manually change country to fit format (i.e. AUSTRALIA to Australia)
write_csv(tally, file = "country.csv")
tally <- read_csv("country edited.csv")

# filter out unknown
tally <- tally |>
  filter(!Country %in% "Unknown")

# heatmap
ggplot(tally) +
  geom_map(dat = world_map, map = world_map, aes(map_id = region),
           fill = "white", colour = "black"
  ) +
  geom_map(map = world_map, aes(map_id = Country, fill = n)) +
  scale_fill_gradient(low = "#D4F1F4", high = "#517891", name = "Number of isolates") +
  expand_limits(x = world_map$long, y = world_map$lat) +
  labs(
    title = "Geographical distribution of isolates"
  ) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 12),
    title = element_text(size = 18),
    legend.text = element_text(size = 10)
  )

##########
# serotype, gpsc, AMR figure
##########

# create a table of number of GPSCs, then filtering out only the top 40
complete$GPSC <- sub(".*?;", "", complete$GPSC)

GPSC_tally <- complete |>
  group_by(GPSC) |>
  tally()

GPSC_40 <- GPSC_tally[order(GPSC_tally$n, decreasing = TRUE),] |>
  filter(!is.na(GPSC)) |>
  slice_head(n = 40) 

# extract top 40 most abundant GPSCs
n_40 <- unique(GPSC_40$GPSC)

# classify not top 40 as others
tidy <- complete |>
  mutate(
    GPSC_tidy = ifelse(GPSC %in% n_40, GPSC, "NA")
  ) |>
  filter(!Serotype %in% "_" & !is.na(GPSC_tidy)) |>
  relocate(GPSC_tidy, .after = GPSC)

# serotype tidying
tidy <- tidy |>
  mutate(
    serotype_modified = ifelse(Serotype %in% c("6A(6A-I)", "6A(6A-II)", "6A(6A-III)", "6A(6A-IV)", "6A(6A-V)", "6E(6A)", "6A(6A-VI)"), "6A",
                               ifelse(Serotype %in% c("6B(6B-I)", "6B(6B-II)", "6E(6B)"), "6B",
                                      ifelse(Serotype %in% "11A(11F_like)", "11A",
                                             ifelse(Serotype %in% c("15B", "15C"), "15B/15C",
                                                    ifelse(Serotype %in% c("19A(19A-I)", "19A(19A-I/19A-II)", "19A(19A-II)", "19A(19A-III)"), "19A",
                                                           ifelse(Serotype %in% c("19F(19AF)", "19F(19F-I)", "19F(19F-II)", "19F(19F-III)", "19F(19F-IV)"), "19F",
                                                                  ifelse(Serotype %in% c("23B1", "23B(23B1)", "23B(23B1)"), "23B",
                                                                         ifelse(Serotype %in% c("24B/24C/24F", "24B/24F"), "24F",
                                                                                ifelse(Serotype %in% c("33F(33F-1a)", "33F(33F-1b)"), "33F",
                                                                                       ifelse(Serotype %in% c("NCC2_aliC_aliD_non_encapsulated", "NCC1_pspK_non_encapsulated", "NCC2_S_mitis_aliC_aliD_non_encapsulated", "NCC3_aliD_non_encapsulated"), "non-encapsulated", 
                                                                                              ifelse(Serotype %in% c("possible 06A	 but wciP gene might not be complete", "possible 06E	 but wciP gene might not be complete", "possible 06D	 but wciP gene might not be complete"), "untypable", Serotype))
                                                                                ))))))))))



# join tidied dataset with colour for annotation
serotype_colours <- read_csv(file = "serotype_colours.csv")

coloured <- left_join(tidy, serotype_colours, by = "serotype_modified") |>
  relocate(In_silico_serotype__colour, .after = )

# vaccine and nvt annotation
coloured$Vaccine <- ifelse(
  coloured$serotype_modified %in% PCV13, "VT", "NVT"
)

# joining tidied dataset with annotated serotype dataframes, then set names for colour and return serotype
graph_dataset <- coloured |>
  group_by(GPSC_tidy, serotype_modified) |>
  tally()

graph_dataset <- left_join(graph_dataset, serotype_colours, by = "serotype_modified")

serotype_colors_vector <- setNames(graph_dataset$In_silico_serotype__colour, graph_dataset$serotype_modified)

# vaccine and nvt annotation
graph_dataset$Vaccine <- ifelse(
  graph_dataset$serotype_modified %in% PCV13, "VT", "NVT"
)

# barchart
bar <- graph_dataset |>
  filter(!GPSC_tidy %in% "NA") |>
  ggplot() +
  (aes(x = as.factor(GPSC_tidy), y = n, fill = serotype_modified, pattern = Vaccine)) +
  geom_col_pattern(
    pattern_spacing = 0.01, 
    pattern_fill = "black", 
    pattern_size = 0.5, 
    pattern_linetype = 0.5,
    pattern_key_scale_factor = 0.5,
    pattern_orientation = "vertical", show.legend = TRUE) +
  scale_fill_manual(name = "Serotype",
                    values = serotype_colors_vector,
                    guide = guide_legend(override.aes = list(pattern = "none"), ncol = 3)) +
  scale_pattern_manual(
    name = "Vaccine",
    values = c("NVT" = "stripe", "VT" = "none")
  )  +
  theme(
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    axis.text.x = element_text(vjust = 0.9),
    plot.margin = unit(c(2,0,-6,0), 'lines')
  ) +
  scale_pattern_fill_manual(values = list) +
  labs(
    title = "Serotype and antibiotic resistance within the 40 most abundant GPSCs",
    x = "GPSC",
    y = "Number of isolates",
    fill = "Serotype"
  )

# geom tile
# empty dataframe 
forloop_output <- data.frame()
csv_output = list()

# for loop to create dataframe for geom tile
for (i in antibiotic) {
  # group by GPSC and looped antibiotic, then pivot wider
  tile <- tidy |>
    group_by(GPSC_tidy, !!sym((i))) |>
    tally() |>
    pivot_wider(names_from = !!sym((i)), values_from = n)
  # sum up rows for each GPSC
  tile$total <- rowSums(test_tile[, !colnames(tile) %in% "GPSC_tidy"], na.rm = TRUE)
  # use previous sum to calculate percentage, then select relevant columns
  row_for_dataframe <- tile |>
    mutate(
      antibiotic = i,
      proportion = R/total
    ) |>
    select(GPSC_tidy, proportion, antibiotic)
  # set NAs as 0
  row_for_dataframe[is.na(row_for_dataframe)] <- 0
  # pivot wider for final output csv
  dataframe_row <- row_for_dataframe |>
    pivot_wider(names_from = GPSC_tidy, values_from = proportion)
  # add antibiotic label to row
  csv_output[[i]] <- dataframe_row
  # rbind
  forloop_output <- rbind(row_for_dataframe, forloop_output)
}
# pivot wider for easier reading
resistance_reference <- forloop_output |>
  pivot_wider(names_from = antibiotic, values_from = proportion)

# write reference file as csv
write_csv(resistance_reference, file = "tile reference.csv")

# relevelling factors
forloop_output$antibiotic <- factor(forloop_output$antibiotic, levels = c("ERY", "MER", "CFT (meningital)", "TAX (meningital)"))

# plot geom tile using for loop output
tile <- forloop_output |>
  filter(!GPSC_tidy %in% "NA") |>
  ggplot(aes(x = GPSC_tidy, y = antibiotic, fill = proportion)) +
  geom_tile(colour = "black") +
  scale_fill_gradientn(colours = c("white", "pink", "red"), values = c("0", "0.25", "0.5", "0.75", "1.00")) +
  coord_fixed() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(-4,0,0,0), 'lines')
  ) +
  labs(
    y = "Antibiotic",
    fill = "Proportion"
  )

plot(tile)

# complete figure generation
figure <- ggarrange(bar, tile + rremove("xlab"),
                    labels = NULL,
                    ncol = 1, align = "v") 
plot(figure)

##########
# year, manifest, region plot
##########

# filter isolates with no region information or no manifest information, group by continent, year, manifest
# tally then mutate another column to join region and manifest information
continent_year_tally <- complete |>
  filter(!is.na(Continent) & !is.na(Manifest_category) & !Manifest_category %in% "") |>
  group_by(Continent, Year, Manifest_category) |>
  tally() |>
  mutate(
    join = paste0(Continent, " ", Manifest_category)
  )

# manually relevel manifest
continent_year_tally$Manifest_category <- factor(continent_year_tally$Manifest_category, levels = c("Carriage", "IPD", "non IPD", "NA"))

# year, manifest, region plot
continent_year_tally |>
  ggplot(aes(x = as.factor(Year), y = n, fill = Manifest_category)) +
  geom_bar(stat = "identity") +
  facet_grid(Continent ~ .) +
  scale_fill_manual(values = c("#BA55D3", "#00BFFF", "#FF8C00")) +
  labs(
    title = "Clinical manifestation of isolates collected in each region",
    x = "Year",
    y = "Number of isolates",
    fill = "Clinical \nManifestation"
  ) +
  geom_vline(xintercept = "2000", linetype = "dashed") +
  geom_vline(xintercept = "2010", linetype = "dashed") +
  geom_vline(xintercept = "2020", linetype = "dashed") +
  theme(
    legend.title = element_text(size = 12),
    title = element_text(size = 14),
    legend.text = element_text(size = 10)
  )


##########
# panels
##########

### TAX (cefotaxime), CFT (ceftriaxone), MER (meropenem), ERY (erythromycin)
antibiotic <- c("TAX (meningital)", "CFT (meningital)", "MER", "ERY")
Continent <- c("AFRICA", "NORTH AMERICA", "CENTRAL AMERICA", "LATIN AMERICA", 
               "EUROPE", "ASIA", "OCEANIA")

world <- map_data("world")
world <- subset(world, region != "Antarctica")

### panel and forest plot for loops
for(i in antibiotic){
  
  #filter out isolates that failed resistance prediction
  complete <- complete |> filter(!!((sym(i))) %in% c("S", "I", "R"))
  
  ### create new directory for each looped variable
  newdir <- paste0(i, " plots")
  
  dir.create(newdir)
  
  ### data tidying, relevel factors
  complete$Age_category <- factor(complete$Age_category, levels = c("< 5", ">=5 - <=65", "> 65"))
  complete$Manifest_category <- factor(complete$Manifest_category, levels = c("Carriage", "IPD", "non IPD", "NA"))
  complete[[i]] <- factor(complete[[i]], levels = c("S", "I", "R"))
  
  ### filter out isolates with no age, manifest information or failed SIR prediction
  # then group by age, SIR and manifest
  age_tally <- complete |>
    filter(!is.na(Age_category)) |>
    filter(!is.na(Manifest_category)) |>
    filter(!!((sym(i))) %in% c("S", "I", "R")) |>
    group_by(Age_category, !!((sym(i))), Manifest_category) |>
    tally()
  
  ### filter out isolates with no age, manifest information or failed SIR prediction
  # calculate total, then group by age
  total_tally <- complete |>
    filter(!is.na(Age_category)) |>
    filter(!is.na(Manifest_category)) |>
    filter(!!((sym(i))) %in% c("S", "I", "R")) |>
    group_by(Age_category) |>
    tally()
  
  # join SIR and total
  join <- left_join(age_tally, total_tally, by = "Age_category")

  ### define path to save csv
  age_manifest_path <- file.path("/path", "to", "directory", paste(newdir, "/age manifest table.csv", sep = ""), fsep = "/")
  
  ### csv file for number of resistant isolates in each age group for each disease manifestation
  write_csv(join, file = age_manifest_path)
  
  ### calculate percentage
  age_tally_manifest <- join |>
    summarise(
      res_manifest = ifelse(!.data[[i]] == "R", "S", as.character(Manifest_category)),
      percentage = n.x/n.y,
      .groups = "keep"
    ) 
  
  ### manual relevelling of factors
  age_tally_manifest$res_manifest <- factor(age_tally_manifest$res_manifest, levels = c("S", "Carriage", "IPD", "non IPD"))
  age_tally_manifest$Age_category <- factor(age_tally_manifest$Age_category, levels = c("< 5", ">=5 - <=65", "> 65"))
  
  ### manifestation in resistant isolates within each age group
  age_manifest <- age_tally_manifest |>
    filter(!res_manifest %in% "S") |>
    ggplot(aes(x = Age_category, y = percentage, fill = res_manifest)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("#BA55D3", "#00BFFF", "#FF8C00")) +
    scale_y_continuous(labels = scales::percent) +
    theme(
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 6),
      legend.key.size = unit(0.2, "cm"),
      plot.margin = margin(0, 0, 0, 5),
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 7),
      axis.text.x = element_text(size = 5)
    ) +
    labs(
      fill = "Clinical \nmanifestation",
      y = "Percentage",
      x = "Age category"
    ) 
  

  # tally looped antibiotic to check if there is "I" 
  grouping <- complete |>
    filter(!!((sym(i))) %in% c("S", "I", "R")) |>
    group_by(!!((sym(i)))) |>
    tally() |>
    pivot_wider(names_from = !!((sym(i))), values_from = n)
  
  ### if looped antibiotic has "I", else
  if("I" %in% names(grouping)){
    # read csv again to prevent exclusion of isolates that were filtered out from previous iteration
    complete <- read_csv(file = "metadata.csv")
    # SIR tally
    pivot_SIR <- complete |>
    filter(!is.na(Continent)) |>
    filter(!!((sym(i))) %in% c("S", "I", "R")) |>
    group_by(Continent, !!((sym(i)))) |>
    tally() |>
    pivot_wider(names_from = !!((sym(i))), values_from = n)
  # set NAs as 0
  pivot_SIR[is.na(pivot_SIR)] = 0
  # add longitude or latitude
  pivot_plot <- pivot_SIR |>
    mutate(
      Longitude = as.numeric(case_when(
        Continent == "AFRICA" ~ 21.890909,
        Continent == "NORTH AMERICA" ~ -100.538118,
        Continent == "CENTRAL AMERICA" ~ -91.573274,
        Continent == "LATIN AMERICA" ~ -53.604523,
        Continent == "EUROPE" ~ 15.738565,
        Continent == "ASIA" ~ 110.308879,
        Continent == "OCEANIA" ~ 150.562787
      )),
      Latitude = as.numeric(case_when(
        Continent == "AFRICA" ~ 6.208920,
        Continent == "NORTH AMERICA" ~ 49.249542,
        Continent == "CENTRAL AMERICA" ~ 18.453068,
        Continent == "LATIN AMERICA" ~ -12.581874,
        Continent == "EUROPE" ~ 50.890969,
        Continent == "ASIA" ~ 41.691447,
        Continent == "OCEANIA" ~ -28.553076
      )))
  # calculate row total
  pivot_plot$Total <- rowSums(pivot_plot[, !colnames(pivot_plot) %in% c("Continent", "Longitude", "Latitude")], na.rm = TRUE)
  } else
    {
    
    complete <- read_csv(file = "metadata.csv")    
      
    pivot_SIR <- complete |>
    filter(!is.na(Continent)) |>
    filter(!!((sym(i))) %in% c("S", "R")) |>
    group_by(Continent, !!((sym(i)))) |>
    tally() |>
    pivot_wider(names_from = !!((sym(i))), values_from = n)
  
  pivot_SIR[is.na(pivot_SIR)] = 0
  
  pivot_plot <- pivot_SIR |>
    mutate(
      Longitude = as.numeric(case_when(
        Continent == "AFRICA" ~ 21.890909,
        Continent == "NORTH AMERICA" ~ -100.538118,
        Continent == "CENTRAL AMERICA" ~ -91.573274,
        Continent == "LATIN AMERICA" ~ -53.604523,
        Continent == "EUROPE" ~ 15.738565,
        Continent == "ASIA" ~ 110.308879,
        Continent == "OCEANIA" ~ 150.562787
      )),
      Latitude = as.numeric(case_when(
        Continent == "AFRICA" ~ 6.208920,
        Continent == "NORTH AMERICA" ~ 49.249542,
        Continent == "CENTRAL AMERICA" ~ 18.453068,
        Continent == "LATIN AMERICA" ~ -12.581874,
        Continent == "EUROPE" ~ 50.890969,
        Continent == "ASIA" ~ 41.691447,
        Continent == "OCEANIA" ~ -28.553076
      )))
  
  pivot_plot$Total <- rowSums(pivot_plot[, !colnames(pivot_plot) %in% c("Continent", "Longitude", "Latitude")], na.rm = TRUE)
    }
  
  # scatterpie with "S" and "R include "I" if looped antibiotic has intermediate
  
  if(!"I" %in% names(pivot_plot)){
    plot <- ggplot() +
      geom_polygon(data = world, aes(x=long, y = lat, group = group)) + 
      coord_fixed(1.3) +
      geom_scatterpie(data = pivot_plot, aes(x = Longitude, y = Latitude, r = 13), cols = c("S", "R"), color = NA) + 
      coord_equal() +
      scale_fill_manual(values = c("#0069EC", "#FF2722")) +
      labs(
        x = "Longitude",
        y = "Latitude",
        fill = "Resistance"
      ) +
      theme(
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.4, "cm"),
        legend.position.inside = c(.1,.35),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        plot.margin = margin(-20, 0, 0, 0)
      ) 
   
  } else
    plot <- ggplot() +
    geom_polygon(data = world, aes(x=long, y = lat, group = group)) + 
    coord_fixed(1.3) +
    geom_scatterpie(data = pivot_plot, aes(x = Longitude, y = Latitude, r = 13), cols = c("S", "I", "R"), color = NA) + 
    coord_equal() +
    scale_fill_manual(values = c("#0069EC", "#F797B1", "#FF2722")) +
    labs(
      x = "Longitude",
      y = "Latitude",
      fill = "Resistance"
    ) +
    theme(
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 6),
      legend.key.size = unit(0.4, "cm"),
      legend.position.inside = c(.1,.35),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid.major = element_blank(),
      plot.margin = margin(-20, 0, 0, 0)
    ) 
  # indicate path to save continent SIR file
  SIR_path <- file.path("/path", "to", "directory", paste(newdir, "/SIR table.csv", sep = ""), fsep = "/")
  write_csv(pivot_plot, file = SIR_path)
  
  # add labels to plots
  map_manifest <- ggarrange(plot, age_manifest, widths = c(1.75, 1.25), heights = c(1, 1), labels = c("(a)", "(b)"), label.y = 1.01, font.label=list(color="black",size=8), hjust = c(0,0,0)) +
    theme(
      plot.margin = margin(10, 0, 0, 0)
    )
  # manual releveling 
  complete[[i]] <- factor(complete[[i]], levels = c("S", "I", "R"))
  
  ### year plot
  # include "I" if looped antibiotic has intermediate, else
  if(!"I" %in% names(pivot_plot)){
  year <- complete |>
    filter(!!((sym(i))) %in% c("S", "I", "R")) |>
    filter(!Year %in% "_" & !is.na(Year)) |>
    ggplot(aes(x = as.factor(Year), fill = !!((sym(i))))) +
    geom_bar() +
    scale_fill_manual(values = c("#0069EC", "#FF2722")) +
    theme(
      axis.text = element_text(size = 7),
      plot.title = element_text(size = 8),
      legend.position = "none",
      axis.title = element_text(size = 7),
      axis.text.x = element_text(size = 5)
    ) +
    labs(
      y = "Number of isolates",
      x = "Year",
      fill = "Resistance"
    ) +
    geom_vline(xintercept = "2000", linetype = "dashed") +
    geom_vline(xintercept = "2010", linetype = "dashed") +
    geom_vline(xintercept = "2020", linetype = "dashed") 
  
  } else
    year <- complete |>
    filter(!!((sym(i))) %in% c("S", "I", "R")) |>
    filter(!Year %in% "_" & !is.na(Year)) |>
    ggplot(aes(x = as.factor(Year), fill = !!((sym(i))))) +
    geom_bar() +
    scale_fill_manual(values = c("#0069EC", "#F797B1", "#FF2722")) +
    theme(
      axis.text = element_text(size = 7),
      legend.position = "none",
      axis.title = element_text(size = 7),
      axis.text.x = element_text(size = 5)
    ) +
    labs(
      y = "Number of isolates",
      x = "Year",
      fill = "Resistance"
    ) +
    geom_vline(xintercept = "2000", linetype = "dashed") +
    geom_vline(xintercept = "2010", linetype = "dashed") +
    geom_vline(xintercept = "2020", linetype = "dashed") 
  
  ### serotype GPSC plot, only show the top 30 most abundant GPSC
  # collapse GPSC
  complete$GPSC <- sub(".*?;", "", complete$GPSC)
  # group resistant isolates by GPSC
  n30 <- complete |>
    filter(!!((sym(i))) %in% "R" & !is.na(GPSC)) |>
    group_by(GPSC) |>
    tally() 
  # extract 30 GPSCs with highest number of resistant isolates
  tally_30 <- n30[order(n30$n, decreasing = TRUE),] |>
    slice_head(n = 30)
  # create vector with 30 GPSCs with highest number of resistant isolates
  GPSC_list <- unique(tally_30$GPSC) 
  # serotype clean up
  complete_modified <- complete |>
    mutate(
      serotype_modified = ifelse(Serotype %in% c("6A(6A-I)", "6A(6A-II)", "6A(6A-III)", "6A(6A-IV)", "6A(6A-V)", "6E(6A)"), "6A",
                                 ifelse(Serotype %in% c("6B(6B-I)", "6B(6B-II)", "6E(6B)"), "6B",
                                        ifelse(Serotype %in% "11A(11F_like)", "11A",
                                               ifelse(Serotype %in% c("15B", "15C"), "15B/15C",
                                                      ifelse(Serotype %in% c("19A(19A-I)", "19A(19A-I/19A-II)", "19A(19A-II)", "19A(19A-III)"), "19A",
                                                             ifelse(Serotype %in% c("19F(19AF)", "19F(19F-I)", "19F(19F-II)", "19F(19F-III)", "19F(19F-IV)"), "19F",
                                                                    ifelse(Serotype %in% c("23B1", "23B(23B1)", "23B(23B1)"), "23B",
                                                                           ifelse(Serotype %in% c("24B/24C/24F", "24B/24F"), "24F",
                                                                                  ifelse(Serotype %in% c("33F(33F-1a)", "33F(33F-1b)"), "33F",
                                                                                         ifelse(Serotype %in% c("NCC2_aliC_aliD_non_encapsulated", "NCC1_pspK_non_encapsulated", "NCC2_S_mitis_aliC_aliD_non_encapsulated", "NCC3_aliD_non_encapsulated"), "non-encapsulated", 
                                                                                                ifelse(Serotype %in% c("possible 06A	 but wciP gene might not be complete", "possible 06E	 but wciP gene might not be complete", "possible 06D	 but wciP gene might not be complete"), "untypable", Serotype))
                                                                                  ))))))))),
      gpsc_tidy = ifelse(GPSC %in% GPSC_list, complete$GPSC, "others")
    )
  # join serotype with colour
  complete_coloured <- left_join(complete_modified, colour, by = "serotype_modified")
  # set names for colour and return serotype
  serotype_colours_vector <- setNames(complete_coloured$In_silico_serotype__colour, complete_coloured$serotype_modified)
  # create separate bar chart for GPSC1
  GPSC_serotype_1 <- complete_coloured |>
    filter(!!((sym(i))) %in% "R" & !gpsc_tidy %in% "others" & gpsc_tidy %in% "1") |>
    ggplot(aes(x = as.factor(gpsc_tidy), fill = serotype_modified)) +
    geom_bar() +
    scale_fill_manual(name = "Serotype",
                      values = serotype_colours_vector) +
    labs(
      x = "GPSC",
      y = "Number of isolates"
    ) +
    theme(
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 6),
      legend.key.size = unit(0.2, "cm"),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 7),
      axis.text.x = element_text(size = 5),
      plot.margin = margin(15, 0, 0, 0)
    )
  # create bar chart for others
  GPSC_serotype_others <- complete_coloured |>
    filter(!!((sym(i))) %in% "R" & !gpsc_tidy %in% "others" & !gpsc_tidy %in% "1") |>
    ggplot(aes(x = as.factor(gpsc_tidy), fill = serotype_modified)) +
    geom_bar() +
    scale_fill_manual(name = "Serotype",
                      values = serotype_colours_vector) +
    labs(
      x = "GPSC",
      y = "Number of isolates"
    ) +
    theme(
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 6),
      legend.key.size = unit(0.2, "cm"),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 7),
      axis.text.x = element_text(size = 5),
      plot.margin = margin(15, 0, 0, 0)
    )
  
  # arrange bar charts
  GPSC_serotype_plot <- ggarrange(GPSC_serotype_1, GPSC_serotype_others, widths = c(1, 5)) 
  # combine bar charts with previous plots
  GPSC_serotype_year <- ggarrange(year, GPSC_serotype_plot, ncol = 1, labels = c("(c)", "(d)"), font.label=list(color="black",size=8), hjust = c(0,0,0))
  # combining to create top half of panel
  full_plot <- ggarrange(map_manifest, GPSC_serotype_year, ncol = 1, heights = c(1.5, 2)) +
    theme(
      plot.margin = margin(0.5,0.1,0.5,0.1, "cm")
    )

  ### forest plots
  # set path to save full plot
  full_plot_path <- file.path("/path", "to", "directory", paste(newdir, "/full plot.png", sep = ""), fsep = "/")
  
  ggsave(full_plot, filename = full_plot_path, dpi = 320, width = 26.58, height = 16.91, units = "cm")
  
  # define and create new directory for each looped antibiotic
  newdir <- paste0(i, " forest plot")
  
  dir.create(newdir)
  
  # define path to save created png
  napath <- file.path("/path", "to", "directory", paste(newdir, "/na forest.png", sep = ""), fsep = "/")
  # filter out failed predictions
  complete <- read_csv(file = "metadata") |>
    filter(!!((sym(i))) %in% c("S", "I", "R"))
  
  #forest plot for resistance in north american countries
  NorthAmerica <- complete |>
    filter(Continent %in% "NORTH AMERICA") |>
    group_by(Country, !!sym((i))) |>
    tally() |>
    pivot_wider(names_from = !!sym((i)), values_from = n)
  # filter out R and perform metaprop analysis, if there are no R isolates set to 0 then run metaprop
  if(!"R" %in% names(NorthAmerica)) {NorthAmerica$R <- 0
  NorthAmerica$total <- rowSums(NorthAmerica[, !colnames(NorthAmerica) %in% "Country"], na.rm = TRUE)
  NorthAmerica[is.na(NorthAmerica)] <- 0
  nameta <- metaprop(R, total, data = NorthAmerica)
  }
  
  else NorthAmerica$total <- rowSums(NorthAmerica[, !colnames(NorthAmerica) %in% "Country"], na.rm = TRUE)
  NorthAmerica[is.na(NorthAmerica)] <- 0
  nameta <- metaprop(R, total, data = NorthAmerica)
  # initiate graphic device, create forest plot then save
  png(file = napath, width = 765, height = 169)
  forest(nameta, comb.fixed = nameta$comb.fixed, comb.random = nameta$comb.random, 
         overall = TRUE, studlab = NorthAmerica$Country, plotwidth=unit(10, "cm"), xlim = c(0, 1), colgap.forest = "2cm")
  dev.off()
  
  #central america
  CentralAmerica <- complete |>
    filter(Continent %in% "CENTRAL AMERICA") |>
    group_by(Country, !!sym((i))) |>
    tally() |>
    pivot_wider(names_from = !!sym((i)), values_from = n)
  # define path to save created png
  capath <- file.path("/path", "to", "directory", paste(newdir, "/ca forest.png", sep = ""), fsep = "/")
  # filter out R and perform metaprop analysis, if there are no R isolates set to 0 then run metaprop
  if(!"R" %in% names(CentralAmerica)) {CentralAmerica$R <- 0
  CentralAmerica$total <- rowSums(CentralAmerica[, !colnames(CentralAmerica) %in% "Country"], na.rm = TRUE)
  CentralAmerica[is.na(CentralAmerica)] <- 0
  cameta <- metaprop(R, total, data = CentralAmerica) } 
  
  else CentralAmerica$total <- rowSums(CentralAmerica[, !colnames(CentralAmerica) %in% "Country"], na.rm = TRUE)
  CentralAmerica[is.na(CentralAmerica)] <- 0
  cameta <- metaprop(R, total, data = CentralAmerica)
  # initiate graphic device, create forest plot then save
  png(file = capath, width = 765, height = 219)
  forest(cameta, comb.fixed = cameta$comb.fixed, comb.random = cameta$comb.random, 
         overall = TRUE, studlab = CentralAmerica$Country, plotwidth=unit(10, "cm"), xlim = c(0, 1), colgap.forest = "2cm")
  dev.off()
  
  #latin america
  LatinAmerica <- complete |>
    filter(Continent %in% "LATIN AMERICA") |>
    group_by(Country, !!sym((i))) |>
    tally() |>
    pivot_wider(names_from = !!sym((i)), values_from = n)
  # define path to save created png
  lapath <- file.path("/path", "to", "directory", paste(newdir, "/la forest.png", sep = ""), fsep = "/")
  # filter out R and perform metaprop analysis, if there are no R isolates set to 0 then run metaprop
  if(!"R" %in% names(LatinAmerica)) {LatinAmerica$R <- 0
  LatinAmerica$total <- rowSums(LatinAmerica[, !colnames(LatinAmerica) %in% "Country"], na.rm = TRUE)
  LatinAmerica[is.na(LatinAmerica)] <- 0
  lameta <- metaprop(R, total, data = LatinAmerica)
  }
  
  else LatinAmerica$total <- rowSums(LatinAmerica[, !colnames(LatinAmerica) %in% "Country"], na.rm = TRUE)
  LatinAmerica[is.na(LatinAmerica)] <- 0
  lameta <- metaprop(R, total, data = LatinAmerica)
  # initiate graphic device, create forest plot then save
  png(file = lapath, width = 765, height = 249)
  forest(lameta, comb.fixed = lameta$comb.fixed, comb.random = lameta$comb.random, 
         overall = TRUE, studlab = LatinAmerica$Country, plotwidth=unit(10, "cm"), xlim = c(0, 1), colgap.forest = "2cm")
  dev.off()
  
  #europe
  Europe <- complete |>
    filter(Continent %in% "EUROPE") |>
    group_by(Country, !!sym((i))) |>
    tally() |>
    pivot_wider(names_from = !!sym((i)), values_from = n)
  # define path to save created png
  eupath <- file.path("/path", "to", "directory", paste(newdir, "/eu forest.png", sep = ""), fsep = "/")
  # filter out R and perform metaprop analysis, if there are no R isolates set to 0 then run metaprop
  if(!"R" %in% names(Europe)){Europe$R <- 0
  Europe$total <- rowSums(Europe[, !colnames(Europe) %in% "Country"], na.rm = TRUE)
  Europe[is.na(Europe)] <- 0
  eumeta <- metaprop(R, total, data = Europe)
  }
  
  else Europe$total <- rowSums(Europe[, !colnames(Europe) %in% "Country"], na.rm = TRUE)
  Europe[is.na(Europe)] <- 0
  eumeta <- metaprop(R, total, data = Europe)
  # initiate graphic device, create forest plot then save
  png(file = eupath, width = 765, height = 409)
  forest(eumeta, comb.fixed = eumeta$comb.fixed, comb.random = eumeta$comb.random, 
         overall = TRUE, studlab = Europe$Country, plotwidth=unit(10, "cm"), xlim = c(0, 1), colgap.forest = "2cm")
  dev.off()
  
  #asia
  Asia <- complete |>
    filter(Continent %in% "ASIA") |>
    group_by(Country, !!sym((i))) |>
    tally() |>
    pivot_wider(names_from = !!sym((i)), values_from = n)
  # define path to save created png
  aspath <- file.path("/path", "to", "directory", paste(newdir, "/as forest.png", sep = ""), fsep = "/")
  # filter out R and perform metaprop analysis, if there are no R isolates set to 0 then run metaprop
  if(!"R" %in% names(Asia)){Asia$R <- 0
  Asia$total <- rowSums(Asia[, !colnames(Asia) %in% "Country"], na.rm = TRUE)
  Asia[is.na(Asia)] <- 0
  asmeta <- metaprop(R, total, data = Asia)
  }
  
  else Asia$total <- rowSums(Asia[, !colnames(Asia) %in% "Country"], na.rm = TRUE)
  Asia[is.na(Asia)] <- 0
  asmeta <- metaprop(R, total, data = Asia)
  # initiate graphic device, create forest plot then save
  png(file = aspath, width = 765, height = 409)
  forest(asmeta, comb.fixed = asmeta$comb.fixed, comb.random = asmeta$comb.random, 
         overall = TRUE, studlab = Asia$Country, plotwidth=unit(10, "cm"), xlim = c(0, 1), colgap.forest = "2cm")
  dev.off()
  
  #africa
  Africa <- complete |>
    filter(Continent %in% "AFRICA") |>
    group_by(Country, !!sym((i))) |>
    tally() |>
    pivot_wider(names_from = !!sym((i)), values_from = n)
  # define path to save created png
  afpath <- file.path("/path", "to", "directory", paste(newdir, "/af forest.png", sep = ""), fsep = "/")
  # filter out R and perform metaprop analysis, if there are no R isolates set to 0 then run metaprop
  if(!"R" %in% names(Africa)){Africa$R <- 0
  Africa$total <- rowSums(Africa[, !colnames(Africa) %in% "Country"], na.rm = TRUE)
  Africa[is.na(Africa)] <- 0
  afmeta <- metaprop(R, total, data = Africa)
  }
  
  else Africa$total <- rowSums(Africa[, !colnames(Africa) %in% "Country"], na.rm = TRUE)
  Africa[is.na(Africa)] <- 0
  afmeta <- metaprop(R, total, data = Africa)
  # initiate graphic device, create forest plot then save
  png(file = afpath, width = 815, height = 409)
  forest(afmeta, comb.fixed = afmeta$comb.fixed, comb.random = afmeta$comb.random, 
         overall = TRUE, studlab = Africa$Country, plotwidth=unit(10, "cm"), xlim = c(0, 1), colgap.forest = "2cm")
  dev.off()
  
  #oceania
  Oceania <- complete |>
    filter(Continent %in% "OCEANIA") |>
    group_by(Country, !!sym((i))) |>
    tally() |>
    pivot_wider(names_from = !!sym((i)), values_from = n)
  # define path to save created png
  ocepath <- file.path("/path", "to", "directory", paste(newdir, "/oce forest.png", sep = ""), fsep = "/")
  # filter out R and perform metaprop analysis, if there are no R isolates set to 0 then run metaprop
  if(!"R" %in% names(Oceania)){Oceania$R <- 0
  Oceania$total <- rowSums(Oceania[, !colnames(Oceania) %in% "Country"], na.rm = TRUE)
  Oceania[is.na(Oceania)] <- 0
  ocemeta <- metaprop(R, total, data = Oceania) 
  }
  
  else Oceania$total <- rowSums(Oceania[, !colnames(Oceania) %in% "Country"], na.rm = TRUE)
  Oceania[is.na(Oceania)] <- 0
  ocemeta <- metaprop(R, total, data = Oceania)
  # initiate graphic device, create forest plot then save
  png(file = ocepath, width = 765, height = 169)
  forest(ocemeta, comb.fixed = ocemeta$comb.fixed, comb.random = ocemeta$comb.random, 
         overall = TRUE, studlab = Oceania$Country, plotwidth=unit(10, "cm"), xlim = c(0, 1), colgap.forest = "2cm")
  dev.off()
  
}