#### Packages laden

library(phyloseq)
library(vegan)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(ggbiplot)
library(devtools)
library(microbiome)
library(factoextra)



#### Handig pakket met veel phyloseq add-ons



install.packages("remotes")
remotes::install_github("cpauvert/psadd")

library(psadd)


#### Alfa-diversiteit

########################################
############# True Rarefaction #########
############# vegan::rarefy ############
########################################


library(psadd)

##########################################################################
########### vergeet HIER niet te starten van de NIET-rarefied versie #####################
##########################################################################

# Link: https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html

sample_sums(physeq_Olig01_Annelida) ##sanity check

## Pick whatever dataset you want to work with

physeq <- physeq_Olig01_Annelida_species

# physeq <- physeq_18S_Annelida

##### cleanup ######

physeq <- subset_taxa(physeq, phylum =="Annelida")

sample_sums(physeq)

## Remove samples according to the dataset you want to work with (e.g. without replicates or without a certain land use type)
## metadata staat voorlopig in het nederlands, eerst vertalen naar engels?

# physeq <- subset_samples_no_zero(physeq, Cmon_PlotID != "VUL_IN")



# Sanity check

physeq

head(sample_data(physeq))

physeq_Olig01_Annelida_species

sample_sums(physeq)



######################################################################
############# Voorbereidende stap #########################################
######################################################################


# convert the otu_table() within a phyloseq object to a vegan compatible data object
psotu2veg <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

veganobject <- psotu2veg(physeq)



###################################################################
############ Perform true rarefaction to estiamte species richness ############
##############################################################################

getwd()

setwd(dir = "G:/Gedeelde drives/PRJ_MBAG/4c_bodembiodiversiteit/data/statistiek/")


write.csv(data.frame(sample_data(physeq)), file = "metadata_physeq_Olig01_Annelida.csv")

sam.new <- read.csv("metadata_physeq_Olig01_Annelida.csv", header = T, sep = ",", stringsAsFactors = T)

# Sanity check

head(sam.new)

sam.new


########################################################
############### True rarefaction using vegan ###################
########################################################


Estimated_species_richness <- vegan::rarefy(veganobject, 30780, se = FALSE, MARGIN = 1)

Estimated_species_richness # Sanity check

write.table(Estimated_species_richness, file = "true_rarefaction_Estimated_species_richness_Olig01_Annelida.txt")

head(Estimated_species_richness) # Sanity check

true_rarefaction_results <- cbind(sam.new, Estimated_species_richness)

getwd()

write.csv(true_rarefaction_results, "true_rarefaction_results_Olig01_Annelida_final.csv")


#################################################################################
################################### Plot results #####################################
#################################################################################


#### Land use type ####

library(ggplot2)
library(ggrepel)

ggplot(true_rarefaction_results, aes(x = Landgebruik_MBAG, y = Estimated_species_richness, color = Landgebruik_MBAG)) + geom_point() + ylab("Expected species richness (Annelida(Olig01))") +
  xlab("Landgebruik") + labs(color = "Landgebruik")

####### Eventueel nog opsmukken met: #################

# namen veranderen naar waarden Landgebruik_MBAG

+ scale_color_discrete(name = "Landgebruik", labels = c("Cropland", "Forest", "Grassland", "Recreation")) + scale_x_discrete(name = "Land use type", labels = c("Cropland", "Forest", "Grassland", "Recreation")) + geom_text_repel(aes(label=Naam), size=2)


###### Diepte #######


ggplot(true_rarefaction_results, aes(x = Diepte, y = Estimated_species_richness, color = Diepte)) + geom_point() + ylab("Expected species richness (Annelida(Olig01))") +
  xlab("Diepte") + labs(color = "Diepte")


#####################################################################
########################   Statistics     ###################################
#####################################################################


shapiro.test(true_rarefaction_results$Estimated_species_richness)

# is p-value van de shapiro test hierboven kleiner dan 0.05, dan is de data niet normaal verdeeld en doen we kruskal-wallis

library(rstatix)

kruskal_test(true_rarefaction_results, Estimated_species_richness ~ Landgebruik_MBAG)
kruskal_test(true_rarefaction_results, Estimated_species_richness ~ Diepte)





#### alfa-diversiteit 2: Kruskal Wallis en Dunn post-hoc test with adjusted p-values using the Benjamini–Hochberg method

# Eerst beslissen voor hoeveel dingen je wil testen. Welke zijn dus de goede variabelen in de metadata

install.packages("dunn.test")
library(dunn.test)
library(rstatix)


shapiro.test(true_rarefaction_results$Estimated_species_richness)

# is p-value van de shapiro test hierboven kleiner dan 0.05, dan is de data niet normaal verdeeld en doen we kruskal-wallis


kruskal_test(true_rarefaction_results, Estimated_species_richness ~ Landgebruik_MBAG)
kruskal_test(true_rarefaction_results, Estimated_species_richness ~ Diepte)


# Define the variable for which you found a significant difference
variable_of_interest <- "Landgebruik_MBAG"  # Replace with the variable name

# Perform the Kruskal-Wallis test
kruskal_result <- kruskal_test(true_rarefaction_results, Estimated_species_richness ~ get(variable_of_interest))



# Perform the Dunn post-hoc test with Benjamini-Hochberg correction
dunn_result <- dunn.test(true_rarefaction_results$Estimated_species_richness, data[[variable_of_interest]], method = "bh")

# Display the Dunn test results
print("Dunn Test Results:")
print(dunn_result)

# Extract adjusted p-values
adjusted_p_values <- dunn_result$P.adjusted

# Print the adjusted p-values for each pairwise comparison
print("Adjusted P-values:")
print(adjusted_p_values)


Identify Groups with Highest/Lowest Alpha Diversity: You can also extract specific values if you want to determine which group has the highest or lowest estimated species richness. For example, you can use the summary function to compute summary statistics and identify the group with the maximum or minimum mean estimated species richness:
  

# Find the group with the highest estimated species richness
max_group <- summary(true_rarefaction_results$Estimated_species_richness[true_rarefaction_results$Landgebruik_MBAG == "group_name"])$mean

# Find the group with the lowest estimated species richness
min_group <- summary(true_rarefaction_results$Estimated_species_richness[true_rarefaction_results$Landgebruik_MBAG == "group_name"])$mean

cat("Group with the highest estimated species richness:", max_group, "\n")
cat("Group with the lowest estimated species richness:", min_group, "\n")



Replace "group_name" with the actual group you want to investigate. This code calculates the mean estimated species richness for the specified group.




#### rarefaction curves



library(vegan)

S <- specnumber(veganobject)

getwd()

write.table(S, file = "Olig01_Annelida_specnumber.txt")

specnumber(veganobject)

min(rowSums(veganobject))

raremax <- min(rowSums(veganobject))

options(scipen=5)

rarecurve(veganobject, step = 100, sample = raremax, col = "blue", cex= 0.6)





#### Phyloseq shiny for quick visualisations

install.packages("shiny")
library(shiny)
shiny::runGitHub("shiny-phyloseq","joey711")

#Als hij vraagt om Rstudio te rebooten, NIET DOEN


#### General results from phyloseq object

OTU counts per phylum / class / order / family / genus


# Load required packages
library(phyloseq)
library(dplyr)

# Count OTUs per Phylum
otu_counts_per_phylum <- tax_table(physeq) %>%
  as.data.frame() %>%
  group_by(Phylum) %>%
  summarise(OTU_Count = n())

# Print the result
print(otu_counts_per_phylum)

write.csv(otu_counts_per_phylum, file = "otu_counts_per_phylum.csv")

## Genus level

# Load required packages
library(phyloseq)
library(dplyr)

# Count OTUs per Genus
otu_counts_per_genus <- tax_table(physeq) %>%
  as.data.frame() %>%
  group_by(Genus) %>%
  summarise(OTU_Count = n())

# Print the result
print(otu_counts_per_genus)

write.csv(otu_counts_per_genus, file = "otu_counts_per_genus.csv")

## Species level

# Load required packages
library(phyloseq)
library(dplyr)

# Count OTUs per Species
otu_counts_per_species <- tax_table(physeq) %>%
  as.data.frame() %>%
  group_by(Species) %>%
  summarise(OTU_Count = n())

# Print the result
print(otu_counts_per_species)

write.csv(otu_counts_per_species, file = "otu_counts_per_species.csv")


##Read counts per phylum   #### werkt miss niet

# Summarize read counts per Phylum
physeq %>%
  tax_glom(taxrank = "Phylum") %>%
  transform_sample_counts(function(x) sum(x)) %>%
  otu_table() %>%
  as.data.frame() %>%
  group_by(Phylum) %>%
  summarise(Read_Count = sum(!!sym(colnames(.))))

# Print the result
print(otu_counts_per_phylum)







#### How to remove unwanted OTUs from phyloseq object

pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(allTaxa, physeq))
}

badTaxa = c("Zotu511")
physeq <- pop_taxa(physeq, badTaxa)

physeq



################################################################
#############   Ordinatie   ####################################



library(phyloseq)
library(vegan)
library(compositions)
library(ggplot2)

# Load phyloseq object
physeq <- physeq_Olig01_Annelida_rar.rarefied.species

physeq

#physeq <- rarefy_even_depth(physeq, sample.size = 30780)

physeq

getwd()

head(sample_data(physeq))

setwd(dir = "G:/Gedeelde drives/PRJ_MBAG/4c_bodembiodiversiteit/data/statistiek/")

getwd()

save(physeq, file = "physeq_MBAG_Olig01_Annelida_30780.RData")





# Load metadata

write.csv(data.frame(sample_data(physeq)), file = "metadata_physeq_Olig01_Annelida.csv")

your_metadata <- read.csv("metadata_physeq_Olig01_Annelida.csv", header = T, sep = ",", stringsAsFactors = T)

## OF ZONDER stringsAsFactors = T ???

# Perform CLR transformation
clr_physeq <- transform(physeq, "clr")

# Calculate Bray distance matrix based on the CLR-transformed data
bray_distance <- distance(physeq, method = "bray")

# Perform PCoA
pcoa_result <- ordinate(physeq, method = "PCoA", distance = "bray")

# Create the ordination plot with color based on land_use and shape based on diepte
plot_ordination(physeq, pcoa_result, type = "samples", color = "Landgebruik_MBAG", shape = "Diepte") +
  theme_minimal()




