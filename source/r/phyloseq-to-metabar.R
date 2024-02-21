# Install/load packages

# Function to install a package if not already installed
install_if_missing <- function(package_name) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    install.packages(package_name, dependencies = TRUE)
  }
}

# List of packages to load
packages_to_load <- c("igraph", "ggplot2", "dplyr", "tidyr",
                      "biomformat", "metabaR", "phyloseq",
                      "microbiome", "biom", "data.table")

# Load packages
for (package_name in packages_to_load) {
  library(package_name, character.only = TRUE)
}

################################################################################
# Function to download and source R script
download_and_source <- function(url, destination) {
  if (!file.exists(destination)) {
    download.file(url, destination, mode = "wb")
  }
  source(destination)
}

# Download and source microfiltR_source_code.R
download_and_source(
  "https://raw.githubusercontent.com/itsmisterbrown/microfiltR/master/microfiltR_source_code.R",
  "microfiltR_source_code.R")

# Download and source check_metabarlist.R
download_and_source(
  "https://raw.githubusercontent.com/metabaRfactory/metabaR/HEAD/R/check_metabarlist.R",
  "check_metabarlist.R")

# Download and source subset_metabarlist.R
download_and_source(
  "https://raw.githubusercontent.com/metabaRfactory/metabaR/HEAD/R/subset_metabarlist.R",
  "subset_metabarlist.R")

# Download and source ggtaxplot.R
download_and_source(
  "https://raw.githubusercontent.com/slambrechts/metabaR/master/R/ggtaxplot.R",
  "ggtaxplot.R")

################################################################################
################################################################################
# Load in Phyloseq Olig01 Annelida MBAG
load("./physeq_Olig01_Annelida_rar_species.Rdata")


fp <- paste(getwd(), "/", sep = "")
write.dataset(
  ps = physeq_Olig01_Annelida_rar.rarefied.species,
  filePATH = fp,
  filePREFIX = "Olig01_Annelida")


################################################################################
# File reads
asv_tab <- as.matrix(
  t(
    phyloseq::otu_table(physeq_Olig01_Annelida_rar.rarefied.species)
    )
  )
file_reads <- as.data.frame(asv_tab)

# Make numeric columns numeric
str(file_reads)

file_reads_num <- file_reads
file_reads_num[] <- as.data.frame(sapply(file_reads, as.numeric))
str(file_reads_num)

# create matrix
file_reads_num_matrix <- as.matrix(file_reads_num)
class(file_reads_num_matrix)

# Add _r1 to the row names
# adding suffix to row names
rownames(file_reads_num_matrix) <- paste(
  rownames(file_reads_num_matrix),
  "r1",
  sep = "_"
  )


# Save file
write.table(
  file_reads_num_matrix,
  file = "Olig01_Annelida_ASV_table_transpose.txt",
  sep = "\t",
  row.names = TRUE)


################################################################################
#file pcrs (moet eerst nog  aangemaakt worden)


# Assuming your data is stored in a matrix or a data frame named 'your_data'
# First, create a dataframe with the first two columns
df_pcr <- data.frame(
  pcr_id = rownames(file_reads_num_matrix),
  sample_id = sub("_r1$", "", rownames(file_reads_num_matrix)),
  stringsAsFactors = FALSE
)

# Add the last two columns
df_pcr$type <- "sample"
df_pcr$control_type <- NA

# create matrix
file_pcrs_matrix <- as.matrix(df_pcr)
class(file_pcrs_matrix)

# Save file
write.table(
  file_pcrs_matrix,
  file = "Olig01_Annelida_PCRs.txt",
  sep = "\t",
  row.names = FALSE)

################################################################################
# file samples

file_samples <- suppressWarnings(
  as.matrix(
    phyloseq::sample_data(physeq_Olig01_Annelida_rar.rarefied.species)
    )
  )

# create matrix
file_samples_matrix <- as.matrix(file_samples)
class(file_samples_matrix)

# Save file
write.table(
  file_samples_matrix,
  file = "Olig01_Annelida_sample_data2.txt",
  sep = "\t",
  row.names = TRUE)

################################################################################
# File motus

# Belangrijk! De sequenties moeten op 1 lijn staan per OTU in de sequentie file (olig01_otus92_mbag.fa)


# Read taxonomy table from physeq
tax.tab <- as.data.frame(
  phyloseq::tax_table(physeq_Olig01_Annelida_rar.rarefied.species)
  )

# Create a taxonomy dictionary
taxonomy_dict <- list()
for (i in 1:nrow(tax.tab)) {
  otu <- rownames(tax.tab)[i]
  taxonomy <- paste(na.omit(tax.tab[i, ]), collapse = ";")
  taxonomy_dict[[otu]] <- taxonomy
}

# Process the input file and create the output file
input_file <- file("olig01_otus92_mbag.fa", "r")
output_file <- file("Olig01_Annelida_motus_seq.txt", "w")

# Add new column names
new_column_names <- c("#ASVID", "taxonomy", "sequence")

# Write column names to the output file
writeLines(paste(new_column_names, collapse = "\t"), output_file)

# Create an empty dataframe to store the data
output_df <- data.frame(matrix(ncol = length(new_column_names), nrow = 0))
colnames(output_df) <- new_column_names

current_otu <- NULL
current_sequence <- NULL

while (length(line <- readLines(input_file, n = 1)) > 0) {
  if (substr(line, 1, 1) == ">") {
    if (!is.null(current_otu)) {
      if (current_otu %in% names(taxonomy_dict) && !is.null(current_sequence)) {
        # Write to file
        writeLines(paste(current_otu, taxonomy_dict[[current_otu]], current_sequence, sep = "\t"), output_file)
        # Add to dataframe
        new_row <- c(current_otu, taxonomy_dict[[current_otu]], current_sequence)
        output_df <- rbind(output_df, new_row)
      }
    }
    current_otu <- gsub(">", "", trimws(line))
    current_sequence <- NULL
  } else {
    current_sequence <- paste(current_sequence, trimws(line), sep = "")
  }
}

# Add the last record if it's in the taxonomy file and has a sequence
if (!is.null(current_otu) &&
    current_otu %in% names(taxonomy_dict) &&
    !is.null(current_sequence)) {
  # Write to file
  writeLines(
    paste(
      current_otu,
      taxonomy_dict[[current_otu]],
      current_sequence,
      sep = "\t"),
    output_file
  )
  # Add to dataframe
  new_row <- c(current_otu, taxonomy_dict[[current_otu]], current_sequence)
  output_df <- rbind(output_df, new_row)
}

close(input_file)
close(output_file)

# Assigning new column names to output_df
colnames(output_df) <- new_column_names


# Now you can use 'output_df' as your dataframe for further processing.

class(output_df)

file_motus <- as.data.frame(output_df)


#redefine row names
rownames(file_motus) <- file_motus[, 1]
file_motus <- file_motus[ ,-1]


# Reorder the OTU_ names
# Extract numeric part from row names
row_numbers <- as.numeric(gsub("OTU_", "", rownames(file_motus)))

# Sort row names based on numeric values
sorted_row_names <- rownames(file_motus)[order(row_numbers)]

# Reorder the data frame using sorted row names
sorted_data_motus <- file_motus[sorted_row_names, ]


sorted_data_motus <- sorted_data_motus %>%
  separate(col = taxonomy, sep = ";", into = paste0('Tax', 1:6), fill = 'right')


# Rename columns
rownames(sorted_data_motus) <- rownames(sorted_data_motus)
colnames(sorted_data_motus) <- c(
  "phylum_name", "class_name", "order_name", "family_name", "genus_name",
  "species_name", "sequence")

# Rename to NA were needed Wanneer een OTU niet geïdentificeerd is op een bepaald niveau (manueel in google sheet)

motus_NA <- sorted_data_motus

# Function to replace names ending with '_otu' followed by digits with NA
replace_names <- function(motus_NA) {
  for (i in 1:ncol(motus_NA)) {
    motus_NA[, i] <- gsub("\\b\\w+_otu\\d+\\b", NA, motus_NA[, i])
  }
  return(motus_NA)
}

# Call the function to replace names
motus_NA <- replace_names(motus_NA)

# Create matrix

motus_NA_matrix <- motus_NA


write.table(
  motus_NA_matrix,
  file = "Olig01_Annelida_motus_NA.txt", sep = "\t",
  row.names = TRUE
)


################################################################################
# Create metabarlist
Olig01_Annelida_metabarlist <- tabfiles_to_metabarlist(
  file_reads = "./Olig01_Annelida_ASV_table_transpose.txt",
  file_motus = "Olig01_Annelida_motus_NA.txt",
  file_pcrs = "./Olig01_Annelida_PCRs.txt",
  file_samples = "./Olig01_Annelida_sample_data2.txt",
  files_sep = "\t")

summary_metabarlist(Olig01_Annelida_metabarlist)

save(Olig01_Annelida_metabarlist,file = "metabarlist_Olig01_Annelida.Rdata")


#################################################################################
#################################################################################
# VISUALISATIE

# Specify what you need
indices <- grepl("Annelida", Olig01_Annelida_metabarlist$motus$phylum_name)

Olig01_Annelida <- subset_metabarlist(
  Olig01_Annelida_metabarlist,
  table = "motus",
  indices = indices
)



# Using a taxonomic table
taxo.col <- c(
  "phylum_name", "class_name", "order_name",
  "family_name", "genus_name", "species_name"
)
plot <- ggtaxplot(Olig01_Annelida, taxo.col)
plot
ggsave("my_plot.pdf", plot = plot, width = 18, height = 10, units = "in")

ggsave("my_plot.png", plot = plot, width = 18, height = 10, units = "in")
