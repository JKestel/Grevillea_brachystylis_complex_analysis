###################################################################################################

# Frankenstein_Ver.1.3 _ A looping script for population-level SNP & Autosomal diversity analysis

###################################################################################################
#Ver.1.0 created by Joshua Kestel with help from Ben Anderson  - April 2024
#Ver.1.0 adapted from material provided by Rachel Binks
#Ver.1.1 modified by Rachel Binks - May 2024
#Ver.1.2 modified by Joshua Kestel - June 2024
#Ver.1.3 modified by Joshua Kestel - September 2024
###################################################################################################
 
                                      ###Function of this script###

# The Frankenstein script subsets data by population, randomly subsets populations by a set number of 
# individuals and calculates both SNP and Autosomal diversity for each of the specified populations
# in your indmetrics file.

# Prior to running Frankenstein clones/dodgy samples need to be removed from each population.

##################################################################################################
##################################################################################################


# load in packages 
{library(dartR)
gl.install.vanilla.dartR()
library(ape)
library(plotly)
library(StAMPP)
library(hierfstat)
library(directlabels)
library(pheatmap)
library(rrBLUP)
library(plotrix)
library(graph4lg)
library(gdata)
library(poppr)
library(RColorBrewer)
library(gplots)
library(adegenet)
library(pbapply)}
library(readr)
library(dplyr)

##################################################################################################
##### 1. REMOVE LOW-QUALITY INDIVIDUALS  #####


# Read the data - Filtered metadata 2 has new pops
Ind_calls <- gl.read.dart(filename="RAW_data.csv", ind.metafile="filtered_metadata2.csv")

Ind_calls
?gl.filter.callrate

# Remove low quality individuals from the dataset
filtered <- gl.filter.callrate(Ind_calls, "ind", threshold=0.20)

# Convert the report to a dataframe
filtered <- as.data.frame(filtered$other$ind.metrics)

# Write the filtered data to a new CSV file
write.csv(filtered, file = "filtered_output.csv")


##################################################################################################
##### 2. SET INPUT VALUES  #####

# load metadata file 

metadata <- read_csv("filtered_output.csv")

# check object summary to determine lowest pop size and define the critical threshold needed

# Input

critical_threshold <- 4

##################################################################################################
##### 3. STANDARDISE THE POPULATIONS #####

# Convert the data frame to a list
data_list <- split(metadata, metadata$pop)

# Initialize an empty list to store the modified populations
standardised_list <- list()

# Initialize an empty list to store the deleted individuals that exceeded the critical_threshold
deleted_individuals <- list()

#Initialize an empty list to store the omitted individuals from populations under the critical_threshold
omitted_individuals <- list()

# Loop over each population
for(i in seq_along(data_list)) {
  # Subset the data for the current population
  gl <- data_list[[i]]
  
  pop_name <- names(data_list)[i]
  
  if(!is.null(pop_name)) {
    # Get the population count
    pop_count <- length(gl$pop)
    
    # get the identifiers of all individuals
    all_individuals <- gl$id
    
    # Print the population count and all individuals
    print(paste("Population count for", pop_name, ":", pop_count))
    print(paste("All individuals in", pop_name, ":", toString(all_individuals)))
    
    # check if the population count exceeds the critical threshold
    if(pop_count > critical_threshold) {
      # calculate the number of individuals to drop
      num_drop <- min(pop_count - critical_threshold, length(all_individuals))
      
      # get the identifiers of the individuals to drop
      ind_drop <- sample(all_individuals, size = num_drop)
      
      # Print the number of individuals to drop and their identifiers
      print(paste("Number of individuals to drop from", pop_name, ":", num_drop))
      print(paste("Identifiers of individuals to drop from", pop_name, ":", toString(ind_drop)))
      
      # Add the identifiers of the deleted individuals to the list
      if(pop_name %in% names(deleted_individuals)) {
        deleted_individuals[[pop_name]] <- c(deleted_individuals[[pop_name]], ind_drop)
      } else {
        deleted_individuals[[pop_name]] <- ind_drop
      }
      
      # use subset to remove the individuals
      gl <- subset(gl, !(id %in% ind_drop))
      
      # Add the modified population data to the standardised list
      standardised_list[[pop_name]] <- gl
    } else if (pop_count == critical_threshold) {
      # If the population count is exactly the critical threshold, keep it as is
      standardised_list[[pop_name]] <- gl
    } else {
      # If the population count is below the critical threshold, add the individuals to the omitted_individuals list
      omitted_individuals[[pop_name]] <- all_individuals
    }
  }
}


# check that each pop has been subset to the critical threshold
for(i in seq_along(standardised_list)) {
  # get the current data frame
  df <- standardised_list[[i]]
  
  # get the population name
  pop_name <- names(standardised_list)[i]
  
  # check if the number of individuals is equal to the critical threshold
  if(nrow(df) == critical_threshold) {
    print(paste(pop_name, "has", critical_threshold, "individuals."))
  } else {
    print(paste(pop_name, "does not have", critical_threshold, "individuals. It has", nrow(df), "individuals."))
  }
}


##################################################################################################
##### 4. SAVE CURATED METADATA FILE AND REMOVED INDIVIDUALS FILE #####


# Save standardised populations as a CSV
new_data_df <- bind_rows(standardised_list)
write_csv(new_data_df, "Curated_populations.csv")

# Create a data frame for deleted individuals
deleted_individuals_df <- data.frame(
  id = unlist(deleted_individuals),
  reason = rep("randomly removed", length(unlist(deleted_individuals)))
)

# Create a data frame for omitted individuals
omitted_individuals_df <- data.frame(
  id = unlist(omitted_individuals),
  reason = rep("insufficient population size", length(unlist(omitted_individuals)))
)

# Combine the two data frames
removed_individuals_df <- rbind(deleted_individuals_df, omitted_individuals_df)

# Save as a CSV
write_csv(removed_individuals_df, "Removed_individuals.csv")


##################################################################################################
##### 5. LOAD IN DATASET WITH CURATED POPULATIONS #####

# load dataset
input <- gl.read.dart(filename = "RAW_data.csv", ind.metafile = "Curated_populations.csv")

##################################################################################################
##### 6. SET UP POPULATIONS INTO SEPARATE OBJECTS

# get the unique population names
pop_names <- unique(input@pop)

# initialize gl_list as an empty list
gl_list <- list()

# loop over each population name
for(pop_name in pop_names) {
  
  # subset the population with gl.keep.pop
  gl <- gl.keep.pop(input, pop = pop_name, recalc = TRUE)
  
  # add the gl object to gl_list with pop_name as its name
  gl_list[[pop_name]] <- gl
}

# check that gl_list contains all the gl objects with their correct names
# you can access each gl object by its name, like gl_list$AUS1, gl_list$AUS2
gl_list

##################################################################################################
##### 7. SNP DIVERSITY - FILTERING #####

# apply filters to each gl object in gl_list
gl_list_filteredSNP <- lapply(gl_list, function(gl) {
  filter0 <- gl.filter.allna(gl, v=0)
  filter1 <- gl.filter.secondaries(gl, method="best", v=3)
  filter2 <- gl.filter.rdepth(filter1, lower=5, upper=100, v=3)
  filter3 <- gl.filter.reproducibility(filter2, t=0.99, v=3)
  filter4 <- gl.filter.callrate(filter3, method="loc", threshold=1, v=3)
  filter5 <- gl.filter.maf(filter4, threshold=0.02, verbose=3)
  filter6 <- gl.filter.monomorphs(filter5, v=3)
  filter7 <- gl.filter.callrate(filter6, method="ind", threshold=0.8, v=3)
  gl.recalc.metrics(filter7)
  return(filter7)
})

gl_list_filteredSNP


# check that each pop has retained the critical threshold number of individuals post-filtering
for(i in seq_along(gl_list_filteredSNP)) {
  # get the current gl object
  gl <- gl_list_filteredSNP[[i]]
  
  # get the population name
  pop_name <- names(gl_list_filteredSNP)[i]
  
  # check if the number of individuals is equal to the critical threshold
  if(length(gl@ind.names) == critical_threshold) {
    print(paste(pop_name, "has", critical_threshold, "individuals."))
  } else {
    print(paste(pop_name, "does not have", critical_threshold, "individuals. It has", length(gl@ind.names), "individuals."))
  }
}


##################################################################################################
##### 8. SNP DIVERSITY - ANALYSIS                
                           

# calculate your SNP diversity measures on all standardised populations
het_stats_pop <- lapply(gl_list_filteredSNP, function(gl) {
  gl.report.heterozygosity(gl)
})

# initialize an empty list to store the extracted data
het_standard <- list()

# extract the desired columns from each dataframe in het_stats_pop
for(i in seq_along(het_stats_pop)) {
  het_standard[[i]] <- het_stats_pop[[i]][, c("pop", "nInd", "nLoc", "polyLoc", "Ho", "HoSD", "He", "HeSD", "FIS")]
}

# combine the extracted data into a single dataframe
het_standard <- do.call(rbind, het_standard)
het_loci <- do.call(rbind, het_stats_pop) #RB: changed from het_stats_loci because that doesn't exist

# calculate metrics (FIS SE and allelic richness will be added in ver_1.1).

poly <- het_standard$polyLoc / het_standard$nLoc * 100 #calculates the percentage of polymorphic loci
het_standard["PCpoly"] <-poly

ho_se <- het_standard$HoSD / sqrt(het_standard$nLoc)  # Calculates the standard error for Ho
het_standard["Ho.SE"] <- ho_se

he_se <- het_standard$HeSD / sqrt(het_standard$nLoc)  # Calculates the standard error for He
het_standard["He.SE"] <- he_se

# re-order the output and save your output

het_standard <- het_standard[, c("pop", "nInd", "nLoc", "polyLoc", "PCpoly", "Ho", "Ho.SE", "He","He.SE", "FIS")]   #re-arrange columns in sensible order
write.csv(het_standard, "SNP_diversity_output.csv")



##################################################################################################
##### 9. AUTOSOMAL DIVERSITY - FILTERING

# apply the filters to each gl object in gl_list_auto
gl_list_filteredAUTO <- lapply(gl_list, function(gl) {
  filter0 <- gl.filter.allna(gl)
  filter2 <- gl.filter.rdepth(filter0, lower=5, upper=100, v=3)
  filter3 <- gl.filter.reproducibility(filter2, t=0.99, v=3)
  filter4 <- gl.filter.callrate(filter3, method="loc", threshold=1, v=3)
  filter5 <- gl.filter.monomorphs(filter4, v=3)
  filter6 <- gl.filter.callrate(filter5, method="ind", threshold=0.80, v=3)
  gl.recalc.metrics(filter6)
  return(filter6)
})

gl_list_filteredAUTO

# check that each pop has been subset to the critical threshold
for(i in seq_along(gl_list_filteredAUTO)) {
  # get the current gl object
  gl <- gl_list_filteredAUTO[[i]]
  
  # get the population name
  pop_name <- names(gl_list_filteredAUTO)[i]
  
  # check if the number of individuals is equal to the critical threshold
  if(length(gl@ind.names) == critical_threshold) {
    print(paste(pop_name, "has", critical_threshold, "individuals."))
  } else {
    print(paste(pop_name, "does not have", critical_threshold, "individuals. It has", length(gl@ind.names), "individuals."))
  }
}


##################################################################################################
##### 10. AUTOSOMAL DIVERSITY - ANALYSIS  

# calculate the number of invariant loci for each population
secondaries <- lapply(gl_list_filteredAUTO, function(gl) {
  gl.report.secondaries(gl)
})

# extract n.invariant for each population
n_invariant_values <- sapply(secondaries, function(stats) {
  stats[7, 2]  # Assuming the value is at row 7, column 2
})

# initialize an empty list to store results
het_adjusted_list <- list()

# loop through each population
for (i in seq_along(gl_list_filteredAUTO)) {
  gl <- gl_list_filteredAUTO[[i]]
  n_variable <- gl@n.loc
  n_invariant <- n_invariant_values[i]  # Population-specific n.invariant
  
  # Calculate heterozygosity using population-specific n.invariant
  het_adjusted <- gl.report.heterozygosity(gl, method = "pop", n.invariant = n_invariant)
  
  # Store the result in the list
  het_adjusted_list[[i]] <- het_adjusted
}

# convert the list to a data frame 
het_adjusted_df <- do.call(rbind, het_adjusted_list)

# Add in the number of invariant sites
het_adjusted_df <- cbind (het_adjusted_df,n_invariant_values)

# calculate the standard error for adjusted observed and expected heterozygosity

ho_adj_se <- het_adjusted_df$Ho.adjSD / sqrt(het_adjusted_df$nLoc + het_adjusted_df$n_invariant_values)  # calculates the standard error for Ho
het_adjusted_df["Ho.adj.se"] <- ho_adj_se

he_adj_se <- het_adjusted_df$He.adjSD/ sqrt(het_adjusted_df$nLoc + het_adjusted_df$n_invariant_values)
het_adjusted_df["He.adj.se"] <- he_adj_se

## removed manual calculation of Wright's measure FIS (F = 1 - (Ho/He)) - ## Now done by DartR

# reorder your columns of interest
het_autosomal <- het_adjusted_df[, c("pop", "nInd", "nLoc", "n_invariant_values", "Ho.adj", "Ho.adj.se", "He.adj","He.adj.se", "FIS")]   #re-arrange columns in sensible order

# save the final output
write.csv(het_autosomal, "AUTO_diversity_output.csv")


##################################################################################################

                                       ### Supplementary Materials ###

                              ### Random subsetting within populations ###

# Populations with a minimum number of individuals (designated the "critical_threshold") will be retained.
# I have set this minimum to 5, but this can be changed depending on your sample size per population.If you are not worried
# about population size minimums, then this value could be set to 1000 and you will keep everything, alternatively, 
# skip this step.

# Note: If you population exceeds the "critical_threshold" of individuals, then each individual is allocated a
# number between 1 and the total number of individuals in that population (designated as "pop_count").
# Then, "pop_count" - "critical_threshold" random numbers are generated, these numbers are designated as "num_drop"
# and this is the number of individuals that have to be removed from this population. Then, the script
# randomly samples "num_drop" individuals from all_individuals and uses gl.drop.ind to remove these individuals and then recalculates
# the metrics for that population. 

# Note: If the population has less than the "critical_threshold" of individuals, this population is 
# dropped by the if(length(gl@ind.names) == critical_threshold) argument and these populations
# will not be carried through for the next part of the analysis.



