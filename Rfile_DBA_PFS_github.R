# Load required libraries
library(reshape2)
library(openxlsx)
library(gemtc)
library(rjags)
library(dplyr)
library(readxl)
library(ggplot2)
library(forcats)
library(multinma)
#upload data
Scenario_PFS_data<-readxl::read_xlsx("Path to the file.xlsx")# Provide the path to your input file here
TPP_PFS_data<-readxl::read_xlsx("Path to the file.xlsx")# Provide the path to your input file here

# Create a new Excel workbook
wb <- createWorkbook()

# Iterate through values of k from 0 to 100
for (i in 1:nrow(TPP_PFS_data)) {

TPP_scenario <- data.frame(diff = TPP_PFS_data$`True target lnhr`[i],std.err=TPP_PFS_data$`Simulated std.error`[i])
  
# Replace the values for Wonderumab-TPP intervention
Scenario_PFS_data[Scenario_PFS_data$treatment == "Wonderumab", c("diff","std.err")] <- TPP_scenario

#################Run the analysis
net_PFS <- set_agd_contrast(Scenario_PFS_data, 
                           study = study,
                           trt = treatment,
                           y = diff, 
                           se = std.err,
                           sample_size = sample_size)


# Fixed effects
set.seed(as.numeric(Sys.time()))#generates a random seed automatically each time
fit_PFS_FE <- nma(net_PFS, 
                 trt_effects = "fixed",
                 prior_intercept = normal(scale = 100),
                 prior_trt = normal(scale = 100),
                 prior_het = half_normal(scale = 5))

rel_PFS_FE=relative_effects(fit_PFS_FE, all_contrasts = TRUE)
rel_PFS_FE=data.frame(rel_PFS_FE)
rel_PFS_FE=rel_PFS_FE[,c(1,2,4,6,10)]
colnames(rel_PFS_FE)=c("first","second","HR","min","max")
rel_PFS_FE$HR <- round(exp(rel_PFS_FE$HR), 2)
rel_PFS_FE$min <- round(exp(rel_PFS_FE$min), 2)
rel_PFS_FE$max <- round(exp(rel_PFS_FE$max), 2)
rel_PFS_FE$label <- paste(rel_PFS_FE$HR, "(", rel_PFS_FE$min, ",", rel_PFS_FE$max, ")")

#################################League table##########################
# Reshape the data
league_table <- dcast(rel_PFS_FE, second ~ first, value.var = "label")
# Set row names
rownames(league_table) <- league_table$second
league_table <- league_table[, -1]
# Set column names
colnames(league_table) <- unique(rel_PFS_FE$first)
# Reorder columns
league_table <- league_table[, rev(colnames(league_table))]

#########################DIC table###################################
dic_result<-dic(fit_PFS_FE)
# Create a data frame with the result
dic_table <- data.frame(
  "Residual deviance" = dic_result$resdev,
  "pD" = dic_result$pd,
  "DIC" = dic_result$dic,
  "data points" = length(dic_result$pointwise$agd_contrast$n_contrast)
)
dic_table

########################sucra ranks###########################
smk_ranks_sucra <- posterior_rank_probs(fit_PFS_FE, lower_better = TRUE, sucra = TRUE, cumulative = TRUE)
sucra_values <- data.frame(smk_ranks_sucra$summary %>% dplyr::select(.trt, sucra))
sucra_values <- sucra_values[order(sucra_values$sucra, decreasing = TRUE),]
sucra_values$sucra <- paste0(round(as.numeric(sucra_values$sucra) * 100, 1), "%")
# Transpose the sucra_values data frame
sucra_table <- t(sucra_values)
sucra_table

######################rank matrix###############################
rank_matrix<-posterior_rank_probs(fit_PFS_FE)$summary
rank_matrix<-rank_matrix %>%
  arrange(desc(.trt))
# Select .trt and probability rank columns
rank_matrix_selected <- rank_matrix %>%
  dplyr::select(.trt, starts_with("p_rank"))
# Rename the columns
colnames(rank_matrix_selected)[-1] <- paste0("Rank_", 1:(ncol(rank_matrix_selected) - 1))
colnames(rank_matrix_selected)[1]<-""


############ Create sheet names#############
sheet_name <- paste("scenario_", i, sep = "") 
# Add a worksheet to the workbook
addWorksheet(wb,sheet_name)
#Write the Wonderumab scenario value
writeData(wb,sheet_name,TPP_scenario,startCol = 10 , rowNames = TRUE)
# Write the data to the worksheet
writeData(wb,sheet_name, league_table, rowNames = TRUE)
# Create a cell style for red font color
red_style <- createStyle(fontColour = "red")
# Create a cell style for red font color
header_style <- createStyle(fgFill = "lightblue")
# Loop through the cells and apply the style if the upper value exceeds 1 or is NA
for (i in 1:nrow(league_table)) {
  for (j in 1:ncol(league_table)) {
    cell_value <- league_table[i, j]
    if (!is.na(cell_value)) {
      # Extract the upper value from the cell value
      upper_value <- as.numeric(gsub(".*\\(\\s*\\d+\\.\\d+\\s*,\\s*(\\d+\\.?\\d*)\\s*\\).*", "\\1", cell_value))
      if (is.na(upper_value) || upper_value >= 1) {
        addStyle(wb, sheet = sheet_name, rows = i + 1, cols = j + 1, style = red_style)
      }
    } else {
      # If cell_value is NA, apply red style
      addStyle(wb, sheet = sheet_name, rows = i + 1, cols = j + 1, style = red_style)
    }
  }
}
# Apply the style to the header cells (first row)
for (i in 1:ncol(league_table)) {
  addStyle(wb, sheet =sheet_name, rows = 1, cols = i+1, style = header_style)
}
# Apply the style to the cells in the first column
for (i in 1:nrow(league_table)) {
  addStyle(wb, sheet =sheet_name, rows = i + 1, cols = 1, style = header_style)
}
# Apply the style to the cells in the first column
for (i in 1:ncol(dic_table)) {
  addStyle(wb, sheet =sheet_name, rows = nrow(league_table) + 5, cols = i+2, style = header_style)
}

# Write the DIC table to the worksheet
writeData(wb,sheet_name, dic_table, startCol = 3, startRow = nrow(league_table) + 5, rowNames = FALSE)
# Apply the style to the cells in the first column
for (i in 1:ncol(sucra_table)) {
  addStyle(wb, sheet =sheet_name, rows = nrow(league_table) + 5 + nrow(dic_table) + 5, cols = i+2, style = header_style)
}

# Write the DIC table to the worksheet
writeData(wb,sheet_name, sucra_table, startCol = 3, startRow = nrow(league_table) + 5+ nrow(dic_table) + 5, rowNames = FALSE, colNames = FALSE)
# Apply the style to the cells in the first column
for (i in 1:ncol(rank_matrix_selected)) {
  addStyle(wb, sheet =sheet_name, rows = nrow(league_table) + 5 + nrow(dic_table) + 5 +nrow(sucra_table)+5, cols = i+3, style = header_style)
}
# Apply the style to the cells in the first column
for (i in 1:nrow(rank_matrix_selected)) {
  addStyle(wb, sheet =sheet_name, rows = nrow(league_table) + 5 + nrow(dic_table) + 5 +nrow(sucra_table)+5+ i, cols = 3, style = header_style)
}
# Write the DIC table to the worksheet
writeData(wb,sheet_name, rank_matrix_selected, startCol = 3, startRow = nrow(league_table) + 5+ nrow(dic_table) + 5+nrow(sucra_table)+5, rowNames = FALSE)
}

# Save the workbook
saveWorkbook(wb,"File path/Result_DBA_PFS.csv", overwrite = TRUE)#Add the folder path to save the DBA result

