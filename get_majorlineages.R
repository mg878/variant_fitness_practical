## This is just a tiny bit of code to convert pangolin calls to major SARS-CoV-2 clades
## Original data is from the ONS-CIS lineage replacement paper (https://royalsocietypublishing.org/doi/10.1098/rspb.2023.1284) 
## and lives here: https://github.com/thomasallanhouse/covid19-lineages/tree/main?tab=MIT-1-ov-file


library(tidyverse)
# Read the dataset from the URL
url <- "https://raw.githubusercontent.com/thomasallanhouse/covid19-lineages/refs/heads/main/ONS_day_lineage.csv"
lineage_data <- read_csv(url)

# Define broader variant categories
alpha_variants <- c("B.1.1.7", "Q.4", "Q.8", "Q.1")
delta_variants <- c("AY.4", "AY.5", "AY.120", "B.1.617.2", "AY.43", "AY.9", "AY.4.2", "AY.6", "AY.7", "AY.126",
                    "AY.4.5", "AY.98", "AY.111", "AY.87", "AY.122", "AY.3", "AY.46.5", "AY.11", "AY.42", "AY.4.1",
                    "AY.46", "AY.98.1", "AY.9.2", "AY.4.13", "AY.90", "AY.46.6", "AY.36", "AY.34", "AY.4.2.1",
                    "AY.43.6", "AY.4.8", "AY.37", "AY.25.1", "AY.4.2.2", "AY.4.3", "AY.80", "AY.4.11", "AY.5.4",
                    "AY.26", "AY.93", "AY.55", "AY.44", "AY.33", "AY.25", "AY.74", "AY.125", "AY.121", "AY.79",
                    "AY.43.3", "AY.4.7", "AY.105", "AY.45", "AY.46.2", "AY.106", "AY.20", "AY.39.1.1", "AY.124",
                    "AY.24", "AY.58", "AY.34.1", "AY.7.2", "AY.4.9", "AY.70", "AY.107", "AY.43.2", "AY.39.1",
                    "AY.129", "AY.5.5", "AY.43.4", "AY.127", "AY.23", "AY.47", "AY.91", "AY.4.2.3", "AY.4.2.5",
                    "AY.34.2", "AY.4.14", "AY.120.1", "AY.103", "AY.121.1", "AY.108", "AY.120.2.1", "AY.109",
                    "AY.7.1", "AY.110", "AY.119", "AY.118", "AY.114", "AY.128", "AY.100", "AY.4.10", "AY.36.1",
                    "AY.124.1.1", "AY.122.5", "AY.95", "AY.122.6", "AY.4.2.4", "AY.39", "AY.43.8", "AY.112.3",
                    "AY.29.1", "AY.116.1", "AY.104", "AY.122.4", "AY.112", "AY.116", "AY.46.4", "AY.112.2",
                    "AY.122.1", "AY.4.6", "AY.85", "AY.61", "AY.94", "AY.39.1.4", "AY.33.2", "AY.20.1", "AY.8",
                    "AY.75", "AY.102", "AY.92", "AY.99", "AY.78", "AY.16", "AY.96", "AY.4.17", "AY.10", "AY.25.3")
ba1_variants <- c("BA.1", "BA.1.17.2", "BA.1.15.1", "BA.1.15", "BA.1.10", "BA.1.16", 
                  "BA.1.9", "BA.1.1", "BA.1.7", "BA.1.1.15", "BA.1.17", "BA.1.1.13", 
                  "BA.1.20", "BA.1.5", "BA.1.1.14", "BA.1.18", "BA.1.8", "BA.1.14", 
                  "BA.1.12", "BA.1.1.1", "BA.1.13", "BA.1.1.12", "BA.1.1.11", 
                  "BA.1.1.18", "BA.1.1.4", "BA.1.21", "BA.1.2", "BA.1.14.1", "BA.1.19", 
                  "BA.1.1.10", "BA.1.3", "BA.1.1.2", "BA.1.21.1", "BA.1.1.7", "BA.1.6", 
                  "BA.1.4", "BA.1.14.2", "BA.1.1.9", "BA.1.17.1", "BA.1.1.16", 
                  "BA.1.22", "BA.1.1.8", "BA.1.1.3", "BD.1", "BA.1.15.2", "BA.1.15.3")
ba2_variants <- c("BA.2.9", "BA.2", "BA.2.10", "BA.2.37", "BA.2.3", "BA.2.38", 
                  "BA.2.1", "BA.2.5", "BA.2.29", "BA.2.7", "BA.2.16", "BA.2.10.1", 
                  "BA.2.36", "BA.2.24", "BA.2.23", "BA.2.31", "BA.2.33", "BA.2.12", 
                  "BA.2.8", "BA.2.2", "BA.2.18", "BA.2.41", "BA.2.22", "BA.2.26", 
                  "BA.2.23.1", "BA.2.17", "BA.2.39", "BA.2.25", "BA.2.9.1", "BA.2.6", 
                  "BA.2.3.1", "BA.2.4", "BA.2.34", "BA.2.32", "BA.2.15", "BA.2.14", 
                  "BA.2.9.5", "BA.2.48", "BA.2.56", "BA.2.58", "BA.2.52", "BA.2.74",
                  "BA.2.53", "BA.2.47", "BA.2.65", "BA.2.51", "BA.2.49", "BA.2.50", 
                  "BA.2.9.2", "BA.2.3.5", "BA.2.69", "BA.2.3.2", "BA.2.62", 
                  "BA.2.3.13", "BA.2.67", "BA.2.27", "BA.2.3.9", "BA.2.9.3", "BA.2.13", 
                  "BA.2.31.1", "BA.2.46", "BA.2.64", "BA.2.55", "BA.2.10.3", "BA.2.72", 
                  "BA.2.3.16", "BA.2.54", "BA.2.12.1", "BA.2.3.11", "BA.2.44", 
                  "BA.2.57", "BA.2.68", "BA.2.79", "BA.2.21", "BA.2.3.15", "BA.2.81", 
                  "BA.2.45", "BA.2.63", "BA.2.11", "BA.2.40.1", "BA.2.25.1", 
                  "BA.2.3.22", "BA.2.60", "BA.2.3.6", "BA.2.30", "BA.2.3.4", "BA.2.61", 
                  "BA.2.70", "BG.2", "BA.2.71", "BA.2.76", "BA.2.38.3", "BA.2.3.8", 
                  "BA.2.13.1", "BA.2.78", "BA.2.35", "BA.2.80", "BA.2.9.7", "BG.5", 
                  "BA.2.28", "BA.2.9.6", "BA.2.83", "BA.2.73", "BA.2.38.2", "BJ.1",
                  "BA.2.42", "CM.2", "BS.1", "BA.2.3.20", "BS.1.1", "CM.4.1", "CM.3", 
                  "CM.12", "CM.7", "CM.8.1", "CM.11")
ba275_variants <- c("BA.2.75", "BA.2.75.1", "BY.1", "BL.1", "BA.2.75.2", "BA.2.75.3", 
                    "BN.1.5", "BL.1.1", "BA.2.75.10", "BA.2.75.4", "CB.1", "BM.4.1.1", 
                    "BN.1.3", "BN.1", "BA.2.75.5", "BM.1.1.3", "BN.6", "BN.5", 
                    "BM.1.1.1", "CA.3", "BY.1.2", "BL.2", "BM.1.1", "BN.1.2.1", "BR.1", 
                    "BN.1.3.1", "BL.1.4", "DS.1", "CA.1", "BM.5", "BR.3", "BN.1.4", 
                    "CA.7", "BN.3", "BN.1.9", "CH.1.1", "CV.1", "BN.1.7", "BN.2", "BM.1", 
                    "BN.1.2", "BM.4.1", "BA.2.75.7", "BN.4", "CA.2", "BA.2.75.6", 
                    "BN.1.4.1", "CH.1.1.2", "BN.1.1", "CH.1.1.5", "BN.2.1", "BR.2.1", 
                    "CJ.1", "BN.3.1", "BY.1.1.1", "BR.2", "CH.1.1.1", "CA.5", "BM.1.1.5", 
                    "CA.6", "BN.1.8", "CV.2", "CH.3.1", "BN.1.1.1", "BR.1.2", "BM.1.1.4", 
                    "CH.1.1.3", "CJ.1.1", "BM.2.1", "CA.3.1", "DV.1", "CH.1.1.6", 
                    "CH.1.1.8", "DV.3", "CH.1.1.9", "BM.2", "DV.2", "CH.1.1.7", "BN.1.10")
ba4_variants <- c("BA.4", "BA.4.1.1", "BA.4.1", "BA.4.6", "BA.4.1.9", "BA.4.4", 
                  "BA.4.2", "BA.4.7", "BA.4.6.5", "BA.4.1.6", "BA.4.5", "BA.4.6.1", 
                  "BA.4.1.8", "BA.4.3", "BA.4.1.4", "BA.4.1.5", "BA.4.6.4", "BA.4.8", 
                  "BA.4.1.10", "DC.1", "BA.4.6.2", "CS.1", "BA.4.6.3")
ba5_variants <- c("BA.5.2.28", "BE.1", "BA.5.2.1", "BF.1", "BA.5.3.1", "BA.5", "BE.5", 
                  "BA.5.1", "BA.5.3.3", "BA.5.1.22", "BA.5.9", "BA.5.2.3", "BA.5.8", 
                  "BA.5.2", "BA.5.1.23", "BA.5.5", "BF.4", "BA.5.3.2", "BA.5.3", 
                  "BA.5.1.19", "BA.5.2.26", "BF.6", "BE.1.1", "BF.28", "BA.5.1.5", 
                  "BA.5.1.25", "BF.10", "BE.1.4", "BA.5.1.10", "BA.5.1.3", "BA.5.1.12", 
                  "BA.5.6", "BF.5", "BA.5.2.33", "BF.18", "BA.5.1.2", "BE.3", 
                  "BA.5.1.4", "BA.5.2.24", "BF.27", "BA.5.2.8", "BA.5.1.30", 
                  "BA.5.1.11", "BA.5.3.4", "BF.31", "BA.5.1.9", "BA.5.1.24", "BE.6", 
                  "BA.5.2.20", "BF.15", "BA.5.1.1", "BA.5.2.9", "BA.5.2.16", "BF.7", 
                  "BA.5.2.2", "BA.5.5.1", "BF.11", "BA.5.2.21", "BE.1.3", "BA.5.1.16", 
                  "BF.9", "BA.5.1.15", "BF.25", "BE.1.4.4", "BF.2", "BE.1.1.2", "DE.1", 
                  "BA.5.5.2", "BE.1.1.1", "BA.5.2.4", "BA.5.2.27", "BA.5.2.19", "BT.1", 
                  "BA.5.2.13", "BE.7", "BF.14", "BA.5.2.6", "BF.12", "BA.5.1.8", 
                  "BF.8", "BA.5.2.7", "BA.5.1.26", "BA.5.1.6", "BA.5.6.2", "BF.17", 
                  "BT.2", "BE.1.2.1", "BF.7.6", "BE.4", "BA.5.1.17", "BK.1", "BF.7.12", 
                  "BA.5.2.18", "BF.20", "BF.26", "BF.21", "BA.5.2.22", "BF.7.4", 
                  "BF.1.1", "BF.32", "BA.5.7", "CG.1", "BF.11.1", "BA.5.1.21", 
                  "BA.5.2.12", "CP.1.3", "BF.29", "BF.3", "BA.5.2.36", "BA.5.1.7", 
                  "BF.23", "BA.5.1.27", "CE.1", "BF.3.1", "CK.1", "CP.1", "BA.5.2.31", 
                  "BA.5.3.5", "BA.5.2.10", "BA.5.2.47", "BA.5.2.37", "CK.2", "BA.5.10", 
                  "BF.11.4", "BA.5.6.1", "BF.16", "BE.9", "BA.5.5.3", "BA.5.1.14", 
                  "BF.7.9", "BV.1", "BA.5.2.39", "CN.1", "BA.5.2.29", "BF.7.5", 
                  "BA.5.2.44", "BF.30", "BE.2", "BF.11.3", "BA.5.6.3", "BF.11.5", 
                  "BA.5.1.28", "BF.13", "BA.5.2.35", "CR.2", "BE.1.2", "BA.5.2.30", 
                  "BV.2", "BF.7.4.2", "BE.8", "CP.1.2", "BU.1", "DF.1", "BF.7.7", 
                  "BE.4.1", "CL.1", "DE.2", "BF.7.13.1", "CK.3", "BA.5.2.48", 
                  "BA.5.2.23", "BA.5.2.14", "CT.1", "BF.34", "BE.1.4.1", "BA.5.2.49", 
                  "BF.7.10", "CQ.1.1", "BA.5.2.34", "CR.1", "BA.5.2.25", "CP.1.1", 
                  "CQ.1", "BF.11.2", "BQ.2", "BA.5.2.41", "CU.1", "BA.5.1.31", 
                  "BA.5.1.18", "CC.1", "BF.7.4.1", "CD.1", "BF.7.15", "BA.5.2.40", 
                  "BE.4.1.1", "BF.7.8", "BA.5.11", "BF.5.1", "BE.4.2", "CP.3", 
                  "BA.5.2.38", "CK.2.1.1", "BF.7.13.2", "BF.7.11", "BF.7.1", "CP.5", 
                  "BF.7.3", "BA.5.1.20", "BA.5.6.4", "BA.5.2.43", "BF.7.5.1", "CK.2.1", 
                  "BW.1", "BW.1.1", "DQ.1", "DA.1", "BA.5.2.42", "BF.31.1", 
                  "BA.5.10.1", "DJ.1.1", "CP.6", "BE.10", "CY.1", "DG.1", "DL.1", 
                  "DB.1", "BA.5.2.56", "BF.7.16", "CL.1.1", "BF.7.19", "BF.7.20", 
                  "BA.5.1.33", "BF.40", "BF.7.23", "BA.5.2.59", "BF.7.24")
bq1_variants <- c("BQ.1.1", "BQ.1.5", "BQ.1", "BQ.1.8", "BQ.1.11", "BQ.1.15", 
                  "BQ.1.2", "BQ.1.1.4", "BQ.1.10.1", "BQ.1.1.1", "BQ.1.13", 
                  "BQ.1.1.16", "BQ.1.1.8", "BQ.1.1.2", "BQ.1.1.3", "BQ.1.16", 
                  "BQ.1.18", "BQ.1.20", "BQ.1.1.5", "BQ.1.1.19", "BQ.1.12", "BQ.1.27", 
                  "BQ.1.8.2", "BQ.1.14", "BQ.1.1.6", "BQ.1.10", "BQ.1.3", "BQ.1.4", 
                  "BQ.1.23", "BQ.1.1.18", "BQ.1.21", "BQ.1.1.15", "BQ.1.1.22", 
                  "BQ.1.1.24", "BQ.1.24", "BQ.1.1.10", "BQ.1.19", "BQ.1.25", 
                  "BQ.1.1.7", "DN.1", "BQ.1.7", "BQ.1.1.11", "BQ.1.1.14", "DN.1.1", 
                  "BQ.1.1.23", "DR.1", "BQ.1.1.17", "BQ.1.1.27", "BQ.1.1.29", 
                  "BQ.1.22", "BQ.1.1.13", "BQ.1.13.1", "BQ.1.1.25", "BQ.1.1.31", 
                  "BQ.1.1.28", "BQ.1.28", "DU.1", "BQ.1.6", "BQ.1.1.9", "BQ.1.1.12", 
                  "BQ.1.1.34", "DT.1", "BQ.1.26", "BQ.1.1.32", "BQ.1.1.26", 
                  "BQ.1.1.20", "DM.1", "BQ.1.1.30", "BQ.1.1.21", "DK.1", "BQ.1.25.1", 
                  "BQ.1.26.1", "BQ.1.17", "BQ.1.1.66", "BQ.1.1.44", "EF.1.1", "EF.1", 
                  "BQ.1.1.45", "EF.2", "ED.2", "BQ.1.1.69", "BQ.1.1.57", "BQ.1.1.56", 
                  "BQ.1.1.35", "DN.1.1.1", "EC.1.1", "BQ.1.1.37", "BQ.1.1.53", 
                  "EF.1.2", "BQ.1.1.36", "BQ.1.1.67", "EC.1", "BQ.1.1.47", "BQ.1.1.38", 
                  "BQ.1.1.46", "BQ.1.1.40", "BQ.1.1.43", "BQ.1.1.39", "BQ.1.1.58", 
                  "BQ.1.1.65", "EA.1", "EF.1.3", "DN.1.1.2", "BQ.1.1.41", "EE.4", 
                  "EF.1.1.1", "EA.2", "BQ.1.31", "BQ.1.1.59", "BQ.1.1.63")
xbb_variants <- c("XBB.3", "XBB", "XBB.1", "XBB.1.2", "XBB.1.9", "XBB.2", "XBB.4",
                  "XBB.6", "XBB.1.4.1", "XBB.2.1", "XBB.4.1", "XBB.1.4", "XBB.1.5", 
                  "XBB.2.2", "XBB.3.1", "XBB.1.1", "XBB.3.2", "XBB.1.9.1", "XBB.1.7", 
                  "XBB.1.5.7", "XBB.2.4", "XBB.1.5.4", "XBB.1.5.1", "XBB.1.5.5", 
                  "XBB.1.9.2", "XBB.1.5.6", "XBB.1.5.9", "XBB.2.3", "XBB.1.5.8", 
                  "XBB.1.5.2", "XBB.1.12", "XBB.1.13", "XBB.1.11.1")

# Classify pangolin_call into broader categories
lineage_data <- lineage_data %>%
  mutate(major_lineage = case_when(
    pangolin_call %in% alpha_variants ~ "B.1.1.7",
    pangolin_call %in% delta_variants ~ "B.1.617.2",
    pangolin_call %in% ba1_variants ~ "BA.1",
    pangolin_call %in% ba2_variants ~ "BA.2",
    pangolin_call %in% ba275_variants ~ "BA.2.75",
    pangolin_call %in% ba4_variants ~ "BA.4",
    pangolin_call %in% ba5_variants ~ "BA.5",
    pangolin_call %in% bq1_variants ~ "BQ.1",
    pangolin_call %in% xbb_variants ~ "XBB",
    TRUE ~ "Other" # Default category for any unlisted variants
  ))

View(lineage_data)
# Ensure the pangolin_call column is treated as character to avoid misinterpretation
lineage_data$pangolin_call <- as.character(lineage_data$pangolin_call)

# Fix potential issues (e.g., non-printable characters or unexpected spaces)
lineage_data$pangolin_call <- gsub("[^[:alnum:].]", "", lineage_data$pangolin_call) # Keep only alphanumeric and "."
lineage_data$pangolin_call <- trimws(lineage_data$pangolin_call) # Remove leading/trailing spaces

# Export to CSV with UTF-8 encoding and quote each field to preserve formatting
write.csv(lineage_data, file = "Downloads/independence/Schedule for computing skills HT 2025/lineage_data.csv", row.names = FALSE, fileEncoding = "UTF-8", quote = TRUE)
