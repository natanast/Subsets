


rm(list = ls())
gc()


#create a file "Samples" with all your samples
#set "Sample" file as working directory
#samples folders should contain the output of IMGT (at least the 1_Summary.txt)



library(data.table)
library(dplyr)



#---------------------Subsets------------------------------

#files names

fls    = "C:/Users/nanastasiadou/Desktop/Carmilla 2.0/Subsets/Samples/"

fls    = fls |> 
           list.files(full.names = TRUE) 




for (file in fls) {

    #load file
    summary_file <- file.path(file, "1_Summary.txt")
    df <- fread(summary_file, sep = "\t", fill = TRUE, select = 1:33)

    
    # --------------- subset 1-99---------------------
    df1 <- df |>
            filter(`V-DOMAIN Functionality` == "productive" |
                   `V-DOMAIN Functionality`== "productive (see comment)",
                   `CDR3-IMGT length` != "X") |>
            mutate(`IGHV.subgroup` = substr(`V-GENE and allele`, 8, 12),
                   `IGHJ.gene` = substr(`J-GENE and allele`, 8, 12),
                   `CDR3.offset` = substr(`AA JUNCTION`, 3, 9),
                   `CDR3.offset.true` = grepl("QWL", `CDR3.offset`),
                   `CDR3.offset.Subset1` = substr(`AA JUNCTION`, 5, 7)) |>
            filter(`IGHV.subgroup` %in% c("IGHV1", "IGHV5", "IGHV7"),
                   `IGHJ.gene` == "IGHJ4", 
                   `CDR3-IMGT length` %in% 11:16,
                   `CDR3.offset.true` == "TRUE")|>
            mutate(`is.Subset1` = ifelse(
                   `CDR3-IMGT length` == 13 &
                   `CDR3.offset.Subset1` == "QWL",
                    TRUE,
                    FALSE)) |>
            mutate(is.Subset99 = ifelse(
                  `CDR3-IMGT length` == 14 &
                  `CDR3.offset.Subset1` == "QWL", 
                   TRUE,
                   FALSE))
    
    
    #Subset1
    Subset_1 <- df1 |>
        filter(is.Subset1 == "TRUE") |>
        select(1:(ncol(df1) - 7))
    
    
    
    #Subset99
    Subset99 <- df1 |>
                    filter(is.Subset99 == "TRUE") |>
                    select(1:(ncol(df1)-7))  
                
    
    #Satellite
    Subset1_99_Sat <- df1 |>
                    filter(is.Subset1 == "FALSE",  is.Subset99 == "FALSE") |>
                    select(1:(ncol(df1) - 7))
    
    
    #save file
    df_list <- list(Subset1 = Subset_1, Subset99 = Subset99,
                    Satellite = Subset1_99_Sat)
    
    # Specify the file path
    excel_file <- paste0(file, "/Subset_1_99.xlsx")
    
    # Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    
    
    
    
    # --------------- subset 2-169---------------------
    df1 <- df |>
            filter(`V-DOMAIN Functionality` == "productive" |
                   `V-DOMAIN Functionality` == "productive (see comment)",
                   `CDR3-IMGT length` != "X") |>
            mutate(`IGHV.subgroup` = substr(`V-GENE and allele`, 8, 12),
                   `IGHV.gene` = substr(`V-GENE and allele`, 8, 15),
                   `IGHJ.gene` = substr(`J-GENE and allele`, 8, 12),
                   `CDR3.offset` = substr(`AA JUNCTION`, 2, 6),
                   `CDR3.offsetD.true` = grepl("D", `CDR3.offset`),
                   `CDR3.offsetE.true` = grepl("E", `CDR3.offset`),
                   `CDR3.offset.Subset2` = substr(`AA JUNCTION`, 4, 4)) |>
            filter(`IGHV.subgroup` %in% c("IGHV3"),
                   #`IGHJ.gene` == "IGHJ6", commend out by Andreas
                   `CDR3-IMGT length` %in% 7:11,
                   `CDR3.offsetD.true` == "TRUE" |
                   `CDR3.offsetE.true` == "TRUE") |>
            mutate(`is.Subset2` = ifelse(
                   `CDR3-IMGT length` == 9 &
                  (`CDR3.offset.Subset2` == "D" | 
                   `CDR3.offset.Subset2` == "E") &
                   `IGHV.gene` == "IGHV3-21",
                    TRUE,
                    FALSE)) |>
            mutate(`is.Subset169` = ifelse(
                   `CDR3-IMGT length` == 9 &
                  (`CDR3.offset.Subset2` == "D" | 
                   `CDR3.offset.Subset2` == "E") &
                   `IGHV.gene` == "IGHV3-48",
                    TRUE,
                    FALSE))
    
    
    #Subset2
    Subset_2 <- df1 |>
        filter(`is.Subset2` == "TRUE") |>
        select(1:(ncol(df1) - 9))
    
    
    #Subset169
    Subset169 <- df1 |>
        filter(`is.Subset169` == "TRUE") |>
        select(1:(ncol(df1) - 9))
    
    
    
    #Satellite
    Subset_2_169_Sat <- df1 |>
        filter(`is.Subset2` == "FALSE", `is.Subset169` == "FALSE") |>
        select(1:(ncol(df1) - 8))
    
    
    #save file
    df_list <- list(Subset2 = Subset_2, Subset169 = Subset169, 
                    Satellite = Subset_2_169_Sat)
    
    # Specify the file path
    excel_file <- paste0(file, "/Subset_2_169_without_J.xlsx")
    
    # Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    
    
    
    
    
    # --------------- subset 8 -----------------------------
    df1 <- df |>
            filter(`V-DOMAIN Functionality` == "productive" |
                   `V-DOMAIN Functionality` == "productive (see comment)",
                   `CDR3-IMGT length` != "X") |>
            mutate(`IGHV.subgroup` = substr(`V-GENE and allele`, 8, 12),
                   `IGHV.gene` = substr(`V-GENE and allele`, 8, 15),
                   `IGHJ.gene` = substr(`J-GENE and allele`, 8, 12),
                   `CDR3.offset` = substr(`AA JUNCTION`, 5, 14),
                   `CDR3.offset.true` = grepl("YSSSWY", `CDR3.offset`),
                   `CDR3.offset.Subset8` = substr(`AA JUNCTION`, 7, 12)) |>
            filter(`IGHV.subgroup` %in% c("IGHV2", "IGHV4", "IGHV6"),
                   `IGHJ.gene` == "IGHJ5", 
                   `CDR3-IMGT length` %in% 16:21,
                   `CDR3.offset.true` == "TRUE") |>
            mutate(`is.Subset.8` = ifelse(
                   `CDR3-IMGT length` == 19 &
                   `CDR3.offset.Subset8` == "YSSSWY",
                    TRUE,
                    FALSE)) |>
            mutate(`is.Subset8B` = ifelse(
                   `CDR3-IMGT length` == 18 &
                   `CDR3.offset.Subset8` == "YSSSWY",
                    TRUE,
                    FALSE))
    
    
    #Subset_8
    Subset_8 <- df1 |>
                     filter(`is.Subset.8` == "TRUE") |>
                     select(1:(ncol(df1) - 8))
            
            
    #Subset_8B
    Subset_8B <- df1 |>
                    filter(`is.Subset8B` == "TRUE") |>
                    select(1:(ncol(df1) - 8))
    
    
    
    #Satellite
    
    Subset_8_8B_Sat <- df1 |>
        filter(`is.Subset.8` == "FALSE",
               `is.Subset8B` == "FALSE") |>
        select(1:(ncol(df1) - 7))
    
    
    
    #save file
    df_list <- list(Subset_8 = Subset_8, Subset_8B = Subset_8B,
                    Satellite = Subset_8_8B_Sat)
    
    # Specify the file path
    excel_file <- paste0(file, "/Subset_8_8Β.xlsx")
    
    # Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    
    
    
    # --------------- subset 16 --------------------------------
    
    #filtering data frame
    #the FYCS motif of subset 16 is at positions 4-7, with an offset
    #relaxation of +-2 the relevant positions are 2-9 
    df1 <- df |>
            filter(`V-DOMAIN Functionality` == "productive" |
                   `V-DOMAIN Functionality` == "productive (see comment)",
                   `CDR3-IMGT length` != "X") |>
            mutate(`IGHV.subgroup` = substr(`V-GENE and allele`, 8, 12),
                   `IGHV.gene` = substr(`V-GENE and allele`, 8, 15),
                   `IGHJ.gene` = substr(`J-GENE and allele`, 8, 12),
                   `CDR3.offset` = substr(`AA JUNCTION`, 3, 10),
                   `CDR3.offset.true` = grepl("FYCS", `CDR3.offset`),
                   `CDR3.offset.Subset16` = substr(`AA JUNCTION`, 5, 8)) |>
            filter(`IGHV.subgroup` %in% c("IGHV2", "IGHV4", "IGHV6"),
                   `IGHJ.gene` == "IGHJ6", 
                   `CDR3-IMGT length` %in% 22:26,
                   `CDR3.offset.true` == "TRUE") |>
            mutate(`is.Subset16` = ifelse(
                   `CDR3.offset.true` == 24 &
                   `CDR3.offset.Subset16` == "FYCS",
                    TRUE,
                    FALSE))
        
    
    #Subset_16
    Subset_16 <- df1 |>
        filter(`is.Subset16` == "TRUE") |>
        select(1:(ncol(df1) - 7))
    
    
    #Satellite
    Subset_16_Sat <- df1 |>
        filter(`is.Subset16` == "FALSE") |>
        select(1:(ncol(df1) - 7))
    
    #save file
    df_list <- list(Subset_16 = Subset_16, Satellite = Subset_16_Sat)
    
    # Specify the file path
    excel_file <- paste0(file, "/Subset_16.xlsx")
    
    # Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    
    
    # --------------- subset 77 ------------------------------------
    
    #filtering data frame
    #the GW motif of subset 77 is at positions 8 to 9, with an offset 
    #relaxation of +-2 the relevant positions are 6-11
    df1 <- df |>
            filter(`V-DOMAIN Functionality` == "productive" |
                   `V-DOMAIN Functionality` == "productive (see comment)",
                   `CDR3-IMGT length` != "X") |>
            mutate(`IGHV.subgroup` = substr(`V-GENE and allele`, 8, 12),
                   `IGHV.gene` = substr(`V-GENE and allele`, 8, 15),
                   `IGHJ.gene` = substr(`J-GENE and allele`, 8, 12),
                   `CDR3.offset` = substr(`AA JUNCTION`, 7, 12),
                   `CDR3.offsetGW.true` = grepl("GW", `CDR3.offset`),
                   `CDR3.offset.Subset77` = substr(`AA JUNCTION`, 9, 10)) |>
            filter(`IGHV.subgroup` %in% c("IGHV2", "IGHV4", "IGHV6"),
                   `IGHJ.gene` == "IGHJ4", 
                   `CDR3-IMGT length` %in%  12:16,
                   `CDR3.offsetGW.true` == "TRUE") |>
            mutate(`is.Subset77` = ifelse(
                   `CDR3-IMGT length` == 14 &
                   `CDR3.offset.Subset77` == "GW",
                    TRUE,
                    FALSE))
        
    
    #Subset_77
    Subset_77 <- df1 |>
        filter(`is.Subset77` == "TRUE") |>
        select(1:(ncol(df1) - 7))
    
    
    #Satellite
    Subset_77_Sat <- df1 |>
        filter(`is.Subset77` == "FALSE") |>
        select(1:(ncol(df1) - 7))
    
    
    #save file
    df_list <- list(Subset_77 = Subset_77, Satellite = Subset_77_Sat)
    
    # Specify the file path
    excel_file <- paste0(file, "/Subset_77.xlsx")
    
    # Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    
    # --------------- subset 3C2,3C3 ------------------------------------
    #filtering data frame
    df1 <- df |>
            filter(`V-DOMAIN Functionality` == "productive" |
                   `V-DOMAIN Functionality` == "productive (see comment)",
                   `CDR3-IMGT length` != "X") |>
            mutate(`IGHV.subgroup` = substr(`V-GENE and allele`, 8, 12),
                   `IGHV.gene` = substr(`V-GENE and allele`, 8, 15),
                   `IGHJ.gene` = substr(`J-GENE and allele`, 8, 12),
                   `CDR3.offset` = substr(`AA JUNCTION`, 4, 16),
                   `CDR3.offset.true` = grepl("DIVVVPAA", `CDR3.offset`),
                   `CDR3.offset.Subset3C2` = substr(`AA JUNCTION`, 6, 13),
                   `CDR3.offset.Subset3C3` = substr(`AA JUNCTION`, 7, 14)) |>
            filter(`IGHV.subgroup` %in% c("IGHV3"),
                   `IGHJ.gene` == "IGHJ6", 
                   `CDR3-IMGT length` %in% 20:24,
                   `CDR3.offset.true` == "TRUE") |>
            mutate(`is.Subset3C2` = ifelse(
                   `CDR3-IMGT length` == 22 &
                   `CDR3.offset.Subset3C2` == "DIVVVPAA",
                    TRUE,
                    FALSE)) |>
            mutate(`is.Subset3C3` = ifelse(
                   `CDR3-IMGT length` == 9 &
                   `CDR3.offset.Subset3C3` == "DIVVVPAA",
                    TRUE,
                    FALSE))
    
    
    
    
    #Subset3C2
    Subset_3C2 <- df1 |>
        filter(`is.Subset3C2` == "TRUE") |>
        select(1:(ncol(df1) - 9))
    
    #Subset3C3
    Subset_3C3 <- df1 |>
        filter(`is.Subset3C3` == "TRUE") |>
        select(1:(ncol(df1) - 9))  
    
    
    
    #Satellite
    
    Subset_3C2_3C3_sat <- df1 |>
        filter(`is.Subset3C2` == "FALSE",
               `is.Subset3C3` == "FALSE") |>
        select(1:(ncol(df1) - 9))
    
    
    #save file
    df_list <- list(Subset_3C2 = Subset_3C2, Subset_3C3 = Subset_3C3, 
                    Satellite = Subset_3C2_3C3_sat)
    

    # Specify the file path
    excel_file <- paste0(file, "/Subset_3C2_3C3.xlsx")
    
    # Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    
    
    
    
    # --------------- subset 4 ---------------------
    
    #filtering data frame
    df1 <- df |>
            filter(`V-DOMAIN Functionality` == "productive" |
                   `V-DOMAIN Functionality`== "productive (see comment)",
                   `CDR3-IMGT length` != "X") |>
            mutate(`IGHV.subgroup` = substr(`V-GENE and allele`, 8, 15),
                   `IGHV.gene` = substr(`V-GENE and allele`, 8, 15),
                   `IGHJ.gene` = substr(`J-GENE and allele`, 8, 12),
                   `CDR3.offset` = substr(`AA JUNCTION`, 10, 19),
                   `CDR3.offsetKR.true` = grepl("KRYYYY", `CDR3.offset`),
                   `CDR3.offsetRR.true` = grepl("RRYYYY", `CDR3.offset`),
                   `CDR3.offset.Subset4` = substr(`AA JUNCTION`, 12, 17)) |>
            filter(`IGHV.subgroup` %in% c("IGHV2", "IGHV4", "IGHV6"),
                   `IGHJ.gene` == "IGHJ6", 
                   `CDR3-IMGT length` %in% 18:22,
                   `CDR3.offsetKR.true` == "TRUE" | 
                   `CDR3.offsetRR.true` == "TRUE") |>
            mutate(`is.Subset4` = ifelse(
                   `CDR3-IMGT length` == 20 &
                   `CDR3.offset.Subset4` == "KRYYYY" | 
                   `CDR3.offset.Subset4` == "RRYYYY",
                    TRUE,
                    FALSE))
        
        
    #Subset_4
    Subset_4 <- df1 |>
        filter(`is.Subset4` == "TRUE") |>
        select(1:(ncol(df1) - 8))
    
    
    
    #Satellite
    Subset_4_Sat <- df1 |>
        filter(`is.Subset4` == "FALSE") |>
        select(1:(ncol(df1) - 8))
    
    
    
    #save file
    df_list <- list(Subset_4 = Subset_4, Satellite = Subset_4_Sat)

    
    # Specify the file path
    excel_file <- paste0(file, "/Subset_4.xlsx")
    
    # Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    
    
    # --------------- subset 5 ---------------------
    
    #filtering data frame
    df1 <- df |>
        filter(`V-DOMAIN Functionality` == "productive" |
               `V-DOMAIN Functionality`== "productive (see comment)",
               `CDR3-IMGT length` != "X") |>
        mutate(`IGHV.subgroup` = substr(`V-GENE and allele`, 8, 15),
               `IGHV.gene` = substr(`V-GENE and allele`, 8, 15),
               `IGHJ.gene` = substr(`J-GENE and allele`, 8, 12),
               `CDR3.offset` = substr(`AA JUNCTION`, 6, 12),
               `CDR3.offset.GVI.true` = grepl("GVI", `CDR3.offset`),
               `CDR3.offset.GVV.true` = grepl("GVV", `CDR3.offset`),
               `CDR3.offset.Subset5` = substr(`AA JUNCTION`, 8, 10)) |>
        filter(`IGHV.subgroup` %in% c("IGHV2", "IGHV4", "IGHV6"),
               `IGHJ.gene` == "IGHJ6", 
               `CDR3-IMGT length` %in% 18:22,
               `CDR3.offset.GVI.true` == "TRUE" |
               `CDR3.offset.GVV.true` == "TRUE") |>
        mutate(`is.Subset5` = ifelse(
               `CDR3-IMGT length` == 20 &
               `CDR3.offset.Subset5` == "GVV" | 
               `CDR3.offset.Subset5` == "GVI",
                TRUE,
                FALSE))
    
    
    
    #Subset5
    Subset_5 <- df1 |>
        filter(`is.Subset5` == "TRUE") |>
        select(1:(ncol(df1) - 8))
    
    
    
    #Satellite
    
    Subset_5_Sat <- df1 |>
        filter(`is.Subset5` == "FALSE") |>
        select(1:(ncol(df1) - 8))
    
    
    
    #save file
    df_list <- list(Subset_5 = Subset_5, Satellite = Subset_5_Sat)
    
    # Specify the file path
    excel_file <- paste0(file, "/Subset_5.xlsx")
    
    # Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    
    
    # --------------- subset 6 ---------------------
    #the WGSY[K/R] motif of subset 6 is at positions 10-14, with an offset 
    #relaxation of +-2 the relevant positions are 8-16
    
    df1 <- df |>
            filter(`V-DOMAIN Functionality` == "productive" |
                   `V-DOMAIN Functionality`== "productive (see comment)",
                   `CDR3-IMGT length` != "X") |>
            mutate(`IGHV.subgroup` = substr(`V-GENE and allele`, 8, 12),
                   `IGHV.gene` = substr(`V-GENE and allele`, 8, 15),
                   `IGHJ.gene` = substr(`J-GENE and allele`, 8, 12),
                   `CDR3.offset` = substr(`AA JUNCTION`, 9, 17),
                   `CDR3.offset.WGSYK.true` = grepl("WGSYK", `CDR3.offset`),
                   `CDR3.offset.WGSYR.true` = grepl("WGSYR", `CDR3.offset`),
                   `CDR3.offset.Subset6` = substr(`AA JUNCTION`, 11, 15)) |>
            filter(`IGHV.subgroup` %in% c("IGHV1", "IGHV5", "IGHV7"),
                   `IGHJ.gene` == "IGHJ3", 
                   `CDR3-IMGT length` %in% 19:23,
                   `CDR3.offset.WGSYK.true` == "TRUE" |
                   `CDR3.offset.WGSYR.true` == "TRUE")|>
            mutate(`is.Subset6` = ifelse(
                   `CDR3-IMGT length` == 21 &
                   `CDR3.offset.Subset6` == "WGSYK" | 
                   `CDR3.offset.Subset6` == "WGSYR",
                    TRUE,
                    FALSE))
    
    
    
    #Subset_6
    Subset_6 <- df1 |>
        filter(`is.Subset6` == "TRUE") |>
        select(1:(ncol(df1) - 9))
    
    
    
    
    #Satellite
    Subset_6_Sat <- df1 |>
        filter(`is.Subset6` == "FALSE") |>
        select(1:(ncol(df1) - 9))
    
    
    #save file
    df_list <- list(Subset_6 = Subset_6, Satellite = Subset_6_Sat)
    
    # Specify the file path
    excel_file <- paste0(file, "/Subset_6.xlsx")
    
    # Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    
    
    
    # --------------- subset 7 -----------------------------------
    #filtering data frame
    #the DFWSGY motif of subsets 7C2 and 7D3 is at positions 7-13, with an offset 
    #relaxation of +-2 the relevant positions are 5-15
    #creates a new column containing the amino acids from CDR3 position 5 to 15 
    df1 <- df |>
        filter(`V-DOMAIN Functionality` == "productive" |
               `V-DOMAIN Functionality` == "productive (see comment)",
               `CDR3-IMGT length` != "X") |>
        mutate(`IGHV.subgroup` = substr(`V-GENE and allele`, 8, 12),
               `IGHV.gene` = substr(`V-GENE and allele`, 8, 15),
               `IGHJ.gene` = substr(`J-GENE and allele`, 8, 12),
               `CDR3.offset` = substr(`AA JUNCTION`, 6, 16),
               `CDR3.offset.true` = grepl("DFWSGY", `CDR3.offset`),
               `CDR3.offset.Subset7C2` = substr(`AA JUNCTION`, 8, 13),
               `CDR3.offset.Subset7D3` = substr(`AA JUNCTION`, 7, 14)) |>
        filter(`IGHV.subgroup` %in% c("IGHV1", "IGHV5", "IGHV7"),
               `IGHJ.gene` == "IGHJ6", 
               `CDR3-IMGT length` %in% 21:26,
               `CDR3.offset.true` == "TRUE") |>
        mutate(`is.Subset.7C2` = ifelse(
               `CDR3-IMGT length` == 23 &
               `CDR3.offset.Subset7C2` == "DFWSGY",
                TRUE,
                FALSE)) |>
        mutate(`is.Subset.7D3` = ifelse(
               `CDR3-IMGT length` == 24 &
               `CDR3.offset.Subset7D3` == "DFWSGY",
                TRUE,
                FALSE))
    
    
    
    
    
    #Subset_7C2
    Subset_7C2 <- df1 |>
        filter(`is.Subset.7C2` == "TRUE") |>
        select(1:(ncol(df1) - 9))
    
    
    
    
    #Subset_7D3
    Subset_7D3 <- df1 |>
        filter(`is.Subset.7D3` == "TRUE") |>
        select(1:(ncol(df1)-9))  
    
    
    
    
    #Satellite
    
    Subset_7C2_7D3_sat <- df1 |>
        filter(`is.Subset.7C2` == "FALSE",
               `is.Subset.7D3` == "FALSE") |>
        select(1:(ncol(df1) - 9))
    
    
    
    #save file
    df_list <- list(Subset_7C2 = Subset_7C2, Subset_7D3 = Subset_7D3,
                    Satellite = Subset_7C2_7D3_sat)
    

    # Specify the file path
    excel_file <- paste0(file, "/Subset_7C2_7D3.xlsx")
    
    # Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    
    
    
    
    # --------------- subset 10 -------------
    
    #filtering data frame
    #the GYCSSTSC motif of subset 10 is at positions 6-13, with an offset 
    #relaxation of +-2 the relevant positions are 4-15
    
    df1 <- df |>
        filter(`V-DOMAIN Functionality` == "productive" |
               `V-DOMAIN Functionality` == "productive (see comment)",
               `CDR3-IMGT length` != "X") |>
        mutate(`IGHV.subgroup` = substr(`V-GENE and allele`, 8, 12),
               `IGHV.gene` = substr(`V-GENE and allele`, 8, 15),
               `IGHJ.gene` = substr(`J-GENE and allele`, 8, 12),
               `CDR3.offset` = substr(`AA JUNCTION`, 5, 16),
               `CDR3.offset.true` = grepl("GYCSSTSC", `CDR3.offset`),
               `CDR3.offset.Subset10` = substr(`AA JUNCTION`, 7, 14)) |>
        filter(`IGHV.subgroup` %in% c("IGHV2", "IGHV4", "IGHV6"),
               `IGHJ.gene` == "IGHJ6", 
               `CDR3-IMGT length` %in% 20:24,
               `CDR3.offset.true` == "TRUE") |>
        mutate(`is.Subset10` = ifelse(
               `CDR3-IMGT length` == 22 &
               `CDR3.offset.Subset10` == "GYCSSTSC",
                TRUE,
                FALSE))
    
    
    
    #Subset_10
    Subset_10 <- df1 |>
        filter(`is.Subset10` == "TRUE") |>
        select(1:(ncol(df1) - 7))
    
    
    
    #Satellite
    Subset_10_Sat <- df1 |>
        filter(`is.Subset10` == "FALSE") |>
        select(1:(ncol(df1) - 7))
    
    
    #save file
    df_list <- list(Subset_10 = Subset_10, Satellite = Subset_10_Sat)
    
    # Specify the file path
    excel_file <- paste0(file, "/Subset_10.xlsx")
    
    # Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    
    
    # --------------- subset 12 ---------------------
    
    #filtering data frame
    #the YYDSSGYY motif of subset 12 is at positions 6-13, with an offset 
    #relaxation of +-2 the relevant positions are 4-15
    df1 <- df |>
            filter(`V-DOMAIN Functionality` == "productive" |
                   `V-DOMAIN Functionality` == "productive (see comment)",
                   `CDR3-IMGT length` != "X") |>
            mutate(`IGHV.subgroup` = substr(`V-GENE and allele`, 8, 12),
                   `IGHV.gene` = substr(`V-GENE and allele`, 8, 15),
                   `IGHJ.gene` = substr(`J-GENE and allele`, 8, 12),
                   `CDR3.offset` = substr(`AA JUNCTION`, 5, 16),
                   `CDR3.offset.true` = grepl("YYDSSGYY", `CDR3.offset`),
                   `CDR3.offset.Subset12` = substr(`AA JUNCTION`, 7, 14)) |>
            filter(`IGHV.subgroup` %in% c("IGHV1", "IGHV5", "IGHV7"),
                   `IGHJ.gene` == "IGHJ4", 
                   `CDR3-IMGT length` %in% 17:21,
                   `CDR3.offset.true` == "TRUE") |>
            mutate(`is.Subset.12` = ifelse(
                   `CDR3-IMGT length` == 19 &
                   `CDR3.offset.Subset12` == "YYDSSGYY",
                    TRUE,
                    FALSE))
        
        
    
    
    #Subset_12
    Subset_12 <- df1 |>
        filter(`is.Subset.12` == "TRUE") |>
        select(1:(ncol(df1) - 7))
    
    
    #Satellite
    Subset_12_Sat <- df1 |>
        filter(`is.Subset.12` == "FALSE") |>
        select(1:(ncol(df1) - 7))
    
    
    
    
    #save file
    df_list <- list(Subset_12 = Subset_12, Satellite = Subset_12_Sat)
    
    #Specify the file path
    excel_file <- paste0(file, "/Subset_12.xlsx")
    
    #Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    
    
    # --------------- subset 14 ---------------------
    
    #filtering data frame
    #the RGG motif of subset 14 is at positions 2-4, with an offset 
    #relaxation of +-2 the relevant positions are 1-6
    df1 <- df |>
            filter(`V-DOMAIN Functionality` == "productive" |
                   `V-DOMAIN Functionality` == "productive (see comment)",
                   `CDR3-IMGT length` != "X") |>
            mutate(`IGHV.subgroup` = substr(`V-GENE and allele`, 8, 12),
                   `IGHV.gene` = substr(`V-GENE and allele`, 8, 15),
                   `IGHJ.gene` = substr(`J-GENE and allele`, 8, 12),
                   `CDR3.offset` = substr(`AA JUNCTION`, 2, 7),
                   `CDR3.offset.true` = grepl("RGG", `CDR3.offset`),
                   `CDR3.offset.Subset14` = substr(`AA JUNCTION`, 3, 5)) |>
            filter(`IGHV.subgroup` %in% c("IGHV2", "IGHV4", "IGHV6"),
                   `IGHJ.gene` == "IGHJ4", 
                   `CDR3-IMGT length` %in% 8:12,
                   `CDR3.offset.true` == "TRUE") |>
            mutate(`is.Subset.14` = ifelse(
                   `CDR3-IMGT length` == 10 &
                   `CDR3.offset.Subset14` == "RGG",
                    TRUE,
                    FALSE))
    
    
    
    
    #Subset_14
    Subset_14 <- df1 |>
        filter(`is.Subset.14` == "TRUE") |>
        select(1:(ncol(df1) - 7))
    
    
    
    
    #Satellite
    Subset_14_Sat <- df1 |>
        filter(`is.Subset.14` == "FALSE") |>
        select(1:(ncol(df1) - 7))
    
    
    
    
    #save file
    df_list <- list(Subset_14 = Subset_14, Satellite = Subset_14_Sat)
    
    #Specify the file path
    excel_file <- paste0(file, "/Subset_14.xlsx")
    
    #Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    
    
    # --------------- subset 28 ---------------------
    
    #filtering data frame
    #the SGS motif of subset 28A is at positions 5-7, with an offset 
    #relaxation of +-2 the relevant positions are 3-9 
    #creates a new column containing the amino acids from CDR3 position 4 to 15 
    df1 <- df |>
        filter(`V-DOMAIN Functionality` == "productive" |
               `V-DOMAIN Functionality` == "productive (see comment)",
               `CDR3-IMGT length` != "X") |>
        mutate(`IGHV.subgroup` = substr(`V-GENE and allele`, 8, 12),
               `IGHV.gene` = substr(`V-GENE and allele`, 8, 15),
               `IGHJ.gene` = substr(`J-GENE and allele`, 8, 12),
               `CDR3.offset` = substr(`AA JUNCTION`, 4, 10),
               `CDR3.offset.true` = grepl("SGS", `CDR3.offset`),
               `CDR3.offset.Subset28A` = substr(`AA JUNCTION`, 6, 8)) |>
        filter(`IGHV.subgroup` %in% c("IGHV1", "IGHV5", "IGHV7"),
               `IGHJ.gene` == "IGHJ6", 
               `CDR3-IMGT length` %in% 15:19,
               `CDR3.offset.true` == "TRUE") |>
        mutate(`is.Subset28A` = ifelse(
               `CDR3-IMGT length` == 17 &
               `CDR3.offset.Subset28A` == "SGS",
                TRUE,
                FALSE))
    
    
    
    #Subset_28A
    Subset_28A <- df1 |>
        filter(`is.Subset28A` == "TRUE") |>
        select(1:(ncol(df1) - 7))
    
    
    
    
    #Satellite
    Subset_28A_Sat <- df1 |>
        filter(`is.Subset28A` == "FALSE") |>
        select(1:(ncol(df1) - 7))
    
    
    #save file
    df_list <- list(Subset_28A = Subset_28A, Satellite = Subset_28A_Sat)
    
    #Specify the file path
    excel_file <- paste0(file, "/Subset_28A.xlsx")
    
    #Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    
    
    
    # --------------- subset 31 ---------------------
    
    #filtering data frame
    #the FWSGY motif of subset 31 is at positions 6-10, with an offset 
    #relaxation of +-2 the relevant positions are 4-12
    df1 <- df |>
        filter(`V-DOMAIN Functionality` == "productive" |
               `V-DOMAIN Functionality` == "productive (see comment)",
               `CDR3-IMGT length` != "X") |>
        mutate(`IGHV.subgroup` = substr(`V-GENE and allele`, 8, 12),
               `IGHV.gene` = substr(`V-GENE and allele`, 8, 15),
               `IGHJ.gene` = substr(`J-GENE and allele`, 8, 12),
               `CDR3.offset` = substr(`AA JUNCTION`, 5, 13),
               `CDR3.offset.true` = grepl("FWSGY", `CDR3.offset`),
               `CDR3.offset.Subset31` = substr(`AA JUNCTION`, 7, 11)) |>
        filter(`IGHV.subgroup` %in% c("IGHV3"),
               `IGHJ.gene` == "IGHJ6", 
               `CDR3-IMGT length` %in% 19:23,
               `CDR3.offset.true` == "TRUE") |>
        mutate(`is.Subset.31` = ifelse(
               `CDR3-IMGT length` == 21 &
               `CDR3.offset.Subset31` == "FWSGY",
                TRUE,
                FALSE))
    
    
    
    
    #Subset_31
    Subset_31 <- df1 |>
        filter(`is.Subset.31` == "TRUE") |>
        select(1:(ncol(df1) - 7))
    
    
    
    #Satellite
    Subset_31_Sat <- df1 |>
        filter(`is.Subset.31` == "FALSE") |>
        select(1:(ncol(df1) - 7))
    
    
    #save file
    df_list <- list(Subset_31 = Subset_31, Satellite = Subset_31_Sat)
    
    #Specify the file path
    excel_file <- paste0(file, "/Subset_31.xlsx")
    
    #Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    
    
    
    
    
    # --------------- subset 59 ---------------------
    
    #filtering data frame
    #the FWSG motif of subset 59 is at positions 6-9, with an offset relaxation 
    #of +-2 the relevant positions are 4-11
    df1 <- df |>
        filter(`V-DOMAIN Functionality` == "productive" |
               `V-DOMAIN Functionality` == "productive (see comment)",
               `CDR3-IMGT length` != "X") |>
        mutate(`IGHV.subgroup` = substr(`V-GENE and allele`, 8, 12),
               `IGHV.gene` = substr(`V-GENE and allele`, 8, 15),
               `IGHJ.gene` = substr(`J-GENE and allele`, 8, 12),
               `CDR3.offset` = substr(`AA JUNCTION`, 5, 12),
               `CDR3.offset.true` = grepl("FWSG", `CDR3.offset`),
               `CDR3.offset.Subset59` = substr(`AA JUNCTION`, 7, 10)) |>
        filter(`IGHV.subgroup` %in% c("IGHV1", "IGHV5", "IGHV7"),
               `IGHJ.gene` == "IGHJ6", 
               `CDR3-IMGT length` %in% 10:14,
               `CDR3.offset.true` == "TRUE") |>
        mutate(`is.Subset.59` = ifelse(
               `CDR3-IMGT length` == 12 &
               `CDR3.offset.Subset59` == "FWSGY",
                TRUE,
                FALSE))
    
    
    
    
    #Subset_59
    Subset_59 <- df1 |>
        filter(`is.Subset.59` == "TRUE") |>
        select(1:(ncol(df1) - 7))
    
    
    
    #Satellite
    Subset_59_Sat <- df1 |>
        filter(`is.Subset.59` == "FALSE") |>
        select(1:(ncol(df1) - 7))
    
    
    
    
    #save file
    df_list <- list(Subset_59 = Subset_59, Satellite = Subset_59_Sat)
    
    #Specify the file path
    excel_file <- paste0(file, "/Subset_59.xlsx")
    
    #Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    
    
    
    
    # --------------- subset 64B ---------------------
    
    #filtering data frame
    
    #the LVV motif of subset 64B is at positions 6-8, with an offset 
    #relaxation of +-2 the relevant positions are 4-10
    df1 <- df |>
        filter(`V-DOMAIN Functionality` == "productive" |
               `V-DOMAIN Functionality` == "productive (see comment)",
               `CDR3-IMGT length` != "X") |>
        mutate(`IGHV.subgroup` = substr(`V-GENE and allele`, 8, 12),
               `IGHV.gene` = substr(`V-GENE and allele`, 8, 15),
               `IGHJ.gene` = substr(`J-GENE and allele`, 8, 12),
               `CDR3.offset` = substr(`AA JUNCTION`, 5, 11),
               `CDR3.offset.true` = grepl("LVV", `CDR3.offset`),
               `CDR3.offset.Subset64B` = substr(`AA JUNCTION`, 7, 9)) |>
        filter(`IGHV.subgroup` %in% c("IGHV3"),
               `IGHJ.gene` == "IGHJ6", 
               `CDR3-IMGT length` %in% 19:23,
               `CDR3.offset.true` == "TRUE") |>
        mutate(`is.Subset.64B` = ifelse(
               `CDR3-IMGT length` == 21 &
               `CDR3.offset.Subset64B` == "LVV",
                TRUE,
                FALSE))
    
    
    
    #Subset_64B
    Subset_64B <- df1 |>
        filter(`is.Subset.64B` == "TRUE") |>
        select(1:(ncol(df1) - 7))
    
    
    
    #Satellite
    Subset_64B_Sat <- df1 |>
        filter(`is.Subset.64B` == "FALSE") |>
        select(1:(ncol(df1) - 7))
    
    
    
    
    #save file
    df_list <- list(Subset_64B = Subset_64B, Satellite = Subset_64B_Sat)
    
    #Specify the file path
    excel_file <- paste0(file, "/Subset_64B.xlsx")
    
    #Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    
    
    
    # --------------- subset 73 ---------------------
    #filtering data frame
    
    #the [K/R]D motif of subset 73 is at positions 2 to 3, with an offset 
    #relaxation of +-2 the relevant positions are 1-5
    
    df1 <- df |>
        filter(`V-DOMAIN Functionality` == "productive" |
               `V-DOMAIN Functionality` == "productive (see comment)",
               `CDR3-IMGT length` != "X")  |>
        mutate(`IGHV.subgroup` = substr(`V-GENE and allele`, 8, 12),
               #IGHV.gene = substr(V.GENE.and.allele, 8, 15),
               `IGHJ.gene` = substr(`J-GENE and allele`, 8, 12),
               `CDR3.offset` = substr(`AA JUNCTION`, 2, 6),
               `CDR3.offsetKD.true` = grepl("KD", `CDR3.offset`),
               `CDR3.offsetRD.true` = grepl("RD", `CDR3.offset`),
               `CDR3.offset2` = substr(`AA JUNCTION`, 10, 14),
               `CDR3.offsetDY.true` = grepl("DY", `CDR3.offset2`),
               `CDR3.offset.Subset73` = substr(`AA JUNCTION`, 3, 4)) |>
        filter(`IGHV.subgroup` %in% c("IGHV3"),
               `IGHJ.gene` == "IGHJ4", 
               `CDR3-IMGT length` %in% 10:14,
               `CDR3.offsetKD.true` == "TRUE" |
               `CDR3.offsetRD.true` == "TRUE",
               `CDR3.offsetDY.true` == "TRUE")|>
        mutate(`is.Subset73` = ifelse(
               `CDR3-IMGT length` == 12 &
               `CDR3.offset.Subset73` == "KD" | 
               `CDR3.offset.Subset73` == "RD",
               TRUE,
               FALSE))
    
    #Subset_73
    Subset_73 <- df1 |>
        filter(`is.Subset73` == "TRUE") |>
        select(1:(ncol(df1) - 9))
    
    
    
    #Satellite
    Subset_73_Sat <- df1 |>
        filter(`is.Subset73` == "FALSE") |>
        select(1:(ncol(df1) - 9))
    
    
    #save file
    df_list <- list(Subset_73 = Subset_73, Satellite = Subset_73_Sat)
    
    #Specify the file path
    excel_file <- paste0(file, "/Subset_73.xlsx")
    
    #Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    
    
    
    
    
    # --------------- subset 77 ---------------------
    
    #filtering data frame
    #the GW motif of subset 77 is at positions 8 to 9, with an offset 
    #relaxation of +-2 the relevant positions are 6-11
    df1 <- df |>
        filter(`V-DOMAIN Functionality` == "productive" |
               `V-DOMAIN Functionality` == "productive (see comment)",
               `CDR3-IMGT length` != "X") |>
        mutate(`IGHV.subgroup` = substr(`V-GENE and allele`, 8, 12),
               `IGHV.gene` = substr(`V-GENE and allele`, 8, 15),
               `IGHJ.gene` = substr(`J-GENE and allele`, 8, 12),
               `CDR3.offset` = substr(`AA JUNCTION`, 7, 12),
               `CDR3.offsetGW.true` = grepl("GW", `CDR3.offset`),
               `CDR3.offset.Subset77` = substr(`AA JUNCTION`, 9, 10)) |>
        filter(`IGHV.subgroup` %in% c("IGHV2", "IGHV4", "IGHV6"),
               `IGHJ.gene` == "IGHJ4", 
               `CDR3-IMGT length` %in%  12:16,
               `CDR3.offsetGW.true` == "TRUE") |>
        mutate(`is.Subset77` = ifelse(
               `CDR3-IMGT length` == 14 &
               `CDR3.offset.Subset77` == "GW",
                TRUE,
                FALSE))
    
    
    #Subset_77
    Subset_77 <- df1 |>
        filter(`is.Subset77` == "TRUE") |>
        select(1:(ncol(df1) - 7))
    
    
    #Satellite
    Subset_77_Sat <- df1 |>
        filter(`is.Subset77` == "FALSE") |>
        select(1:(ncol(df1) - 7))
    
    
    #save file
    df_list <- list(Subset_77 = Subset_77, Satellite = Subset_77_Sat)
    
    #Specify the file path
    excel_file <- paste0(file, "/Subset_77.xlsx")
    
    #Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    
    
    
    
    # --------------- subset 111 ---------------------
    
    #filtering data frame
    
    #the SGG motif of subset 111 is at positions 4 to 6, with an offset
    #relaxation of +-2 the relevant positions are 2-8
    df1 <- df |>
        filter(`V-DOMAIN Functionality` == "productive" |
               `V-DOMAIN Functionality` == "productive (see comment)",
               `CDR3-IMGT length` != "X") |>
        mutate(`IGHV.subgroup` = substr(`V-GENE and allele`, 8, 12),
               `IGHV.gene` = substr(`V-GENE and allele`, 8, 15),
               `IGHJ.gene` = substr(`J-GENE and allele`, 8, 12),
               `CDR3.offset` = substr(`AA JUNCTION`, 3, 9),
               `CDR3.offsetSGG.true` = grepl("SGG", `CDR3.offset`),
               `CDR3.offset.Subset111` = substr(`AA JUNCTION`, 5, 7)) |>
        filter(`IGHV.subgroup` %in% c("IGHV1", "IGHV5", "IGHV7"),
               `IGHJ.gene` == "IGHJ6", 
               `CDR3-IMGT length` %in% 8:12,
               `CDR3.offsetSGG.true` == "TRUE") |>
        mutate(`is.Subset111` = ifelse(
               `CDR3-IMGT length` == 10 &
               `CDR3.offset.Subset111` == "SGG",
                TRUE,
                FALSE))
    
    
    
    
    #Subset_111
    Subset_111 <- df1 |>
        filter(`is.Subset111` == "TRUE") |>
        select(1:(ncol(df1) - 7))
    
    
    
    #Satellite
    Subset_111_Sat <- df1 |>
        filter(`is.Subset111` == "FALSE") |>
        select(1:(ncol(df1) - 7))
    
    
    #save file
    df_list <- list(Subset_111 = Subset_111, Satellite = Subset_111_Sat)
    
    #Specify the file path
    excel_file <- paste0(file, "/Subset_111.xlsx")
    
    #Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    
    
    
    
    
    
    
    # --------------- subset 148Β ---------------------
    
    #filtering data frame
    #the HR motif of subset 148B is at positions 2-3, with an offset 
    #relaxation of +-2 the relevant positions are 1-5
    df1 <- df |>
        filter(`V-DOMAIN Functionality` == "productive" |
               `V-DOMAIN Functionality` == "productive (see comment)",
               `CDR3-IMGT length` != "X") |>
        mutate(`IGHV.subgroup` = substr(`V-GENE and allele`, 8, 12),
               `IGHV.gene` = substr(`V-GENE and allele`, 8, 15),
               `IGHJ.gene` = substr(`J-GENE and allele`, 8, 12),
               `CDR3.offset` = substr(`AA JUNCTION`, 2, 6),
               `CDR3.offset.true` = grepl("HR", `CDR3.offset`),
               `CDR3.offset2` = substr(`AA JUNCTION`, 9, 13),
               `CDR3.offset2.true` = grepl("W", `CDR3.offset2`),
               `CDR3.offset.Subset148B.A` = substr(`AA JUNCTION`, 3, 4),
               `CDR3.offset.Subset148B.B` = substr(`AA JUNCTION`, 11, 11)) |>
        filter(`IGHV.subgroup` %in% c("IGHV2", "IGHV4", "IGHV6"),
               #IGHJ.gene == "IGHJ4", #commend out from Andreas
               `CDR3-IMGT length` %in% 15:19,
               `CDR3.offset.true` == "TRUE",
               `CDR3.offset2.true` == "TRUE") |>
        mutate(`is.Subset.148B` = ifelse(
               `CDR3-IMGT length` == 17 &
               `CDR3.offset.Subset148B.A` == "HR" &
               `CDR3.offset.Subset148B.B` == "W",
                TRUE,
                FALSE))
    
    
    
    #Subset_148B
    Subset_148B <- df1 |>
        filter(`is.Subset.148B` == "TRUE") |>
        select(1:(ncol(df1) - 10))
    
    
    
    #Satellite
    Subset_148B_Sat <- df1 |>
        filter(`is.Subset.148B` == "FALSE") |>
        select(1:(ncol(df1) - 10))
    
    
    #save file
    df_list <- list(Subset_148B = Subset_148B, Satellite = Subset_148B_Sat)
    
    #Specify the file path
    excel_file <- paste0(file, "/Subset_148B.xlsx")
    
    #Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    
    
    
    
    
    # --------------- subset 188 ---------------------
    
    #filtering data frame
    #the [Y/F]C motif of subset 188 is at positions 5 to 6, with an offset
    #relaxation of +-2 the relevant positions are 3-8
    df1 <- df |>
        filter(`V-DOMAIN Functionality` == "productive" |
               `V-DOMAIN Functionality` == "productive (see comment)",
               `CDR3-IMGT length` != "X") |>
        mutate(`IGHV.subgroup` = substr(`V-GENE and allele`, 8, 12),
               `IGHV.gene` = substr(`V-GENE and allele`, 8, 15),
               `IGHJ.gene` = substr(`J-GENE and allele`, 8, 12),
               `CDR3.offset` = substr(`AA JUNCTION`, 4, 9),
               `CDR3.offsetYC.true` = grepl("YC", `CDR3.offset`),
               `CDR3.offsetFC.true` = grepl("FC", `CDR3.offset`),
               `CDR3.offset2` = substr(`AA JUNCTION`, 10, 15),
               `CDR3.offsetCR.true` = grepl("CR", `CDR3.offset2`),
               `CDR3.offsetCY.true` = grepl("CY", `CDR3.offset2`),
               `CDR3.offset.Subset188` = substr(`AA JUNCTION`, 6, 7)) |>
        filter(`IGHV.subgroup` %in% c("IGHV3"),
               `IGHJ.gene` == "IGHJ6", 
               `CDR3-IMGT length` %in% 15:19,
               `CDR3.offsetYC.true` == "TRUE" |
               `CDR3.offsetFC.true` == "TRUE",
               `CDR3.offsetCR.true` == "TRUE" |
               `CDR3.offsetCY.true` == "TRUE") |>
        mutate(`is.Subset188` = ifelse(
               `CDR3-IMGT length` == 17 &
               `CDR3.offset.Subset188` == "YC",
               TRUE,
               FALSE))
    
    
    
    
    #Subset_188
    Subset_188 <- df1 |>
        filter(`is.Subset188` == "TRUE") |>
        select(1:(ncol(df1) - 11))
    
    
    
    
    #Satellite
    Subset_188_Sat <- df1 |>
        filter(`is.Subset188` == "FALSE") |>
        select(1:(ncol(df1) - 11))
    
    
    
    
    #save file
    df_list <- list(Subset_188 = Subset_188, Satellite = Subset_188_Sat)
    
    #Specify the file path
    excel_file <- paste0(file, "/Subset_188.xlsx")
    
    #Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    
    
    
    
    
    
    # --------------- subset 201 ---------------------
    
    #filtering data frame
    #the ARR motif of subset 201 is at positions 1-3, with an offset 
    #relaxation of +-2 the relevant positions are 1-5
    
    df1 <- df |>
        filter(`V-DOMAIN Functionality` == "productive" |
               `V-DOMAIN Functionality` == "productive (see comment)",
               `CDR3-IMGT length` != "X") |>
        mutate(`IGHV.subgroup` = substr(`V-GENE and allele`, 8, 12),
               `IGHV.gene` = substr(`V-GENE and allele`, 8, 15),
               `IGHJ.gene` = substr(`J-GENE and allele`, 8, 12),
               `CDR3.offset` = substr(`AA JUNCTION`, 2, 6),
               `CDR3.offset.true` = grepl("ARR", `CDR3.offset`),
               `CDR3.offset2` = substr(`AA JUNCTION`, 6, 10),
               `CDR3.offset2.true` = grepl("W", `CDR3.offset2`),
               `CDR3.offset.Subset201A` = substr(`AA JUNCTION`, 2, 4),
               `CDR3.offset.Subset201B` = substr(`AA JUNCTION`, 8, 8)) |>
        filter(`IGHV.subgroup` %in% c("IGHV2", "IGHV4", "IGHV6"),
               `IGHJ.gene` == "IGHJ3", 
               `CDR3-IMGT length` %in% 15:19,
               `CDR3.offset.true` == "TRUE",
               `CDR3.offset2.true` == "TRUE") |>
        mutate(`is.Subset201` = ifelse(
               `CDR3-IMGT length` == 17 &
               `CDR3.offset.Subset201A` == "ARR" &
               `CDR3.offset.Subset201B` == "W",
                TRUE,
                FALSE))
    
    
    
    #Subset_201
    Subset_201 <- df1 |>
        filter(`is.Subset201` == "TRUE") |>
        select(1:(ncol(df1) - 10))
    
    
    
    #Satellite
    Subset_201_Sat <- df1 |>
        filter(`is.Subset201` == "FALSE") |>
        select(1:(ncol(df1) - 10))
    
    
    
    
    #save file
    df_list <- list(Subset_201 = Subset_201, Satellite = Subset_201_Sat)
    
    #Specify the file path
    excel_file <- paste0(file, "/Subset_201.xlsx")
    
    #Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    
    
    
    # --------------- subset 202 ---------------------
    
    #filtering data frame
    #the GDY motif of subset 202 is at positions 6-8, with an offset 
    #relaxation of +-2 the relevant positions are 4-10
    
    df1 <- df |>
        filter(`V-DOMAIN Functionality` == "productive" |
               `V-DOMAIN Functionality` == "productive (see comment)",
               `CDR3-IMGT length` != "X") |>
        mutate(`IGHV.subgroup` = substr(`V-GENE and allele`, 8, 12),
               `IGHV.gene` = substr(`V-GENE and allele`, 8, 15),
               `IGHJ.gene` = substr(`J-GENE and allele`, 8, 12),
               `CDR3.offset` = substr(`AA JUNCTION`, 5, 11),
               `CDR3.offset.GDY.true` = grepl("GDY", `CDR3.offset`),
               `CDR3.offset.Subset202` = substr(`AA JUNCTION`, 7, 9)) |>
        filter(`IGHV.subgroup` %in% c("IGHV3"),
               #IGHJ.gene == "IGHJ3", commend out by Andreas
               `CDR3-IMGT length` %in% 12:16,
               `CDR3.offset.GDY.true` == "TRUE") |>
        mutate(`is.Subset.202` = ifelse(
               `CDR3-IMGT length` == 14 &
               `CDR3.offset.Subset202` == "GDY",
               TRUE,
               FALSE))
    
    
    
    #Subset_202
    Subset_202 <- df1 |>
        filter(`is.Subset.202` == "TRUE") |>
        select(1:(ncol(df1) - 7))
    
    
    
    #Satellite
    Subset_202_Sat <- df1 |>
        filter(`is.Subset.202` == "FALSE") |>
        select(1:(ncol(df1) - 7))
    
    
    
    
    #save file
    df_list <- list(Subset_202 = Subset_202, Satellite = Subset_202_Sat)
    
    #Specify the file path
    excel_file <- paste0(file, "/Subset_202.xlsx")
    
    #Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    
    
    
    # --------------- subset 252 ---------------------
    
    #filtering data frame
    #the G landmark of subset 252 is at position 7, with an offset 
    #relaxation of +-2 the relevant positions are 5-9
    
    df1 <- df |>
            filter(`V-DOMAIN Functionality` == "productive" |
                   `V-DOMAIN Functionality` == "productive (see comment)",
                   `CDR3-IMGT length` != "X") |>
            mutate(`IGHV.subgroup` = substr(`V-GENE and allele`, 8, 12),
                   `IGHV.gene` = substr(`V-GENE and allele`, 8, 15),
                   `IGHJ.gene` = substr(`J-GENE and allele`, 8, 12),
                   `CDR3.offset` = substr(`AA JUNCTION`, 6, 10),
                   `CDR3.offset.true` = grepl("G", `CDR3.offset`),
                   `CDR3.offset2` = substr(`AA JUNCTION`, 8, 12),
                   `CDR3.offset2.true` = grepl("F", `CDR3.offset2`),
                   `CDR3.offset.Subset252A` = substr(`AA JUNCTION`, 8, 8),
                   `CDR3.offset.Subset252B` = substr(`AA JUNCTION`, 10, 10)) |>
            filter(`IGHV.subgroup` %in% c("IGHV3"),
                   `IGHJ.gene` == "IGHJ6", 
                   `CDR3-IMGT length` %in% 17:21,
                   `CDR3.offset.true` == "TRUE",
                   `CDR3.offset2.true` == "TRUE") |>
            mutate(`is.Subset.252` = ifelse(
                   `CDR3-IMGT length` == 19 &
                   `CDR3.offset.Subset252A` == "G" &
                   `CDR3.offset.Subset252B` == "F",
                    TRUE,
                    FALSE))
    
    
    
    #Subset_252
    Subset_252 <- df1 |>
        filter(`is.Subset.252` == "TRUE") |>
        select(1:(ncol(df1) - 10))
    
    
    
    #Satellite
    Subset_252_Sat <- df1 |>
        filter(`is.Subset.252` == "FALSE") |>
        select(1:(ncol(df1) - 10))
    
    
    
    
    #save file
    df_list <- list(Subset_252 = Subset_252, Satellite = Subset_252_Sat)
    
    #Specify the file path
    excel_file <- paste0(file, "/Subset_252.xlsx")
    
    #Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    
    
    
    
    
    
    
    
    
    
}
