


rm(list = ls())
gc()


library(data.table)
library(dplyr)




#function-----------------------------------------
carmilla <- function(
        
    q, IGHV.subgroup = FALSE, vgene = FALSE, jgene = FALSE, 
    CDR3.offset.x = NULL, CDR3.offset.y = NULL, 
    CDR3.offset.a = "[A-Z]+", CDR3.offset.b = "[A-Z]+",
    CDR3.offset.Subset.1.x = NULL, CDR3.offset.Subset.1.y = NULL,
    CDR3.offset.Subset.2.x = NULL, CDR3.offset.Subset.2.y = NULL,
    subgroup = c("[A-Z]+[0-9]", "[A-Z]+[0-9]", "[A-Z]+[0-9]"),
    ighj.filter = "[A-Z]+", CDR3.IMGT.length.filter = NULL, 
    CDR3.offset.a.true.filter = FALSE, CDR3.offset.b.true.filter = FALSE
    
) {
    
    
    q1 <- q |>
        filter(
            `V-DOMAIN Functionality` == "productive" |
            `V-DOMAIN Functionality`== "productive (see comment)",
            `CDR3-IMGT length` != "X"
        ) 
    
    
    if(IGHV.subgroup) {
        
        q1 <- q1 |> mutate(`IGHV.subgroup` = substr(`V-GENE and allele`, 8, 12))
        
    }
    
    
    if(vgene) {
        
        q1 <- q1 |> mutate(`IGHV.gene` = substr(`V-GENE and allele`, 8, 15))
        
    }
    
    
    if(jgene) {
        
        q1 <- q1 |> mutate(`IGHJ.gene` = substr(`J-GENE and allele`, 8, 12))
        
    }
    
    
    if(!is.null(CDR3.offset.x)  & !is.null(CDR3.offset.y)) {
        
        q1 <- q1 |> mutate(
            `CDR3.offset` = substr(`AA JUNCTION`, CDR3.offset.x, CDR3.offset.y)
        )
        
    }
    
    
    if(!is.null(CDR3.offset.a)) {
        
        q1 <- q1 |> mutate(
            `CDR3.offset.a.true` = grepl(CDR3.offset.a, `CDR3.offset`)
        )
        
    }
    
    
    if(!is.null(CDR3.offset.b)) {
        
        q1 <- q1 |> mutate(
            `CDR3.offset.b.true` = grepl(CDR3.offset.b, `CDR3.offset`)
        )
        
    }
    
    
    if(CDR3.offset.Subset.1.x & CDR3.offset.Subset.1.y) {
        
        q1 <- q1 |> mutate(
                `CDR3.offset.Subset.1` = substr(
                `AA JUNCTION`, CDR3.offset.Subset.1.x, CDR3.offset.Subset.1.y
                )
             )
        
    }
    
    if(CDR3.offset.Subset.2.x & CDR3.offset.Subset.2.y) {
        
        q1 <- q1 |> mutate(
                `CDR3.offset.Subset.2` = substr(
                `AA JUNCTION`, CDR3.offset.Subset.2.x, CDR3.offset.Subset.2.y
                )
        )
        
    }
    
    
    if(length(subgroup > 0)) {
        
        q1 <- q1 |> filter(
            `IGHV.subgroup` %in% subgroup
        )
        
    }
    
    
    if(ighj.filter != "") {
        
        q1 <- q1 |> filter(
            `IGHJ.gene` == ighj.filter
        )
        
    }
    
    
    if(length(CDR3.IMGT.length.filter > 0)) {
        
        q1 <- q1 |> filter(`CDR3-IMGT length` %in% CDR3.IMGT.length.filter)
        
    }
    
    
    if(CDR3.offset.a.true.filter | CDR3.offset.b.true.filter) {
        
        q1 <- q1 |> filter(`CDR3.offset.a.true` == CDR3.offset.a.true.filter |
                           `CDR3.offset.b.true` == CDR3.offset.b.true.filter  )
        
    }
    
    
    
}



#---------------------Subsets------------------------------

#load path-----------------------------------------
fls    = "C:/Users/nanastasiadou/Desktop/Carmilla 2.0/Subsets_1/Samples/"

fls    = fls |> 
    list.files(full.names = TRUE) 





for (file in fls) {
    
    #load file
    summary_file <- file.path(file, "1_Summary.txt")
    df <- fread(summary_file, sep = "\t", fill = TRUE, select = 1:33)
    
    
    
    # ---------------- subset 1-99---------------------
    
    data <- carmilla(df, IGHV.subgroup = TRUE, vgene = FALSE, jgene = TRUE,
                     CDR3.offset.x = 3,   CDR3.offset.y = 9, 
                     CDR3.offset.a = "QWL", CDR3.offset.b = "",
                     CDR3.offset.Subset.1.x = 5, CDR3.offset.Subset.1.y = 7,
                     CDR3.offset.Subset.2.x = 0, CDR3.offset.Subset.2.y = 0,
                     subgroup =  c("IGHV1", "IGHV5", "IGHV7"), ighj.filter = "IGHJ4",
                     CDR3.IMGT.length.filter = 11:16, 
                     CDR3.offset.a.true.filter = TRUE, 
                     CDR3.offset.b.true.filter = FALSE)
    
    
    
    df1 <- data |>
              mutate(`is.Subset1` = ifelse(
                     `CDR3-IMGT length` == 13 &
                     `CDR3.offset.Subset.1` == "QWL",
                      TRUE,
                      FALSE)) |>
              mutate(`is.Subset99` = ifelse(
                     `CDR3-IMGT length` == 14 &
                     `CDR3.offset.Subset.1` == "QWL", 
                      TRUE,
                      FALSE))
    
    
    
    #Subset1
    Subset_1 <- df1 |>
        filter(`is.Subset1` == "TRUE") |>
        select(1:(ncol(df1) - 10))
    
    
    
    #Subset99
    Subset99 <- df1 |>
        filter(`is.Subset99` == "TRUE") |>
        select(1:(ncol(df1)-10))  
    
    
    #Satellite
    Subset1_99_Sat <- df1 |>
        filter(is.Subset1 == "FALSE",  is.Subset99 == "FALSE") |>
        select(1:(ncol(df1) - 10))
    
    
    #save file
    df_list <- list(Subset1 = Subset_1, Subset99 = Subset99,
                    Satellite = Subset1_99_Sat)
    
    # Specify the file path
    excel_file <- paste0(file, "/Subset_1_99_f.xlsx")
    
    # Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    
    
    
    
    # ---------------- subset 2-169---------------------
    data <- carmilla(df, IGHV.subgroup = TRUE, vgene = TRUE, jgene = TRUE,
                     CDR3.offset.x = 2,   CDR3.offset.y = 6, 
                     CDR3.offset.a = "D", CDR3.offset.b = "E",
                     CDR3.offset.Subset.1.x = 4, CDR3.offset.Subset.1.y = 4,
                     CDR3.offset.Subset.2.x = 0, CDR3.offset.Subset.2.y = 0,
                     subgroup = c("IGHV3"), ighj.filter = "",
                     CDR3.IMGT.length.filter = 7:11, 
                     CDR3.offset.a.true.filter = TRUE, 
                     CDR3.offset.b.true.filter = TRUE) 
    

    df1 <- data |>
             mutate(`is.Subset2` = ifelse(
                    `CDR3-IMGT length` == 9 &
                   (`CDR3.offset.Subset.1` == "D" | 
                    `CDR3.offset.Subset.1` == "E") &
                    `IGHV.gene` == "IGHV3-21",
                     TRUE,
                     FALSE)) |>
             mutate(`is.Subset169` = ifelse(
                    `CDR3-IMGT length` == 9 &
                   (`CDR3.offset.Subset.1` == "D" | 
                    `CDR3.offset.Subset.1` == "E") &
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
        select(1:(ncol(df1) - 9))
    
    
    #save file
    df_list <- list(Subset2 = Subset_2, Subset169 = Subset169, 
                    Satellite = Subset_2_169_Sat)
    
    # Specify the file path
    excel_file <- paste0(file, "/Subset_2_169_without_J_f.xlsx")
    
    # Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    
    
    
    
    
    
    # --------------- subset 8 -----------------------------
    data <- carmilla(df, IGHV.subgroup = TRUE, vgene = TRUE, jgene = TRUE,
                     CDR3.offset.x = 5,   CDR3.offset.y = 14, 
                     CDR3.offset.a = "YSSSWY", CDR3.offset.b = "",
                     CDR3.offset.Subset.1.x = 7, CDR3.offset.Subset.1.y = 12,
                     CDR3.offset.Subset.2.x = 0, CDR3.offset.Subset.2.y = 0,
                     subgroup = c("IGHV2", "IGHV4", "IGHV6"), ighj.filter = "IGHJ5",
                     CDR3.IMGT.length.filter = 16:21, 
                     CDR3.offset.a.true.filter = TRUE, 
                     CDR3.offset.b.true.filter = FALSE) 
    
    
    
    df1 <- data |>
        mutate(`is.Subset.8` = ifelse(
               `CDR3-IMGT length` == 19 &
               `CDR3.offset.Subset.1` == "YSSSWY",
                TRUE,
                FALSE)) |>
        mutate(`is.Subset8B` = ifelse(
               `CDR3-IMGT length` == 18 &
               `CDR3.offset.Subset.1` == "YSSSWY",
                TRUE,
                FALSE))
    
    
    #Subset_8
    Subset_8 <- df1 |>
        filter(`is.Subset.8` == "TRUE") |>
        select(1:(ncol(df1) - 9))
    
    
    #Subset_8B
    Subset_8B <- df1 |>
        filter(`is.Subset8B` == "TRUE") |>
        select(1:(ncol(df1) - 9))
    
    
    
    #Satellite
    
    Subset_8_8B_Sat <- df1 |>
        filter(`is.Subset.8` == "FALSE",
               `is.Subset8B` == "FALSE") |>
        select(1:(ncol(df1) - 9))
    
    
    
    #save file
    df_list <- list(Subset_8 = Subset_8, Subset_8B = Subset_8B,
                    Satellite = Subset_8_8B_Sat)
    
    # Specify the file path
    excel_file <- paste0(file, "/Subset_8_8Β_f.xlsx")
    
    # Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    
    
    
    
    
    
    # --------------- subset 16 --------------------------------
    
    #filtering data frame
    #the FYCS motif of subset 16 is at positions 4-7, with an offset
    #relaxation of +-2 the relevant positions are 2-9 
    data <- carmilla(df, IGHV.subgroup = TRUE, vgene = TRUE, jgene = TRUE,
                     CDR3.offset.x = 3,   CDR3.offset.y = 10, 
                     CDR3.offset.a = "FYCS", CDR3.offset.b = "",
                     CDR3.offset.Subset.1.x = 5, CDR3.offset.Subset.1.y = 8,
                     CDR3.offset.Subset.2.x = 0, CDR3.offset.Subset.2.y = 0,
                     subgroup = c("IGHV2", "IGHV4", "IGHV6"), ighj.filter = "IGHJ6",
                     CDR3.IMGT.length.filter = 22:26, 
                     CDR3.offset.a.true.filter = TRUE, 
                     CDR3.offset.b.true.filter = FALSE) 
    
    
    
    df1 <- data |>
        mutate(`is.Subset16` = ifelse(
               `CDR3-IMGT length` == 24 &
               `CDR3.offset.Subset.1` == "FYCS",
                TRUE,
                FALSE))
    
    
    #Subset_16
    Subset_16 <- df1 |>
        filter(`is.Subset16` == "TRUE") |>
        select(1:(ncol(df1) - 8))
    
    
    #Satellite
    Subset_16_Sat <- df1 |>
        filter(`is.Subset16` == "FALSE") |>
        select(1:(ncol(df1) - 8))
    
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
    data <- carmilla(df, IGHV.subgroup = TRUE, vgene = TRUE, jgene = TRUE,
                     CDR3.offset.x = 7,   CDR3.offset.y = 12, 
                     CDR3.offset.a = "GW", CDR3.offset.b = "",
                     CDR3.offset.Subset.1.x = 9, CDR3.offset.Subset.1.y = 10,
                     CDR3.offset.Subset.2.x = 0, CDR3.offset.Subset.2.y = 0,
                     subgroup = c("IGHV2", "IGHV4", "IGHV6"), ighj.filter = "IGHJ4",
                     CDR3.IMGT.length.filter = 12:16, 
                     CDR3.offset.a.true.filter = TRUE, 
                     CDR3.offset.b.true.filter = FALSE) 
    
    
    
    df1 <- data |>
        mutate(`is.Subset77` = ifelse(
               `CDR3-IMGT length` == 14 &
               `CDR3.offset.Subset.1` == "GW",
                TRUE,
                FALSE))
    
    
    #Subset_77
    Subset_77 <- df1 |>
        filter(`is.Subset77` == "TRUE") |>
        select(1:(ncol(df1) - 8))
    
    
    #Satellite
    Subset_77_Sat <- df1 |>
        filter(`is.Subset77` == "FALSE") |>
        select(1:(ncol(df1) - 8))
    
    
    #save file
    df_list <- list(Subset_77 = Subset_77, Satellite = Subset_77_Sat)
    
    # Specify the file path
    excel_file <- paste0(file, "/Subset_77.xlsx")
    
    # Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    # --------------- subset 3C2,3C3 ------------------------------------
    #filtering data frame
    data <- carmilla(df, IGHV.subgroup = TRUE, vgene = TRUE, jgene = TRUE,
                     CDR3.offset.x = 4,   CDR3.offset.y = 16, 
                     CDR3.offset.a = "DIVVVPAA", CDR3.offset.b = "",
                     CDR3.offset.Subset.1.x = 6, CDR3.offset.Subset.1.y = 13,
                     CDR3.offset.Subset.2.x = 7, CDR3.offset.Subset.2.y = 14,
                     subgroup = c("IGHV3"), ighj.filter = "IGHJ6",
                     CDR3.IMGT.length.filter = 20:24, 
                     CDR3.offset.a.true.filter = TRUE, 
                     CDR3.offset.b.true.filter = FALSE) 
    
    
    
    df1 <- data |>
             mutate(`is.Subset3C2` = ifelse(
                   `CDR3-IMGT length` == 22 &
                   `CDR3.offset.Subset.1` == "DIVVVPAA",
                    TRUE,
                    FALSE)) |>
             mutate(`is.Subset3C3` = ifelse(
                   `CDR3-IMGT length` == 9 &
                   `CDR3.offset.Subset.2` == "DIVVVPAA",
                    TRUE,
                    FALSE))
    
    
    
    
    #Subset3C2
    Subset_3C2 <- df1 |>
        filter(`is.Subset3C2` == "TRUE") |>
        select(1:(ncol(df1) - 10))
    
    #Subset3C3
    Subset_3C3 <- df1 |>
        filter(`is.Subset3C3` == "TRUE") |>
        select(1:(ncol(df1) - 10))  
    
    
    
    #Satellite
    Subset_3C2_3C3_sat <- df1 |>
        filter(`is.Subset3C2` == "FALSE",
               `is.Subset3C3` == "FALSE") |>
        select(1:(ncol(df1) - 10))
    
    
    #save file
    df_list <- list(Subset_3C2 = Subset_3C2, Subset_3C3 = Subset_3C3, 
                    Satellite = Subset_3C2_3C3_sat)
    
    
    # Specify the file path
    excel_file <- paste0(file, "/Subset_3C2_3C3.xlsx")
    
    # Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    
    
    
    
    
    
    
    # --------------- subset 4 ---------------------
    
    #filtering data frame
    data <- carmilla(df, IGHV.subgroup = TRUE, vgene = TRUE, jgene = TRUE,
                     CDR3.offset.x = 10,   CDR3.offset.y = 19, 
                     CDR3.offset.a = "KRYYYY", CDR3.offset.b = "RRYYYY",
                     CDR3.offset.Subset.1.x = 12, CDR3.offset.Subset.1.y = 17,
                     CDR3.offset.Subset.2.x = 0, CDR3.offset.Subset.2.y = 0,
                     subgroup = c("IGHV2", "IGHV4", "IGHV6"), ighj.filter = "IGHJ6",
                     CDR3.IMGT.length.filter = 18:22, 
                     CDR3.offset.a.true.filter = TRUE, 
                     CDR3.offset.b.true.filter = TRUE) 
    
    df1 <- data |>
        mutate(`is.Subset4` = ifelse(
               `CDR3-IMGT length` == 20 &
               `CDR3.offset.Subset.1` == "KRYYYY" | 
               `CDR3.offset.Subset.1` == "RRYYYY",
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
    excel_file <- paste0(file, "/Subset_4_f.xlsx")
    
    # Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    
    
    
}
