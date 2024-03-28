


rm(list = ls())
gc()


library(data.table)
library(dplyr)




#function-----------------------------------------
carmilla <- function(
        
    q, IGHV.subgroup = FALSE, vgene = FALSE, jgene = FALSE, 
    CDR3.offset.x = NULL, CDR3.offset.y = NULL, 
    CDR3.offset.a = "[A-Z]+", CDR3.offset.b = "[A-Z]+",
    CDR3.offset.Subset.x = NULL, CDR3.offset.Subset.y = NULL,
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
    
    
    if(CDR3.offset.Subset.x & CDR3.offset.Subset.y) {
        
        q1 <- q1 |> mutate(
                `CDR3.offset.Subset` = substr(
                `AA JUNCTION`, CDR3.offset.Subset.x, CDR3.offset.Subset.y)
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
                     CDR3.offset.Subset.x = 5, CDR3.offset.Subset.y = 7,
                     subgroup =  c("IGHV1", "IGHV5", "IGHV7"), ighj.filter = "IGHJ4",
                     CDR3.IMGT.length.filter = 11:16, 
                     CDR3.offset.a.true.filter = TRUE, 
                     CDR3.offset.b.true.filter = FALSE)
    
    
    
    df1 <- data |>
              mutate(`is.Subset1` = ifelse(
                     `CDR3-IMGT length` == 13 &
                     `CDR3.offset.a.true` == "QWL",
                      TRUE,
                      FALSE)) |>
              mutate(`is.Subset99` = ifelse(
                     `CDR3-IMGT length` == 14 &
                     `CDR3.offset.a.true` == "QWL", 
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
    excel_file <- paste0(file, "/Subset_1_99_f.xlsx")
    
    # Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    
    
    
    
    # ---------------- subset 2-169---------------------
    data <- carmilla(df, IGHV.subgroup = TRUE, vgene = TRUE, jgene = TRUE,
                     CDR3.offset.x = 2,   CDR3.offset.y = 6, 
                     CDR3.offset.a = "D", CDR3.offset.b = "E",
                     CDR3.offset.Subset.x = 4, CDR3.offset.Subset.y = 4,
                     subgroup = c("IGHV3"), ighj.filter = "",
                     CDR3.IMGT.length.filter = 7:11, 
                     CDR3.offset.a.true.filter = TRUE, 
                     CDR3.offset.b.true.filter = TRUE) 
    
    
    df1 <- data |>
             mutate(`is.Subset2` = ifelse(
                    `CDR3-IMGT length` == 9 &
                   (`CDR3.offset.Subset` == "D" | 
                    `CDR3.offset.Subset` == "E") &
                    `IGHV.gene` == "IGHV3-21",
                     TRUE,
                     FALSE)) |>
             mutate(`is.Subset169` = ifelse(
                    `CDR3-IMGT length` == 9 &
                   (`CDR3.offset.Subset` == "D" | 
                    `CDR3.offset.Subset` == "E") &
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
    excel_file <- paste0(file, "/Subset_2_169_without_J_f.xlsx")
    
    # Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    
    
    
    
    
    
    # --------------- subset 8 -----------------------------
    data <- carmilla(df, IGHV.subgroup = TRUE, vgene = TRUE, jgene = TRUE,
                     CDR3.offset.x = 5,   CDR3.offset.y = 14, 
                     CDR3.offset.a = "YSSSWY", CDR3.offset.b = "",
                     CDR3.offset.Subset.x = 7, CDR3.offset.Subset.y = 12,
                     subgroup = c("IGHV2", "IGHV4", "IGHV6"), ighj.filter = "IGHJ5",
                     CDR3.IMGT.length.filter = 16:21, 
                     CDR3.offset.a.true.filter = TRUE, 
                     CDR3.offset.b.true.filter = FALSE) 
    
    
    
    df1 <- data |>
        mutate(`is.Subset.8` = ifelse(
               `CDR3-IMGT length` == 19 &
               `CDR3.offset.Subset` == "YSSSWY",
                TRUE,
                FALSE)) |>
        mutate(`is.Subset8B` = ifelse(
               `CDR3-IMGT length` == 18 &
               `CDR3.offset.Subset` == "YSSSWY",
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
    excel_file <- paste0(file, "/Subset_8_8Β_f.xlsx")
    
    # Write the list of data frames to an Excel file with different sheets
    writexl::write_xlsx(df_list, path = excel_file)
    
    
    
    
    
    
    
    
    
}
