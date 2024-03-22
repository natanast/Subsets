

carmilla <- function(
        q, vgene = FALSE, jgene = FALSE, 
        CDR3.offset.x = 8, CDR3.offset.y = 12, CDR3.offset.lb = "QWL",
        CDR3.offset2.x = NULL, CDR3.offset2.y = NULL, CDR3.offset2.lb = NULL
) {
    
    
    q1 <- q |>
            filter(
                
                `V-DOMAIN Functionality` == "productive" |
                `V-DOMAIN Functionality`== "productive (see comment)",
                `CDR3-IMGT length` != "X"
                
              ) 
    
    if(vgene) {
        
        q1 <- q1 |> mutate(`IGHV.gene` = substr(`V-GENE and allele`, 8, 15))
        
    }
    
    
    if(jgene) {
        
        q1 <- q1 |> mutate(`IGHJ.gene` = substr(`J-GENE and allele`, 8, 12))
        
    }
    
    if(CDR3.offset.x) {
        q1 <- q1 |> mutate(`CDR3.offset` = substr(`AA JUNCTION`, 3, 9))
    }
    
    if(CDR3.offset.y) {
        q1 <- q1 |> mutate(`CDR3.offset2` = substr(`AA JUNCTION`, 3, 9))
    }
    
    if(CDR3.offset.y) {
        q1 <- q1 |> mutate(`CDR3.offset2` = substr(`AA JUNCTION`, 3, 9))
    }
    
    
    
    df1 <-  |>
        
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
    
    
}