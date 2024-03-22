#load path-----------------------------------------
fls    = "C:/Users/nanastasiadou/Desktop/Carmilla 2.0/Subsets/Samples/"

fls    = fls |> 
    list.files(full.names = TRUE) 


#function----------------------
carmilla <- function(
        
        q, IGHV.subgroup = FALSE, vgene = FALSE, jgene = FALSE, CDR3.offset = FALSE,
        CDR3.offset.x = NULL, CDR3.offset.y = NULL, CDR3.offset.lb = "[A-Z]+",
        CDR3.offset.Subset.x = NULL, CDR3.offset.Subset.y = NULL
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
    
    
    if(CDR3.offset.x & CDR3.offset.y) {

        q1 <- q1 |> mutate(
            `CDR3.offset` = substr(`AA JUNCTION`, CDR3.offset.x, CDR3.offset.y)
        )

    }
    

    
    if(CDR3.offset.lb != "") {

        q1 <- q1 |> mutate(
            `CDR3.offset.true` = grepl(CDR3.offset.lb, `CDR3.offset`)
        )

    }

    if(CDR3.offset.Subset.x & CDR3.offset.Subset.y) {
        
        q1 <- q1 |> mutate(
            `CDR3.offset.Subset1` = substr(
                      `AA JUNCTION`, CDR3.offset.Subset.x, CDR3.offset.Subset.y
                    )
        )
        
    }
    
}





for (file in fls) {
    
    #load file
    summary_file <- file.path(file, "1_Summary.txt")
    df <- fread(summary_file, sep = "\t", fill = TRUE, select = 1:33)
    
    
    # --------------- subset 1-99---------------------
    data <- carmilla(df, IGHV.subgroup = TRUE, vgene = TRUE, jgene = TRUE,
                     CDR3.offset.x = 3,   CDR3.offset.y = 9, 
                     CDR3.offset.lb = "QWL",
                     CDR3.offset.Subset.x = 5, CDR3.offset.Subset.y = 7)
    
    
    df1 <- df |>
        filter(`V-DOMAIN Functionality` == "productive" |
               `V-DOMAIN Functionality`== "productive (see comment)",
               `CDR3-IMGT length` != "X") |>
        mutate(`IGHV.subgroup` = substr(`V-GENE and allele`, 8, 12),
               `IGHV.gene` = substr(`V-GENE and allele`, 8, 15),
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
