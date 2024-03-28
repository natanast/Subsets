


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
fls    = "C:/Users/natan/Υπολογιστής/Carmilla 2.0/Samples"

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
                  (`CDR3.offset.a.true` == "D" | 
                   `CDR3.offset.b.true` == "E") &
                   `IGHV.gene` == "IGHV3-21",
                    TRUE,
                    FALSE)) |>
            mutate(`is.Subset169` = ifelse(
                   `CDR3-IMGT length` == 9 &
                  (`CDR3.offset.a.true` == "D" | 
                   `CDR3.offset.b.true` == "E") &
                   `IGHV.gene` == "IGHV3-48",
                    TRUE,
                    FALSE))
    
    
    
    
    
    
    # --------------- subset 8 -----------------------------
    data <- carmilla(df, IGHV.subgroup = TRUE, vgene = TRUE, jgene = TRUE,
                     CDR3.offset.x = 5,   CDR3.offset.y = 14, 
                     CDR3.offset.a = "YSSSWY", CDR3.offset.b = "",
                     CDR3.offset.Subset.x = 5, CDR3.offset.Subset.y = 14,
                     subgroup = c("IGHV2", "IGHV4", "IGHV6"), ighj.filter = "IGHJ5",
                     CDR3.IMGT.length.filter = 16:21, 
                     CDR3.offset.a.true.filter = TRUE, 
                     CDR3.offset.b.true.filter = FALSE) 
    
    
    
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
    
    
    
}