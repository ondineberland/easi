## Ondine's notes ##

## Step by Step

### 0. Data
#Load data 

### 1. Modules
#Run easi.R
    ## get 'est' the results of the estimation
#Run intermediate blocs.R 
    ## In summary, the intermediate.blocs function processes an input object to recover and calculate various coefficients and sums related to equations. These calculations are based on the "EASI made EASIER" methodology by PENDAKUR (2008). The function returns a list containing all the recovered and calculated values.##
#Run elastic.R
    ## then run the following code to export elasticities into latex Tables
## ELASTICITIES 
elasticities_result <- elastic(est, type = c("price", "income", "demographics"))

  # Convert the matrix to LaTeX format

selected_names <- c("EP", "EP_SE", "EPS", "EPQ", "ELASTPRICE", "ELASTINCOME")
  for(name in selected_names) {
  
   latex_code <- xtable(elasticities_result[[name]])
   file_name <- paste0("table_", name, ".tex")
   print(latex_code, type = "latex", file = file_name)
  }

  # one matrix per separate object
  for(name in selected_names) {
    assign(name, elasticities_result[[name]])
  }


### NEXT TEST



 # A function that adds stars to the coefficients
add_stars <- function(coef, se) {
  z <- abs(coef/se)
  
  # Calculate p-value for two-tailed z-test
  p_value <- 2 * (1 - pnorm(z))
  
  # Add stars based on p-value
  if (p_value < 0.01) {
    return(paste0(coef, "***"))
  } else if (p_value < 0.05) {
    return(paste0(coef, "**"))
  } else if (p_value < 0.1) {
    return(paste0(coef, "*"))
  } else {
    return(coef)
  }
}

EP_starred <- EP
for (i in 2:nrow(EP)) {
  for (j in 2:ncol(EP)) {
    EP_starred[i, j] <- add_stars(as.numeric(EP[i, j]), as.numeric(EP_SE[i, j]))
  }
}

### THEN I HAVE A PROBLEM BECAUSE EP_starred is a 9x9 matrix but appears weird -> find someone to help me 
merged_matrix <- data.frame()
