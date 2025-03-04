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
  

selected_names <- c("EP", "EP_SE")
for(name in selected_names) {
  latex_code <- xtable(elasticities_result[[name]])
  file_name <- paste0("table_", name, ".tex")
  print(latex_code, type = "latex", file = file_name)
}

# one matrix per separate object
for(name in selected_names) {
  assign(name, elasticities_result[[name]])
}

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

format_two_digits <- function(number) {
  return(sprintf("%.2f", number))
}

for (i in 1:nrow(EP)) {
  for (j in 1:ncol(EP)) {
    EP[i, j] <- format_two_digits(as.numeric(EP[i, j]))
    EP_SE[i, j] <- format_two_digits(as.numeric(EP_SE[i, j]))
  }
}

EP_starred <- EP
for (i in 2:nrow(EP)) {
  for (j in 2:ncol(EP)) {
    EP_starred[i, j] <- add_stars(as.numeric(EP[i, j]), as.numeric(EP_SE[i, j]))
  }
}

fusionner_dataframes <- function(df1, df2) {
  rows <- nrow(df1)
  result <- data.frame(matrix(nrow = rows * 2, ncol = ncol(df1), 
                              dimnames = list(NULL, colnames(df1))))
  for (i in 1:rows) {
    result[i * 2 - 1, ] <- df1[i, ]
    result[i * 2, ] <- paste0("(", sprintf("%.2f", as.numeric(df2[i, ])), ")")
  }
  return(result)
}

# Utilisez la fonction pour fusionner les dataframes
EP_concat <- fusionner_dataframes(EP_starred, EP_SE)
latex_code <- xtable(EP_concat)
# Affichez le dataframe résultant
print(EP_concat)

