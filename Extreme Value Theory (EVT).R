# Load the required libraries
library(xts)
library(readxl)
library(extRemes)

# Read the data from the Excel file (assuming the data contains annual maxima rainfall)
data <- read_excel("usman.xlsx", sheet = 1)
data <- data[, -c(2)] 

# Define the number of years in each time slice (30 years in this case)
moving_window <- 30

# Define the total number of years in the data
total_years <- nrow(data)

# Define a vector to store the results
return_levels <- numeric()


# Loop through the data and fit GEV distributions to each time slice
for (start_year in 1921:(1921 + total_years - moving_window)) {
  last_year <- start_year + moving_window - 1
  subset_years <- data[(start_year - 1921 + 1):(last_year - 1921 + 1), ]
  
  # Fit GEV distribution based on L-moments estimation for the current time slice
  fit_lmom <- fevd(as.vector(subset_years$B028), method = "Lmoments", type = "GEV")
  
  # Calculate the 100-year return level for the current time slice
  rl_lmom <- return.level(fit_lmom, conf = 0.05, return.period = c(100))
  
  # Store the return level in the results vector
  return_levels <- c(return_levels, rl_lmom)
}

write.csv(return_levels, file = "30-moving window RL.csv")

#Read Trend and Senslope from 30 Year Moving Average 

dataset <- read.csv("30year moving window-RLs.csv")
colnames(dataset)<- c("Year","B006","B008")
db<- dataset$B006

# Mann-Kendall Test
mk_result <- mk.test(db, alternative = "two.sided", continuity = TRUE)
print("Mann-Kendall Test:")
print(mk_result)

# Sen's Slope Estimation
sens_slope_result <- sens.slope(db)
print("Sen's Slope:")
print(sens_slope_result)