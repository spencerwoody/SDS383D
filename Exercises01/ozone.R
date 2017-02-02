# Taken from Prof James Scott

# Load the library
# you might have to install this the first time
library(mlbench)

# Load the data
ozone = data(Ozone, package='mlbench')

# Look at the help file for details
?Ozone

#  1  Month: 1 = January, ..., 12 = December
#  2  Day of month
#  3  Day of week: 1 = Monday, ..., 7 = Sunday
#  4  Daily maximum one-hour-average ozone reading
#  5  500 millibar pressure height (m) measured at Vandenberg AFB
#  6  Wind speed (mph) at Los Angeles International Airport (LAX)
#  7  Humidity (%) at LAX
#  8  Temperature (degrees F) measured at Sandburg, CA
#  9  Temperature (degrees F) measured at El Monte, CA
# 10  Inversion base height (feet) at LAX
# 11  Pressure gradient (mm Hg) from LAX to Daggett, CA
# 12  Inversion base temperature (degrees F) at LAX
# 13  Visibility (miles) measured at LAX

# Scrub the missing values
# Extract the relevant columns 
ozone = na.omit(Ozone)[,4:13]

y = ozone[,1]
x = as.matrix(ozone[,2:10])

# add an intercept
x = cbind(1,x)

# compute the estimator
betahat = solve(t(x) %*% x) %*% t(x) %*% y

# Fill in the blank
# betacov = ?

# Now compare to lm
# the 'minus 1' notation says not to fit an intercept (we've already hard-coded it as an extra column)
lm1 = lm(y~x-1)

summary(lm1)
betacovlm = vcov(lm1)
sqrt(diag(betacovlm))