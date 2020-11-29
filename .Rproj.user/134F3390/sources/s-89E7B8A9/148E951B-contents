#######################
## Analysis of       ##
##  Covid            ##
##  recovery rates   ##
##    and            ##
##  Population       ##
##                   ##
##                   ##
##                   ##
##      NO. 3        ##
##                   ##
## Analysis of       ##
#       the data     ##
##                   ##
#######################


# Clear memory
rm(list=ls())

# Packages to use
library(tidyverse)
# For scaling ggplots
require(scales)
# Estimate piecewise linear splines
#install.packages("lspline")
library(lspline)
# Estimate robust SE
#install.packages("estimatr")
library(lmtest)
## Check 

library(estimatr)
# Compare models with robust SE
#install.packages("texreg")
library(texreg)
# For different themes
#install.packages(ggthemes)
library(ggthemes)


# Call the data from local file
data_path <- "data/clean/"
df <- read_csv(paste0(data_path,'covid_pop_10_09_2020_clean.csv'))


#### Aim of the Analysis #######
### Research Question: 
### What is the impact of the population size on my covid recovery rates?

### Variables Description
### recovery_rate (explanatory): describes the ratio of recovered cases divided by the total confirmed cases on a given day.
### population_mi (dependent): population size in million units by country - 2019
### Sample: Recovery Rates of a single day / Population: Total recovery rates whole Covid season
### Possible Quality Issues: Reliability (Data integrity - Governmental Interests) , Content (Caused by Covid or simply with Covid at death moment?)


####
# 
# Quick check on HISTOGRAMS and Boxplots

######## Recovery Rates ###############

summary(df$recovery_rate)

hist(df$recovery_rate, main = "Covid Recovery rates 10-09-20", xlab="Recovery rates(%)", 
     ylab = "Number of Countries", border="blue", col="lightgreen", xlim=c(0,1), las=1, 
     breaks=10)

boxplot(df$recovery_rate, main = "Covid Recovery rates (%)", xlab = "Recovery Rates(%)",
        col = "orange", border = "brown", horizontal = TRUE, notch = TRUE)

## Left tail Distribution: median > mean, outliers distributed on the left part of the graph 

############ Population ##########################

summary(df$population_mi)

hist(df$population_mi, main = "Population Size", xlab="Number of Citizens (mi)", 
     ylab = "Number of Countries", border="blue", col="lightgreen", las=1, 
     breaks=10)

boxplot(df$population_mi, main = "Population Size", xlab = "Number of Citizens (mi)",
        col = "grey", border = "brown", horizontal = TRUE, notch = TRUE)

## Right tail Distribution: median < mean, outliers distributed on the right part of the graph


## Remove all countries with low confirmed cases and really low population - 
## not relevant for recover analysis virus not sufficiently spread

df <- df %>% filter(confirmed > 5000,
                            population > 100000)

######
# Check basic scatter-plots!
# 
# Where to use log-transformation? - level-level vs level-log vs log-level vs log-log
#
# 1) lifeexp - gdptot: level-level model without scaling
ggplot( df , aes(x = population_mi , y = recovery_rate)) +
  geom_point() +
  geom_smooth(method="loess")+
  labs(x = "Population Size (mi)", y = "Recovery Rates (%)") +
  ylim(0,1)
  
# You can change the scale for Total GDP for checking log-transformation
ggplot( df , aes(x = population_mi , y = recovery_rate)) +
  geom_point() +
  geom_smooth(method="loess") +
  labs(x = " Population Size (mi) , ln scale )",y = "Recovery Rates (%)") + ylim(0,1) +
  scale_x_continuous( trans = log_trans(),  breaks = c(1,2,5,10,20,50,100,200,500,1000))

# You can change the scale for Total GDP and life-expectancy for checking log-transformation
ggplot( df , aes(x = population_mi, y = recovery_rate ))  +
  geom_point() +
  geom_smooth(method="loess")+
  labs(x = " Population Size (mi) , ln scale )",y = "Recovery Rates (%), ln scale)") +
  ylim(0,1) +
  scale_x_continuous( trans = log_trans(),  breaks = c(1,2,5,10,20,50,100,200,500,1000) )+
  scale_y_continuous( trans = log_trans())


####
# Conclusions:
#   1) taking log of population_mi is needed, it's a highly right-tailed skewed distribution:
#       - Better fit, easier to interpret. 
#   2) using log of recovery_rate is possible, however the level are easier to interpret 
#       - Substantive: Log changes is harder to interpret and our aim is to get a direct % based comparison
#       - Statistical: The patterns are better visible on level


df <- df %>% mutate(ln_population_mi = log(population_mi))

######
# Make some models:
#   w ln_population_mi:
#     reg1: recovery_rate = alpha + beta * ln_population_mi
#     reg2: recovery_rate = alpha + beta_1 * ln_population_mi + beta_2 * ln_population_mi^2
#     reg3: recovery_rate = alpha + beta_1 * ln_population_mi + beta_2 * ln_population_mi^2 + beta_3 * ln_population_mi^3


###
# Two ways to handle polynomials: 
#
# 1) Add powers of the variable(s) to the dataframe:
df <- df %>% mutate( ln_population_mi_sq = ln_population_mi^2,
                     ln_population_mi_cb = ln_population_mi^3)

# Regressions
#
# Built in regression in R
reg_b <- lm(recovery_rate ~ ln_population_mi , data = df )
summary( reg_b )

### Check homoskedastic - Breusch-Pagan Test
bptest( reg_b )
### As showed into the loess graph there is a high chance of Heteroscedasticity
### p-value is really big, so we reject the null hypothesis (homoscedasticity). 
### So we should choose the robust model

# First model:
reg1 <- lm_robust( recovery_rate ~ ln_population_mi , data = df , se_type = "HC2" )
summary( reg1 )

# Visual inspection:
ggplot( data = df, aes( x = ln_population_mi , y = recovery_rate ) ) + 
  geom_point( color='blue') +
  geom_smooth( method = lm , color = 'red' )

# Second and third model with gdptot
reg2 <- lm_robust( recovery_rate ~ ln_population_mi + ln_population_mi_sq , data = df )
summary( reg2 )
ggplot( data = df, aes( x = ln_population_mi , y = recovery_rate ) ) + 
  geom_point( color='blue') +
  geom_smooth( formula = y ~ poly(x,2) , method = lm , color = 'red' )

reg3 <- lm_robust( recovery_rate ~ ln_population_mi + ln_population_mi_sq + ln_population_mi_cb , data = df )
summary ( reg3 )
ggplot( data = df, aes( x = ln_population_mi , y = recovery_rate ) ) + 
  geom_point( color='blue') +
  geom_smooth( formula = y ~ poly(x,3) , method = lm , color = 'red' )


# Regression with piecewise linear spline:
# 1st define the cutoff for 
cutoff <- 10
# 2nd we use a log transformation -> cutoff needs to be transformed as well
cutoff_ln <- log( cutoff )
# Use simple regression with the lspline function
reg6 <- lm_robust(recovery_rate ~ lspline( ln_population_mi , cutoff_ln ), data = df )
summary( reg6 )
ggplot( data = df, aes( x = ln_population_mi, y = recovery_rate ) ) + 
  geom_point( color='blue') +
  geom_smooth( formula = y ~ lspline(x,cutoff_ln) , method = lm , color = 'red' )

#####
# Creating model summary with texreg
data_out <- "out/"
htmlreg( list(reg1 , reg2 , reg3 , reg6),
         type = 'html',
         custom.model.names = c("Population Size - linear","Population Size - quadratic","Population Size - cubic",
                                "GDP/capita - PLS"),
         caption = "Modelling recovery_rates of countries",
         file = paste0( data_out ,'model_comparison.html'), include.ci = FALSE)

######
# Based on model comparison our chosen model is reg2 - recovery_rate ~ ln_gdppc +  
#   Substantive: - level-log interpretation works properly for countries population
#                - magnitude of population coefficient is meaningful
#                - Lower R-Adjusted square, model is significant however improper for prediction
#
#   Statistical: - simple model, easy to interpret
#                - Comparatively high R2 and captures variation well


#################################
## Testing hypothesis
#

##
# 1) Coefficient is equal to 0:
# Implemented by default...
summary( reg2 )

# Let test: H0:  = 0, HA: ln_population_mi neq 5
linearHypothesis( reg2 , "ln_population_mi = 0")
# Can reject null hypothesis:  equal to 0 

######
# Residual analysis.

# Get the predicted y values from the model
df$reg2_y_pred <- reg2$fitted.values
# Calculate the errors of the model
df$reg2_res <- df$recovery_rate - df$reg2_y_pred 

# Find countries with largest negative errors
df %>% top_n( -5 , reg2_res ) %>% 
  select( country , recovery_rate , reg2_y_pred , reg2_res )

# Find countries with largest positive errors
df %>% top_n( 5 , reg2_res ) %>% 
  select( country , recovery_rate , reg2_y_pred , reg2_res )


#################################
## Prediction uncertainty
#

# CI of predicted value/regression line is implemented in ggplot
ggplot( data = df, aes( x = ln_population_mi , y = recovery_rate ) ) + 
  geom_point( color='blue') +
  geom_smooth( method = lm , color = 'red' , se = T )


# CI of regression line
pred2_CI <- predict( reg2, newdata = df , interval ="confidence" , alpha = 0.05 )
pred2_CI


# Hand made CI for regression line
# 1) Add to datatset:
df <- df %>% mutate( CI_reg2_lower = pred2_CI$fit[,2],
                     CI_reg2_upper = pred2_CI$fit[,3] )
# 2) Plot
ggplot(  ) + 
  geom_point( data = df, aes( x = ln_population_mi , y = recovery_rate ) , color='blue') +
  geom_line( data = df, aes( x = ln_population_mi, y = reg2_y_pred ) , color = 'red' , size = 1 ) +
  geom_line( data = df, aes( x = ln_population_mi, y = CI_reg2_lower ) , color = 'green' ,
             size = 1 , linetype = "dashed" ) +
  geom_line( data = df, aes( x = ln_population_mi, y = CI_reg2_upper ) , color = 'black' ,
             size = 1 , linetype = "dashed" ) +
  labs(x = "ln( Population Size (mi) )",y = "Recovery Rates (%)") 


# As expected, due the high R2 considered in the model, the predicted results weren´t a good fit for the model


##
# Prediction intervals!
#
pred2_PI <- predict( reg2, newdata = df , interval ="prediction" , alpha = 0.05 )

# Hand made Prediction Interval for regression line

df <- df %>% mutate( PI_reg2_lower = pred2_PI$fit[,2],
                     PI_reg2_upper = pred2_PI$fit[,3] )
# 2) Plot
ggplot(  ) + 
  geom_point( data = df, aes( x = ln_population_mi, y = recovery_rate ) , color='blue') +
  geom_line( data = df, aes( x = ln_population_mi, y = reg2_y_pred ) , color = 'red' , size = 1 ) +
  geom_line( data = df, aes( x = ln_population_mi, y = PI_reg2_lower ) , color = 'green' ,
             size = 1 , linetype = "dotted" ) +
  geom_line( data = df, aes( x = ln_population_mi, y = PI_reg2_upper ) , color = 'black' ,
             size = 1 , linetype = "dotted" ) +
  labs(x = "ln( GDP/capita, 2018 int. const. $, PPP)",y = "Life expectancy  (years)") 

##### Conclusion 

### The population size demonstrated itself significant to the recovery rates on the following dataset. 
### However, the impact alone of the variable isn´t not sufficient to describe the complexity of the recovery rates,
### other variables, should be added to the model for a better understanding of the recovery_rates

