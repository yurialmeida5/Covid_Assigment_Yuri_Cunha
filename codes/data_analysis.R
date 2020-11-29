#######################
## Analysis of       ##
## Number            ##
##    of Registered  ##
##  Deaths           ##
##                   ##    
##    and            ##
##  Registered Cases ##
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
#Linear Hypothesis
library(car)


# Call the data from local file
data_path <- "data/clean/"
df <- read_csv(paste0(data_path,'covid_pop_09_29_2020_clean.csv'))


#### Aim of the Analysis #######
### Research Question: 
### Is there any association between number of deaths and confirmed cases on this dataset? If yes, what is it?
### Variables Description
### death (Explanatory): number of accumulated covid deaths until a given day by country    
### confirmed (Dependent): number of accumulated confirmed covid cases until a given day by country 
### Sample: number of cases and deaths (until 29/09)  / Population: Total number of cases and deaths whole Covid season
### Possible Quality Issues: Reliability (Data integrity - Governmental Interests) , Content (Caused by Covid or simply with Covid at death moment?)


####
#### 
# Quick check on HISTOGRAMS and Boxplots

######## Deaths ###############

summary(df$death)

hist(df$death, main = "Covid Deaths (until 09-29-20)", xlab="Number of Deaths", 
     ylab = "Frequency", border="blue", col="lightgreen", las=1, 
     breaks=10)

## Right tail Distribution: median < mean, outliers distributed on the right part of the graph

############ Confirmed Cases ##########################

summary(df$confirmed)

hist(df$confirmed, main = "Covid Confirmerd Cases (until 09-29-20)", xlab="Number of Confirmed Cases", 
     ylab = "Frequency", border="blue", col="lightgreen", las=1, 
     breaks=10)

## Right tail Distribution: median < mean, outliers distributed on the right part of the graph


## Remove all countries with low confirmed cases and really low population - not relevant for analysis. 
## Virus not sufficiently spread, and also helps to eliminate countries with specific geographic characteristics.
## Ex: Brunei, Comoros, Dominica, Fiji and some other islands.

df <- df %>% filter(confirmed > 5000,
                            population > 100000)


## Log Transformations

### Creating the new variables

df <- df %>% mutate(ln_death = log(death),
                    ln_confirmed = log(confirmed))

### Check homoskedastic - Breusch-Pagan Test
lvl_lvl <- lm( death ~ confirmed , data = df)
log_lvl <- lm( death ~ ln_confirmed , data = df )
lvl_log <- lm( ln_death ~ confirmed , data = df )
log_log <- lm( ln_death ~ ln_confirmed , data = df )

bptest( lvl_lvl ) # small p-value can´t reject the null hypothesis (homoscedasticity)
bptest( log_lvl ) # small p-value can´t reject the null hypothesis (homoscedasticity)
bptest( lvl_log ) # big p-value can t reject the null hypothesis (homoscedasticity)  
bptest( log_log ) # big p-value can reject the null hypothesis (homoscedasticity)

# Use lm_robust lvl_log / log_log
lvl_log <- lm_robust( ln_death ~ confirmed , data = df , se_type = "HC2" )
log_log <- lm_robust( ln_death ~ ln_confirmed , data = df , se_type = "HC2" )

#####
# Log Transformation Summary texreg
data_out <- "out/"
htmlreg( list(lvl_lvl , log_lvl , lvl_log,  log_log ),
         type = 'html',
         custom.model.names = c("lvl_lvl", "log_lvl","lvl_log",
                                "log_log"),
         caption = "Log Transformations Comparison",
         file = paste0( data_out ,'logtransformation_comparison.html'), include.ci = FALSE)
####
# Conclusions:
#   Comparing the p-values (significance) and the R^2 adjusted (highness) for the four types of possible transformations,
#   it´s possible to check that the lvl_lvl and log_log represented the best possible models to choose.

# Graphical Check
# Death - confirmed : level-level 

ggplot( df , aes(x = confirmed , y = death)) +
  geom_point() +
  geom_smooth(method="loess")+
  labs(x = " Number of Confirmed Cases", y = "Number of Deaths") 

# Death - confirmed : log-log
ggplot( df , aes(x = confirmed, y = death ))  +
  geom_point() +
  geom_smooth(method="loess")+
  labs(x = " Number of Confirmed Cases , ln scale )",y = "Number of Deaths, ln scale)") +
  scale_x_continuous( trans = log_trans(),  breaks = c(10000,50000,1000000))+
  scale_y_continuous( trans = log_trans())

####
# Conclusions:
#  Although the models lvl_lvl and log_log showed good p and R^2-Ajusted values,
#  checking the graphs its possible to determine that taking logs will be needed for both cases.
#       Both deaths and confirmed cases are highly right-tailed skewed distributions with a few outliers.
#       Although, level is easier to interpret in general, the logs in this case produce a better fit, 
#       making the association pattern between them better to identify and analyze. 
#       Absolute values are not that relevant in both cases.

# Remove not necessary variables
rm(log_log,log_lvl,lvl_lvl,lvl_log)

## Estimating Different Models

reg1 <- lm_robust( ln_death ~ ln_confirmed , data = df , se_type = "HC2" )
summary( reg1 )

# Visual inspection:
ggplot( data = df, aes( x = ln_death, y = ln_confirmed ) ) + 
  geom_point( color='blue') +
  geom_smooth( method = lm , color = 'red' )


### Quadratic Linear Regression

### Add quadratic variable to the dataframe

df <- df %>% mutate(ln_confirmed_sq = ln_confirmed^2)

# Second and third model with gdptot
reg2 <- lm_robust( ln_death ~ ln_confirmed + ln_confirmed_sq , data = df )
summary( reg2 )

# Visual inspection:
ggplot( data = df, aes( x = ln_confirmed , y = ln_death ) ) + 
  geom_point( color='blue') +
  geom_smooth( formula = y ~ poly(x,2) , method = lm , color = 'red' )

# Piecewise linear spline:

cutoff <- 400000
cutoff_ln <- log( cutoff )

reg3 <- lm_robust(ln_death ~ lspline( ln_confirmed , cutoff_ln ), data = df )
summary( reg3 )

ggplot( data = df, aes( x = ln_confirmed, y = ln_death ) ) + 
  geom_point( color='blue') +
  geom_smooth( formula = y ~ lspline(x,cutoff_ln) , method = lm , color = 'red' )

# Weighted-OLS

reg4 <- lm_robust(ln_death ~ ln_confirmed, data = df , weights = population)
summary( reg4 )

ggplot(data = df, aes(x = ln_confirmed, y = ln_death)) +
  geom_point(data = df, aes(size=population),  color = 'blue', shape = 16, alpha = 0.6,  show.legend=F) +
  geom_smooth(aes(weight = population), method = "lm", color='red')+
  scale_size(range = c(1, 15)) +
  labs(x = "ln(Number of Confirmed Cases) ",y = "ln(Number of Deaths)")
  
#####
# Model summary with texreg
data_out <- "out/"
htmlreg( list(reg1 , reg2 , reg3 , reg4),
         type = 'html',
         custom.model.names = c("Linear Regression", "Quadratic (Linear) Regression","Piecewise Linear Spline",
                                "Weighted-OLS"),
         caption = "Modelling number of deaths by covid based on number of confirmed cases",
         file = paste0( data_out ,'model_comparison.html'), include.ci = FALSE)

######
# Based on model comparison our chosen model is reg4  - Weighted-OLS   
#   Substantive: - The model demonstrated a really great fit for the dataset, and can be used for prediction 
#                - The beta parameter states that for the following dataset an increase of 10% on the number of confirmed cases
#                  will increase, on average, the number of deaths in 9.4%
#                - The alfa in that case represents the mean effect on number of deaths of variables not included in the model.
#
#           
#
#   Statistical: - The models showed really similar statistical results (p-values and r-adjusted), 
#                   capturing the variation really well. 
#                - As expected, the weighted one by population size could demonstrate a even better fit.
#                  the number of 


# Remove non-used variables

rm(cutoff,cutoff_ln,reg1,reg2,reg3)


#################################
## Testing hypothesis
#

##
# 1) Coefficient is equal to 0:
# Implemented by default...
summary( reg4 )

# Let test: H0:  = 0, HA: ln_population_mi neq 5
linearHypothesis( reg4 , "ln_confirmed = 0")
# Can reject null hypothesis:  equal to 0 

######
# Residual analysis.

# Get the predicted y values from the model
df$reg4_y_pred <- reg4$fitted.values
# Calculate the errors of the model
df$reg4_res <- df$ln_death - df$reg4_y_pred 

# Find countries with largest negative errors
df %>% top_n( -5 , reg4_res ) %>% 
  select( country , ln_death , reg4_y_pred , reg4_res )

# Find countries with largest positive errors
df %>% top_n( 5 , reg4_res ) %>% 
  select( country , ln_death , reg4_y_pred , reg4_res )


#################################
## Prediction uncertainty
#

# CI of predicted value/regression line is implemented in ggplot
ggplot( data = df, aes( x = ln_confirmed , y = ln_death ) ) + 
  geom_point( color='blue') +
  geom_smooth( method = lm , color = 'red' , se = T )


# CI of regression line
pred4_CI <- predict( reg4, newdata = df , interval ="confidence" , alpha = 0.05 )
pred4_CI


# Hand made CI for regression line
# 1) Add to datatset:
df <- df %>% mutate( CI_reg4_lower = pred4_CI$fit[,2],
                     CI_reg4_upper = pred4_CI$fit[,3] )
# 2) Plot
ggplot(  ) + 
  geom_point( data = df, aes( x = ln_confirmed , y = ln_death ) , color='blue') +
  geom_line( data = df, aes( x = ln_confirmed , y = reg4_y_pred ) , color = 'red' , size = 1 ) +
  geom_line( data = df, aes( x = ln_confirmed, y = CI_reg4_lower ) , color = 'green' ,
             size = 1 , linetype = "dashed" ) +
  geom_line( data = df, aes( x = ln_confirmed, y = CI_reg4_upper ) , color = 'black' ,
             size = 1 , linetype = "dashed" ) +
  labs(x = "ln( Number of Confirmed Cases )",y = "ln (number of deaths)") 


# As expected, due the high R2 considered in the model, the predicted results were a good fit for the model




