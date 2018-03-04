library(glmnet)
library(corrplot)
library(MASS)
library(psych) 
library(polycor) 
library(readr)

loan <- read_csv("C:/.../lending_club_loans_Cleaned.csv", 
           col_types = cols(int_rate = col_number(), 
           revol_util = col_number()))
View(loan)

#More Data Cleaning
loan = subset(loan, purpose != 'other')
loan = subset(loan, loan$annual_inc<6000000)

#Shorter names for the boxplots graphs and new variables
loan$purpose[loan$purpose=='credit_card'] <- "crca"
loan$purpose[loan$purpose=='small_business'] <- "smb"
loan$purpose[loan$purpose=='debt_consolidation'] <- "dbtc"
loan$purpose[loan$purpose=='home_improvement'] <- "hmim"
loan$purpose[loan$purpose=='major_purchase'] <- "mapu"
loan$purpose[loan$purpose=='medical'] <- "med"
loan$purpose[loan$purpose=='moving'] <- "mov"
loan$purpose[loan$purpose=='vacation'] <- "vac"
loan$purpose[loan$purpose=='educational'] <- "edu"
loan$purpose[loan$purpose=='renewable_energy'] <- "reen"
loan$log_annual_inc = log(loan$annual_inc/1000)
loan$log_annual_inc = ifelse(log(loan$annual_inc/1000) < 0, 0, log(loan$annual_inc/1000))
loan$fico = (loan$fico_range_low + loan$fico_range_high) / 2
loan$amnt1000 = loan$loan_amnt/1000
loan$sqrt_amnt = sqrt(loan$loan_amnt/1000)
loan$log_fico = log(loan$fico)

#Histograms
hist(loan$loan_amnt)
hist(loan$sqrt_amnt)
hist(loan$annual_inc)
hist(loan$log_annual_inc)
hist(loan$fico)
hist(loan$log_fico)

#Boxplot, Quartiles and Anova for different variables and purpose
BoxPlotAnova <- function(df, predictor_col, cat_col, pred_col_name, cat_col_name)
  { aggregate(df[[predictor_col]], by=list(Category=df[[cat_col]]), na.rm = TRUE, FUN=quantile)
    boxplot(df[[predictor_col]]  ~ df[[cat_col]], main="", xlab=cat_col_name, ylab=pred_col_name)
    anova= lm(df[[predictor_col]] ~ 0+factor(df[[cat_col]]))
    summary(anova)
    int_anova = round(confint(anova),2)
    row.names(int_anova) <- levels(factor(df[[cat_col]]))
    print(int_anova)
    return(int_anova)
}

a1 = BoxPlotAnova(loan, 'amnt1000', 'purpose', "Loan amount (Thousands of $)", "Purpose")
a2 = BoxPlotAnova(loan, 'sqrt_amnt', 'purpose', "Sqrt. Loan amount", "Purpose")
a3 = BoxPlotAnova(loan, 'int_rate', 'purpose', "Interest Rate", "Purpose")
a4 = BoxPlotAnova(loan, 'log_annual_inc', 'purpose', "Log Annual Income", "Purpose")
a5 = BoxPlotAnova(loan, 'annual_inc', 'purpose', "Annual Income", "Purpose")
a6 = BoxPlotAnova(loan, 'fico', 'purpose', "FICO", "Purpose")

#Loan Dataframe
loandf = as.data.frame(loan)
loandf$grade = as.factor(loandf$grade)
loandf$term = as.factor(loandf$term)
loandf$home_ownership = as.factor(loandf$home_ownership)
loandf$purpose = as.factor(loandf$purpose)
loandf$annual_inc = as.integer(loandf$annual_inc)

#MODELLING TO KNOW MAIN PREDICTORS
fit = lm(int_rate ~ sqrt_amnt + term + emp_length + factor(home_ownership) + verification_status + factor(purpose) + log_fico + dti + revol_bal + revol_util + log_annual_inc, data=loan)
summary(fit)

fitp = lm(int_rate ~ factor(purpose), data=loan)
summary(fitp)

fitf = lm(int_rate ~ fico + amnt1000, data=loan)
summary(fitf)

fitficoamnt1 = lm(int_rate ~ log_fico + sqrt_amnt, data=loan)
summary(fitficoamnt1)

fitficoamnt2 = lm(int_rate ~ fico + amnt1000, data=loan)
summary(fitficoamnt2)


#######################
######### MDS #########
#######################

# Normalizing and merging Interest Rate, Loan Amount and FICO score by purpose category
ratea = as.matrix(a3)
rateamed = (ratea[,1]+ratea[,2]) / 2
rateaval = (rateamed - mean(rateamed)) / sd(rateamed)
print(cbind(ratea, rateaval))

ficoa = as.matrix(a6)
ficoamed = (fitoa[,1]+fitoa[,2]) / 2
ficoaval = (fitoamed - mean(fitoamed)) / sd(ficoamed)
ficoStdz = round((ficoamed - mean(ficoamed)) / sd(ficoamed),2)
print(cbind(ficoa, ficoStdz))

amnt1000a = as.matrix(a3)
amnt1000amed = (amnt1000a[,1]+amnt1000a[,2]) / 2
amnt1000aval = (amnt1000amed - mean(amnt1000amed)) / sd(amnt1000amed)
amnt1000stdz = round((amnt1000amed - mean(amnt1000amed)) / sd(amnt1000amed),2)
print(cbind(amnt1000a, amnt1000stdz))

#MDS columns (FICO + Loan Amount)
mdscols = cbind(amnt1000aval, ficoaval)
fplusa = ficoaval+amnt1000aval
feuca = (ficoaval^2 +amnt1000aval^2)^(1/2)
rplusf = ficoaval+rateaval
mdscols = cbind(mdscols, fplusa, feuca, rplusf)
print(mdscols)

#Calculate Distances between categories by Interest Rate, Loan Amount and FICO score
# - Interest Rate
# - FICO score
# - Loan Amount
ratedist = dist(rateaval, method = "minkowski", diag = FALSE, upper = FALSE, p = 2)
print(round(ratedist,2))
ficodist = dist(ficoaval, method = "minkowski", diag = FALSE, upper = FALSE, p = 2)
print(round(ficodist,2))
amnt1000dist = dist(amnt1000aval, method = "minkowski", diag = FALSE, upper = FALSE, p = 2)
print(round(amnt1000dist,2))

# Total distances (composed)
# - Euclidean distance (FICO - Loan Amount)
# - Manhattan distance (FICO - Loan Amount)
# - Euclidean distance (FICO - Int. rate)
# - Manhattan distance (FICO - Int. rate)
distfae = sqrt(ficodist**2 + amnt1000dist**2)
distfam = ficodist + amnt1000dist
distrfe = sqrt(ficodist**2 + ratedist**2)
distrfm = ficodist + ratedist

# Fit MDS and plot
# - Euclidean distance (FICO - Loan Amount)
# - Manhattan distance (FICO - Loan Amount)
# - Euclidean distance (FICO - Int. rate)
# - Manhattan distance (FICO - Int. rate)
fitMDS <- function(distances, numofdim, main_text, Coord1, Coord2)
{ fitmds <- isoMDS(distances, k=numofdim) # k is the number of dim
  fitmds # view results
  x <- fitmds$points[,1]
  y <- fitmds$points[,2]
  plot(x, y, xlab=Coord1, ylab=Coord2, 
     main=main_text, type="n")
  text(x, y, labels = row.names(ficoa), cex=.7)
  return(fitmds)
}

fitmds1 = fitMDS(distfae, 2, "Euclidean distance (FICO - Loan Amount) MDS", "Coordinate 1", "Coordinate 2")
fitmds2 = fitMDS(distfam, 2, "Manhattan distance (FICO - Loan Amount) MDS", "Coordinate 1", "Coordinate 2")
fitmds3 = fitMDS(distrfe, 2, "Euclidean distance (FICO - Int. rate) MDS", "Coordinate 1", "Coordinate 2")
fitmds4 = fitMDS(distrfm, 2, "Manhattan distance (FICO - Int. rate) MDS", "Coordinate 1", "Coordinate 2")

#Finding an explanation for the axis
# - Example: Euclidean distance (FICO - Loan Amount) MDS
# - Correlation helps us to know if there is any relationship between the axis
#   and variables of the dataset.
mdspts1 = as.data.frame(fitmds1$points)
mdspoints1 = cbind(mdspts1, fplusa, feuca, rateaval, fitoaval, amnt1000aval)
names(mdspoints1) <- c("edsx", "edsy", "fico+amnt", 'ficoEucAmnt', 'rate','fico', 'amnt')
corrmds1 = corr.test(mdspoints1,adjust="none")
corrmds1
corrplot(cor(mdspoints1), method="ellipse")


###########################################
# REORDERING AND CHANGING NAMES
x <- -fitmds1$points[,1]
y <- -fitmds1$points[,2]
plot(x, y, xlab="Loan Amount + FICO", ylab="Interest rate", 
     main="Purpose MDS", type="n")
text(x, y, labels = row.names(fitoa), cex=.7)
