\name{measureQ}
\alias{measureQ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating Construct Reliability, Convergent and Discriminant Validity with Sampling Errors using Bootstrapping
}
\description{
Bootstrap percentile and bias-corrected confidence intervals for: construct reliability, standardized factor loadings, Average Variance Extracted (AVE), comparing AVE with square correlation coefficients, and disattenuated correlation.
}
\usage{
measureQ(model, data.source, b.no=1000, cluster="NULL", CI="PCI", omega="omegaH", R.decimal=3,
         HTMT="FALSE")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
User specified CFA model
}
  \item{data.source}{
A data frame containing the observed variables used in the model.
}
  \item{b.no}{
Number of bootstrap samples, must be between 200 and 5,000. Default is 1,000.
}
  \item{cluster}{
Cluster variable for nested data. Parametric bootstrapping instead of non-parametric bootstrapping will be conducted for nested data.
}
  \item{CI}{
Confidence intervals will be used in the outputs. If "PCI", percentile confidence intervals will be used (default). If "BCCI", biased-corrected confidence intervals will be used.
}
  \item{omega}{
Type of omega used for construct reliability of second-order construct. If "omegaH", omega Hierarchical (default)  will be used. If "omegaT", omega Total will be used.
}
  \item{R.decimal}{
Number of decimal places in reports. Default is 3.
}
  \item{HTMT}{
Heterotrait-monotrait ratio, which is the disattenuated correlation coefficient. Henseler, Ringle, and Sarstedt’s (2015) HTMT equation is incorrect, of which composite is formed by the simple average of standardized item scores instead of the simple average of item scores.
}
}
\details{
1.	The significance tests are based on the type of confidence intervals chosen.
2.	Missing values are handled with full-information maximum likelihood.
3.	No cross-loading or correlated indicator residuals is allowed.
4.  Second-order factors should be defined after all first-order factors in the model.
5.	For second-order factor, omegaT shows the overall reliability of the first-order factors, omegaH shows the reliability of the second-order factor.
6.  Before generating the outputs, bootstrapped samples with improper solutions (non-convergence, negative variance, and correlation greater than 1) are removed.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Cheung, G. W., Cooper-Thomas, H. D., Lau, R.S., & Wang, L. C. (2024).  Reporting reliability, convergent and discriminant validity with structural equation modeling: A review and best-practice recommendations.  Asia Pacific Journal of Management, 41(2), 745-783.  DOI: 10.1007/s10490-023-09871-y
}
\author{
Gordon W. Cheung
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

# ===== Example A (Basic Model)
# Data simulated from Yu, M., Lin, H., Wang, G. G., Liu, Y., & Zheng, X. (2021). Is too much as
# bad as too little? The S-curve relationship between corporate philanthropy and employee
# performance. Asia Pacific Journal of Management. https://doi.org/10.1007/s10490-021-09775-9.

# Set working directory to "C:/research/validity"
# setwd("c:/research/validity")


# Load data file (Example_A.csv)
# Data_A <- read.csv(file = "Example_A.csv")

# Specify the model (Model.A) for standardized coefficients
Model.A <- '
              EB =~ EB1 + EB2 + EB3 +EB4 + EB5 + EB6 + EB7 + EB8
              IV =~ IV1 + IV2 + IV3 + IV4 + IV5
              NVG =~ NVG1 + NVG2 + NVG3
              NVA =~ NVA1 + NVA2 + NVA3 + NVA4
             '

# Call the function "measureQ"
measureQ(Model.A, Data_A, b.no=2000)



# ===== Example B (Model with a Single-Item Factor):
# Data simulated from Zahoor, N., Khan, H., Khan, Z., & Akhtar, P. (2021). Responsible innovation
# in emerging markets’ SMEs: The role of alliance learning and absorptive capacity. Asia Pacific
# Journal of Management. doi:10.1007/s10490-022-09843-8.

# IAPD is a single indicator variable in this example
# For single indicator: factor loading = 1; residuals = (1 - reliability)*observed variance

# Set working directory to "C:/research/validity"
# setwd("c:/research/validity")


# Load data file (Example_B.csv)
# Data_B <- read.csv(file = "Example_B.csv")

# Specify the model (Model.B) for standardized coefficients
Model.B <- '
              MD =~ MD1 + MD2 + MD3 + MD4 + MD5
              TD =~ TD1 + TD2 + TD3 + TD4
              CI =~ CI1 + CI2 + CI3 + CI4
              CCKA =~ CCKA1 + CCKA2 + CCKA3 + CCKA4
              RE =~ RE1 + RE2 + RE3 + RE4
              DP =~ DP1 + DP2 + DP3
              IAPD =~ 1*IAPD1
              IAPD1 ~~ .0097*IAPD1
             '

# Call the function "measureQ"
measureQ(Model.B, Data_B, R.decimal=4)



# ===== Example C (Model with Higher-Order Factors):
# Data simulated from Lythreatis, S., El-Kassar, A., Smart, P., & Ferraris, A. (2022).
# Participative leadership, ethical climate and responsible innovation perceptions: evidence
# from South Korea. Asia-Pacific Journal of Management.  doi: 10.1007/s10490-022-09856-3

# All second-order factors should be defined after all first-order factors in the model.

# Set working directory to "C:/research/validity"
# setwd("c:/research/validity")


# load data file (Example_C.csv)
# Data_C <- read.csv(file = "Example_C.csv")

# Specify the model (Model.C) for standardized coefficients
Model.C <- '
              PL =~ PL1 + PL2 + PL3 + PL4 + PL5 + PL6
              CRNG =~ ETH1 + ETH2 + ETH3 + ETH4 + ETH5 + ETH6
              RULES =~ ETH7 + ETH8 + ETH9
              LC =~ ETH10 + ETH11
              ANTC =~ RI1 + RI2 + RI3
              RFLX =~ RI4 + RI5 + RI6
              INCL =~ RI7 + RI8 + RI9
              RSPN =~ RI10 + RI11 + RI12
              ETH =~ CRNG + RULES + LC
              RI =~ ANTC + RFLX + INCL + RSPN
             '

# Call the function "measureQ"
measureQ(Model.C, Data_C, HTMT="TRUE")



# ===== Example D (Model with High-Order Factor):
# Data simulated from Way, S. A., Tracey, J. B., Fay, C. H., Wright, P. M., Snell, S. A., Chang,
# S., & Gong, Y. (2015). Validation of a multidimensional HR flexibility measure.
# Journal of Management, 41, 1098-1131.

# All second-order factors should be defined after all first-order factors in the model.

# Set working directory to "C:/research/validity"
# setwd("c:/research/validity")


# Load data file (Example_D.csv)
# Data_D <- read.csv(file = "Example_D.csv")

# Specify the model (Model.D) for standardized coefficients
Model.D <- '
            RFHRP =~ RFHRP1+RFHRP2+RFHRP3+RFHRP4+RFHRP5
            RFE =~ RFE1+RFE2+RFE3+RFE4
            CFHRP =~ CFHRP1+CFHRP2+CFHRP3+CFHRP4
            CFCW =~ CFCW1+CFCW2+CFCW3+CFCW4
            CFE =~ CFE1+CFE2+CFE3+CFE4
            HRF =~ RFHRP+RFE+CFHRP+CFCW+CFE
             '

# Call the function "measureQ"
measureQ(Model.D, Data_D)



# ===== Example E (Model with Higher-Order Factor and Nested Data)
# Data simulated from Wayne, S. J., Lemmon, G., Hoobler, J. M., Cheung, G. W., & Wilson, M. S.
# (2017). The ripple effect: A spillover model of the detrimental impact of work-family conflict
# on job success. Journal of Organizational Behavior, 38, 876-894.

# All second-order factors should be defined after all first-order factors in the model.
# The data are nested within the cluster variable = "leader"

# Set working directory to "C:/research/validity"
# setwd("c:/research/validity")


# load data file (Example_E.csv)
# Data_E <- read.csv(file = "Example_E.csv")

# Specify the model (Model.E) for standardized coefficients
Model.E <- '
              X =~ x1 + x2 + x3
              M =~ m1 + m2 + m3
              W =~ w1 + w2 + w3
              YV =~ y1 + y2 + y3
              YD =~ y4 + y5 + y6
              YA =~ y7 + y8 + y9
              sprom =~ sprom1 + sprom2 + sprom3
              Salary =~ 1*SAL
              SAL ~~ 0*SAL
              Perform =~ 1* perfev
              perfev ~~ 0*perfev
              Y =~ YV + YD + YA
             '

# Call the function "measureQ"
measureQ(Model.E, Data_E, cluster="leader")



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
