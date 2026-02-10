## ======================================================================================== ##
## R Code for Estimating Construct Reliability and Testing Convergent/Discriminant Validity ##
## ======================================================================================== ##



# ==================== Creating Function "measureQ" ==================== #
measureQ <- function(model, data.source, b.no=1000, cluster="NULL", CI="PCI", omega="omegaH", R.decimal=3, HTMT="FALSE") {

options("width"=210)

match.arg(omega, c("omegaH","omegaT"))
match.arg(HTMT, c("FALSE","TRUE"))
match.arg(CI, c("PCI","BCCI"))
if (CI == "BCCI") {
  PCI <- "FALSE"
  BCCI <- "TRUE"
} else {
  PCI <- "TRUE"
  BCCI <- "FALSE"
}

## Test Program Begin

# load data file (Example_E.csv)
#Data_E <- read.csv(file = "Example_E.csv")

# specify the model (Model.E) for standardized coefficients
#Model.E <- '
#              X =~ x1 + x2 + x3
#              M =~ m1 + m2 + m3
#              W =~ w1 + w2 + w3
#              YV =~ y1 + y2 + y3
#              YD =~ y4 + y5 + y6
#              YA =~ y7 + y8 + y9
#              sprom =~ sprom1 + sprom2 + sprom3
#              Extra =~ X+M
#              Salary =~ 1*SAL
#              SAL ~~ 0*SAL
#              Perform =~ 1* perfev
#              perfev ~~ 0*perfev
#              Y =~ YV + YD + YA
#             '

#model <- Model.E
#data.source <- Data_E
#b.no <- 1000
#cluster = "NULL"
#omega <- "omegaH"
#HTMT = "TRUE"
#PCI = "TRUE"
#BCCI = "TRUE"

## Test Program End


# Check for Report decimal places (R.decimal)
R.integer <- R.decimal == round(R.decimal)
if (R.integer == "FALSE") stop("R.decimal must be an integer")
if (R.decimal > 6) stop("More than 6 decimal places in reports is not recommended")
if (R.decimal < 1) stop("Less than 1 decimal place in reports is not recommended")


# Check for bootstrap sample number (b.no)
b.no.integer <- b.no == round(b.no)
if (b.no.integer == "FALSE") stop("Bootstrap sample number must be an integer")
if (b.no > 5000) stop("Bootstrap sample number greater than 5,000 is not recommended")
if (b.no < 500) stop("Bootstrap sample number smaller than 500 is not recommended")

# Estimate model with MLR
if (cluster == "NULL") {
  Model.EST1 <- lavaan::cfa(model, data.source, meanstructure=TRUE, estimator="MLR", missing="fiml")
} else {
  Model.EST1 <- lavaan::cfa(model, data.source, meanstructure=TRUE, estimator="MLR", cluster=cluster, missing="fiml")
}


# Check for cross-loading
DL <- lavaan::inspect(Model.EST1, what="partable")$lambda
if (sum(DL != 0) > nrow(DL))
  stop("There is cross-loading in the model, which has violated discriminant validity")

# Check for correlated residuals
CE <- lavaan::inspect(Model.EST1, what="partable")$theta
if (sum(CE[lower.tri(CE)] != 0) > 0) {
  cat(rep("\n", 3))
  cat("******************************************************", "\n")
  cat("*  WARNING: Some indicator residuals are correlated  *", "\n")
  cat("******************************************************", rep("\n",3))
}

# Search for first-order lv for second-order factors
names.lv.nox <- lavaan::lavNames(model, type = "lv.nox")
no.lv.nox <- length(names.lv.nox)

# extract estimated parameters
SDS <- lavaan::standardizedSolution(Model.EST1) # standardized solution (SDS)
SDS <- SDS[SDS[,"op"] == "=~",]
lv.cor <- lavaan::inspect(Model.EST1, what = "cor.lv") # latent variable correlation
lv.cov <- lavaan::inspect(Model.EST1, what = "cov.lv") # latent variable covariance

# Check for negative standardized factor loading
if (sum(SDS[,4] < 0) > 0) {
  # ===== Print summary statistics and fit indices of unstandardized parameters
  cat(rep("\n", 3))
  cat("\n", "Estimates from Original Sample", rep("\n", 2))
  print(lavaan::summary(Model.EST1, fit.measures=TRUE))
  cat("\n")
  stop("All factor loadings must be positive for estimation of reliability. Please recode your observed variables.")
}


# latent variable and indicator names
names.ov <- lavaan::lavNames(model, type = "ov.ind")

xxx <- c(rep("", no.lv.nox))
for (r in 1: no.lv.nox) {
  xxx[r] <- SDS[which(SDS["rhs"] == names.lv.nox[r] & SDS["op"] == "=~"), "lhs"]
}
so.lv <- xxx[!duplicated(xxx)] # second-order latent variable (so.lv)
names.xxx <-  lavaan::lavNames(Model.EST1, type = "lv.x")
names.lv <- c(names.xxx[!names.xxx %in% so.lv], lavaan::lavNames(Model.EST1, type = "lv.nox"), so.lv)
no.factor <- length(names.lv) # no. of factors

### Check for second-order factor requesting HTMT
##if (length(so.lv) > 0 & HTMT == "TRUE")
##  stop("HTMT is not defined for second-order factor")


# total number of correlation coefficients among latent variables
tot.lv.cor <- no.factor*(no.factor + 1)/2
COR.lv <- matrix(nrow = tot.lv.cor, ncol = 4)

# Replace diagonal elements of latent correlation with standard deviations
diag(lv.cor) <- sqrt(diag(lv.cov))

n <- 0
for (factor.a in 1:no.factor) {
  for (factor.b in factor.a:no.factor) {
    n <- n + 1
    COR.lv[n, 1] <- names.lv[factor.a]
    COR.lv[n, 2] <- "~~"
    COR.lv[n, 3] <- names.lv[factor.b]
    COR.lv[n, 4] <- lv.cor[which(rownames(lv.cor) == names.lv[factor.a]), which(rownames(lv.cor) == names.lv[factor.b])]
  }
}
colnames(COR.lv) <- colnames(SDS)[1:4]

# total number of factor loadings
tot.fl <- sum(SDS[,"op"] == "=~")

# Combining standardized factor loadings with latent variance/correlation
SDS <- rbind(SDS[1:tot.fl,1:4], COR.lv)
SDS[, 4] <- as.numeric(SDS[, 4])
par.est <- as.numeric(SDS[, 4])

# Reorder data.source according to observed variable names
if (cluster == "NULL") {
  data.source <- data.source[, names.ov]
  } else {
  data.source <- data.source[, c(names.ov, cluster)]
}


# name of items (name.items[]) and number of items (no.items[]) of each factor
no.items <- matrix(1:no.factor, nrow = 1)
for (factor.no in 1:no.factor) {
  no.items[factor.no] <- sum(SDS[,"lhs"] == names.lv[factor.no] & SDS[,"op"] == "=~")
}
name.items <- matrix(nrow = no.factor, ncol = length(names.ov))
rownames(name.items) <- names.lv
if (length(so.lv) > 0) {
  for (factor.no in 1:no.factor) {
    if (names.lv[(factor.no)] %in% so.lv) {
      yyy <- SDS[which(SDS["lhs"] == names.lv[factor.no] & SDS["op"] == "=~"), "rhs"]
      for (r in 1: length(yyy)) {
        if (r == 1) {
          zzz <- SDS[which(SDS["lhs"] == yyy[r] & SDS["op"] == "=~"), "rhs"]
        } else {
          zzz <- c(zzz, SDS[which(SDS["lhs"] == yyy[r] & SDS["op"] == "=~"), "rhs"])
        }
      }
      name.items[factor.no, 1:length(zzz)] <- zzz
      no.items[factor.no] <- length(zzz)
    } else {
      name.items[factor.no, 1:no.items[factor.no]] <- SDS[which(SDS["lhs"] == names.lv[factor.no] & SDS["op"] == "=~"), "rhs"]
    }
  }
} else {
  for (factor.no in 1:no.factor) {
    name.items[factor.no, 1:no.items[factor.no]] <- SDS[which(SDS["lhs"] == names.lv[factor.no] & SDS["op"] == "=~"), "rhs"]
  }
}


# ==================== Bootstrapping ==================== #

BCOR.lv <- matrix(nrow = tot.lv.cor, ncol = 1)

if (HTMT == "TRUE"){
  htmtbs <- matrix(nrow = tot.lv.cor, ncol = 1)
  htmt.mat <- matrix(nrow=no.factor, ncol=no.factor)
  colnames(htmt.mat) <- names.lv
  rownames(htmt.mat) <- names.lv
#  begin.item <- matrix(no.factor)  # beginning item number (h)
#  begin.item[1] <- 1
#  for (i in 2:no.factor) {
#    begin.item[i] <- begin.item[i-1] + no.items[i-1]
#  }

  myFUN <- function(x) {
    BDS <- lavaan::standardizedSolution(x)
    BDS <- BDS[BDS[,"op"] == "=~",]
    lv.cor <- lavaan::inspect(x, what = "cor.lv")
    lv.cov <- lavaan::inspect(x, what = "cov.lv")
    diag(lv.cor) <- sqrt(diag(lv.cov))
    m <- lavaan::inspect(x, what = "sampstat")
    m <- m$cov

    n <- 0
    for (factor.a in 1:no.factor) {
      for (factor.b in factor.a:no.factor) {
        n <- n + 1
        BCOR.lv[n] <- lv.cor[which(rownames(lv.cor) == names.lv[factor.a]), which(rownames(lv.cor) == names.lv[factor.b])]
     }
    }
    for (i in 1:no.factor) {
      for (j in i:no.factor) {
        if (no.items[i] == 1) {
          MTi <- SDS[which(SDS["lhs"] == names.lv[i] & SDS["op"] == "=~"), "est.std"]^2*m[name.items[i],name.items[i]]
        } else {
          n <- m[name.items[i,1:no.items[i]], name.items[i,1:no.items[i]]]
          MTi <- mean(n[lower.tri(n)])
        }
        if (no.items[j] == 1) {
          MTj <- SDS[which(SDS["lhs"] == names.lv[j] & SDS["op"] == "=~"), "est.std"]^2*m[name.items[j,1:no.items[j]], name.items[j,1:no.items[j]]]
        } else {
          n <- m[name.items[j,1:no.items[j]], name.items[j,1:no.items[j]]]
          MTj <- mean(n[lower.tri(n)])
        }
        htmt.mat[j,i] <- mean(m[name.items[j,1:no.items[j]], name.items[i,1:no.items[i]]])/(MTi*MTj)^0.5
      }
    }
    for (j in 1:no.factor) {
      htmt.mat[j,j] <- 1
    }
    htmtbs <- htmt.mat[lower.tri(htmt.mat, diag="TRUE")]
    ADS <- c(as.numeric(BDS[,4]), BCOR.lv, htmtbs)
  }

} else {
  myFUN <- function(x) {
    BDS <- lavaan::standardizedSolution(x)
    BDS <- BDS[BDS[,"op"] == "=~",]
    lv.cor <- lavaan::inspect(x, what = "cor.lv")
    lv.cov <- lavaan::inspect(x, what = "cov.lv")
    diag(lv.cor) <- sqrt(diag(lv.cov))
    n <- 0
    for (factor.a in 1:no.factor) {
      for (factor.b in factor.a:no.factor) {
        n <- n + 1
        BCOR.lv[n] <- lv.cor[which(rownames(lv.cor) == names.lv[factor.a]), which(rownames(lv.cor) == names.lv[factor.b])]
      }
    }
    ADS <- c(as.numeric(BDS[,4]), BCOR.lv)
  }
}

# ===== Simplified bootstrapping model
Model.EST2 <- lavaan::cfa(model, data = data.source, missing = "fiml", se = "none", test = "none", check.start = FALSE, check.post = FALSE, check.gradient = FALSE)

# ===== Bootstrapping
if (cluster == "NULL") {
  # Nonparametric Bootstrapping
  bootcoef <- lavaan::bootstrapLavaan(Model.EST2, R = b.no, FUN = myFUN, parallel="snow")
} else {
  # Parametric Bootstrapping when cluster != "NULL"
  bootcoef <- lavaan::bootstrapLavaan(Model.EST2, R = b.no, FUN = myFUN, parallel="snow", type = "parametric")
}


# ===== Remove bootstrap samples with standardized factor loading or correlation larger than 1
for (r in 1: tot.fl) {
  bootcoef <- bootcoef[bootcoef[,r] <= 1,]
}

for (i in 1: factor.no) {
  for (j in 1: i) {
    r <- (j-1)*(no.factor) - sum(1:j-1) + factor.no
    if (i == j) next
    bootcoef <- bootcoef[which(abs(bootcoef[,tot.fl+r]) <= 1),]
    bootcoef <- bootcoef[abs(bootcoef[,tot.fl+r]) <= 1,]
  }
}

r.b.no <- b.no
b.no <- nrow(bootcoef)

# ===== Print number of bootstrap samples
cat("\n", "   Number of requested bootstrap sample = ", r.b.no)
cat("\n", "   Number of completed bootstrap sample = ", b.no, rep("\n",2))

# ==================================================================== #


# =========================== Sample Estimates ======================================#


# ===== Construct Reliability & Average Variance Extracted
SCR <- 1:no.factor
SAVE <- 1: no.factor
names(SCR) <- names.lv
names(SAVE) <- names.lv
for (factor.no in 1:no.factor) {
  if (names.lv[(factor.no)] %in% so.lv) {
    no.items.so.lv <- sum(SDS[,"lhs"] == names.lv[factor.no] & SDS[,"op"] == "=~")
    yyy <- SDS[which(SDS["lhs"] == names.lv[factor.no] & SDS["op"] == "=~"), "rhs"]
    G <- c(1:no.items.so.lv)
    E <- c(1:no.items.so.lv)
    L <- c(1:no.items.so.lv)
    A <- c(1:no.items.so.lv)
    NI <- c(1:no.items.so.lv) # number of items in first order factor
    for (r in 1: no.items.so.lv) {
      rows <- which(SDS["lhs"] == yyy[r] & SDS["op"] == "=~")
      L[r] <- sum(SDS[rows,"est.std"])
      E[r] <- sum(1 - SDS[rows,"est.std"]^2)
      G[r] <- SDS[which(SDS["rhs"] == yyy[r] & SDS["op"] == "=~"), "est.std"]
      A[r] <- sum(SDS[rows,"est.std"]^2)*G[r]^2
      NI[r] <- NROW(rows)
    }
    omegaH <- sum(L*G)^2/(sum(L*G)^2+sum((1-G^2)*L^2)+sum(E))
    omegaT <- (sum(L*G)^2+sum((1-G^2)*L^2))/(sum(L*G)^2+sum((1-G^2)*L^2)+sum(E))
    if (omega == "omegaH") {
       SCR[factor.no] <- omegaH
    } else {
       SCR[factor.no] <- omegaT
    }
    SAVE[factor.no] <- sum(A)/sum(NI)
  } else {
    SCR[factor.no] <- sum(as.numeric(SDS[which(SDS[,"lhs"] == names.lv[factor.no] & SDS[,"op"] == "=~"), "est.std"]))^2 /
                   (sum(as.numeric(SDS[which(SDS[,"lhs"] == names.lv[factor.no] & SDS[,"op"] == "=~"), "est.std"]))^2 +
                    no.items[factor.no] -
                    sum(as.numeric(SDS[which(SDS[,"lhs"] == names.lv[factor.no] & SDS[,"op"] == "=~"), "est.std"])^2))
    SAVE[factor.no] <- sum(as.numeric(SDS[which(SDS[,"lhs"] == names.lv[factor.no] & SDS[,"op"] == "=~"), "est.std"])^2)/no.items[factor.no]
  }
}

# ============= Generate simple average score and Cronbach's alpha of summated scales (SCA) ================ #

SCA <- matrix(0, no.factor)

for (factor.no in 1:no.factor) {
  if (no.items[factor.no] == 1) {
    cols <- which(SDS["lhs"] == names.lv[factor.no] & SDS["op"] == "=~")
    data.source[, names.lv[factor.no]] <- data.source[, SDS[cols,"rhs"]]
    SCA[factor.no] <- SCR[factor.no]
  } else {
    cols <- which(SDS["lhs"] == names.lv[factor.no] & SDS["op"] == "=~")
    data.source[, names.lv[factor.no]] <- rowMeans(data.source[, SDS[cols,"rhs"]], na.rm = TRUE)
    var.item.total <- sum(diag(var(data.source[, cols], na.rm = TRUE)))/(no.items[factor.no]^2)
    SCA[factor.no] <- no.items[factor.no]/(no.items[factor.no] - 1) * (1 - var.item.total/var(data.source[,names.lv[factor.no]], na.rm = TRUE))
  }
}

if (length(so.lv) > 0) {
  for (factor.no in 1:no.factor) {
    if (no.items[factor.no] == 1) {
      cols <- which(SDS["lhs"] == names.lv[factor.no] & SDS["op"] == "=~")
      data.source[, names.lv[factor.no]] <- data.source[, SDS[cols,"rhs"]]
      SCA[factor.no] <- SCR[factor.no]
    } else {
      if (names.lv[(factor.no)] %in% so.lv) {
        no.items.so.lv <- sum(SDS[,"lhs"] == names.lv[factor.no] & SDS[,"op"] == "=~")
        NI <- c(1:no.items.so.lv) # number of items in first order factor
        yyy <- SDS[which(SDS["lhs"] == names.lv[factor.no] & SDS["op"] == "=~"), "rhs"]
        for (r in 1: no.items.so.lv) {
          if (r == 1) {
            zzz <- which(SDS["lhs"] == yyy[r] & SDS["op"] == "=~")
          } else {
            zzz <- c(zzz, which(SDS["lhs"] == yyy[r] & SDS["op"] == "=~"))
          }
        }
        data.source[, names.lv[factor.no]] <- rowMeans(data.source[, SDS[zzz,"rhs"]], na.rm = TRUE)
        var.item.total <- sum(diag(var(data.source[, zzz], na.rm = TRUE)))/(length(zzz)^2)
        SCA[factor.no] <- length(zzz)/(length(zzz) - 1) * (1 - var.item.total/var(data.source[,names.lv[factor.no]], na.rm = TRUE))
      } else {
        cols <- which(SDS["lhs"] == names.lv[factor.no] & SDS["op"] == "=~")
        data.source[, names.lv[factor.no]] <- rowMeans(data.source[, SDS[cols,"rhs"]], na.rm = TRUE)
        var.item.total <- sum(diag(var(data.source[, cols], na.rm = TRUE)))/(no.items[factor.no]^2)
        SCA[factor.no] <- no.items[factor.no]/(no.items[factor.no] - 1) * (1 - var.item.total/var(data.source[,names.lv[factor.no]], na.rm = TRUE))
      }
    }
  }
}


# ======================================================================= #



# ===== Discriminant Validity - AVE compare with rXY-square
SDV <- matrix(nrow = tot.lv.cor, ncol = 4)
for (r in 1:tot.lv.cor) {
  factor.a <- COR.lv[r,1]
  factor.b <- COR.lv[r,3]
  SDV[r, 1] <- factor.a
  SDV[r, 2] <- factor.b
  SDV[r, 3] <- SAVE[factor.a] - as.numeric(COR.lv[r, 4])^2
  SDV[r, 4] <- SAVE[factor.b] - as.numeric(COR.lv[r, 4])^2
}


# ==================================================================================== #



# =========================== Bootstrap Estimates ==================================== #

# ===== Average Variance Extracted (BAVE) and Construct Reliability (BCR)
BAVE <- matrix(0, nrow = b.no, ncol = no.factor)
colnames(BAVE) <- names.lv
BCR <- matrix(0, nrow = b.no, ncol = no.factor)
colnames(BCR) <- names.lv
omegaH <- matrix(0, nrow = b.no, ncol = 1)
omegaT <- matrix(0, nrow = b.no, ncol = 1)

for (boot.no in 1:b.no) {
  for (factor.no in 1:no.factor) {
    if (names.lv[(factor.no)] %in% so.lv) {
      no.items.so.lv <- sum(SDS[,"lhs"] == names.lv[factor.no] & SDS[,"op"] == "=~")
      yyy <- SDS[which(SDS["lhs"] == names.lv[factor.no] & SDS["op"] == "=~"), "rhs"]
      BG <- matrix(0, nrow = b.no, ncol = no.items.so.lv)
      BE <- matrix(0, nrow = b.no, ncol = no.items.so.lv)
      BL <- matrix(0, nrow = b.no, ncol = no.items.so.lv)
      BA <- matrix(0, nrow = b.no, ncol = no.items.so.lv)  ## Average variance extracted for first-ordr factor r
      BNI <- matrix(0, nrow = b.no, ncol = no.items.so.lv) # number of items in first order factor
      for (r in 1: no.items.so.lv) {
        cols <- which(SDS["lhs"] == yyy[r] & SDS["op"] == "=~")
        BL[boot.no, r] <- sum(bootcoef[boot.no, cols])
        BE[boot.no, r] <- sum(1 - bootcoef[boot.no, cols]^2)
        BG[boot.no, r] <- bootcoef[boot.no, which(SDS["rhs"] == yyy[r] & SDS["op"] == "=~")]
        BA[boot.no, r] <- sum(bootcoef[boot.no, cols]^2)*BG[boot.no, r]^2
        BNI[boot.no, r] <- NROW(cols)
      }
      omegaH[boot.no] <- sum(BL[boot.no,]*BG[boot.no,])^2/(sum(BL[boot.no,]*BG[boot.no,])^2+sum((1-BG[boot.no,]^2)*BL[boot.no,]^2)+sum(BE[boot.no,]))
      omegaT[boot.no] <- (sum(BL[boot.no,]*BG[boot.no,])^2+sum((1-BG[boot.no]^2)*BL[boot.no,]^2))/(sum(BL[boot.no,]*BG[boot.no,])^2+sum((1-BG[boot.no,]^2)*BL[boot.no,]^2)+sum(BE[boot.no,]))
      if (omega == "omegaH") {
        BCR[boot.no, factor.no] <- omegaH[boot.no]
      } else {
        BCR[boot.no, factor.no] <- omegaT[boot.no]
      }
      BAVE[boot.no, factor.no] <- sum(BA[boot.no,])/sum(BNI[boot.no,])
    } else {
    BAVE[boot.no, factor.no] <- sum(bootcoef[boot.no, which(SDS[,"lhs"] == names.lv[factor.no] & SDS[,"op"] == "=~")]^2)/no.items[factor.no]
    BCR[boot.no, factor.no] <- sum(bootcoef[boot.no, which(SDS[,"lhs"] == names.lv[factor.no] & SDS[,"op"] == "=~")])^2 /
                   (sum(bootcoef[boot.no, which(SDS[,"lhs"] == names.lv[factor.no] & SDS[,"op"] == "=~")])^2 +
                    no.items[factor.no] -
                    sum(bootcoef[boot.no, which(SDS[,"lhs"] == names.lv[factor.no] & SDS[,"op"] == "=~")]^2))
    }
  }
}

# ===== Discriminant Validity - AVE compare with rXY-square
BDV <- array(0, dim = c(b.no, tot.lv.cor, 2))
for (boot.no in 1:b.no) {
  for (r in 1:tot.lv.cor) {
    factor.a <- COR.lv[r, 1]
    factor.b <- COR.lv[r, 3]
    BDV[boot.no, r, 1] <- BAVE[boot.no, factor.a] - bootcoef[boot.no, which(SDS[,"lhs"] == factor.a & SDS[,"rhs"] == factor.b & SDS[,"op"] == "~~")]^2
    BDV[boot.no, r, 2] <- BAVE[boot.no, factor.b] - bootcoef[boot.no, which(SDS[,"lhs"] == factor.a & SDS[,"rhs"] == factor.b & SDS[,"op"] == "~~")]^2
  }
}

# ==================================================================================== #


# ========================== Construct Reliability      ============================== #

# ===== Percentile Probability (CR.pp)
CR.pp <- matrix(0, no.factor)
for (r in 1:no.factor) {
  if (quantile(BCR[, r],probs=0.5)>0) {
      CR.pp[r] = 2*(sum(BCR[, r]<0)/b.no)
  } else {
    CR.pp[r] = 2*(sum(BCR[, r]>0)/b.no)
  }
}

# ===== Percentile Confidence Intervals of Construct Reliability (CR.PCI)
CR.PCI <- matrix(0, nrow = no.factor, ncol = 9)
colnames(CR.PCI) <- c(" Factor","    0.5%","    2.5%","      5%"," Estimate","     95%","    97.5%","    99.5%"," p-value")

for (r in 1:no.factor) {
  CR.PCI[r, 1] <- names.lv[r]
  CR.PCI[r, 2] <- format(round(quantile(BCR[, r], c(0.005)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  CR.PCI[r, 3] <- format(round(quantile(BCR[, r], c(0.025)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  CR.PCI[r, 4] <- format(round(quantile(BCR[, r], c(0.050)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  CR.PCI[r, 5] <- format(round(SCR[r], digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  CR.PCI[r, 6] <- format(round(quantile(BCR[, r], c(0.950)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  CR.PCI[r, 7] <- format(round(quantile(BCR[, r], c(0.975)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  CR.PCI[r, 8] <- format(round(quantile(BCR[, r], c(0.995)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  CR.PCI[r, 9] <- format(round(CR.pp[r], digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
}

# ===== Print percentile confidence intervals
if (PCI == "TRUE") {
  cat("\n", "   Percentile Confidence Intervals for Construct Reliability", rep("\n", 2))
  rownames(CR.PCI) <- rep("    ", nrow(CR.PCI))
  print(CR.PCI, quote=FALSE, right=TRUE)
  cat("\n")
}

# ===== Bias-Corrected Factor
CR.z <- matrix(0, no.factor)
for (r in 1:no.factor) {
  CR.z[r] <- qnorm(sum(BCR[, r]<SCR[r])/b.no)
}

# ===== Bias-Corrected Probability (CR.bcp)
CR.bcp <- matrix(0, no.factor)  # Bias-Corrected Probability #

for (r in 1:no.factor) {
  if ((SCR[r]>0 & min(BCR[, r])>0) | (SCR[r]<0 & max(BCR[, r])<0)) {
      CR.bcp[r] = 0
    } else if (qnorm(sum(BCR[, r]>0)/b.no)+2*CR.z[r]<0) {
      CR.bcp[r] = 2*pnorm((qnorm(sum(BCR[, r]>0)/b.no)+2*CR.z[r]))
    } else {
      CR.bcp[r] = 2*pnorm(-1*(qnorm(sum(BCR[, r]>0)/b.no)+2*CR.z[r]))
  }
}

# ===== Bias-Corrected Confidence Intervals for Construct Reliability
CR.BCCI <- matrix(0, nrow = no.factor, ncol=9)
colnames(CR.BCCI) <- c(" Factor","    0.5%","    2.5%","      5%"," Estimate","     95%","    97.5%","    99.5%"," p-value")

for (r in 1:no.factor) {
  CR.BCCI[r, 1] <- names.lv[r]
  CR.BCCI[r, 2] <- format(round(quantile(BCR[, r],probs = pnorm(2*CR.z[r]+qnorm(0.005))), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  CR.BCCI[r, 3] <- format(round(quantile(BCR[, r],probs = pnorm(2*CR.z[r]+qnorm(0.025))), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  CR.BCCI[r, 4] <- format(round(quantile(BCR[, r],probs = pnorm(2*CR.z[r]+qnorm(0.050))), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  CR.BCCI[r, 5] <- format(round(SCR[r], digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  CR.BCCI[r, 6] <- format(round(quantile(BCR[, r],probs = pnorm(2*CR.z[r]+qnorm(0.950))), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  CR.BCCI[r, 7] <- format(round(quantile(BCR[, r],probs = pnorm(2*CR.z[r]+qnorm(0.975))), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  CR.BCCI[r, 8] <- format(round(quantile(BCR[, r],probs = pnorm(2*CR.z[r]+qnorm(0.995))), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  CR.BCCI[r, 9] <- format(round(CR.bcp[r], digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
}

# ===== Print bias-corrected confidence intervals
if (BCCI == "TRUE") {
  cat("\n", "   Bias-Corrected Confidence Intervals for Construct Reliability", rep("\n", 2))
  rownames(CR.BCCI) <- rep("    ", nrow(CR.BCCI))
  print(CR.BCCI, quote=FALSE, right=TRUE)
  cat("\n")
}

# ==================================================================================== #


# ========================== Average Variance Extracted ============================== #

# ===== Calculate Percentile Probability (AVE.pp)
AVE.pp <- matrix(0, no.factor)
for (r in 1:no.factor) {
  if (quantile(BAVE[, r],probs=0.5)>0) {
      AVE.pp[r] = 2*(sum(BAVE[, r]<0)/b.no)
  } else {
    AVE.pp[r] = 2*(sum(BAVE[, r]>0)/b.no)
  }
}

# ===== Percentile Confidence Intervals of Average Variance Extracted (AVE.PCI)
AVE.PCI <- matrix(0, nrow = no.factor, ncol = 9)
colnames(AVE.PCI) <- c(" Factor","    0.5%","    2.5%","      5%"," Estimate","     95%","   97.5%","   99.5%","  p-value")

for (r in 1:no.factor) {
  AVE.PCI[r, 1] <- names.lv[r]
  AVE.PCI[r, 2] <- format(round(quantile(BAVE[, r], c(0.005)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  AVE.PCI[r, 3] <- format(round(quantile(BAVE[, r], c(0.025)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  AVE.PCI[r, 4] <- format(round(quantile(BAVE[, r], c(0.050)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  AVE.PCI[r, 5] <- format(round(SAVE[r], digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  AVE.PCI[r, 6] <- format(round(quantile(BAVE[, r], c(0.950)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  AVE.PCI[r, 7] <- format(round(quantile(BAVE[, r], c(0.975)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  AVE.PCI[r, 8] <- format(round(quantile(BAVE[, r], c(0.995)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  AVE.PCI[r, 9] <- format(round(AVE.pp[r], digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
}

# ===== Print percentile confidence intervals
if (PCI == "TRUE") {
  cat("\n", "   Percentile Confidence Intervals for Average Variance Extracted", rep("\n", 2))
  rownames(AVE.PCI) <- rep("    ", nrow(AVE.PCI))
  print(AVE.PCI, quote=FALSE, right = TRUE)
  cat("\n")
}

# ===== Bias-Correction Factor for AVE (AVE.z)
AVE.z <- matrix(0, no.factor)  # Bias-Corrected Factor #
for (r in 1:no.factor) {
  AVE.z[r] <- qnorm(sum(BAVE[, r]<SAVE[r])/b.no)
}

# ===== Bias-Corrected Probability (AVE.bcp)
AVE.bcp <- matrix(0, no.factor)

for (r in 1:no.factor) {
  if ((SAVE[r]>0 & min(BAVE[, r])>0) | (SAVE[r]<0 & max(BAVE[, r])<0)) {
      AVE.bcp[r] = 0
    } else if (qnorm(sum(BAVE[, r]>0)/b.no)+2*AVE.z[r]<0) {
      AVE.bcp[r] = 2*pnorm((qnorm(sum(BAVE[, r]>0)/b.no)+2*AVE.z[r]))
    } else {
      AVE.bcp[r] = 2*pnorm(-1*(qnorm(sum(BAVE[, r]>0)/b.no)+2*AVE.z[r]))
  }
}

# ===== Bias-Corrected confidence intervals (AVE.BCCI)
AVE.BCCI <- matrix(0, nrow = no.factor, ncol=9)
colnames(AVE.BCCI) <- c(" Factor","    0.5%","    2.5%","      5%"," Estimate","     95%","   97.5%","   99.5%","  p-value")

for (r in 1:no.factor) {
  AVE.BCCI[r, 1] <- names.lv[r]
  AVE.BCCI[r, 2] <- format(round(quantile(BAVE[, r],probs = pnorm(2*AVE.z[r]+qnorm(0.005))), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  AVE.BCCI[r, 3] <- format(round(quantile(BAVE[, r],probs = pnorm(2*AVE.z[r]+qnorm(0.025))), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  AVE.BCCI[r, 4] <- format(round(quantile(BAVE[, r],probs = pnorm(2*AVE.z[r]+qnorm(0.050))), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  AVE.BCCI[r, 5] <- format(round(SAVE[r], digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  AVE.BCCI[r, 6] <- format(round(quantile(BAVE[, r],probs = pnorm(2*AVE.z[r]+qnorm(0.950))), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  AVE.BCCI[r, 7] <- format(round(quantile(BAVE[, r],probs = pnorm(2*AVE.z[r]+qnorm(0.975))), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  AVE.BCCI[r, 8] <- format(round(quantile(BAVE[, r],probs = pnorm(2*AVE.z[r]+qnorm(0.995))), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  AVE.BCCI[r, 9] <- format(round(AVE.bcp[r], digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
}

# ===== Print bias-corrected confidence intervals
if (BCCI == "TRUE") {
  cat("\n", "   Bias-Corrected Confidence Intervals for Average Variance Extracted", rep("\n", 2))
  rownames(AVE.BCCI) <- rep("    ", nrow(AVE.BCCI))
  print(AVE.BCCI, quote=FALSE, right=TRUE)
  cat("\n")
}

# ==================================================================================== #


# ========================= Standardized Factor Loadings ============================= #

# ===== Percentile Probability for Standardized Factor Loadings (SFL.pp)
SFL.pp <- matrix(0, tot.fl)
for (r in 1:tot.fl) {
  if (quantile(bootcoef[, r],probs=0.5)>0) {
      SFL.pp[r] = 2*(sum(bootcoef[, r]<0)/b.no)
  } else {
    SFL.pp[r] = 2*(sum(bootcoef[, r]>0)/b.no)
  }
}

# ===== Percentile Confidence Intervals of Standardized Factor Loadings (SFL.PCI)
SFL.PCI <- matrix(0, nrow = tot.fl, ncol = 9)
colnames(SFL.PCI) <- c(" Factor Loading","    0.5%","    2.5%","      5%","  Estimate","     95%","   97.5%","   99.5%","  p-value")

for (r in 1:tot.fl) {
  SFL.PCI[r, 1] <- paste(SDS[r,1], SDS[r,2], SDS[r,3], sep = " ")
  SFL.PCI[r, 2] <- format(round(quantile(bootcoef[, r], c(0.005)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  SFL.PCI[r, 3] <- format(round(quantile(bootcoef[, r], c(0.025)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  SFL.PCI[r, 4] <- format(round(quantile(bootcoef[, r], c(0.050)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  SFL.PCI[r, 5] <- format(round(par.est[r], digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  SFL.PCI[r, 6] <- format(round(quantile(bootcoef[, r], c(0.950)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  SFL.PCI[r, 7] <- format(round(quantile(bootcoef[, r], c(0.975)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  SFL.PCI[r, 8] <- format(round(quantile(bootcoef[, r], c(0.995)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  SFL.PCI[r, 9] <- format(round(SFL.pp[r], digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
}

# ===== Print percentile confidence intervals
if (PCI == "TRUE") {
  cat("\n", "   Percentile Confidence Intervals for Standardized Factor Loadings", rep("\n",2))
  rownames(SFL.PCI) <- rep("    ", nrow(SFL.PCI))
  print(SFL.PCI, quote=FALSE, right=TRUE)
  cat("\n")
}

# ===== Bias-Correction Factor for Standardized Factor Loadings (SFL.z)
SFL.z <- matrix(0, tot.fl)  # Bias-Corrected Factor #
for (r in 1:tot.fl) {
  SFL.z[r] <- qnorm(sum(bootcoef[, r]<par.est[r])/b.no)
}

# ===== Bias-Corrected Probability for Standardized Factor Loadings (SFL.bcp)
SFL.bcp <- matrix(0, tot.fl)  # Bias-Corrected Probability #
for (r in 1:tot.fl) {
  if ((par.est[r]>0 & min(bootcoef[, r])>0) | (par.est[r]<0 & max(bootcoef[, r])<0)) {
      SFL.bcp[r] = 0
    } else if (qnorm(sum(bootcoef[, r]>0)/b.no)+2*SFL.z[r]<0) {
      SFL.bcp[r] = 2*pnorm((qnorm(sum(bootcoef[, r]>0)/b.no)+2*SFL.z[r]))
    } else {
      SFL.bcp[r] = 2*pnorm(-1*(qnorm(sum(bootcoef[, r]>0)/b.no)+2*SFL.z[r]))
  }
}

# ===== Bias-Corrected Confidence Intervlas for Standardized Factor Loadings (SFL.BCI)
SFL.BCCI <- matrix(0, nrow = tot.fl, ncol=9)
colnames(SFL.BCCI) <- c(" Factor Loading","    0.5%","    2.5%","      5%","  Estimate","     95%","   97.5%","   99.5%","  p-value")

for (r in 1:tot.fl) {
  SFL.BCCI[r, 1] <- paste(SDS[r,1], SDS[r,2], SDS[r,3], sep = " ")
  SFL.BCCI[r, 2] <- format(round(quantile(bootcoef[, r],probs = pnorm(2*SFL.z[r]+qnorm(0.005))), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  SFL.BCCI[r, 3] <- format(round(quantile(bootcoef[, r],probs = pnorm(2*SFL.z[r]+qnorm(0.025))), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  SFL.BCCI[r, 4] <- format(round(quantile(bootcoef[, r],probs = pnorm(2*SFL.z[r]+qnorm(0.050))), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  SFL.BCCI[r, 5] <- format(round(par.est[r], digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  SFL.BCCI[r, 6] <- format(round(quantile(bootcoef[, r],probs = pnorm(2*SFL.z[r]+qnorm(0.950))), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  SFL.BCCI[r, 7] <- format(round(quantile(bootcoef[, r],probs = pnorm(2*SFL.z[r]+qnorm(0.975))), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  SFL.BCCI[r, 8] <- format(round(quantile(bootcoef[, r],probs = pnorm(2*SFL.z[r]+qnorm(0.995))), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  SFL.BCCI[r, 9] <- format(round(SFL.bcp[r], digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
}

# ===== Print bias-corrected confidence intervals
if (BCCI == "TRUE") {
  cat("\n", "   Bias-Corrected Confidence Intervals for Standardized Factor Loadings", rep("\n", 2))
  rownames(SFL.BCCI) <- rep("    ", nrow(SFL.BCCI))
  print(SFL.BCCI, quote=FALSE, right=TRUE)
  cat("\n")
}

# ==================================================================================== #


# ============================ Correlation Coefficient =============================== #

# ===== Percentile Probability for Latent Correlation (COR.pp)
COR.pp <- matrix(0, tot.lv.cor)  # Percentile Probability #
for (r in 1:tot.lv.cor) {
  if (quantile(bootcoef[, tot.fl+r],probs=0.5)>0) {
      COR.pp[r] = 2*(sum(bootcoef[, tot.fl+r]<0)/b.no)
  } else {
    COR.pp[r] = 2*(sum(bootcoef[, tot.fl+r]>0)/b.no)
  }
}

# ===== Percentile Confidence Intervals of Correlation Coefficients (COR.PCI)
COR.PCI <- matrix(0, nrow = (tot.lv.cor - no.factor - length(xxx)), ncol = 9)
colnames(COR.PCI) <- c(" Correlation","     0.5%","     2.5%","       5%","   Estimate","      95%","    97.5%","    99.5%","  p-value")

r <- 0
for (j in 1:tot.lv.cor) {
  if (SDS[tot.fl+j, 1] == SDS[tot.fl+j, 3]) {
    next
  }
  if (SDS[tot.fl+j, 1] %in% so.lv) {
    yyy <- SDS[which(SDS["lhs"] == SDS[tot.fl+j, 1] & SDS["op"] == "=~"), "rhs"]
    if (SDS[tot.fl+j, 3] %in% yyy) {
      next
    }
  }
  if (SDS[tot.fl+j, 3] %in% so.lv) {
    yyy <- SDS[which(SDS["lhs"] == SDS[tot.fl+j, 3] & SDS["op"] == "=~"), "rhs"]
    if (SDS[tot.fl+j, 1] %in% yyy) {
      next
    }
  }
  r <- r + 1
  COR.PCI[r, 1] <- paste(SDS[tot.fl+j,1], SDS[tot.fl+j,2], SDS[tot.fl+j,3], sep = " ")
  COR.PCI[r, 2] <- format(round(quantile(bootcoef[, tot.fl+j], c(0.005)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  COR.PCI[r, 3] <- format(round(quantile(bootcoef[, tot.fl+j], c(0.025)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  COR.PCI[r, 4] <- format(round(quantile(bootcoef[, tot.fl+j], c(0.050)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  COR.PCI[r, 5] <- format(round(SDS[tot.fl+j,4], digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  COR.PCI[r, 6] <- format(round(quantile(bootcoef[, tot.fl+j], c(0.950)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  COR.PCI[r, 7] <- format(round(quantile(bootcoef[, tot.fl+j], c(0.975)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  COR.PCI[r, 8] <- format(round(quantile(bootcoef[, tot.fl+j], c(0.995)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  COR.PCI[r, 9] <- format(round(COR.pp[j], digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
}

# ===== Print percentile confidence intervals
if (PCI == "TRUE") {
  cat("\n", "   Percentile Confidence Intervals for Correlation Coefficients", rep("\n", 2))
  rownames(COR.PCI) <- rep("    ", nrow(COR.PCI))
  print(COR.PCI, quote=FALSE, right=TRUE)
  cat("\n")
}

# ===== Bias-Correction Factor for Correlation (COR.z)
COR.z <- matrix(0, tot.lv.cor)  # Bias-Corrected Factor #
for (r in 1:tot.lv.cor) {
  COR.z[r] <- qnorm(sum(bootcoef[, tot.fl+r]<SDS[tot.fl+r,4])/b.no)
}

# ===== Bias-Corrected Probability for Correlation
COR.bcp <- matrix(0, tot.lv.cor)
for (r in 1:tot.lv.cor) {
  if ((par.est[tot.fl+r]>0 & min(bootcoef[, tot.fl+r])>0) | (par.est[tot.fl+r]<0 & max(bootcoef[, tot.fl+r])<0)) {
      COR.bcp[r] = 0
    } else if (qnorm(sum(bootcoef[, tot.fl+r]>0)/b.no)+2*COR.z[r]<0) {
      COR.bcp[r] = 2*pnorm((qnorm(sum(bootcoef[, tot.fl+r]>0)/b.no)+2*COR.z[r]))
    } else {
      COR.bcp[r] = 2*pnorm(-1*(qnorm(sum(bootcoef[, tot.fl+r]>0)/b.no)+2*COR.z[r]))
  }
}

# ===== Bias-Corrected Confidence Intervals for Correlation (COR.BCCI)
COR.BCCI <- matrix(0, nrow = tot.lv.cor - no.factor - length(xxx), ncol=9)
colnames(COR.BCCI) <- c(" Correlation","     0.5%","     2.5%","       5%","   Estimate","      95%","    97.5%","    99.5%","  p-value")

r <- 0
for (j in 1:tot.lv.cor) {
  if (SDS[tot.fl+j, 1] == SDS[tot.fl+j, 3]) {
    next
  }
  if (SDS[tot.fl+j, 1] %in% so.lv) {
    yyy <- SDS[which(SDS["lhs"] == SDS[tot.fl+j, 1] & SDS["op"] == "=~"), "rhs"]
    if (SDS[tot.fl+j, 3] %in% yyy) {
      next
    }
  }
  if (SDS[tot.fl+j, 3] %in% so.lv) {
    yyy <- SDS[which(SDS["lhs"] == SDS[tot.fl+j, 3] & SDS["op"] == "=~"), "rhs"]
    if (SDS[tot.fl+j, 1] %in% yyy) {
      next
    }
  }
  r <- r + 1
  COR.BCCI[r, 1] <- paste(SDS[tot.fl+j,1], SDS[tot.fl+j,2], SDS[tot.fl+j,3], sep = " ")
  COR.BCCI[r, 2] <- format(round(quantile(bootcoef[, tot.fl+j],probs = pnorm(2*COR.z[j]+qnorm(0.005))), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  COR.BCCI[r, 3] <- format(round(quantile(bootcoef[, tot.fl+j],probs = pnorm(2*COR.z[j]+qnorm(0.025))), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  COR.BCCI[r, 4] <- format(round(quantile(bootcoef[, tot.fl+j],probs = pnorm(2*COR.z[j]+qnorm(0.050))), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  COR.BCCI[r, 5] <- format(round(SDS[tot.fl+j,4], digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  COR.BCCI[r, 6] <- format(round(quantile(bootcoef[, tot.fl+j],probs = pnorm(2*COR.z[j]+qnorm(0.950))), digits = R.decimal), nsmall = R.decimal, scientific = FALSE)
  COR.BCCI[r, 7] <- format(round(quantile(bootcoef[, tot.fl+j],probs = pnorm(2*COR.z[j]+qnorm(0.975))), digits = R.decimal), nsmall = R.decimal, scientific = FALSE)
  COR.BCCI[r, 8] <- format(round(quantile(bootcoef[, tot.fl+j],probs = pnorm(2*COR.z[j]+qnorm(0.995))), digits = R.decimal), nsmall = R.decimal, scientific = FALSE)
  COR.BCCI[r, 9] <- format(round(COR.bcp[j], digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
}

# ===== Print bias-corrected confidence intervals
if (BCCI == "TRUE") {
  cat("\n", "   Bias-Corrected Confidence Intervals for Correlation Coefficients", rep("\n", 2))
  cat("\n")
  rownames(COR.BCCI) <- rep("    ", nrow(COR.BCCI))
  print(COR.BCCI, quote=FALSE, right=TRUE)
  cat("\n")
}

# ==================================================================================== #


# ======================== Dircriminant Validity ===================================== #

# ===== Percentile Probability for Discriminant Validity
DV.pp <- matrix(0, tot.lv.cor, 2)

for (r in 1:tot.lv.cor) {
  if (quantile(BDV[, r, 1],probs=0.5)>0) {
    DV.pp[r, 1] = 2*(sum(BDV[, r, 1]<0)/b.no)
  } else {
    DV.pp[r, 1] = 2*(sum(BDV[, r, 1]>0)/b.no)
  }
  if (quantile(BDV[, r, 2],probs=0.5)>0) {
    DV.pp[r, 2] = 2*(sum(BDV[, r, 2]<0)/b.no)
  } else {
    DV.pp[r, 2] = 2*(sum(BDV[, r, 2]>0)/b.no)
  }
}

# ===== Percentile Confidence Intervals for Discriminant Validity (DV.PCI)
tot.corX2 <- (tot.lv.cor - no.factor - length(xxx))*2
DV.PCI <- matrix(0, nrow = tot.corX2, ncol = 10)
colnames(DV.PCI) <- c(" Correlation","      Factor","     0.5%","     2.5%","       5%","   Estimate","      95%","    97.5%","    99.5%","  p-value")

k <- 0
for (j in 1:tot.lv.cor) {
  if (SDS[tot.fl+j, 1] == SDS[tot.fl+j, 3]) {
    next
  }
  if (SDS[tot.fl+j, 1] %in% so.lv) {
    yyy <- SDS[which(SDS["lhs"] == SDS[tot.fl+j, 1] & SDS["op"] == "=~"), "rhs"]
    if (SDS[tot.fl+j, 3] %in% yyy) {
      next
    }
  }
  if (SDS[tot.fl+j, 3] %in% so.lv) {
    yyy <- SDS[which(SDS["lhs"] == SDS[tot.fl+j, 3] & SDS["op"] == "=~"), "rhs"]
    if (SDS[tot.fl+j, 1] %in% yyy) {
      next
    }
  }
  i <- k + 1
  DV.PCI[i, 1] <- paste(SDS[tot.fl+j,1], SDS[tot.fl+j,2], SDS[tot.fl+j,3], sep = " ")
  DV.PCI[i, 2] <- SDV[j,1]
  DV.PCI[i, 3] <- format(round(quantile(BDV[, j, 1], c(0.005)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  DV.PCI[i, 4] <- format(round(quantile(BDV[, j, 1], c(0.025)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  DV.PCI[i, 5] <- format(round(quantile(BDV[, j, 1], c(0.050)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  DV.PCI[i, 6] <- format(round(as.numeric(SDV[j, 3]), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  DV.PCI[i, 7] <- format(round(quantile(BDV[, j, 1], c(0.950)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  DV.PCI[i, 8] <- format(round(quantile(BDV[, j, 1], c(0.975)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  DV.PCI[i, 9] <- format(round(quantile(BDV[, j, 1], c(0.995)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  DV.PCI[i,10] <- format(round(DV.pp[j, 1], digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  k <- i + 1
  DV.PCI[k, 1] <- paste(SDS[tot.fl+j,1], SDS[tot.fl+j,2], SDS[tot.fl+j,3], sep = " ")
  DV.PCI[k, 2] <- SDV[j, 2]
  DV.PCI[k, 3] <- format(round(quantile(BDV[, j, 2], c(0.005)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  DV.PCI[k, 4] <- format(round(quantile(BDV[, j, 2], c(0.025)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  DV.PCI[k, 5] <- format(round(quantile(BDV[, j, 2], c(0.050)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  DV.PCI[k, 6] <- format(round(as.numeric(SDV[j, 4]), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  DV.PCI[k, 7] <- format(round(quantile(BDV[, j, 2], c(0.950)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  DV.PCI[k, 8] <- format(round(quantile(BDV[, j, 2], c(0.975)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  DV.PCI[k, 9] <- format(round(quantile(BDV[, j, 2], c(0.995)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  DV.PCI[k,10] <- format(round(DV.pp[j, 2], digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
}

# ===== Print percentile confidence intervals
if (PCI == "TRUE") {
  cat("\n", "   Percentile Confidence Intervals for Discriminant Validity", rep("\n",2))
  rownames(DV.PCI) <- rep("    ", nrow(DV.PCI))
  print(DV.PCI, quote=FALSE, right=TRUE)
  cat("\n")
}

# ===== Bias-Correction Factor for Discriminant Validity
DV.z <- matrix(0, tot.lv.cor, 2)
for (r in 1:tot.lv.cor) {
  DV.z[r, 1] <- qnorm(sum(BDV[, r, 1] < SDV[r, 3])/b.no)
  DV.z[r, 2] <- qnorm(sum(BDV[, r, 2] < SDV[r, 4])/b.no)
}

# ===== Bias-Corrected Probability for Discriminant Validity (DV.bcp)
DV.bcp <- matrix(0, tot.lv.cor, 2)

for (r in 1:tot.lv.cor) {
  if ((SDV[r, 1]>0 & min(BDV[, r,1]) > 0) | (SDV[r, 1] < 0 & max(BDV[, r, 1]) < 0)) {
      DV.bcp[r, 1] = 0
    } else if (qnorm(sum(BDV[, r, 1] > 0)/b.no)+2*DV.z[r, 1]<0) {
      DV.bcp[r, 1] = 2*pnorm((qnorm(sum(BDV[, r, 1]>0)/b.no)+2*DV.z[r, 1]))
    } else {
      DV.bcp[r, 1] = 2*pnorm(-1*(qnorm(sum(BDV[, r, 1]>0)/b.no)+2*DV.z[r, 1]))
  }
  if ((SDV[r, 2]>0 & min(BDV[, r, 2]) > 0) | (SDV[r, 2] < 0 & max(BDV[, r, 2]) < 0)) {
      DV.bcp[r, 2] = 0
    } else if (qnorm(sum(BDV[, r, 2] > 0)/b.no)+2*DV.z[r, 2]<0) {
      DV.bcp[r, 2] = 2*pnorm((qnorm(sum(BDV[, r, 2]>0)/b.no)+2*DV.z[r, 2]))
    } else {
      DV.bcp[r, 2] = 2*pnorm(-1*(qnorm(sum(BDV[, r, 2]>0)/b.no)+2*DV.z[r, 2]))
  }
}

# ===== Bias-Corrected Confidence Intervals for Discriminant Validity (DV.BCCI)
DV.BCCI <- matrix(0, nrow = tot.corX2, ncol=10)
colnames(DV.BCCI) <- c(" Correlation","      Factor","     0.5%","     2.5%","       5%","   Estimate","      95%","    97.5%","    99.5%","  p-value")

k <- 0
for (j in 1:tot.lv.cor) {
  if (SDS[tot.fl+j, 1] == SDS[tot.fl+j, 3]) {
    next
  }
  if (SDS[tot.fl+j, 1] %in% so.lv) {
    yyy <- SDS[which(SDS["lhs"] == SDS[tot.fl+j, 1] & SDS["op"] == "=~"), "rhs"]
    if (SDS[tot.fl+j, 3] %in% yyy) {
      next
    }
  }
  if (SDS[tot.fl+j, 3] %in% so.lv) {
    yyy <- SDS[which(SDS["lhs"] == SDS[tot.fl+j, 3] & SDS["op"] == "=~"), "rhs"]
    if (SDS[tot.fl+j, 1] %in% yyy) {
      next
    }
  }
  i <- k + 1
  DV.BCCI[i, 1] <- paste(SDS[tot.fl+j,1], SDS[tot.fl+j,2], SDS[tot.fl+j,3], sep = " ")
  DV.BCCI[i, 2] <- SDV[j, 1]
  DV.BCCI[i, 3] <- format(round(quantile(BDV[, j, 1],probs = pnorm(2*DV.z[j, 1]+qnorm(0.005))), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  DV.BCCI[i, 4] <- format(round(quantile(BDV[, j, 1],probs = pnorm(2*DV.z[j, 1]+qnorm(0.025))), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  DV.BCCI[i, 5] <- format(round(quantile(BDV[, j, 1],probs = pnorm(2*DV.z[j, 1]+qnorm(0.050))), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  DV.BCCI[i, 6] <- format(round(as.numeric(SDV[j, 3]), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  DV.BCCI[i, 7] <- format(round(quantile(BDV[, j, 1],probs = pnorm(2*DV.z[j, 1]+qnorm(0.950))), digits = R.decimal), nsmall = R.decimal, scientific = FALSE)
  DV.BCCI[i, 8] <- format(round(quantile(BDV[, j, 1],probs = pnorm(2*DV.z[j, 1]+qnorm(0.975))), digits = R.decimal), nsmall = R.decimal, scientific = FALSE)
  DV.BCCI[i, 9] <- format(round(quantile(BDV[, j, 1],probs = pnorm(2*DV.z[j, 1]+qnorm(0.995))), digits = R.decimal), nsmall = R.decimal, scientific = FALSE)
  DV.BCCI[i,10] <- format(round(DV.bcp[j, 1], digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  k <- i + 1
  DV.BCCI[k, 1] <- paste(SDS[tot.fl+j,1], SDS[tot.fl+j,2], SDS[tot.fl+j,3], sep = " ")
  DV.BCCI[k, 2] <- SDV[j, 2]
  DV.BCCI[k, 3] <- format(round(quantile(BDV[, j, 2],probs = pnorm(2*DV.z[j, 2]+qnorm(0.005))), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  DV.BCCI[k, 4] <- format(round(quantile(BDV[, j, 2],probs = pnorm(2*DV.z[j, 2]+qnorm(0.025))), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  DV.BCCI[k, 5] <- format(round(quantile(BDV[, j, 2],probs = pnorm(2*DV.z[j, 2]+qnorm(0.050))), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  DV.BCCI[k, 6] <- format(round(as.numeric(SDV[j, 4]), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
  DV.BCCI[k, 7] <- format(round(quantile(BDV[, j, 2],probs = pnorm(2*DV.z[j, 2]+qnorm(0.950))), digits = R.decimal), nsmall = R.decimal, scientific = FALSE)
  DV.BCCI[k, 8] <- format(round(quantile(BDV[, j, 2],probs = pnorm(2*DV.z[j, 2]+qnorm(0.975))), digits = R.decimal), nsmall = R.decimal, scientific = FALSE)
  DV.BCCI[k, 9] <- format(round(quantile(BDV[, j, 2],probs = pnorm(2*DV.z[j, 2]+qnorm(0.995))), digits = R.decimal), nsmall = R.decimal, scientific = FALSE)
  DV.BCCI[k,10] <- format(round(DV.bcp[j, 2], digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
}

# ===== Print bias-corrected confidence intervals
if (BCCI == "TRUE") {
  cat("\n", "   Bias-Corrected Confidence Intervals for Discriminant Validity", rep("\n",2))
  rownames(DV.BCCI) <- rep("    ", nrow(DV.BCCI))
  print(DV.BCCI, quote=FALSE, right=TRUE)
  cat("\n")
}

# ==================================================================================== #

# ============================ HTMT =============================== #
if (HTMT == "TRUE") {
  mHTMT <- matrix(nrow=no.factor, ncol=no.factor)
  colnames(mHTMT) <- names.lv
  rownames(mHTMT) <- names.lv
  m <- lavaan::inspect(Model.EST1, what = "sampstat")
  m <- m$cov

  for (i in 1:no.factor) {
    for (j in i:no.factor) {
      if (no.items[i] == 1) {
        MTi <- SDS[which(SDS["lhs"] == names.lv[i] & SDS["op"] == "=~"), "est.std"]^2*m[name.items[i],name.items[i]]
      } else {
        n <- m[name.items[i,1:no.items[i]], name.items[i,1:no.items[i]]]
        MTi <- mean(n[lower.tri(n)], na.rm = TRUE)
      }
      if (no.items[j] == 1) {
        MTj <- SDS[which(SDS["lhs"] == names.lv[j] & SDS["op"] == "=~"), "est.std"]^2*m[name.items[j,1:no.items[j]], name.items[j,1:no.items[j]]]
      } else {
        n <- m[name.items[j,1:no.items[j]], name.items[j,1:no.items[j]]]
        MTj <- mean(n[lower.tri(n)], na.rm = TRUE)
      }
      mHTMT[j,i] <- mean(m[name.items[j,1:no.items[j]], name.items[i,1:no.items[i]]], na.rm = TRUE)/(MTi*MTj)^0.5
    }
  }
  for (j in 1:no.factor) {
    mHTMT[j,j] <- 1
  }


  # ===== Percentile Probability for HTMT (HTMT.pp)
  HTMT.pp <- matrix(0, tot.lv.cor)  # Percentile Probability #
  for (r in 1:tot.lv.cor) {
    if (quantile(bootcoef[, tot.fl+tot.lv.cor+r],probs=0.5)>0) {
      HTMT.pp[r] = 2*(sum(bootcoef[, tot.fl+tot.lv.cor+r]<0)/b.no)
    } else {
      HTMT.pp[r] = 2*(sum(bootcoef[, tot.fl+tot.lv.cor+r]>0)/b.no)
    }
  }

  # ===== Percentile Confidence Intervals of HTMT (HTMT.PCI)
  HTMT.PCI <- matrix(0, nrow = (tot.lv.cor - no.factor - length(xxx)), ncol = 9)
  colnames(HTMT.PCI) <- c("         HTMT","     0.5%","     2.5%","       5%","   Estimate","      95%","    97.5%","    99.5%","  p-value")

  r <- 0
  k <- 0
  for (i in 1:no.factor) {
    for (j in i:no.factor) {
      k <- k + 1
      if (i == j) {
        next
      }
      if (colnames(mHTMT)[i] %in% so.lv) {
        yyy <- SDS[which(SDS["lhs"] == SDS[tot.fl+j, 1] & SDS["op"] == "=~"), "rhs"]
        if (colnames(mHTMT)[j] %in% yyy) {
          next
        }
      }
      if (colnames(mHTMT)[j] %in% so.lv) {
        yyy <- SDS[which(SDS["lhs"] == SDS[tot.fl+j, 3] & SDS["op"] == "=~"), "rhs"]
        if (colnames(mHTMT)[i] %in% yyy) {
          next
        }
      }
      r <- r + 1
      HTMT.PCI[r, 1] <- paste(colnames(mHTMT)[i], " ~~ ", colnames(mHTMT)[j], sep = " ")
      HTMT.PCI[r, 2] <- format(round(quantile(bootcoef[, tot.fl+tot.lv.cor+k], c(0.005)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
      HTMT.PCI[r, 3] <- format(round(quantile(bootcoef[, tot.fl+tot.lv.cor+k], c(0.025)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
      HTMT.PCI[r, 4] <- format(round(quantile(bootcoef[, tot.fl+tot.lv.cor+k], c(0.050)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
      HTMT.PCI[r, 5] <- format(round(mHTMT[j,i], digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
      HTMT.PCI[r, 6] <- format(round(quantile(bootcoef[, tot.fl+tot.lv.cor+k], c(0.950)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
      HTMT.PCI[r, 7] <- format(round(quantile(bootcoef[, tot.fl+tot.lv.cor+k], c(0.975)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
      HTMT.PCI[r, 8] <- format(round(quantile(bootcoef[, tot.fl+tot.lv.cor+k], c(0.995)), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
      HTMT.PCI[r, 9] <- format(round(HTMT.pp[k], digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
    }
  }

  # ===== Print percentile confidence intervals
  if (PCI == "TRUE") {
    cat("\n", "   Percentile Confidence Intervals for Disattenuated Correlation (Corrected HTMT)", rep("\n", 2))
    rownames(HTMT.PCI) <- rep("    ", nrow(HTMT.PCI))
    print(HTMT.PCI, quote=FALSE, right=TRUE)
    cat("\n")
  }

  # ===== Bias-Correction Factor for HTMT (HTMT.z)
  HTMT.z <- matrix(0, tot.lv.cor)  # Bias-Corrected Factor #
  r <- 0
  for (i in 1:no.factor) {
    for (j in i:no.factor) {
      r <- r + 1
      HTMT.z[r] <- qnorm(sum(bootcoef[, tot.fl+tot.lv.cor+r]<mHTMT[j,i])/b.no)
    }
  }

  # ===== Bias-Corrected Probability for HTMT
  HTMT.bcp <- matrix(0, tot.lv.cor)
  r <- 0
  for (i in 1:no.factor){
    for (j in i:no.factor) {
      r <- r + 1
      if ((mHTMT[j,i]>0 & min(bootcoef[, tot.fl+tot.lv.cor+r])>0) | (mHTMT[j,i]<0 & max(bootcoef[, tot.fl+tot.lv.cor+r])<0)) {
        HTMT.bcp[r] = 0
      } else if (qnorm(sum(bootcoef[, tot.fl+tot.lv.cor+r]>0)/b.no)+2*HTMT.z[r]<0) {
        HTMT.bcp[r] = 2*pnorm((qnorm(sum(bootcoef[, tot.fl+tot.lv.cor+r]>0)/b.no)+2*HTMT.z[r]))
      } else {
        HTMT.bcp[r] = 2*pnorm(-1*(qnorm(sum(bootcoef[, tot.fl+tot.lv.cor+r]>0)/b.no)+2*HTMT.z[r]))
      }
    }
  }

  # ===== Bias-Corrected Confidence Intervals for Correlation (HTMT.BCCI)
  HTMT.BCCI <- matrix(0, nrow = tot.lv.cor - no.factor - length(xxx), ncol=9)
  colnames(HTMT.BCCI) <- c("        HTMT","     0.5%","     2.5%","       5%","   Estimate","      95%","    97.5%","    99.5%","  p-value")

  r <- 0
  k <- 0
  for (i in 1:no.factor) {
    for (j in i:no.factor) {
      k <- k + 1
      if (i == j) {
        next
      }
      if (colnames(mHTMT)[i] %in% so.lv) {
        yyy <- SDS[which(SDS["lhs"] == SDS[tot.fl+j, 1] & SDS["op"] == "=~"), "rhs"]
        if (colnames(mHTMT)[j] %in% yyy) {
          next
        }
      }
      if (colnames(mHTMT)[j] %in% so.lv) {
        yyy <- SDS[which(SDS["lhs"] == SDS[tot.fl+j, 3] & SDS["op"] == "=~"), "rhs"]
        if (colnames(mHTMT)[i] %in% yyy) {
          next
        }
      }
      r <- r + 1
      HTMT.BCCI[r, 1] <- paste(colnames(mHTMT)[i], " ~~ ", colnames(mHTMT)[j], sep = " ")
      HTMT.BCCI[r, 2] <- format(round(quantile(bootcoef[, tot.fl+tot.lv.cor+k],probs = pnorm(2*HTMT.z[j]+qnorm(0.005))), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
      HTMT.BCCI[r, 3] <- format(round(quantile(bootcoef[, tot.fl+tot.lv.cor+k],probs = pnorm(2*HTMT.z[j]+qnorm(0.025))), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
      HTMT.BCCI[r, 4] <- format(round(quantile(bootcoef[, tot.fl+tot.lv.cor+k],probs = pnorm(2*HTMT.z[j]+qnorm(0.050))), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
      HTMT.BCCI[r, 5] <- format(round(mHTMT[j,i], digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
      HTMT.BCCI[r, 6] <- format(round(quantile(bootcoef[, tot.fl+tot.lv.cor+k],probs = pnorm(2*HTMT.z[j]+qnorm(0.950))), digits = R.decimal), nsmall = R.decimal, scientific = FALSE)
      HTMT.BCCI[r, 7] <- format(round(quantile(bootcoef[, tot.fl+tot.lv.cor+k],probs = pnorm(2*HTMT.z[j]+qnorm(0.975))), digits = R.decimal), nsmall = R.decimal, scientific = FALSE)
      HTMT.BCCI[r, 8] <- format(round(quantile(bootcoef[, tot.fl+tot.lv.cor+k],probs = pnorm(2*HTMT.z[j]+qnorm(0.995))), digits = R.decimal), nsmall = R.decimal, scientific = FALSE)
      HTMT.BCCI[r, 9] <- format(round(HTMT.bcp[k], digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
    }
  }

  # ===== Print bias-corrected confidence intervals
  if (BCCI == "TRUE") {
    cat("\n", "   Bias-Corrected Confidence Intervals for Disattenuated Correlation (Corrected HTMT)", rep("\n", 2))
    cat("\n")
    rownames(HTMT.BCCI) <- rep("    ", nrow(HTMT.BCCI))
    print(HTMT.BCCI, quote=FALSE, right=TRUE)
    cat("\n")
  }
}
# ==================================================================================== #

# ===== Print summary statistics and fit indices of unstandardized parameters

cat(rep("\n", 3))
cat("\n", "Estimates from Original Sample", rep("\n", 2))
print(lavaan::summary(Model.EST1, fit.measures=TRUE, standardized=TRUE, rsq=TRUE))
cat("\n")

# ========================= Generate Table using Bias-Corrected Confidence Intervals ====================== #

if (BCCI == "TRUE") {

  # ====================== Table 1. Standardized Factor Loadings ======================= #

  table.1 <- matrix(" ", nrow = tot.fl, ncol=3)
  colnames(table.1) <- c(" Factor Loading","  Estimate", "  ")

  for (r in 1:tot.fl) {
    table.1[r, 1] <- paste(SDS[r,1], SDS[r,2], SDS[r,3], sep = " ")
    table.1[r, 2] <- format(round(par.est[r], digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
    table.1[r, 3] <- " "
    if (quantile(bootcoef[, r],probs = pnorm(2*SFL.z[r]+qnorm(0.950))) < 0.4) {
      table.1[r, 3] <- "a"
    } else if (quantile(bootcoef[, r],probs = pnorm(2*SFL.z[r]+qnorm(0.950))) < 0.5) {
      table.1[r, 3] <- "b"
    } else if (quantile(bootcoef[, r],probs = pnorm(2*SFL.z[r]+qnorm(0.950))) < 0.7) {
      table.1[r, 3] <- "c"
    }
  }

  # ====================================================================================== #


  # ========================== Table 2 - Latent Variables ================================ #

  no.row = no.factor
  if (length(so.lv) != 0) no.row = no.factor + 1
  table.2 <- matrix(" ", nrow = no.row, ncol = no.factor + 4)
  colnames(table.2) <- rep("       ", no.factor + 4)

  if (no.lv.nox > 0) {
    colnames(table.2)[1] <- " First-order Factor"
  } else {
    colnames(table.2)[1] <- "  Factor"
  }
  colnames(table.2)[2] <- "    Mean"
  colnames(table.2)[3] <- "     s.d."
  colnames(table.2)[4] <- "      AVE   "
  for (factor.no in 1:no.factor) {
    colnames(table.2)[factor.no + 4] <- paste0("   ", names.lv[factor.no], "    ")
  }
  rownames(table.2) <- rep("    ", nrow(table.2))

  sh <- 0
  for (factor.no in 1:no.row) {
    if (sh == 0) {
      if (names.lv[factor.no] %in% so.lv) {
        table.2[factor.no, 1] <- "Second-order factor"
        sh <- 1
        next
      }
    }
    table.2[factor.no, 1] <- names.lv[factor.no-sh]
    table.2[factor.no, 2] <- format(round(mean(data.source[, names.lv[factor.no-sh]], na.rm = TRUE), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
    table.2[factor.no, 3] <- format(round(lv.cov[names.lv[factor.no-sh], names.lv[factor.no-sh]]^0.5, digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
    if (quantile(BAVE[, factor.no-sh], probs = pnorm(2*AVE.z[factor.no-sh]+qnorm(0.950))) < 0.49999) {
      table.2[factor.no, 4] <- paste0(format(round(SAVE[factor.no-sh], digits = R.decimal), nsmall = R.decimal, scientific=FALSE),"#  ")
    } else {
      table.2[factor.no, 4] <- paste0(format(round(SAVE[factor.no-sh], digits = R.decimal), nsmall = R.decimal, scientific=FALSE),"   ")
    }
    for (j in 1:(factor.no-sh)) {
      r <- (j-1)*(no.factor) - sum(1:j-1) + factor.no - sh
      if (j == (factor.no-sh)) {
        if (quantile(BCR[, factor.no-sh], probs = pnorm(2*CR.z[factor.no-sh]+qnorm(0.950))) < 0.69999) {
          table.2[factor.no, j+4] <- paste0("(",format(round(SCR[factor.no-sh], digits = R.decimal), nsmall = R.decimal, scientific=FALSE),")A  ")
        } else if (quantile(BCR[, factor.no-sh], probs = pnorm(2*CR.z[factor.no-sh]+qnorm(0.950))) < 0.79999) {
          table.2[factor.no, j+4] <- paste0("(",format(round(SCR[factor.no-sh], digits = R.decimal), nsmall = R.decimal, scientific=FALSE),")\u1D2E  ")
        } else {
        table.2[factor.no, j+4] <- paste0("(",format(round(SCR[factor.no-sh], digits = R.decimal), nsmall = R.decimal, scientific=FALSE),")   ")
        }
        next
      }
      if (names.lv[(factor.no-sh)] %in% so.lv) {
         if (names.lv[j] %in% names.lv.nox) {
          next
        }
      }
      table.2[factor.no, j+4] <- paste0(format(round(SDS[tot.fl+r, 4], digits = R.decimal), nsmall = R.decimal, scientific=FALSE), "    ")
      if (factor.no-sh != j) {
        k <- 0
        if (quantile(BDV[, r, 1],probs = pnorm(2*DV.z[r, 1]+qnorm(0.950))) < 0) {
          table.2[factor.no, j+4] <- paste0(format(round(SDS[tot.fl+r, 4], digits = R.decimal), nsmall = R.decimal, scientific=FALSE), "$   ")
          k <- 1
        }
        if (quantile(BDV[, r, 2],probs = pnorm(2*DV.z[r, 2]+qnorm(0.950))) < 0) {
          table.2[factor.no, j+4] <- paste0(format(round(SDS[tot.fl+r, 4], digits = R.decimal), nsmall = R.decimal, scientific=FALSE), "$   ")
          k <- 1
        }
        if (k == 1) {
          if (SDS[tot.fl+r , 4] > 0) {
            if (quantile(bootcoef[, tot.fl+r],probs = pnorm(2*COR.z[r]+qnorm(0.050))) > 0.850001) {
              table.2[factor.no, j+4] <- paste0(trimws(table.2[factor.no, j+4], "r"), "a  ")
            } else if (abs(quantile(bootcoef[, tot.fl+r],probs = pnorm(2*COR.z[r]+qnorm(0.050)))) > 0.800001) {
              table.2[factor.no, j+4] <- paste0(trimws(table.2[factor.no, j+4], "r"), "b  ")
            } else if (abs(quantile(bootcoef[, tot.fl+r],probs = pnorm(2*COR.z[r]+qnorm(0.050)))) > 0.700001) {
              table.2[factor.no, j+4] <- paste0(trimws(table.2[factor.no, j+4], "r"), "c  ")
            }
          } else {
            if (quantile(bootcoef[, tot.fl+r],probs = pnorm(2*COR.z[r]+qnorm(0.950))) < -0.849999) {
              table.2[factor.no, j+4] <- paste0(trimws(table.2[factor.no, j+4], "r"), "a  ")
            } else if ((quantile(bootcoef[, tot.fl+r],probs = pnorm(2*COR.z[r]+qnorm(0.950)))) < -0.799999) {
              table.2[factor.no, j+4] <- paste0(trimws(table.2[factor.no, j+4], "r"), "b  ")
            } else if ((quantile(bootcoef[, tot.fl+r],probs = pnorm(2*COR.z[r]+qnorm(0.950)))) < -0.699999) {
              table.2[factor.no, j+4] <- paste0(trimws(table.2[factor.no, j+4], "r"), "c  ")
            }
          }
        } else {
          if (SDS[tot.fl+r, 4] > 0) {
            if (quantile(bootcoef[, tot.fl+r],probs = pnorm(2*COR.z[r]+qnorm(0.050))) > 0.850001) {
              table.2[factor.no, j+4] <- paste0(trimws(table.2[factor.no, j+4], "r"), "a   ")
            } else if (abs(quantile(bootcoef[, tot.fl+r],probs = pnorm(2*COR.z[r]+qnorm(0.050)))) > 0.800001) {
              table.2[factor.no, j+4] <- paste0(trimws(table.2[factor.no, j+4], "r"), "b   ")
            } else if (abs(quantile(bootcoef[, tot.fl+r],probs = pnorm(2*COR.z[r]+qnorm(0.050)))) > 0.700001) {
              table.2[factor.no, j+4] <- paste0(trimws(table.2[factor.no, j+4], "r"), "c   ")
            }
          } else {
            if (quantile(bootcoef[, tot.fl+r],probs = pnorm(2*COR.z[r]+qnorm(0.950))) < -0.849999) {
              table.2[factor.no, j+4] <- paste0(trimws(table.2[factor.no, j+4], "r"), "a   ")
            } else if ((quantile(bootcoef[, tot.fl+r],probs = pnorm(2*COR.z[r]+qnorm(0.950)))) < -0.799999) {
              table.2[factor.no, j+4] <- paste0(trimws(table.2[factor.no, j+4], "r"), "b   ")
            } else if ((quantile(bootcoef[, tot.fl+r],probs = pnorm(2*COR.z[r]+qnorm(0.950)))) < -0.699999) {
              table.2[factor.no, j+4] <- paste0(trimws(table.2[factor.no, j+4], "r"), "c   ")
            }
          }
        }
      }
    }
  }
  # ====================================================================================== #


  # ========================== Table 3 - Observed Variables ============================== #

  # ===== Cronbach's alpha of summated scales (SCA)

  no.row = no.factor
  if (length(so.lv) != 0) no.row = no.factor + 1

  table.3 <- matrix(" ", nrow = no.row, ncol = no.factor + 4)
  colnames(table.3) <- rep("     ", no.factor + 4)
  if (no.lv.nox > 0) {
    colnames(table.3)[1] <- " First-Order Factor"
  } else {
    colnames(table.3)[1] <- "  Factor"
  }
  colnames(table.3)[2] <- "    Mean"
  colnames(table.3)[3] <- "     s.d."
  colnames(table.3)[4] <- "      AVE   "
  for (factor.no in 1:no.factor) {
    colnames(table.3)[factor.no + 4] <- paste0(names.lv[factor.no], "    ")
  }
  rownames(table.3) <- rep("    ", nrow(table.3))

  sh <- 0
  for (factor.no in 1:(no.row)) {
    if (sh == 0) {
      if (names.lv[factor.no] %in% so.lv) {
        table.3[factor.no, 1] <- "Second-order factor"
        sh <- 1
        next
      }
    }
    table.3[factor.no, 1] <- names.lv[factor.no-sh]
    table.3[factor.no, 2] <- format(round(mean(data.source[, names.lv[factor.no-sh]], na.rm = TRUE), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
    table.3[factor.no, 3] <- format(round(sd(data.source[, names.lv[factor.no-sh]], na.rm = TRUE), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
    if (quantile(BAVE[, factor.no-sh], probs = pnorm(2*AVE.z[factor.no-sh]+qnorm(0.950))) < 0.49999) {
      table.3[factor.no, 4] <- paste0(format(round(SAVE[factor.no-sh], digits = R.decimal), nsmall = R.decimal, scientific=FALSE),"#  ")
    } else {
      table.3[factor.no, 4] <- paste0(format(round(SAVE[factor.no-sh], digits = R.decimal), nsmall = R.decimal, scientific=FALSE),"   ")
    }
    for (j in 1:(factor.no-sh)) {
      if (j == (factor.no-sh)) {
        table.3[factor.no, j+4] <- paste0("(", format(round(SCA[factor.no-sh], digits = R.decimal), nsmall = R.decimal, scientific=FALSE), ")   ")
        next
      }
      if (names.lv[(factor.no-sh)] %in% so.lv) {
        if (names.lv[j] %in% names.lv.nox) {
          next
        }
      }
      table.3[factor.no, j+4] <- paste0(format(round(cor(data.source[, names.lv[factor.no-sh]], data.source[, names.lv[j]], use = "na.or.complete"), digits = R.decimal), nsmall = R.decimal, scientific=FALSE), "    ")
    }
  }

  # ========================== Table 4 - Disattenuated Correlation (Corrected HTMT) ============================== #
  if (HTMT == "TRUE") {
    no.row = no.factor
    if (length(so.lv) != 0) no.row = no.factor + 1
    table.4 <- matrix(" ", nrow = no.row, ncol = no.factor + 4)
    colnames(table.4) <- rep("     ", no.factor + 4)

    if (no.lv.nox > 0) {
      colnames(table.4)[1] <- " First-Order Factor"
    } else {
      colnames(table.4)[1] <- "  Factor"
    }
    colnames(table.4)[2] <- "    Mean"
    colnames(table.4)[3] <- "     s.d."
    colnames(table.4)[4] <- "      AVE   "
    for (factor.no in 1:no.factor) {
      colnames(table.4)[factor.no + 4] <- paste0(names.lv[factor.no], "    ")
    }
    rownames(table.4) <- rep("    ", nrow(table.4))

    sh <- 0
    for (factor.no in 1:(no.row)) {
      if (sh == 0) {
        if (names.lv[factor.no] %in% so.lv) {
          table.4[factor.no, 1] <- "Second-order factor"
          sh <- 1
          next
        }
      }
      table.4[factor.no, 1] <- names.lv[factor.no-sh]
      table.4[factor.no, 2] <- format(round(mean(data.source[, names.lv[factor.no-sh]], na.rm = TRUE), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
      table.4[factor.no, 3] <- format(round(sd(data.source[, names.lv[factor.no-sh]], na.rm = TRUE), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
      if (quantile(BAVE[, factor.no-sh], c(0.950)) < 0.49999) {
        table.4[factor.no, 4] <- paste0(format(round(SAVE[factor.no-sh], digits = R.decimal), nsmall = R.decimal, scientific=FALSE),"#  ")
      } else {
        table.4[factor.no, 4] <- paste0(format(round(SAVE[factor.no-sh], digits = R.decimal), nsmall = R.decimal, scientific=FALSE),"   ")
      }

      for (j in 1:(factor.no-sh)) {
        k <- (j-1)*(no.factor) - sum(1:j-1) + factor.no - sh
        if (names.lv[(factor.no-sh)] %in% so.lv) {
           if (names.lv[j] %in% names.lv.nox) {
            next
          }
        }
        table.4[factor.no, j+4] <- paste0(format(round(mHTMT[names.lv[factor.no-sh], names.lv[j]], digits = R.decimal), nsmall = R.decimal, scientific=FALSE), "    ")
        if (j == (factor.no-sh)) {
          table.4[factor.no, j+4] <- paste0("(", format(round(SCA[factor.no-sh], digits = R.decimal), nsmall = R.decimal, scientific=FALSE), ")   ")
          next
        }
        if (mHTMT[names.lv[factor.no-sh], names.lv[j]] > 0) {
          if (quantile(bootcoef[, tot.fl+tot.lv.cor+k], probs = pnorm(2*HTMT.z[k]+qnorm(0.050))) > 0.850001) {
            table.4[factor.no, j+4] <- paste0(trimws(table.4[factor.no, j+4], "r"), "a   ")
          } else if (abs(quantile(bootcoef[, tot.fl+tot.lv.cor+k], probs = pnorm(2*HTMT.z[k]+qnorm(0.050)))) > 0.800001) {
            table.4[factor.no, j+4] <- paste0(trimws(table.4[factor.no, j+4], "r"), "b   ")
          } else if (abs(quantile(bootcoef[, tot.fl+tot.lv.cor+k], probs = pnorm(2*HTMT.z[k]+qnorm(0.050)))) > 0.700001) {
            table.4[factor.no, j+4] <- paste0(trimws(table.4[factor.no, j+4], "r"), "c   ")
          }
        } else {
          if (quantile(bootcoef[, tot.fl+tot.lv.cor+k], probs = pnorm(2*HTMT.z[k]+qnorm(0.950))) < -0.849999) {
            table.4[factor.no, j+4] <- paste0(trimws(table.4[factor.no, j+4], "r"), "a   ")
          } else if (quantile(bootcoef[, tot.fl+tot.lv.cor+k], probs = pnorm(2*HTMT.z[k]+qnorm(0.950))) < -0.799999) {
            table.4[factor.no, j+4] <- paste0(trimws(table.4[factor.no, j+4], "r"), "b   ")
          } else if (quantile(bootcoef[, tot.fl+tot.lv.cor+k], probs = pnorm(2*HTMT.z[k]+qnorm(0.950))) < -0.699999) {
            table.4[factor.no, j+4] <- paste0(trimws(table.4[factor.no, j+4], "r"), "c   ")
          }
        }
      }
    }
  }
}
# ====================================================================================== #


# ========================= Generate Table using Percentile Confidence Intervals ====================== #

if (PCI == "TRUE") {

  # ====================== Table 1. Standardized Factor Loadings ======================= #

  table.1 <- matrix(" ", nrow = tot.fl, ncol=3)
  colnames(table.1) <- c(" Factor Loading","  Estimate", "  ")

  for (r in 1:tot.fl) {
    table.1[r, 1] <- paste(SDS[r,1], SDS[r,2], SDS[r,3], sep = " ")
    table.1[r, 2] <- format(round(par.est[r], digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
    table.1[r, 3] <- " "
    if (quantile(bootcoef[, r], c(0.950)) < 0.4) {
      table.1[r, 3] <- "a"
    } else if (quantile(bootcoef[, r], c(0.950)) < 0.5) {
      table.1[r, 3] <- "b"
    } else if (quantile(bootcoef[, r], c(0.950)) < 0.7) {
      table.1[r, 3] <- "c"
    }
  }

  # ====================================================================================== #


  # ========================== Table 2 - Latent Variables ================================ #

  no.row = no.factor
  if (length(so.lv) != 0) no.row = no.factor + 1
  table.2 <- matrix(" ", nrow = no.row, ncol = no.factor + 4)
  colnames(table.2) <- rep("       ", no.factor + 4)

  if (no.lv.nox > 0) {
    colnames(table.2)[1] <- "First-order Factor"
  } else {
    colnames(table.2)[1] <- "  Factor"
  }
  colnames(table.2)[2] <- "    Mean"
  colnames(table.2)[3] <- "     s.d."
  colnames(table.2)[4] <- "      AVE   "
  for (factor.no in 1:no.factor) {
    colnames(table.2)[factor.no + 4] <- paste0("   ", names.lv[factor.no], "    ")
  }
  rownames(table.2) <- rep("    ", nrow(table.2))

  sh <- 0
  for (factor.no in 1:no.row) {
    if (sh == 0) {
      if (names.lv[factor.no] %in% so.lv) {
        table.2[factor.no, 1] <- "Second-order factor"
        sh <- 1
        next
      }
    }
    table.2[factor.no, 1] <- names.lv[factor.no-sh]
    table.2[factor.no, 2] <- format(round(mean(data.source[, names.lv[factor.no-sh]], na.rm = TRUE), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
    table.2[factor.no, 3] <- format(round(lv.cov[names.lv[factor.no-sh], names.lv[factor.no-sh]]^0.5, digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
    if (quantile(BAVE[, factor.no-sh], c(0.950)) < 0.49999) {
      table.2[factor.no, 4] <- paste0(format(round(SAVE[factor.no-sh], digits = R.decimal), nsmall = R.decimal, scientific=FALSE),"#  ")
    } else {
      table.2[factor.no, 4] <- paste0(format(round(SAVE[factor.no-sh], digits = R.decimal), nsmall = R.decimal, scientific=FALSE),"   ")
    }
    for (j in 1:(factor.no-sh)) {
      r <- (j-1)*(no.factor) - sum(1:j-1) + factor.no - sh
      if (j == (factor.no-sh)) {
        if (quantile(BCR[, factor.no-sh], c(0.950)) < 0.69999) {
          table.2[factor.no, j+4] <- paste0("(",format(round(SCR[factor.no-sh], digits = R.decimal), nsmall = R.decimal, scientific=FALSE),")A  ")
        } else if (quantile(BCR[, factor.no-sh], c(0.950)) < 0.79999) {
          table.2[factor.no, j+4] <- paste0("(",format(round(SCR[factor.no-sh], digits = R.decimal), nsmall = R.decimal, scientific=FALSE),")B  ")
        } else {
        table.2[factor.no, j+4] <- paste0("(",format(round(SCR[factor.no-sh], digits = R.decimal), nsmall = R.decimal, scientific=FALSE),")   ")
        }
        next
      }
      if (names.lv[(factor.no-sh)] %in% so.lv) {
         if (names.lv[j] %in% names.lv.nox) {
          next
        }
      }
      table.2[factor.no, j+4] <- paste0(format(round(SDS[tot.fl+r, 4], digits = R.decimal), nsmall = R.decimal, scientific=FALSE), "    ")
      if (factor.no-sh != j) {
        k <- 0
        if (quantile(BDV[, r, 1],c(0.950)) < 0) {
          table.2[factor.no, j+4] <- paste0(format(round(SDS[tot.fl+r, 4], digits = R.decimal), nsmall = R.decimal, scientific=FALSE), "$   ")
          k <- 1
        }
        if (quantile(BDV[, r, 2], c(0.950)) < 0) {
          table.2[factor.no, j+4] <- paste0(format(round(SDS[tot.fl+r, 4], digits = R.decimal), nsmall = R.decimal, scientific=FALSE), "$   ")
          k <- 1
        }
        if (k == 1) {
          if (SDS[tot.fl+r , 4] > 0) {
            if (quantile(bootcoef[, tot.fl+r], c(0.050)) > 0.850001) {
              table.2[factor.no, j+4] <- paste0(trimws(table.2[factor.no, j+4], "r"), "a  ")
            } else if (abs(quantile(bootcoef[, tot.fl+r], c(0.050))) > 0.800001) {
              table.2[factor.no, j+4] <- paste0(trimws(table.2[factor.no, j+4], "r"), "b  ")
            } else if (abs(quantile(bootcoef[, tot.fl+r], c(0.050))) > 0.700001) {
              table.2[factor.no, j+4] <- paste0(trimws(table.2[factor.no, j+4], "r"), "c  ")
            }
          } else {
            if (quantile(bootcoef[, tot.fl+r], c(0.950)) < -0.849999) {
              table.2[factor.no, j+4] <- paste0(trimws(table.2[factor.no, j+4], "r"), "a  ")
            } else if (quantile(bootcoef[, tot.fl+r], c(0.950)) < -0.799999) {
              table.2[factor.no, j+4] <- paste0(trimws(table.2[factor.no, j+4], "r"), "b  ")
            } else if (quantile(bootcoef[, tot.fl+r], c(0.950)) < -0.699999) {
              table.2[factor.no, j+4] <- paste0(trimws(table.2[factor.no, j+4], "r"), "c  ")
            }
          }
        } else {
          if (SDS[tot.fl+r, 4] > 0) {
            if (quantile(bootcoef[, tot.fl+r], c(0.050)) > 0.850001) {
              table.2[factor.no, j+4] <- paste0(trimws(table.2[factor.no, j+4], "r"), "a   ")
            } else if (abs(quantile(bootcoef[, tot.fl+r], c(0.050))) > 0.800001) {
              table.2[factor.no, j+4] <- paste0(trimws(table.2[factor.no, j+4], "r"), "b   ")
            } else if (abs(quantile(bootcoef[, tot.fl+r], c(0.050))) > 0.700001) {
              table.2[factor.no, j+4] <- paste0(trimws(table.2[factor.no, j+4], "r"), "c   ")
            }
          } else {
            if (quantile(bootcoef[, tot.fl+r], c(0.950)) < -0.849999) {
              table.2[factor.no, j+4] <- paste0(trimws(table.2[factor.no, j+4], "r"), "a   ")
            } else if (quantile(bootcoef[, tot.fl+r], c(0.950)) < -0.799999) {
              table.2[factor.no, j+4] <- paste0(trimws(table.2[factor.no, j+4], "r"), "b   ")
            } else if (quantile(bootcoef[, tot.fl+r], c(0.950)) < -0.699999) {
              table.2[factor.no, j+4] <- paste0(trimws(table.2[factor.no, j+4], "r"), "c   ")
            }
          }
        }
      }
    }
  }
  # ====================================================================================== #


  # ========================== Table 3 - Observed Variables ============================== #

  no.row = no.factor
  if (length(so.lv) != 0) no.row = no.factor + 1

  table.3 <- matrix(" ", nrow = no.row, ncol = no.factor + 4)
  colnames(table.3) <- rep("     ", no.factor + 4)
  if (no.lv.nox > 0) {
    colnames(table.3)[1] <- "First-Order Factor"
  } else {
    colnames(table.3)[1] <- "  Factor"
  }
  colnames(table.3)[2] <- "    Mean"
  colnames(table.3)[3] <- "     s.d."
  colnames(table.3)[4] <- "      AVE   "
  for (factor.no in 1:no.factor) {
    colnames(table.3)[factor.no + 4] <- paste0(names.lv[factor.no], "    ")
  }
  rownames(table.3) <- rep("    ", nrow(table.3))

  sh <- 0
  for (factor.no in 1:(no.row)) {
    if (sh == 0) {
      if (names.lv[factor.no] %in% so.lv) {
        table.3[factor.no, 1] <- "Second-order factor"
        sh <- 1
        next
      }
    }
    table.3[factor.no, 1] <- names.lv[factor.no-sh]
    table.3[factor.no, 2] <- format(round(mean(data.source[, names.lv[factor.no-sh]], na.rm = TRUE), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
    table.3[factor.no, 3] <- format(round(sd(data.source[, names.lv[factor.no-sh]], na.rm = TRUE), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
    if (quantile(BAVE[, factor.no-sh], c(0.950)) < 0.49999) {
      table.3[factor.no, 4] <- paste0(format(round(SAVE[factor.no-sh], digits = R.decimal), nsmall = R.decimal, scientific=FALSE),"#  ")
    } else {
      table.3[factor.no, 4] <- paste0(format(round(SAVE[factor.no-sh], digits = R.decimal), nsmall = R.decimal, scientific=FALSE),"   ")
    }
    for (j in 1:(factor.no-sh)) {
      if (j == (factor.no-sh)) {
        table.3[factor.no, j+4] <- paste0("(", format(round(SCA[factor.no-sh], digits = R.decimal), nsmall = R.decimal, scientific=FALSE), ")   ")
        next
      }
      if (names.lv[(factor.no-sh)] %in% so.lv) {
        if (names.lv[j] %in% names.lv.nox) {
          next
        }
      }
      testcor <- cor.test(data.source[, names.lv[factor.no-sh]], data.source[, names.lv[j]])
      if (testcor$p.value < 0.01) {
        table.3[factor.no, j+4] <- paste0(format(round(cor(data.source[, names.lv[factor.no-sh]], data.source[, names.lv[j]], use = "na.or.complete"), digits = R.decimal), nsmall = R.decimal, scientific=FALSE), "**  ")
      } else if (testcor$p.value < 0.05) {
        table.3[factor.no, j+4] <- paste0(format(round(cor(data.source[, names.lv[factor.no-sh]], data.source[, names.lv[j]], use = "na.or.complete"), digits = R.decimal), nsmall = R.decimal, scientific=FALSE), "*   ")
      } else {
        table.3[factor.no, j+4] <- paste0(format(round(cor(data.source[, names.lv[factor.no-sh]], data.source[, names.lv[j]], use = "na.or.complete"), digits = R.decimal), nsmall = R.decimal, scientific=FALSE), "    ")
      }
    }
  }

  # ========================== Table 4 - Disattenuated Correlation (Corrected HTMT) ============================== #

  if (HTMT == "TRUE") {
    no.row = no.factor
    if (length(so.lv) != 0) no.row = no.factor + 1

    table.4 <- matrix(" ", nrow = no.row, ncol = no.factor + 4)
    colnames(table.4) <- rep("     ", no.factor + 4)
    if (no.lv.nox > 0) {
      colnames(table.4)[1] <- "First-Order Factor"
    } else {
      colnames(table.4)[1] <- "  Factor"
    }
    colnames(table.4)[2] <- "    Mean"
    colnames(table.4)[3] <- "     s.d."
    colnames(table.4)[4] <- "      AVE   "
    for (factor.no in 1:no.factor) {
      colnames(table.4)[factor.no + 4] <- paste0(names.lv[factor.no], "    ")
    }
    rownames(table.4) <- rep("    ", nrow(table.4))

    sh <- 0
    for (factor.no in 1:(no.row)) {
      if (sh == 0) {
        if (names.lv[factor.no] %in% so.lv) {
          table.4[factor.no, 1] <- "Second-order factor"
          sh <- 1
          next
        }
      }
      table.4[factor.no, 1] <- names.lv[factor.no-sh]
      table.4[factor.no, 2] <- format(round(mean(data.source[, names.lv[factor.no-sh]], na.rm = TRUE), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
      table.4[factor.no, 3] <- format(round(sd(data.source[, names.lv[factor.no-sh]], na.rm = TRUE), digits = R.decimal), nsmall = R.decimal, scientific=FALSE)
      if (quantile(BAVE[, factor.no-sh], c(0.950)) < 0.49999) {
        table.4[factor.no, 4] <- paste0(format(round(SAVE[factor.no-sh], digits = R.decimal), nsmall = R.decimal, scientific=FALSE),"#  ")
      } else {
        table.4[factor.no, 4] <- paste0(format(round(SAVE[factor.no-sh], digits = R.decimal), nsmall = R.decimal, scientific=FALSE),"   ")
      }
      for (j in 1:(factor.no-sh)) {
        k <- (j-1)*(no.factor) - sum(1:j-1) + factor.no - sh
        if (names.lv[(factor.no-sh)] %in% so.lv) {
           if (names.lv[j] %in% names.lv.nox) {
            next
          }
        }
        table.4[factor.no, j+4] <- paste0(format(round(mHTMT[names.lv[factor.no-sh], names.lv[j]], digits = R.decimal), nsmall = R.decimal, scientific=FALSE), "    ")
        if (j == (factor.no-sh)) {
          table.4[factor.no, j+4] <- paste0("(", format(round(SCA[factor.no-sh], digits = R.decimal), nsmall = R.decimal, scientific=FALSE), ")   ")
          next
        }
        if (mHTMT[names.lv[factor.no-sh], names.lv[j]] > 0) {
          if (quantile(bootcoef[, tot.fl+tot.lv.cor+k], c(0.050)) > 0.850001) {
            table.4[factor.no, j+4] <- paste0(trimws(table.4[factor.no, j+4], "r"), "a   ")
          } else if (abs(quantile(bootcoef[, tot.fl+tot.lv.cor+k], c(0.050))) > 0.800001) {
            table.4[factor.no, j+4] <- paste0(trimws(table.4[factor.no, j+4], "r"), "b   ")
          } else if (abs(quantile(bootcoef[, tot.fl+tot.lv.cor+k], c(0.050))) > 0.700001) {
            table.4[factor.no, j+4] <- paste0(trimws(table.4[factor.no, j+4], "r"), "c   ")
          }
        } else {
          if (quantile(bootcoef[, tot.fl+tot.lv.cor+k], c(0.950)) < -0.849999) {
            table.4[factor.no, j+4] <- paste0(trimws(table.4[factor.no, j+4], "r"), "a   ")
          } else if (quantile(bootcoef[, tot.fl+tot.lv.cor+k], c(0.950)) < -0.799999) {
            table.4[factor.no, j+4] <- paste0(trimws(table.4[factor.no, j+4], "r"), "b   ")
          } else if (quantile(bootcoef[, tot.fl+tot.lv.cor+k], c(0.950)) < -0.699999) {
            table.4[factor.no, j+4] <- paste0(trimws(table.4[factor.no, j+4], "r"), "c   ")
          }
        }
        if (names.lv[(factor.no-sh)] %in% so.lv) {
          if (names.lv[j] %in% names.lv.nox) {
            next
          }
        }
      }
    }
  }
}
# ====================================================================================== #




# ========================== Print Tables ================================ #

cat("\n", "   Table 1.  Standardized Factor Loadings", rep("\n", 2))
rownames(table.1) <- rep("    ", nrow(table.1))
print(table.1, quote=FALSE, right=TRUE)
cat(rep("\n", 2), "  Note: a = standardized factor loading significantly less than 0.4;","\n")
cat("         b = standardized factor loading significantly less than 0.5;","\n")
cat("         c = standardized factor loading significantly less than 0.7 (p < .05)", rep("\n", 3))

cat("\n", "  Table 2. Descriptive Statistics (Observed Mean, Latent s.d., AVE, Construct Reliability, Latent Correlation)", rep("\n", 2))
rownames(table.2) <- rep("    ", nrow(table.2))
print(table.2, quote=FALSE, right=TRUE)
cat("\n", "  Note: AVE = Average Variance Extracted; # = AVE significantly lower than 0.5 (p < .05)","\n")
if (length(so.lv) == 0) {
  cat("         diagonal elements in brackets = Construct Reliability","\n")
} else {
  cat(paste0("         diagonal elements in brackets = Construct Reliability for first-order factor and ", omega, " for second-order factor"), "\n")
}
cat("         A = Construct Reliability significantly lower than 0.7; B = Construct Reliability significantly lower than 0.8 (p < .05)", "\n")
cat("         Correlation coefficient: a = significantly larger than 0.85; b = significantly larger than 0.8; c = significantly larger than 0.7 (p < .05)", "\n")
if (length(so.lv) == 0) {
  cat("         $ = AVE is significantly less than squared-correlation (p < .05)", rep("\n", 3))
} else {
  cat("         $ = AVE is significantly less than squared-correlation (p < .05)", rep("\n"))
  cat("         Observed mean of second-order factor is based on all items ignoring the second-order structure", rep("\n", 3))
}



cat("\n", "  Table 3.  Descriptive Statistics (Observed Mean, Observed s.d., AVE, Reliability, Observed Correlation)", rep("\n", 2))
print(table.3, quote=FALSE, right=TRUE)
cat("\n", "  Note: AVE = Average Variance Extracted; # = AVE significantly lower than 0.5 (p < .05)","\n")
cat("         * p < .05; ** p < .01", "\n")
if (length(so.lv) == 0) {
  cat("         Diagonal elements in brackets = Cronbach's Alpha", rep("\n", 3))
} else {
  cat("         Diagonal elements in brackets = Cronbach's Alpha", rep("\n"))
  cat("         Observed mean, s.d., correlation and Cronbach's Alpha of second-order factor are based on all items ignoring the second-order structure", rep("\n", 3))
}

if (HTMT == "TRUE") {
  cat("\n", "  Table 4.  Descriptive Statistics (Observed Mean, Observed s.d., AVE, Reliability, Disattenuated Correlation)", rep("\n", 2))
  print(table.4, quote=FALSE, right=TRUE)
  cat("\n", "  Note: AVE = Average Variance Extracted; # = AVE significantly lower than 0.5 (p < .05)","\n")
  cat("         Disattenuated correlation: a = significantly larger than 0.85; b = significantly larger than 0.8; c = significantly larger than 0.7 (p < .05)", "\n")
  if (length(so.lv) == 0) {
    cat("         Diagonal elements in brackets = Cronbach's Alpha", rep("\n", 3))
  } else {
    cat("         Diagonal elements in brackets = Cronbach's Alpha", rep("\n"))
    cat("         Observed mean, s.d., correlation and Cronbach's Alpha of second-order factor are based on all items ignoring the second-order structure", rep("\n", 3))
  }
}

# ========================== Create Excel File ================================ #
## create workbook
wb <- openxlsx::createWorkbook()

## Add 3 worksheets
openxlsx::addWorksheet(wb, "Table 1")
openxlsx::addWorksheet(wb, "Table 2")
openxlsx::addWorksheet(wb, "Table 3")

## Write Table 1 ##
H.Table.1 <- c("Table 1.  Standardized Factor Loadings")
openxlsx::writeData(wb, "Table 1", H.Table.1, startCol = 1, startRow = 1, rowNames = FALSE)
openxlsx::writeData(wb, "Table 1", table.1, startCol = 1, startRow = 3, rowNames = FALSE)

## Right justify and set width of first column ##
openxlsx::addStyle(wb = wb, sheet = 'Table 1', rows = 1:nrow(table.1)+3, cols = 2L, style = createStyle(halign = 'right'))
openxlsx::setColWidths(wb, sheet = 'Table 1', cols = c(1), widths = (max(nchar(table.1[,1]))+3))

## Create remarks
Note.Table.1 <- matrix(3,)
Note.Table.1[1] <- c("Note: a = standardized factor loading significantly less than 0.4;")
Note.Table.1[2] <- c("           b = standardized factor loading significantly less than 0.5;")
Note.Table.1[3] <- c("           c = standardized factor loading significantly less than 0.7 (p < .05)")
openxlsx::writeData(wb, "Table 1", Note.Table.1, startCol = 1, startRow = nrow(table.1)+5, rowNames = FALSE)

## Write Table 2 ##
H.Table.2 <- c("Table 2. Descriptive Statistics (Observed Mean, Latent s.d., AVE, Construct Reliability, Latent Correlation)")
openxlsx::writeData(wb, "Table 2", H.Table.2, startCol = 1, startRow = 1, rowNames = FALSE)
openxlsx::writeData(wb, "Table 2", table.2, startCol = 1, startRow = 3, rowNames = FALSE)

## Right justify and set width of first column ##
openxlsx::addStyle(wb = wb, sheet = 'Table 2', rows = 3:(nrow(table.2)+3), cols = 2:(ncol(table.2)+1),
         style = createStyle(halign = 'right'), gridExpand = T)
openxlsx::setColWidths(wb, sheet = 'Table 2', cols = c(1), widths = (max(nchar(table.2[,1]))+3))

## Create remarks
Note.Table.2 <- matrix(5,)
Note.Table.2[1] <- c("Note: AVE = Average Variance Extracted; # = AVE significantly lower than 0.5 (p < .05)")
Note.Table.2[2] <- c("           diagonal elements in brackets = Construct Reliability")
Note.Table.2[3] <- c("           A = Construct Reliability significantly lower than 0.7; B = Construct Reliability significantly lower than 0.8 (p < .05)")
Note.Table.2[4] <- c("           Correlation coefficient: a = significantly larger than 0.85; b = significantly larger than 0.8; c = significantly larger than 0.7 (p < .05)")
Note.Table.2[5] <- c("           $ = AVE is significantly less than squared-correlation (p < .05)")
openxlsx::writeData(wb, "Table 2", Note.Table.2, startCol = 1, startRow = nrow(table.2)+5, rowNames = FALSE)


## Write Table 3 ##
H.Table.3 <- c("Table 3.  Descriptive Statistics (Observed Mean, Observed s.d., AVE, Reliability, Observed Correlation)")
openxlsx::writeData(wb, "Table 3", H.Table.3, startCol = 1, startRow = 1, rowNames = FALSE)
openxlsx::writeData(wb, "Table 3", table.3, startCol = 1, startRow = 3, rowNames = FALSE)

## Right justify and set width of first column ##
openxlsx::addStyle(wb = wb, sheet = 'Table 3', rows = 3:(nrow(table.3)+3), cols = 2:(ncol(table.3)+1),
         style = createStyle(halign = 'right'), gridExpand = T)
openxlsx::setColWidths(wb, sheet = 'Table 3', cols = c(1), widths = (max(nchar(table.3[,1]))+3))

## Create remarks
if(sh == 0) {
  Note.Table.3 <- matrix(3,)
  Note.Table.3[1] <- c("Note: AVE = Average Variance Extracted; # = AVE significantly lower than 0.5 (p < .05)")
  Note.Table.3[2] <- c("           * p < .05; ** p < .01")
  Note.Table.3[3] <- c("           Diagonal elements in brackets = Cronbach's Alpha")
} else {
  Note.Table.3 <- matrix(4,)
  Note.Table.3[1] <- c("Note: AVE = Average Variance Extracted; # = AVE significantly lower than 0.5 (p < .05)")
  Note.Table.3[2] <- c("             * p < .05; ** p < .01")
  Note.Table.3[3] <- c("           Diagonal elements in brackets = Cronbach's Alpha")
  Note.Table.3[4] <- c("           Observed mean, s.d., correlation and Cronbach's Alpha of second-order factor are based on all items ignoring the second-order structure")
}
openxlsx::writeData(wb, "Table 3", Note.Table.3, startCol = 1, startRow = nrow(table.3)+5, rowNames = FALSE)


## Write Table 4 ##
if (HTMT == "TRUE") {
  openxlsx::addWorksheet(wb, "Table 4")
  H.Table.4 <- c("Table 4.  Descriptive Statistics (Observed Mean, Observed s.d., AVE, Reliability, Disattenuated Correlation)")
  openxlsx::writeData(wb, "Table 4", H.Table.4, startCol = 1, startRow = 1, rowNames = FALSE)
  openxlsx::writeData(wb, "Table 4", table.4, startCol = 1, startRow = 3, rowNames = FALSE)

  ## Right justify and set width of first column ##
  openxlsx::addStyle(wb = wb, sheet = 'Table 4', rows = 3:(nrow(table.4)+3), cols = 2:(ncol(table.4)+1),
           style = createStyle(halign = 'right'), gridExpand = T)
  openxlsx::setColWidths(wb, sheet = 'Table 4', cols = c(1), widths = (max(nchar(table.4[,1]))+3))

  ## Create remarks
  if(sh == 0){
    Note.Table.4 <- matrix(3,)
    Note.Table.4[1] <- c("Note: AVE = Average Variance Extracted; # = AVE significantly lower than 0.5 (p < .05)")
    Note.Table.4[2] <- c("           Disattenuated correlation: a = significantly larger than 0.85; b = significantly larger than 0.8; c = significantly larger than 0.7 (p < .05)")
    Note.Table.4[3] <- c("           Diagonal elements in brackets = Cronbach's Alpha")
  } else {
    Note.Table.4 <- matrix(4,)
    Note.Table.4[1] <- c("Note: AVE = Average Variance Extracted; # = AVE significantly lower than 0.5 (p < .05)")
    Note.Table.4[2] <- c("           Disattenuated correlation: a = significantly larger than 0.85; b = significantly larger than 0.8; c = significantly larger than 0.7 (p < .05)")
    Note.Table.4[3] <- c("           Diagonal elements in brackets = Cronbach's Alpha")
    Note.Table.4[4] <- c("           Observed mean, s.d., correlation and Cronbach's Alpha of second-order factor are based on all items ignoring the second-order structure")
  }
  openxlsx::writeData(wb, "Table 4", Note.Table.4, startCol = 1, startRow = nrow(table.4)+5, rowNames = FALSE)
}

## Save Workbook ##
openxlsx::saveWorkbook(wb, "measureQ.xlsx", overwrite = TRUE)

cat("\n", "## == Output tables are saved in MeasureQ.xlsx == ##", rep("\n", 2))

}






