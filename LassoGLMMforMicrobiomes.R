### LassoGLMM for Microbiome Studies ###
### Written by Laura Tipton          ###
### Last edited: August 2, 2017      ###

## Data should be in the following formats:
# MBdat: 16S/ITS relative abundance data in 1 matrix, samples in rows and OTUs/species in columns with identifiable names
# dat: continuous response varibles in 1 matrix, samples in rows and variables in columns with identifiable names
# demos: categorical explanatory variables in 1 matrix, samples in rows and variables in columns with identifiable names
# ids: identification random effect variables in 1 matrix, samples in rows and variables in columns with identifiable names
# all rows in the above 4 matrices should be in the same order

## Data cleanup, skip if data is already clean
# relabun = function to calculate relative abundance, if not already in this format
relabun <- function(x){
  sums <- apply(x, 1, sum, na.rm=TRUE)
  y <- x/sums
  return(y)
}

# gt0 = function to count samples that contain OTUs 
gt0 <- function(vec){
  v <- as.numeric(vec)
  s <- sum(v>0)
  return(s)
}

# remove OTUs in less than 2 samples by applying gt0 function
MBdat.gt0 <- as.matrix(apply(MBdat, 2, gt0))
MBdat2 <- MBdat[,-which(dat.gt0<2)]

## Variable screening step based on Correlations
# UPDATE 5/16: added option to specify correation test and changed default to Spearman
# corrpairs = function to calculate correlations between all OTU-response variable pairs
corrpairs <- function(ys, xs, fName="Correlations.csv", useQ=FALSE, met="spearman"){
  sums <- apply(xs, 2, gt0)
  res <- vector("list", length(ncol(ys)))
  #names(res) <- colnames(ys)
  for(i in 1:ncol(ys)){
    res[[i]] <- vector("list")
    for(j in 1:ncol(xs)){
      c <- cor.test(ys[,i], xs[,j], na.rm=TRUE, method=met)
      p <- c$p.value
      q <- c$p.value*(ncol(ys)*ncol(xs))
      if (!is.na(p)){
        c2 <- c(colnames(ys)[i], colnames(xs)[j], as.numeric(c$estimate), c$conf.in, p, q, sums[[j]])
        write(c2, file=fName, ncolumns=8, append=TRUE, sep=",") 
        if (useQ){
          if (q < 0.05){
            res[[i]] <- append(res[[i]], j)
          }
        }
        else {
          if (p < 0.05){
            res[[i]] <- append(res[[i]], j)
          }
        }
      }
    }
    print(paste("Completed y variable ", i, ", ", colnames(ys)[i]))
  }
  return(res)
}

# calculate correlations, in order to move on, assign this to a variable (corrs)
corr <- corrpairs(dat, MBdat2, fName="Correlations.csv")

## Perform LassoGLMM
# penGLMM = function to regress MBdat on dat accounting for demos and ids
penGLMM <- function(ys, xs, corrs, randE, fName='Regression.txt', lam=seq(0,200,1), strat=NULL, rtrnmod=FALSE){
  require(glmmLasso)
  for (i in 1:ncol(ys)){
    vars <- unlist(corrs[[i]])
    vars2 <- colnames(xs)[vars]
    vars3 <- ''
    for (j in 1:length(vars2)){
      vars3 <- paste(vars3, vars2[j], sep="+")
    }
    tmp <- data.frame(na.omit(cbind(ys[,i], xs, randE)))
    colnames(tmp) <- c(colnames(ys)[i], colnames(xs), colnames(randE))
    ranEf <- list()
    for(k in 1:ncol(randE)){ ranEf <- append(ranEf, as.formula(paste(colnames(randE)[k], "=~1"))); names(ranEf)[[k]] <- colnames(randE)[k]}
    if (!is.null(strat)){
      tmp <- data.frame(na.omit(cbind(ys[,i], xs, randE, strat)))
      colnames(tmp) <- c(colnames(ys)[i], colnames(xs), colnames(randE), colnames(strat))
      for (k in 1:ncol(strat)){
        vars3 <- paste(vars3, paste0("as.factor(",colnames(strat)[k],")"), sep="+")
      }
    }
    tmp[,c(1:ncol(xs)+1)] <- apply(tmp[,c(1:ncol(xs)+1)],2, as.numeric)
    tmp[,c((ncol(xs)+2):(ncol(xs)+ncol(randE)+1))] <- apply(tmp[,c((ncol(xs)+2):(ncol(xs)+ncol(randE)+1))], 2, as.factor)
    min <- Inf
    lamb <- 0
    minmod <- list()
    minmod$coefficients <- 0
    minmod$ranef <- 0
    tmp[,1] <- as.numeric(as.character(tmp[,1]))
    for (l in lam){
      try({
        mod <- glmmLasso(fix=as.formula(paste("tmp[,1]~", substr(vars3,2,nchar(vars3)))), rnd=ranEf, data=data.frame(tmp), lambda=l, control=list(q_start=diag(0.1, ncol(randE))))
        if (mod$bic < min){
          minmod <- mod
          min <- mod$bic
          lamb <- l
        }
      }, silent=T)
    }
    write(paste("Y =", colnames(ys)[i]), file=fName, append=T)
    write("Fixed Effects:", file=fName, append=T)
    write.table(as.matrix(minmod$coefficients[abs(minmod$coefficients)>0]), file=fName, append=T)
    write("Random Effects:", file=fName, append=T)
    write.table(as.matrix(minmod$ranef), file=fName, append=T)
    write(paste("Optimal Lamba: ", lamb), file=fName, append=T)
    write(paste("Minimum BIC: ", min), file=fName, append=T)
    write(paste(""), file=fName, append=T)
    print(paste("Completed y variable ", i, ", ", colnames(ys)[i]))
    if (rtrnmod){ return(minmod) }
  }
}

# apply LassoGLMM, this does not need to be assigned to a variable
penGLMM(dat, MBdat2, corr, ids, strat=demos)

## Plotting example
require(car)
par(family="sans")

# create a temporary dataset sorted by response variable of interest (using 1 in this example and assuming OTU-1 is strongly associated)
tmpplot <- cbind(MBdat2[order(dat[,1]),], dat[order(dat[,1]),1])

# plot a "none" plot to set axes and labels
matplot(tmpplot[,1], tmpplot[,ncol(tmpplot)], pch=19, type="n", ylab="Response Variable", xlab="log relative abundance", main="OTU-1")

# plot grey dashed lines for all responses
for(i in 1:nrow(tmpplot)){ lines(c(-20,20), c(tmpplot[i,ncol(tmpplot)], tmpplot[i,ncol(tmpplot)]), lty=2, col="grey")}

# finally plot abundances in red
matplot(tmpplot[,1], tmpplot[,ncol(tmpplot)], pch=19, type="o", add=TRUE, col="red")

## Examine R^2 for overfitting
# Added 8/17
require(piecewiseSEM)
require(lme4)

# create temporary dataset with all relevant variables, including random effects; exclude observations with missing response variables
tmpmod <- data.frame(na.omit(cbind(MBdat2[,1], dat, ids, demos)))

# run linear mixed model with only stongly associated OTUs (again assuming OTU-1 is strongly associated with response variable 1)
# I found it easier to use variable names in this step so have used DEMO1 and ID1 for variables in demos and ids respectively
mod1 <- lmer(tmpmod[,1]~OTU_1+as.factor(DEMO1)+(1|ID1), data=tmpmod)

# calculate marginal and conditional R^2, columns 5 and 6, respectively
rsq <- sem.model.fits(mod1)

## Residual Plots; must be run after examining R^2 so that you have mod1 in the environment
# Added 8/17

plot(tmpmod[,1], predict(mod1), main="Variable 1", xlab="Observed", ylab="Predicted")

# add red line where prediction equals observation
abline(0,1, col="red", lty=2)

# add correlation value; x and y coordinates will need to be updated to an empty spot on the plot
text(x=0, y=0, paste("corr=", round(cor(tmpmod[,1], predict(mod1)), 3)))