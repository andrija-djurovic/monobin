#' Four-stage monotonic binning procedure including regression with nested dummies
#'
#' \code{ndr.bin} implements extension of three-stage monotonic binning procedure (\code{\link{iso.bin}})
#' with step of regression with nested dummies as fourth stage. 
#' The first stage is isotonic regression used to achieve the monotonicity. The next two stages are possible corrections for
#' minimum percentage of observations and target rate, while the last regression stage is used to identify
#' statistically significant cut points.
#'
#'@seealso \code{\link{iso.bin}} for three-stage monotonic binning procedure.
#'
#'@param x Numeric vector to be binned.
#'@param y Numeric target vector (binary or continuous).
#'@param sc Numeric vector with special case elements. Default values are \code{c(NA, NaN, Inf)}.
#' Recommendation is to keep the default values always and add new ones if needed. Otherwise, if these values exist
#' in \code{x} and are not defined in the \code{sc} vector, function will report the error. 
#'@param sc.method Define how special cases will be treated, all together or separately.
#' Possible values are \code{"together", "separately"}.
#'@param y.type Type of \code{y}, possible options are \code{"bina"} (binary) and \code{"cont"} (continuous).
#' If default value is passed, then algorithm will identify if y is 0/1 or continuous variable.
#'@param min.pct.obs Minimum percentage of observations per bin. Default is 0.05 or 30 observations.
#'@param min.avg.rate Minimum \code{y} average rate. Default is 0.05 or 30 observations.
#'@param p.val Threshold for p-value of regression coefficients. Default is 0.05.
#' For a binary target binary logistic regression is estimated, whereas for a continuous target, linear regression is used.
#'@param force.trend If the expected trend should be forced. Possible values: \code{"i"} for
#' increasing trend (\code{y} increases with increase of \code{x}), \code{"d"} for decreasing trend 
#' (\code{y} decreases with decrease of \code{x}). Default value is \code{NA}.
#' If the default value is passed, then trend will be identified based on the sign of the Spearman correlation 
#' coefficient between \code{x} and \code{y} on complete cases.
#'
#'@return The command \code{ndr.bin} generates a list of two objects. The first object, data frame \code{summary.tbl}
#' presents a summary table of final binning, while \code{x.trans} is a vector of discretized values.
#' In case of single unique value for \code{x} or \code{y} of complete cases (cases different than special cases), 
#' it will return data frame with info. 
#'
#'@examples
#' suppressMessages(library(monobin))
#' data(gcd)
#' age.bin <- ndr.bin(x = gcd$age, y = gcd$qual)
#' age.bin[[1]]
#' table(age.bin[[2]])
#' #linear regression example
#' amount.bin <- ndr.bin(x = gcd$amount, y = gcd$qual, y.type = "cont", p.val = 0.05)
#' #create nested dummies
#' db.reg <- gcd[, c("qual", "amount")]
#' db.reg$amount.bin <- amount.bin[[2]]
#' amt.s <- db.reg %>% 
#'	      group_by(amount.bin) %>% 
#'	      summarise(qual.mean = mean(qual),
#'			    amt.min = min(amount))
#' mins <- amt.s$amt.min
#' for (i in 2:length(mins)) {
#' 	 level.l <- mins[i]
#'	 nd <- ifelse(db.reg$amount < level.l, 0, 1)
#'	 db.reg <- cbind.data.frame(db.reg, nd)
#'	 names(db.reg)[ncol(db.reg)] <- paste0("dv_", i)
#'	 }
#' reg.f <- paste0("qual ~ dv_2 + dv_3")
#' lrm <- lm(as.formula(reg.f), data = db.reg)
#' lr.coef <- data.frame(summary(lrm)$coefficients)
#' lr.coef
#' cumsum(lr.coef$Estimate)
#' #check
#' as.data.frame(amt.s)
#' diff(amt.s$qual.mean)
#'
#'@importFrom stats cor isoreg pt sd weighted.mean
#'@importFrom Hmisc cut2
#'@import dplyr
#'@export

ndr.bin <- function(x, y, sc = c(NA, NaN, Inf), sc.method = "together", y.type = NA, 
			 min.pct.obs = 0.05, min.avg.rate = 0.01, p.val = 0.05, force.trend = NA) {
	options(scipen = 20)

	checks.init(x = x, y = y, sc = sc, sc.method = sc.method, 
			y.type = y.type, force.trend = force.trend)

	d <- data.frame(y, x)
	d <- d[!is.na(y), ]
	d.sc <- d[d$x%in%sc, ]
	d.cc <- d[!d$x%in%sc, ]

	checks.res <- checks.iter(d = d, d.cc = d.cc, y.type = y.type)
	if	(checks.res[[1]] > 0) {
		return(eval(parse(text = checks.res[[2]])))
		} 
	y.check <- checks.res[[3]]

	nr <- nrow(d)
	min.obs <- ceiling(ifelse(nr * min.pct.obs < 30, 30, nr * min.pct.obs))
	if	(y.check == "bina") {
		nd <- sum(d$y)
		min.rate <- ceiling(ifelse(nd * min.avg.rate < 1, 1, nd * min.avg.rate))		
		} else {
		min.rate <- min.avg.rate
		}
	ds <- iso(tbl.sc = d.sc, tbl.cc = d.cc, method = sc.method, min.obs = min.obs, 
		    min.rate = min.rate, y.check = y.check, force.trend = force.trend)
	ds.sc <- ds[ds$type%in%"special cases", ]
	ds.cc <- ds[ds$type%in%"complete cases", ]
	ds.cc <- ndr(tbl = ds.cc, mdb = d.cc, p.val = p.val, y.check = y.check)
	ds <- bind_rows(ds.sc, ds.cc)
	ds <- woe.calc(tbl = ds, y.check = y.check)
	sc.u <- unique(sc)
	sc.g <- ds$bin[ds$type%in%"special cases"]
	x.mins <- ds$x.min[!ds$bin%in%sc.u & !ds$bin%in%"SC"]
	x.maxs <- ds$x.max[!ds$bin%in%sc.u & !ds$bin%in%"SC"]
	x.trans <- slice.variable(x.orig = d$x, x.lb = x.mins, x.ub = x.maxs, 
					  sc.u = sc.u, sc.g = sc.g) 
return(list(summary.tbl = ds, x.trans = x.trans))
}
#regression on nested dummies
ndr <- function(tbl, mdb, p.val, y.check) {
	if	(nrow(tbl) == 1) {
		return(tbl)
		}
	level.u <- sort(tbl$x.min)
	for	(i in 2:length(level.u)) {
		level.l <- level.u[i]
		nd <- ifelse(mdb$x < level.l, 0, 1)
		mdb <- cbind.data.frame(mdb, nd)
		names(mdb)[ncol(mdb)] <- paste0("dv_", i)
		}
	nst.d <- names(mdb)[grepl("dv_", names(mdb))]
	lrm <- stepwise.reg(mdb = mdb, nst.d = nst.d, p.val = p.val, y.check = y.check)
	if	(length(lrm) > 0) {
		indx <- sapply(strsplit(lrm, "_"), "[", 2)
		cp <- level.u[as.numeric(indx)]
		mdb$bin <- cut2(mdb$x, cuts = cp)
		} else {
		mdb$bin <- "NS"
		}
	res <- iso.summary(tbl = mdb[, c("y", "x", "bin")], bin = "bin")
	res$bin <- format.bin(x.lb = res$x.min, x.ub = res$x.max)	
	res$type <- "complete cases"		
return(as.data.frame(res))
}
#stepwise regressions
stepwise.reg <- function(mdb, nst.d, p.val, y.check) {
	if	(y.check == "bina") {
		reg.exp <- "glm(as.formula(reg.f), family = binomial(link = logit), data = mdb)"
		} else {
		reg.exp <- "lm(as.formula(reg.f), data = mdb)"
		}
	repeat {
		 reg.f <- paste0("y ~ ", paste(nst.d , collapse = " + "))
		 lr <- eval(parse(text = reg.exp))
		 lr.coef <- data.frame(summary(lr)$coefficients)
		 p.vals <- lr.coef[-1, 4]
		 cond <- all(p.vals < p.val)
		 if	(cond) {break}
		 p.max <- which.max(p.vals)
		 nst.d <- nst.d[-p.max]
		 if	(length(nst.d) == 0) {break}
		 }
return(nst.d)
}





