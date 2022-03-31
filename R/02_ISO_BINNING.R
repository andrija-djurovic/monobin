#' Three-stage monotonic binning procedure
#'
#' \code{iso.bin} implements three-stage monotonic binning procedure. The first stage is isotonic regression
#' used to achieve the monotonicity, while the remaining two stages are possible corrections for
#' minimum percentage of observations and target rate.
#' 
#' The corrections of isotonic regression results present an important step in credit rating model development.
#' The minimum percentage of observation is capped to minimum 30 observations per bin, while target rate for 
#' binary target is capped to 1 bad case. 
#'
#'@param x Numeric vector to be binned.
#'@param y Numeric target vector (binary or continuous).
#'@param sc Numeric vector with special case elements. Default values are \code{c(NA, NaN, Inf, -Inf)}.
#' Recommendation is to keep the default values always and add new ones if needed. Otherwise, if these values exist
#' in \code{x} and are not defined in the \code{sc} vector, function will report the error.  
#'@param sc.method Define how special cases will be treated, all together or in separate bins.
#' Possible values are \code{"together", "separately"}.
#'@param y.type Type of \code{y}, possible options are \code{"bina"} (binary) and \code{"cont"} (continuous).
#' If default value (\code{NA}) is passed, then algorithm will identify if \code{y} is 0/1 or continuous variable.
#'@param min.pct.obs Minimum percentage of observations per bin. Default is 0.05 or minimum 30 observations.
#'@param min.avg.rate Minimum \code{y} average rate. Default is 0.01 or minimum 1 bad case for y 0/1.
#'@param force.trend If the expected trend should be forced. Possible values: \code{"i"} for
#' increasing trend (\code{y} increases with increase of \code{x}), \code{"d"} for decreasing trend 
#' (\code{y} decreases with decrease of \code{x}). Default value is \code{NA}.
#' If the default value is passed, then trend will be identified based on the sign of the Spearman correlation 
#' coefficient between \code{x} and \code{y} on complete cases.
#' 
#'@return The command \code{iso.bin} generates a list of two objects. The first object, data frame \code{summary.tbl}
#' presents a summary table of final binning, while \code{x.trans} is a vector of discretized values. 
#' In case of single unique value for \code{x} or \code{y} of complete cases (cases different than special cases), 
#' it will return data frame with info.
#'
#'@examples
#' suppressMessages(library(monobin))
#' data(gcd)
#' age.bin <- iso.bin(x = gcd$age, y = gcd$qual)
#' age.bin[[1]]
#' table(age.bin[[2]])
#' # force increasing trend
#' iso.bin(x = gcd$age, y = gcd$qual, force.trend = "i")[[1]]
#'
#' #stage by stage example
#' #inputs
#' x <- gcd$age		#risk factor
#' y <- gcd$qual	#binary dependent variable
#' min.pct.obs <- 0.05	#minimum percentage of observations per bin
#' min.avg.rate <- 0.01	#minimum percentage of defaults per bin
#' #stage 1: isotonic regression
#' db <- data.frame(x, y)
#' db <- db[order(db$x), ]
#' cc.sign <- sign(cor(db$y, db$x, method = "spearman", use = "complete.obs"))
#' iso.r <- isoreg(x = db$x, y = cc.sign * db$y)
#' db$y.hat <- iso.r$yf
#' db.s0 <- db %>%
#' 	   group_by(bin = y.hat) %>%
#' 	   summarise(no = n(),
#' 			 y.sum = sum(y),
#' 			 y.avg = mean(y),
#' 			 x.avg = mean(x),
#' 			 x.min = min(x),
#' 			 x.max = max(x))
#' db.s0 
#' #stage 2: merging based on minimum percentage of observations
#' db.s1 <- db.s0
#' thr.no <- ceiling(ifelse(nrow(db) * min.pct.obs < 30, 30, nrow(db) * min.pct.obs))
#' thr.no #threshold for minimum number of observations per bin
#' repeat {
#' 		 if	(nrow(db.s1) == 1) {break}
#' 		 values <- db.s1[, "no"]
#' 		 if	(all(values >= thr.no)) {break}
#' 		 gap <- min(which(values < thr.no))
#' 		 if	(gap == nrow(db.s1)) {
#' 			db.s1$bin[(gap - 1):gap] <- db.s1$bin[(gap - 1)]
#' 			} else {
#' 			db.s1$bin[gap:(gap + 1)] <- db.s1$bin[gap + 1]
#' 			}	
#' 		 db.s1 <- db.s1 %>%
#' 			    group_by(bin) %>%
#' 			    mutate(
#' 				y.avg = weighted.mean(y.avg, no),
#' 				x.avg = weighted.mean(x.avg, no)) %>% 
#' 			    summarise(
#' 				no = sum(no),
#' 				y.sum = sum(y.sum),
#' 				y.avg = unique(y.avg),
#' 				x.avg = unique(x.avg),
#' 				x.min = min(x.min),
#' 				x.max = max(x.max))
#' 		} 
#' db.s1
#' #stage 3: merging based on minimum percentage of bad cases
#' db.s2 <- db.s1
#' thr.nb <- ceiling(ifelse(nrow(db) * min.avg.rate < 1, 1, nrow(db) * min.avg.rate))
#' thr.nb #threshold for minimum number of observations per bin
#' #already each bin has more bad cases than selected threshold hence no need for further merging
#' all(db.s2$y.sum > thr.nb)
#' #final result
#' db.s2
#' #result of the iso.bin function (formatting and certain metrics has been added)
#' iso.bin(x = gcd$age, y = gcd$qual)[[1]]
#'
#'@importFrom stats cor isoreg pt sd weighted.mean
#'@importFrom Hmisc cut2
#'@import dplyr
#'@export
iso.bin <- function(x, y, sc = c(NA, NaN, Inf, -Inf), sc.method = "together", y.type = NA, 
			 min.pct.obs = 0.05, min.avg.rate = 0.01, force.trend = NA) {
	ops <- options(scipen = 20)
	on.exit(options(ops)) 

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
	ds <- woe.calc(tbl = ds, y.check = y.check)
	sc.u <- unique(sc)
	sc.g <- ds$bin[ds$type%in%"special cases"]
	x.mins <- ds$x.min[!ds$bin%in%sc.u & !ds$bin%in%"SC"]
	x.maxs <- ds$x.max[!ds$bin%in%sc.u & !ds$bin%in%"SC"]
	x.trans <- slice.variable(x.orig = d$x, x.lb = x.mins, x.ub = x.maxs, 
					  sc.u = sc.u, sc.g = sc.g) 
return(list(summary.tbl = ds, x.trans = x.trans))
}

#iso summary
iso.summary <- function(tbl, bin) {
	tbl %>%
	group_by_at(bin) %>%
	summarise(
		no = n(),
		y.sum = sum(y),
		y.avg = mean(y),
		x.avg = mean(x),
		x.min = min(x),
		x.max = max(x))
}
#iso binning plus correction for min obs and rate
iso <- function(tbl.sc, tbl.cc, method, min.obs, min.rate, y.check, force.trend) {	
	#special cases
	if	(nrow(tbl.sc) > 0) {	
		if	(method == "together") {
			tbl.sc$bin <- "SC"
			} else {
			tbl.sc$bin <- as.character(tbl.sc$x)
			}
		tbl.sc.s <- iso.summary(tbl = tbl.sc, bin = "bin")
		tbl.sc.s$type <- "special cases"
		} else {
		tbl.sc.s <- data.frame()
		}
	#complete cases
	tbl.cc <- tbl.cc[order(tbl.cc$x), ]
	if	(is.na(force.trend)) {
		cor.coef <- cor(tbl.cc$y, tbl.cc$x, method = "spearman", use = "complete.obs")
		cc.sign <- sign(cor.coef)			
		} else {
		cc.sign <- ifelse(force.trend == "i", 1, -1)
		}
	iso.r <- isoreg(x = tbl.cc$x, y = cc.sign * tbl.cc$y)
	tbl.cc$y.hat <- iso.r$yf
	tbl.cc.s <- iso.summary(tbl = tbl.cc, bin = "y.hat")
	tbl.cc.s <- tbl.cc.s[order(cc.sign * tbl.cc.s$y.hat), ]
	tbl.cc.s <- cbind.data.frame(bin = 1:nrow(tbl.cc.s), tbl.cc.s[, -1])
	tbl.cc.s <- tbl.correction.00(tbl = tbl.cc.s, mno = min.obs, mrate = min.rate,
					      what = "obs", y.check = y.check)
	tbl.cc.s <- tbl.correction.00(tbl = tbl.cc.s, mno = min.obs, mrate = min.rate,
					      what = "rate", y.check = y.check)
	tbl.cc.s <- tbl.cc.s[order(tbl.cc.s$x.avg), ]
	tbl.cc.s$bin <- format.bin(x.lb = tbl.cc.s$x.min, x.ub = tbl.cc.s$x.max)
	tbl.cc.s$type <- "complete cases"
	tbl.s <- bind_rows(tbl.sc.s, tbl.cc.s)
return(as.data.frame(tbl.s))
}
#correction for num of obs and min rate
tbl.correction.00 <- function(tbl, mno, mrate, what, y.check) {
	if	(what == "obs") {
		cn <- "no"; thr <- mno
		} else {
		if	(y.check == "bina") {cn <- "y.sum"} else {cn <- "y.avg"}
		thr <- mrate}
	repeat {
		 if	(nrow(tbl) == 1) {break}
		 values <- tbl[, cn]
		 if	(all(values >= thr)) {break}
		 gap <- min(which(values < thr))
		 if	(gap == nrow(tbl)) {
			tbl$bin[(gap - 1):gap] <- tbl$bin[(gap - 1)]
			} else {
			tbl$bin[gap:(gap + 1)] <- tbl$bin[gap + 1]
			}	
		 tbl <- tbl %>%
			  group_by(bin) %>%
			  mutate(
				y.avg = weighted.mean(y.avg, no),
				x.avg = weighted.mean(x.avg, no)) %>% 
			  summarise(
				no = sum(no),
				y.sum = sum(y.sum),
				y.avg = unique(y.avg),
				x.avg = unique(x.avg),
				x.min = min(x.min),
				x.max = max(x.max))
		} 
return(tbl)
}
#woe calculation
woe.calc <- function(tbl, y.check) {
	if	(y.check == "bina") {
		res <- tbl %>% 
			 mutate(
				so = sum(no),
				sg = sum(no) - sum(y.sum),
				sb = sum(y.sum), 
				dist.g = (no - y.sum) / sg,
				dist.b = y.sum / sb,
				woe = log(dist.g / dist.b),
				iv.b = (dist.g - dist.b) * woe)
		} else {
		res <- tbl %>% 
			 mutate(
				so = sum(no),
				sy = sum(y.sum),
				pct.obs = no / so,
				pct.y.sum = y.sum / sy,
				woe = log(pct.y.sum / pct.obs),
				iv.b = (pct.y.sum - pct.obs) * woe)
		}
return(res)
}



