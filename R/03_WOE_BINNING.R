#' Four-stage monotonic binning procedure with WoE threshold
#'
#' \code{woe.bin} implements extension of the three-stage monotonic binning procedure (\code{\link{iso.bin}})
#' with weights of evidence (WoE) threshold. 
#' The first stage is isotonic regression used to achieve the monotonicity. The next two stages are possible corrections for
#' minimum percentage of observations and target rate, while the last stage is iterative merging of  
#' bins until WoE threshold is exceeded.
#'
#'@seealso \code{\link{iso.bin}} for three-stage monotonic binning procedure.
#'
#'@param x Numeric vector to be binned.
#'@param y Numeric target vector (binary or continuous).
#'@param sc Numeric vector with special case elements. Default values are \code{c(NA, NaN, Inf)}.
#' Recommendation is to keep the default values always and add new ones if needed. Otherwise, if these values exist
#' in \code{x} and are not defined in the \code{sc} vector, function will report the error. 
#'@param sc.method Define how special cases will be treated, all together or in separate bins.
#' Possible values are \code{"together", "separately"}.
#'@param y.type Type of \code{y}, possible options are \code{"bina"} (binary) and \code{"cont"} (continuous).
#' If default value (\code{NA}) is passed, then algorithm will identify if \code{y} is 0/1 or continuous variable.
#'@param min.pct.obs Minimum percentage of observations per bin. Default is 0.05 or minimum 30 observations.
#'@param min.avg.rate Minimum \code{y} average rate. Default is 0.01 or minimum 1 bad case for y 0/1.
#'@param woe.gap Minimum WoE gap between bins. Default is 0.1.
#'@param force.trend If the expected trend should be forced. Possible values: \code{"i"} for
#' increasing trend (\code{y} increases with increase of \code{x}), \code{"d"} for decreasing trend 
#' (\code{y} decreases with decrease of \code{x}). Default value is \code{NA}.
#' If the default value is passed, then trend will be identified based on the sign of the Spearman correlation 
#' coefficient between \code{x} and \code{y} on complete cases.
#' 
#'@return The command \code{woe.bin} generates a list of two objects. The first object, data frame \code{summary.tbl}
#' presents a summary table of final binning, while \code{x.trans} is a vector of discretized values. 
#' In case of single unique value for \code{x} or \code{y} of complete cases (cases different than special cases), 
#' it will return data frame with info.
#'
#'@examples
#' suppressMessages(library(monobin))
#' data(gcd)
#' amount.bin <- woe.bin(x = gcd$amount, y = gcd$qual)
#' amount.bin[[1]]
#' diff(amount.bin[[1]]$woe)
#' tapply(gcd$amount, amount.bin[[2]], function(x) c(length(x), mean(x)))
#' woe.bin(x = gcd$maturity, y = gcd$qual)[[1]]
#' woe.bin(x = gcd$maturity, y = gcd$qual, woe.gap = 0.5)[[1]]
#'
#'@importFrom stats cor isoreg pt sd weighted.mean
#'@importFrom Hmisc cut2
#'@import dplyr
#'@export

woe.bin <- function(x, y, sc = c(NA, NaN, Inf), sc.method = "together", y.type = NA, 
			 min.pct.obs = 0.05, min.avg.rate = 0.01, woe.gap = 0.1, force.trend = NA) {
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
	ds <- tbl.correction.01(tbl = ds, woe.gap = woe.gap, y.check = y.check)
	sc.u <- unique(sc)
	sc.g <- ds$bin[ds$type%in%"special cases"]
	x.mins <- ds$x.min[!ds$bin%in%sc.u & !ds$bin%in%"SC"]
	x.maxs <- ds$x.max[!ds$bin%in%sc.u & !ds$bin%in%"SC"]
	x.trans <- slice.variable(x.orig = d$x, x.lb = x.mins, x.ub = x.maxs, 
					  sc.u = sc.u, sc.g = sc.g) 
return(list(summary.tbl = ds, x.trans = x.trans))
}
#correction for min woe gap
tbl.correction.01 <- function(tbl, woe.gap, y.check) {
	sc <- tbl[tbl$type%in%"special cases", ]
	cc <- tbl[tbl$type%in%"complete cases", ]
	cc <- cc[order(cc$y.avg), ]
	cc$bin <- sprintf("%02d", 1:nrow(cc))
	woe.gap <- ifelse(y.check == "bina", -woe.gap, woe.gap)
	cond <- ifelse(y.check == "bina", "all(woe.diff < woe.gap)", "all(woe.diff > woe.gap)")
	gap.indx <- ifelse(y.check == "bina", "which.max(woe.diff)[1] + 1", 
				 "which.min(woe.diff)[1] + 1")
	repeat {
		 if	(nrow(cc) == 1) {break}
	 	 woe.diff <- diff(cc$woe)
		 if	(eval(parse(text = cond))) {break}
		 gap <- eval(parse(text = gap.indx))
 		 cc$bin[(gap - 1):gap] <- cc$bin[(gap - 1)]
		 cc <- cc %>%
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
				x.max = max(x.max),
				type = unique(type))
 		 scc <- woe.calc(tbl = bind_rows(sc[, names(cc)], cc), y.check = y.check)
		 cc <- scc[scc$type%in%"complete cases", ]
		}
	cc <- cc[order(cc$x.avg), ]
	cc$bin <- format.bin(x.lb = cc$x.min, x.ub = cc$x.max)
	res <- bind_rows(sc, cc)
return(res)
}



