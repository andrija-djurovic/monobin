#' Monotonic binning based on maximum cumulative target rate (MAPA)
#'
#' \code{cum.bin} implements monotonic binning based on maximum cumulative target rate. 
#' This algorithm is known as MAPA (Monotone Adjacent Pooling Algorithm).
#'@param x Numeric vector to be binned.
#'@param y Numeric target vector (binary or continuous).
#'@param sc Numeric vector with special case elements. Default values are \code{c(NA, NaN, Inf)}.
#' Recommendation is to keep the default values always and add new ones if needed. Otherwise, if these values exist
#' in \code{x} and are not defined in the \code{sc} vector, function will report the error.  
#'@param sc.method Define how special cases will be treated, all together or in separate bins.
#' Possible values are \code{"together", "separately"}.
#'@param g Number of starting groups. Default is 15.
#'@param y.type Type of \code{y}, possible options are \code{"bina"} (binary) and  \code{"cont"} (continuous).
#' If default value (\code{NA}) is passed, then algorithm will identify if \code{y} is 0/1 or continuous variable.
#'@param force.trend If the expected trend should be forced. Possible values: \code{"i"} for
#' increasing trend (\code{y} increases with increase of \code{x}), \code{"d"} for decreasing trend 
#' (\code{y} decreases with decrease of \code{x}). Default value is \code{NA}.
#' If the default value is passed, then trend will be identified based on the sign of the Spearman correlation 
#' coefficient between \code{x} and \code{y} on complete cases.
#'
#'@return The command \code{cum.bin} generates a list of two objects. The first object, data frame \code{summary.tbl}
#' presents a summary table of final binning, while \code{x.trans} is a vector of discretized values. 
#' In case of single unique value for \code{x} or \code{y} in complete cases (cases different than special cases), 
#' it will return data frame with info.
#'
#'@examples
#' suppressMessages(library(monobin))
#' data(gcd)
#' amount.bin <- cum.bin(x = gcd$amount, y = gcd$qual)
#' amount.bin[[1]]
#' gcd$amount.bin <- amount.bin[[2]]
#' gcd %>% group_by(amount.bin) %>% summarise(n = n(), y.avg = mean(qual))
#' #increase default number of groups (g = 20)
#' amount.bin.1 <- cum.bin(x = gcd$amount, y = gcd$qual, g = 20)
#' amount.bin.1[[1]]
#' #force trend to decreasing
#' cum.bin(x = gcd$amount, y = gcd$qual, g = 20, force.trend = "d")[[1]]
#'
#'@importFrom stats cor isoreg pt sd weighted.mean
#'@importFrom Hmisc cut2
#'@importFrom dplyr group_by summarise ungroup mutate
#'@export

#iso.binning
cum.bin <- function(x, y, sc = c(NA, NaN, Inf), sc.method = "together", g = 15, y.type = NA,
			  force.trend = NA) {
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

	#special cases
	if	(nrow(d.sc) > 0) {	
		if	(sc.method[1] == "together") {
			d.sc$bin <- "SC"
			} else {
			d.sc$bin <- as.character(d.sc$x)
			}
		d.sc.s <- iso.summary(tbl = d.sc, bin = "bin")
		d.sc.s$type <- "special cases"
		} else {
		d.sc.s <- data.frame()
		}
	#complete cases
	tbl <- cbind.data.frame(d.cc, bin = cut2(d.cc$x, g = g))
	d.cc.s <- cum.bin.aux(tbl = tbl, force.trend = force.trend, y.check = y.check)
	ds <- as.data.frame(bind_rows(d.sc.s, d.cc.s))
	ds <- woe.calc(tbl = ds, y.check = y.check)
	sc.u <- unique(sc)
	sc.g <- ds$bin[ds$type%in%"special cases"]
	x.mins <- ds$x.min[!ds$bin%in%sc.u & !ds$bin%in%"SC"]
	x.maxs <- ds$x.max[!ds$bin%in%sc.u & !ds$bin%in%"SC"]
	x.trans <- slice.variable(x.orig = d$x, x.lb = x.mins, x.ub = x.maxs, 
					  sc.u = sc.u, sc.g = sc.g) 
return(list(summary.tbl = ds, x.trans = x.trans))
}
#cumulative binning function
cum.bin.aux <- function(tbl, force.trend, y.check){
	tbl <- tbl[order(tbl$x), ]
	if	(is.na(force.trend)) {
		sort.d <- sign(cor(tbl$y, tbl$x, method = "spearman", use = "complete.obs"))
		} else {
		sort.d <- ifelse(force.trend == "i", 1, -1)
		}
	tbl <- iso.summary(tbl = tbl, "bin")
	if	(nrow(tbl) == 1) {
		tbl$type <- "complete cases"
		return(as.data.frame(tbl))
		}
	tbl <- tbl[order(tbl$bin, decreasing = ifelse(sort.d == 1, TRUE, FALSE)), ]
	if	(all(diff(tbl$y.avg) > 0)) {
		return(tbl)
		}
	tbl.i <- tbl
	cp <- c()
	repeat {
		 cs <- cumsum(tbl.i$y.sum) / cumsum(tbl.i$no)
		 indx <- which(cs == max(cs))[1]
		 if	(indx == nrow(tbl.i)) {cp <- c(cp, indx); break}
		 cp <- c(cp, indx)
		 tbl.i <- tbl.i[-(1:indx), ]
		 }
	if	(length(cp) == 1) {
		res <- cbind.data.frame(tbl, mod = 1, stringsAsFactors = FALSE)
		} else {
		mod <- rep(cumsum(cp), times = cp)
		res <- cbind.data.frame(tbl, mod)
		}	
	res.s <- res %>%
		   group_by(bin = mod) %>%
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
	res.s <- res.s[order(res.s$x.avg), ]
	res.s$bin <- format.bin(x.lb = res.s$x.min, x.ub = res.s$x.max)
	res.s$type <- "complete cases"	
return(as.data.frame(res.s))
}


