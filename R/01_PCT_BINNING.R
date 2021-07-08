#' Monotonic binning based on percentiles
#'
#' \code{pct.bin} implements percentile-based monotonic binning by the iterative discretization. 
#'@param x Numeric vector to be binned.
#'@param y Numeric target vector (binary or continuous).
#'@param sc Numeric vector with special case elements. Default values are \code{c(NA, NaN, Inf)}.
#' Recommendation is to keep the default values always and add new ones if needed. Otherwise, if these values exist
#' in \code{x} and are not defined in the \code{sc} list some statistics cannot be calculated properly. 
#'@param sc.method Define how special cases will be treated, all together or in separate bins.
#' Possible values are \code{"together", "separately"}.
#'@param g Number of starting groups. Default is 15.
#'@param y.type Type of \code{y}, possible options are \code{"bina"} (binary) and \code{"cont"} (continuous). 
#' If default value is passed, then algorithm will identify if y is 0/1 or continuous variable.
#'@param woe.trend Applied only for a continuous target (\code{y}) as weights of evidence (WoE) trend check. Default is TRUE.
#'@param force.trend If the expected trend should be forced. Possible values: \code{"i"} for
#' increasing trend (\code{y} increases with increase of \code{x}), \code{"d"} for decreasing trend 
#' (\code{y} decreases with decrease of \code{x}). Default value is \code{NA}.
#' If the default value is passed, algorithm will stop if perfect negative or positive correlation (Spearman) is achieved
#' between average \code{y} and average \code{x} per bin. Otherwise, it will stop only if the forced trend is achieved. 
#' 
#'@return The command \code{pct.bin} generates a list of two objects. The first object, data frame \code{summary.tbl}
#' presents a summary table of final binning, while \code{x.trans} is a vector of discretized values. 
#' In case of single unique value for \code{x} or \code{y} of complete cases (cases different than special cases), 
#' it will return data frame with info.
#'
#'@examples
#' suppressMessages(library(monobin))
#' data(gcd)
#' #binary target
#' mat.bin <- pct.bin(x = gcd$maturity, y = gcd$qual)
#' mat.bin[[1]]
#' table(mat.bin[[2]])
#' #continuous target, separate groups for special cases
#' set.seed(123)
#' gcd$age.d <- gcd$age
#' gcd$age.d[sample(1:nrow(gcd), 10)] <- NA
#' gcd$age.d[sample(1:nrow(gcd), 3)] <- 9999999999
#' age.d.bin <- pct.bin(x = gcd$age.d, 
#' 			   	y = gcd$qual, 
#' 			   	sc = c(NA, NaN, Inf, 9999999999), 
#' 			  	sc.method = "separately",
#' 			   	force.trend = "d")
#' age.d.bin[[1]]
#' gcd$age.d.bin <- age.d.bin[[2]]
#' gcd %>% group_by(age.d.bin) %>% summarise(n = n(), y.avg = mean(qual))
#'
#'@importFrom stats cor isoreg pt sd weighted.mean
#'@importFrom Hmisc cut2
#'@import dplyr
#'@export

#pct.binning
pct.bin <- function(x, y, sc = c(NA, NaN, Inf), sc.method = "together", g = 15, 
			  y.type = NA, woe.trend = TRUE, force.trend = NA) {
	ops <- options(scipen = 20)
	on.exit(options(ops)) 
	
	if	(!is.logical(woe.trend)) {
		stop("woe.trend has to be logical: TRUE or FALSE")
		}
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

	if	(y.check == "bina") {
		ds <- pct.bin.bina(tbl.sc = d.sc, tbl.cc = d.cc, 
					 method = sc.method, g = g, force.trend = force.trend)
		} else {
		ds <- pct.bin.cont(tbl.sc = d.sc, tbl.cc = d.cc, 
					 method = sc.method, g = g, woe.trend = woe.trend,
					 force.trend = force.trend)
		}
	sc.u <- unique(sc)
	sc.g <- ds$bin[ds$type%in%"special cases"]
	x.mins <- ds$x.min[!ds$bin%in%sc.u & !ds$bin%in%"SC"]
	x.maxs <- ds$x.max[!ds$bin%in%sc.u & !ds$bin%in%"SC"]
	x.trans <- slice.variable(x.orig = d$x, x.lb = x.mins, x.ub = x.maxs, 
					  sc.u = sc.u, sc.g = sc.g) 	
return(list(summary.tbl = ds, x.trans = x.trans))
}

#formatting bins
format.bin <- function(x.lb, x.ub) {
	x.lb.lag <- c(x.lb[-1], Inf)
	bin.n <- sprintf("%02d", 1:length(x.lb))  
	bin.f <- ifelse(abs(x.lb - x.ub) < 1e-8, 
			    paste0(bin.n, " [", round(x.lb, 4), "]"), 
			    paste0(bin.n, " [", round(x.lb, 4), ",", round(x.lb.lag, 4), ")"))
return(bin.f)
}

#summary binary 
tbl.summary.bina <- function(tbl, g.tot, b.tot) {
	tbl %>% 
	group_by(bin) %>%
	summarise(
		no = n(),
		y.sum = sum(y),
		y.avg = mean(y),
		x.avg = mean(x),
		x.min = min(x),
		x.max = max(x)) %>%
	ungroup() %>%
	mutate(
		so = g.tot + b.tot,
		sg = g.tot,
		sb = b.tot, 
		dist.g = (no - y.sum) / sg,
		dist.b = y.sum / sb,
		woe = log(dist.g / dist.b),
		iv.b = (dist.g - dist.b) * woe)
}

#binning binary
pct.bin.bina <- function(tbl.sc, tbl.cc, method, g, force.trend) {	
	y.tot <- nrow(tbl.sc) + nrow(tbl.cc)
	b.tot <- sum(tbl.sc$y) + sum(tbl.cc$y)
	g.tot <- y.tot - b.tot
	#special cases
	if	(nrow(tbl.sc) > 0) {	
		if	(method == "together") {
			tbl.sc$bin <- "SC"
			} else {
			tbl.sc$bin <- as.character(tbl.sc$x)
			}
		tbl.sc.s <- tbl.summary.bina(tbl = tbl.sc, g.tot = g.tot, b.tot = b.tot)
		tbl.sc.s$type <- "special cases"
		} else {
		tbl.sc.s <- data.frame()
		}
	#complete cases
	if	(!is.na(force.trend) & force.trend == "i") {
		cond.exp <- "isTRUE(all.equal(cor.coef, 1)) | g == 1"
		}
	if	(!is.na(force.trend) & force.trend == "d") {
		cond.exp <- "isTRUE(all.equal(cor.coef, -1)) | g == 1"
		}
	if	(is.na(force.trend)) {
		cond.exp <- "isTRUE(all.equal(cor.coef, 1)) | isTRUE(all.equal(cor.coef, -1)) | g == 1"
		}	
	repeat {
		 tbl.cc$bin = cut2(tbl.cc$x, g = g)
		 tbl.cc.s <- tbl.summary.bina(tbl = tbl.cc, g.tot = g.tot, b.tot = b.tot)
		 if	(nrow(tbl.cc.s) == 1) {break}
		 cor.coef <- cor(tbl.cc.s$y.avg, tbl.cc.s$x.avg, method = "spearman", use = "complete.obs")
		 if	(eval(parse(text = cond.exp))) {break}
		 g <- g - 1
		 }
	tbl.cc.s$bin <- format.bin(x.lb = tbl.cc.s$x.min, x.ub = tbl.cc.s$x.max)
	tbl.cc.s$type <- "complete cases"
	tbl.s <- bind_rows(tbl.sc.s, tbl.cc.s)
return(as.data.frame(tbl.s))
}

#summary continuous
tbl.summary.cont <- function(tbl, n.tot, y.tot) {
	tbl %>% 
	group_by(bin) %>%
	summarise(
		no = n(),
		y.sum = sum(y),
		y.avg = mean(y),
		x.avg = mean(x),
		x.min = min(x),
		x.max = max(x)) %>%
	ungroup() %>%
	mutate(
		so = n.tot,
		sy = y.tot,
		pct.obs = no / n.tot,
		pct.y.sum = y.sum / y.tot,
		woe = log(pct.y.sum / pct.obs),
		iv.b = (pct.y.sum - pct.obs) * woe)
}

#binning continuous
pct.bin.cont <- function(tbl.sc, tbl.cc, method, g, woe.trend, force.trend) {
	n.tot <- nrow(tbl.sc) + nrow(tbl.cc)	
	y.tot <- sum(tbl.sc$y) + sum(tbl.cc$y) 
	#special cases
	if	(nrow(tbl.sc) > 0) {	
		if	(method == "together") {
			tbl.sc$bin <- "SC"
			} else {
			tbl.sc$bin <- as.character(tbl.sc$x)
			}
		tbl.sc.s <- tbl.summary.cont(tbl = tbl.sc, n.tot = n.tot, y.tot = y.tot)
		tbl.sc.s$type <- "special cases"
		} else {
		tbl.sc.s <- data.frame()
		}
	#complete cases
	if	(!is.na(force.trend) & force.trend == "i") {
		cond.exp.1 <- "isTRUE(all.equal(cor.coef, 1)) | g == 1"
		cond.exp.2 <- "all(diff(tbl.cc.s$woe) > 0)"
		}
	if	(!is.na(force.trend) & force.trend == "d") {
		cond.exp.1 <- "isTRUE(all.equal(cor.coef, -1)) | g == 1"
		cond.exp.2 <- "all(diff(tbl.cc.s$woe) < 0)"
		}
	if	(is.na(force.trend)) {
		cond.exp.1 <- "isTRUE(all.equal(cor.coef, 1)) | isTRUE(all.equal(cor.coef, -1)) | g == 1"
		cond.exp.2 <- "all(diff(tbl.cc.s$woe) > 0) | all(diff(tbl.cc.s$woe) < 0)"		
		}
	repeat {
		 tbl.cc$bin = cut2(tbl.cc$x, g = g)
		 tbl.cc.s <- tbl.summary.cont(tbl = tbl.cc, n.tot = n.tot, y.tot = y.tot)
		 if	(nrow(tbl.cc.s) == 1) {break}
		 if	(woe.trend) {
			monocheck <- eval(parse(text = cond.exp.2))
			} else {
		 	cor.coef <- cor(tbl.cc.s$y.avg, tbl.cc.s$x.avg, method = "spearman", use = "complete.obs")
			monocheck <- eval(parse(text = cond.exp.1))
			}
		 if(monocheck | g == 1) {break}
		 g <- g - 1
		 }
	tbl.cc.s$bin <- format.bin(tbl.cc.s$x.min, tbl.cc.s$x.max)
	tbl.cc.s$type <- "complete cases"
	tbl.s <- bind_rows(tbl.sc.s, tbl.cc.s)
return(as.data.frame(tbl.s))
}

#transform x vector
slice.variable <- function(x.orig, x.lb, x.ub, sc.u, sc.g) {
	lx <- length(x.orig)
	lg <- length(x.lb)
	x.trans <- x.orig
	if	(length(sc.g) > 0) {
			if	("SC"%in%sc.g) {
				x.trans[x.trans%in%sc.u] <- "SC"
				}
		}
	x.lb.lag <- c(x.lb[-1], Inf)
	for	(i in 1:lg) {
		x.lb.l <- x.lb[i]
		x.lb.lag.l <- x.lb.lag[i]
		x.ub.l <- x.ub[i]
		bin.n <- sprintf("%02d", i)
		bin.f <- ifelse(x.lb.l == x.ub.l, 
				    paste0(bin.n, " [", x.lb.l, "]"), 
				    paste0(bin.n, " [", x.lb.l, ",", x.lb.lag.l, ")"))
		rep.indx <- which(!x.orig%in%sc.u & !x.orig%in%sc.g & x.orig >= x.lb.l & x.orig <= x.ub.l)
		x.trans[rep.indx] <- bin.f
		}	
return(x.trans)
}
          

