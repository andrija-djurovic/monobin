#' Four-stage monotonic binning procedure with statistical test correction
#'
#' \code{sts.bin} implements extension of the three-stage monotonic binning procedure (\code{\link{iso.bin}})
#' with final step of iterative merging of adjacent bins based on
#' statistical test.
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
#'@param p.val Threshold for p-value of statistical test. Default is 0.05. For binary target test of two proportion
#' is applied, while for continuous two samples independent t-test.
#'@param force.trend If the expected trend should be forced. Possible values: \code{"i"} for
#' increasing trend (\code{y} increases with increase of \code{x}), \code{"d"} for decreasing trend 
#' (\code{y} decreases with decrease of \code{x}). Default value is \code{NA}.
#' If the default value is passed, then trend will be identified based on the sign of the Spearman correlation 
#' coefficient between \code{x} and \code{y} on complete cases. 
#'
#'@return The command \code{sts.bin} generates a list of two objects. The first object, data frame \code{summary.tbl}
#' presents a summary table of final binning, while \code{x.trans} is a vector of discretized values.
#' In case of single unique value for \code{x} or \code{y} of complete cases (cases different than special cases), 
#' it will return data frame with info. 
#'
#'@examples
#' suppressMessages(library(monobin))
#' data(gcd)
#' #binary target
#' maturity.bin <- sts.bin(x = gcd$maturity, y = gcd$qual)
#' maturity.bin[[1]]
#' tapply(gcd$qual, maturity.bin[[2]], function(x) c(length(x), sum(x), mean(x)))
#' prop.test(x = c(sum(gcd$qual[maturity.bin[[2]]%in%"01 [4,8)"]), 
#'		       sum(gcd$qual[maturity.bin[[2]]%in%"02 [8,16)"])), 
#'	       n = c(length(gcd$qual[maturity.bin[[2]]%in%"01 [4,8)"]),
#'		       length(gcd$qual[maturity.bin[[2]]%in%"02 [8,16)"])), 
#'	       alternative = "less", 
#'	       correct = FALSE)$p.value
#' #continuous target
#' age.bin <- sts.bin(x = gcd$age, y = gcd$qual, y.type = "cont")
#' age.bin[[1]]
#' t.test(x = gcd$qual[age.bin[[2]]%in%"01 [19,26)"], 
#'	    y = gcd$qual[age.bin[[2]]%in%"02 [26,35)"],
#'	    alternative = "greater")$p.value
#' 
#'@importFrom stats cor isoreg pt sd weighted.mean
#'@importFrom Hmisc cut2
#'@import dplyr
#'@export

#statistical test binning
sts.bin <- function(x, y, sc = c(NA, NaN, Inf), sc.method = "together", y.type = NA, 
			 min.pct.obs = 0.05, min.avg.rate = 0.01, p.val = 0.05, force.trend = NA) {
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
	ds.sc <- ds[ds$type%in%"special cases", ]
	ds.cc <- ds[ds$type%in%"complete cases", ]

	if	(y.check == "bina") {
		ds.cc <- t2p.merge(tbl = ds.cc, sig = p.val)
		} else {
		ds.cc$y.sd <- add.sd(tbl = ds.cc, x = d$x, y = d$y, sc = sc)
		ds.cc <- ttg.merge(tbl = ds.cc, sig = p.val, ds = d.cc)
		}
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
#add standard deviation for ds table
add.sd <- function(tbl, x, y, sc) {
	ngr <- nrow(tbl)
	sd.gr <- rep(NA, ngr)	
	for	(i in 1:ngr) {
		x.lb.l <- tbl$x.min[i]
		x.ub.l <- tbl$x.max[i]
		y.l <- y[!x%in%sc & !x%in%"SC" & x >= x.lb.l & x <= x.ub.l]
		sd.gr[i] <- sd(y.l)
		}
return(sd.gr)
}
#bin merging based on test of 2 proportions
t2p.merge <- function(tbl, sig) {
	if	(nrow(tbl) <= 1) {
		tbl$p.val <- 99
		return(tbl)
		}
	sts <- ifelse(sign(cor(tbl$y.avg, tbl$x.avg, method = "spearman")) == 1, "less", "greater")
	test.exp <- "prop.test(x = c(y.sum1, y.sum2), 
				     n = c(no1, no2), alternative = sts, 
				     correct = FALSE)$p.value"
	tbl$p.val <- NA
	for	(i in 2:nrow(tbl)) {
		y.sum1 <- tbl$y.sum[i - 1]
		y.sum2 <- tbl$y.sum[i]
		no1 <- tbl$no[i - 1]
		no2 <- tbl$no[i]
		p.val.l <- eval(parse(text = test.exp))
		tbl$p.val[i] <- p.val.l
		}
	tbl$mod <- 1:nrow(tbl)
	repeat {
		#criteria to exit the loop
		if	(all(tbl$p.val[!is.na(tbl$p.val)] < sig)) {break}
		if	(nrow(tbl) == 1) {break}
		p.max <- which.max(tbl$p.val)[1]
		tbl$mod[(p.max -1):p.max] <- tbl$mod[p.max]
		bm <- c(as.character(tbl$bin[p.max]), as.character(tbl$bin[p.max - 1]))
		#find previous and next bins after merging bins in previous step
		if	((p.max - 2) < 1) {
			bp <- NULL
			} else {
			bp <- as.character(tbl$bin[p.max - 2])
			}
		if	((p.max + 1) > nrow(tbl)) {
			bn <- NULL
			} else {
			bn <- as.character(tbl$bin[p.max + 1])
			}
		#recalculate p-values for new group and its neighbors
		if	(!is.null(bp)) {
			y.sum1 <- tbl$y.sum[tbl$bin%in%bp]; no1 <- tbl$no[tbl$bin%in%bp]
			y.sum2 <- sum(tbl$y.sum[tbl$bin%in%bm]); no2 <- sum(tbl$no[tbl$bin%in%bm])
			p.val <- eval(parse(text = test.exp))
			tbl$p.val[tbl$mod%in%tbl$mod[p.max]] <- p.val
			} else {
			tbl$p.val[p.max - 1] <- tbl$p.val[p.max]
			}
		if	(!is.null(bn)) {
			y.sum1 <- sum(tbl$y.sum[tbl$bin%in%bm]); no1 <- sum(tbl$no[tbl$bin%in%bm])
			y.sum2 <- tbl$y.sum[tbl$bin%in%bn]; no2 <- tbl$no[tbl$bin%in%bn]
			p.val <- eval(parse(text = test.exp))
			tbl$p.val[tbl$mod%in%tbl$mod[p.max + 1]] <- p.val
			} else {
			tbl$p.val[p.max - 1] <- tbl$p.val[p.max]
			}
		#condition for only 2 remaining groups with insignificant split
		if	(is.null(bp) & is.null(bn)) {
			tbl$p.val[(p.max -1):p.max] <- 1
			}
		#summarize data based on group correction
		 tbl <- tbl %>%
			  group_by(mod) %>%
			  mutate(
				y.avg = weighted.mean(y.avg, no),
				x.avg = weighted.mean(x.avg, no)) %>% 
			  summarise(
				bin = paste0(bin, collapse = " + "),
				no = sum(no),
				y.sum = sum(y.sum),
				y.avg = unique(y.avg),
				x.avg = unique(x.avg),
				x.min = min(x.min),
				x.max = max(x.max),
				p.val = unique(p.val))
		tbl$p.val[1] <- NA
		}
	res <- as.data.frame(tbl)[, -1]
	res <- res[order(res$x.avg), ]
	res$bin <- format.bin(x.lb = res$x.min, x.ub = res$x.max)
	res$type <- "complete cases"
return(res)
}
#bin merging based on t test
t.test.g <- function(x1, x2, sd1, sd2, no1, no2, alternative) {
	std1 <- sd1/sqrt(no1)
	std2 <- sd2/sqrt(no2)
	std <- sqrt(std1^2 + std2^2)
	t.stat <- (x1 - x2) / std
	df <-  std^4/(std1^4 / (no1 - 1) + std2^4 / (no2 - 1))
	if	(alternative == "less") {
		 p.val <- pt(t.stat, df)
		 } else {
		 p.val <- pt(t.stat, df, lower.tail = FALSE)
		 }
return(p.val)
}
ttg.merge <- function(tbl, sig, ds) {
	if	(nrow(tbl) == 1) {
		tbl$p.val <- 99
		return(tbl)
		}
	sts <- ifelse(sign(cor(tbl$y.avg, tbl$x.avg, method = "spearman")) == 1, "less", "greater")

	test.exp.init <- "t.test.g(x1 = x1, x2 = x2, sd1 = sd1, sd2 = sd2,
				    no1 = no1, no2 = no2, alternative = sts)"
	test.exp.iter <- "t.test(x = y1, y = y2, alternative = sts)"

	tbl$p.val <- NA
	for	(i in 2:nrow(tbl)) {
		x1 <- tbl$y.avg[i - 1]
		x2 <- tbl$y.avg[i]
		sd1 <- tbl$y.sd[i - 1]
		sd2 <- tbl$y.sd[i]
		no1 <- tbl$no[i - 1]
		no2 <- tbl$no[i]
		p.val.l <- eval(parse(text = test.exp.init))
		tbl$p.val[i] <- p.val.l
		}
	tbl$mod <- 1:nrow(tbl)
	repeat {
		#criteria to exit the loop
		if	(all(tbl$p.val[!is.na(tbl$p.val)] < sig)) {break}
		if	(nrow(tbl) == 1) {break}
		p.max <- which.max(tbl$p.val)[1]
		tbl$mod[(p.max -1):p.max] <- tbl$mod[p.max]
		bm <- c(as.character(tbl$bin[p.max]), as.character(tbl$bin[p.max - 1]))
		#find previous and next bins after merging 2 in previous step
		if	((p.max - 2) < 1) {
			bp <- NULL
			} else {
			bp <- as.character(tbl$bin[p.max - 2])
			}
		if	((p.max + 1) > nrow(tbl)) {
			bn <- NULL
			} else {
			bn <- as.character(tbl$bin[p.max + 1])
			}
		#recalculate p-values for new group and its neighbors
		if	(!is.null(bp)) {
			y1 <- ds$y[ds$x >= tbl$x.min[tbl$bin%in%bp] & 
				     ds$x <= tbl$x.max[tbl$bin%in%bp]]
			y2 <- ds$y[ds$x >= min(tbl$x.min[tbl$bin%in%bm]) & 
				     ds$x <= max(tbl$x.max[tbl$bin%in%bm])]
			p.val <- eval(parse(text = test.exp.iter))$p.value
			tbl$p.val[tbl$mod%in%tbl$mod[p.max]] <- p.val
			} else {
			tbl$p.val[p.max - 1] <- tbl$p.val[p.max]
			}
		if	(!is.null(bn)) {
			y1 <- ds$y[ds$x >= min(tbl$x.min[tbl$bin%in%bm]) & 
				     ds$x <= max(tbl$x.max[tbl$bin%in%bm])]
			y2 <- ds$y[ds$x >= tbl$x.min[tbl$bin%in%bn] & 
				     ds$x <= tbl$x.max[tbl$bin%in%bn]]
			p.val <- eval(parse(text = test.exp.iter))$p.value
			tbl$p.val[tbl$mod%in%tbl$mod[p.max + 1]] <- p.val
			} else {
			tbl$p.val[p.max - 1] <- tbl$p.val[p.max]
			}
		#condition for only 2 remaining groups with insignificant split
		if	(is.null(bp) & is.null(bn)) {
			tbl$p.val[(p.max -1):p.max] <- 1
			}
		#summarize data based on group correction
		 tbl <- tbl %>%
			  group_by(mod) %>%
			  mutate(
				y.avg = weighted.mean(y.avg, no),
				x.avg = weighted.mean(x.avg, no)) %>% 
			  summarise(
				bin = paste0(bin, collapse = " + "),
				no = sum(no),
				y.sum = sum(y.sum),
				y.avg = unique(y.avg),
				x.avg = unique(x.avg),
				x.min = min(x.min),
				x.max = max(x.max),
				p.val = unique(p.val))
		tbl$p.val[1] <- NA
		}
	res <- as.data.frame(tbl)[, -1]
	res <- res[order(res$x.avg), ]
	res$bin <- format.bin(x.lb = res$x.min, x.ub = res$x.max)
	res$type <- "complete cases"
return(res)
}




