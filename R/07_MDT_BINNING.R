#' Monotonic binning driven by decision tree
#'
#' \code{mdt.bin} implements monotonic binning driven by decision tree. As a splitting metric for continuous target
#' algorithm uses sum of squared errors, while for the binary target Gini index is used.
#'@param x Numeric vector to be binned.
#'@param y Numeric target vector (binary or continuous).
#'@param g Number of splitting groups for each node. Default is 50.
#'@param sc Numeric vector with special case elements. Default values are \code{c(NA, NaN, Inf, -Inf)}.
#' Recommendation is to keep the default values always and add new ones if needed. Otherwise, if these values exist
#' in \code{x} and are not defined in the \code{sc} vector, function will report the error.  
#'@param sc.method Define how special cases will be treated, all together or in separate bins.
#' Possible values are \code{"together", "separately"}.
#'@param y.type Type of \code{y}, possible options are \code{"bina"} (binary) and  \code{"cont"} (continuous).
#' If default value (\code{NA}) is passed, then algorithm will identify if \code{y} is 0/1 or continuous variable.
#'@param min.pct.obs Minimum percentage of observations per bin. Default is 0.05 or minimum 30 observations.
#'@param min.avg.rate Minimum \code{y} average rate. Default is 0.01 or minimum 1 bad case for y 0/1.
#'@param force.trend If the expected trend should be forced. Possible values: \code{"i"} for
#' increasing trend (\code{y} increases with increase of \code{x}), \code{"d"} for decreasing trend 
#' (\code{y} decreases with decrease of \code{x}). Default value is \code{NA}.
#' If the default value is passed, then trend will be identified based on the sign of the Spearman correlation 
#' coefficient between \code{x} and \code{y} on complete cases.
#'
#'@return The command \code{mdt.bin} generates a list of two objects. The first object, data frame \code{summary.tbl}
#' presents a summary table of final binning, while \code{x.trans} is a vector of discretized values. 
#' In case of single unique value for \code{x} or \code{y} in complete cases (cases different than special cases), 
#' it will return data frame with info.
#'
#'@examples
#'suppressMessages(library(monobin))
#'data(gcd)
#'amt.bin <- mdt.bin(x = gcd$amount, y = gcd$qual)
#'amt.bin[[1]]
#'table(amt.bin[[2]])
#'#force decreasing trend
#'mdt.bin(x = gcd$amount, y = gcd$qual, force.trend = "d")[[1]]
#'
#'@importFrom Hmisc cut2
#'@import dplyr
#'@export
mdt.bin <- function(x, y, g = 50, sc = c(NA, NaN, Inf, -Inf), sc.method = "together", y.type = NA, 
			 min.pct.obs = 0.05, min.avg.rate = 0.01, force.trend = NA) {
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
	#special cases
	if	(nrow(d.sc) > 0) {	
		if	(sc.method[1] == "together") {
			d.sc$bin <- "SC"
			} else {
			d.sc$bin <- as.character(d.sc$x)
			}
		ds.sc <- iso.summary(tbl = d.sc, bin = "bin")
		ds.sc$type <- "special cases"
		} else {
		ds.sc <- data.frame()
		}
	#complete cases	
	ds.cc <- mdt(tbl = d.cc, g = g, min.obs = min.obs, min.rate = min.rate, 
			 y.check = y.check, force.trend = force.trend)
	ds <- as.data.frame(bind_rows(ds.sc, ds.cc))
	ds <- woe.calc(tbl = ds, y.check = y.check)
	sc.u <- unique(sc)
	sc.g <- ds$bin[ds$type%in%"special cases"]
	x.mins <- ds$x.min[!ds$bin%in%sc.u & !ds$bin%in%"SC"]
	x.maxs <- ds$x.max[!ds$bin%in%sc.u & !ds$bin%in%"SC"]
	x.trans <- slice.variable(x.orig = d$x, x.lb = x.mins, x.ub = x.maxs, 
					  sc.u = sc.u, sc.g = sc.g) 
return(list(summary.tbl = ds, x.trans = x.trans))
}

#monotone decision tree
mdt <- function(tbl, g, min.obs, min.rate, y.check, force.trend) {
	if	(is.na(force.trend)) {
		cor.obs <- sign(cor(tbl$y, tbl$x, method = "spearman", use = "complete.obs"))
		} else {
		cor.obs <- ifelse(force.trend == "i", 1, -1)
		}
	cor.dir <- ifelse(cor.obs == 1, "<", ">")
	#initialize the process
	do.splits <- TRUE
	tree.info <- data.frame(NODE = 1, 
					NOBS = nrow(tbl), 
					FILTER = NA,
					TERMINAL = "SPLIT",
					BREAKS = NA,
					LR.NODE = NA,
					THRESHOLD = ifelse(cor.obs == 1, -Inf, Inf), 
					stringsAsFactors = FALSE)
	#start splitting
	while	(do.splits) {
		to.calculate <- which(tree.info$TERMINAL == "SPLIT")
		for	(i in to.calculate) {
      		if (!is.na(tree.info[i, "FILTER"])) {
        		tbl.l <- subset(tbl, eval(parse(text = tree.info[i, "FILTER"])))
      		} else {
        		tbl.l <- tbl
      		}
			tbl.l <- cbind.data.frame(tbl.l, bin = cut2(tbl.l$x, g = g))
      		splitting <- node.split(tbl = tbl.l, 
							min.obs = min.obs, 
							min.rate = min.rate, 
							y.check = y.check,
							node.thr = tree.info[i, "THRESHOLD"],
							node.lr = tree.info[i, "LR.NODE"],
							monodir = cor.dir)   
			if	(is.na(splitting[1]))  {
				leaf.adj <- mean(tbl.l$y)
				tree.info[tree.info$TERMINA%in%"SPLIT", "THRESHOLD"] <- leaf.adj
				tree.info[i, "TERMINAL"] <- "LEAF"
				next
				} 
      		tmp.splitter <- unname(splitting[2])
			mn <- max(tree.info$NODE) 
			tmp.filter <- c(paste("x", "<", tmp.splitter),
            	          	    paste("x", ">=", tmp.splitter))
			if	(!is.na(tree.info[i, "FILTER"])) {
				tmp.filter <- paste(tree.info[i, "FILTER"], tmp.filter, sep = " & ")
				} 
			tmp.nobs <- sapply(tmp.filter, FUN = function(i, x) {		
									 nrow(subset(x = x, subset = eval(parse(text = i))))
									 }, x = tbl.l) 
			children <- data.frame(NODE = c(mn + 1, mn + 2),
						     NOBS = tmp.nobs,
						     FILTER = tmp.filter,
						     TERMINAL = rep("SPLIT", 2),
						     BREAKS = tmp.splitter, 
						     LR.NODE = c("L", "R"),
						     THRESHOLD = splitting["y.r"],
						     row.names = NULL)    
			tree.info[i, "TERMINAL"] <- "PARENT"
			tree.info$THRESHOLD[tree.info$TERMINAL%in%"SPLIT"] <- splitting["y.r"]
			tree.info <- rbind(tree.info, children)
    			}
		do.splits <- !all(tree.info$TERMINAL != "SPLIT")
  		}
	cp <- sort(unique(tree.info$BREAKS[tree.info$TERMINAL%in%"LEAF"]))
	tbl$bin <- cut2(tbl$x, cuts = c(-Inf, cp, Inf))
	res <- iso.summary(tbl = tbl[, c("y", "x", "bin")], bin = "bin")
	res <- res[order(res$x.avg), ]
	res$bin <- format.bin(x.lb = res$x.min, x.ub = res$x.max)
	res$type <- "complete cases"
return(as.data.frame(res))
}
#node split function
node.split <- function(tbl, min.obs, min.rate, y.check, node.thr, node.lr, monodir) {
	tbl.s <- tbl %>%
		   group_by(bin) %>%
		   summarise(n = n(),
		   y.sum = sum(y),
		   y.avg = mean(y), 
		   x.min = min(x),
		   x.max = max(x)) %>%
		   ungroup() %>%
		   mutate(n.cs = cumsum(n),
			    y.cs = cumsum(y.sum),
			    y.cs.a = y.cs / n.cs,
			    n.cs.rev = cumsum(rev(n)),
			    y.cs.rev = cumsum(rev(y.sum)),
			    y.cs.a.rev = y.cs.rev / n.cs.rev)
	if	(y.check == "bina") { 
		lt <- min(which(tbl.s$n.cs >= min.obs & tbl.s$y.cs >= min.rate)) + 1
		ut <- nrow(tbl.s) - min(which(tbl.s$n.cs.rev >= min.obs & tbl.s$y.cs.rev >= min.rate)) - 1
		} else {
		lt <- min(which(tbl.s$n.cs >= min.obs & tbl.s$y.cs.a >= min.rate)) + 1
		ut <- nrow(tbl.s) - min(which(tbl.s$n.cs.rev >= min.obs & tbl.s$y.cs.a.rev >= min.rate)) - 1
		}
	if	(lt > ut) {
		return(NA)
		}
	nr <- nrow(tbl.s)
	eval.exp <- paste("lt.avg", monodir, "ut.avg")
	if	(is.na(node.lr)) {
		rt.dir <- ifelse(monodir == ">", "<", ">")
		} else {
		if	(node.lr == "L") {
			rt.dir <- monodir
			} else {
			rt.dir <- ifelse(monodir == ">", "<", ">")
			}
		}
	eval.exp <- paste0(eval.exp, " & (lt.avg ", rt.dir, " ", node.thr, " & ut.avg ", rt.dir, " ", node.thr, ")") 
	sp <- c()
	for	(i in lt:ut) {
		lt.avg <- tbl.s[i - 1, "y.cs.a"]
		ut.avg <- tbl.s[nr - i + 1, "y.cs.a.rev"]
		cond <- eval(parse(text = eval.exp))
		if	(cond) {sp <- c(sp, i)}
		}
	if	(length(sp) == 0) {
		return(NA)
		}
	if	(y.check == "bina") {
		sp.metric <- "gini(y = tbl$y, x = tbl$x, sp = split.l)"
		sp.indx <- "which.max(ssv)"
		} else {
		sp.metric <- "sse(y = tbl$y, x = tbl$x, sp = split.l)"
		sp.indx <- "which.min(ssv)"
		}
	splits <- tbl.s$x.min[sp]
	sl <- length(splits)
	ssv <- rep(NA, sl)
	for	(i in 1:sl) {
		split.l <- splits[i]
		ssv[i] <- eval(parse(text = sp.metric))
		}
	split.at <- splits[eval(parse(text = sp.indx))]
	y.l <- mean(tbl$y[tbl$x < split.at])
	y.r <- mean(tbl$y[tbl$x >= split.at])
return(c(ssv = min(ssv), split = split.at, y.l = y.l, y.r = y.r))
}
#sse
sse <- function(y, x, sp) {
	res <- sum((y[x < sp] - mean(y[x < sp]))^2) +
		 sum((y[x >= sp] - mean(y[x >= sp]))^2) 
return(res)	
}
#gini
gini <- function(y, x, sp) {
	ct <- table(y, x < sp)
	nx <-  apply(ct, 2, sum)
	n <- sum(ct)
	pxy <- ct / matrix(rep(nx, each = 2), nrow = 2)
	omega <- matrix(rep(nx, each = 2), nrow = 2) / n
	res <- -sum(omega * pxy * (1 - pxy))
return(res)
}



