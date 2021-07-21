checks.init <- function(x, y, sc, sc.method, y.type, force.trend) {

	scm.opts <- c("together", "separately")
	yt.opts <- c(NA, "bina", "cont")
	ft.opts <- c(NA, "i", "d")

	cond.01 <- !is.numeric(x) | !is.numeric(y) | ifelse(length(sc) == 1, ifelse(is.numeric(sc) | is.na(sc), FALSE,  TRUE), 
							    ifelse(is.numeric(sc), FALSE, TRUE))
	cond.02 <- !sc.method[1]%in%scm.opts
	cond.03 <- !y.type[1]%in%yt.opts
	cond.04 <- !force.trend[1]%in%ft.opts
	cond.all <- c(cond.01, cond.02, cond.03, cond.04)	
	if	(sum(cond.all) > 0) {
		which.cond <- min(which(cond.all))
		} else {
		which.cond <- 0
		}
	error <- switch(as.character(which.cond), 
			    "0" = "NULL", 
			    "1" = "stop('x, y & sc have to be numeric vectors')",
			    "2" = "stop(paste0('sc.method  argument has to be one of: ', 
				        paste0(scm.opts, collapse = ', ')))",	
			    "3" = "stop(paste0('y.type  argument has to be one of: ', 
				        paste0(yt.opts, collapse = ', ')))",
			    "4" = "stop(paste0('force.trend  argument has to be one of: ', 
				        paste0(ft.opts, collapse = ', ')))")
return(eval(parse(text = error)))
}

checks.iter <- function(d, d.cc, y.type) {
	cond.01 <- length(unique(d.cc$y)) == 1
	cond.02 <- length(unique(d.cc$x)) == 1
	if	(is.na(y.type)) {
		y.check <- ifelse(sum(d$y[!is.na(d$y)]%in%c(0, 1)) == length(d$y[!is.na(d$y)]), 
					"bina", "cont")
		cond.03 <- FALSE
		} else {
		if	(y.type[1] == "bina") {
			cond.03 <- !sum(d$y[!is.na(d$y)]%in%c(0, 1)) == length(d$y[!is.na(d$y)])
			} else {
			cond.03 <- FALSE
			}
		y.check <- y.type[1]
		}
	cond4 <- nrow(d.cc) == 0
	cond.all <- c(cond.01, cond.02, cond.03, cond4)
		if	(sum(cond.all) > 0) {
		which.cond <- min(which(cond.all))
		} else {
		which.cond <- 0
		}
	msger <- switch(as.character(which.cond), 
			    "1" = "data.frame(bin = 'y has single unique value for complete cases')",
			    "2" = "data.frame(bin = 'x has single unique value for complete cases')",	
			    "3" = "stop('y is not 0/1 varibale')",
			    "4" = "data.frame(bin = 'no complete cases')")
return(list(which.cond, msger, y.check))
}
