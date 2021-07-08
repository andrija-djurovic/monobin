# monobin 0.1.0
This is a small release focusing on fixing the bin indexing for categorization of numeric risk factors. 
The change affects only the second element of functions' output and aligns it with summary table indexing of the bins. 

# monobin 0.1.1
Planned actions for the next CRAN realise:<br/>
1. change check for sc (function check.init)<br/>
  replace: ```!is.numeric(sc)``` with ```ifelse(length(sc) == 1, ifelse(is.numeric(sc) | is.na(sc), FALSE,  TRUE), ifelse(is.numeric(sc), FALSE, TRUE))```

2. check of y.type<br/>
   delete this line of the code form check.iter function ```cond.03 <- FALSE```<br/>
   additionally quatation marked corrected: instead of ``"3" = "stop('y is not 0/1 varibale)'",`` use ```"3" = "stop('y is not 0/1 varibale')",```<br/>
   
 Until new CRAN realise, github version 0.1.1 will be updated.
