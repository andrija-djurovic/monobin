# monobin 0.1.0
This is a small release focusing on fixing the bin indexing for categorization of numeric risk factors. 
The change affects only the second element of functions' output and aligns it with summary table indexing of the bins. 

# monobin 0.1.1
Changes:
1. special case value functions' argument now accepts NA only, not producing the error that numeric vector is needed (check for sc argument in the function check.init):<br/>
   ```!is.numeric(sc)``` replaced with ```ifelse(length(sc) == 1, ifelse(is.numeric(sc) | is.na(sc), FALSE,  TRUE), ifelse(is.numeric(sc), FALSE, TRUE))```

2. check if binary target contains only 0/1 values (y.type argument check in the function check.iter):<br/>
   a) ```cond.03 <- FALSE``` moved to the proper place in if clause when checking for binary target type <br/>
   b) quotation marks corrected:  ``"3" = "stop('y is not 0/1 varibale)'",`` corrected to ```"3" = "stop('y is not 0/1 varibale')",```<br/>
   
3. github page now available: www.github.com/andrija-djurovic/monobin
   
