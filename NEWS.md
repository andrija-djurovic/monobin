# monobin 0.1.0
This is a small release focusing on fixing the bin indexing for categorization of numeric risk factors. 
The change affects only the second element of functions' output and aligns it with summary table indexing of the bins. 

# monobin 0.1.1
Changes:
1. check for sc in the function check.init<br/>
   ```!is.numeric(sc)``` replaced with ```ifelse(length(sc) == 1, ifelse(is.numeric(sc) | is.na(sc), FALSE,  TRUE), ifelse(is.numeric(sc), FALSE, TRUE))```

2. check of y.type in the fucntion check.iter<br/>
   a) deleted the following line of the code form check.iter function: ```cond.03 <- FALSE```<br/>
   b) additionally quatation marks corrected:  ``"3" = "stop('y is not 0/1 varibale)'",`` corrected to ```"3" = "stop('y is not 0/1 varibale')",```<br/>
   
3. github page now available: www.github.com/andrija-djurovic/monobin
   
 Until new CRAN realise, github version 0.1.1 will be updated.
