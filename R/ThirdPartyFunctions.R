### v1.0.1 7FEB2024
#This function was obtained from Reddit and written by Josh Bredeweg, user jbraids1421. 

################
##countSpaces
#'Count the Spaces in a Character String
#'
#'Counts the number of spaces in a character string.
#'
#'@param x A character string.
#'
#'@return A numeric value representing the number of spaces (' ') in x.
#'
#'@author Josh Bredeweg, Reddit user jbraids1421
#'
#'@export
#'
#'@source https://www.reddit.com/r/rstats/comments/2th8ic/function_to_count_the_number_of_white_spaces/
#'
#'@examples
#'countSpaces("abc def")
#'
countSpaces <- function(x){
  counter <- 0
  coll <- numeric()
  vec <- strsplit(x," ")[[1]]
  for(i in 1:length(vec)){
    if (vec[i]==""){
      counter <- counter+1
    }
    else{
      if (counter!=0) coll <- c(coll,counter)
      counter <- 1
    }
  }
  coll
}
