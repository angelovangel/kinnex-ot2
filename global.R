# 

make_s <- function(samples, plex) {
# samples is a vector with len 16, i.e. "s1" "s2" "."  "."  "." "."
# should be len 6 (each is two columns, so 6 double columns x 16 wells = 96)
  c(
    sapply(
      samples,
      function(x){
        c(rep(x, plex), rep('.', 16-plex))  
      }
    )
  )
}
