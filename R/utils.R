# To cut data into batches
chunk <- function(x,n){
  if (n==1){
    list(x)
  }
  else{
    split(x, cut(seq_along(x), n, labels = FALSE))
  }
}

# Regularisation soft-threshold function
S.z.g <- function(S.z,S.g,b.a){
  to_return <- array(0,dim=dim(S.z))
  to_return <-
    (S.z - S.g) * ((S.z > 0) * (S.g < abs(S.z))) +
    (S.z + S.g) * ((S.z < 0) * (S.g < abs(S.z)))
  return(to_return)
}
