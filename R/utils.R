# To cut data into batches
chunk <- function(x,n){
  if (n==1){
    list(x)
  }
  else{
    split(x, cut(seq_along(x), n, labels = FALSE))
  }
}



procrustes <- function(A, B){
  # Rotation matrix T
  svd.results <- svd(B %*% t(A))
  U <- svd.results$u
  V <- svd.results$v
  TT <- V %*% t(U)

  # B transformed
  B.transformed <- TT %*% B

  # Error after superimposition
  RSS <- norm(A - B.transformed,  type = "F")

  # Return
  return(list(A.normalized = A, B.normalized = B, rotation.mtx = TT, B.transformed = B.transformed, RSS = RSS))
}
