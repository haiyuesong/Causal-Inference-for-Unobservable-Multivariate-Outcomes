# define a function to find high-correlated regions
high_corr <- function(roi1, roi2, std_fisherz, threshold=0.8,std_err){  
  std_fisherz_threshold <- fisherz(threshold)/std_err
  z_diff <- std_fisherz[roi1,roi2]-std_fisherz_threshold
  p_value <- 1 - pnorm(z_diff) # one-tailed test with null hypothesis: correlation < threshold
  return(p_value)
}

# intersection
intersection_conditional_granger <- function(roi_1,roi_2,ts,order_p,high_corr){
  if(roi_1==roi_2) return(NA)
  else{
    cond_roi <- intersect(which(high_corr[,roi_1]),which(high_corr[,roi_2]))
    if(length(cond_roi)>0){
      return(condGranger(cbind(ts[,roi_1],ts[,roi_2],ts[,cond_roi]),nx=1,ny=1,order=order_p,perm = T)$orig)
    }else{
      return(grangertest(y~x,order=order_p,data=cbind(y=ts[,roi_2],x=ts[,roi_1]))$`F`[2])
    }
  }
}
aal2_txt <- read.table("aal2_ROI_names.txt")
vec_to_square_matrix <- function(vec, p) {
  mat <- matrix(NA, nrow = p, ncol = p)
  for (i in 1:p) {
    start_index <- (i - 1) * (p - 1) + 1
    mid_index   <- start_index + (i - 1) - 1
    end_index   <- start_index + (p - 2)
    if (i > 1) mat[i, 1:(i - 1)] <- vec[start_index:mid_index]
    if (i < p) mat[i, (i + 1):p] <- vec[(mid_index + 1):end_index]
  }
  colnames(mat) <- rownames(mat) <- aal2_txt$V2
  return(mat)
}
mat_to_df <- function(mat, label = "value") {
  as.data.frame(as.table(mat)) %>%
    rename(Row = Var1, Col = Var2, !!label := Freq)
}
