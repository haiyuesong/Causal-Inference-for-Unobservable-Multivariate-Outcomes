// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
using namespace Rcpp;
using arma::mat; using arma::vec; using arma::uword;

inline double fisher_z(double r) {
  if (r >  0.999999) r =  0.999999;
  if (r < -0.999999) r = -0.999999;
  return 0.5 * std::log((1.0 + r) / (1.0 - r));
}

// Build conditioning sets from split-1 correlations (top-k by p-value)
 // [[Rcpp::export]]
List condset_topk_cpp(const arma::mat& X,
                      const double r_thresh = 0.1,
                      const double alpha = 0.05,
                      const int cap_k = 12) {
  const uword T = X.n_rows, p = X.n_cols;
  arma::rowvec mu = arma::mean(X, 0);
  arma::mat Xc = X.each_row() - mu;
  arma::rowvec ss = arma::sqrt(arma::sum(arma::square(Xc), 0));
  for (uword j=0;j<p;++j) if (ss(j) > 0) Xc.col(j) /= ss(j);
  arma::mat R = (Xc.t() * Xc) / (T - 1);

  const double zthr = fisher_z(r_thresh);
  const double se   = 1.0 / std::sqrt((double)T - 3.0);

  List out(p);
  for (uword k=0;k<p;++k) {
    std::vector<std::pair<double,int>> cand; cand.reserve(p);
    for (uword j=0;j<p;++j) if (j!=k) {
      double z  = fisher_z(R(k,j)) / se;
      double pv = 1.0 - R::pnorm(std::fabs(z) - zthr, 0.0, 1.0, 1, 0);
      if (pv <= alpha && std::isfinite(pv)) cand.emplace_back(pv, (int)j+1);
    }
    std::sort(cand.begin(), cand.end(),
      [](const std::pair<double,int>& a, const std::pair<double,int>& b){
        return a.first < b.first;
      });
    int take = std::min((int)cand.size(), cap_k);
    IntegerVector idx(take);
    for (int t=0;t<take;++t) idx[t] = cand[t].second;
    out[k] = idx;
  }
  return out;
}

// Batched F-statistics for all ordered pairs (i->j and j->i)
 // [[Rcpp::export]]
arma::mat batched_F_cpp(const arma::mat& Y, const List& cond_sets) {
  const uword T = Y.n_rows, p = Y.n_cols;
  const uword TL = T - 1;
  arma::mat RY = Y.rows(1, T-1);
  arma::mat LY = Y.rows(0, T-2);
  arma::vec ones = arma::ones<arma::vec>(TL);
  arma::mat F(p, p, arma::fill::zeros);

  for (uword i=0; i<p; ++i) {
    IntegerVector ci = cond_sets[i];
    std::vector<uword> ni; ni.reserve(ci.size());
    for (int t=0;t<ci.size();++t) { int v=ci[t]; if (v>=1 && (uword)v-1!=i) ni.push_back((uword)v-1); }
    std::sort(ni.begin(), ni.end());

    for (uword j=i+1; j<p; ++j) {
      IntegerVector cj = cond_sets[j];
      std::vector<uword> nj; nj.reserve(cj.size());
      for (int t=0;t<cj.size();++t) { int v=cj[t]; if (v>=1 && (uword)v-1!=j) nj.push_back((uword)v-1); }
      std::sort(nj.begin(), nj.end());

      std::vector<uword> cond_idx;
      std::set_intersection(ni.begin(), ni.end(), nj.begin(), nj.end(), std::back_inserter(cond_idx));

      // j as target: test i->j
      arma::vec yj = RY.col(j);
      uword k_red_j = 2 + cond_idx.size();
      arma::mat Xrj(TL, k_red_j);
      Xrj.col(0)=ones; Xrj.col(1)=LY.col(j);
      for (uword c=0;c<cond_idx.size();++c) Xrj.col(2+c)=LY.col(cond_idx[c]);
      arma::mat Xfj(TL, k_red_j+1); Xfj.cols(0,k_red_j-1)=Xrj; Xfj.col(k_red_j)=LY.col(i);

      arma::vec bj  = arma::solve(Xrj, yj);
      arma::vec rj  = yj - Xrj*bj; double ssr_rj = arma::dot(rj,rj);
      arma::vec bjf = arma::solve(Xfj, yj);
      arma::vec rjf = yj - Xfj*bjf; double ssr_fj = arma::dot(rjf,rjf);
      int dfj = (int)TL - (int)Xfj.n_cols;
      double Fj = 0.0; if (dfj>0 && ssr_fj>0.0) Fj = ((ssr_rj-ssr_fj)/1.0)/(ssr_fj/(double)dfj);

      // i as target: test j->i
      arma::vec yi = RY.col(i);
      uword k_red_i = 2 + cond_idx.size();
      arma::mat Xri(TL, k_red_i);
      Xri.col(0)=ones; Xri.col(1)=LY.col(i);
      for (uword c=0;c<cond_idx.size();++c) Xri.col(2+c)=LY.col(cond_idx[c]);
      arma::mat Xfi(TL, k_red_i+1); Xfi.cols(0,k_red_i-1)=Xri; Xfi.col(k_red_i)=LY.col(j);

      arma::vec bi  = arma::solve(Xri, yi);
      arma::vec ri  = yi - Xri*bi; double ssr_ri = arma::dot(ri,ri);
      arma::vec bif = arma::solve(Xfi, yi);
      arma::vec rif = yi - Xfi*bif; double ssr_fi = arma::dot(rif,rif);
      int dfi = (int)TL - (int)Xfi.n_cols;
      double Fi = 0.0; if (dfi>0 && ssr_fi>0.0) Fi = ((ssr_ri-ssr_fi)/1.0)/(ssr_fi/(double)dfi);

      F(i,j)=Fj; F(j,i)=Fi;
    }
  }
  return F;
}
