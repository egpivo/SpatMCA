// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppParallel.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;
using namespace arma;

arma::mat cubicmatrix(const vec z1)
{
  
  vec h(z1.n_elem - 1);
  mat Q, R;
  int p = z1.n_elem;
  Q.zeros(p - 2, p);
  R.zeros(p - 2, p - 2);
  
  for (unsigned i = 0; i < z1.n_elem - 1; i++)
    h[i] = z1[i + 1] - z1[i];
  
  for (unsigned j = 0; j < p - 2; ++j)
  {
    Q(j, j) = 1 / h[j];
    Q(j, j + 1) = -1 / h[j] - 1 / h[j + 1];
    Q(j, j + 2) = 1 / h[j + 1];
    R(j, j) = (h[j] + h[j + 1]) / 3;
  }
  for (unsigned j = 0; j < p - 3; ++j)
    R(j, j + 1) = R(j + 1, j) = h[j + 1] / 6;
  
  return (Q.t() * solve(R, Q));
}

struct tpm : public RcppParallel::Worker
{
  const mat &P;
  mat &L;
  int p;
  int d;
  tpm(const mat &P, mat &L, int p, int d) : P(P), L(L), p(p), d(d) {}
  void operator()(std::size_t begin, std::size_t end)
  {
    for (std::size_t i = begin; i < end; i++)
    {
      for (unsigned j = 0; j < p; j++)
      {
        if (j > i)
        {
          if (d == 2)
          {
            double r = sqrt(pow(P(i, 0) - P(j, 0), 2) + (pow(P(i, 1) - P(j, 1), 2)));
            L(i, j) = r * r * log(r) / (8.0 * datum::pi);
          }
          else if (d == 1)
          {
            double r = sqrt(pow(P(i, 0) - P(j, 0), 2));
            L(i, j) = pow(r, 3) / 12;
          }
          else if (d == 3)
          {
            double r = sqrt(pow(P(i, 0) - P(j, 0), 2) +
                            pow(P(i, 1) - P(j, 1), 2) +
                            pow(P(i, 2) - P(j, 2), 2));
            L(i, j) = -r / (8.0 * datum::pi);
          }
        }
      }
      
      L(i, p) = 1;
      for (unsigned k = 0; k < d; ++k)
      {
        L(i, p + k + 1) = P(i, k);
      }
    }
  }
};

arma::mat tpmatrix(const arma::mat P)
{
  mat L, Lp, Ip;
  int p = P.n_rows, d = P.n_cols;
  
  L.zeros(p + d + 1, p + d + 1);
  Ip.eye(p + d + 1, p + d + 1);
  tpm tpm(P, L, p, d);
  parallelFor(0, p, tpm);
  L = symmatl(L);
  Lp = inv_sympd(L + 1e-8 * Ip);
  Lp.shed_cols(p, p + d);
  Lp.shed_rows(p, p + d);
  L.shed_cols(p, p + d);
  L.shed_rows(p, p + d);
  mat result = Lp.t() * (L * Lp);
  return (result);
}

//' Internal function: thin-plane spline matrix
//' @keywords internal
//' @param z A new location matrix
//' @param P A location matrix
//' @param Phi An eigenvector matrix
//' @return A thin-plane spline matrix
//'
// [[Rcpp::export]]
arma::mat tpm2(const arma::mat z, const arma::mat P, const arma::mat Phi)
{
  mat L;
  int p = P.n_rows, d = P.n_cols, K = Phi.n_cols;
  L.zeros(p + d + 1, p + d + 1);
  tpm tpm(P, L, p, d);
  parallelFor(0, p, tpm);
  L = L + L.t();
  mat Phi_star, para(p + d + 1, K);
  Phi_star.zeros(p + d + 1, K);
  Phi_star.rows(0, p - 1) = Phi;
  para = solve(L, Phi_star);
  int pnew = z.n_rows;
  mat eigen_fn(pnew, K);
  double psum, r;
  for (unsigned newi = 0; newi < pnew; newi++)
  {
    for (unsigned i = 0; i < K; i++)
    {
      psum = 0;
      for (unsigned j = 0; j < p; j++)
      {
        if (d == 2)
        {
          r = sqrt(pow(z(newi, 0) - P(j, 0), 2) + (pow(z(newi, 1) - P(j, 1), 2)));
          if (r != 0)
            psum += para(j, i) * r * r * log(r) / (8.0 * datum::pi);
        }
        else if (d == 1)
        {
          r = norm(z.row(newi) - P.row(j), 'f');
          if (r != 0)
            psum += para(j, i) * pow(r, 3) / 12;
        }
        else if (d == 3)
        {
          double r = sqrt(pow(z(newi, 0) - P(j, 0), 2) +
                          pow(z(newi, 1) - P(j, 1), 2) +
                          pow(z(newi, 2) - P(j, 2), 2));
          if (r != 0)
            psum += -para(j, i) * r / (8.0 * arma::datum::pi);
        }
      }
      if (d == 1)
        eigen_fn(newi, i) = psum + para(p + 1, i) * z(newi, 0) + para(p, i);
      else if (d == 2)
        eigen_fn(newi, i) = psum + para(p + 1, i) * z(newi, 0) + para(p + 2, i) * z(newi, 1) + para(p, i);
      else if (d == 3)
        eigen_fn(newi, i) = psum + para(p + 1, i) * z(newi, 0) + para(p + 2, i) * z(newi, 1) + para(p + 3, i) * z(newi, 2) + para(p, i);
    }
  }
  return (eigen_fn);
}

void spatmca_tau1(mat &G, mat &C, mat &Lambda2, const mat matrixinv, const int p1, const int p2, const double zeta, const int maxit, const double tol)
{
  
  int p = G.n_rows;
  mat temp, U1, U2, V1, V2, R = G, Cold = C, Lambda2old = Lambda2;
  vec er(2), S1, S2;
  
  for (unsigned iter = 0; iter < maxit; iter++)
  {
    G = matrixinv * (zeta * (R + Cold) - Lambda2old);
    R = G;
    temp = zeta * G + Lambda2old;
    svd_econ(U1, S1, V1, temp.rows(0, p1 - 1));
    
    C.rows(0, p1 - 1) = U1.cols(0, V1.n_cols - 1) * V1.t();
    svd_econ(U2, S2, V2, temp.rows(p1, p1 + p2 - 1));
    C.rows(p1, p1 + p2 - 1) = U2.cols(0, V1.n_cols - 1) * V2.t();
    Lambda2 = Lambda2old + zeta * (G - C);
    
    er[0] = norm(G - C, "fro") / sqrt(p / 1.0);
    er[1] = norm((C - Cold), "fro") / sqrt(p / 1.0);
    
    if (max(er) <= tol)
      break;
    Cold = C;
    Lambda2old = Lambda2;
  }
}

void spatmca_tau2(mat &G, mat &R, mat &C, mat &Lambda1, mat &Lambda2,
                  const mat matrixinv, double tau2u, double tau2v, const int p1,
                  const int p2, const double zeta, const int maxit, const double tol)
{
  
  int p = G.n_rows;
  int K = G.n_cols;
  
  mat temp, zero, one, U1, U2, V1, V2;
  vec er(4), S1, S2;
  mat Rold = R, Cold = C, Lambda1old = Lambda1, Lambda2old = Lambda2;
  zero.zeros(p, K);
  one.ones(p, K);
  for (unsigned iter = 0; iter < maxit; iter++)
  {
    G = matrixinv * (zeta * (Rold + Cold) - Lambda1old - Lambda2old);
    R.rows(0, p1 - 1) = sign((Lambda1old.rows(0, p1 - 1) / zeta + G.rows(0, p1 - 1))) % max(zero.rows(0, p1 - 1), abs((Lambda1old.rows(0, p1 - 1) / zeta + G.rows(0, p1 - 1))) - tau2u * one.rows(0, p1 - 1) / zeta);
    R.rows(p1, p1 + p2 - 1) = sign((Lambda1old.rows(p1, p1 + p2 - 1) / zeta + G.rows(p1, p1 + p2 - 1))) % max(zero.rows(p1, p1 + p2 - 1), abs((Lambda1old.rows(p1, p1 + p2 - 1) / zeta + G.rows(p1, p1 + p2 - 1))) - tau2v * one.rows(p1, p1 + p2 - 1) / zeta);
    temp = zeta * G + Lambda2old;
    svd_econ(U1, S1, V1, temp.rows(0, p1 - 1));
    C.rows(0, p1 - 1) = U1.cols(0, V1.n_cols - 1) * V1.t();
    svd_econ(U2, S2, V2, temp.rows(p1, p1 + p2 - 1));
    C.rows(p1, p1 + p2 - 1) = U2.cols(0, V1.n_cols - 1) * V2.t();
    
    Lambda1 = Lambda1old + zeta * (G - R);
    Lambda2 = Lambda2old + zeta * (G - C);
    
    er[0] = norm(G - R, "fro") / sqrt(p / 1.0);
    er[1] = norm((R - Rold), "fro") / sqrt(p / 1.0);
    er[2] = norm(G - C, "fro") / sqrt(p / 1.0);
    er[3] = norm((C - Cold), "fro") / sqrt(p / 1.0);
    
    if (max(er) <= tol)
      break;
    Rold = R;
    Cold = C;
    Lambda1old = Lambda1;
    Lambda2old = Lambda2;
  }
}

struct spatmcacv_p : public RcppParallel::Worker
{
  const mat &X;
  const mat &Y;
  const int K;
  const mat &Omega1;
  const mat &Omega2;
  const vec &tau1u;
  const vec &tau1v;
  const vec &nk;
  const int p1;
  const int p2;
  const int maxit;
  const double tol;
  cube &output;
  cube &S12train;
  cube &S12traint;
  cube &S12valid;
  cube &S1221;
  cube &S2112;
  vec &zetatemp;
  cube &Ucv;
  cube &Vcv;
  cube &Lmbd12cv;
  cube &Lmbd22cv;
  spatmcacv_p(const mat &X, const mat &Y, const int K, const mat &Omega1, const mat &Omega2, const vec &tau1u,
              const vec &tau1v, const vec &nk, const int p1, const int p2, const int maxit,
              const double tol, cube &output, cube &S12train, cube &S12traint, cube &S12valid, cube &S1221, cube &S2112,
              vec &zetatemp, cube &Ucv, cube &Vcv, cube &Lmbd12cv, cube &Lmbd22cv) : X(X), Y(Y), K(K), Omega1(Omega1), Omega2(Omega2), tau1u(tau1u), tau1v(tau1v),
              nk(nk), p1(p1), p2(p2), maxit(maxit), tol(tol), output(output), S12train(S12train), S12traint(S12traint), S12valid(S12valid), S1221(S1221), S2112(S2112), zetatemp(zetatemp),
              Ucv(Ucv), Vcv(Vcv), Lmbd12cv(Lmbd12cv), Lmbd22cv(Lmbd22cv) {}
  void operator()(std::size_t begin, std::size_t end)
  {
    mat Ip;
    Ip.eye(Y.n_cols, Y.n_cols);
    for (std::size_t k = begin; k < end; k++)
    {
      mat Uoldtemp, Voldtemp, Goldtemp, G(p1 + p2, K), C(p1 + p2, K), Gamma2(p1 + p2, K), Gold(p1 + p2, K), Cold(p1 + p2, K), Gamma2old(p1 + p2, K), matrixinv(p1 + p2, p1 + p2), Ip;
      vec svdtemp, D;
      mat Ytrain = Y.rows(find(nk != (k + 1)));
      mat Yvalid = Y.rows(find(nk == (k + 1)));
      mat Xtrain = X.rows(find(nk != (k + 1)));
      mat Xvalid = X.rows(find(nk == (k + 1)));
      mat Theta(p1 + p2, p1 + p2);
      S12train.slice(k) = (Xtrain.t()) * Ytrain / Xtrain.n_rows;
      S12valid.slice(k) = Xvalid.t() * Yvalid / Xvalid.n_rows;
      svd_econ(Uoldtemp, svdtemp, Voldtemp, S12train.slice(k));
      zetatemp[k] = 10 * max(svdtemp);
      Ucv.slice(k) = Uoldtemp.cols(0, K - 1);
      Gold.rows(0, p1 - 1) = Uoldtemp.cols(0, K - 1);
      Cold.rows(0, p1 - 1) = Uoldtemp.cols(0, K - 1);
      Vcv.slice(k) = Voldtemp.cols(0, K - 1);
      Gold.rows(p1, p1 + p2 - 1) = Voldtemp.cols(0, K - 1);
      Cold.rows(p1, p1 + p2 - 1) = Voldtemp.cols(0, K - 1);
      Gamma2old.rows(0, p1 - 1) = Lmbd12cv.slice(k) = S12train.slice(k) * Gold.rows(p1, p1 + p2 - 1);
      S12traint.slice(k) = S12train.slice(k).t();
      Gamma2old.rows(p1, p1 + p2 - 1) = Lmbd22cv.slice(k) = S12traint.slice(k) * Gold.rows(0, p1 - 1);
      Ip.eye(p1 + p2, p1 + p2);
      Theta.submat(0, p1, p1 - 1, p1 + p2 - 1) = S12train.slice(k) / 2;
      Theta.submat(p1, 0, p1 + p2 - 1, p1 - 1) = S12traint.slice(k) / 2;
      S1221.slice(k) = S12train.slice(k) * S12traint.slice(k);
      S2112.slice(k) = S12traint.slice(k) * S12train.slice(k);
      vec zero;
      zero.zeros(K);
      for (uword i = 0; i < tau1u.n_elem; i++)
      {
        G = Gold;
        C = Cold;
        Gamma2 = Gamma2old;
        Theta.submat(0, 0, p1 - 1, p1 - 1) = -tau1u[i] * Omega1;
        for (uword j = 0; j < tau1v.n_elem; j++)
        {
          Theta.submat(p1, p1, p1 + p2 - 1, p1 + p2 - 1) = -tau1v[j] * Omega2;
          if (j == 0 && i == 0)
          {
            output(i, j, k) = norm((S12valid.slice(k) - G.rows(0, p1 - 1) * diagmat(max(zero, svdtemp.subvec(0, K - 1))) * G.rows(p1, p1 + p2 - 1).t()), "fro");
          }
          else
          {
            matrixinv = 0.5 * inv_sympd(zetatemp[k] * Ip - Theta);
            spatmca_tau1(G, C, Gamma2, matrixinv, p1, p2, zetatemp[k], maxit, tol);
            D = max(zero, diagvec(G.rows(0, p1 - 1).t() * S12train.slice(k) * G.rows(p1, p1 + p2 - 1)));
            output(i, j, k) = norm((S12valid.slice(k) - G.rows(0, p1 - 1) * diagmat(D) * G.rows(p1, p1 + p2 - 1).t()), "fro");
          }
          if (j == 0)
          {
            Gold = G;
            Cold = C;
            Gamma2old = Gamma2;
          }
        }
      }
    }
  }
};

struct spatmcacv_pp : public RcppParallel::Worker
{
  const mat &X;
  const mat &Y;
  const int K;
  const mat &Omega1;
  const mat &Omega2;
  const vec &tau1u;
  const vec &tau1v;
  const vec &nk;
  const int p1;
  const int p2;
  const int maxit;
  const double tol;
  cube &output;
  cube &S12train;
  cube &S12traint;
  cube &S12valid;
  cube &S1221;
  cube &S2112;
  vec &zetatemp;
  cube &Ucv;
  cube &Vcv;
  cube &Lmbd12cv;
  cube &Lmbd22cv;
  spatmcacv_pp(const mat &X, const mat &Y, const int K, const mat &Omega1, const mat &Omega2, const vec &tau1u,
               const vec &tau1v, const vec &nk, const int p1, const int p2, const int maxit,
               const double tol, cube &output, cube &S12train, cube &S12traint, cube &S12valid, cube &S1221,
               cube &S2112, vec &zetatemp, cube &Ucv, cube &Vcv, cube &Lmbd12cv, cube &Lmbd22cv) : X(X), Y(Y), K(K), Omega1(Omega1), Omega2(Omega2), tau1u(tau1u), tau1v(tau1v),
               nk(nk), p1(p1), p2(p2), maxit(maxit), tol(tol), output(output), S12train(S12train), S12traint(S12traint), S12valid(S12valid), S1221(S1221), S2112(S2112), zetatemp(zetatemp),
               Ucv(Ucv), Vcv(Vcv), Lmbd12cv(Lmbd12cv), Lmbd22cv(Lmbd22cv) {}
  void operator()(std::size_t begin, std::size_t end)
  {
    mat Ip;
    Ip.eye(Y.n_cols, Y.n_cols);
    for (std::size_t k = begin; k < end; k++)
    {
      mat Uoldtemp, Voldtemp, Goldtemp, G(p1 + p2, K), C(p1 + p2, K), Gamma2(p1 + p2, K), Gold(p1 + p2, K), Cold(p1 + p2, K), Gamma2old(p1 + p2, K), matrixinv(p1 + p2, p1 + p2), Ip;
      vec svdtemp;
      mat Ytrain = Y.rows(find(nk != (k + 1)));
      mat Yvalid = Y.rows(find(nk == (k + 1)));
      mat Xtrain = X.rows(find(nk != (k + 1)));
      mat Xvalid = X.rows(find(nk == (k + 1)));
      mat Theta(p1 + p2, p1 + p2);
      S12train.slice(k) = (Xtrain.t()) * Ytrain / Xtrain.n_rows;
      S12valid.slice(k) = Xvalid.t() * Yvalid / Xvalid.n_rows;
      svd_econ(Uoldtemp, svdtemp, Voldtemp, S12train.slice(k));
      zetatemp[k] = 5 * max(svdtemp);
      Ucv.slice(k) = Uoldtemp.cols(0, K - 1);
      Gold.rows(0, p1 - 1) = Uoldtemp.cols(0, K - 1);
      Cold.rows(0, p1 - 1) = Uoldtemp.cols(0, K - 1);
      Vcv.slice(k) = Voldtemp.cols(0, K - 1);
      Gold.rows(p1, p1 + p2 - 1) = Voldtemp.cols(0, K - 1);
      Cold.rows(p1, p1 + p2 - 1) = Voldtemp.cols(0, K - 1);
      Gamma2old.rows(0, p1 - 1) = Lmbd12cv.slice(k) = S12train.slice(k) * Gold.rows(p1, p1 + p2 - 1);
      S12traint.slice(k) = S12train.slice(k).t();
      Gamma2old.rows(p1, p1 + p2 - 1) = Lmbd22cv.slice(k) = S12traint.slice(k) * Gold.rows(0, p1 - 1);
      Ip.eye(p1 + p2, p1 + p2);
      Theta.submat(0, p1, p1 - 1, p1 + p2 - 1) = S12train.slice(k) / 2;
      Theta.submat(p1, 0, p1 + p2 - 1, p1 - 1) = S12traint.slice(k) / 2;
      S1221.slice(k) = S12train.slice(k) * S12traint.slice(k);
      S2112.slice(k) = S12traint.slice(k) * S12train.slice(k);
    }
  }
};

struct spatmcacv_p2 : public RcppParallel::Worker
{
  const cube &S12train;
  const cube &S12traint;
  const cube &S12valid;
  const cube &S1221;
  const cube &S2112;
  const vec &zetatemp;
  const cube &Ucv;
  const cube &Vcv;
  const cube &Lmbd12cv;
  const cube &Lmbd22cv;
  const mat &Omega1;
  const mat &Omega2;
  int K;
  double tau1u;
  double tau1v;
  const vec &tau2u;
  const vec &tau2v;
  const vec &nk;
  int p1;
  int p2;
  int maxit;
  double tol;
  cube &output;
  spatmcacv_p2(const cube &S12train, const cube &S12traint, const cube &S12valid, const cube &S1221, const cube &S2112, const vec &zetatemp, const cube &Ucv, const cube &Vcv, const cube &Lmbd12cv,
               const cube &Lmbd22cv, const mat &Omega1, const mat &Omega2, int K, double tau1u, double tau1v, const vec &tau2u, const vec &tau2v,
               const vec &nk, int p1, int p2, int maxit, double tol, cube &output) : S12train(S12train), S12traint(S12traint), S12valid(S12valid), S1221(S1221), S2112(S2112), zetatemp(zetatemp), Ucv(Ucv), Vcv(Vcv), Lmbd12cv(Lmbd12cv),
               Lmbd22cv(Lmbd22cv), Omega1(Omega1), Omega2(Omega2), K(K), tau1u(tau1u), tau1v(tau1v), tau2u(tau2u), tau2v(tau2v), nk(nk), p1(p1), p2(p2), maxit(maxit), tol(tol), output(output) {}
  void operator()(std::size_t begin, std::size_t end)
  {
    for (std::size_t k = begin; k < end; k++)
    {
      mat G(p1 + p2, K), C(p1 + p2, K), R(p1 + p2, K), Gamma1(p1 + p2, K), Gamma2(p1 + p2, K), Gold(p1 + p2, K), Cold(p1 + p2, K), Rold(p1 + p2, K), Gamma1old(p1 + p2, K), Gamma2old(p1 + p2, K), Ip;
      mat Theta(p1 + p2, p1 + p2), matrixinv(p1 + p2, p1 + p2);
      vec D;
      Gold.rows(0, p1 - 1) = Ucv.slice(k);
      Cold.rows(0, p1 - 1) = Ucv.slice(k);
      Gold.rows(p1, p1 + p2 - 1) = Vcv.slice(k);
      Cold.rows(p1, p1 + p2 - 1) = Vcv.slice(k);
      Gamma2old.rows(0, p1 - 1) = S12train.slice(k) * Gold.rows(p1, p1 + p2 - 1);
      Gamma2old.rows(p1, p1 + p2 - 1) = S12traint.slice(k) * Gold.rows(0, p1 - 1);
      Ip.eye(p1 + p2, p1 + p2);
      Theta.submat(0, p1, p1 - 1, p1 + p2 - 1) = S12train.slice(k) / 2;
      Theta.submat(p1, 0, p1 + p2 - 1, p1 - 1) = S12traint.slice(k) / 2;
      Theta.submat(0, 0, p1 - 1, p1 - 1) = -tau1u * Omega1;
      Theta.submat(p1, p1, p1 + p2 - 1, p1 + p2 - 1) = -tau1v * Omega2;
      if (tau1u != 0)
      {
        mat M1 = inv_sympd(2 * zetatemp[k] * Ip.submat(0, 0, p1 - 1, p1 - 1) - 2 * Theta.submat(0, 0, p1 - 1, p1 - 1) - S1221.slice(k) / (2 * zetatemp[k]));
        matrixinv.submat(0, 0, p1 - 1, p1 - 1) = M1;
        matrixinv.submat(0, p1, p1 - 1, p1 + p2 - 1) = M1 * S12train.slice(k) / (2 * zetatemp[k]);
        matrixinv.submat(p1, 0, p1 + p2 - 1, p1 - 1) = matrixinv.submat(0, p1, p1 - 1, p1 + p2 - 1).t();
        matrixinv.submat(p1, p1, p1 + p2 - 1, p1 + p2 - 1) = (S12traint.slice(k) * matrixinv.submat(0, p1, p1 - 1, p1 + p2 - 1) + Ip.submat(p1, p1, p1 + p2 - 1, p1 + p2 - 1)) / (2 * zetatemp[k]);
        spatmca_tau1(Gold, Cold, Gamma2old, matrixinv, p1, p2, zetatemp[k], maxit, tol);
      }
      if (tau1v != 0)
      {
        if (tau1u == 0)
        {
          mat M2 = inv_sympd(2 * zetatemp[k] * Ip.submat(p1, p1, p1 + p2 - 1, p1 + p2 - 1) - 2 * Theta.submat(p1, p1, p1 + p2 - 1, p1 + p2 - 1) - S2112.slice(k) / (2 * zetatemp[k]));
          matrixinv.submat(p1, p1, p1 + p2 - 1, p1 + p2 - 1) = M2;
          matrixinv.submat(p1, 0, p1 + p2 - 1, p1 - 1) = M2 * S12traint.slice(k) / (2 * zetatemp[k]);
          matrixinv.submat(0, p1, p1 - 1, p1 + p2 - 1) = matrixinv.submat(p1, 0, p1 + p2 - 1, p1 - 1).t();
          matrixinv.submat(0, 0, p1 - 1, p1 - 1) = (S12train.slice(k) * matrixinv.submat(p1, 0, p1 + p2 - 1, p1 - 1) + Ip.submat(0, 0, p1 - 1, p1 - 1)) / (2 * zetatemp[k]);
          spatmca_tau1(Gold, Cold, Gamma2old, matrixinv, p1, p2, zetatemp[k], maxit, tol);
        }
        else
        {
          matrixinv = 0.5 * inv_sympd(zetatemp[k] * Ip - Theta);
          spatmca_tau1(Gold, Cold, Gamma2old, matrixinv, p1, p2, zetatemp[k], maxit, tol);
        }
      }
      if (tau1u == 0 && tau1v == 0)
      {
        matrixinv = 0.5 * inv_sympd(zetatemp[k] * Ip - Theta);
      }
      Rold = Gold;
      Gamma1old = 0 * Gamma2old;
      vec zero;
      zero.zeros(K);
      for (uword i = 0; i < tau2u.n_elem; i++)
      {
        G = Gold;
        R = Rold;
        C = Cold;
        Gamma1 = Gamma1old;
        Gamma2 = Gamma2old;
        for (uword j = 0; j < tau2v.n_elem; j++)
        {
          spatmca_tau2(G, R, C, Gamma1, Gamma2, matrixinv, tau2u[i], tau2v[j], p1, p2, zetatemp[k], maxit, tol);
          D = max(zero, diagvec(G.rows(0, p1 - 1).t() * S12train.slice(k) * G.rows(p1, p1 + p2 - 1)));
          output(i, j, k) = norm((S12valid.slice(k) - G.rows(0, p1 - 1) * diagmat(D) * G.rows(p1, p1 + p2 - 1).t()), "fro");
          if (j == 0)
          {
            Gold = G;
            Rold = R;
            Cold = C;
            Gamma1old = Gamma1;
            Gamma2old = Gamma2;
          }
        }
      }
    }
  }
};

//' Internal function: M-fold Cross-validation of SpatMCA
//' @keywords internal
//' @param sxr A location matrix for a variable X
//' @param sxr A location matrix for a variable Y
//' @param Yr A data matrix of Y
//' @param Yr A data matrix of X
//' @param M A number of folds
//' @param K The number of estimated eigenfunctions
//' @param tau1ur A range of tau1u
//' @param tau2ur A range of tau2u
//' @param tau1vr A range of tau1v
//' @param tau2vr A range of tau2v
//' @param nkr A vector of fold numbers
//' @param maxit A maximum number of iteration
//' @param tol A tolerance rate
//' @param l2ur A given tau2u
//' @param l2vr A given tau2v
//' @return A list of selected parameters
// [[Rcpp::export]]
List spatmcacv_rcpp(NumericMatrix sxr, NumericMatrix syr, NumericMatrix Xr, NumericMatrix Yr, int M, int K,
                    NumericVector tau1ur, NumericVector tau2ur, NumericVector tau1vr, NumericVector tau2vr,
                    NumericVector nkr, int maxit, double tol, NumericVector l2ur, NumericVector l2vr)
{
  int n = Yr.nrow(), p = Xr.ncol(), q = Yr.ncol(), d = sxr.ncol();
  mat X(Xr.begin(), n, p, false);
  mat Y(Yr.begin(), n, q, false);
  mat sx(sxr.begin(), p, d, false);
  mat sy(syr.begin(), q, d, false);
  colvec tau1u(tau1ur.begin(), tau1ur.size(), false);
  colvec tau2u(tau2ur.begin(), tau2ur.size(), false);
  colvec tau1v(tau1vr.begin(), tau1vr.size(), false);
  colvec tau2v(tau2vr.begin(), tau2vr.size(), false);
  colvec nk(nkr.begin(), nkr.size(), false);
  colvec l2u(l2ur.begin(), l2ur.size(), false);
  colvec l2v(l2vr.begin(), l2vr.size(), false);
  cube cv(tau1u.n_elem, tau1v.n_elem, M);
  mat Omega1, Omega2, out, out2;
  out.zeros(tau1u.n_elem, tau1v.n_elem);
  
  if (max(tau1u) == 0)
  {
    Omega1.eye(p, p);
  }
  else
  {
    if (d == 1)
      Omega1 = cubicmatrix(sx);
    else
      Omega1 = tpmatrix(sx);
  }
  if (max(tau1v) == 0)
  {
    Omega2.eye(q, q);
  }
  else
  {
    if (d == 1)
      Omega2 = cubicmatrix(sy);
    else
      Omega2 = tpmatrix(sy);
  }
  mat S12est = X.t() * Y / n;
  mat S12estt = S12est.t();
  mat Utempest, Vtempest;
  vec SPhiest;
  svd_econ(Utempest, SPhiest, Vtempest, S12est);
  // Need to set a reasonable upper bound
  double zeta = n * max(SPhiest);
  mat Gest(p + q, K), Cest(p + q, K), Gamma2est(p + q, K), Thetaest(p + q, p + q), matrixinv(p + q, p + q);
  Gest.rows(0, p - 1) = Utempest.cols(0, K - 1);
  Cest.rows(0, p - 1) = Utempest.cols(0, K - 1);
  Gest.rows(p, p + q - 1) = Vtempest.cols(0, K - 1);
  Cest.rows(p, p + q - 1) = Vtempest.cols(0, K - 1);
  Gamma2est.rows(0, p - 1) = S12est * Vtempest.cols(0, K - 1);
  Gamma2est.rows(p, p + q - 1) = S12estt * Utempest.cols(0, K - 1);
  Thetaest.submat(0, p, p - 1, p + q - 1) = S12est / 2;
  Thetaest.submat(p, 0, p + q - 1, p - 1) = S12estt / 2;
  mat S1221est = S12est * S12estt, S2112est = S12estt * S12est;
  cube S12train(p, q, M), S12traint(q, p, M), S12valid(p, q, M), S1221(p, p, M), S2112(q, q, M), Ucv(p, K, M), Vcv(q, K, M), Lmbd12cv(p, K, M), Lmbd22cv(q, K, M);
  double cvtau1u, cvtau1v, cvtau2u, cvtau2v, tempout;
  mat Ip;
  vec zetatemp(M);
  Ip.eye(p + q, p + q);
  
  if (tau1u.n_elem == 1 && tau1v.n_elem == 1)
  {
    cvtau1u = max(tau1u);
    cvtau1v = max(tau1v);
    Thetaest.submat(0, 0, p - 1, p - 1) = -cvtau1u * Omega1;
    Thetaest.submat(p, p, p + q - 1, p + q - 1) = -cvtau1v * Omega2;
    
    if (cvtau1u != 0)
    {
      mat M1 = inv_sympd(2 * zeta * Ip.submat(0, 0, p - 1, p - 1) - 2 * Thetaest.submat(0, 0, p - 1, p - 1) - S1221est / (2 * zeta));
      matrixinv.submat(0, 0, p - 1, p - 1) = M1;
      matrixinv.submat(0, p, p - 1, p + q - 1) = M1 * S12est / (2 * zeta);
      matrixinv.submat(p, 0, p + q - 1, p - 1) = matrixinv.submat(0, p, p - 1, p + q - 1).t();
      matrixinv.submat(p, p, p + q - 1, p + q - 1) = (S12estt * matrixinv.submat(0, p, p - 1, p + q - 1) + Ip.submat(p, p, p + q - 1, p + q - 1)) / (2 * zeta);
      spatmca_tau1(Gest, Cest, Gamma2est, matrixinv, p, q, zeta, maxit, tol);
    }
    if (cvtau1u == 0 && cvtau1v == 0)
    {
      matrixinv = 0.5 * inv_sympd(zeta * Ip - Thetaest);
    }
    else if (cvtau1v != 0)
    {
      if (cvtau1u == 0)
      {
        mat M2 = inv_sympd(2 * zeta * Ip.submat(p, p, p + q - 1, p + q - 1) - 2 * Thetaest.submat(p, p, p + q - 1, p + q - 1) - S2112est / (2 * zeta));
        matrixinv.submat(p, p, p + q - 1, p + q - 1) = M2;
        matrixinv.submat(p, 0, p + q - 1, p - 1) = M2 * S12estt / (2 * zeta);
        matrixinv.submat(0, p, p - 1, p + q - 1) = matrixinv.submat(p, 0, p + q - 1, p - 1).t();
        matrixinv.submat(0, 0, p - 1, p - 1) = (S12est * matrixinv.submat(p, 0, p + q - 1, p - 1) + Ip.submat(0, 0, p - 1, p - 1)) / (2 * zeta);
        spatmca_tau1(Gest, Cest, Gamma2est, matrixinv, p, q, zeta, maxit, tol);
      }
      else
      {
        matrixinv = 0.5 * inv_sympd(zeta * Ip - Thetaest);
        spatmca_tau1(Gest, Cest, Gamma2est, matrixinv, p, q, zeta, maxit, tol);
      }
    }
    if (tau2u.n_elem != 1 || tau2v.n_elem != 1)
    {
      spatmcacv_pp spatmcacv_pp(X, Y, K, Omega1, Omega2, tau1u, tau1v, nk, p, q, maxit,
                                tol, cv, S12train, S12traint, S12valid, S1221, S2112, zetatemp, Ucv, Vcv, Lmbd12cv, Lmbd22cv);
      RcppParallel::parallelFor(0, M, spatmcacv_pp);
    }
    out.zeros(1);
  }
  else
  {
    spatmcacv_p spatmcacv_p(X, Y, K, Omega1, Omega2, tau1u, tau1v, nk, p, q, maxit,
                            tol, cv, S12train, S12traint, S12valid, S1221, S2112, zetatemp, Ucv, Vcv, Lmbd12cv, Lmbd22cv);
    RcppParallel::parallelFor(0, M, spatmcacv_p);
    uword row1, col1;
    for (uword m = 0; m < M; m++)
      out += cv.slice(m) / M;
    tempout = out.min(row1, col1);
    cvtau1u = tau1u[row1];
    cvtau1v = tau1v[col1];
    Thetaest.submat(0, 0, p - 1, p - 1) = -cvtau1u * Omega1;
    Thetaest.submat(p, p, p + q - 1, p + q - 1) = -cvtau1v * Omega2;
    if (cvtau1u != 0)
    {
      mat M1 = inv_sympd(2 * zeta * Ip.submat(0, 0, p - 1, p - 1) - 2 * Thetaest.submat(0, 0, p - 1, p - 1) - S1221est / (2 * zeta));
      matrixinv.submat(0, 0, p - 1, p - 1) = M1;
      matrixinv.submat(0, p, p - 1, p + q - 1) = M1 * S12est / (2 * zeta);
      matrixinv.submat(p, 0, p + q - 1, p - 1) = matrixinv.submat(0, p, p - 1, p + q - 1).t();
      matrixinv.submat(p, p, p + q - 1, p + q - 1) = (S12estt * matrixinv.submat(0, p, p - 1, p + q - 1) + Ip.submat(p, p, p + q - 1, p + q - 1)) / (2 * zeta);
      spatmca_tau1(Gest, Cest, Gamma2est, matrixinv, p, q, zeta, maxit, tol);
    }
    if (cvtau1u == 0 && cvtau1v == 0)
    {
      matrixinv = 0.5 * inv_sympd(zeta * Ip - Thetaest);
    }
    else if (cvtau1v != 0)
    {
      if (cvtau1u == 0)
      {
        mat M2 = inv_sympd(2 * zeta * Ip.submat(p, p, p + q - 1, p + q - 1) - 2 * Thetaest.submat(p, p, p + q - 1, p + q - 1) - S2112est / (2 * zeta));
        matrixinv.submat(p, p, p + q - 1, p + q - 1) = M2;
        matrixinv.submat(p, 0, p + q - 1, p - 1) = M2 * S12estt / (2 * zeta);
        matrixinv.submat(0, p, p - 1, p + q - 1) = matrixinv.submat(p, 0, p + q - 1, p - 1).t();
        matrixinv.submat(0, 0, p - 1, p - 1) = (S12est * matrixinv.submat(p, 0, p + q - 1, p - 1) + Ip.submat(0, 0, p - 1, p - 1)) / (2 * zeta);
        spatmca_tau1(Gest, Cest, Gamma2est, matrixinv, p, q, zeta, maxit, tol);
      }
      else
      {
        matrixinv = inv_sympd(2 * zeta * Ip - 2 * Thetaest);
        spatmca_tau1(Gest, Cest, Gamma2est, matrixinv, p, q, zeta, maxit, tol);
      }
    }
  }
  if (tau2u.n_elem == 1 && tau2v.n_elem == 1)
  {
    cvtau2u = max(tau2u);
    cvtau2v = max(tau2v);
    if (cvtau2u > 0 || cvtau2v > 0)
    {
      mat Rest = Gest;
      mat Gamma1est = 0 * Gamma2est;
      for (uword i = 0; i < l2u.n_elem; i++)
        spatmca_tau2(Gest, Rest, Cest, Gamma1est, Gamma2est, matrixinv, l2u[i], 0, p, q, zeta, maxit, tol);
      for (uword i = 0; i < l2v.n_elem; i++)
        spatmca_tau2(Gest, Rest, Cest, Gamma1est, Gamma2est, matrixinv, max(l2u), l2v[i], p, q, zeta, maxit, tol);
    }
    out2 = tempout;
  }
  else
  {
    out2.zeros(tau2u.n_elem, tau2v.n_elem);
    cube cv2(tau2u.n_elem, tau2v.n_elem, M);
    spatmcacv_p2 spatmcacv_p2(S12train, S12traint, S12valid, S1221, S2112, zetatemp, Ucv, Vcv, Lmbd12cv, Lmbd22cv, Omega1, Omega2,
                              K, cvtau1u, cvtau1v, tau2u, tau2v, nk, p, q, maxit, tol, cv2);
    RcppParallel::parallelFor(0, M, spatmcacv_p2);
    uword row2, col2;
    for (uword m = 0; m < M; m++)
      out2 += cv2.slice(m) / M;
    out2.min(row2, col2);
    cvtau2u = tau2u[row2];
    cvtau2v = tau2v[col2];
    
    mat Rest = Gest;
    mat Gamma1est = 0 * Gamma2est;
    for (uword i = 0; i <= row2; i++)
      spatmca_tau2(Gest, Rest, Cest, Gamma1est, Gamma2est, matrixinv, tau2u[i], 0, p, q, zeta, maxit, tol);
    for (uword i = 0; i <= col2; i++)
      spatmca_tau2(Gest, Rest, Cest, Gamma1est, Gamma2est, matrixinv, cvtau2u, tau2v[i], p, q, zeta, maxit, tol);
  }
  vec zeros;
  zeros.zeros(K);
  vec D = max(zeros, diagvec(Gest.rows(0, p - 1).t() * S12est * Gest.rows(p, p + q - 1)));
  return List::create(Named("cv1") = out, Named("cv2") = out2, Named("Uest") = Gest.rows(0, p - 1), Named("Vest") = Gest.rows(p, p + q - 1), Named("Dest") = D, Named("cvtau1u") = cvtau1u, Named("cvtau2u") = cvtau2u, Named("cvtau1v") = cvtau1v, Named("cvtau2v") = cvtau2v);
}

struct spatmcacv_pall : public RcppParallel::Worker
{
  const mat &X;
  const mat &Y;
  const int K;
  const mat &Omega1;
  const mat &Omega2;
  const vec &tau1u;
  const vec &tau1v;
  const vec &tau2u;
  const vec &tau2v;
  const vec &nk;
  const int p1;
  const int p2;
  const int maxit;
  const double tol;
  cube &output;
  spatmcacv_pall(const mat &X, const mat &Y, const int K, const mat &Omega1, const mat &Omega2,
                 const vec &tau1u, const vec &tau1v, const vec &tau2u, const vec &tau2v,
                 const vec &nk, const int p1, const int p2, const int maxit, const double tol, cube &output) : X(X), Y(Y), K(K), Omega1(Omega1), Omega2(Omega2), tau1u(tau1u), tau1v(tau1v),
                 tau2u(tau2u), tau2v(tau2v), nk(nk), p1(p1), p2(p2), maxit(maxit), tol(tol), output(output) {}
  void operator()(std::size_t begin, std::size_t end)
  {
    mat Ip;
    Ip.eye(Y.n_cols, Y.n_cols);
    for (std::size_t k = begin; k < end; k++)
    {
      mat Uoldtemp, Voldtemp, Goldtemp, G(p1 + p2, K), R(p1 + p2, K), C(p1 + p2, K), Gamma1(p1 + p2, K),
      Gamma2(p1 + p2, K), Gold(p1 + p2, K), Rold(p1 + p2, K), Cold(p1 + p2, K), Gamma1old(p1 + p2, K), Gamma2old(p1 + p2, K),
      Gold0(p1 + p2, K), Rold0(p1 + p2, K), Cold0(p1 + p2, K), Gamma1old0(p1 + p2, K), Gamma2old0(p1 + p2, K),
      Gold2(p1 + p2, K), Rold2(p1 + p2, K), Cold2(p1 + p2, K), Gamma1old2(p1 + p2, K), Gamma2old2(p1 + p2, K),
      Gold3(p1 + p2, K), Rold3(p1 + p2, K), Cold3(p1 + p2, K), Gamma1old3(p1 + p2, K), Gamma2old3(p1 + p2, K),
      S12train(p1, p2), S12traint(p2, p1), S12valid(p1, p2), matrixinv(p1 + p2, p1 + p2), Ip;
      
      vec svdtemp, D;
      mat Ytrain = Y.rows(find(nk != (k + 1)));
      mat Yvalid = Y.rows(find(nk == (k + 1)));
      mat Xtrain = X.rows(find(nk != (k + 1)));
      mat Xvalid = X.rows(find(nk == (k + 1)));
      mat Theta(p1 + p2, p1 + p2);
      S12train = (Xtrain.t()) * Ytrain / Xtrain.n_rows;
      S12valid = Xvalid.t() * Yvalid / Xvalid.n_rows;
      svd_econ(Uoldtemp, svdtemp, Voldtemp, S12train);
      double zetatemp = 10 * max(svdtemp);
      Gold.rows(0, p1 - 1) = Uoldtemp.cols(0, K - 1);
      Rold.rows(0, p1 - 1) = Uoldtemp.cols(0, K - 1);
      Cold.rows(0, p1 - 1) = Uoldtemp.cols(0, K - 1);
      Gold.rows(p1, p1 + p2 - 1) = Voldtemp.cols(0, K - 1);
      Rold.rows(p1, p1 + p2 - 1) = Voldtemp.cols(0, K - 1);
      Cold.rows(p1, p1 + p2 - 1) = Voldtemp.cols(0, K - 1);
      Gamma2old.rows(0, p1 - 1) = S12train * Gold.rows(p1, p1 + p2 - 1);
      S12traint = S12train.t();
      Gamma2old.rows(p1, p1 + p2 - 1) = S12traint * Gold.rows(0, p1 - 1);
      
      Gold0 = Gold2 = Gold3 = Gold;
      Rold0 = Rold2 = Rold3 = Rold;
      Cold0 = Cold2 = Cold3 = Cold;
      Gamma2old0 = Gamma2old2 = Gamma2old3 = Gamma2old;
      
      Ip.eye(p1 + p2, p1 + p2);
      Gamma1old.zeros(p1 + p2, K);
      Gamma1old0 = Gamma1old2 = Gamma1old3 = Gamma1old;
      Theta.submat(0, p1, p1 - 1, p1 + p2 - 1) = S12train / 2;
      Theta.submat(p1, 0, p1 + p2 - 1, p1 - 1) = S12traint / 2;
      vec zero;
      zero.zeros(K);
      for (uword i = 0; i < tau1u.n_elem; i++)
      {
        G = Gold0;
        R = Rold0;
        C = Cold0;
        Gamma1 = Gamma1old0;
        Gamma2 = Gamma2old0;
        Theta.submat(0, 0, p1 - 1, p1 - 1) = -tau1u[i] * Omega1;
        for (uword j = 0; j < tau1v.n_elem; j++)
        {
          Theta.submat(p1, p1, p1 + p2 - 1, p1 + p2 - 1) = -tau1v[j] * Omega2;
          
          matrixinv = 0.5 * inv_sympd(zetatemp * Ip - Theta);
          vec zero;
          zero.zeros(K);
          for (uword l = 0; l < tau2u.n_elem; l++)
          {
            G = Gold3;
            R = Rold3;
            C = Cold3;
            Gamma1 = Gamma1old3;
            Gamma2 = Gamma2old3;
            for (uword m = 0; m < tau2v.n_elem; m++)
            {
              spatmca_tau2(G, R, C, Gamma1, Gamma2, matrixinv, tau2u[l], tau2v[m], p1, p2, zetatemp, maxit, tol);
              D = max(zero, diagvec(G.rows(0, p1 - 1).t() * S12train * G.rows(p1, p1 + p2 - 1)));
              output(tau2u.n_elem * i + l, tau2v.n_elem * j + m, k) =
                norm((S12valid - G.rows(0, p1 - 1) * diagmat(D) * G.rows(p1, p1 + p2 - 1).t()), "fro");
              if (l == 0 && m == 0)
              {
                Gold2 = G;
                Rold2 = R;
                Cold2 = C;
                Gamma1old2 = Gamma1;
                Gamma2old2 = Gamma2;
              }
              if (m == 0)
              {
                Gold3 = G;
                Rold3 = R;
                Cold3 = C;
                Gamma1old3 = Gamma1;
                Gamma2old3 = Gamma2;
              }
            }
          }
          
          if (j == 0)
          {
            Gold0 = Gold2;
            Rold0 = Rold2;
            Cold0 = Cold2;
            Gamma1old0 = Gamma1old2;
            Gamma2old0 = Gamma2old2;
          }
          Gold = Gold2;
          Rold = Rold2;
          Cold = Cold2;
          Gamma1old = Gamma1old2;
          Gamma2old = Gamma2old2;
        }
      }
    }
  }
};

//' Internal function: Overall M-fold Cross-validation of SpatMCA
//' @keywords internal
//' @param sxr A location matrix for a variable X
//' @param sxr A location matrix for a variable Y
//' @param Yr A data matrix of Y
//' @param Yr A data matrix of X
//' @param M A number of folds
//' @param K The number of estimated eigenfunctions
//' @param tau1ur A range of tau1u
//' @param tau2ur A range of tau2u
//' @param tau1vr A range of tau1v
//' @param tau2vr A range of tau2v
//' @param nkr A vector of fold numbers
//' @param maxit A maximum number of iteration
//' @param tol A tolerance rate
//' @param l2ur A given tau2u
//' @param l2vr A given tau2v
//' @return A list of selected parameters
// [[Rcpp::export]]
List spatmcacvall_rcpp(NumericMatrix sxr, NumericMatrix syr, NumericMatrix Xr, NumericMatrix Yr, int M, int K,
                       NumericVector tau1ur, NumericVector tau2ur, NumericVector tau1vr, NumericVector tau2vr,
                       NumericVector nkr, int maxit, double tol, NumericVector l2ur, NumericVector l2vr)
{
  int n = Yr.nrow(), p = Xr.ncol(), q = Yr.ncol(), d = sxr.ncol();
  mat X(Xr.begin(), n, p, false);
  mat Y(Yr.begin(), n, q, false);
  mat sx(sxr.begin(), p, d, false);
  mat sy(syr.begin(), q, d, false);
  colvec tau1u(tau1ur.begin(), tau1ur.size(), false);
  colvec tau2u(tau2ur.begin(), tau2ur.size(), false);
  colvec tau1v(tau1vr.begin(), tau1vr.size(), false);
  colvec tau2v(tau2vr.begin(), tau2vr.size(), false);
  colvec nk(nkr.begin(), nkr.size(), false);
  colvec l2u(l2ur.begin(), l2ur.size(), false);
  colvec l2v(l2vr.begin(), l2vr.size(), false);
  cube cv(tau1u.n_elem * tau2u.n_elem, tau1v.n_elem * tau2v.n_elem, M);
  mat Omega1, Omega2, out, out2;
  out.zeros(tau1u.n_elem, tau1v.n_elem);
  
  if (d == 1)
  {
    Omega1 = cubicmatrix(sx);
    Omega2 = cubicmatrix(sy);
  }
  else
  {
    Omega1 = tpmatrix(sx);
    Omega2 = tpmatrix(sy);
  }
  
  mat S12est = X.t() * Y / n;
  mat S12estt = S12est.t();
  mat Utempest, Vtempest;
  vec SPhiest;
  svd_econ(Utempest, SPhiest, Vtempest, S12est);
  double zeta = 10 * max(SPhiest);
  mat Gest(p + q, K), Cest(p + q, K), Gamma2est(p + q, K), Thetaest(p + q, p + q), matrixinv(p + q, p + q);
  Gest.rows(0, p - 1) = Utempest.cols(0, K - 1);
  Cest.rows(0, p - 1) = Utempest.cols(0, K - 1);
  Gest.rows(p, p + q - 1) = Vtempest.cols(0, K - 1);
  Cest.rows(p, p + q - 1) = Vtempest.cols(0, K - 1);
  Gamma2est.rows(0, p - 1) = S12est * Vtempest.cols(0, K - 1);
  Gamma2est.rows(p, p + q - 1) = S12estt * Utempest.cols(0, K - 1);
  Thetaest.submat(0, p, p - 1, p + q - 1) = S12est / 2;
  Thetaest.submat(p, 0, p + q - 1, p - 1) = S12estt / 2;
  mat S1221est = S12est * S12estt, S2112est = S12estt * S12est;
  double cvtau1u, cvtau1v, cvtau2u, cvtau2v;
  mat Ip;
  Ip.eye(p + q, p + q);
  
  out.zeros(tau1u.n_elem * tau2u.n_elem, tau1v.n_elem * tau2v.n_elem);
  spatmcacv_pall spatmcacv_pall(X, Y, K, Omega1, Omega2, tau1u,
                                tau1v, tau2u, tau2v, nk, p, q, maxit,
                                tol, cv);
  RcppParallel::parallelFor(0, M, spatmcacv_pall);
  
  uword row, col;
  for (uword m = 0; m < M; m++)
    out += cv.slice(m) / M;
  out.min(row, col);
  
  cvtau1u = tau1u[floor(row / tau2u.n_elem)];
  cvtau2u = tau2u[row % tau2u.n_elem];
  cvtau1v = tau1v[floor(col / tau2v.n_elem)];
  cvtau2v = tau2v[col % tau2v.n_elem];
  
  Thetaest.submat(0, 0, p - 1, p - 1) = -cvtau1u * Omega1;
  Thetaest.submat(p, p, p + q - 1, p + q - 1) = -cvtau1v * Omega2;
  
  matrixinv = inv_sympd(2 * zeta * Ip - 2 * Thetaest);
  spatmca_tau1(Gest, Cest, Gamma2est, matrixinv, p, q, zeta, maxit, tol);
  
  mat Rest = Gest;
  mat Gamma1est = 0 * Gamma2est;
  
  for (uword i = 0; i <= row % tau2u.n_elem; i++)
    spatmca_tau2(Gest, Rest, Cest, Gamma1est, Gamma2est, matrixinv, tau2u[i], 0, p, q, zeta, maxit, tol);
  for (uword i = 0; i <= col % tau2v.n_elem; i++)
    spatmca_tau2(Gest, Rest, Cest, Gamma1est, Gamma2est, matrixinv, cvtau2u, tau2v[i], p, q, zeta, maxit, tol);
  vec zeros;
  zeros.zeros(K);
  vec D = max(zeros, diagvec(Gest.rows(0, p - 1).t() * S12est * Gest.rows(p, p + q - 1)));
  return List::create(Named("cvall") = out, Named("Uest") = Gest.rows(0, p - 1), Named("Vest") = Gest.rows(p, p + q - 1), Named("Dest") = D, Named("cvtau1u") = cvtau1u, Named("cvtau2u") = cvtau2u, Named("cvtau1v") = cvtau1v, Named("cvtau2v") = cvtau2v);
}