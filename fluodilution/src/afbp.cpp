// Copyright (c) 2015-2018 Hadrien Chauvin
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//////////////////////////////////////////

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
//#include <RcppParallel.h>
#include <string>

using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::as;
using Rcpp::List;
using Rcpp::max;

double half_sqrt_2pi = sqrt(2 * M_PI) * 0.5;

struct AfBpControl {
  int kSmallMax;
  double tol;
  double sumTol;
  double pRightTail;
  double mmTransition;
  double minExp;
  double kernelScale;
  bool parallelize;
  int nThreads;
  int minParallelGrain;
};

const AfBpControl defaultAfBpControl = {
  // From R core: 30 is somewhat arbitrary: it is on the *safe* side:
  // both speed and precision are clearly improved for k < 30.
  30,
  1e-10,
  1e-8,
  1e-9,
  7,
  -700,
  4,
  false, // true,
  10,
  100
};

inline double fastPow(double base, int exp) {
  double result = 1;
  while (exp) {
    if (exp & 1) {
      result *= base;
    }
    exp >>= 1;
    base *= base;
  }
  return result;
}

inline double choose(double n, double k) {
  if (n-k < k && n >= 0) k = n-k; /* <- Symmetry */
  if (k == 0) return 1.;
  /* else: k >= 1 */
  double r = n;
  for (int j = 2; j <= k; j++) {
    r *= (n-j+1)/j;
  }
  return r;
}

arma::mat probas(17, 10000);

const arma::mat *fdBinomialPartitioning(
  double mm,
  double s,
  double mgen,
  const AfBpControl ctl
) {
  // Maximum number of molecules for which to find probability
  double nMax = R::qlnorm(
    1 - ctl.pRightTail,  // p
    mm,                  // meanlog
    s,                   // sdlog
    1,                   // lower.tail
    0);                  // log.p

  // This will contain the result.
  //arma::mat probas(mgen + 1, nMax, arma::fill::zeros);
  probas.zeros(mgen + 1, nMax);

  // First iterations: log-norm centered on mm/2^gen with SD s
  int m = 0;
  do {
    double prev = 0;
    for (int i = 1; i < probas.n_cols; ++i) {
      double cur = R::plnorm(i, mm, s, 1, 0);
      double diff = cur - prev;
      probas(m, i) = diff;
      prev = cur;
    }
    mm -= log(2);
    ++m;
  } while (mm > ctl.mmTransition && m < probas.n_rows);

  for (; m < probas.n_rows; ++m) {
    // Only take into account in the sum probas that are greater
    // than a given tolerance.  As the distribution is unimodal,
    // we will thus have cutoff levels on the left and on the right.
    int minn = -1;  // Left cutoff, inclusive
    for (int i = 0; i < probas.n_cols; ++i) {
      if (probas(m - 1, i) >= ctl.tol) {
        minn = i;
        break;
      }
    }

    int maxx = -1;  // Right cutoff, exclusive
    if (minn >= 0) {
      for (int i = probas.n_cols - 1; i >= minn; --i) {
        if (probas(m - 1, i) >= ctl.tol) {
          maxx = i + 1;
          break;
        }
      }
    }

    if (minn == -1 || maxx == -1) {
      // Almost only 0 molecules, all the way down the bottom of the
      // `probas` matrix.
      for (int k = m; k < probas.n_rows; ++k) {
        probas(k, 0) = 1;
      }
      break;
    }

    double aSeed = fastPow(2.0, minn);

    int rowSum = 0;
    for (int i = 0; i < maxx; ++i) {
      int curMinn = minn > i ? minn : i;  // max(minn, i)
      double sum = 0;
      double a = aSeed;
      if (i > minn) {
        a = aSeed *= 2;
      }
      bool notFinite = false;
      for (int j = curMinn; j < maxx; ++j) {
        if ((i < ctl.kSmallMax || (j - i) < ctl.kSmallMax) && !notFinite) {
          double cur = probas(m - 1, j) * choose(j, i) / a;
          if (!std::isfinite(cur)) {
            notFinite = true;
            --j;
            continue;
          } else {
            sum += cur;
            a *= 2;
          }
        } else {
          // we use the normality approximation
          double b = i - j * 0.5;
          double operand = -b * b / (j / 2);
          if (operand < ctl.minExp) {
            break;
          }
          sum += probas(m - 1, j) *
            exp(operand) /
            (half_sqrt_2pi * sqrt(j));
        }
      }
      probas(m, i) = sum;
      rowSum += sum;
    }
  }

  return &probas;
}

struct AfBpParams {
  double ftorExp;
  double sdaf;
  const arma::mat &calc;
  NumericVector a;
  NumericVector b;
  NumericVector gen;
};

class AfBpFunctor {
 private:
  int m_i;
  int m_j;
  const AfBpParams m_params;
  double m_middle;
  double m_diff;
  double m_dnormCoef;

 public:
  AfBpFunctor(int i, int j, const AfBpParams params)
      : m_i(i), m_j(j), m_params(params) {
    double lower = m_params.a[m_i];
    double upper = m_params.b[m_i];
    m_middle = (upper + lower) / 2;
    m_diff = upper - lower;
    m_dnormCoef = M_1_SQRT_2PI / m_params.sdaf;
  }

  inline double operator()(const double x) const {
    int n = floor((m_middle - x) * m_params.ftorExp);
    if (n >= 0 && n < m_params.calc.n_cols) {
      return m_params.calc(m_params.gen[m_j], n) * m_params.ftorExp * m_diff;
    }
    return 0;
  }
};

class GaussKernel {
 private:
  arma::vec m_kernel;

 public:
  double m_lowerBound;
  double m_step;
  int m_nSamples;

  GaussKernel(int samples, double sigma)
      : m_nSamples(samples), m_kernel(samples) {
    m_lowerBound = -2 * sigma;
    double upperBound = 2 * sigma;

    m_step = (upperBound - m_lowerBound) / samples;
    double step_sigma = m_step / sigma;
    double cst = M_1_SQRT_2PI / sigma;

    double a = m_lowerBound / sigma;
    for (int i = 0; i < samples; ++i) {
      m_kernel[i] = cst * exp(-0.5 * a * a);
      a += step_sigma;
    }
  }

  inline double operator[](int i) const {
    return m_kernel[i];
  }
};

class AfBpWorker //: public RcppParallel::Worker
{
 private:
  const AfBpParams m_params;
  //RcppParallel::RMatrix<double> m_ans;
  arma::mat &m_ans;
  const GaussKernel &m_kernel;

 public:
  AfBpWorker(
    const AfBpParams params,
    const GaussKernel &kernel,
    arma::mat &ans)
      : m_params(params), m_ans(ans), m_kernel(kernel) {}

  void operator()(std::size_t begin, std::size_t end) {
    double lowerBound = -2 * m_params.sdaf;
    double upperBound = 2 * m_params.sdaf;

    std::size_t i = begin / m_ans.n_cols;
    std::size_t j = begin % m_ans.n_cols;

    std::size_t cur = begin;
    while (cur != end) {
      AfBpFunctor f(i, j, m_params);
      double errEst;
      int errCode;

      double curAns = 0;
      double x = m_kernel.m_lowerBound;
      curAns += m_kernel[0] * f(x);
      x += m_kernel.m_step;
      for (int k = 1; k < m_kernel.m_nSamples - 1; ++k) {
        curAns += 2.0 * m_kernel[k] * f(x);
        x += m_kernel.m_step;
      }
      curAns += m_kernel[m_kernel.m_nSamples - 1] * f(x);
      curAns *= m_kernel.m_step / 2.0;
      m_ans(i, j) = curAns;

      if (++j == m_ans.n_cols) {
        j = 0;
        ++i;
      }
      ++cur;
    }
  }
};

arma::mat ans(10000, 17);

// [[Rcpp::export]]
SEXP afBpDiffCdf(
  NumericVector a,
  NumericVector b,
  NumericVector gen,
  List paramsList
) {
  double ftor = as<double>(paramsList["ftor"]);

  AfBpControl ctl = defaultAfBpControl;

  const arma::mat *calc = fdBinomialPartitioning(
    as<double>(paramsList["m0"]) - ftor,
    paramsList["sd0"],
    Rcpp::max(gen),
    ctl);

  AfBpParams params = {
    exp(-ftor),
    as<double>(paramsList["sdaf"]),
    *calc,
    a,
    b,
    gen
  };

  GaussKernel kernel(ctl.kernelScale * params.sdaf, params.sdaf);
  //NumericMatrix ans(a.length(), gen.length());
  ans.zeros(a.length(), gen.length());
  AfBpWorker worker(params, kernel, ans);
  // if (ctl.parallelize) {
  //   int grain = ans.length() / ctl.nThreads;
  //   if (grain < ctl.minParallelGrain) grain = ctl.minParallelGrain;
  //   RcppParallel::parallelFor(0, ans.length(), worker, grain);
  // } else {
    worker(0, ans.n_elem);
  // }

  return Rcpp::wrap(ans);
}


