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

#include <vector>

#include <Rcpp.h>

using Rcpp::NumericVector;
using Rcpp::Function;
using Rcpp::NumericMatrix;
using Rcpp::Environment;
using Rcpp::stop;
using Rcpp::_;
using Rcpp::List;
using Rcpp::as;
using Rcpp::XPtr;

// #define COUT

class Dist {
 private:
  NumericVector m_a, m_b, m_loc, m_mul;
  Function pgamma;

 public:
  explicit Dist(NumericVector dst, int m = 1)
    : m_a(1), m_b(1), m_loc(1), m_mul(1),
      pgamma(Environment("package:stats")["pgamma"]) {
    if (dst["delta"] == 0) {
        // Dirac
        m_a[0] = INFINITY;
        m_b[0] = INFINITY;
        m_loc[0] = dst["mm"];
    } else if (dst["delta"] == INFINITY) {
        // "One" (heavyside because over positive support)
        m_a[0] = INFINITY;
        m_b[0] = 0;
        m_loc[0] = dst["mm"];
    } else {
        double d2 = dst["delta"]*dst["delta"];
        m_a[0] = 4 * dst["ss"]*dst["ss"] / d2;
        m_b[0] = 2 * dst["ss"] / (dst["mm"] * d2);
        m_loc[0] = dst["mm"] * (1 - 2 * dst["ss"]);
    }
    m_mul[0] = m;
  }

  Dist(NumericVector a, NumericVector b, NumericVector loc,
       NumericVector mul) :
    m_a(a), m_b(b), m_loc(loc), m_mul(mul),
    pgamma(Environment("package:stats")["pgamma"]) {
  }

  static Dist dirac(int m = 1, double n = 0) {
    NumericVector vec;
    vec["delta"] = 0;
    vec["mm"] = n;
    return Dist(vec, m);
  }

  static Dist one(int m = 1, double n = 0) {
    NumericVector vec;
    vec["delta"] = INFINITY;
    vec["mm"] = n;
    return Dist(vec, m);
  }

  NumericVector unpack(NumericVector times) {
    NumericVector res(times.length());
    for (int i = 0; i < m_a.length(); ++i) {
      if (std::isinf(m_mul[i])) stop("multiplicative factor not finite");
      if (std::isinf(m_a[i])) {  // Dirac or "one"
        if (std::isinf(m_b[i])) stop("Dirac: error");
        res = res + m_mul[i];
      } else {
        NumericVector p = pgamma(
          times - m_loc[i],
          _["shape"] = m_a[i],
          _["rate"] = m_b[i]);
        if (any(!is_finite(p)).is_true()) {
          stop("pgamma not finite");
        }
        res = res + m_mul[i] * p;
      }
    }
    return pmax(0, res);
  }

  Dist &operator+=(Dist rhs) {
    for (int i = 0; i < rhs.m_a.length(); ++i) {
      if (rhs.m_mul[i] != 0) {
        m_a.push_front(rhs.m_a[i]);
        m_b.push_front(rhs.m_b[i]);
        m_loc.push_front(rhs.m_loc[i]);
        m_mul.push_front(rhs.m_mul[i]);
      }
    }
    return (*this);
  }

  friend Dist operator+(Dist x, Dist y);
  friend Dist operator-(Dist x, Dist y);

  friend Dist operator*(double m, Dist x);
  friend Dist operator*(Dist x, Dist y);

  Dist nconv(int n) {
    if (n == 0) return dirac();
    else
      return Dist(n * m_a, m_b, n * m_loc, pow(m_mul, n));
  }

#ifdef COUT
    friend std::ostream &operator<<(std::ostream &os, Dist x);
#endif
};

Dist operator+(Dist x, Dist y) {
  NumericVector a, b, loc, mul;
  for (int i = 0; i < x.m_a.length(); ++i) {
    if (x.m_mul[i] != 0) {
      a.push_front(x.m_a[i]);
      b.push_front(x.m_b[i]);
      loc.push_front(x.m_loc[i]);
      mul.push_front(x.m_mul[i]);
    }
  }
  for (int i = 0; i < y.m_a.length(); ++i) {
    if (y.m_mul[i] != 0) {
      a.push_front(y.m_a[i]);
      b.push_front(y.m_b[i]);
      loc.push_front(y.m_loc[i]);
      mul.push_front(y.m_mul[i]);
    }
  }
  return Dist(a, b, loc, mul);
}

Dist operator-(Dist x, Dist y) {
  NumericVector a, b, loc, mul;
  for (int i = 0; i < x.m_a.length(); ++i) {
    if (x.m_mul[i] != 0) {
      a.push_front(x.m_a[i]);
      b.push_front(x.m_b[i]);
      loc.push_front(x.m_loc[i]);
      mul.push_front(x.m_mul[i]);
    }
  }
  for (int i = 0; i < y.m_a.length(); ++i) {
    if (y.m_mul[i] != 0) {
      a.push_front(y.m_a[i]);
      b.push_front(y.m_b[i]);
      loc.push_front(y.m_loc[i]);
      mul.push_front(-y.m_mul[i]);
    }
  }
  return Dist(a, b, loc, mul);
}

Dist operator*(double m, Dist x) {
  return Dist(x.m_a, x.m_b, x.m_loc, x.m_mul * m);
}

Dist operator*(Dist x, Dist y) {
  NumericVector a, b, loc, mul;
  for (int i = 0; i < x.m_a.length(); ++i) {
    for (int j = 0; j < y.m_a.length(); ++j) {
      mul.push_back(x.m_mul[i] * y.m_mul[j]);
      loc.push_back(x.m_loc[i] + y.m_loc[j]);
      if (std::isinf(x.m_a[i])) {
        a.push_back(y.m_a[j]);
        b.push_back(y.m_b[j]);
      } else if (std::isinf(y.m_a[j])) {
        a.push_back(x.m_a[i]);
        b.push_back(x.m_b[i]);
      } else {
        double A = x.m_a[i] / x.m_b[i] + y.m_a[j] / y.m_b[j];
        double B = x.m_a[i] / (x.m_b[i] * x.m_b[i]) +
            y.m_a[j] / (y.m_b[j] * y.m_b[j]);
        a.push_back(A*A / B);
        b.push_back(A / B);
      }
    }
  }
  return Dist(a, b, loc, mul);
}

#ifdef COUT
std::ostream &operator<<(std::ostream &os, Dist x) {
    os << "a: "; os << x.m_a;
    os << "\nb: "; os << x.m_b;
    os << "\nloc: "; os << x.m_loc;
    os << "\nmul: "; os << x.m_mul;
    return os;
}
#endif

NumericVector rowSums(NumericMatrix mat) {
  NumericVector ret(mat.nrow());
  for (int i = 0; i < mat.nrow(); ++i)
    ret[i] = sum(mat(i, _));
  return ret;
}

class ModelBranching {
 public:
  const int startPop;
  std::vector<std::vector<Dist> > lstLostPop;
  List lstLostPopDisp, lstLivePop, lstN, lstNLost, lstCumInflux;
  const NumericMatrix initial;
  const bool recDSM;
  const NumericVector times;
  const int Gm;

  ModelBranching(const List params,
                 int startPop,
                 int Gm,
                 NumericVector times,
                 const NumericMatrix initial,
                 bool recDSM)
      : startPop(startPop),
        initial(initial),
        recDSM(recDSM),
        times(times),
        Gm(Gm) {
    lstLivePop = List::create();
    lstN = List::create();
    lstNLost = List::create();
    lstCumInflux = List::create();

    if (startPop > 0) {
      model(params, 0);
    }
  }

  void model(const List params, int startPop);
};

void ModelBranching::model(const List params, int startPop) {
  for (int curPop = startPop; curPop < params.length(); ++curPop) {
    List it = params[curPop];

    // Probabilities:
    double p0 = it["p0"], res0 = it["res0"];
    double p = it["p"], res = it["res"];
    // Could be less than 0 owing to rounding errors:
    double d0 = fmax(0, 1 - p0 - res0), d = fmax(0, 1 - p - res);

    // Distributions:
    Dist G0(as<NumericVector>(it["g0"]));
    Dist G(as<NumericVector>(it["g"]));
    Dist F0(as<NumericVector>(it["f0"]));
    Dist F(as<NumericVector>(it["f"]));

    std::vector<Dist> H;

    if (curPop == 0) {
      for (int n = 0; n <= Gm; ++n)
        H.push_back(Dist::one(this->initial(curPop, n)));
    } else {
      NumericVector mrates0 = as<NumericVector>(it["mrates0"]);
      NumericVector mrates = as<NumericVector>(it["mrates"]);

      if (mrates0.length() > this->lstLostPop.size() ||
          mrates.length() > this->lstLostPop.size()) {
        stop("migration influxes format error");
      }

      // Get cumulative influx H
      Dist h = Dist::one(this->initial(curPop, 0));
      for (int i = 0; i < this->lstLostPop.size(); ++i) {
        h += mrates0[i] * this->lstLostPop[i][0];
      }
      H.push_back(h);
      if (Gm > 0) {
        for (int n = 1; n <= Gm; ++n) {
          h = Dist::one(this->initial(curPop, n));
          for (int i = 0; i < this->lstLostPop.size(); ++i) {
            h += mrates[i] * this->lstLostPop[i][n];
          }
          H.push_back(h);
        }
      }
    }

    // Get live and lost populations
    Dist opL0 = Dist::dirac() - p0 * F0 - d0 * G0,
        opL = Dist::dirac() - p * F - d * G,
        opM0 = 2 * p0 * F0,
        opM = 2 * p * F,
        opB0 = d0 * G0,
        opB = d * G;
    NumericMatrix matLivePop(times.length(), Gm + 1);

    matLivePop(_, 0) = (opL0 * H[0]).unpack(times);
    std::vector<Dist> lostPop;
    Dist h = opB0 * H[0];
    lostPop.push_back(h);
    NumericMatrix matLostPopDisp(times.length(), Gm + 1);
    if (this->recDSM) matLostPopDisp(_, 0) = h.unpack(times);
    for (int n = 1; n <= Gm; ++n) {
      Dist fact = opM0 * opM.nconv(n-1) * H[0];
      for (int k = 1; k <= n; ++k)
          fact += opM.nconv(n-k) * H[k];
      matLivePop(_, n) = (opL * fact).unpack(times);

      h = opB * fact;
      if (this->recDSM) {
        if (it.containsElementNamed("h")) {
          Dist opH = Dist::dirac() - Dist(as<NumericVector>(it["h"]));
          h = h * opH;
        }
        matLostPopDisp(_, n) = h.unpack(times);
      }
      lostPop.push_back(h);
    }

    NumericVector s = rowSums(matLivePop);
    this->lstN.push_back(s);

    NumericMatrix matProp(matLivePop.nrow(), matLivePop.ncol());
    for (int i = 0; i < s.length(); ++i) {
      if (s[i] > 0)
        matProp(i, _) = matLivePop(i, _) / s[i];
    }
    this->lstLivePop.push_back(matProp);

    this->lstLostPop.push_back(lostPop);

    if (this->recDSM) {
      s = rowSums(matLostPopDisp);
      this->lstNLost.push_back(s);
      NumericMatrix matPropLost(matLostPopDisp.nrow(),
                          matLostPopDisp.ncol());
      for (int i = 0; i < s.length(); ++i) {
        if (s[i] > 0)
          matPropLost(i, _) = matLostPopDisp(i, _) / s[i];
      }
      this->lstLostPopDisp.push_back(matPropLost);

      NumericMatrix matCumInflux(times.length(), Gm + 1);
      for (int n = 0; n <= Gm; ++n)
        matCumInflux(_, n) = H[n].unpack(times);
      this->lstCumInflux.push_back(matCumInflux);
    }
  }
}

// [[Rcpp::export]]
SEXP modelBranchingInit(const List params, int startPop,
                   int Gm, NumericVector times,
                           const NumericMatrix initial,
                           bool verbose = true) {
  return wrap(XPtr<ModelBranching>(
    new ModelBranching(params, startPop, Gm, times, initial, verbose),
    true));
}

// [[Rcpp::export]]
void modelBranchingRelease(SEXP obj) {
  as<XPtr<ModelBranching> >(obj).release();
}

// [[Rcpp::export]]
List modelBranching(const List params, SEXP obj) {
  ModelBranching objnew = ModelBranching(*as<XPtr<ModelBranching> >(obj));
  objnew.model(params, objnew.startPop);

  List ret;
  ret["live_pop"] = objnew.lstLivePop;
  ret["Ns"] = objnew.lstN;
  if (objnew.recDSM) {
    ret["lost_pop"] = objnew.lstLostPopDisp;
    ret["Ns_lost"] = objnew.lstNLost;
    ret["cum_influx"] = objnew.lstCumInflux;
  }
  return ret;
}
