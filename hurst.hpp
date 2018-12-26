#ifndef __HURST_H__
#define __HURST_H__

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <deque>
#include <vector>
#include <iostream>

class Hurst
{
public:
  Hurst(size_t cap=1000, size_t n_lags=20, size_t start=10) :
    cap_(cap), n_lags_(n_lags), start_(start) {
    ConstructLags();
    mean_.resize(cap);
    sub_.resize(cap);
    max_.resize(cap);
    min_.resize(cap);
    std_.resize(cap);
    cum_.resize(cap);
    log_RS_.resize(n_lags_);
  }

  bool compute() {
    for (size_t i = 0; i < n_lags_; i++) {
      size_t lag = lags_[i];
      size_t n = cap_ / lag;
      // compute means first
      for (size_t j = 0; j < lag; j++) {
        mean_[j] = 0;
        std_[j] = 0;
      }
      size_t index = 0;
      for (size_t j = 0; j < cap_; j++) {
        mean_[index] += ts_[j];
        index++;
        if (index == lag)
          index = 0;
      }
      for (size_t j = 0; j < lag; j++)
        mean_[j] /= n;
      // demean
      index = 0;
      for (size_t j = 0; j < cap_; j++) {
        sub_[j] = ts_[j] - mean_[index];
        cum_[j] = (j < lag) ? sub_[j] : sub_[j] + cum_[j - lag];
        // record max and min
        max_[index] = (j < lag) ? cum_[j] : std::max(max_[index], cum_[j]);
        min_[index] = (j < lag) ? cum_[j] : std::min(min_[index], cum_[j]);
        // compute std
        std_[index] += sub_[j] * sub_[j];

        index++;
        if (index == lag)
          index = 0;
      }
      double RS_sum = 0;
      for (size_t j = 0; j < lag; j++) {
        std_[j] = sqrt(std_[j] / n);
        RS_sum += (max_[j] - min_[j]) / std_[j];
      }
      RS_sum /= lag;
      std::cout << RS_sum << std::endl;
      log_RS_[i] = log(RS_sum);
    }

    // regression on log_RS_ and log_lags_
    // x: log_lags_
    // y: log_RS_
    double cov = Covariance(log_lags_, mean_log_lags_, log_RS_, Mean(log_RS_));
    std::cout << cov << " " << var_log_lags_ << " " << cov / var_log_lags_ << std::endl;
    h_ = cov / var_log_lags_;
    return h_;
  }

  void push_back(double x) {
    ts_.push_back(x);
    if (ts_.size() > cap_)
      ts_.pop_front();
  }

private:
  size_t cap_, n_lags_, start_;
  double var_log_lags_, mean_log_lags_, h_;
  std::vector<size_t> lags_;
  std::deque<double> ts_;
  std::vector<double> sub_, mean_, max_, min_, std_, cum_, log_RS_, log_lags_;

  void ConstructLags() {
    for (size_t n = start_; n <= cap_ / start_; n++) {
      if (cap_ % n == 0) {
        lags_.push_back(cap_ / n);
        log_lags_.push_back(log(n));
        if (lags_.size() > n_lags_)
          break;
      }
    }

    if (n_lags_ > lags_.size())
      n_lags_ = lags_.size();

    mean_log_lags_ = Mean(log_lags_);
    var_log_lags_ = Variance(log_lags_, mean_log_lags_);
  }

  double Mean(const std::vector<double>& array) {
    double sum = 0;
    for (auto& x : array)
      sum += x;
    return sum / array.size();
  }

  double Variance(const std::vector<double>& array, const double mean) {
    double sum = 0;
    for (auto& x : array)
      sum += (x - mean) * (x - mean);
    return sum / array.size();
  }

  double Covariance(const std::vector<double>& x, const double mean_x, 
                    const std::vector<double>& y, const double mean_y) {
    if (x.size() != y.size())
      return 0;
    double sum = 0;
    for (size_t i = 0; i < x.size(); i++)
      sum += (x[i] - mean_x) * (y[i] - mean_y);
    return sum / x.size();
  }
};

#endif // __HURST_H__