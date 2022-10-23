#include <Eigen/Dense>
#include <fstream>
#include <Eigen/Sparse>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>

using namespace Eigen;

static constexpr auto get_potential(double x, double y) { return 0.0; }

void set_matrix(SparseMatrix<std::complex<double>>& A, const VectorXd& xs,
                const VectorXd& ys, double dt) {
  auto dx = xs[1] - xs[0];
  auto dy = ys[1] - ys[0];
  std::complex<double> ct{0.0, dt};
  auto cx = 0.5 * ct / (dx * dx);
  auto cy = 0.5 * ct / (dy * dy);
  auto cxy = 1.0 + ct * (1.0 / (dx * dx) + 1.0 / (dy * dy));

  std::vector<Triplet<std::complex<double>>> nnz;
  const auto get_id = [n = xs.size()](auto j, auto i) { return j + i * n; };

  for (auto i = 0; i < ys.size(); ++i) {
    for (auto j = 0; j < xs.size(); ++j) {
      auto v = get_potential(xs[j], ys[i]);
      auto id = get_id(j, i);

      nnz.emplace_back(id, id, cxy - 0.5 * ct * v);
      if (j - 1 >= 0) nnz.emplace_back(id, get_id(j - 1, i), -cx);
      if (j + 1 < xs.size()) nnz.emplace_back(id, get_id(j + 1, i), -cx);
      if (i - 1 >= 0) nnz.emplace_back(id, get_id(j, i - 1), -cy);
      if (i + 1 < ys.size()) nnz.emplace_back(id, get_id(j, i + 1), -cy);
    }
  }

  A.setFromTriplets(nnz.begin(), nnz.end());
}

void set_rhs(VectorXcd& b, const VectorXcd& u0, const VectorXd& xs,
             const VectorXd& ys, double dt) {
  auto dx = xs[1] - xs[0];
  auto dy = ys[1] - ys[0];
  std::complex<double> ct{0.0, dt};
  auto cx = 0.5 * ct / (dx * dx);
  auto cy = 0.5 * ct / (dy * dy);
  auto cxy = 1.0 - ct * (1.0 / (dx * dx) + 1.0 / (dy * dy));

  const auto get_id = [n = xs.size()](auto j, auto i) { return j + i * n; };

  for (auto i = 0; i < ys.size(); ++i) {
    for (auto j = 0; j < xs.size(); ++j) {
      auto v = get_potential(xs[j], ys[i]);
      auto id = get_id(j, i);

      b[id] = u0[id] * (cxy + 0.5 * ct * v);
      if (j - 1 >= 0) b[id] += u0[get_id(j - 1, i)] * cx;
      if (j + 1 < xs.size()) b[id] += u0[get_id(j + 1, i)] * cx;
      if (i - 1 >= 0) b[id] += u0[get_id(j, i - 1)] * cy;
      if (i + 1 < ys.size()) b[id] = u0[get_id(j, i + 1)] * cy;
    }
  }
}

struct GaussianWavePacket {
  double x0 = 0.0;
  double y0 = 0.0;
  double sigma = 1.0;
  double kx = 0.0;
  double ky = 0.0;

  auto operator()(double x, double y) const {
    auto dx = x - x0;
    auto dy = y - y0;
    auto r = std::exp(-0.5 * (dx * dx + dy * dy) / sigma / sigma) /
             (std::numbers::pi * sigma * sigma);
    auto theta = -(kx * dx + ky * dy);
    return std::polar(r, theta);
  }
};

void set_initial_condition(VectorXcd& u0, const VectorXd& xs,
                           const VectorXd& ys) {
  GaussianWavePacket gwp{0.2, 0.5, 0.05, 0.0, 0.0};

  for (auto i = 0; i < ys.size(); ++i)
    for (auto j = 0; j < xs.size(); ++j)
      u0[j + i * xs.size()] = gwp(xs[j], ys[i]);
}

void output_solution(std::ostream& os, const VectorXcd& u, const VectorXd& xs,
                     const VectorXd& ys) {
  for (auto i = 0; i < ys.size(); ++i) {
    for (auto j = 0; j < xs.size(); ++j) {
      auto id = j + i * xs.size();
      os << std::setprecision(15) << xs[j] << ' ' << ys[i] << ' '
         << std::norm(u[id]) << ' ' << u[id].real() << ' ' << u[id].imag()
         << '\n';
    }
    os << '\n';
  }
}

int main() {
  auto nx = 512;
  auto ny = nx;
  auto sz = nx * ny;

  VectorXd xs = VectorXd::LinSpaced(nx, 0.0, 1.0);
  VectorXd ys = VectorXd::LinSpaced(ny, 0.0, 1.0);

  auto dx = xs[1] - xs[0];
  auto dy = ys[1] - ys[0];
  auto dt = nx * ny * 1.0e-2;

  SparseMatrix<std::complex<double>> A(sz, sz);
  VectorXcd b;
  VectorXcd u;
  VectorXcd u0;
  b.setZero(sz);
  u.setZero(sz);
  u0.setZero(sz);

  set_matrix(A, xs, ys, dt);
  set_rhs(b, u0, xs, ys, dt);
  set_initial_condition(u0, xs, ys);
  SimplicialCholesky<decltype(A)> chol(A);

  std::ofstream ofs{"./data/u0.dat"};
  output_solution(ofs, u0, xs, ys);
  ofs.close();

  auto t = 0.0;
  for (auto i = 1; i <= 128; ++i, t += dt) {
    u = chol.solve(b);
    u0 = u;
    set_rhs(b, u0, xs, ys, dt);

    std::cout << u0.norm() << '\n';

    ofs.open("./data/u" + std::to_string(i) + ".dat");
    output_solution(ofs, u0, xs, ys);
    ofs.close();

    if (i == 2) break;
  }
}