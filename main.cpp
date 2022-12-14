#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>

#include "Timer.hpp"

using namespace Eigen;

static constexpr auto get_potential(double x, double y) { return 0.0; }

template <typename SparseMatrixType>
void set_matrix(SparseMatrixType& A, const VectorXd& xs, const VectorXd& ys,
                double dt) {
  auto dx = xs[1] - xs[0];
  auto dy = ys[1] - ys[0];
  std::complex<double> ct{0.0, dt};
  auto cx = ct / (4.0 * dx * dx);
  auto cy = ct / (4.0 * dy * dy);

  std::vector<Triplet<std::complex<double>>> nz;
  const auto get_id = [n = xs.size()](auto j, auto i) { return j + i * n; };

  for (auto i = 0; i < ys.size(); ++i) {
    for (auto j = 0; j < xs.size(); ++j) {
      auto v = get_potential(xs[j], ys[i]);
      auto cxy = 1.0 + ct / 2.0 * (v + 1.0 / (dx * dx) + 1.0 / (dy * dy));
      auto id = get_id(j, i);

      nz.emplace_back(id, id, cxy);
      if (j - 1 >= 0) nz.emplace_back(id, get_id(j - 1, i), -cx);
      if (j + 1 < xs.size()) nz.emplace_back(id, get_id(j + 1, i), -cx);
      if (i - 1 >= 0) nz.emplace_back(id, get_id(j, i - 1), -cy);
      if (i + 1 < ys.size()) nz.emplace_back(id, get_id(j, i + 1), -cy);
    }
  }

  A.setFromTriplets(nz.begin(), nz.end());
}

void set_rhs(VectorXcd& b, const VectorXcd& u0, const VectorXd& xs,
             const VectorXd& ys, double dt) {
  auto dx = xs[1] - xs[0];
  auto dy = ys[1] - ys[0];
  std::complex<double> ct{0.0, dt};
  auto cx = ct / (4.0 * dx * dx);
  auto cy = ct / (4.0 * dy * dy);

  const auto get_id = [n = xs.size()](auto j, auto i) { return j + i * n; };

  for (auto i = 0; i < ys.size(); ++i) {
    for (auto j = 0; j < xs.size(); ++j) {
      auto v = get_potential(xs[j], ys[i]);
      auto cxy = 1.0 - ct / 2.0 * (v + 1.0 / (dx * dx) + 1.0 / (dy * dy));
      auto id = get_id(j, i);

      b[id] = u0[id] * cxy;
      if (j - 1 >= 0) b[id] += u0[get_id(j - 1, i)] * cx;
      if (j + 1 < xs.size()) b[id] += u0[get_id(j + 1, i)] * cx;
      if (i - 1 >= 0) b[id] += u0[get_id(j, i - 1)] * cy;
      if (i + 1 < ys.size()) b[id] += u0[get_id(j, i + 1)] * cy;
    }
  }
}

struct GaussianWavePacket {
  double x0 = 0.5;
  double y0 = 0.5;
  double sigma = 1.0;
  double kx = 0.0;
  double ky = 0.0;
  double sigma2 = sigma * sigma;

  auto operator()(double x, double y) const {
    auto dx = x - x0;
    auto dy = y - y0;
    auto r = std::exp(-0.5 * (dx * dx + dy * dy) / sigma2) /
             (std::numbers::pi * sigma2);
    auto theta = kx * dx + ky * dy;
    return std::polar(r, theta);
  }
};

void set_initial_condition(VectorXcd& u0, const VectorXd& xs,
                           const VectorXd& ys) {
  auto kx = std::numbers::pi * 100;
  auto ky = std::numbers::pi * 50;
  GaussianWavePacket gwp{0.2, 0.4, 0.08, kx, ky};

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
  Timer timer(true);

  auto nx = 256;
  auto ny = nx;
  auto sz = nx * ny;
  VectorXd xs = VectorXd::LinSpaced(nx, 0.0, 1.0);
  VectorXd ys = xs;
  auto dx = xs[1] - xs[0];
  auto dy = dx;
  auto dt = dx * dy / 4;
  auto nt = 4096;
  auto fps = 64;

  SparseMatrix<std::complex<double>, RowMajor> A(sz, sz);
  VectorXcd b;
  VectorXcd u;
  VectorXcd u0;
  b.setZero(sz);
  u.setZero(sz);
  u0.setZero(sz);

  set_initial_condition(u0, xs, ys);
  set_matrix(A, xs, ys, dt);
  set_rhs(b, u0, xs, ys, dt);
  BiCGSTAB<decltype(A)> solver;
  // SparseLU<decltype(A)> solver;
  solver.compute(A);

  if (solver.info() != Success) {
    std::cerr << "Error\n";
    exit(EXIT_FAILURE);
  }

  std::ofstream ofs;
  std::cout << "Solving dofs: " << sz << ", threads: " << nbThreads() << '\n';

  auto t = 0.0;
  auto print_index = 0;

  for (auto i = 0; i < nt; ++i, t += dt) {
    if (i % fps == 0) {
      ofs.open("./data/u" + std::to_string(print_index++) + ".dat");
      output_solution(ofs, u0, xs, ys);
      ofs.close();
    }

    u = solver.solveWithGuess(b, u0);
    u0 = u;
    b.setZero();
    set_rhs(b, u0, xs, ys, dt);
  }

  ofs.open("./data/u" + std::to_string(print_index++) + ".dat");
  output_solution(ofs, u0, xs, ys);
  ofs.close();

  std::cout << "Iterations: " << nt << '\n';
  std::cout << "Printed: " << print_index << '\n';
  std::cout << "End time: " << t << '\n';
}