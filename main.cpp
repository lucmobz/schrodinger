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

using namespace Eigen;

static constexpr auto get_potential(double x, double y) { return 0.0; }

void set_matrix(SparseMatrix<std::complex<double>>& A, const VectorXd& xs,
                const VectorXd& ys, double dt) {
  using namespace std::complex_literals;

  auto dx = xs[1] - xs[0];
  auto dy = ys[1] - ys[0];
  auto ct = 1.0i * dt;
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
  using namespace std::complex_literals;

  auto dx = xs[1] - xs[0];
  auto dy = ys[1] - ys[0];
  auto ct = 1.0i * dt;
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
  double kx = 1.0;
  double ky = 1.0;
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
  GaussianWavePacket gwp{0.2, 0.4, 0.05, std::numbers::pi * 15.0, 0.0};

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
         << std::norm(u[id]) << ' '
         << u[id].real() << ' ' << u[id].imag() << '\n';
    }
    os << '\n';
  }
}

int main() {
  auto nx = 256;
  auto ny = nx;
  auto sz = nx * ny;
  VectorXd xs = VectorXd::LinSpaced(nx, 0.0, 1.0);
  VectorXd ys = xs;
  auto dx = xs[1] - xs[0];
  auto dy = dx;
  auto dt = dx * dy / 4;

  SparseMatrix<std::complex<double>> A(sz, sz);
  VectorXcd b;
  VectorXcd u;
  VectorXcd u0;
  b.setZero(sz);
  u.setZero(sz);
  u0.setZero(sz);

  set_initial_condition(u0, xs, ys);
  set_matrix(A, xs, ys, dt);
  set_rhs(b, u0, xs, ys, dt);
  SparseLU<decltype(A)> chol;
  chol.compute(A);

  if (chol.info() != Success) {
    std::cerr << "Error\n";
    exit(EXIT_FAILURE);
  }

  std::ofstream ofs;
  std::cout << "Solving dofs: " << sz << "\n";

  auto i = 0;
  auto print = 0;
  auto t = 0.0;
  for (; i < 2048; ++i, t += dt) {
    if (i % 8 == 0) {
      ofs.open("./data/u" + std::to_string(print++) + ".dat");
      output_solution(ofs, u0, xs, ys);
      ofs.close();
    }

    u = chol.solve(b);
    u0 = u;
    b.setZero();
    set_rhs(b, u0, xs, ys, dt);
  }

  ofs.open("./data/u" + std::to_string(print++) + ".dat");
  output_solution(ofs, u0, xs, ys);
  ofs.close();

  std::cout << i << '\n';
  std::cout << print << '\n';
  std::cout << t << '\n';
}