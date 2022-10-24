#include <chrono>
#include <string_view>

struct Timer {
  using Clock = std::chrono::high_resolution_clock;

  Timer(bool flag = false) : t0{Clock::now()}, print_flag{flag} {}

  ~Timer() {
    if (print_flag) print();
  }

  auto start() { t0 = Clock::now(); }

  auto stop() const {
    return std::chrono::duration<double>(Clock::now() - t0).count();
  }

  void print(std::string_view msg = "") const {
    if (!msg.empty()) std::cout << msg;
    std::cout << stop() << "s\n";
  }

  std::chrono::time_point<Clock> t0;
  bool print_flag = false;
};
