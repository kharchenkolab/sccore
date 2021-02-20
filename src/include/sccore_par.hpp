#ifndef _sccore_parallel_HPP
#define _sccore_parallel_HPP

#include <Rcpp.h>
#include <progress.hpp>

#include <mutex>
#include <thread>
#include <atomic>

namespace sccore {
using namespace Rcpp;

//' @md
//' @name ThreadProgress
//' @title Thread-safe Progress Bar
//' @description Alternative to RcppProgress that is safe for C++11 multithreading
//' @field ThreadProgress Constructor. Parameters:
//'   - `total_steps`: total number of steps for the progress
//'   - `verbose`: whether to show progress
//' @field increment Increment the counter.
//' @field isInterrupted was the process interrupted by user from R?
class ThreadProgress {
public:
  ThreadProgress(size_t total_steps, bool verbose)
    : total_steps(total_steps), verbose(verbose),
      master_thread_id(std::this_thread::get_id()), step(0), is_interrupted(false), n_printed(0)
  {
    if (verbose) {
      Rcout << "0%   10   20   30   40   50   60   70   80   90   100%\n";
      Rcout << "[----|----|----|----|----|----|----|----|----|----|\n";
    }
  }

  void increment() {
    ++this->step;
    if (!this->verbose || !this->isMasterThread())
      return;

    int need_printed = (this->total_chars * this->step) / this->total_steps;
    if (this->n_printed < need_printed) {
      Rcout << std::string(need_printed - this->n_printed, '*') << std::flush;
      this->n_printed = need_printed;
    }
  }

  bool isInterrupted() {
    this->updateInterrupted();
    return this->is_interrupted;
  }

  virtual ~ThreadProgress() {
    int need_printed = (this->total_chars * this->step) / this->total_steps;
    if (this->verbose && (this->n_printed < need_printed)) {
      Rcout << std::string(this->total_chars - this->n_printed, '*') << std::endl;
    }
  }

protected:
  bool isMasterThread() const {
    return std::this_thread::get_id() == this->master_thread_id;
  }

  bool updateInterrupted() {
    if (!this->isMasterThread())
      return false;

    if (this->is_interrupted)
      return false;

    this->is_interrupted = checkInterrupt();
    return true;
  }

private:
  const size_t total_steps;
  const bool verbose;
  const std::thread::id master_thread_id;
  const int total_chars = 51;
  std::atomic<unsigned long> step;
  std::atomic<bool> is_interrupted;
  int n_printed;
};


//' Run Task Parallel
//'
//' @param task task to run. It's a lambda that accepts ThreadProgress object and
//' is responsible for tracking the progress and user interruptions
//' @param n_cores number of cores
//' @param n_steps total number of steps, forwarded to ThreadProgress
//' @param verbose whether to show progress, forwarded to ThreadProgress
void runTaskParallel(const std::function<void(ThreadProgress&)> &task, int n_cores, int n_steps=1, bool verbose=false) {
  ThreadProgress p(n_steps, verbose);
  std::vector<std::thread> tasks;
  for (int thread_num = 0; thread_num < n_cores - 1; ++thread_num) {
    tasks.emplace_back(task, std::ref(p));
  }

  task(p);

  for (auto &t: tasks) {
    t.join();
  }

  if (p.isInterrupted())
    stop("Interrupted by user");
}

//' Run Task in a Parallel For Loop
//'
//' @param start int index from which the loop starts
//' @param end int index at which the loop starts
//' @param task task to run. It's a lambda that accepts integer index from the for loop
//' @param n_cores number of cores
//' @param n_steps total number of steps, forwarded to ThreadProgress
//' @param verbose whether to show progress, forwarded to ThreadProgress
void runTaskParallelFor(const int start, const int end, const std::function<void(int)> &task,
                        const int n_cores, const bool verbose=false) {
  std::atomic<int> glob_i(start);
  auto loop_task = [&task, &end, &glob_i](ThreadProgress &p) {
    while (true) {
      int i = glob_i++;
      if (i >= end)
        break;

      task(i);

      p.increment();
      if (p.isInterrupted())
        break;
    }
  };

  runTaskParallel(loop_task, n_cores, end - start, verbose);
}

}

#endif
