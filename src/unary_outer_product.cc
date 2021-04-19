#include "unary_outer_product.h"
#include <condition_variable>
#include <thread>
#include <iostream>


template <Mode mode>
std::vector<Share<mode>> populate_seeds(const ShareCSpan<mode>& x, std::size_t& missing) {
  const auto n = x.rows();

  // We maintain the seed buffer by putting seeds into appropriate tree locations.
  // The buffer only has to be large enough for the final layer as we only
  // store intermediate seeds temporarily.
  std::vector<Share<mode>> seeds(1 << n);

  // E keeps track of the missing tree node.
  missing = 0;

  // As a base case, we can derive the first two seeds from the possible labels for x[0].
  if constexpr (mode == Mode::G) {
    if (x[n-1].color()) {
      seeds[0] = x[n-1].H(); // S00
      seeds[1] = (~x[n-1]).H(); // S01
    } else {
      seeds[1] = x[n-1].H(); // S00
      seeds[0] = (~x[n-1]).H(); // S01
    }
  } else {
    seeds[!x[n-1].color()] = x[n-1].H();
  }
  missing |= x[n-1].color();
  ++Share<mode>::nonce;

  const auto one = Share<mode>::bit(true);
  const auto zero = Share<mode>::bit(false);

  // Now, iterate over the levels of the tree.
  for (std::size_t i = 1; i < n; ++i) {

    const auto key0 = x[n-i-1] ^ (x[n-i-1].color() ? zero : one);
    const auto key1 = key0 ^ one;

    if constexpr (mode == Mode::G) {
      // Maintain xor sums of all odd seeds/all even seeds
      Share<mode> odds = std::bitset<128> { 0 };
      Share<mode> evens = std::bitset<128> { 0 };
      // Work backwards across the level so as to not overwrite the parent seed
      // until it is no longer needed.
      for (int j = (1 << (i-1)); j >= 0; --j) {
        seeds[j*2 + 1] = seeds[j].H(0);
        seeds[j*2] = seeds[j].H(1);
        evens ^= seeds[j*2];
        odds ^= seeds[j*2 + 1];
      }

      (evens ^ (key0.H())).send();
      (odds ^ (key1.H())).send();

      const auto bit = x[n-i-1].color();
      missing = (missing << 1) | bit;
    } else {
      const auto g_evens = Share<mode>::recv();
      const auto g_odds = Share<mode>::recv();

      Share<mode> e_evens = std::bitset<128> { 0 };
      Share<mode> e_odds = std::bitset<128> { 0 };

      for (int j = (1 << (i-1)); j >= 0; --j) {
        if (j != missing) {
          seeds[j*2 + 1] = seeds[j].H(0);
          seeds[j*2] = seeds[j].H(1);
          e_evens ^= seeds[j*2];
          e_odds ^= seeds[j*2 + 1];
        }
      }

      // use the color of the `i`th share to figure out which element is missing at the next level.
      const auto bit = x[n-i-1].color();
      missing = (missing << 1) | bit;

      // assign the sibling of the missing node by (1) decrypting the appropriate row given by
      seeds[missing ^ 1] = x[n-i-1].H() ^ (bit ? (g_evens ^ e_evens) : (g_odds ^ e_odds));
    }

    ++Share<mode>::nonce;
  }

  return seeds;
}


// multithreading coordiation
std::vector<std::thread> threads;
std::vector<int> ready;
std::size_t njobs;
std::atomic<int> finished_job_counter;
std::condition_variable cv;
std::mutex mutex;
bool done;


struct GCtxt {
  std::size_t n;
  std::size_t l;
  std::span<const Share<Mode::G>> seeds;
  std::span<Share<Mode::G>> messages;
  const ShareSpan<Mode::G>* out;
  const ShareCSpan<Mode::G>* y;
  const Table* f;
};


GCtxt gctxt;

struct GJob {
  std::size_t start;
  std::size_t stop;

  void operator()() const {
    const auto n = gctxt.n;
    const auto l = gctxt.l;
    for (std::size_t j = start; j < stop; ++j) {
      Share<Mode::G> sum = std::bitset<128> { 0 };
      for (std::size_t i = 0; i < (1 << n); ++i) {
        const auto s = gctxt.seeds[i].H(Share<Mode::G>::nonce + (1 << n)*j + i);
        sum ^= s;
        std::size_t frow = (*gctxt.f)(i);
        for (std::size_t k = 0; k < l; ++k) {
          if (frow & (1 << k)) { (*gctxt.out)(k, j) ^= s; }
        }
      }
      sum ^= (*gctxt.y)[j];
      gctxt.messages[j] = sum;
    }
  }
};

std::vector<GJob> gjobs;


void gjob(std::size_t job_number) {
  ++finished_job_counter;
  while (true) {
    std::unique_lock<std::mutex> lock(mutex);
    cv.wait(lock, [job_number] { return ready[job_number] > 0 || done; });
    if (done) { return; }
    ready[job_number] = 0;
    lock.unlock();
    gjobs[job_number]();
    ++finished_job_counter;
  }
}


void initialize_gjobs() {
  done = false;
  /* njobs = std::thread::hardware_concurrency(); */
  njobs = 4;
  gjobs.resize(njobs);
  ready.resize(njobs);
  threads.resize(0);
  for (std::size_t i = 0; i < njobs; ++i) {
    threads.emplace_back([i] { gjob(i); });
  }

  int expected = njobs;
  while(!atomic_compare_exchange_strong(&finished_job_counter, &expected, 0)) {
    expected = njobs;
    // wait until all jobs start
  }
}

void finalize_jobs() {
  done = true;
  cv.notify_all();
  for (auto& th: threads) { th.join(); }
}



struct ECtxt {
  std::size_t missing;
  std::size_t n;
  std::size_t l;
  std::span<const Share<Mode::E>> seeds;
  std::span<Share<Mode::E>> messages;
  const ShareSpan<Mode::E>* out;
  const ShareCSpan<Mode::E>* y;
  const Table* f;
};


ECtxt ectxt;


struct EJob {
  std::size_t start;
  std::size_t stop;

  void operator()() const {
    const auto n = ectxt.n;
    const auto l = ectxt.l;
    for (std::size_t j = start; j < stop; ++j) {
      const Share<Mode::E> g_sum = ectxt.messages[j];
      Share<Mode::E> e_sum = std::bitset<128> { 0 };
      for (std::size_t i = 0; i < (1 << n); ++i) {
        if (i != ectxt.missing) {
          const auto s = ectxt.seeds[i].H(Share<Mode::E>::nonce + (1 << n)*j + i);
          e_sum ^= s;
          std::size_t frow = (*ectxt.f)(i);
          for (std::size_t k = 0; k < l; ++k) {
            if (frow & (1 << k)) { (*ectxt.out)(k, j) ^= s; }
          }
        }
      }
      const auto s = e_sum ^ g_sum ^ (*ectxt.y)[j];
      std::size_t frow = (*ectxt.f)(ectxt.missing);
      for (std::size_t k = 0; k < l; ++k) {
        if (frow & (1 << k)) { (*ectxt.out)(k, j) ^= s; }
      }
    }
  }
};

std::vector<EJob> ejobs;


void ejob(std::size_t job_number) {
  ++finished_job_counter;
  while (true) {
    std::unique_lock<std::mutex> lock(mutex);
    cv.wait(lock, [job_number] { return ready[job_number] > 0 || done; });
    if (done) { return; }
    ready[job_number] = 0;
    lock.unlock();
    ejobs[job_number]();
    ++finished_job_counter;
  }
}


void initialize_ejobs() {
  done = false;
  /* njobs = std::thread::hardware_concurrency(); */
  njobs = 4;
  ejobs.resize(njobs);
  ready.resize(njobs);
  threads.resize(0);
  for (std::size_t i = 0; i < njobs; ++i) {
    threads.emplace_back([i] { ejob(i); });
  }

  int expected = njobs;
  while(!atomic_compare_exchange_strong(&finished_job_counter, &expected, 0)) {
    expected = njobs;
    // wait until all jobs start
  }
}


template <Mode mode>
void unary_outer_product(
    const Table& f, const ShareCSpan<mode>& x, const ShareCSpan<mode>& y, const ShareSpan<mode>& out) {

  assert(x.cols() == 1);
  assert(y.cols() == 1);

  const auto n = x.rows();
  const auto m = y.rows();
  const auto l = out.rows();
  assert(out.cols() == m);

  std::size_t missing = 0;
  const auto seeds = populate_seeds<mode>(x, missing);

  // Now we are ready to compute the outer product.
  // For each share (B, B + bDelta)
  // G sends the sum (XOR_i A_i) + B, which allows E to obtain A_{x + gamma} + bDelta
  if constexpr (mode == Mode::G) {
    std::vector<Share<mode>> messages(m);
    gctxt = { n, l, seeds, messages, &out, &y, &f };

    {
      std::unique_lock<std::mutex> lock(mutex);
      for (std::size_t jb = 0; jb < njobs; ++jb) {
        GJob job { jb*(m/njobs), std::min((jb+1)*(m/njobs), m) };
        gjobs[jb] = job;
        ready[jb] = 1;
      }
    }

    // dispatch jobs
    cv.notify_all();
    int expected = njobs;
    while(!atomic_compare_exchange_strong(&finished_job_counter, &expected, 0)) {
      expected = njobs;
      // wait until all jobs finish
    }

    for (auto& m: messages) { m.send(); }
  } else {
    std::vector<Share<mode>> messages(m);
    for (auto& m: messages) { m = Share<mode>::recv(); }
    ectxt = { missing, n, l, seeds, messages, &out, &y, &f };

    {
      std::unique_lock<std::mutex> lock(mutex);
      for (std::size_t jb = 0; jb < njobs; ++jb) {
        EJob job { jb*(m/njobs), std::min((jb+1)*(m/njobs), m) };
        ejobs[jb] = job;
        ready[jb] = 1;
      }
    }

    // dispatch jobs
    cv.notify_all();
    int expected = njobs;
    while(!atomic_compare_exchange_strong(&finished_job_counter, &expected, 0)) {
      expected = njobs;
      // wait until all jobs finish
    }
  }
  Share<mode>::nonce += (1<<n)*m;
}



template void unary_outer_product(
    const Table&,
    const ShareCSpan<Mode::G>&,
    const ShareCSpan<Mode::G>&,
    const ShareSpan<Mode::G>&);
template void unary_outer_product(
    const Table&,
    const ShareCSpan<Mode::E>&,
    const ShareCSpan<Mode::E>&,
    const ShareSpan<Mode::E>&);
