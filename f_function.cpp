#include<omp.h>
#include<bits/stdc++.h>
#include <boost/algorithm/string.hpp>
using namespace std;
#define int long long

int *biggest_prime_divisor,
    *divided,
    *shared_F,
    *shared_f;

int N, n, threads;
vector<int> primes;

void allocate() {
    biggest_prime_divisor = new int[N + 1];
    memset(biggest_prime_divisor, 0, sizeof(int) * (N + 1));

    divided = new int[N + 1];

    shared_F = new int[(threads + 1) * (N + 1)];
    memset(shared_F, 0, sizeof(int) * (threads + 1) * (N + 1));

    shared_f = new int[(threads + 1) * (N + 1)];
}

void sieve() {
    for (int i = 2; i <= N; i++) if (biggest_prime_divisor[i] == 0) {
        primes.push_back(i);
        for (int j = i; j <= N; j += i) biggest_prime_divisor[j] = i;
    }
    for (int i = 2; i <= N; i++) divided[i] = i / biggest_prime_divisor[i];
    assert(N >= 2 * primes[n - 1] - 1);
    N = 2 * primes[n - 1] - 1;
    for (int i = primes[n - 1] + 1; i <= N; i++) {
        if (biggest_prime_divisor[i] == i) {
            divided[i] = 0;
        }
    }
}

inline void calculate(int mask) {
    int id = omp_get_thread_num();
    int *f  = shared_f + (id * (N + 1));
    int *F  = shared_F + (id * (N + 1));

    f[1] = 0;
    int pos = 0;
    while (pos < n) {
        if (mask & 1) {
            f[primes[pos]] = 0;
        } else {
            f[primes[pos]] = 1;
        }
        pos++;
        mask >>= 1;
    }

    for (int i = 4; i <= N; i++) {
        f[i] = (
            divided[i] == 0 ? 1 : f[biggest_prime_divisor[i]] ^ f[divided[i]]
        );
    }

    f[1] = 1;
    for (int i = 2; i <= N; i++) {
        f[i] = f[i - 1] + (f[i] ? -1 : 1);
        F[i] = min(F[i], f[i]);
    }
}

void join_F() {
    shared_F[1] = 1;
    for (int i = 2; i <= N; i++) {
        for (int thr = 0; thr < threads; thr++) {
            shared_F[i] = min(shared_F[i], shared_F[i + (N + 1) * thr]);
        }
    }
}

void print_for_a_graph() {
    std::ofstream plik;

    plik.open ("values_of_f_function.out", std::ofstream::out | std::ofstream::binary | std::ofstream::trunc);

    for (int i = 1; i <= N; i ++) {
        plik << shared_F[i] << "\n";
    }
    plik.close();
}

void jebaj() {
    allocate();
    sieve();
    int max_mask = (1LL << n) - 1LL;
    int mask;
    #pragma omp parallel for shared(biggest_prime_divisor, shared_f, primes, N) private(mask) num_threads(threads)
    for (mask = 0; mask <= max_mask; mask++) {
        calculate(mask);
    }

    join_F();
    print_for_a_graph();
}

#undef int
int main(int argc, char *argv[]) {
    for (int i = 1; i < argc; i++) {
        auto s = string(argv[i]);
        std::vector<std::string> strs;
        boost::split(strs, s, boost::is_any_of("="));
        if (strs[0] == "-n") n = stoi(strs[1]);
        if (strs[0] == "-N") N = stoi(strs[1]);
        if (strs[0] == "-threads") threads = stoi(strs[1]);

    }
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    jebaj();
    end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "finished computation at " << std::ctime(&end_time)
             << "elapsed time: " << elapsed_seconds.count() << "s\n";
}
