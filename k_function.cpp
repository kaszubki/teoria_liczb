#include<omp.h>
#include<bits/stdc++.h>
#include <boost/algorithm/string.hpp>
using namespace std;

int *biggest_prime_divisor,
    *K;
bool *k_function;

long long N;
int threads;
int print_every;

bool starting_primes[100];
string TITLE;

void allocate() {
    biggest_prime_divisor = new int[N + 1];
    fill_n(biggest_prime_divisor, (N + 1), 0);
    k_function = new bool[N + 1];
}

void sieve() {
    long long sqrt_of_N = sqrt(N) + 1;
    for (long long  i = 2; i <= sqrt_of_N; i++) {
        if (biggest_prime_divisor[i] == 0) {
            #pragma omp parallel for shared(biggest_prime_divisor) num_threads(threads)
            for (long long j = i; j <= N; j += i) {
                biggest_prime_divisor[j] = i;
            }
        }
    }
}

inline void calculate(long long pos) {
    if (pos == 1) {
        k_function[pos] = 0;
        return;
    }
    if (pos < 100 && starting_primes[pos] == true) {
        k_function[pos] = 0;
        return;
    }

    long long bpd = biggest_prime_divisor[pos];
    if (bpd == pos || bpd == 0) {
        k_function[pos] = 1;
        return;
    }

    k_function[pos] = k_function[bpd] ^ k_function[pos / bpd];
    return;
}

void summation() {
    K = biggest_prime_divisor; // substituting
    K[1] = 1;
    for (long long i = 2; i <= N; i++) {
        K[i] = K[i - 1] + (k_function[i] ? -1 : 1);
    }
}

void print() {
    std::ofstream plik;
    plik.open(TITLE+".out",
        std::ofstream::out | std::ofstream::binary | std::ofstream::trunc);
    for (long long i = 1; i <= N; i += print_every) {
        plik << K[i] << "\n";
    }
    plik.close();
}

void jebaj() {
    allocate();
    sieve();
    long long kon;
    for (long long i = 1; i <= N; i *= 2) {
        kon = min(N, 2 * i - 1);
        #pragma omp parallel for shared(biggest_prime_divisor, k_function) num_threads(threads)
        for (long long j = i; j <= kon; j++) {
            calculate(j);
        }
    }
    summation();
    print();
}

int main(int argc, char *argv[]) {
    TITLE = "values_of_summatory_function_of_L";
    for (int i = 1; i < argc; i++) {
        auto s = string(argv[i]);
        std::vector<std::string> strs;
        boost::split(strs, s, boost::is_any_of("="));
        if (strs[0] == "-N") N = stoll(strs[1]);
        if (strs[0] == "-threads") threads = stoi(strs[1]);
        if (strs[0] == "-print_every") print_every = stoi(strs[1]);
        if (strs[0] == "-primes") {
            TITLE += s;
            vector<string> strs_primes;
            boost::split(strs_primes, strs[1], boost::is_any_of(","));
            for (auto sp : strs_primes) {
                starting_primes[stoi(sp)] = 1;
            }
        }
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
