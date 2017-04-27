#include<omp.h>
#include<bits/stdc++.h>
using namespace std;

int *biggest_prime_divisor,
    *K,
    *L;

bool *k_function,
     *l_function;

int N;

int PRINT_EVERY;

void allocate() {
    biggest_prime_divisor = new int[N + 1];
    k_function = new bool[N + 1];
    l_function = new bool[N + 1];
    K = new int[N + 1];
}

void sieve() {
    // make parallel
    #pragma omp parallel for
    for (int i = 2; i <= N; i++) biggest_prime_divisor[i] = 0;

    int sqrt_of_N = sqrt(N) + 1;

    for (int i = 2; i <= sqrt_of_N; i++) if (biggest_prime_divisor[i] == 0) {
        // make parallel
        #pragma omp parallel for
        for (int j = i; j <= N; j += i) biggest_prime_divisor[j] = i;
    }
}

inline void calculate(int pos) {
    if (pos == 1) {
        k_function[pos] = l_function[pos] = 0;
        return;
    }
    if (pos == 2) {
        k_function[pos] = 0;
        l_function[pos] = 1;
        return;
    }

    int bpd = biggest_prime_divisor[pos];
    if (bpd == pos || bpd == 0) {
        k_function[pos] = l_function[pos] = 1;
        return;
    }

    k_function[pos] = k_function[bpd] ^ k_function[pos / bpd];
    l_function[pos] = l_function[bpd] ^ l_function[pos / bpd];
    return;
}

void summation() {
    L = biggest_prime_divisor; // substituting
    K[1] = L[1] = 1;
    for (int i = 2; i <= N; i++) {
        K[i] = K[i - 1] + (k_function[i] ? -1 : 1);
        L[i] = L[i - 1] + (l_function[i] ? -1 : 1);
    }
}

void print() {
    for (int i = 1; i <= N; i++) {
        cout << i << ":  " << K[i] << " " << L[i] << "   " << L[i] - K[i] << "\n";
    }
}

void print_for_a_graph() {
    std::ofstream plik;

    plik.open ("liczby.out", std::ofstream::out | std::ofstream::binary | std::ofstream::trunc);

    for (int i = 1; i <= N; i += PRINT_EVERY) {
        plik << K[i] << " " << L[i] << "\n";
    }
    plik.close();
}

void jebaj() {
    N = 200000000;
    allocate();
    sieve();
    int kon;
    for (int i = 1; i <= N; i *= 2) {
        kon = min(N, 2 * i - 1);
        #pragma omp parallel for
        for (int j = i; j <= kon; j++) {
            calculate(j);
        }
    }
    summation();
    // print();
    PRINT_EVERY = 500;
    print_for_a_graph();
}

int main() {
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    jebaj();
    end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "finished computation at " << std::ctime(&end_time)
             << "elapsed time: " << elapsed_seconds.count() << "s\n";
}
