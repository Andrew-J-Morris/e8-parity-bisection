/**
 * Cumulative lattice points in E_8 with norm squared <= target_r2.
 * Parameterization: Set target_r2 in main() to evaluate a single domain bound.
 */

#include <iostream>
#include <vector>
#include <cstdint>
#include <chrono>

using namespace std;
using namespace std::chrono;

// Bitwise integer square root executing entirely on the integer ALU
static inline uint64_t pure_isqrt(uint64_t value) {
    uint64_t res = 0;
    uint64_t bit = 1ULL << 62;
    while (bit > value) bit >>= 2;
    while (bit != 0) {
        if (value >= res + bit) {
            value -= res + bit;
            res = (res >> 1) + bit;
        } else {
            res >>= 1;
        }
        bit >>= 2;
    }
    return res;
}

// Phase 1: 2D Base Profile Construction (O(r^2))
void build_E8_2D_profiles(uint64_t scaled_r2,
                          vector<uint64_t>& f2_even_s0, vector<uint64_t>& f2_even_s2,
                          vector<uint64_t>& f2_odd_s0,  vector<uint64_t>& f2_odd_s2) {
    int64_t max_coord = static_cast<int64_t>(pure_isqrt(scaled_r2));

    // All-Even 2D Coordinate Sweep (Steps by 2)
    int64_t min_even = (max_coord % 2 != 0) ? -(max_coord - 1) : -max_coord;
    for (int64_t x = min_even; x <= max_coord; x += 2) {
        uint64_t x2 = static_cast<uint64_t>(x * x);
        for (int64_t y = min_even; y <= max_coord; y += 2) {
            uint64_t m = x2 + static_cast<uint64_t>(y * y);
            if (m <= scaled_r2) {
                if ((x + y) % 4 == 0) f2_even_s0[m]++;
                else f2_even_s2[m]++;
            }
        }
    }

    // All-Odd 2D Coordinate Sweep (Steps by 2)
    int64_t min_odd = (max_coord % 2 == 0) ? -(max_coord - 1) : -max_coord;
    for (int64_t x = min_odd; x <= max_coord; x += 2) {
        uint64_t x2 = static_cast<uint64_t>(x * x);
        for (int64_t y = min_odd; y <= max_coord; y += 2) {
            uint64_t m = x2 + static_cast<uint64_t>(y * y);
            if (m <= scaled_r2) {
                if ((x + y) % 4 == 0) f2_odd_s0[m]++;
                else f2_odd_s2[m]++;
            }
        }
    }
}

// Phase 2: 4D Intermediate Convolutions with Sparse-Index Acceleration
void build_E8_4D_profiles(uint64_t scaled_r2,
                          const vector<uint64_t>& f2_even_s0, const vector<uint64_t>& f2_even_s2,
                          const vector<uint64_t>& f2_odd_s0,  const vector<uint64_t>& f2_odd_s2,
                          vector<uint64_t>& f4_even_s0, vector<uint64_t>& f4_even_s2,
                          vector<uint64_t>& f4_odd_s0,  vector<uint64_t>& f4_odd_s2) {
    // Collect non-zero entries for sparse traversal
    vector<size_t> nz_even, nz_odd;
    nz_even.reserve(scaled_r2 / 8);
    nz_odd.reserve(scaled_r2 / 8);
    for (size_t m = 0; m <= scaled_r2; ++m) {
        if (f2_even_s0[m] != 0 || f2_even_s2[m] != 0) nz_even.push_back(m);
        if (f2_odd_s0[m] != 0  || f2_odd_s2[m] != 0)  nz_odd.push_back(m);
    }

    // Even + Even Crosses
    for (size_t i = 0; i < nz_even.size(); ++i) {
        size_t u = nz_even[i];
        uint64_t u_s0 = f2_even_s0[u], u_s2 = f2_even_s2[u];
        for (size_t j = 0; j < nz_even.size(); ++j) {
            size_t v = nz_even[j];
            if (u + v > scaled_r2) break; // Exact monotonic bound
            uint64_t v_s0 = f2_even_s0[v], v_s2 = f2_even_s2[v];

            if (u_s0 && v_s0) f4_even_s0[u + v] += u_s0 * v_s0;
            if (u_s2 && v_s2) f4_even_s0[u + v] += u_s2 * v_s2;
            if (u_s0 && v_s2) f4_even_s2[u + v] += u_s0 * v_s2;
            if (u_s2 && v_s0) f4_even_s2[u + v] += u_s2 * v_s0;
        }
    }

    // Odd + Odd Crosses
    for (size_t i = 0; i < nz_odd.size(); ++i) {
        size_t u = nz_odd[i];
        uint64_t u_s0 = f2_odd_s0[u], u_s2 = f2_odd_s2[u];
        for (size_t j = 0; j < nz_odd.size(); ++j) {
            size_t v = nz_odd[j];
            if (u + v > scaled_r2) break; // Exact monotonic bound
            uint64_t v_s0 = f2_odd_s0[v], v_s2 = f2_odd_s2[v];

            if (u_s0 && v_s0) f4_odd_s0[u + v] += u_s0 * v_s0;
            if (u_s2 && v_s2) f4_odd_s2[u + v] += u_s2 * v_s2;
            if (u_s0 && v_s2) f4_odd_s2[u + v] += u_s0 * v_s2;
            if (u_s2 && v_s0) f4_odd_s2[u + v] += u_s2 * v_s0;
        }
    }
}

// Phase 3: Cumulative Prefix Sums & Terminal Single-Pass Dot Product (Strictly O(r^2))
uint64_t execute_E8_volume(uint64_t scaled_r2,
                           const vector<uint64_t>& f4_even_s0, const vector<uint64_t>& f4_even_s2,
                           const vector<uint64_t>& f4_odd_s0,  const vector<uint64_t>& f4_odd_s2) {
    vector<uint64_t> F4_even_s0(scaled_r2 + 1, 0), F4_even_s2(scaled_r2 + 1, 0);
    vector<uint64_t> F4_odd_s0(scaled_r2 + 1, 0),  F4_odd_s2(scaled_r2 + 1, 0);

    uint64_t sum_e0 = 0, sum_e2 = 0, sum_o0 = 0, sum_o2 = 0;
    for (size_t m = 0; m <= scaled_r2; ++m) {
        sum_e0 += f4_even_s0[m]; F4_even_s0[m] = sum_e0;
        sum_e2 += f4_even_s2[m]; F4_even_s2[m] = sum_e2;
        sum_o0 += f4_odd_s0[m];  F4_odd_s0[m]  = sum_o0;
        sum_o2 += f4_odd_s2[m];  F4_odd_s2[m]  = sum_o2;
    }

    // Terminal single-pass 8D dot product
    uint64_t total_E8 = 0;
    for (size_t v = 0; v <= scaled_r2; ++v) {
        size_t rem = scaled_r2 - v;
        if (f4_even_s0[v]) total_E8 += f4_even_s0[v] * F4_even_s0[rem];
        if (f4_even_s2[v]) total_E8 += f4_even_s2[v] * F4_even_s2[rem];
        if (f4_odd_s0[v])  total_E8 += f4_odd_s0[v]  * F4_odd_s0[rem];
        if (f4_odd_s2[v])  total_E8 += f4_odd_s2[v]  * F4_odd_s2[rem];
    }
    return total_E8;
}

int main() {
    // ========================================================================
    // PARAMETER SELECTION ZONE
    // Target squared radius in standard E_8 coordinates (||x||^2 <= target_r2)
    // Preset: target_r2 = 10 matches terminal entry of Table I (count = 56,881)
    // ========================================================================
    uint64_t target_r2 = 10;
    uint64_t scaled_r2 = 4 * target_r2; // 2E_8 scalar dilation

    cout << "====================================================================\n";
    cout << "  E_8 LATTICE HARDWARE-NATIVE BISECTION BENCHMARK ENGINE            \n";
    cout << "====================================================================\n";
    cout << "Target Squared Radius (r^2) : " << target_r2 << "\n";
    cout << "Scaled Domain Bound ((2r)^2): " << scaled_r2 << "\n";
    cout << "--------------------------------------------------------------------\n";

    auto start_time = high_resolution_clock::now();

    // Phase 1: 2D Profiling
    vector<uint64_t> f2_even_s0(scaled_r2 + 1, 0), f2_even_s2(scaled_r2 + 1, 0);
    vector<uint64_t> f2_odd_s0(scaled_r2 + 1, 0),  f2_odd_s2(scaled_r2 + 1, 0);
    build_E8_2D_profiles(scaled_r2, f2_even_s0, f2_even_s2, f2_odd_s0, f2_odd_s2);

    // Phase 2: 4D Sparse Convolution
    vector<uint64_t> f4_even_s0(scaled_r2 + 1, 0), f4_even_s2(scaled_r2 + 1, 0);
    vector<uint64_t> f4_odd_s0(scaled_r2 + 1, 0),  f4_odd_s2(scaled_r2 + 1, 0);
    build_E8_4D_profiles(scaled_r2, f2_even_s0, f2_even_s2, f2_odd_s0, f2_odd_s2,
                         f4_even_s0, f4_even_s2, f4_odd_s0, f4_odd_s2);

    // Phase 3: Terminal Single-Pass 8D Volume
    uint64_t total_points = execute_E8_volume(scaled_r2, f4_even_s0, f4_even_s2,
                                              f4_odd_s0, f4_odd_s2);

    auto stop_time = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop_time - start_time);

    cout << "Total E_8 Lattice Points    : " << total_points << "\n";
    cout << "Execution Time              : " << duration.count() / 1000.0 << " ms\n";
    if (target_r2 == 10) {
        cout << "Verification Status         : " 
             << (total_points == 56881 ? "PASS (Exact match with OEIS A004010)" : "FAIL") << "\n";
    }
    cout << "====================================================================\n";

    return 0;
}