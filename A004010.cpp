/**
 * OEIS A004010: Number of lattice points in 8-dimensional E_8 lattice with norm squared <= 2n.
 *
 * Logic derived from:
 * Paper V: Morris (2026), "Parity-Filtered Bisection of the E_8 Lattice: Hardware-Native
 *          Dimension-Paired Enumeration via Diophantine Scalar Transformation", Appendix.
 *
 * Architecture:
 * - 2E_8 scalar mapping: y = 2x eliminates half-integers, mapping norm <= 2n to scaled_norm <= 8n.
 * - Residue collapse: modulo-4 sum constraint bifurcates into All-Even and All-Odd parity sublattices.
 * - Hardware-Native: Executes strictly on the integer ALU with zero floating-point emulation.
 * - Double-Check Verification: Asserts each calculated term against the divisor-sum formula
 *   r_{E8}(2k) = 240 * \sigma_3(k) in real time.
 *
 * Visual Studio Usage: Press Ctrl+F5 to run. Automatically creates "b004010.txt" in project folder.
 */

#include <iostream>
#include <vector>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <string>
#include <chrono>
#include <cmath>

using namespace std;
using namespace std::chrono;

// Bitwise integer square root to execute purely on the ALU
static inline uint64_t pure_isqrt(uint64_t value) {
    uint64_t res = 0;
    uint64_t bit = 1ULL << 62;
    while (bit > value) bit >>= 2;
    while (bit != 0) {
        if (value >= res + bit) {
            value -= res + bit;
            res = (res >> 1) + bit;
        }
        else {
            res >>= 1;
        }
        bit >>= 2;
    }
    return res;
}

int main(int argc, char* argv[]) {
    size_t max_n = 10000;
    string out_filename = "b004010.txt";
    if (argc > 1) {
        max_n = static_cast<size_t>(atoll(argv[1]));
    }

    // 2E_8 scalar mapping: ||y||^2 = 4 * ||x||^2 <= 4 * (2 * max_n) = 8 * max_n
    const size_t scaled_max_m = 8 * max_n;

    cout << "====================================================================\n";
    cout << "  OEIS A004010 EXTENDED B-FILE GENERATOR (E_8 LATTICE BISECTION)   \n";
    cout << "====================================================================\n";
    cout << "Target Terms           : n = 0 to " << max_n << " (" << (max_n + 1) << " total terms)\n";
    cout << "Scaled Domain Bound    : ||y||^2 <= " << scaled_max_m << " (via 2E_8 scalar dilation)\n";
    cout << "Output File            : " << out_filename << "\n";
    cout << "--------------------------------------------------------------------\n";

    auto start_total = high_resolution_clock::now();

    // ========================================================================
    // PHASE 1: 2D BASE PROFILE GENERATION (O(r^2))
    // ========================================================================
    cout << "[1/4] Constructing 2D parity-filtered base profiles...\n";
    auto t0 = high_resolution_clock::now();

    vector<uint64_t> f2_even_s0(scaled_max_m + 1, 0), f2_even_s2(scaled_max_m + 1, 0);
    vector<uint64_t> f2_odd_s0(scaled_max_m + 1, 0), f2_odd_s2(scaled_max_m + 1, 0);

    int64_t max_coord = static_cast<int64_t>(pure_isqrt(scaled_max_m));

    // All-Even 2D Coordinate Sweep (Steps by 2)
    int64_t min_even = (max_coord % 2 != 0) ? -(max_coord - 1) : -max_coord;
    for (int64_t x = min_even; x <= max_coord; x += 2) {
        uint64_t x2 = static_cast<uint64_t>(x * x);
        for (int64_t y = min_even; y <= max_coord; y += 2) {
            uint64_t m = x2 + static_cast<uint64_t>(y * y);
            if (m <= scaled_max_m) {
                if (abs(x + y) % 4 == 0) f2_even_s0[m]++;
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
            if (m <= scaled_max_m) {
                if (abs(x + y) % 4 == 0) f2_odd_s0[m]++;
                else f2_odd_s2[m]++;
            }
        }
    }

    auto t1 = high_resolution_clock::now();
    cout << "      Completed in " << duration_cast<microseconds>(t1 - t0).count() / 1000.0 << " ms.\n";

    // ========================================================================
    // PHASE 2: 4D INTERMEDIATE CONVOLUTION (Sparse Index-Accelerated)
    // ========================================================================
    cout << "[2/4] Convolving 2D profiles into 4D intermediate submanifolds...\n";
    t0 = high_resolution_clock::now();

    // Collect non-zero entries for sparse iteration
    vector<size_t> nz_even, nz_odd;
    nz_even.reserve(scaled_max_m / 10);
    nz_odd.reserve(scaled_max_m / 10);
    for (size_t m = 0; m <= scaled_max_m; ++m) {
        if (f2_even_s0[m] != 0 || f2_even_s2[m] != 0) nz_even.push_back(m);
        if (f2_odd_s0[m] != 0 || f2_odd_s2[m] != 0)   nz_odd.push_back(m);
    }

    vector<uint64_t> f4_even_s0(scaled_max_m + 1, 0), f4_even_s2(scaled_max_m + 1, 0);
    vector<uint64_t> f4_odd_s0(scaled_max_m + 1, 0), f4_odd_s2(scaled_max_m + 1, 0);

    // Even + Even crosses
    for (size_t i = 0; i < nz_even.size(); ++i) {
        size_t u = nz_even[i];
        uint64_t u_s0 = f2_even_s0[u];
        uint64_t u_s2 = f2_even_s2[u];
        for (size_t j = 0; j < nz_even.size(); ++j) {
            size_t v = nz_even[j];
            if (u + v > scaled_max_m) break;
            uint64_t v_s0 = f2_even_s0[v];
            uint64_t v_s2 = f2_even_s2[v];

            if (u_s0 != 0 && v_s0 != 0) f4_even_s0[u + v] += u_s0 * v_s0;
            if (u_s2 != 0 && v_s2 != 0) f4_even_s0[u + v] += u_s2 * v_s2;
            if (u_s0 != 0 && v_s2 != 0) f4_even_s2[u + v] += u_s0 * v_s2;
            if (u_s2 != 0 && v_s0 != 0) f4_even_s2[u + v] += u_s2 * v_s0;
        }
    }

    // Odd + Odd crosses
    for (size_t i = 0; i < nz_odd.size(); ++i) {
        size_t u = nz_odd[i];
        uint64_t u_s0 = f2_odd_s0[u];
        uint64_t u_s2 = f2_odd_s2[u];
        for (size_t j = 0; j < nz_odd.size(); ++j) {
            size_t v = nz_odd[j];
            if (u + v > scaled_max_m) break;
            uint64_t v_s0 = f2_odd_s0[v];
            uint64_t v_s2 = f2_odd_s2[v];

            if (u_s0 != 0 && v_s0 != 0) f4_odd_s0[u + v] += u_s0 * v_s0;
            if (u_s2 != 0 && v_s2 != 0) f4_odd_s0[u + v] += u_s2 * v_s2;
            if (u_s0 != 0 && v_s2 != 0) f4_odd_s2[u + v] += u_s0 * v_s2;
            if (u_s2 != 0 && v_s0 != 0) f4_odd_s2[u + v] += u_s2 * v_s0;
        }
    }

    t1 = high_resolution_clock::now();
    cout << "      Completed in " << duration_cast<milliseconds>(t1 - t0).count() << " ms.\n";

    // ========================================================================
    // PHASE 3: PARALLEL VERIFICATION SIEVE (Divisor-Sum \sigma_3(k))
    // ========================================================================
    cout << "[3/4] Initializing ground-truth verification sieve (\\sigma_3)...\n";
    vector<uint64_t> sigma3(max_n + 1, 0);
    for (size_t d = 1; d <= max_n; ++d) {
        uint64_t d3 = static_cast<uint64_t>(d) * d * d;
        for (size_t k = d; k <= max_n; k += d) {
            sigma3[k] += d3;
        }
    }

    // ========================================================================
    // PHASE 4: 8D SPATIAL ENUMERATION & FILE GENERATION
    // ========================================================================
    cout << "[4/4] Resolving 8D lattice bounds and writing " << out_filename << "...\n";
    t0 = high_resolution_clock::now();

    ofstream outfile(out_filename);
    if (!outfile.is_open()) {
        cerr << "Error: Could not open file " << out_filename << " for writing.\n";
        return 1;
    }

    uint64_t total_points = 0;
    bool all_verified = true;

    for (size_t n = 0; n <= max_n; ++n) {
        size_t T = 8 * n; // Scaled shell target
        uint64_t shell_points = 0;
        size_t half = T / 2;

        // Symmetric convolution of 4D submanifolds
        for (size_t u = 0; u < half; ++u) {
            size_t v = T - u;
            uint64_t e_cross = f4_even_s0[u] * f4_even_s0[v] + f4_even_s2[u] * f4_even_s2[v];
            uint64_t o_cross = f4_odd_s0[u] * f4_odd_s0[v] + f4_odd_s2[u] * f4_odd_s2[v];
            shell_points += 2 * (e_cross + o_cross);
        }

        // Center reflection point (u == v == half)
        uint64_t e_mid = f4_even_s0[half] * f4_even_s0[half] + f4_even_s2[half] * f4_even_s2[half];
        uint64_t o_mid = f4_odd_s0[half] * f4_odd_s0[half] + f4_odd_s2[half] * f4_odd_s2[half];
        shell_points += e_mid + o_mid;

        // Accumulate cumulative total
        total_points += shell_points;

        // On-the-fly mathematical verification against theta series
        uint64_t expected_shell = (n == 0) ? 1 : (240 * sigma3[n]);
        if (shell_points != expected_shell) {
            all_verified = false;
            cerr << "DISCREPANCY at n = " << n << "! Spatial: " << shell_points
                << " vs Theta Series: " << expected_shell << "\n";
        }

        outfile << n << " " << total_points << "\n";
    }

    outfile.close();
    t1 = high_resolution_clock::now();
    cout << "      Completed in " << duration_cast<milliseconds>(t1 - t0).count() << " ms.\n";

    auto stop_total = high_resolution_clock::now();
    auto total_ms = duration_cast<milliseconds>(stop_total - start_total).count();

    cout << "--------------------------------------------------------------------\n";
    cout << "Final Total at n = " << max_n << " : " << total_points << "\n";
    cout << "Verification Status      : "
        << (all_verified ? "PASS (100% Bit-Exact across all terms)" : "FAIL (Discrepancy detected)") << "\n";
    cout << "Total Elapsed Time       : " << total_ms << " ms ("
        << total_ms / 1000.0 << " seconds)\n";
    cout << "Generated File           : " << out_filename << " (" << (max_n + 1) << " entries)\n";
    cout << "====================================================================\n";

    return 0;
}