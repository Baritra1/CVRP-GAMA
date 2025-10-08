#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <sstream>
#include <mutex>
#include <unordered_set>
#include <omp.h>
#include <chrono>
#include <algorithm>
#include <execution>



using namespace std;

/**
 * I ported compute_graver_basis.py as found in https://doi.org/10.1287/ijoc.2024.0574.cd from
 * Python to C++, which already makes it significantly faster. Then, I optimized some slow numpy
 * functions to be significantly faster (some stuff I optimized included reducing the number of O(n^2)
 * copies and using some compiler options like -O3 and -ffast-math), and used the fact that the solutions are all either 0 1
 * or -1 to make the functions work with int8_t, which is significantly faster than numpy working
 * with doubles. Finally, I made another very significant optimization by parallelizing the entirety
 * of compute_kernel_feas_sols() and a large part of compute_gbasis_par() to significantly increase
 * speed. For a feas_sols_sorted.txt that started off as having 95 solutions, compute_graver_basis.py
 * was able to get through 200 iterations of compute_gbasis_par() in 3 hours, whereas compute_graver_basis.cpp
 * gets through ~24k iterations (and solves the problem) in ~6 minutes (I didn't time this I'm estimating)
 */


// Vector functions

double sum_abs_diff(const vector<int8_t>& a, const vector<int8_t>& b) {
    double s = 0;
    for (size_t i = 0; i < a.size(); ++i)
        s += abs(a[i] - b[i]);
    return s;
}

vector<int8_t> vec_sub(const vector<int8_t>& a, const vector<int8_t>& b) {
    vector<int8_t> r(a.size());
    for (size_t i = 0; i < a.size(); ++i)
        r[i] = a[i] - b[i];
    return r;
}

vector<int8_t> vec_abs(const vector<int8_t>& a) {
    vector<int8_t> r(a.size());
    for (size_t i = 0; i < a.size(); ++i)
        r[i] = abs(a[i]);
    return r;
}

bool vec_all_ge0(const vector<int8_t>& a) {
    for (double v : a)
        if (v < 0) return false;
    return true;
}

bool vec_leq_abs(const vector<int8_t>& a, const vector<int8_t>& b) {
    for (size_t i = 0; i < a.size(); ++i)
        if (abs(a[i]) > abs(b[i]))
            return false;
    return true;
}

bool vec_equal(const vector<int8_t>& a, const vector<int8_t>& b, double eps = 1e-9) {
    for (size_t i = 0; i < a.size(); ++i)
        if (abs(a[i] - b[i]) > eps)
            return false;
    return true;
}

bool vec_zero(const vector<int8_t>& a, double eps = 1e-9) {
    for (double v : a)
        if (abs(v) > eps)
            return false;
    return true;
}

std::vector<std::vector<int8_t>> dedup(const std::vector<std::vector<int8_t>>& data) {
    auto data_copy = data;
    auto start_time = std::chrono::steady_clock::now();
    sort(execution::par, data_copy.begin(), data_copy.end());
    auto last = unique(execution::par, data_copy.begin(), data_copy.end());
    data_copy.erase(last, data_copy.end());
    return data_copy;
}

vector<int8_t> normal_form(vector<int8_t> r_i, const vector<vector<int8_t>>& g_basis, size_t skip_index) {
    vector<int8_t> res = r_i;
    for (int i=0;i<g_basis.size();i++) {
        if(i==skip_index) continue;
        auto& g_i_chk =g_basis[i];
        bool f_chk1 = true;
        bool f_chk2 = true;

        // Only iterate over non-zero entries of g_i_chk
        for (size_t k = 0; k < g_i_chk.size(); ++k) {
            int8_t gi = g_i_chk[k];
            if (gi == 0) continue;

            int8_t ri = r_i[k];

            if (gi * ri < 0) f_chk1 = false;
            if (abs(gi) > abs(ri)) f_chk2 = false; 

            if (!f_chk1 && !f_chk2) break; 
        }

        if (f_chk1 && f_chk2) {
            for (size_t k = 0; k < r_i.size(); ++k) {
                r_i[k] -= g_i_chk[k];
            }
        }
    }
    return r_i;
}

vector<vector<int8_t>> compute_kernel_feas_sols(const vector<vector<int8_t>>& feas_sols) {
    size_t n = feas_sols.size();
    vector<vector<vector<int8_t>>> thread_local_results;

    int num_threads = 1;
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        #pragma omp single
        num_threads = omp_get_num_threads();
        vector<vector<int8_t>> local_results;

        #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i + 1; j < n; ++j) {
                vector<int8_t> diff1 = vec_sub(feas_sols[i], feas_sols[j]);
                vector<int8_t> diff2 = vec_sub(feas_sols[j], feas_sols[i]);
                local_results.push_back(diff1);
                local_results.push_back(diff2);

                if (sum_abs_diff(feas_sols[i], feas_sols[j]) == 1.0)
                    #pragma omp critical
                    cout << i << " " << j << endl;
            }
        }

        #pragma omp critical
        thread_local_results.push_back(std::move(local_results));
    }

    // merge results
    vector<vector<int8_t>> ker_sols;
    for (auto& tvec : thread_local_results) {
        ker_sols.insert(ker_sols.end(),
                        make_move_iterator(tvec.begin()),
                        make_move_iterator(tvec.end()));
    }

    cout << "done1" << endl;
    return ker_sols;
}

vector<vector<int8_t>> compute_gbasis_par(vector<vector<int8_t>> ker_par) {
    vector<vector<int8_t>> ker_el_ones_copy = ker_par;
    bool tst1 = true;
    int cc = 0;

    while (tst1) {
        tst1 = false;
        size_t n_rows = ker_el_ones_copy.size();

        // store updated rows
        vector<vector<int8_t>> updated_rows(n_rows);

        // mark what's being deleted
        vector<bool> del_mask(n_rows, false);

        // Parallel computation of normal_form for each row
        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < n_rows; ++i) {
            // cout<<i<<"\n";
            updated_rows[i] = normal_form(ker_el_ones_copy[i],ker_el_ones_copy, i);
        }

        // Sequential update changed rows
        for (size_t i = 0; i < n_rows; ++i) {
            const auto& new_row = updated_rows[i];
            auto& old_row = ker_el_ones_copy[i];

            if (!vec_equal(new_row, old_row)) {
                if (!vec_zero(new_row)) {
                    old_row = new_row;
                    tst1 = true;
                    cout << "i=" << i << " cc=" << cc << endl;
                } else {
                    del_mask[i] = true;
                }
            }
        }

        // Delete marked rows with boolean mask
        vector<vector<int8_t>> new_mat;
        new_mat.reserve(n_rows - count(del_mask.begin(), del_mask.end(), true));
        for (size_t i = 0; i < n_rows; ++i)
            if (!del_mask[i])
                new_mat.push_back(std::move(ker_el_ones_copy[i]));

        ker_el_ones_copy.swap(new_mat);
        cc++;
        cout << "done2" << endl;
    }
    cout<<ker_el_ones_copy.size()<<"\n";
    return ker_el_ones_copy;
}
vector<vector<int8_t>> ker_sols;
vector<vector<int8_t>> ker_sols_uniq;
auto get_graver_basis(const vector<vector<int8_t>>& fsol_list_reduced, long long max_gbasis_size) {

    ker_sols = compute_kernel_feas_sols(fsol_list_reduced);

    // Deduplicate
    ker_sols_uniq = dedup(ker_sols);
    vector<vector<int8_t>> slice_ker_sols_uniq;


    // Shuffle
    random_device rd;
    mt19937 g(rd());
    shuffle(ker_sols_uniq.begin(), ker_sols_uniq.end(), g);

    for (size_t i = 0; i < min((long long)ker_sols_uniq.size(),max_gbasis_size); ++i) {
        slice_ker_sols_uniq.push_back(ker_sols_uniq[i]);
    }
    
    vector<vector<int8_t>> g_basis_par_i = compute_gbasis_par(slice_ker_sols_uniq);
    return g_basis_par_i;
}

int main(int argc, char* argv[]) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    ifstream fin("data/feas_sols_sorted_2.txt");
    vector<vector<int8_t>> feas_solns;
    string line;

    while (getline(fin, line)) {
        stringstream ss(line);
        double val;
        vector<int8_t> row;
        while (ss >> val)
            row.push_back(val);
        if (!row.empty())
            feas_solns.push_back(row);
    }
    int n=stoi(argv[1]);
    long long max_gbasis_size=stoll(argv[2]);
    auto gbasis = get_graver_basis(feas_solns,max_gbasis_size);
    cout << "Final g_basis size: " << gbasis.size() << " x " 
         << (gbasis.empty() ? 0 : gbasis[0].size()) << endl;

    ofstream fout("data/graver_basis_test_2.txt");
    for (const auto& row : gbasis) {
        fout.write(reinterpret_cast<const char*>(row.data()), row.size());
    }
    cout<<"you can end this program now (if there are a lot of feas_sols the destructor for feas_sols takes very long)"<<"\n";
    cout<<gbasis.size()<<" "<< (gbasis.empty() ? 0 : gbasis[0].size()) << endl;
    fout.close();

    return 0;
}
