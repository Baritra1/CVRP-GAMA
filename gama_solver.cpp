#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdint>
#include <random>
#include <limits>

using namespace std;


vector<int8_t> vec_add(const vector<int8_t>& a, const vector<int8_t>& b) {
    vector<int8_t> r(a.size());
    for (size_t i = 0; i < a.size(); ++i)
        r[i] = a[i] + b[i];
    return r;
}

bool vec_all_bin(const vector<int8_t>& a) {
    for (double v : a)
        if (v > 1 || v < 0) return false;
    return true;
}

template<typename T>
vector<vector<T>> read_matrix(const string& filename) {
    ifstream fin(filename);
    vector<vector<T>> matrix;
    string line;

    if (!fin.is_open()) {
        cerr << "Error: could not open file " << filename << "\n";
        return matrix;
    }

    while (getline(fin, line)) {
        stringstream ss(line);
        double value;
        vector<T> row;
        while (ss >> value)
            row.push_back(value);
        if (!row.empty())
            matrix.push_back(row);
    }

    return matrix;
}

double objective_function(vector<vector<double>>& cost, vector<int8_t>& candidate, int n){
    double cst=0;
    for (int i=0;i<n*n+n;i++) {
        if(candidate[i]) {
            if(i%(n+1)==n){
                cst+=cost[i/(n+1)][0];
            } else {
                cst+=cost[i/(n+1)][i%(n+1)];
            }
        }
    }
    return cst;
}


int main() {
    vector<vector<int8_t>> feas_sols = read_matrix<int8_t>("data/feas_sols_sorted_2.txt");
    vector<vector<int8_t>> g_basis = read_matrix<int8_t>("data/graver_basis_test_2.txt");
    vector<vector<double>> cost = read_matrix<double>("data/cost_matrix_2.txt");
    cout << "Loaded feas_sols: " << feas_sols.size() << "x" << (feas_sols.empty() ? 0 : feas_sols[0].size()) << "\n";
    cout << "Loaded g_basis: " << g_basis.size() << "x" << (g_basis.empty() ? 0 : g_basis[0].size()) << "\n";
    cout << "Loaded cost: " << cost.size() << "x" << (cost.empty() ? 0 : cost[0].size()) << "\n";
    int n=3;
    vector<int> seeds;
    if (feas_sols.size()<1000) {
        for (int i=0;i<feas_sols.size();i++){
            seeds.push_back(i);
        }
    } else {
        int m=1000;
        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<int> dist(0, feas_sols.size() - 1);
        for (size_t i = 0; i < m; ++i)
            seeds.push_back(static_cast<int>(dist(gen)));
    }
    //ASSUMPTION: there is a path with cost <1.8^10^308
    double global_min = numeric_limits<double>::max();
    vector<int8_t> min_path;
    for(int i=0;i<seeds.size();i++){
        vector<int8_t> sol=feas_sols[seeds[i]];
        double first_cst=objective_function(cost,sol,n);
        if(first_cst<global_min){
            global_min=first_cst;
            min_path=sol;
        }
        bool tst=true;
        while(tst){
            tst=false;
            for (int j=0;j<g_basis.size();j++){
                vector<int8_t> new_sol = vec_add(sol,g_basis[j]);
                if (!vec_all_bin(new_sol)) continue;
                double cst=objective_function(cost,new_sol,n);
                if(cst<global_min) {
                    tst=true;
                    global_min=cst;
                    min_path=new_sol;
                }
            }
        }
    }
    cout << global_min << "\n";
    for (int i=0;i<n;i++) {
        for (int j=0;j<n+1;j++) {
            cout << (int)min_path[i*(n+1)+j] << " ";
        }
        cout<<"\n";
    }
}