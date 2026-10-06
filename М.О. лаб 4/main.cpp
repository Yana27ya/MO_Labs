#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>
#ifdef _WIN32
#include <windows.h>
#endif

using namespace std;

// варіант 7
double f(double x) {
    return 1.75 * (x + 2.3) * (x - 1.05) * (x - 1.7) * (x - 3.0) + 1.0;
}

int main() {
     #ifdef _WIN32
    SetConsoleCP(CP_UTF8);
    SetConsoleOutputCP(CP_UTF8);
#endif
    double a = -1.35, b = -1.05, eps = 0.05;
    double L0 = b - a;
    double R = L0 / eps;

    // 1. таблиця чисел фібоначчі
    vector<double> F = {0, 1, 1}; 
    while (F.back() < R) {
        F.push_back(F[F.size() - 1] + F[F.size() - 2]);
    }
    int N = F.size() - 1;

    double x1, x2, f1, f2;
    int nf = 0; 
    bool x1_inherited = false, x2_inherited = false;
    double delta = eps / 10.0; 

    cout << fixed << setprecision(5);
    cout << "R = " << R << ", N = " << N << ", F_N = " << F[N] << "\n\n";
    cout << "k\tm\tF_m-2/F_m\tF_m-1/F_m\ta\t\tb\t\tx1\t\tx2\t\tf(x1)\t\tf(x2)\t\tL_k\tРішення\n";

    // 2.основний цикл методу фібоначчі
    for (int k = 0; k <= N - 3; ++k) {
        int m = N - k;
        double frac1 = F[m - 2] / F[m];
        double frac2 = F[m - 1] / F[m];

        if (m > 3) {
            if (!x1_inherited) { x1 = a + frac1 * (b - a); f1 = f(x1); nf++; }
            if (!x2_inherited) { x2 = a + frac2 * (b - a); f2 = f(x2); nf++; }
        } else { 
    
            if (!x1_inherited) { x1 = a + (b - a) / 2.0; f1 = f(x1); nf++; }
            x2 = x1 + delta; 
            f2 = f(x2); nf++;
            frac1 = 0.5; frac2 = 0.5;
        }

        cout << k << "\t" << m << "\t" << F[m-2] << "/" << F[m] << "\t\t" << F[m-1] << "/" << F[m] << "\t\t"
             << a << "\t" << b << "\t" << x1 << "\t" << x2 << "\t" << f1 << "\t" << f2 << "\t" << (b - a) << "\t";

        //правило виключення
        if (f1 <= f2) {
            cout << "f1<=f2 -> b=x2\n";
            b = x2; x2 = x1; f2 = f1;
            x2_inherited = true; x1_inherited = false;
        } else {
            cout << "f1>f2 -> a=x1\n";
            a = x1; x1 = x2; f1 = f2;
            x1_inherited = true; x2_inherited = false;
        }
    }

    double L_kin = b - a;
    double x_star = (a + b) / 2.0;

    cout << "------------------------------------------------------------------------------------------------------------------------\n";
    cout << "оптимальна точка x* = " << x_star << "\n";
    cout << "кількість обчислень функції N_f = " << nf << "\n";
    cout << "кінцева довжина L_kin = " << L_kin << "\n";
    cout << "теоретична межа L0/F_N = " << L0 / F[N] << "\n";

    return 0;
}