#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#ifdef _WIN32
#include <windows.h>
#endif

using namespace std;

double F(double x1, double x2) {
    return 3.0 * x1 * x1 + x2 * x2 + 2.0 * x1 * x2 + 16.0 * x1 + 4.0 * x2;
}

void grad(double x1, double x2, double &g1, double &g2) {
    g1 = 6.0 * x1 + 2.0 * x2 + 16.0;
    g2 = 2.0 * x2 + 2.0 * x1 + 4.0;
}

double sign_val(double v) {
    if (v > 0) return 1.0;
    if (v < 0) return -1.0;
    return 0.0;
}

void direction(double g1, double g2, const string &metric, double &d1, double &d2) {
    if (metric == "sphere") {
        double norm = sqrt(g1 * g1 + g2 * g2);
        d1 = -g1 / norm;
        d2 = -g2 / norm;
    } else if (metric == "cube") {
        d1 = -sign_val(g1);
        d2 = -sign_val(g2);
    } else { // octa
        if (abs(g1) >= abs(g2)) {
            d1 = -sign_val(g1);
            d2 = 0.0;
        } else {
            d1 = 0.0;
            d2 = -sign_val(g2);
        }
    }
}

void check_derivatives(double x1, double x2, double h = 1e-5) {
    double g1, g2;
    grad(x1, x2, g1, g2);
    double num_g1 = (F(x1 + h, x2) - F(x1 - h, x2)) / (2.0 * h);
    double num_g2 = (F(x1, x2 + h) - F(x1, x2 - h)) / (2.0 * h);

    cout << "=== ПЕРЕВІРКА ПОХІДНИХ У x0 = (" << x1 << "; " << x2 << ") ===" << endl;
    cout << fixed << setprecision(8);
    cout << "dF/dx1: аналіт. = " << g1 << ", чисельн. = " << num_g1
         << ", похибка = " << scientific << setprecision(2) << abs(g1 - num_g1) << endl;
    cout << fixed << setprecision(8);
    cout << "dF/dx2: аналіт. = " << g2 << ", чисельн. = " << num_g2
         << ", похибка = " << scientific << setprecision(2) << abs(g2 - num_g2) << endl << endl;
}

struct Result {
    string metric;
    double h0;
    int n;
    int nF;
    int nG;
    double x1, x2, Fx, err;
};

Result run_gradient(double x1_0, double x2_0, double h0, double eps, const string &metric, bool verbose = true) {
    double x1 = x1_0, x2 = x2_0;
    double Fx = F(x1, x2);
    double h = h0;
    int nF = 1, nG = 0, k = 0;

    if (verbose) {
        cout << "=== МЕТРИКА: " << metric << " (h0 = " << fixed << setprecision(1) << h0 << ") ===" << endl;
        cout << left << setw(4) << "k" << setw(11) << "x1" << setw(11) << "x2"
             << setw(12) << "F(x^k)" << setw(10) << "g1" << setw(10) << "g2"
             << setw(10) << "||g||2" << setw(9) << "d1" << setw(9) << "d2"
             << setw(9) << "h" << setw(7) << "дробл" << "F(x^(k+1))" << endl;
        cout << string(112, '-') << endl;
    }

    for (k = 1; k <= 5000; k++) {
        double g1, g2;
        grad(x1, x2, g1, g2);
        nG++;
        double gnorm = sqrt(g1 * g1 + g2 * g2);
        if (gnorm <= eps) {
            k--; // остання ітерація не відбулася
            break;
        }

        double d1, d2;
        direction(g1, g2, metric, d1, d2);

        int halvings = 0;
        double xb1, xb2, Fb;
        while (true) {
            xb1 = x1 + h * d1;
            xb2 = x2 + h * d2;
            Fb = F(xb1, xb2);
            nF++;
            if (Fb < Fx) break;
            h /= 2.0;
            halvings++;
            if (h < 1e-10) break;
        }

        if (verbose) {
            cout << left << setw(4) << k << fixed << setprecision(5)
                 << setw(11) << x1 << setw(11) << x2 << setw(12) << Fx
                 << setprecision(3) << setw(10) << g1 << setw(10) << g2 << setw(10) << gnorm
                 << setw(9) << d1 << setw(9) << d2 << setprecision(4) << setw(9) << h
                 << setw(7) << halvings << setprecision(5) << Fb << endl;
        }

        x1 = xb1; x2 = xb2; Fx = Fb;
    }

    double err = sqrt(pow(x1 - (-3.0), 2) + pow(x2 - 1.0, 2));
    if (verbose) {
        cout << string(112, '-') << endl;
        cout << fixed << setprecision(6);
        cout << "x* = (" << x1 << "; " << x2 << "), F(x*) = " << Fx
             << ", n = " << k << ", NF = " << nF << ", N_grad = " << nG
             << ", Разом = " << (nF + nG) << ", ||x* - x_точн|| = "
             << scientific << setprecision(2) << err << endl << endl;
    }

    return {metric, h0, k, nF, nG, x1, x2, Fx, err};
}

int main() {
#ifdef _WIN32
    SetConsoleCP(CP_UTF8);
    SetConsoleOutputCP(CP_UTF8);
#endif

    double x1_0 = 3.0, x2_0 = -3.0, eps = 0.01;
    check_derivatives(x1_0, x2_0);


    Result r1 = run_gradient(x1_0, x2_0, 1.0, eps, "sphere", true);
    Result r2 = run_gradient(x1_0, x2_0, 1.0, eps, "cube", true);
    Result r3 = run_gradient(x1_0, x2_0, 1.0, eps, "octa", true);

    Result r4 = run_gradient(x1_0, x2_0, 0.7, eps, "sphere", false);
    Result r5 = run_gradient(x1_0, x2_0, 0.7, eps, "cube", false);
    Result r6 = run_gradient(x1_0, x2_0, 0.7, eps, "octa", false);

    // зведена таблиця
    cout << "=== ЗВЕДЕНА ТАБЛИЦЯ ===" << endl;
    cout << left << setw(10) << "Метрика" << setw(6) << "h0" << setw(5) << "n"
         << setw(6) << "NF" << setw(7) << "N_g" << setw(8) << "Разом"
         << setw(26) << "x*" << "||x* - x_точн||" << endl;
    cout << string(82, '-') << endl;

    Result all[6] = {r1, r2, r3, r4, r5, r6};
    for (int i = 0; i < 6; i++) {
        cout << left << setw(10) << all[i].metric
             << fixed << setprecision(1) << setw(6) << all[i].h0
             << setw(5) << all[i].n << setw(6) << all[i].nF
             << setw(7) << all[i].nG << setw(8) << (all[i].nF + all[i].nG)
             << "(" << setprecision(4) << all[i].x1 << "; " << all[i].x2 << ")   "
             << scientific << setprecision(2) << all[i].err << endl;
    }

    return 0;
}