#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#ifdef _WIN32
#include <windows.h>
#endif
using namespace std;

double f(double x) {
    return 1.75 * (x + 2.3) * (x - 1.05) * (x - 1.7) * (x - 3.0) + 1.0;
}

double df(double x) {
    return 7.0 * pow(x, 3) - 18.1125 * pow(x, 2) - 11.165 * x + 31.019625;
}

double d2f(double x) {
    return 21.0 * pow(x, 2) - 36.225 * x - 11.165;
}

// чисельна перевірка похідних 
void check_derivatives(double x_test, double h = 1e-5) {
    double df_exact = df(x_test);
    double df_num = (f(x_test + h) - f(x_test - h)) / (2.0 * h);

    double d2f_exact = d2f(x_test);
    double d2f_num = (df(x_test + h) - df(x_test - h)) / (2.0 * h);

    cout << "=== перевірка похідних у точці x0 = " << fixed << setprecision(4) << x_test << " ===" << endl;
    cout << left << setw(10) << "похідна" 
         << setw(16) << "аналітично" 
         << setw(16) << "чисельно" 
         << "похибка" << endl;
    cout << string(56, '-') << endl;
    
    cout << left << setw(10) << "f'(x0)" 
         << setw(16) << fixed << setprecision(8) << df_exact 
         << setw(16) << df_num 
         << scientific << setprecision(2) << abs(df_exact - df_num) << endl;
         
    cout << left << setw(10) << "f''(x0)" 
         << setw(16) << fixed << setprecision(8) << d2f_exact 
         << setw(16) << d2f_num 
         << scientific << setprecision(2) << abs(d2f_exact - d2f_num) << endl;
    cout << endl;
}

// 2. частина А: класичний метод Ньютона 
void newton_classic(double x0, double eps = 0.001, int max_iter = 50, bool verbose = true) {
    double x = x0;
    int n_df = 0, n_d2f = 0;
    int k;

    if (verbose) {
        cout << "=== класичний метод ньютона (x0 = " << fixed << setprecision(4) << x0 << ") ===" << endl;
        cout << left << setw(4) << "k" 
             << setw(12) << "xk" 
             << setw(14) << "f(xk)" 
             << setw(14) << "f'(xk)" 
             << setw(14) << "|f'(xk)|" 
             << "f''(xk)" << endl;
        cout << string(68, '-') << endl;
    }

    for (k = 0; k < max_iter; k++) {
        double fx = f(x);
        double d1 = df(x);
        n_df++;

        if (abs(d1) <= eps) {
            if (verbose) {
                cout << left << setw(4) << k 
                     << setw(12) << fixed << setprecision(6) << x 
                     << setw(14) << fx 
                     << setw(14) << d1 
                     << setw(14) << abs(d1) 
                     << "-" << endl;
            }
            break;
        }

        double d2 = d2f(x);
        n_d2f++;

        if (d2 <= 1e-12) {
            if (verbose) {
                cout << string(68, '-') << endl;
                cout << "зупинка: у точці x = " << x << " f''(x) = " << d2 << " <= 0 (напрямок до максимуму!)" << endl << endl;
            }
            return;
        }

        if (verbose) {
            cout << left << setw(4) << k 
                 << setw(12) << fixed << setprecision(6) << x 
                 << setw(14) << fx 
                 << setw(14) << d1 
                 << setw(14) << abs(d1) 
                 << d2 << endl;
        }

        x = x - d1 / d2;
    }

    if (verbose) {
        cout << string(68, '-') << endl;
        cout << "оптимальна точка x* = " << fixed << setprecision(6) << x << endl;
        cout << "значення функції f(x*) = " << f(x) << endl;
        cout << "кількість ітерацій n = " << k << endl;
        cout << "кількість обчислень: " << n_df << " f' + " << n_d2f << " f''\n" << endl;
    }
}

// 3. частина Б: метод Ньютона-Рафсона з регулюванням кроку 
void newton_raphson_regulated(double x0, double eps = 0.001, double alpha = 0.5, int max_iter = 50, bool verbose = true) {
    double x = x0;
    double fx = f(x), d1 = df(x), d2 = d2f(x);
    int n_f = 1, n_df = 1, n_d2f = 1;
    int k;

    if (verbose) {
        cout << "=== НЬЮТОН-РАФСОН ІЗ РЕГ. КРОКУ (x0 = " << fixed << setprecision(4) << x0 << ") ===" << endl;
        cout << left << setw(4) << "k" 
             << setw(11) << "xk" 
             << setw(12) << "f(xk)" 
             << setw(11) << "f'(xk)" 
             << setw(11) << "f''(xk)" 
             << setw(7)  << "t" 
             << setw(11) << "x_tilde" 
             << setw(12) << "f(x_til)" 
             << setw(12) << "f'(x_til)" 
             << "(3.10)" << endl;
        cout << string(104, '-') << endl;
        
        cout << left << setw(4) << 0 
             << setw(11) << fixed << setprecision(6) << x 
             << setw(12) << fx 
             << setw(11) << d1 
             << setw(11) << d2 
             << setw(7)  << "-" 
             << setw(11) << "-" 
             << setw(12) << "-" 
             << setw(12) << "-" 
             << "Старт" << endl;
    }

    for (k = 1; k <= max_iter; k++) {
        if (abs(d1) <= eps) break;

        if (d2 <= 1e-12) {
            if (verbose) {
                cout << string(104, '-') << endl;
                cout << "зупинка: f''(x) = " << d2 << " <= 0 — напрямок Ньютона не вниз!" << endl << endl;
            }
            return;
        }

        double t = 1.0;
        while (true) {
            double xb = x - t * d1 / d2;
            double d1b = df(xb);
            n_df++;

            if (abs(d1b) <= eps) {
                double fb = f(xb);
                if (verbose) {
                    cout << left << setw(4) << k 
                         << setw(11) << fixed << setprecision(6) << x 
                         << setw(12) << fx 
                         << setw(11) << d1 
                         << setw(11) << d2 
                         << setw(7)  << setprecision(3) << t 
                         << setw(11) << setprecision(6) << xb 
                         << setw(12) << fb 
                         << setw(12) << d1b 
                         << "Стоп" << endl;
                    cout << string(104, '-') << endl;
                    cout << "оптимальна точка x* = " << xb << endl;
                    cout << "значення функції f(x*) = " << fb << endl;
                    cout << "кількість ітерацій n = " << k << endl;
                    cout << "кількість обчислень: " << n_f << " f + " << n_df << " f' + " << n_d2f << " f''\n" << endl;
                }
                return;
            }

            double fb = f(xb);
            n_f++;
            bool cond = (fb - fx) <= -alpha * t * (d1 * d1) / d2;

            if (verbose) {
                string status = cond ? "Так(v)" : "Ні(x)";
                cout << left << setw(4) << k 
                     << setw(11) << fixed << setprecision(6) << x 
                     << setw(12) << fx 
                     << setw(11) << d1 
                     << setw(11) << d2 
                     << setw(7)  << setprecision(3) << t 
                     << setw(11) << setprecision(6) << xb 
                     << setw(12) << fb 
                     << setw(12) << d1b 
                     << status << endl;
            }

            if (cond) {
                x = xb;
                fx = fb;
                d1 = d1b;
                d2 = d2f(x);
                n_d2f++;
                break;
            }

            t /= 2.0;
            if (t < 1e-8) {
                cout << "помилка: дроблення кроку не допомогло (t < 1e-8)." << endl;
                return;
            }
        }
    }
}

int main() {
     #ifdef _WIN32
    SetConsoleCP(CP_UTF8);
    SetConsoleOutputCP(CP_UTF8);
#endif
    double a0 = -1.35;
    double b0 = -1.05;
    double eps = 0.001;
    double alpha = 0.5;

    // пункт 1: перевірка похідних
    check_derivatives(a0);

    //пункт 3: частина А — класичний метод Ньютона
    newton_classic(a0, eps);

    //пункт 4: частина Б — метод Ньютона-Рафсона з регулюванням кроку
    newton_raphson_regulated(a0, eps, alpha);

    //пункт 7: експеримент зі стартами
    cout << "=== експеримент зі стартами (пункт 7) ===" << endl;
    
    cout << "1) старт з правої межі інтервалу Свенна x0 = b0 = -1.05:" << endl;
    newton_classic(b0, eps);

    cout << "2) старт біля точки перегину x0 = -0.35 (класичний Ньютон):" << endl;
    newton_classic(-0.35, eps);

    cout << "3) старт з x0 = -0.35 (Ньютон-Рафсон із регулюванням кроку):" << endl;
    newton_raphson_regulated(-0.35, eps, alpha);

    cout << "4) старт у точці з від'ємною кривиною x0 = 0.0:" << endl;
    newton_classic(0.0, eps);

    return 0;
}