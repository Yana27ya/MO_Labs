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

// похідна
double df(double x) {
    return 7.0 * pow(x, 3) - 18.1125 * pow(x, 2) - 11.165 * x + 31.019625;
}

void secant_method_standard(double a0, double b0, double eps = 0.001) {
    double xm = a0, xk = b0;
    double pm = df(xm), pk = df(xk);
    int evals = 2; // 2 обчислення на старті

    cout << " МЕТОД СІЧНИХ " << endl;
    if (pm * pk >= 0) {
        cout << "помилка: f'(a0) та f'(b0) мають однаковий знак!" << endl;
        return;
    }

    cout << left << setw(4) << "k" << setw(12) << "x_{k-1}" << setw(12) << "f'(x_{k-1})"
         << setw(12) << "x_k" << setw(12) << "f'(x_k)" << setw(12) << "x_{k+1}"
         << setw(14) << "f'(x_{k+1})" << setw(14) << "f(x_{k+1})" << "Дія" << endl;
    cout << string(105, '-') << endl;

    for (int k = 1; k <= 50; k++) {
        if (abs(pk - pm) < 1e-12) {
            cout << "помилка: ділення на нуль (січна горизонтальна)." << endl;
            break;
        }

        double xn = (pk * xm - pm * xk) / (pk - pm); // (3.16)
        double pn = df(xn);
        double fn = f(xn);
        evals++;

        cout << left << setw(4) << k << fixed << setprecision(6) 
             << setw(12) << xm << setw(12) << pm << setw(12) << xk << setw(12) << pk 
             << setw(12) << xn << setw(14) << pn << setw(14) << fn;

        if (abs(pn) <= eps) {
            cout << "СТОП (|f'|<=eps)" << endl;
            cout << string(105, '-') << endl;
            cout << "оптимальна точка x* = " << xn << ", f(x*) = " << fn << endl;
            cout << "ітерацій: " << k << ", звернень до f': " << evals << endl << endl;
            break;
        }

        if (pn > 0) {
            cout << "xk <- x_{k+1}" << endl;
            xk = xn; pk = pn; //замінюємо праву
        } else {
            cout << "x_{k-1} <- x_{k+1}" << endl;
            xm = xn; pm = pn; //замінюємо ліву
        }
    }
}

void secant_method_clean(double a0, double b0, double eps = 0.001) {
    double xm = a0, xk = b0;
    double pm = df(xm), pk = df(xk);
    int evals = 2;

    cout << " ЕКСПЕРИМЕНТ: «ЧИСТІ» СІЧНІ " << endl;
    cout << left << setw(4) << "k" << setw(12) << "x_{k-1}" << setw(12) << "x_k"
         << setw(12) << "x_{k+1}" << setw(14) << "f'(x_{k+1})" << endl;
    cout << string(60, '-') << endl;

    for (int k = 1; k <= 50; k++) {
        if (abs(pk - pm) < 1e-12) break;

        double xn = (pk * xm - pm * xk) / (pk - pm);
        double pn = df(xn);
        evals++;

        cout << left << setw(4) << k << fixed << setprecision(6) 
             << setw(12) << xm << setw(12) << xk << setw(12) << xn << setw(14) << pn << endl;

        if (abs(pn) <= eps) {
            cout << string(60, '-') << endl;
            cout << "оптимальна точка x* = " << xn << endl;
            cout << "ітерацій: " << k << ", звернень до f': " << evals << endl << endl;
            break;
        }

        //зсуваємо пару
        xm = xk; pm = pk;
        xk = xn; pk = pn;
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

    //  (Таблиця 4.11)
    secant_method_standard(a0, b0, eps);

    //  (Пункт 7)
    secant_method_clean(a0, b0, eps);

    return 0;
}