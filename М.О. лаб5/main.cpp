#include <iostream>
#include <math.h>
#include <iomanip>
#ifdef _WIN32
#include <windows.h>
#endif

using namespace std;
double f(double x) {
    return 1.75 * (x + 2.3) * (x - 1.05) * (x - 1.7) * (x - 3.0) + 1.0;
}

int main() {
    #ifdef _WIN32
    SetConsoleCP(CP_UTF8);
    SetConsoleOutputCP(CP_UTF8);
#endif
    double x1 = -1.35; 
    double dx = 0.15;  
    double sf = 0.01, sx = 0.001; 

    double x2 = x1 + dx;
    double f1 = f(x1), f2 = f(x2);
    double x3 = (f1 > f2) ? (x1 + 2 * dx) : (x1 - dx);
    double f3 = f(x3);
    
    int nf = 3; 

    cout << fixed << setprecision(6);
    cout << "k\tx1\t\tx2\t\tx3\t\tf(x1)\t\tf(x2)\t\tf(x3)\t\ta1\t\ta2\t\txt\t\tf(xt)\n";

    double xs[3] = {x1, x2, x3};
    double fs[3] = {f1, f2, f3};

    for (int k = 1; k <= 50; k++) {
    
        int i_min = 0, i_max = 0;
        for (int i = 1; i < 3; i++) {
            if (fs[i] < fs[i_min]) i_min = i;
            if (fs[i] > fs[i_max]) i_max = i;
        }

        double a1 = (fs[1] - fs[0]) / (xs[1] - xs[0]);
        double a2 = ((fs[2] - fs[0]) / (xs[2] - xs[0]) - a1) / (xs[2] - xs[1]);

        double xt;
        if (a2 > 0) {
            xt = -a1 / (2 * a2) + (xs[0] + xs[1]) / 2.0;
        } else {
            xt = 2 * xs[i_min] - xs[i_max]; // Защита от горба
        }

        double ft = f(xt);
        nf++; 

        cout << k << "\t" << xs[0] << "\t" << xs[1] << "\t" << xs[2] << "\t" 
             << fs[0] << "\t" << fs[1] << "\t" << fs[2] << "\t" 
             << a1 << "\t" << a2 << "\t" << xt << "\t" << ft << "\n";

        if (abs(fs[i_min] - ft) <= sf && abs(xs[i_min] - xt) <= sx) {
            cout << "--------------------------------------------------------\n";
            cout << "оптимальна точка x* = " << xt << "\n";
            cout << "значення функції f(x*) = " << ft << "\n";
            cout << "кількість ітерацій n = " << k << "\n";
            cout << "кількість обчислень Nf = " << nf << "\n";
            break;
        }

        xs[i_max] = xt;
        fs[i_max] = ft;
    }

    return 0;
}