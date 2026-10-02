#include <iostream>
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
    //інтервал з методу Свенна
    double a = -1.3500; 
    double b = -1.0500;  
    double L0 = b - a;   
    double eps = 0.05;   

    //обчислення сталої золотого перерізу
    double tau = (sqrt(5.0) - 1.0) / 2.0; 
    
    //початкові точки
    double x1 = a + (1.0 - tau) * (b - a);
    double x2 = a + tau * (b - a);
    
    double f1 = f(x1);
    double f2 = f(x2);
    
    int k = 1;
    int nf = 2; //лічильник обчислень функції

    cout << "k\ta\t\tb\t\tx1\t\tx2\t\tf(x1)\t\tf(x2)\t\tL\n";
    cout << fixed << setprecision(5);

    //цикл методу золотого перерізу
    while (b - a > eps) {
        cout << k << "\t" << a << "\t" << b << "\t" 
             << x1 << "\t" << x2 << "\t" << f1 << "\t" << f2 << "\t" << (b - a) << "\n";
             
        if (f1 <= f2) {
            b = x2;
            x2 = x1;
            f2 = f1; 
            x1 = a + (1.0 - tau) * (b - a);
            f1 = f(x1); 
        } else {
            a = x1;
            x1 = x2;
            f1 = f2; 
            x2 = a + tau * (b - a);
            f2 = f(x2); 
        }
        k++;
        nf++;
    }
    
    double x_star = (a + b) / 2.0;
    double L_final = b - a;
    
    // коефіцієнт стиснення
    double eta = pow((L_final / L0), 1.0 / nf);

    cout << "--------------------------------------------------------\n";
    cout << "оптимальна точка x* = " << x_star << "\n";
    cout << "значення функції f(x*) = " << f(x_star) << "\n";
    cout << "кількість ітерацій: " << k - 1 << "\n";
    cout << "кількість обчислень функції Nf: " << nf << "\n";
    cout << "коефіцієнт стиснення (eta): " << eta << "\n";

    return 0;
}