#include <iostream>
#include <math.h>
#include <iomanip>

#ifdef _WIN32
#include <windows.h>
#endif

using namespace std;

// ф-ія для 7 варіанту
double f(double x) {
    // f(x) = 1.75 * (x + 2.3) * (x - 1.05) * (x - 1.7) * (x - 3) + 1
    return 1.75 * (x + 2.3) * (x - 1.05) * (x - 1.7) * (x - 3.0) + 1.0;
}

// алгоритм Свенна
void methodSvenn(double x0, double delta, double &a_out, double &b_out) {
    cout << "алгоритм Свенна" << endl;
    cout << "початкова точка x0 = " << x0 << ", крок delta = " << delta << endl;
    
    double x_k = x0;
    double f_k = f(x_k);
    double f_right = f(x_k + delta);
    double f_left = f(x_k - delta);
    
    int evals = 3; 
    
    // напрямок
    if (f_left >= f_k && f_k >= f_right) {
        // праворуч
    } else if (f_left <= f_k && f_k <= f_right) {
        delta = -delta; //ліворуч
    } else if (f_left >= f_k && f_k <= f_right) {
        a_out = x_k - delta;
        b_out = x_k + delta;
        cout << "мінімум локалізовано на першому кроці: [" << a_out << ", " << b_out << "]" << endl;
        return;
    }

    double x_next = x_k + delta;
    double f_next = f(x_next);
    evals++;
    
    int k = 1;
    double a = x_k - delta; 
    
    cout << left << setw(5) << "k" << setw(15) << "x_k" << setw(15) << "f(x_k)" << setw(15) << "крок" << endl;
    cout << left << setw(5) << 0 << setw(15) << x_k << setw(15) << f_k << setw(15) << delta << endl;
    
    while (f_next < f_k) {
        cout << left << setw(5) << k << setw(15) << x_next << setw(15) << f_next << setw(15) << pow(2, k) * delta << endl;
        a = x_k;
        x_k = x_next;
        f_k = f_next;
        
        double step = pow(2, k) * delta;
        x_next = x_k + step;
        f_next = f(x_next);
        evals++;
        k++;
    }
    
    cout << left << setw(5) << k << setw(15) << x_next << setw(15) << f_next << setw(15) << "-" << endl;
    
    a_out = min(a, x_next);
    b_out = max(a, x_next);
    cout << "\nзнайдений інтервал: [" << a_out << ", " << b_out << "]" << endl;
    cout << "обчислень функції: " << evals << "\n\n";
}

// метод дихотомії
void methodDichotomy(double a, double b, double sigma, double epsilon) {
    cout << "метод дихотомії " << endl;
    cout << "інтервал: [" << a << ", " << b << "], sigma = " << sigma << ", eps = " << epsilon << endl;
    
    int k = 0;
    int evals = 0;
    double L = b - a;
    
    cout << left << setw(5) << "k" << setw(12) << "a" << setw(12) << "b" << setw(12) << "x1" << setw(12) << "x2" << setw(12) << "L_k" << endl;
    
    while (L > sigma) {
        double x1 = (a + b) / 2.0 - epsilon / 2.0;
        double x2 = x1 + epsilon;
        
        double f1 = f(x1);
        double f2 = f(x2);
        evals += 2;
        
        cout << left << setw(5) << k << setw(12) << a << setw(12) << b << setw(12) << x1 << setw(12) << x2 << setw(12) << L << endl;
        
        if (f1 < f2) {
            b = x2;
        } else {
            a = x1;
        }
        L = b - a;
        k++;
    }
    
    cout << "\nмінімум x* = " << (a + b) / 2.0 << endl;
    cout << "обчислень функції: " << evals << "\n\n";
}

// метод половинного поділу
void methodHalfDivision(double a, double b, double sigma) {
    cout << "метод половинного поділу" << endl;
    cout << "інтервал: [" << a << ", " << b << "], sigma = " << sigma << endl;
    
    int k = 0;
    int evals = 0;
    double L = b - a;
    
    double xm = (a + b) / 2.0;
    double fm = f(xm);
    evals++;
    
    cout << left << setw(5) << "k" << setw(12) << "x1" << setw(12) << "xm" << setw(12) << "x2" << setw(12) << "L_k" << endl;
    
    while (L > sigma) {
        double x1 = a + L / 4.0;
        double x2 = b - L / 4.0;
        
        double f1 = f(x1);
        double f2 = f(x2);
        evals += 2;
        
        cout << left << setw(5) << k << setw(12) << x1 << setw(12) << xm << setw(12) << x2 << setw(12) << L << endl;
        
        if (f1 < fm) {
            b = xm;
            xm = x1;
            fm = f1;
        } else if (f2 < fm) {
            a = xm;
            xm = x2;
            fm = f2;
        } else {
            a = x1;
            b = x2;
        }
        
        L = b - a;
        k++;
    }
    
    cout << "\nмінімум x* = " << (a + b) / 2.0 << endl;
    cout << "обчислень функції: " << evals << endl;
}

int main() {
#ifdef _WIN32
    SetConsoleCP(CP_UTF8);
    SetConsoleOutputCP(CP_UTF8);
#endif

    cout << fixed << setprecision(4);
    
    //дані варіанту 7
    double x0 = -1.35;
    double delta = 0.1;
    double sigma = 0.05;
    double epsilon = 0.01; 
    
    double a, b; // межі інтервалу
    
    //виконання методів
    methodSvenn(x0, delta, a, b);
    methodDichotomy(a, b, sigma, epsilon);
    methodHalfDivision(a, b, sigma);
    
    return 0;
}