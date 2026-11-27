#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <functional>
#ifdef _WIN32
#include <windows.h>
#endif

using namespace std;

struct Vec2 {
    double x1, x2;
    Vec2 operator+(const Vec2 &o) const { return {x1 + o.x1, x2 + o.x2}; }
    Vec2 operator-(const Vec2 &o) const { return {x1 - o.x1, x2 - o.x2}; }
    Vec2 operator*(double s) const { return {x1 * s, x2 * s}; }
};

double dot(const Vec2 &a, const Vec2 &b) { return a.x1 * b.x1 + a.x2 * b.x2; }
double norm2(const Vec2 &a) { return sqrt(dot(a, a)); }

double F7(const Vec2 &x) {
    return 3.0 * x.x1 * x.x1 + x.x2 * x.x2 + 2.0 * x.x1 * x.x2 + 16.0 * x.x1 + 4.0 * x.x2;
}
Vec2 grad7(const Vec2 &x) {
    return {6.0 * x.x1 + 2.0 * x.x2 + 16.0, 2.0 * x.x1 + 2.0 * x.x2 + 4.0};
}

//функція Розенброка
double F_ros(const Vec2 &x) {
    double t1 = 1.0 - x.x1;
    double t2 = x.x2 - x.x1 * x.x1;
    return t1 * t1 + 100.0 * t2 * t2;
}
Vec2 grad_ros(const Vec2 &x) {
    double t2 = x.x2 - x.x1 * x.x1;
    return {-2.0 * (1.0 - x.x1) - 400.0 * x.x1 * t2, 200.0 * t2};
}

//  лінійний пошук 
double line_search(function<double(const Vec2&)> F, const Vec2 &x, const Vec2 &d,
                   int &nF_count, double delta0 = 0.01, double eps = 1e-6) {
    auto phi = [&](double h) {
        nF_count++;
        return F({x.x1 + h * d.x1, x.x2 + h * d.x2});
    };

    double h_prev = 0.0, f_prev = phi(0.0);
    double step = delta0;
    double h_curr = step, f_curr = phi(h_curr);
    double a = 0.0, b = 0.0;

    if (f_curr >= f_prev) {
        while (f_curr >= f_prev && step > 1e-12) {
            step /= 2.0;
            h_curr = step;
            f_curr = phi(h_curr);
        }
        a = 0.0;
        b = h_curr * 2.0;
    } else {
        while (true) {
            step *= 2.0;
            double h_next = h_curr + step;
            double f_next = phi(h_next);
            if (f_next >= f_curr || h_next > 1e6) {
                a = h_prev;
                b = h_next;
                break;
            }
            h_prev = h_curr; f_prev = f_curr;
            h_curr = h_next; f_curr = f_next;
        }
    }

    const double tau = (sqrt(5.0) - 1.0) / 2.0;
    double x1 = b - tau * (b - a);
    double x2 = a + tau * (b - a);
    double f1 = phi(x1), f2 = phi(x2);

    while ((b - a) > eps) {
        if (f1 < f2) {
            b = x2; x2 = x1; f2 = f1;
            x1 = b - tau * (b - a);
            f1 = phi(x1);
        } else {
            a = x1; x1 = x2; f1 = f2;
            x2 = a + tau * (b - a);
            f2 = phi(x2);
        }
    }
    return 0.5 * (a + b);
}

struct Res {
    string name;
    int n, nF, nG;
    Vec2 x;
    double Fx;
};

// метод Флетчера-Рівса 
Res fletcher_reeves(const string &name, function<double(const Vec2&)> F,
                    function<Vec2(const Vec2&)> grad, Vec2 x0,
                    double eps = 0.01, bool restart = true, bool steepest = false,
                    bool verbose = false, int max_iter = 5000) {
    Vec2 x = x0;
    Vec2 g = grad(x);
    Vec2 d = g * (-1.0);
    int nF = 0, nG = 1, k = 0;

    Vec2 d0_saved = d, g1_saved = {0, 0}, d1_saved = {0, 0};

    if (verbose) {
        cout << "=== " << name << " ===" << endl;
        cout << left << setw(3) << "k" << setw(22) << "x_(k)" << setw(12) << "F(x_(k))"
             << setw(20) << "grad F" << setw(11) << "||g||2"
             << setw(22) << "d_(k)" << setw(10) << "h_k" << "beta_(k+1)" << endl;
        cout << string(110, '-') << endl;
    }

    for (k = 0; k < max_iter; k++) {
        double gn = norm2(g);
        double Fx = F(x);
        if (gn <= eps) {
            if (verbose) {
                cout << left << setw(3) << k << fixed << setprecision(6)
                     << "(" << x.x1 << "; " << x.x2 << ")   "
                     << setw(12) << Fx << "СТОП (||g||2 = " << scientific << setprecision(2) << gn << " <= eps)" << endl;
            }
            break;
        }

        double h = line_search(F, x, d, nF, 0.01, 1e-6);
        Vec2 x_new = x + d * h;
        Vec2 g_new = grad(x_new);
        nG++;

        double beta = 0.0;
        if (!steepest) {
            beta = dot(g_new, g_new) / dot(g, g);
            if (restart && ((k + 1) % 2 == 0)) beta = 0.0;
        }

        if (verbose) {
            cout << left << setw(3) << k << fixed << setprecision(5)
                 << "(" << setw(8) << x.x1 << ";" << setw(8) << x.x2 << ")   "
                 << setw(12) << Fx
                 << "(" << setprecision(3) << setw(6) << g.x1 << ";" << setw(7) << g.x2 << ")   "
                 << setprecision(5) << setw(11) << gn
                 << "(" << setprecision(3) << setw(7) << d.x1 << ";" << setw(7) << d.x2 << ")   "
                 << setprecision(6) << setw(10) << h << beta << endl;
        }

        Vec2 d_new = g_new * (-1.0) + d * beta;
        if (k == 0) { g1_saved = g_new; d1_saved = d_new; }
        x = x_new; g = g_new; d = d_new;
    }

    if (verbose) {
        cout << string(110, '-') << endl;
        cout << fixed << setprecision(6) << "x* = (" << x.x1 << "; " << x.x2
             << "), F(x*) = " << F(x) << ", n = " << k << ", NF = " << nF << ", N_g = " << nG << endl;
        if (!steepest && name.find("F7") != string::npos) {
            // Пункт 4: Три перевірки (H = [[6, 2], [2, 2]])
            Vec2 Hd0 = {6.0 * d0_saved.x1 + 2.0 * d0_saved.x2, 2.0 * d0_saved.x1 + 2.0 * d0_saved.x2};
            cout << "Перевірка 1: <grad(x1), d0> = " << scientific << setprecision(2) << dot(g1_saved, d0_saved) << " (~ 0)" << endl;
            cout << "Перевірка 2: d1^T * H * d0  = " << dot(d1_saved, Hd0) << " (~ 0)" << endl;
            cout << "Перевірка 3: <grad(x1), d1> = " << fixed << setprecision(4) << dot(g1_saved, d1_saved) << " (< 0)" << endl;
        }
        cout << endl;
    }
    return {name, k, nF, nG, x, F(x)};
}

int main() {
#ifdef _WIN32
    SetConsoleCP(CP_UTF8);
    SetConsoleOutputCP(CP_UTF8);
#endif
    Vec2 x0_7 = {3.0, -3.0};
    Vec2 g0 = grad7(x0_7), d0 = g0 * (-1.0);
    Vec2 Hd0 = {6.0 * d0.x1 + 2.0 * d0.x2, 2.0 * d0.x1 + 2.0 * d0.x2};
    int tmp_nF = 0;
    double h_exact = -dot(g0, d0) / dot(d0, Hd0);
    double h_num = line_search(F7, x0_7, d0, tmp_nF);

    cout << "=== ПЕРЕВІРКА ЛІНІЙНОГО ПОШУКУ НА КРОЦІ 0 (F7) ===" << endl;
    cout << fixed << setprecision(8) << "Точне h* = " << h_exact << ", чисельне h = " << h_num
         << ", похибка = " << scientific << setprecision(2) << abs(h_exact - h_num) << endl << endl;

    // 1. Флетчер-Рівс на F7
    Res r_fr7 = fletcher_reeves("Флетчер-Рівс (F7)", F7, grad7, x0_7, 0.01, true, false, true);
    // 2. Антиградієнт + точний пошук на F7
    Res r_sd7 = fletcher_reeves("Антиградієнт + точний пошук (F7)", F7, grad7, x0_7, 0.01, false, true, true);
    // 3. Розенброк з (-1.2; 1)
    Vec2 x0_ros = {-1.2, 1.0};
    Res r_ros_res = fletcher_reeves("Розенброк (перезапуск кожні 2)", F_ros, grad_ros, x0_ros, 0.01, true, false, false);
    Res r_ros_no  = fletcher_reeves("Розенброк (без перезапуску)", F_ros, grad_ros, x0_ros, 0.01, false, false, false);

    cout << "=== ЗВЕДЕНА ТАБЛИЦЯ ===" << endl;
    cout << left << setw(34) << "Метод" << setw(6) << "n" << setw(7) << "NF"
         << setw(7) << "N_g" << setw(8) << "Разом" << "F(x*)" << endl;
    cout << string(72, '-') << endl;
    Res all[4] = {r_fr7, r_sd7, r_ros_res, r_ros_no};
    for (int i = 0; i < 4; i++) {
        cout << left << setw(34) << all[i].name << setw(6) << all[i].n
             << setw(7) << all[i].nF << setw(7) << all[i].nG << setw(8) << (all[i].nF + all[i].nG)
             << scientific << setprecision(2) << all[i].Fx << endl;
    }
    return 0;
}