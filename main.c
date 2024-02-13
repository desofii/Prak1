#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define s 5
#define alpha 1
#define tmin -1

double f (double t, double x);
double f (double t, double x) {
    
    double func;
    func = 1 - alpha * t * t - x * x;
    return func;
}

int sign (double a);
int sign (double a) {
    
    if (a > 0) return 1;
    else if (a < 0) return -1;
    else return 0;
}

double min (double a, double b);
double min (double a, double b) {
    
    if (a < b) return a;
    else return b;
}

double hnew (double eps, double err, double h);
double hnew (double eps, double err, double h) {
    
    double h_new;
    double chi;
    
    if (fabs(err) < 1e-15)
        chi = 10;
    else {
        chi = pow(eps / err, 1.0 / (s + 1));
        if (chi > 10) chi = 10;
        if (chi < 0.1) chi = 0.1;
    }
    h_new = 0.95 * h * chi;
    
    return h_new;
}

double RungeKutta (double h, double x0, double *x1, double t);
double RungeKutta (double h, double x0, double *x1, double t) {

    double k1, k2, k3, k4, k5, k6;
    double err = 0;
    static const double a27 = 1.0 / 27.0;
    static const double a6 = 1.0 / 6.0;
    static const double a336 = 1.0 / 336.0;
    
    k1 = h * f(t, x0);
    k2 = h * f(t + h * 0.5, x0 + k1 * 0.5);
    k3 = h * f(t + h * 0.5, x0 + (k1 + k2) * 0.25);
    k4 = h * f(t + h, x0 - k2 + 2 * k3);
    k5 = h * f(t + 4 * h * a6, x0 + (7 * k1 + 10 * k2 + k4) * a27);
    k6 = h * f(t + h * 0.2, x0 + (28 * k1 - 125 * k2 + 546 * k3 + 54 * k4 - 378 * k5) * 0.0016);
    *x1 = x0 + (k1 + 4 * k3 + k4) * a6;
    err += fabs((42 * k1 + 224 * k3 + 21 * k4 - 162 * k5 - 125 * k6) * a336);
    
    return err;
}

double GlobalError (double h, double t_max, double x0, double eps, double *delta, int ind);
double GlobalError (double h, double t_max, double x0, double eps, double *delta, int ind) {
    
    double err;
    double delta1;
    double t;
    double h_new;
    double x_prev;
    double x;
    FILE *file;
    file = fopen("output.txt","w");
    x = x0;
    h_new = h;
    t = tmin;
    delta1 = 0;
    while (t < t_max) {
        x_prev = x;
        if (ind == 0) {
            fprintf(file, "%lf %lf", t, x_prev);
            fprintf(file, "\n");
        }
        err = RungeKutta(h_new, x_prev, &x, t);
        while (err > eps) {
            h_new = min(hnew(eps, err, h_new), t_max - t);
            err = RungeKutta(h_new, x_prev, &x, t);
        }
        delta1 = err + delta1 * exp(h_new * (-2 * x_prev)); // интеграл вычисляем по формуле левых прямоугольников
        t += h_new;
        h_new = min(hnew(eps, err, h_new), t_max - t);
    }
    if (ind == 0) {
        fprintf(file, "%lf %lf", t_max, x);
    }
    *delta = delta1;
    fclose(file);
    return x;
}

int Chord_method (double h, double *x, double eps, double a, double b, double eps0);
int Chord_method (double h, double *x, double eps, double a, double b, double eps0) {
    
    double global_err = 0;
    double res[2];
    double fa, fb, fc;
    int i = 0;
    int sa, sb, sc;
    
    res[0] = a;
    res[1] = b;
    
    fa = GlobalError(h, 1.0, a, eps, &global_err, 1) - a;
    fb = GlobalError(h, 1.0, b, eps, &global_err, 1) - b;
    
    if(fabs(fa) < eps) {
        *x = a;
        return 0;
    }
    if(fabs(fb) < eps) {
        *x = b;
        return 0;
    }
    sa = sign(fa);
    sb = sign(fb);
    if(sa * sb != -1) {
        return -1;
    }
    while(fabs(res[1] - res[0]) > eps) {
        i =! i;
        res[i] = a - ((b - a) * fa) / (fb - fa);
        fc = GlobalError(h, 1.0, res[i], eps, &global_err, 1) - res[i];
        sc = sign(fc);
        if(fabs(fc) < eps) {
            *x = res[i];
            return 0;
        }
        if(sc * sa == 1) {
            a = res[i];
            fa = fc;
        }
        else {
            b = res[i];
            fb = fc;
        }
    }
    *x = res[i];
    return 0;
}
 
int main() {
    double h = 0.1;
    double x0, x;
    double global_err = 0;
    
    if(Chord_method(h, &x0, 1e-11, 0.0, 2.0, 1e-15) == -1)
        printf("!!!!!");
    x = GlobalError(h, 1.0, x0, 1e-11, &global_err, 0);
    
    printf("x = %.12f\n", x);
    printf("global_error = %.12f\n", global_err);
    
    return 0;
}
