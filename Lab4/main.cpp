#include <iostream>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_monte_miser.h>
#include <sys/param.h>

double exact = 2;

double function(double x) {
    return 1 / sqrt(x);
}

double function2(double x) {
    return x * x;
}

double randomPoint(double a, double b) {
    return a + drand48() * (b - a);
}

int funcIn(double x, double y) {
    if (y > 0 && y <= function2(x) || y <= 0 && y >= function2(x))
        return 1;
    return 0;
}

double hit_and_miss(double start, double end, double up, double down, int samples) {
    double x;
    double y;
    if (end > 1.0) end = 1.0;
    int hits = 0;
    for (int i = 0; i < samples; ++i) {
        x = randomPoint(start, end);
        y = randomPoint(down, up);
        hits += funcIn(x, y);
    }
    return (double) hits / samples * abs(down - up) * abs(end - start);
}

void display_results(char *title, double result, double error) {
    printf("%s ==================\n", title);
    printf("result = % .6f\n", result);
    printf("sigma  = % .6f\n", error);
    printf("exact  = % .6f\n", exact);
    printf("error  = % .6f = %.2g sigma\n", result - exact,
           fabs(result - exact) / error);

}

void montecarlo_plain(gsl_monte_function G, double *xl, double *xu, size_t calls, gsl_rng *r) {
    double res, err;
    gsl_monte_plain_state *s = gsl_monte_plain_alloc(2);
    gsl_monte_plain_integrate(&G, xl, xu, 2, calls, r, s,
                              &res, &err);
    gsl_monte_plain_free(s);

    display_results("Monte Carlo Plain", res, err);
}

void montecarlo_vegas(gsl_monte_function G, double *xl, double *xu, size_t calls, gsl_rng *r) {
    double res, err;

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(2);

    gsl_monte_vegas_integrate(&G, xl, xu, 2, calls, r, s,
                              &res, &err);
    display_results("Monte Carlo Vegas", res, err);


    gsl_monte_vegas_free(s);
}

void montecarlo_miser(gsl_monte_function G, double *xl, double *xu, size_t calls, gsl_rng *r) {
    double res, err;

    gsl_monte_miser_state *s = gsl_monte_miser_alloc(2);
    gsl_monte_miser_integrate(&G, xl, xu, 2, calls, r, s,
                              &res, &err);
    display_results("Monte Carlo Miser", res, err);

    gsl_monte_miser_free(s);
}

int main() {
    double xl[2] = {0.0, 0.0};
    double xu[2] = {1.0, 1.0};
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_monte_function G = {reinterpret_cast<double (*)(double *, size_t, void *)>(&function), 2, 0};
    size_t calls = 5000;
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    montecarlo_vegas(G, xl, xu, calls, r);

    for (int i = 100; i < 50000; i += 100) {
        montecarlo_vegas(G, xl, xu, static_cast<size_t>(i), r);
    }
    
    return 0;
}