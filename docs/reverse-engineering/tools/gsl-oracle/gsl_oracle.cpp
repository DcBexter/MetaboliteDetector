// Black-box oracle for the GSL linear algebra statically linked into the old MetaboliteDetector 3.3
// (Linux). Loaded with LD_PRELOAD into the old binary like the golden driver: it replaces main() and
// calls functions of the binary by address (syms.inc from `nm`, see run.sh). It only feeds inputs and
// prints outputs; it contains no GSL code. Prototypes are the public API from the GSL manual.
//
// stdin: one command per line, numbers as C literals (decimal or hex float); stdout: one line per
// command with every output double as 16 hex digits (the IEEE bits). Commands:
//   fit  x0 x1 x2 y0 y1 y2          tools::Statistics::polynomialRegression(X, Y, 2) -> c0 c1 c2
//   mfl  n p X(n*p) y(n)            gsl_multifit_linear -> c(p) cov(p*p) chisq
//   bal  m n A(m*n)                 gsl_linalg_balance_columns -> A D
//   svd  m n A(m*n)                 gsl_linalg_SV_decomp -> U S V
//   svdm m n A(m*n)                 gsl_linalg_SV_decomp_mod -> U S V
//   bid  m n A(m*n)                 gsl_linalg_bidiag_decomp -> A tau_U tau_V
//   hh   n v(n)                     gsl_linalg_householder_transform -> tau v
//   hm   tau n v(n) m k A(m*k)      gsl_linalg_householder_hm -> A
//   mh   tau n v(n) m k A(m*k)      gsl_linalg_householder_mh -> A
//   sol  m n U(m*n) S(n) V(n*n) b(m) gsl_linalg_SV_solve -> x
//   nrm2 n v(n)                     cblas_dnrm2 -> r
//   dot  n x(n) y(n)                cblas_ddot -> r
//   xhyp k a b                      k-th internal xhypot helper of the binary -> r
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dlfcn.h>
#include <string>
#include <vector>

namespace {

uintptr_t base;
template <class F> F fn(uintptr_t addr) { return reinterpret_cast<F>(base + addr); }

#define SYM(name, addr) constexpr uintptr_t name = addr;
#include "syms.inc"
#undef SYM

using M = void*;  // gsl_matrix*
using V = void*;  // gsl_vector*

M m_alloc(size_t a, size_t b) { return fn<M (*)(size_t, size_t)>(MATRIX_ALLOC)(a, b); }
V v_alloc(size_t n) { return fn<V (*)(size_t)>(VECTOR_ALLOC)(n); }
void m_free(M m) { fn<void (*)(M)>(MATRIX_FREE)(m); }
void v_free(V v) { fn<void (*)(V)>(VECTOR_FREE)(v); }
double mg(M m, size_t i, size_t j) { return fn<double (*)(M, size_t, size_t)>(MATRIX_GET)(m, i, j); }
void ms(M m, size_t i, size_t j, double x) { fn<void (*)(M, size_t, size_t, double)>(MATRIX_SET)(m, i, j, x); }
double vg(V v, size_t i) { return fn<double (*)(V, size_t)>(VECTOR_GET)(v, i); }
void vs(V v, size_t i, double x) { fn<void (*)(V, size_t, double)>(VECTOR_SET)(v, i, x); }

std::vector<double> nums;
size_t pos;
double num() { return pos < nums.size() ? nums[pos++] : 0.0; }
size_t count() { return static_cast<size_t>(num()); }

void put(double x) {
    uint64_t u;
    std::memcpy(&u, &x, 8);
    std::printf(" %016llx", static_cast<unsigned long long>(u));
}
M read_m(size_t a, size_t b) {
    M m = m_alloc(a, b);
    for (size_t i = 0; i < a; i++)
        for (size_t j = 0; j < b; j++) ms(m, i, j, num());
    return m;
}
V read_v(size_t n) {
    V v = v_alloc(n);
    for (size_t i = 0; i < n; i++) vs(v, i, num());
    return v;
}
void put_m(M m, size_t a, size_t b) {
    for (size_t i = 0; i < a; i++)
        for (size_t j = 0; j < b; j++) put(mg(m, i, j));
}
void put_v(V v, size_t n) {
    for (size_t i = 0; i < n; i++) put(vg(v, i));
}

void command(const std::string& op) {
    if (op == "fit") {
        std::vector<double> x(3), y(3);
        for (auto& v : x) v = num();
        for (auto& v : y) v = num();
        using Fit = std::vector<double> (*)(std::vector<double>&, std::vector<double>&, unsigned);
        for (double c : fn<Fit>(POLY_REGRESSION)(x, y, 2)) put(c);
    } else if (op == "mfl") {
        size_t n = count(), p = count();
        M x = read_m(n, p);
        V y = read_v(n), c = v_alloc(p);
        M cov = m_alloc(p, p);
        void* w = fn<void* (*)(size_t, size_t)>(MFL_ALLOC)(n, p);
        double chisq = 0;
        fn<int (*)(M, V, V, M, double*, void*)>(MFL)(x, y, c, cov, &chisq, w);
        put_v(c, p);
        put_m(cov, p, p);
        put(chisq);
        fn<void (*)(void*)>(MFL_FREE)(w);
        m_free(x), m_free(cov), v_free(y), v_free(c);
    } else if (op == "bal") {
        size_t m = count(), n = count();
        M a = read_m(m, n);
        V d = v_alloc(n);
        fn<int (*)(M, V)>(BALANCE_COLUMNS)(a, d);
        put_m(a, m, n);
        put_v(d, n);
        m_free(a), v_free(d);
    } else if (op == "svd" || op == "svdm") {
        size_t m = count(), n = count();
        M a = read_m(m, n), v = m_alloc(n, n);
        V s = v_alloc(n), work = v_alloc(n);
        if (op == "svd") {
            fn<int (*)(M, M, V, V)>(SV_DECOMP)(a, v, s, work);
        } else {
            M x = m_alloc(n, n);
            fn<int (*)(M, M, M, V, V)>(SV_DECOMP_MOD)(a, x, v, s, work);
            m_free(x);
        }
        put_m(a, m, n);
        put_v(s, n);
        put_m(v, n, n);
        m_free(a), m_free(v), v_free(s), v_free(work);
    } else if (op == "bid") {
        size_t m = count(), n = count();
        M a = read_m(m, n);
        V tu = v_alloc(n), tv = v_alloc(n - 1);  // the manual: tau_V has length N - 1
        fn<int (*)(M, V, V)>(BIDIAG_DECOMP)(a, tu, tv);
        put_m(a, m, n);
        put_v(tu, n);
        put_v(tv, n - 1);
        m_free(a), v_free(tu), v_free(tv);
    } else if (op == "hh") {
        size_t n = count();
        V v = read_v(n);
        put(fn<double (*)(V)>(HOUSEHOLDER_TRANSFORM)(v));
        put_v(v, n);
        v_free(v);
    } else if (op == "hm" || op == "mh") {
        double tau = num();
        size_t n = count();
        V v = read_v(n);
        size_t a = count(), b = count();
        M m = read_m(a, b);
        fn<int (*)(double, V, M)>(op == "hm" ? HOUSEHOLDER_HM : HOUSEHOLDER_MH)(tau, v, m);
        put_m(m, a, b);
        v_free(v), m_free(m);
    } else if (op == "sol") {
        size_t m = count(), n = count();
        M u = read_m(m, n);
        V s = read_v(n);
        M v = read_m(n, n);
        V b = read_v(m), x = v_alloc(n);
        fn<int (*)(M, M, V, V, V)>(SV_SOLVE)(u, v, s, b, x);
        put_v(x, n);
        m_free(u), m_free(v), v_free(s), v_free(b), v_free(x);
    } else if (op == "nrm2" || op == "dot") {
        size_t n = count();
        std::vector<double> x(n), y(n);
        for (auto& v : x) v = num();
        if (op == "dot") {
            for (auto& v : y) v = num();
            put(fn<double (*)(int, const double*, int, const double*, int)>(CBLAS_DDOT)(int(n), x.data(), 1, y.data(), 1));
        } else {
            put(fn<double (*)(int, const double*, int)>(CBLAS_DNRM2)(int(n), x.data(), 1));
        }
    } else if (op == "xhyp") {  // the binary's internal xhypot helpers, called as a black box
        static const uintptr_t xh[] = {XHYPOT0, XHYPOT1, XHYPOT2, XHYPOT3};
        size_t k = count();
        double a = num(), b = num();
        put(fn<double (*)(double, double)>(xh[k])(a, b));
    } else {
        std::printf(" unknown");
    }
}

int oracle_main(int, char**, char**) {
    fn<void* (*)()>(ERROR_HANDLER_OFF)();
    static char line[1 << 20];
    while (std::fgets(line, sizeof line, stdin)) {
        char* p = line;
        while (*p == ' ' || *p == '\t') p++;
        char* end = p;
        while (*end && *end != ' ' && *end != '\t' && *end != '\n') end++;
        std::string op(p, end);
        if (op.empty()) continue;
        nums.clear();
        pos = 0;
        for (char* q = end;;) {
            char* r;
            double v = std::strtod(q, &r);
            if (r == q) break;
            nums.push_back(v);
            q = r;
        }
        std::printf("%s", op.c_str());
        command(op);
        std::printf("\n");
    }
    return 0;
}

}  // namespace

extern "C" int __libc_start_main(int (*main)(int, char**, char**), int argc, char** argv, void (*init)(),
                                 void (*fini)(), void (*rtld_fini)(), void* stack_end) {
    base = reinterpret_cast<uintptr_t>(main) - MAIN;
    using Start = int (*)(int (*)(int, char**, char**), int, char**, void (*)(), void (*)(), void (*)(), void*);
    auto real = reinterpret_cast<Start>(dlsym(RTLD_NEXT, "__libc_start_main"));
    return real(oracle_main, argc, argv, init, fini, rtld_fini, stack_end);
}
