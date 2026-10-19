// LD_PRELOAD into the old binary: replaces main() and prints tools::NaturalMIDCalculator results.
// usage: metabolitedetector_app ISOTOPES.tsv N FORMULA...
#include <dlfcn.h>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
namespace {
uintptr_t base;
constexpr uintptr_t MAIN = 0x3b1f2d, CTOR = 0x6d931c, MID = 0x6d941e, EXACT = 0x6d9c6a;
template <class F> F fn(uintptr_t a) { return reinterpret_cast<F>(base + a); }
int my_main(int argc, char** argv, char**) {
    std::ifstream f(argv[1]);
    std::stringstream ss; ss << f.rdbuf();
    std::string table = ss.str();
    alignas(16) static char calc[4096];
    fn<void (*)(void*, const std::string&)>(CTOR)(calc, table);
    unsigned long n = std::stoul(argv[2]);
    for (int i = 3; i < argc; i++) {
        std::string formula = argv[i];
        try {
            double m = fn<double (*)(const void*, const std::string&)>(EXACT)(calc, formula);
            std::vector<double> mid = fn<std::vector<double> (*)(const void*, const std::string&, unsigned long)>(MID)(calc, formula, n);
            std::printf("%s\t%.10g", formula.c_str(), m);
            for (double v : mid) std::printf("\t%.12g", v);
            std::printf("\n");
        } catch (...) {
            std::printf("%s\terror\n", formula.c_str());
        }
    }
    return 0;
}
}
extern "C" int __libc_start_main(int (*main)(int, char**, char**), int argc, char** argv, void (*init)(),
                                 void (*fini)(), void (*rtld_fini)(), void* stack_end) {
    base = reinterpret_cast<uintptr_t>(main) - MAIN;
    using Start = int (*)(int (*)(int, char**, char**), int, char**, void (*)(), void (*)(), void (*)(), void*);
    auto real = reinterpret_cast<Start>(dlsym(RTLD_NEXT, "__libc_start_main"));
    return real(my_main, argc, argv, init, fini, rtld_fini, stack_end);
}
