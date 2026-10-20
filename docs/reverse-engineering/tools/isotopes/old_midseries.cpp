// LD_PRELOAD into the old binary (sha256 7883d70c…4794): replaces main() and runs
// gcms::MIDCompoundSeries, the per-run MID of the wizard's "library" and "sum formula" methods, on
// integrals given in a case file. See make-midseries.sh.
// usage: metabolitedetector_app ISOTOPES.tsv LIBRARY.lbr CASES
// CASES, one block per case:
//   case NAME COMPOUND_INDEX TRACER FLAG RUNS   (FLAG 1 = library spectrum, 0 = sum formula; TRACER - = none)
//   frag LO HI FORMULA SIZE                     (FORMULA - = none)
//   int RUN ION VALUE
//   end
// Output: "norm ION VALUE" (normalised library intensity), then per fragment, ion and run
// "mid NAME FRAG ION RUN MID CI R2", or "error NAME".
#include <dlfcn.h>
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
namespace {
uintptr_t base;
constexpr uintptr_t MAIN = 0x3b1f2d, LC_FROM_DISK = 0x620a30, LC_ADD_FRAGMENT = 0x61f232,
                    SET_QUANT_IONS = 0x477296, CS_CTOR = 0x62fef6, CS_NULL = 0x63130a, CS_SET = 0x632cfe,
                    CALC_CTOR = 0x6d931c, MCS_CTOR = 0x6edb74;
template <class F> F fn(uintptr_t a) { return reinterpret_cast<F>(base + a); }

using Map = std::map<float, std::vector<double>>;
struct Fragment { Map mid, ci, r2; };  // MIDCompoundSeries::MIDFragment, ctor 0x6ef794

struct Case {
    std::string name, tracer;
    size_t compound = 0, runs = 0;
    bool flag = false;
    struct Frag { unsigned long lo, hi, size; std::string formula; };
    std::vector<Frag> frags;
    std::vector<std::tuple<size_t, float, double>> ints;
};

void run(const Case& k, const char* lbr, const void* calc) {
    std::vector<void*> lib;
    fn<void (*)(std::vector<void*>*, const char*)>(LC_FROM_DISK)(&lib, lbr);
    void* c = lib.at(k.compound);
    std::set<float> quant;
    for (auto& f : k.frags) {
        fn<void (*)(void*, unsigned long, unsigned long, const std::string&, unsigned long)>(LC_ADD_FRAGMENT)(
            c, f.lo, f.hi, f.formula, f.size);
        for (unsigned long i = std::min(f.lo, f.hi); i <= std::max(f.lo, f.hi); i++) quant.insert(float(i));
    }
    // As MIDWizard::createMIDLibrary: the quant ions are all fragment ions.
    fn<void (*)(void*, const std::set<float>&)>(SET_QUANT_IONS)(c, quant);
    auto* norm = reinterpret_cast<const std::map<float, int>*>(static_cast<char*>(c) + 0x98);
    for (float ion : quant) {
        auto it = norm->find(ion);
        std::printf("norm\t%s\t%g\t%d\n", k.name.c_str(), ion, it == norm->end() ? 0 : it->second);
    }
    alignas(16) static char cs[1 << 14];
    fn<void (*)(void*, const void*, bool)>(CS_CTOR)(cs, c, false);
    fn<void (*)(void*, unsigned long)>(CS_NULL)(cs, k.runs);
    for (auto& [r, ion, v] : k.ints) fn<void (*)(void*, float, unsigned long, double)>(CS_SET)(cs, ion, r, v);
    alignas(16) static char mcs[0xd8];
    fn<void (*)(void*, const void*, const void*, const void*, const std::string&, bool)>(MCS_CTOR)(
        mcs, c, cs, calc, k.tracer, k.flag);
    auto* frags = reinterpret_cast<const std::vector<Fragment>*>(mcs + 0x10);
    for (size_t f = 0; f < frags->size(); f++) {
        const Fragment& g = (*frags)[f];
        for (auto& [ion, mids] : g.mid)
            for (size_t r = 0; r < mids.size(); r++)
                std::printf("mid\t%s\t%zu\t%g\t%zu\t%.17g\t%.17g\t%.17g\n", k.name.c_str(), f, ion, r, mids[r],
                            g.ci.at(ion).at(r), g.r2.at(ion).at(r));
    }
}

int my_main(int argc, char** argv, char**) {
    if (argc < 4) { std::fprintf(stderr, "usage: %s ISOTOPES.tsv LIBRARY.lbr CASES\n", argv[0]); return 2; }
    std::ifstream t(argv[1]);
    std::stringstream ss; ss << t.rdbuf();
    alignas(16) static char calc[4096];
    fn<void (*)(void*, const std::string&)>(CALC_CTOR)(calc, ss.str());
    std::ifstream in(argv[3]);
    std::string line;
    Case k;
    while (std::getline(in, line)) {
        std::istringstream w(line);
        std::string tag;
        w >> tag;
        if (tag == "case") {
            k = Case();
            int flag;
            w >> k.name >> k.compound >> k.tracer >> flag >> k.runs;
            if (k.tracer == "-") k.tracer.clear();
            k.flag = flag != 0;
        } else if (tag == "frag") {
            Case::Frag f;
            w >> f.lo >> f.hi >> f.formula >> f.size;
            if (f.formula == "-") f.formula.clear();
            k.frags.push_back(f);
        } else if (tag == "int") {
            size_t r; float ion; double v;
            w >> r >> ion >> v;
            k.ints.emplace_back(r, ion, v);
        } else if (tag == "end") {
            try { run(k, argv[2], calc); } catch (...) { std::printf("error\t%s\n", k.name.c_str()); }
            std::fflush(stdout);
        }
    }
    return 0;
}
}  // namespace
extern "C" int __libc_start_main(int (*main)(int, char**, char**), int argc, char** argv, void (*init)(),
                                 void (*fini)(), void (*rtld_fini)(), void* stack_end) {
    base = reinterpret_cast<uintptr_t>(main) - MAIN;
    using Start = int (*)(int (*)(int, char**, char**), int, char**, void (*)(), void (*)(), void (*)(), void*);
    auto real = reinterpret_cast<Start>(dlsym(RTLD_NEXT, "__libc_start_main"));
    return real(my_main, argc, argv, init, fini, rtld_fini, stack_end);
}
