// LD_PRELOAD into the old binary (sha256 7883d70c…4794): replaces main() and runs the MID result
// window's group statistics and text exports (gcms::ReplicateMIDData, MIDSingleResultWidget) on
// the per-run MIDs of gcms::MIDCompoundSeries. Same case file as old_midseries.cpp; see make-midstats.sh.
// usage: QT_QPA_PLATFORM=offscreen metabolitedetector_app ISOTOPES.tsv LIBRARY.lbr CASES
// Output per case: "stat NAME MI MEAN STDERR" for every mass isotopomer of group GL (all runs),
// then the two text exports between "begin data|stats NAME" and "end" lines, or "error NAME".
#include <dlfcn.h>
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <list>
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
                    CALC_CTOR = 0x6d931c, MCS_CTOR = 0x6edb74, RMD_CTOR = 0x6f586a, SET_NAMES = 0x6f37f8,
                    SET_GROUPS = 0x6f58d4, GROUP_MEAN = 0x6f6052, GROUP_STDERR = 0x6f64ee,
                    WIDGET_CTOR = 0x4ed8c6, OUT_DATA = 0x4ef548, OUT_STATS = 0x4f0090;
template <class F> F fn(uintptr_t a) { return reinterpret_cast<F>(base + a); }
template <class F> F qt(const char* name) { return reinterpret_cast<F>(dlsym(RTLD_DEFAULT, name)); }

struct Case {
    std::string name, tracer;
    size_t compound = 0, runs = 0;
    bool flag = false;
    struct Frag { unsigned long lo, hi, size; std::string formula; };
    std::vector<Frag> frags;
    std::vector<std::tuple<size_t, float, double>> ints;
};

// Writes one text export of the result widget to stdout through a QTextStream(FILE*).
void dump(void* widget, uintptr_t out, const std::string& what, const std::string& name) {
    std::printf("begin\t%s\t%s\n", what.c_str(), name.c_str());
    std::fflush(stdout);
    alignas(16) char ts[64];
    qt<void (*)(void*, FILE*, int)>("_ZN11QTextStreamC1EP8_IO_FILE6QFlagsIN9QIODevice12OpenModeFlagEE")(ts, stdout, 2);
    fn<void (*)(void*, void*)>(out)(widget, ts);
    qt<void (*)(void*)>("_ZN11QTextStream5flushEv")(ts);
    qt<void (*)(void*)>("_ZN11QTextStreamD1Ev")(ts);
    std::printf("\nend\n");
    std::fflush(stdout);
}

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
    fn<void (*)(void*, const std::set<float>&)>(SET_QUANT_IONS)(c, quant);
    // The widget keeps pointers into the series, so they live until the process ends.
    auto* cs = new (std::align_val_t(16)) char[1 << 14];
    fn<void (*)(void*, const void*, bool)>(CS_CTOR)(cs, c, false);
    fn<void (*)(void*, unsigned long)>(CS_NULL)(cs, k.runs);
    for (auto& [r, ion, v] : k.ints) fn<void (*)(void*, float, unsigned long, double)>(CS_SET)(cs, ion, r, v);
    auto* mcs = new (std::align_val_t(16)) char[0xd8];
    fn<void (*)(void*, const void*, const void*, const void*, const std::string&, bool)>(MCS_CTOR)(
        mcs, c, cs, calc, k.tracer, k.flag);
    // As MIDWizard::midCalcFinished: one ReplicateMIDData, sample names, the wizard's group GL.
    auto* rmd = new (std::align_val_t(16)) char[0x60];
    fn<void (*)(void*, void*)>(RMD_CTOR)(rmd, mcs);
    std::vector<std::string> names;
    for (size_t r = 0; r < k.runs; r++) names.push_back("/data/run" + std::to_string(r) + ".cdf");
    fn<void (*)(void*, const std::vector<std::string>&)>(SET_NAMES)(rmd, names);
    std::map<std::string, std::list<std::string>> groups{{"GL", {names.begin(), names.end()}}};
    fn<void (*)(void*, const std::map<std::string, std::list<std::string>>&)>(SET_GROUPS)(rmd, groups);
    // MIDCompoundSeries+0x90: one MID vector per mass isotopomer, fragments in order.
    auto* mis = reinterpret_cast<const std::vector<const void*>*>(mcs + 0x90);
    const std::string gl = "GL";
    for (size_t i = 0; i < mis->size(); i++) {
        double mean = fn<double (*)(void*, const std::string&, const void*, unsigned long, bool)>(GROUP_MEAN)(
            rmd, gl, mcs, i, true);
        double se = fn<double (*)(void*, const std::string&, const void*, unsigned long, bool)>(GROUP_STDERR)(
            rmd, gl, mcs, i, true);
        std::printf("stat\t%s\t%zu\t%.17g\t%.17g\n", k.name.c_str(), i, mean, se);
    }
    auto* widget = new (std::align_val_t(16)) char[0x90];
    fn<void (*)(void*, void*, void*)>(WIDGET_CTOR)(widget, rmd, nullptr);
    dump(widget, OUT_DATA, "data", k.name);
    dump(widget, OUT_STATS, "stats", k.name);
}

int my_main(int argc, char** argv, char**) {
    if (argc < 4) { std::fprintf(stderr, "usage: %s ISOTOPES.tsv LIBRARY.lbr CASES\n", argv[0]); return 2; }
    static int qargc = 1;
    alignas(16) static char app[256];
    // Version argument as in the old main (0x3b1f75).
    qt<void (*)(void*, int&, char**, int)>("_ZN12QApplicationC1ERiPPci")(app, qargc, argv, 0x50c08);
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
