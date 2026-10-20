// LD_PRELOAD into the old binary (sha256 7883d70c…4794): replaces main() and runs the MID wizard's
// group methods (MIDWizard::Calculate, 0 uncorrected, 1 by unlabelled reference runs) on integrals
// given in a case file, then the group result window's CSV export (MIDResultWidget). See make-midseries.sh.
// usage: QT_QPA_PLATFORM=offscreen metabolitedetector_app ISOTOPES.tsv LIBRARY.lbr CASES
// CASES, one block per case (GL, GUL: comma-separated runs of the labelled and reference group, - = none):
//   case NAME COMPOUND_INDEX METHOD M1 GL GUL RUNS
//   frag LO HI FORMULA SIZE                     (FORMULA - = none)
//   int RUN ION VALUE
//   end
// Output per case: "gmid NAME FRAG K MID CI" and "gr2 NAME FRAG R2" per fragment key in map order,
// "gerror NAME MESSAGE" for the feature ERROR, then the export between "begin group NAME" and "end".
#include <cxxabi.h>
#include <dlfcn.h>
#include <unistd.h>
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
                    CALC_CTOR = 0x6d931c, RRS_CTOR = 0x639018, GROUP_CTOR = 0x4fbb28, GROUP_RUN = 0x4fa57c,
                    GET_FEATURE = 0x60fa5c, LIST_CTOR = 0x4f40b6, LIST_APPEND = 0x4fee54,
                    WIDGET_CTOR = 0x4e9104, EXPORT = 0x4e9212;
template <class F> F fn(uintptr_t a) { return reinterpret_cast<F>(base + a); }
template <class F> F qt(const char* name) { return reinterpret_cast<F>(dlsym(RTLD_DEFAULT, name)); }

std::string export_path;

struct Case {
    std::string name;
    size_t compound = 0, runs = 0;
    int method = 0;
    double m1 = 0;
    std::list<size_t> gl, gul;
    struct Frag { unsigned long lo, hi, size; std::string formula; };
    std::vector<Frag> frags;
    std::vector<std::tuple<size_t, float, double>> ints;
};

std::list<size_t> runs_of(const std::string& s) {
    std::list<size_t> out;
    std::istringstream w(s);
    for (std::string r; s != "-" && std::getline(w, r, ',');) out.push_back(std::stoul(r));
    return out;
}

// ResultSeries<int,float>: +0x20 map<const Compound*, CompoundSeries> (value right after the key),
// +0x80 sample names; ReplicateResultSeries +0x150 map<group, list<sample name>>. The binary only
// looks series up, so a value larger than CompoundSeries (at least 0x158 bytes) is fine.
struct alignas(8) Series { char b[0x1000]; };

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
    // The MID compounds and the widget keep pointers into the series, so nothing is freed.
    auto* rrs = new (std::align_val_t(16)) char[0x400]();
    fn<void (*)(void*, void*, double, bool, bool)>(RRS_CTOR)(rrs, nullptr, 0.0, false, false);
    auto& names = *reinterpret_cast<std::vector<std::string>*>(rrs + 0x80);
    for (size_t r = 0; r < k.runs; r++) names.push_back("/data/run" + std::to_string(r) + ".cdf");
    auto& groups = *reinterpret_cast<std::map<std::string, std::list<std::string>>*>(rrs + 0x150);
    for (auto& [g, idx] : {std::pair{"GL", k.gl}, std::pair{"GUL", k.gul}})
        for (size_t r : idx) groups[g].push_back(names.at(r));
    auto& series = *reinterpret_cast<std::map<const void*, Series>*>(rrs + 0x20);
    void* cs = &series[c];
    fn<void (*)(void*, const void*, bool)>(CS_CTOR)(cs, c, false);
    fn<void (*)(void*, unsigned long)>(CS_NULL)(cs, k.runs);
    for (auto& [r, ion, v] : k.ints) fn<void (*)(void*, float, unsigned long, double)>(CS_SET)(cs, ion, r, v);
    // As MIDWizard::accept: Calculate(series, method, M1 / 100, tracer, calculator) per substance.
    alignas(8) char tracer[8], group[0x28];
    qt<void (*)(void*, const char*, int)>("_ZN7QString15fromUtf8_helperEPKci")(tracer, "C", 1);
    fn<void (*)(void*, const void*, int, double, const void*, const void*)>(GROUP_CTOR)(
        group, rrs, k.method, k.m1, tracer, calc);
    char* mc = fn<char* (*)(void*, const void*)>(GROUP_RUN)(group, c);
    // MIDCompound+0x290: QList<MIDSolver*>, one per fragment key in map order. MIDSolver: R² at +0,
    // MIDs at +0x8, CIs at +0x20.
    auto* d = *reinterpret_cast<int**>(mc + 0x290);
    auto** solvers = reinterpret_cast<char**>(d + 4) + d[2];
    for (int f = 0; f < d[3] - d[2]; f++) {
        auto& mids = *reinterpret_cast<const std::vector<double>*>(solvers[f] + 0x8);
        auto& cis = *reinterpret_cast<const std::vector<double>*>(solvers[f] + 0x20);
        for (size_t i = 0; i < mids.size(); i++)
            std::printf("gmid\t%s\t%d\t%zu\t%.17g\t%.17g\n", k.name.c_str(), f, i, mids[i],
                        i < cis.size() ? cis[i] : -2.0);
        std::printf("gr2\t%s\t%d\t%.17g\n", k.name.c_str(), f, *reinterpret_cast<double*>(solvers[f]));
    }
    std::string err = fn<std::string (*)(const void*, std::string)>(GET_FEATURE)(mc, "ERROR");
    if (!err.empty()) std::printf("gerror\t%s\t%s\n", k.name.c_str(), err.c_str());
    alignas(8) char list[8];
    fn<void (*)(void*)>(LIST_CTOR)(list);
    fn<void (*)(void*, char* const&)>(LIST_APPEND)(list, mc);
    auto* widget = new (std::align_val_t(16)) char[0x400]();
    fn<void (*)(void*, void*, void*)>(WIDGET_CTOR)(widget, list, nullptr);
    std::remove(export_path.c_str());
    fn<void (*)(void*)>(EXPORT)(widget);
    std::ifstream in(export_path);
    std::stringstream text; text << in.rdbuf();
    std::printf("begin\tgroup\t%s\n%s\nend\n", k.name.c_str(), text.str().c_str());
}

int my_main(int argc, char** argv, char**) {
    if (argc < 4) { std::fprintf(stderr, "usage: %s ISOTOPES.tsv LIBRARY.lbr CASES\n", argv[0]); return 2; }
    export_path = "/tmp/old_midgroups_" + std::to_string(getpid()) + ".csv";
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
            std::string gl, gul;
            w >> k.name >> k.compound >> k.method >> k.m1 >> gl >> gul >> k.runs;
            k.gl = runs_of(gl);
            k.gul = runs_of(gul);
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
            try {
                run(k, argv[2], calc);
            } catch (...) {
                std::printf("error\t%s\t%s\n", k.name.c_str(), abi::__cxa_current_exception_type()->name());
            }
            std::fflush(stdout);
        }
    }
    std::remove(export_path.c_str());
    std::fflush(stdout);
    _exit(0);  // the leaked widgets crash Qt's teardown
}
}  // namespace
// The export asks for a file name (QFileDialog::getSaveFileName, returns a QString); answer with
// export_path instead of a dialog.
extern "C" void* save_file_name(void* ret, void*, const void*, const void*, const void*, void*, int) __asm__(
    "_ZN11QFileDialog15getSaveFileNameEP7QWidgetRK7QStringS4_S4_PS2_6QFlagsINS_6OptionEE");
extern "C" void* save_file_name(void* ret, void*, const void*, const void*, const void*, void*, int) {
    qt<void (*)(void*, const char*, int)>("_ZN7QString15fromUtf8_helperEPKci")(ret, export_path.c_str(),
                                                                              int(export_path.size()));
    return ret;
}
extern "C" int __libc_start_main(int (*main)(int, char**, char**), int argc, char** argv, void (*init)(),
                                 void (*fini)(), void (*rtld_fini)(), void* stack_end) {
    base = reinterpret_cast<uintptr_t>(main) - MAIN;
    using Start = int (*)(int (*)(int, char**, char**), int, char**, void (*)(), void (*)(), void (*)(), void*);
    auto real = reinterpret_cast<Start>(dlsym(RTLD_NEXT, "__libc_start_main"));
    return real(my_main, argc, argv, init, fini, rtld_fini, stack_end);
}
