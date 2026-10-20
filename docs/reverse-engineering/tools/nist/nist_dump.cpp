// Black-box checks of the NIST text import and library search of the old MetaboliteDetector 3.3
// (Linux). Loaded with LD_PRELOAD like old_md_dump.cpp: it replaces main() and calls the binary's own
// NISTParser, SQLLibrary, Spectrum and Compound functions by address (syms.inc, generated from nm by
// make-nist-dump.sh). The input files are made by make-nist-dump.sh; no real library data is used.
//
// usage: LD_PRELOAD=nist_dump.so metabolitedetector_app parse NIST.txt
//        LD_PRELOAD=nist_dump.so metabolitedetector_app library LIB.sqlite NIST.txt QUERY.txt
//
// parse:   NISTParser over every record. Prints "field<TAB>RECORD<TAB>NAME<TAB>VALUE" per field;
//          newlines in VALUE are printed as \n (values are the trimmed lines joined with "\n").
// library: creates SQLLibrary(LIB.sqlite), runs importNIST(NIST.txt), then getPotentialHits and
//          getLibraryHits for the query spectrum in QUERY.txt. QUERY.txt: a line "RT <rt>", then one
//          "<mz> <intensity>" line per peak. Prints "potential<TAB>ID<TAB>NAME" in list order and
//          "hit<TAB>ID<TAB>NAME<TAB>S<TAB>RI<TAB>O" in the order getLibraryHits returns.
//
// Struct layouts used (all Observed in the disassembly, see issue-21-nist-re notes):
//   LibraryHit (0x28 bytes, vector stride from the push_back code): +0x00 Compound*, +0x08 LibraryCompound*,
//   +0x10 spectrum similarity S, +0x18 RI similarity, +0x20 overall score O (sort key).
//   std::list node: {next, prev, value}; the list head object is {next, prev, size}.
#include <dlfcn.h>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <list>
#include <utility>
#include <sstream>
#include <string>
#include <vector>

namespace {

uintptr_t base;

template <class F> F fn(uintptr_t addr) { return reinterpret_cast<F>(base + addr); }

#define SYM(name, addr) constexpr uintptr_t name = addr;
#include "syms.inc"
#undef SYM

// Objects are allocated generously; the sizes are not needed, the class layouts stay inside them.
constexpr size_t kObj = 0x400;

std::string esc(const std::string& s) {
    std::string out;
    for (char c : s) {
        if (c == '\n') out += "\\n";
        else if (c == '\t') out += "\\t";
        else out += c;
    }
    return out;
}

void dump_parser(const char* path) {
    alignas(16) static char parser[kObj];
    fn<void (*)(void*, const char*)>(NISTPARSER_CTOR)(parser, path);
    int record = 0;
    while (fn<bool (*)(void*)>(NEXT_RECORD)(parser)) {
        record++;
        while (fn<bool (*)(void*)>(NEXT_FIELD)(parser)) {
            std::string name, value;
            fn<void (*)(std::string*, void*)>(GET_FIELD_NAME)(&name, parser);
            fn<void (*)(std::string*, void*)>(GET_FIELD_VALUE)(&value, parser);
            std::printf("field\t%d\t%s\t%s\n", record, esc(name).c_str(), esc(value).c_str());
        }
    }
    std::printf("records\t%d\n", record);
    fn<void (*)(void*)>(NISTPARSER_DTOR)(parser);
}

// Query spectrum: "RT x", optional "RI x", then "mz intensity" lines. Masses are floats and intensities ints: that is the
// layout the binary's own importer uses (Spectrum<int,float>: intensity vector<int>, mass vector<float>).
bool read_query(const char* path, double& rt, double& ri, std::vector<float>& masses, std::vector<int>& ints) {
    std::ifstream in(path);
    if (!in) return false;
    std::string line;
    if (!std::getline(in, line)) return false;
    std::istringstream head(line);
    std::string key;
    head >> key >> rt;
    while (std::getline(in, line)) {
        std::istringstream row(line);
        double mz, intensity;
        if (line.compare(0, 3, "RI ") == 0) {
            row >> key >> ri;
        } else if (row >> mz >> intensity) {
            masses.push_back(static_cast<float>(mz));
            ints.push_back(static_cast<int>(intensity));
        }
    }
    return true;
}

void dump_library(const char* db, const char* nist, const char* query) {
    std::remove(db);
    alignas(16) static char lib[kObj];
    fn<void (*)(void*, const char*)>(SQL_CTOR)(lib, db);
    fn<void (*)(void*, const char*)>(SQL_IMPORT_NIST)(lib, nist);
    std::printf("imported\t%s\n", nist);

    double rt = 0, ri = -1;
    std::vector<float> masses;
    std::vector<int> ints;
    if (!read_query(query, rt, ri, masses, ints)) {
        std::fprintf(stderr, "cannot read query %s\n", query);
        return;
    }
    // The query is a Compound (RT, peaks, the last float is passed as 1.0; its meaning is not yet decoded).
    // The Compound constructor modifies the passed vectors in place (sorts ions, drops duplicate masses).
    alignas(16) static char comp[kObj];
    fn<void (*)(void*, double, std::vector<int>*, std::vector<float>*, float)>(COMP_CTOR)(comp, rt, &ints, &masses, 1.0f);
    fn<void (*)(void*, double)>(SPEC_SET_RI)(comp, ri);

    // Fragments of the query: Spectrum::getSpectrumFragments(list<list<pair<float,int>>>& groups, double, int)
    // returns std::list<Fragment> by hidden pointer (the same call importNIST makes). Printed for each int argument.
    using Inner = std::list<std::pair<float, int>>;
    for (int k = 0; k <= 2; k++) {
        alignas(8) static char fragments[24];
        std::list<Inner> groups;
        fn<void (*)(void*, const void*, void*, double, int)>(SPEC_FRAGMENTS)(fragments, comp, &groups, 0.1, k);
        std::printf("fragments\tint=%d\tgroups=%zu\tfragments=%ld\n", k, groups.size(),
                    *reinterpret_cast<long*>(fragments + 16));
    }

    // getPotentialHits returns std::list<LibraryCompound*> by hidden pointer; the double argument is
    // never read by the binary (the caller passes 0.5).
    alignas(8) static char pot[24];
    try {
        fn<void (*)(void*, const void*, const void*, double)>(SQL_POTENTIAL)(pot, lib, comp, 0.5);
    } catch (...) {
        std::printf("exception\tgetPotentialHits\n");
        return;
    }
    void* const head = pot;
    for (void* n = *static_cast<void**>(head); n != head; n = *static_cast<void**>(n)) {
        void* c = static_cast<void**>(n)[2];
        long long id = fn<long long (*)(const void*)>(COMP_ID)(c);
        std::string name;
        fn<void (*)(std::string*, const void*)>(COMP_NAME)(&name, c);
        std::printf("potential\t%lld\t%s\n", id, esc(name).c_str());
    }

    // getLibraryHits(Compound const&, float): the bin width argument 1.0 is not read by the binary.
    alignas(8) static char hits[24];
    try {
        fn<void (*)(void*, const void*, const void*, float)>(SQL_LIB_HITS)(hits, lib, comp, 1.0f);
    } catch (...) {
        std::printf("exception\tgetLibraryHits\n");
        return;
    }
    char* begin = *reinterpret_cast<char**>(hits);
    char* end = *reinterpret_cast<char**>(hits + 8);
    for (char* h = begin; h != end; h += 0x28) {
        void* lc = *reinterpret_cast<void**>(h + 8);
        long long id = fn<long long (*)(const void*)>(COMP_ID)(lc);
        std::string name;
        fn<void (*)(std::string*, const void*)>(COMP_NAME)(&name, lc);
        std::printf("hit\t%lld\t%s\t%.9g\t%.9g\t%.9g\n", id, esc(name).c_str(), *reinterpret_cast<double*>(h + 0x10),
                    *reinterpret_cast<double*>(h + 0x18), *reinterpret_cast<double*>(h + 0x20));
    }
    std::printf("hits\t%zu\n", static_cast<size_t>((end - begin) / 0x28));
    std::fflush(stdout);
    // Objects are leaked on purpose: the process ends after this run.
}

}  // namespace

// usage: see the comment at the top of this file.
namespace {
int dump_main(int argc, char** argv, char**) {
    if (argc >= 3 && std::strcmp(argv[1], "parse") == 0) {
        dump_parser(argv[2]);
        return 0;
    }
    if (argc >= 5 && std::strcmp(argv[1], "library") == 0) {
        dump_library(argv[2], argv[3], argv[4]);
        return 0;
    }
    std::fprintf(stderr, "usage: %s parse NIST.txt | library LIB.sqlite NIST.txt QUERY.txt\n", argv[0]);
    return 2;
}
}  // namespace

extern "C" int __libc_start_main(int (*main)(int, char**, char**), int argc, char** argv, void (*init)(),
                                 void (*fini)(), void (*rtld_fini)(), void* stack_end) {
    base = reinterpret_cast<uintptr_t>(main) - MAIN;
    using Start = int (*)(int (*)(int, char**, char**), int, char**, void (*)(), void (*)(), void (*)(), void*);
    auto real = reinterpret_cast<Start>(dlsym(RTLD_NEXT, "__libc_start_main"));
    return real(dump_main, argc, argv, init, fini, rtld_fini, stack_end);
}
