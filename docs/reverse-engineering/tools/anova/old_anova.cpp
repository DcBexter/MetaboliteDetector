// LD_PRELOAD into the old binary (sha256 7883d70c…4794): replaces main() and calls
// tools::Statistics::oneWayANOVA (0x6c3af6) on the cases of a file. See make-anova.sh.
// usage: metabolitedetector_app CASES      (cases: "case NAME", one "group V V ..." line per group, "end")
// Output per case: "anova NAME P" (17 significant digits, nan/inf as printf writes them) or "error NAME" where the old code aborts.
#include <dlfcn.h>
#include <sys/wait.h>
#include <unistd.h>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <list>
#include <sstream>
#include <string>
#include <vector>
namespace {
uintptr_t base;
constexpr uintptr_t MAIN = 0x3b1f2d, ONE_WAY_ANOVA = 0x6c3af6;
using Groups = std::list<std::vector<double>>;

int my_main(int argc, char** argv, char**) {
    if (argc < 2) return 2;
    auto anova = reinterpret_cast<double (*)(const Groups&)>(base + ONE_WAY_ANOVA);
    std::ifstream in(argv[1]);
    std::string line, name;
    Groups groups;
    while (std::getline(in, line)) {
        std::istringstream w(line);
        std::string tag;
        w >> tag;
        if (tag == "case") {
            w >> name;
            groups.clear();
        } else if (tag == "group") {
            groups.emplace_back();
            double v;
            while (w >> v) groups.back().push_back(v);
        } else if (tag == "end") {
            // One process per case: the GSL error handler aborts (e.g. one group: domain error).
            std::fflush(stdout);
            pid_t pid = fork();
            if (pid == 0) {
                std::printf("anova\t%s\t%.17g\n", name.c_str(), anova(groups));
                std::fflush(stdout);
                _exit(0);
            }
            int status = 0;
            waitpid(pid, &status, 0);
            if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) std::printf("error\t%s\n", name.c_str());
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
