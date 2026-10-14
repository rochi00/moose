// Kokkos Tools library that counts device allocations, fences and kernel launches, grouped by the
// first MOOSE (non-Kokkos) frames of the call stack that triggered them.
//
// Build: c++ -O2 -std=c++17 -shared -fPIC kokkos_event_counter.cpp -o libkokkos_event_counter.so
// Use:   KOKKOS_TOOLS_LIBS=/path/libkokkos_event_counter.so KEC_OUT=counts.txt <moose app> ...
//
// Per-assembly counts come from two runs with different numbers of assemblies (difference
// method). Copies made with Kokkos::Impl::DeepCopy (as Moose::Kokkos::Array does) do not reach
// the tools interface; those copies always fence, so they show up in the fence counts.

#include <cxxabi.h>
#include <dlfcn.h>
#include <execinfo.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <mutex>
#include <string>
#include <vector>

namespace
{
struct SpaceHandle
{
  char name[64];
};

std::mutex mutex;
std::map<std::string, std::uint64_t> counts;
std::uint64_t next_id = 0;

std::string
demangle(const char * name)
{
  int status = 0;
  char * out = abi::__cxa_demangle(name, nullptr, nullptr, &status);
  std::string result = status == 0 && out ? out : name;
  std::free(out);
  return result;
}

bool
isInfrastructure(const std::string & name)
{
  static const char * const skip[] = {"Kokkos::",
                                      "kokkosp_",
                                      "(anonymous namespace)::",
                                      "Moose::Kokkos::ArrayBase",
                                      "Moose::Kokkos::Array<",
                                      "std::",
                                      "void Kokkos::",
                                      "auto Kokkos::"};
  for (const auto * const prefix : skip)
    if (name.rfind(prefix, 0) == 0)
      return true;
  // Moose::Kokkos::Array internals, whose demangled names may start with a return type
  const auto head = name.substr(0, name.find('('));
  return head.find("Moose::Kokkos::Array") != std::string::npos ||
         head.find("Kokkos::Impl::") != std::string::npos;
}

// Up to three first non-infrastructure frames, outermost last, without template arguments
std::string
origin()
{
  void * frames[48];
  const int n = backtrace(frames, 48);
  std::string result;
  int kept = 0;

  for (int i = 1; i < n && kept < 3; ++i)
  {
    Dl_info info;
    if (!dladdr(frames[i], &info) || !info.dli_sname)
      continue;

    auto name = demangle(info.dli_sname);
    if (isInfrastructure(name))
      continue;

    // Keep the qualified function name only
    const auto paren = name.find('(');
    if (paren != std::string::npos)
      name.resize(paren);
    if (name.size() > 110)
      name = name.substr(0, 110) + "...";

    result += (kept ? " <- " : "") + name;
    ++kept;
  }

  return result.empty() ? "?" : result;
}

void
record(const std::string & kind, const std::string & detail = "")
{
  auto key = kind + " | " + origin();
  if (!detail.empty())
    key += " | " + detail;

  std::lock_guard<std::mutex> lock(mutex);
  ++counts[key];
}
}

extern "C"
{
struct Kokkos_Tools_ToolSettings
{
  bool requires_global_fencing;
  bool padding[255];
};

// Without this, Kokkos fences around every kernel while a tool is loaded, which would be counted
void
kokkosp_request_tool_settings(const std::uint32_t, Kokkos_Tools_ToolSettings * settings)
{
  settings->requires_global_fencing = false;
}

void
kokkosp_init_library(int, std::uint64_t, std::uint32_t, void *)
{
}

void
kokkosp_finalize_library()
{
  const char * const path = std::getenv("KEC_OUT");
  FILE * out = path ? std::fopen(path, "w") : stderr;

  std::map<std::string, std::uint64_t> totals;
  for (const auto & [key, count] : counts)
    totals[key.substr(0, key.find(" |"))] += count;

  for (const auto & [kind, count] : totals)
    std::fprintf(out, "TOTAL %s %llu\n", kind.c_str(), static_cast<unsigned long long>(count));

  std::vector<std::pair<std::uint64_t, std::string>> sorted;
  for (const auto & [key, count] : counts)
    sorted.emplace_back(count, key);
  std::sort(sorted.rbegin(), sorted.rend());

  for (const auto & [count, key] : sorted)
    std::fprintf(out, "%8llu  %s\n", static_cast<unsigned long long>(count), key.c_str());

  if (path)
    std::fclose(out);
}

void
kokkosp_begin_parallel_for(const char *, std::uint32_t, std::uint64_t * id)
{
  record("launch");
  *id = next_id++;
}

void
kokkosp_end_parallel_for(std::uint64_t)
{
}

void
kokkosp_begin_parallel_reduce(const char *, std::uint32_t, std::uint64_t * id)
{
  record("launch");
  *id = next_id++;
}

void
kokkosp_end_parallel_reduce(std::uint64_t)
{
}

void
kokkosp_begin_fence(const char * name, std::uint32_t, std::uint64_t * id)
{
  record("fence", name);
  *id = next_id++;
}

void
kokkosp_end_fence(std::uint64_t)
{
}

void
kokkosp_allocate_data(SpaceHandle, const char *, const void *, std::uint64_t)
{
  record("allocate");
}

void
kokkosp_deallocate_data(SpaceHandle, const char *, const void *, std::uint64_t)
{
  record("free");
}

void
kokkosp_begin_deep_copy(SpaceHandle dst, const char *, const void *, SpaceHandle src, const char *,
                        const void *, std::uint64_t)
{
  record("deep_copy", std::string(src.name) + "->" + dst.name);
}

void
kokkosp_end_deep_copy()
{
}
}
