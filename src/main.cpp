#include "CLI11.hpp"
#include "kseq++/seqio.hpp"
#include "libsais.h"
#include <iostream>
#include <string>
#include <vector>

auto main(int argc, char *argv[]) -> int {
  CLI::App app{"libsais driver"};
  argv = app.ensure_utf8(argv);

  std::string filename = "default";
  app.add_option("-f,--file", filename, "A help string");

  CLI11_PARSE(app, argc, argv);

  std::cout << "file: " << filename << "\n";

  std::string genome;
  using klibpp::KSeq;
  using klibpp::SeqStreamIn;
  KSeq record;
  SeqStreamIn iss(filename.c_str());
  while (iss >> record) {
    std::cout << record.name << std::endl;
    if (!record.comment.empty())
      std::cout << record.comment << std::endl;
    std::cout << record.seq << std::endl;
    if (!record.qual.empty())
      std::cout << record.qual << std::endl;

    genome = record.seq;
  }

  std::vector<int32_t> sa(genome.size(), 0);
  int32_t ret = libsais(reinterpret_cast<const uint8_t *>(genome.c_str()),
                        sa.data(), genome.size(), 0, nullptr);

  if (ret != 0) {
    std::cerr << "error: return code " << ret << "\n";
  }
  std::cout << "ret = " << ret << "\n";
  std::cout << "SA = [";
  for (size_t i = 0; i < sa.size(); ++i) {
    std::cout << sa[i];
    if (i < sa.size() - 1) {
      std::cout << ", ";
    } else {
      std::cout << "]\n";
    }
  }
  return 0;
}

/**
 * Constructs the suffix array of a given string.
 * @param T [0..n-1] The input string.
 * @param SA [0..n-1+fs] The output array of suffixes.
 * @param n The length of the given string.
 * @param fs The extra space available at the end of SA array (0 should be
 * enough for most cases).
 * @param freq [0..255] The output symbol frequency table (can be NULL).
 * @return 0 if no error occurred, -1 or -2 otherwise.
 */
//    LIBSAIS_API int32_t libsais(const uint8_t * T, int32_t * SA, int32_t n,
//    int32_t fs, int32_t * freq);
