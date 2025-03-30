#include "CLI11.hpp"
#include "kseq++/seqio.hpp"
#include "libsais.h"
#include "libsais64.h"

// logger
#include "quill/Backend.h"
#include "quill/Frontend.h"
#include "quill/LogMacros.h"
#include "quill/Logger.h"
#include "quill/sinks/ConsoleSink.h"
#include "quill/std/Array.h"

// std
#include <utility>
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <sstream>
#include <cstdint>
#include <limits>
#include <filesystem>

enum class InputType : uint8_t { DNA, Text, Integer };


template <typename IdxT>
IdxT build_text_sa(const uint8_t* text, std::vector<IdxT>& sa, size_t len, int32_t nthreads, quill::Logger* logger) {
  IdxT ret = 0;
  if constexpr (std::is_same_v<IdxT, std::int32_t>) {
    LOG_INFO(logger, "using 32-bit (int32_t) indices");
    if (nthreads == 1) {
      ret = libsais(text, sa.data(), len, 0, nullptr);
    } else {
      ret = libsais_omp(text, sa.data(), len, 0, nullptr, nthreads);
    }
  } else {
    LOG_INFO(logger, "using 64-bit (int32_t) indices");
    if (nthreads == 1) {
      ret = libsais64(text, sa.data(), len, 0, nullptr);
    } else {
      ret = libsais64_omp(text, sa.data(), len, 0, nullptr, nthreads);
    }
  }

  if (ret != 0) {
    LOG_ERROR(logger, "error: libsais return code", ret);
  }
  LOG_INFO(logger, "ret :{}", ret);

  return ret;
}

template <typename IdxT>
IdxT build_int_sa(IdxT* text, std::vector<IdxT>& sa, size_t len, int32_t nthreads, quill::Logger* logger) {
  IdxT ret = 0;
  if constexpr (std::is_same_v<IdxT, std::int32_t>) {
    LOG_INFO(logger, "int alphabet using 32-bit (int32_t) indices");
    int64_t m = 0;
    for (size_t i = 0; i < len; ++i) {
      if (text[i] > m) { m = text[i]; }
    }
    if (nthreads == 1) {
      ret = libsais64_long(text, sa.data(), len, m, 0);
    } else {
      ret = libsais64_long_omp(text, sa.data(), len, m, 0, nthreads);
    }
  } else {
    LOG_INFO(logger, "int alphabet using 64-bit (int64_t) indices");
    int64_t m = 0;
    for (size_t i = 0; i < len; ++i) {
      if (text[i] > m) { m = text[i]; }
    }
    if (nthreads == 1) {
      ret = libsais64_long(text, sa.data(), len, m, 0);
    } else {
      ret = libsais64_long_omp(text, sa.data(), len, m, 0, nthreads);
    }
    /*
    LOG_INFO(logger, "using 64-bit (int32_t) indices");
    if (nthreads == 1) {
      ret = libsais64(text, sa.data(), len, 0, nullptr);
    } else {
      ret = libsais64_omp(text, sa.data(), len, 0, nullptr, nthreads);
    }
    */
  }

  if (ret != 0) {
    LOG_ERROR(logger, "error: libsais return code", ret);
  }
  LOG_INFO(logger, "ret :{}", ret);

  return ret;
}

template <typename IdxT>
int write_output(const std::string& output, const std::vector<IdxT>& sa) {
  std::ofstream out(output, std::ios::binary);
  uint64_t nelem = static_cast<uint64_t>(sa.size());
  uint8_t elem_size = sizeof(decltype(sa[0]));
  out.write(reinterpret_cast<char*>(&nelem), sizeof(nelem));
  out.write(reinterpret_cast<char*>(&elem_size), sizeof(elem_size));
  out.write(reinterpret_cast<const char*>(sa.data()), elem_size * nelem);
  out.close();
  return 0;
}

auto main(int argc, char *argv[]) -> int {
  quill::Backend::start();
  // Frontend
  auto console_sink = quill::Frontend::create_or_get_sink<quill::ConsoleSink>("sink_id_1");
  quill::Logger* logger = quill::Frontend::create_or_get_logger("root", std::move(console_sink));


  CLI::App app{"libsais driver"};
  argv = app.ensure_utf8(argv);

  std::string filename = "default";
  std::string output = "default";
  size_t nthreads = 4;
  app.add_option("-f,--file", filename, "input filename");
  app.add_option("-o,--output", output, "output filename")->required();

  // specify string->value mappings
  InputType in_ty{InputType::DNA};
  std::map<std::string, InputType> map{{"dna", InputType::DNA}, {"text", InputType::Text}, {"integer", InputType::Integer}};
  app.add_option("--input-type", in_ty, "input type")->required()->transform(CLI::CheckedTransformer(map, CLI::ignore_case));


  app.add_option("-t,--threads", nthreads, "number of threads");
  CLI11_PARSE(app, argc, argv);
    
  std::string type_str = (in_ty == InputType::DNA) ? "dna" : 
    (in_ty == InputType::Text) ? "text" : "integer";
  LOG_INFO(logger, "input_types: {}", type_str);
  LOG_INFO(logger, "file :{}", filename);

  std::string genome;
  std::vector<int64_t> int_genome;

  switch (in_ty) {
    case InputType::DNA:
      {
        using klibpp::KSeq;
        using klibpp::SeqStreamIn;
        KSeq record;
        SeqStreamIn iss(filename.c_str());
        while (iss >> record) {
          LOG_INFO(logger, "genome size is : {}", record.seq.size());
          genome = record.seq;
        }
      }
      break;
    case InputType::Text:
      {
        std::filesystem::path p{filename};
        size_t fsize = std::filesystem::file_size(p);
        genome.resize(fsize);
        std::ifstream ifile(filename, std::ios::binary);
        ifile.read(reinterpret_cast<char*>(genome.data()), genome.size());
      }
      break;
    case InputType::Integer:
      {
        std::ifstream ifile(filename, std::ios::binary);
        uint64_t len{0};
        ifile.read(reinterpret_cast<char*>(&len), sizeof(len));
        int_genome.resize(len);
        ifile.read(reinterpret_cast<char*>(int_genome.data()), sizeof(len)*len);
      }
  }

  size_t input_len = genome.size();

  // this can certainly be made cleaner, but I'm just getting it working now; this dispatches
  // based on the size of the input.
  if (in_ty == InputType::DNA or in_ty == InputType::Text) {
    if (input_len >= std::numeric_limits<int32_t>::max()) {
      std::vector<int64_t> sa(genome.size(), 0);
      int64_t ret = build_text_sa(reinterpret_cast<const uint8_t*>(genome.c_str()), sa, genome.size(), static_cast<int32_t>(nthreads), logger);
      (void)ret;
      write_output(output, sa);
    } else {
      std::vector<int32_t> sa(genome.size(), 0);
      int32_t ret = build_text_sa(reinterpret_cast<const uint8_t*>(genome.c_str()), sa, genome.size(), static_cast<int32_t>(nthreads), logger);
      (void)ret;
      write_output(output, sa);
    }
  } else {
      std::vector<int64_t> sa(int_genome.size(), 0);
      int64_t ret = build_int_sa(int_genome.data(), sa, int_genome.size(), static_cast<int32_t>(nthreads), logger);
      (void)ret;
      write_output(output, sa);
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
