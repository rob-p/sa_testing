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

// build a suffix array over input that is textual (i.e. "DNA" or "Text"), where the input characters
// are encoded as `uint8_t`.
template <typename IdxT>
IdxT build_text_sa(const uint8_t* text, std::vector<IdxT>& sa, size_t len, int32_t nthreads, quill::Logger* logger) {
  IdxT ret = 0;
  // if we are using 32-bit indices
  if constexpr (std::is_same_v<IdxT, std::int32_t>) {
    LOG_INFO(logger, "using 32-bit (int32_t) indices");
    if (nthreads == 1) {
      ret = libsais(text, sa.data(), len, 0, nullptr);
    } else {
      ret = libsais_omp(text, sa.data(), len, 0, nullptr, nthreads);
    }
  } else {
    // if we are using 64-bit indices
    LOG_INFO(logger, "using 64-bit (int32_t) indices");
    if (nthreads == 1) {
      ret = libsais64(text, sa.data(), len, 0, nullptr);
    } else {
      ret = libsais64_omp(text, sa.data(), len, 0, nullptr, nthreads);
    }
  }
  
  // log an error if the return code was non-zero
  if (ret != 0) {
    LOG_ERROR(logger, "error: libsais return code", ret);
  }
  LOG_INFO(logger, "ret :{}", ret);

  return ret;
}

// build a suffix array over an integer alphabet input. In addition to taking the input (which can be 
// an int32_t or int64_t), it also requires the maximum token, which is given as input to this function. 
// This function performs no alphabet remapping, so if you want to remap the alphabet to a minimal/compacted
// space, you should do that prior to calling this function.
template <typename IdxT>
IdxT build_int_sa(IdxT* text, std::vector<IdxT>& sa, size_t len, IdxT max_token, int32_t nthreads, quill::Logger* logger) {
  IdxT ret = 0;
  // the integer alphabet (and text size) both fit in int32_t width
  if constexpr (std::is_same_v<IdxT, std::int32_t>) {
    LOG_INFO(logger, "int alphabet using 32-bit (int32_t) indices");
    if (nthreads == 1) {
      ret = libsais_int(text, sa.data(), len, max_token, 0);
    } else {
      ret = libsais_int_omp(text, sa.data(), len, max_token, 0, nthreads);
    }
  } else {
    // the integer alphabet (and text size) require int64_t width
    LOG_INFO(logger, "int alphabet using 64-bit (int64_t) indices");
    if (nthreads == 1) {
      ret = libsais64_long(text, sa.data(), len, max_token, 0);
    } else {
      ret = libsais64_long_omp(text, sa.data(), len, max_token, 0, nthreads);
    }
  }

  // log and error if the return code was non-zero
  if (ret != 0) {
    LOG_ERROR(logger, "error: libsais return code", ret);
  }
  LOG_INFO(logger, "ret :{}", ret);

  return ret;
}

// write the output suffix array to a file named `output`.  Here the suffix array 
// will be either of int32_t type or int64_t.  No information about the size of the 
// input alphabet (e.g. DNA, Text, Integer) is currently encoded in the output.
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
  quill::Logger* logger = quill::Frontend::create_or_get_logger("SA driver", std::move(console_sink));

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
    
  // decode the choice to a string; fancier solutions are available but we don't need them right now.
  std::string type_str = (in_ty == InputType::DNA) ? "dna" : 
    (in_ty == InputType::Text) ? "text" : "integer";
  LOG_INFO(logger, "input_types: {}", type_str);
  LOG_INFO(logger, "file :{}", filename);

  std::string genome;
  std::vector<int32_t> int32_genome;
  std::vector<int64_t> int64_genome;
  uint64_t max_token{0};

  switch (in_ty) {
    case InputType::DNA:
      {
        // read the input from a FASTA file
        using klibpp::KSeq;
        using klibpp::SeqStreamIn;
        bool first_record = true;
        KSeq record;
        SeqStreamIn iss(filename.c_str());
        while (iss >> record) {
          if (!first_record) {
            LOG_WARNING(logger, "There was more than one record in the FASTA file, but generalized suffix arrays are not yet supported; just taking the first record.");
            break;
          }
          LOG_INFO(logger, "genome size is : {}", record.seq.size());
          genome = record.seq;
          first_record = false;
        }
      }
      break;
    case InputType::Text:
      {
        // read the input as a text string
        std::filesystem::path p{filename};
        size_t fsize = std::filesystem::file_size(p);
        genome.resize(fsize);
        std::ifstream ifile(filename, std::ios::binary);
        ifile.read(reinterpret_cast<char*>(genome.data()), genome.size());
      }
      break;
    case InputType::Integer:
      {
        // read the input as an integer vector
        std::ifstream ifile(filename, std::ios::binary);
        uint64_t len{0};
        ifile.read(reinterpret_cast<char*>(&len), sizeof(len));
        ifile.read(reinterpret_cast<char*>(&max_token), sizeof(max_token));

        // if the length or the max token is too large for 32-bits, then the file is 64-bits
        if (len >= std::numeric_limits<int32_t>::max() or max_token >= std::numeric_limits<int32_t>::max()) {
          int64_genome.resize(len);
          ifile.read(reinterpret_cast<char*>(int64_genome.data()), sizeof(len)*len);
        } else {
          int32_genome.resize(len);
          ifile.read(reinterpret_cast<char*>(int32_genome.data()), sizeof(len)*len);
        }
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
    // 32-bit
    if (int64_genome.empty()) {
      std::vector<int32_t> sa(int32_genome.size(), 0);
      int32_t ret = build_int_sa(int32_genome.data(), sa, int32_genome.size(), static_cast<int32_t>(max_token), static_cast<int32_t>(nthreads), logger);
      (void)ret;
      write_output(output, sa);
    } else {
      std::vector<int64_t> sa(int64_genome.size(), 0);
      int64_t ret = build_int_sa(int64_genome.data(), sa, int64_genome.size(), static_cast<int64_t>(max_token), static_cast<int32_t>(nthreads), logger);
      (void)ret;
      write_output(output, sa);
    }
  }

  return 0;
}
