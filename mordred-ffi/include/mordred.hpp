#ifndef MORDRED_HPP_
#define MORDRED_HPP_

#include <cstddef>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

extern "C" {
#include "mordred.h"
}

namespace mordred {

// Custom deleters following Google C++ style (structs with operator()).
struct CalculatorDeleter {
  void operator()(MordredCalculator* p) const { mordred_calculator_free(p); }
};

struct MoleculeDeleter {
  void operator()(MordredMolecule* p) const { mordred_molecule_free(p); }
};

struct ResultDeleter {
  void operator()(MordredResult* p) const { mordred_result_free(p); }
};

// Forward declarations.
class Molecule;
class Result;

/// RAII wrapper around a mordred descriptor calculator.
class Calculator {
 public:
  Calculator() : handle_(mordred_calculator_new()) {
    if (!handle_) {
      throw std::runtime_error("failed to create mordred calculator");
    }
  }

  /// Return the number of descriptors in the calculator.
  size_t descriptor_count() const {
    return mordred_calculator_descriptor_count(handle_.get());
  }

  /// Return the name of the descriptor at the given index.
  ///
  /// Throws std::out_of_range if the index is invalid.
  std::string descriptor_name(size_t index) const {
    const char* name =
        mordred_calculator_descriptor_name(handle_.get(), index);
    if (name == nullptr) {
      throw std::out_of_range("descriptor index out of range");
    }
    return std::string(name);
  }

  /// Return the names of all descriptors.
  std::vector<std::string> descriptor_names() const {
    std::vector<std::string> names;
    size_t count = descriptor_count();
    names.reserve(count);
    for (size_t i = 0; i < count; ++i) {
      names.push_back(descriptor_name(i));
    }
    return names;
  }

  /// Calculate all descriptors for a molecule and return the results.
  Result calculate(const Molecule& mol) const;

 private:
  std::unique_ptr<MordredCalculator, CalculatorDeleter> handle_;

  // Allow Result's factory to access the raw pointer.
  friend class Result;
};

/// RAII wrapper around a parsed molecule.
class Molecule {
 public:
  /// Parse a SMILES string into a Molecule.
  ///
  /// Throws std::runtime_error on parse failure.
  static Molecule FromSmiles(const std::string& smiles) {
    char* error = nullptr;
    MordredMolecule* raw =
        mordred_molecule_from_smiles(smiles.c_str(), &error);
    if (raw == nullptr) {
      std::string msg = (error != nullptr) ? std::string(error)
                                           : "unknown SMILES parse error";
      if (error != nullptr) {
        mordred_string_free(error);
      }
      throw std::runtime_error(msg);
    }
    return Molecule(raw);
  }

 private:
  explicit Molecule(MordredMolecule* raw) : handle_(raw) {}

  std::unique_ptr<MordredMolecule, MoleculeDeleter> handle_;

  // Calculator needs access to the raw handle.
  friend class Calculator;
};

/// RAII wrapper around a calculation result set.
class Result {
 public:
  /// Return the number of descriptor results.
  size_t size() const { return mordred_result_len(handle_.get()); }

  /// Return the descriptor name at the given index.
  ///
  /// Throws std::out_of_range if the index is invalid.
  std::string name(size_t index) const {
    const char* n = mordred_result_name(handle_.get(), index);
    if (n == nullptr) {
      throw std::out_of_range("result index out of range");
    }
    return std::string(n);
  }

  /// Return whether the descriptor at the given index was calculated
  /// successfully.
  bool is_ok(size_t index) const {
    return mordred_result_is_ok(handle_.get(), index);
  }

  /// Return the descriptor value at the given index, or std::nullopt if the
  /// calculation failed.
  ///
  /// Throws std::out_of_range if the index is invalid.
  std::optional<double> value(size_t index) const {
    if (index >= size()) {
      throw std::out_of_range("result index out of range");
    }
    if (!is_ok(index)) {
      return std::nullopt;
    }
    return mordred_result_value(handle_.get(), index);
  }

 private:
  explicit Result(MordredResult* raw) : handle_(raw) {}

  std::unique_ptr<MordredResult, ResultDeleter> handle_;

  friend class Calculator;
};

// Inline definition (needs Result to be complete).
inline Result Calculator::calculate(const Molecule& mol) const {
  MordredResult* raw = mordred_calculate(handle_.get(), mol.handle_.get());
  if (raw == nullptr) {
    throw std::runtime_error("mordred_calculate returned null");
  }
  return Result(raw);
}

}  // namespace mordred

#endif  // MORDRED_HPP_
