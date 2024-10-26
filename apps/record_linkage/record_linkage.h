#ifndef PARTY_RL
#define PARTY_RL

#include "../../core/Party.h"
#include "../../utils/connection.h"
#include "../../booleancore/core.h"
#include <memory>
#include <numeric>
#include "../../core/core.h"
#include <cryptopp/blake2.h>
#include <algorithm>

// uncomment to disable assert()
// #define NDEBUG
#include <cassert>
#include <utility>
// Use (void) to silence unused warnings.
#define assertm(exp, msg) assert(((void)msg, exp))

#define ARRAY_SIZE 900
#define BIT_PRECISION 8 // this gives precision to roughly the 6th decimal point

/**
* @brief A boolean array with its associated size and hamming weight (number of true values).
*
*/
struct BooleanArray{
  double hamming_weight;
  size_t size;
  std::unique_ptr<bool[]> values;

  /**
  * @brief Constructs a boolean array.
  *
  * @param size p_size: The size of the array.
  */
  explicit BooleanArray(size_t size = ARRAY_SIZE) {
    hamming_weight = 0;
    this->size = size;
    values = std::make_unique<bool[]>(size); // this does zero-initialisation
  }
};

/**
* @brief A secret shared IDAT field that allows partial matches.
*
* String-encoded IDAT fields are converted to a boolean array of bigrams, similar to a Bloom filter.
* Since we don't want to leak whether a field exists in a record, an empty boolean array is passed for non-existing fields, which is why there is a value {0, 1} to denote whether the field exists. For fields that combine multiple record attributes, this is a fractional value (e.g. 0.666 if two out of three attributes have a value).
* Lastly, the Hamming Weight (i.e. the number of true values in the array) is also stored to be used for Dice coefficient computation.
*/
struct FuzzyField{
  uint64_t has_value{};
  uint64_t hamming_weight{};
  std::unique_ptr<uint8_t[]> boolean_array;

  /**
  * @brief Constructs a FuzzyField.
  *
  * If a helper party constructs a FuzzyField, only a nullptr is stored instead of an array.
  * @param proxy_role p_proxy_role: The Role of the Party which creates the FuzzyField.
  * @param size p_size: The number of elements stored in the array.
  */
  explicit FuzzyField(size_t boolean_array_byte_size, bool initialise_array = true) {
    if (initialise_array) {
      boolean_array = std::make_unique<uint8_t[]>(boolean_array_byte_size); // each bit of the uint8_t is one value
    } else {
      boolean_array = nullptr;
    }
  }
};

/**
* @brief An array of secret shared IDAT fields allowing partial matches.
*
* To make operations with CECILIA easier,  each component of the FuzzyField gets its own array.
* String-encoded IDAT fields are converted to a boolean array of bigrams, similar to a Bloom filter.
* Since we don't want to leak whether a field exists in a record, an empty boolean array is passed for non-existing fields, which is why there is a value {0, 1} to denote whether the field exists.
* Lastly, the Hamming Weight (i.e. the number of true values in the array) is also stored to be used for Dice coefficient computation.
*/
struct FuzzyFields{
  size_t size;
  std::unique_ptr<uint64_t[]> has_value;
  std::unique_ptr<uint64_t[]> hamming_weights;
  std::unique_ptr<uint8_t[]> boolean_arrays;

  /**
  * @brief Constructs FuzzyFields.
  *
  * If a helper party constructs FuzzyFields, only nullptrs are stored instead of arrays.
  * @param proxy_role p_proxy_role: The Role of the Party which creates the FuzzyFields.
  * @param size p_size: The number of FuzzyFields to store.
  * @param boolean_array_size p_boolean_array_size: The number of elements in each boolean array.
  */
  explicit FuzzyFields(size_t size, size_t byte_size_per_boolean_array, bool initialise_arrays = true) {
    this->size = size;
    if (initialise_arrays) {
      has_value = std::make_unique<uint64_t[]>(size);
      hamming_weights = std::make_unique<uint64_t[]>(size);
      boolean_arrays = std::make_unique<uint8_t[]>(size*byte_size_per_boolean_array);
    } else {
      has_value = nullptr;
      hamming_weights = nullptr;
      boolean_arrays = nullptr;
    }
  }
};

/**
 * @brief A simple pair of secret shared values where one indicates whether the value exists.
 *
 * To avoid information leaking we need to store values even if they do not exist. has_value is either 0 (the value does not exist) or 1 (the value exists).
 */
struct OptionalValue{
  uint64_t has_value;
  uint64_t value;
};

/**
* @brief An array of OptionalValues.
*
* To make operations with CECILIA easier, each component gets its own array.
*/
struct OptionalValues{
  size_t size;
  std::unique_ptr<uint64_t[]> has_value;
  std::unique_ptr<uint64_t[]> values;

  /**
  * @brief Constructs OptionalValues.
  *
  * If a helper party constructs FuzzyFields, only nullptrs are stored instead of arrays.
  * @param proxy_role p_proxy_role: The Role of the party creating the OptionalValues.
  * @param size p_size: The number of OptionalValues to store.
  */
  explicit OptionalValues(size_t size, bool initialise_arrays = true) {
    this->size = size;
    if (initialise_arrays) {
      has_value = std::make_unique<uint64_t[]>(size);
      values = std::make_unique<uint64_t[]>(size);
    } else {
      has_value = nullptr;
      values = nullptr;
    }
  }

  OptionalValues(OptionalValues &&values) noexcept
  : size(values.size)
  , has_value(std::move(values.has_value))
  , values(std::move(values.values)) {}

  OptionalValues& operator= (OptionalValues &&) = default;
};

/**
 * @brief The secret shared IDAT fields of a record. Consists of fields allowing inexact matches and fields only allowing exact matches.
 *
 */
struct Record{
  FuzzyFields fuzzy_fields;
  OptionalValues exact_fields;
  /**
  * @brief Constructs a Record.
  *
  * If a helper party constructs a Record, only nullptrs are stored instead of arrays.
  * @param proxy_role p_proxy_role: The Role of the party creating the Record.
  * @param fuzzy_field_size p_fuzzy_field_size: The number of FuzzyFields of the record.
  * @param exact_field_size p_exact_field_size: The number of exact fields of the record.
  */
  Record(size_t fuzzy_field_size, size_t exact_field_size, size_t byte_size_per_boolean_array, bool initialise_array = true)
  : fuzzy_fields(fuzzy_field_size, byte_size_per_boolean_array, initialise_array)
  , exact_fields(exact_field_size, initialise_array)
  {}
};

/**
* @brief An array of secret shared records. Each field gets its own array to make working with CECILIA easier.
*
*/
struct Records{
  size_t size;
  FuzzyFields fuzzy_fields;
  OptionalValues exact_fields;

  /**
  * @brief Constructs Records.
  *
  * If a helper party constructs Records, only nullptrs are stored instead of arrays.
  * @param proxy_role p_proxy_role: The Role of the party creating the Records.
  * @param fuzzy_field_size p_fuzzy_field_size: The number of FuzzyFields per record.
  * @param exact_field_size p_exact_field_size: The number of exact fields per record.
  * @param size p_size: The number of records.
  */
  Records(
    size_t fuzzy_fields_per_record,
    size_t exact_fields_per_record,
    size_t byte_size_per_boolean_array,
    size_t size,
    bool initialise_arrays = true
  ) : fuzzy_fields(fuzzy_fields_per_record*size, byte_size_per_boolean_array,initialise_arrays)
    , exact_fields(exact_fields_per_record*size, initialise_arrays)
  {
      this->size = size;
  }
};

/**
* @brief A secret shared match score between two records.
*
* To avoid having to compute too many expensive divisions, the numerator and denominator are stored separately.
* The only time the actual score has to be computed is to determine whether the match threshold was passed.
*/
struct Score{
  uint64_t numerator;
  uint64_t denominator;
};

/**
* @brief An array of secret shared match scores.
*
*/
struct Scores{
  size_t size;
  std::unique_ptr<uint64_t[]> numerators;
  std::unique_ptr<uint64_t[]> denominators;

  /**
  * @brief Constructs Scores.
  *
  * If a helper party constructs Scores, only nullptrs are stored instead of arrays.
  * @param proxy_role p_proxy_role: The Role of the party creating the Scores.
  * @param size p_size: The number of scores.
  */
  explicit Scores(size_t size, bool initialise_array = true, bool initialise_to_zero = false) {
    this->size = size;
    if (initialise_array) {
      if (initialise_to_zero) {
        numerators = std::make_unique<uint64_t[]>(size);
        denominators = std::make_unique<uint64_t[]>(size);
      } else {
        numerators = std::make_unique<uint64_t[]>(size);
        denominators = std::make_unique<uint64_t[]>(size);
      }
    } else {
      numerators = nullptr;
      denominators = nullptr;
    }
  }
};

/**
 * @brief The secret shared maximum score and its (secret shared) index.
 *
 */
struct MaxScore : Score {
  uint64_t index;
};

/**
* @brief An array of MaxScores.
*
*/
struct MaxScores : Scores {
  std::unique_ptr<uint64_t[]> indices;

  /**
  * @brief Constructs MaxScores.
  *
  * If a helper party constructs MaxScores, only nullptrs are stored instead of arrays.
  * @param proxy_role p_proxy_role: The Role of the party creating the MaxScores.
  * @param size p_size: The number of MaxScores.
  */
  explicit MaxScores(size_t size, bool initialise_arrays = true) : Scores(size, initialise_arrays) {
    if (initialise_arrays) {
      indices = std::make_unique<uint64_t[]>(size);
    } else {
      indices = nullptr;
    }
  }
};

/**
 * @brief An extension of Party implementing the required methods to compute record linkage.
 *
 * Since most of the functions required for record linkage are very specific, they are simply member methods.
 */
class RecordLinkage {
  friend class RecordLinkageTest; // allows RecordLinkageTest to access private methods of PartyRL for testing
public:
  /**
  *
  * @param Role p_role: The Role of the party (proxy1, proxy2 or helper)
  * @param helper_port p_helper_port: The port on which the helper is listening
  * @param helper_ip p_helper_ip: The IP address of the helper
  * @param p1_port p_p1_port: The port on which proxy1 is listening
  * @param p1_ip p_p1_ip: The IP address of proxy1
  * @param exact_field_count p_exact_field_count: The number of ExactFields each record has that are not in exchange groups
  * @param fuzzy_field_count p_fuzzy_field_count: The number of FuzzyFields each record has
  * @param exact_field_weights p_exact_field_weights: The weights of each ExactField (minus those of the exchange groups) for computing the similarity score
  * @param fuzzy_field_weights p_fuzzy_field_weights: The weights of each FuzzyField for computing the similarity score
  * @param boolean_array_size p_boolean_array_size: The number of items in each boolean array of the fuzzy fields
  * @param match_threshold p_match_threshold: The score threshold at which similarity two records are considered matches
  * @param bit_precision p_bit_precision: How many bits of the uint64_t should be preserved for floating points
  */
  RecordLinkage(
    Role role,
    uint16_t helper_port,
    const string &helper_ip,
    uint16_t p1_port,
    const string &p1_ip,
    std::vector<double> &exact_field_weights,
    std::vector<double> &fuzzy_field_weights,
    double match_threshold,
    size_t boolean_array_size = ARRAY_SIZE,
    int bit_precision = BIT_PRECISION
  ) : hasher((unsigned int) 8) {
    this->exact_fields_per_record_ = exact_field_weights.size();
    this->fuzzy_fields_per_record_ = fuzzy_field_weights.size();
    this->boolean_array_size_ = boolean_array_size;
    this->bit_precision_ = bit_precision;
    this->boolean_array_byte_size_ = boolean_array_size/8;
    if (boolean_array_size % 8 != 0) {
      this->boolean_array_byte_size_++;
    }
    this->proxy = std::make_shared<Party>(role, helper_port, helper_ip, p1_port, p1_ip, bit_precision_);
    this->exact_field_weights_ = std::vector<uint64_t>();
    exact_field_weights_.reserve(exact_fields_per_record_);
    this->fuzzy_field_weights_ = std::vector<uint64_t>();
    fuzzy_field_weights_.reserve(fuzzy_fields_per_record_);
    for (size_t i = 0; i < exact_fields_per_record_; i++) {
      this->exact_field_weights_[i] = ConvertToUint64(exact_field_weights[i], bit_precision);
    }
    for (size_t i = 0; i < fuzzy_fields_per_record_; i++) {
      this->fuzzy_field_weights_[i] = ConvertToUint64(fuzzy_field_weights[i], bit_precision);
    }
    this->match_threshold_ = ConvertToUint64(match_threshold, BIT_PRECISION);
  }

  /**
   *
   * @param proxy p_proxy: An already initialised Party object
   * @param exact_field_count p_exact_field_count: The number of ExactFields each record has that are not in exchange groups
   * @param fuzzy_field_count p_fuzzy_field_count: The number of FuzzyFields each record has
   * @param exact_field_weights p_exact_field_weights: The weights of each ExactField (minus those of the exchange groups) for computing the similarity score
   * @param fuzzy_field_weights p_fuzzy_field_weights: The weights of each FuzzyField for computing the similarity score
   * @param boolean_array_size p_boolean_array_size: The number of items in each boolean array of the fuzzy fields
   * @param match_threshold p_match_threshold: The score threshold at which similarity two records are considered matches
   * @param bit_precision p_bit_precision: How many bits of the uint64_t should be preserved for floating points
   */
  RecordLinkage(
  const std::shared_ptr<Party>& proxy,
  size_t exact_field_count,
  size_t fuzzy_field_count,
  std::vector<double> &exact_field_weights,
  std::vector<double> &fuzzy_field_weights,
  double match_threshold,
  size_t boolean_array_size = ARRAY_SIZE,
  int bit_precision = BIT_PRECISION
  ) {
    this->exact_fields_per_record_ = exact_field_count;
    this->fuzzy_fields_per_record_ = fuzzy_field_count;
    this->boolean_array_size_ = boolean_array_size;
    this->bit_precision_ = bit_precision;
    this->boolean_array_byte_size_ = boolean_array_size/8;
    if (boolean_array_size % 8 != 0) {
      this->boolean_array_byte_size_++;
    }
    this->proxy = proxy;
    this->exact_field_weights_ = std::vector<uint64_t>();
    exact_field_weights_.reserve(exact_fields_per_record_);
    this->fuzzy_field_weights_ = std::vector<uint64_t>();
    fuzzy_field_weights_.reserve(fuzzy_fields_per_record_);
    for (size_t i = 0; i < exact_fields_per_record_; i++) {
      this->exact_field_weights_[i] = ConvertToUint64(exact_field_weights[i], bit_precision);
    }
    for (size_t i = 0; i < fuzzy_fields_per_record_; i++) {
      this->fuzzy_field_weights_[i] = ConvertToUint64(fuzzy_field_weights[i], bit_precision);
    }
    this->match_threshold_ = ConvertToUint64(match_threshold, BIT_PRECISION);
  }

  Role GetPRole() {
    return proxy->GetPRole();
  }

  /**
  * @brief creates a share with the set bit precision of the Party.
  *
  * Since the original CreateShare of Party does not consider bit precision,  this shadows that method.
  * @param value p_value: The value to create a share for
  * @return uint64_t The secret share of the value
  */
  uint64_t CreateShare(double value) {
    uint64_t share;
    if (GetPRole() == proxy1) {
      share = proxy->GenerateCommonRandom();
    }
    else if (GetPRole() == proxy2) {
      share = ConvertToUint64(value, bit_precision_) - proxy->GenerateCommonRandom();
    }
    return share;
  }

  /**
  * @brief creates shares with the set bit precision of the Party.
  *
  * Since the original CreateShare of Party does not consider bit precision,  this shadows that method.
  * @param value p_value: The values to create shares for
  * @return uint64_t* The secret shares of the value
  */
  uint64_t* CreateShare(double *value, uint32_t size){ // retains the naming to shadow the original function of Party
    uint64_t *converted = ConvertToUint64(value,size, bit_precision_);
    auto *share = new uint64_t[size];
    for (uint32_t i=0;i<size;i++){
      if (proxy->GetPRole() == proxy1) {
        share[i] = proxy->GenerateCommonRandom();
      }
      else{
        share[i] = converted[i] - proxy->GenerateCommonRandom();
      }
    }
    delete[] converted;
    return share;
  }

  /**
   * @brief sends the number of records that have been read in to the other proxy and the helper.
   * @param value the number of records
   */
  void SendSize(size_t value) {
    if (proxy->GetPRole() != helper) {
      int* socket;
      if (proxy->GetPRole() == proxy1) {
        socket = proxy->GetSocketP2();
      } else { // proxy2
        socket = proxy->GetSocketP1();
      }
      std::memcpy(proxy->GetBuffer1(), &value, 8);
      Send(socket, proxy->GetBuffer1(), 8);
      Send(proxy->GetSocketHelper(), proxy->GetBuffer1(), 8);
    }
  }

  /**
   * @brief receives the number of records that have been read in by the other proxy.
   * @param sender the role of the other proxy. Must be either proxy1 or proxy2.
   * @return
   */
  size_t ReceiveSize(Role sender) {
    int* socket;
    if (sender == proxy1) {
      socket = proxy->GetSocketP1();
    } else if (sender == proxy2) {
      socket = proxy->GetSocketP2();
    } else { // helper
      std::cerr << "Tried to receive size from helper" << std::endl;
      exit(EXIT_FAILURE);
    }
    Receive(socket, proxy->GetBuffer1(), 8);
    return ((size_t*) proxy->GetBuffer1())[0];
  }

  uint64_t ShareValue(double value) {
    uint64_t converted = ConvertToUint64(value, bit_precision_);
    return converted - proxy->GenerateCommonRandom();
  }

  uint64_t ReceiveValue() {
    return proxy->GenerateCommonRandom();
  }

  [[nodiscard]] size_t GetBooleanArrayByteSize() const {
    return boolean_array_byte_size_;
  }

  [[nodiscard]] size_t GetBooleanArraySize() const {
    return boolean_array_size_;
  }

  [[nodiscard]] size_t GetFuzzyFieldsPerRecord() const {
    return fuzzy_fields_per_record_;
  }

  [[nodiscard]] size_t GetExactFieldsPerRecord() const {
    return exact_fields_per_record_;
  }

  /**
  * @brief Reconstructs a secret share and converts it to double.
  *
  * @param value p_value: The secret share
  * @return double The reconstructed value
  */
  double ReconstructDouble(uint64_t value) {
    uint64_t reconstructed = Reconstruct(proxy.get(), value);
    return ConvertToDouble(reconstructed, bit_precision_);
  }

  /**
  * @brief Reconstructs secret shares and converts them to double.
  *
  * @param values p_values: The secret shares
  * @param size p_size: the number of shares in values
  * @return double The reconstructed values
  */
  double* ReconstructDouble(const uint64_t *const values, size_t size) {
    uint64_t *reconstructed = Reconstruct(proxy.get(), values, size);
    double *converted = ConvertToDouble(reconstructed, size, bit_precision_);
    delete[] reconstructed;
    return converted;
  }

  /**
    * @brief Returns a unique integer ranging from 0 to 899 for each possible bigram of the alphabet considered.
    *
    * We only consider the letters a-z, "-", " " and "'". Every other symbol is considered identical.
    * Unicode characters should be converted to ASCII in advance using unidecode.
    * @param first_char p_first_char: The first character of the bigram
    * @param second_char p_second_char: The second character of the bigram
    * @return int the mapping of the bigram
    */
  static int MapBigramToInt(char first_char,  char second_char) {
    return MapCharToInt(first_char) * 30 + MapCharToInt(second_char);
  }

  /**
   * @brief computes a uint64_t hash for a given string.
   * This is necessary because even for fields that contain integer values, human error or data corruption mean that they have to be represented as strings.
   * @param value the string to convert the hash for
   * @return the hash
   */
  uint64_t Hash(std::string &value) {
    hasher.Update((const CryptoPP::byte*)value.data(), value.size());
    uint64_t hash;
    hasher.Final((CryptoPP::byte*) &hash);
    return hash;
  }

  /**
  * @brief Computes a boolean array similar to a Bloom filter from the given string.
  *
  * A uint8_t array of size ARRAY_SIZE is zero-initialised.
  * Then, a predefined integer i is computed for each bigram contained in the string and the value at the ith index of the array is set to 1.
  * The bigrams also include "space with the first character" and "the last character with space".
  * @param string p_string: The string to be converted
  * @return BooleanArray The boolean array including hamming weight.
  */
  static BooleanArray ConvertStringToBooleanArray(std::string string) {
    BooleanArray array;
    int mapped_bigram;
    array.hamming_weight = 1;
    array.values[MapBigramToInt(' ', string[0])] = true;
    for (size_t i = 1; i<string.length(); i++) {
      mapped_bigram = MapBigramToInt(string[i-1],  string[i]);
      //if (not array.values[mapped_bigram]) {
        array.hamming_weight += 1;
        array.values[mapped_bigram] = true;
     // }
    }
    mapped_bigram = MapBigramToInt(string[string.length()-1],  ' ');
    //if (not array.values[mapped_bigram]) {
        array.hamming_weight += 1;
        array.values[mapped_bigram] = true;
    //  }
    return array;
  }

  /**
   * @brief Share a fuzzy field with the other proxy.
   *
   * @param array p_array: the boolean array of the fuzzy field
   * @param has_value p_has_value: A value between 0 and 1 denoting what ratio of contained attributes actually exists (e.g. 0.666 if 2 out of 3 attributes exist)
   * @return FuzzyField
   */
  FuzzyField ShareFuzzyField(const BooleanArray &array, double has_value = 1.0) {
    bool is_not_helper = proxy->GetPRole() != helper;
    FuzzyField field(boolean_array_byte_size_, is_not_helper);
    if (is_not_helper) {
      field.has_value = ConvertToUint64(has_value, bit_precision_) - proxy->GenerateCommonRandom();
      if (has_value > 0) {
        field.hamming_weight = ConvertToUint64(array.hamming_weight, bit_precision_) - proxy->GenerateCommonRandom();
        std::unique_ptr<uint8_t[]> array_values(ConvertToUint8(array.values.get(), array.size));
        for (size_t i = 0; i < boolean_array_byte_size_; i++) {
          field.boolean_array[i] = array_values[i] ^ proxy->GenerateCommonRandomByte();
        }
      } else { // we need to make sure that the same number of values has been randomly generated
        proxy->GenerateCommonRandom();
        for (size_t i = 0; i < boolean_array_byte_size_; i++) {
          proxy->GenerateCommonRandomByte();
        }
      }
    }
    return field;
  }

  /**
   * @brief Share a fuzzy field with the other proxy and store the shared values in an existing FuzzyFields object.
   *
   * @param fields p_fields: The FuzzyFields object in which to store the shared values.
   * @param index p_index: The index to store the shared values at.
   * @param array p_array: the boolean array of the fuzzy field
   * @param has_value p_has_value: A value between 0 and 1 denoting what ratio of contained attributes actually exists (e.g. 0.666 if 2 out of 3 attributes exist)
   */
  void ShareFuzzyField(FuzzyFields &fields, size_t index, const BooleanArray &array, double has_value = 1.0) {
    if (proxy->GetPRole() != helper) {
      assertm(array.size == boolean_array_size_, "The number of items in the boolean array does not match FuzzyFields");
      assertm(index < fields.size, "Index out of bound for fields");
      FuzzyField field = ShareFuzzyField(array, has_value);
      fields.has_value[index] = field.has_value;
      if (has_value > 0) {
        fields.hamming_weights[index] = field.hamming_weight;
        std::copy(
          field.boolean_array.get(),
          field.boolean_array.get() + boolean_array_byte_size_,
          fields.boolean_arrays.get() + GetArrayIndex(index)
          );
      }
    }
  }

  /**
   * @brief Receive a fuzzy field from the other proxy.
   *
   * @return FuzzyField
   */
  FuzzyField ReceiveFuzzyField() {
    bool is_not_helper = proxy->GetPRole() != helper;
    FuzzyField field(boolean_array_byte_size_, is_not_helper);
    if (is_not_helper) {
      field.has_value = proxy->GenerateCommonRandom();
      field.hamming_weight = proxy->GenerateCommonRandom();
      for (size_t i = 0; i < boolean_array_byte_size_; i++) {
        field.boolean_array[i] = proxy->GenerateCommonRandomByte();
      }
    }
    return field;
  }

  /**
   * @brief Receive a fuzzy field from the other proxy and store the shared values in an existing FuzzyFields object.
   *
   * @param fields p_fields: The FuzzyFields object in which to store the shared values.
   * @param index p_index: The index to store the shared values at.
   * @param array_size p_array_size: The size of the fuzzy field's boolean array.
   */
  void ReceiveFuzzyField(FuzzyFields &fields, size_t index) {
    if (proxy->GetPRole() != helper) {
      assertm(index < fields.size, "Index within fields size");
      FuzzyField field = ReceiveFuzzyField();
      fields.has_value[index] = field.has_value;
      fields.hamming_weights[index] = field.hamming_weight;
      std::copy(
        field.boolean_array.get(),
        field.boolean_array.get() + boolean_array_byte_size_,
        fields.boolean_arrays.get() + GetArrayIndex(index)
      );
    }
  }

  /**
   * @brief shares an exact field with the other proxy.
   * @param value the value to share
   * @param has_value whether a value actually exists for the given field
   * @return An OptionalValue object containing the shares for the value and has_value
   */
  OptionalValue ShareExactField(std::string &value, bool has_value = true) {
    if (proxy->GetPRole() != helper) {
      uint64_t converted_boolean = ConvertToUint64(has_value, bit_precision_);
      OptionalValue field{converted_boolean - proxy->GenerateCommonRandom(), Hash(value) - proxy->GenerateCommonRandom()};
      return field;
    } else {
      return {};
    }
  }

  /**
   * @brief shares an exact field with the other proxy.
   * The resulting shares are inserted into an existing OptionalValues object.
   * @param fields where the resulting shares should be stored
   * @param index the index at which the shares will be inserted
   * @param value the value to share
   * @param has_value whether a value actually exists for the given field
   */
  void ShareExactField(OptionalValues &fields, size_t index, std::string &value, bool has_value = true) {
    if (proxy->GetPRole() != helper) {
      assertm(index < fields.size, "Index within bound for fields");
      uint64_t converted_boolean = ConvertToUint64(has_value, bit_precision_);
      fields.has_value[index] = converted_boolean - proxy->GenerateCommonRandom();
      fields.values[index] = Hash(value) - proxy->GenerateCommonRandom();
    }
  }

  /**
   * @brief receives an exact field from the other proxy.
   * @return the shares of the exact field.
   */
  OptionalValue ReceiveExactField() {
    if (proxy->GetPRole() != helper) {
      OptionalValue field{proxy->GenerateCommonRandom(), proxy->GenerateCommonRandom()};
      return field;
    } else {
      return {};
    }
  }

  /**
   * @brief receives an exact field from the other proxy.
   * The resulting shares are inserted into an existing OptionalValues object.
   * @param fields where the resulting shares should be stored
   * @param index the index at which the shares will be inserted
   */
  void ReceiveExactField(OptionalValues &fields, size_t index) {
    if (proxy->GetPRole() != helper) {
      assertm(index < fields.size, "Index out of bound for fields");
      fields.has_value[index] = proxy->GenerateCommonRandom();
      fields.values[index] = proxy->GenerateCommonRandom();
    }
  }

  [[nodiscard]] size_t GetExactIndex(size_t index) const {
    return index * exact_fields_per_record_;
  }

  [[nodiscard]] size_t GetFuzzyIndex(size_t index) const {
    return index * fuzzy_fields_per_record_;
  }

  [[nodiscard]] Records ConstructRecords(size_t size, bool initialise_arrays = true) const {
    return {fuzzy_fields_per_record_, exact_fields_per_record_, boolean_array_byte_size_, size, initialise_arrays};
  }


  /**
  * @brief Convenience method to return the pointer to the start of a secret shared boolean array.
  *
  * For performance reasons, all boolean arrays are stored sequentially in a single array. This method allows for the easy retrieval of the start of a boolean array given its index.
  * @param index p_index: The index of the boolean array to retrieve.
  * @return uint8_t* A pointer to the start of the boolean array,  or nullptr if the array does not exist/the index is out of bounds.
  */
  [[nodiscard]] size_t GetArrayIndex(size_t index) const {
    return index*boolean_array_byte_size_;
  }

    // TODO test
  /**
  * @brief Compute the pairwise similarity score between two sets of records.
  *
  * This is doing a 1 by 1 comparison so that index 0 is compared to index 0,  index 1 to index 1 and so on.
  * @param records1 p_records1: The first set of records
  * @param records2 p_records2: The second set of records
  * @return Scores
  */
  Scores ComputeRecordSimilarities(const Records &records1, const Records &records2) {
    bool is_not_helper = proxy->GetPRole() != helper;
    Scores total_scores(records1.size,is_not_helper,true);
    if (is_not_helper) {
      Scores binary_scores = ComputeBinarySimilarities(records1.exact_fields,  records2.exact_fields);
      Scores fuzzy_field_scores = ComputeDiceScores(records1.fuzzy_fields,  records2.fuzzy_fields); 
      size_t exact_index, fuzzy_index;
      for (size_t record_index = 0; record_index < records1.size; record_index++) {
        exact_index = GetExactIndex(record_index);
        fuzzy_index = GetFuzzyIndex(record_index);
        for (size_t i = 0; i < exact_fields_per_record_; i++) {
          total_scores.numerators[record_index]
          += LocalMultiply(binary_scores.numerators[exact_index + i], exact_field_weights_[i], bit_precision_);
          total_scores.denominators[record_index]
          += LocalMultiply(binary_scores.denominators[exact_index + i], exact_field_weights_[i], bit_precision_);
        }
        for (size_t i = 0; i < fuzzy_fields_per_record_; i++) {
          total_scores.numerators[record_index]
          += LocalMultiply(fuzzy_field_scores.numerators[fuzzy_index + i], fuzzy_field_weights_[i], bit_precision_);
          total_scores.denominators[record_index]
          += LocalMultiply(fuzzy_field_scores.denominators[fuzzy_index + i], fuzzy_field_weights_[i], bit_precision_);
        }
      }
    } else {                                                //  helper
      ComputeBinarySimilarities(records1.exact_fields,  records2.exact_fields);
      ComputeDiceScores(records1.fuzzy_fields,  records2.fuzzy_fields);
    }
    return total_scores;
  }

  /**
  * @brief Computes whether the maximum score of a record passes the match threshold.
  *
  * If no match was found,  the match index is set to zero to avoid leaking information.
  * @param max_scores p_max_scores: Maximum scores per record
  * @return OptionalValues the index of the matched record
  */
  OptionalValues GetMatches(MaxScores &max_scores) {
    if (proxy->GetPRole() !=  helper) {
      std::unique_ptr<uint64_t[]> threshold = std::make_unique<uint64_t[]>(max_scores.size);
      for (size_t i = 0; i < max_scores.size; i++) {
        threshold[i] = LocalMultiply(match_threshold_, max_scores.denominators[i], bit_precision_);
      }

      OptionalValues matches(max_scores.size, false);
      matches.has_value = std::unique_ptr<uint64_t[]>(Compare(
        proxy.get(),
        max_scores.numerators.get(),
        threshold.get(),
        max_scores.size,
        bit_precision_
      ));

    matches.values = std::unique_ptr<uint64_t[]>(Multiply(
        proxy.get(),
        matches.has_value.get(),
        max_scores.indices.get(),
        max_scores.size,
        bit_precision_
      ));
/** This is for testing the correctness of the scores
 * If you want to get exact scores put below back into work
 * Also the corresponding line in helper
 * has_values will have exact scores instead of the comparison result with threshold
 * */
//      matches.has_value = std::unique_ptr<uint64_t[]>(Normalise(
//              proxy.get(),
//                max_scores.numerators.get(),
//                max_scores.denominators.get(),
//                max_scores.size,
//                bit_precision_
//        ));

      return matches;
    } else {                                                // helper
      Compare(proxy.get(),  nullptr,  nullptr,  max_scores.size,  bit_precision_);
      Multiply(proxy.get(),  nullptr,  nullptr,  max_scores.size,  bit_precision_);
      /**Turn below on to get the exact scores for testing*/
      //Normalise(proxy.get(),  nullptr,  nullptr,  max_scores.size,  bit_precision_);
      return OptionalValues(max_scores.size, false);
    }
  }

  /**
  * @brief For each record in records1, return whether a match is found in records2 and the index of the matching record.
  *
  * @param records1 p_records1: the first set of records
  * @param records2 p_records2: the second set of records
  * @param parallelism p_parallelism: how many records should be processed at the same time. This has memory implications.
  * @return OptionalValues
  */
  OptionalValues FindAllMatches(Records &records1,  Records &records2,  size_t parallelism) {
    parallelism = std::min<size_t>(records1.size, parallelism);
    size_t length = records2.size*parallelism;
    size_t remainder = records1.size % parallelism;
    bool is_not_helper = proxy->GetPRole() != helper;
    Records single_records = ConstructRecords(length, is_not_helper);
    Records all_records = ConstructRecords(length, is_not_helper);
    MaxScores max_scores(records1.size, is_not_helper);
    if (is_not_helper) {
      // copy over records2
      for (size_t i = 0; i < parallelism; i++) {
        CopyRecords(records2,  all_records,  records2.size, 0,  records2.size*i);
      }
      // generate indices
      std::unique_ptr<uint64_t[]> indices = std::make_unique<uint64_t[]>(length);
      for (size_t i = 0; i < parallelism; i++) {
        for (size_t ii = 0; ii < records2.size; ii++) {
          indices[i*records2.size+ii] = CreateShare(ii);
        }
      }
    for (size_t block = 0; block < (records1.size/parallelism); block++) {
        // duplicate relevant records for vectorised operation
        for (size_t i = 0; i < parallelism; i++) {
          DuplicateRecord(records1,  single_records, records2.size, block*parallelism+i,  i*records2.size);
        }

        Scores scores = ComputeRecordSimilarities(single_records,  all_records);

        MaxScores block_max = ComputeMaxScores(scores, records2.size, parallelism, indices);

        // copy the results of the block into the max_scores vector:
        std::copy(
          block_max.indices.get(),
          block_max.indices.get()+parallelism,
          max_scores.indices.get()+block*parallelism
        );
        std::copy(
          block_max.numerators.get(),
          block_max.numerators.get()+parallelism,
          max_scores.numerators.get()+block*parallelism
        );
        std::copy(
          block_max.denominators.get(),
          block_max.denominators.get()+parallelism,
          max_scores.denominators.get()+block*parallelism
        );
      }
      if (remainder != 0) {
        // duplicate relevant records for vectorised operation
        for (size_t i = 0; i < remainder; i++) {
          DuplicateRecord(records1,  single_records, records2.size, records1.size-remainder+i,  i*records2.size);
        }
        single_records.size = remainder*records2.size;
        single_records.exact_fields.size = single_records.size * exact_fields_per_record_;
        single_records.fuzzy_fields.size = single_records.size * fuzzy_fields_per_record_;
        Scores scores = ComputeRecordSimilarities(single_records,  all_records);
        MaxScores block_max = ComputeMaxScores(scores, records2.size, remainder, indices);
        // copy the results of the block into the max_scores vector:
        std::copy(
          block_max.indices.get(),
          block_max.indices.get()+remainder,
          max_scores.indices.get()+records1.size-remainder
        );
        std::copy(
          block_max.numerators.get(),
          block_max.numerators.get()+remainder,
          max_scores.numerators.get()+records1.size-remainder
        );
        std::copy(
          block_max.denominators.get(),
          block_max.denominators.get()+remainder,
          max_scores.denominators.get()+records1.size-remainder
        );
      }


      auto matches = GetMatches(max_scores);

      return matches;
    } else {                                                //  helper
      std::unique_ptr<uint64_t[]> indices(nullptr);
      for (size_t block = 0; block < (records1.size/parallelism); block++) {
        Scores scores = ComputeRecordSimilarities(single_records,  all_records);
        ComputeMaxScores(scores,  records2.size,  parallelism, indices);
      }
      if (remainder != 0) {
        single_records.size = remainder*records2.size;
        single_records.exact_fields.size = single_records.size * exact_fields_per_record_;
        single_records.fuzzy_fields.size = single_records.size * fuzzy_fields_per_record_;
        Scores scores = ComputeRecordSimilarities(single_records,  all_records);
        MaxScores block_max = ComputeMaxScores(scores, records2.size, remainder, indices);
      }
      return GetMatches(max_scores);
    }
  }

private:
  size_t exact_fields_per_record_;
  size_t fuzzy_fields_per_record_;
  size_t boolean_array_size_;
  size_t boolean_array_byte_size_;
  int bit_precision_;
  uint64_t match_threshold_;
  std::vector<uint64_t> exact_field_weights_;
  std::vector<uint64_t> fuzzy_field_weights_;
  std::shared_ptr<Party> proxy;
  CryptoPP::BLAKE2b hasher;

  static int MapCharToInt(char character) {
    char lowerCase = tolower(character);
    if ((97 <= lowerCase) & (lowerCase <= 122)) {
      return lowerCase - 97;
    }
    switch(lowerCase) {
      case ' ' :
        return 26;
      case '-' :
        return 27;
      case '.' :
        return 28;
      default:
        return 29;
    }
  }

  /**
  * @brief Copy a number of records from one Records object into another one.
  *
  * @param source p_source: The Records object to copy records from.
  * @param destination p_destination: The Records objects to copy records into.
  * @param source_index p_source_index: The first index of the records to copy in source.
  * @param dest_index p_dest_index: The first index of where to store the records in destination.
  * @param count p_count: How many records to copy over.
  */
  void CopyRecords(
    Records &source,
    Records &destination,
    size_t count,
    size_t source_index = 0,
    size_t dest_index = 0
  ) const {
    assertm(dest_index+count <= destination.size, "Records can't be copied because destination is too small");
    size_t source_start, source_end, destination_start;
    source_start = source_index * exact_fields_per_record_;
    source_end = (source_index+count) * exact_fields_per_record_;
    destination_start = dest_index * exact_fields_per_record_;
    std::copy(
      source.exact_fields.has_value.get()+source_start,
      source.exact_fields.has_value.get()+source_end,
      destination.exact_fields.has_value.get()+destination_start
    );
    std::copy(
      source.exact_fields.values.get()+source_start,
      source.exact_fields.values.get()+source_end,
      destination.exact_fields.values.get()+destination_start
    );
    source_start = source_index*fuzzy_fields_per_record_;
    source_end = (source_index+count)*fuzzy_fields_per_record_;
    destination_start = dest_index*fuzzy_fields_per_record_;
    std::copy(
      source.fuzzy_fields.boolean_arrays.get()+source_start*boolean_array_byte_size_,
      source.fuzzy_fields.boolean_arrays.get()+source_end*boolean_array_byte_size_,
      destination.fuzzy_fields.boolean_arrays.get()+destination_start*boolean_array_byte_size_
    );
    std::copy(
      source.fuzzy_fields.has_value.get() + source_start,
      source.fuzzy_fields.has_value.get() + source_end,
      destination.fuzzy_fields.has_value.get() + destination_start
    );
    std::copy(
      source.fuzzy_fields.hamming_weights.get() + source_start,
      source.fuzzy_fields.hamming_weights.get() + source_end,
      destination.fuzzy_fields.hamming_weights.get() + destination_start
    );
  }

  /**
  * @brief Copy a record from one Records to another multiple times.
  *
  * @param source p_source: The Records object to copy the record from.
  * @param destination p_destination: The Records objects to copy the record into.
  * @param count p_count: How often to duplicate the record.
  * @param source_index p_source_index: The index of the record to copy in source.
  * @param dest_index p_dest_index: The first index of where to store the duplicated records in destination.
  */
  void DuplicateRecord(
    Records &source,
    Records &destination,
    size_t count,
    size_t source_index = 0,
    size_t dest_index = 0
  ) {
    assertm(dest_index + count <= destination.size, "Record can't be duplicated because destination is too small");
    // copy source record into destination
    if (count > 0) {
      CopyRecords(source,  destination,  1, source_index,  dest_index);
    }
    // duplicate existing values until count is reached
    size_t copied_count = 1;
    size_t length = 1;
    while (copied_count < count) {
      CopyRecords(destination,  destination,  length, dest_index, dest_index + copied_count);
      copied_count += length;
      length = std::min(copied_count, count - copied_count);
    }
  }

    /**
   * @brief Compute the Dice coefficient (i.e. similarity score) between two sets of FuzzyFields.
   *
   * Only fields1.size coefficients are computed.
   * @param fields1 p_fields1: The first FuzzyFields set
   * @param fields2 p_fields2: The second FuzzyFields set
   * @return Scores
   */
  Scores ComputeDiceScores(const FuzzyFields &fields1, const FuzzyFields &fields2) {
    bool is_not_helper = proxy->GetPRole() != helper;
    Scores scores(fields1.size, is_not_helper);
    if (is_not_helper) {
      std::unique_ptr<uint64_t[]> hamming_weight_sum = std::make_unique<uint64_t[]>(fields1.size);
      for (size_t i = 0; i < fields1.size; i++) {
        hamming_weight_sum[i] = fields1.hamming_weights[i] + fields2.hamming_weights[i];
      }
      // since has_value can be a fraction (e.g. 2 out of 3 attributes have a value: has_value=0.66), we need to select the smaller of the two values here
      std::unique_ptr<uint64_t[]> selection(Compare(
        proxy.get(),
        fields1.has_value.get(),
        fields2.has_value.get(),
        fields1.size,
        bit_precision_
      ));

      std::unique_ptr<uint64_t[]> has_value(Multiplex(
        proxy.get(),
        fields1.has_value.get(),
        fields2.has_value.get(),
        selection.get(),
        fields1.size,
        bit_precision_
      ));

       selection = nullptr; // clean memory
      scores.denominators = std::unique_ptr<uint64_t[]>(Multiply(
        proxy.get(),
        hamming_weight_sum.get(),
        has_value.get(),
        fields1.size,
        bit_precision_
      ));

      std::unique_ptr<uint8_t[]> shared_array(And(
        proxy.get(),
        fields1.boolean_arrays.get(),
        fields2.boolean_arrays.get(),
        boolean_array_byte_size_ * fields1.size
      ));

      std::unique_ptr<uint64_t[]> arithmetic_array(XorToArithmetic2(
        proxy.get(),
        shared_array.get(),
        boolean_array_byte_size_*fields1.size*8,
        bit_precision_
      ));
      shared_array = nullptr;

      std::unique_ptr<uint64_t[]> twice_bits_in_common = std::make_unique<uint64_t[]>(fields1.size);
      for (size_t i = 0; i < fields1.size; i++) {
        twice_bits_in_common[i] = 2 * std::accumulate(
          arithmetic_array.get() + 8*i*boolean_array_byte_size_,
          arithmetic_array.get() + 8*(i+1)*boolean_array_byte_size_,
          (uint64_t) 0
        );
      }
      arithmetic_array = nullptr; // clean memory
      scores.numerators = std::unique_ptr<uint64_t[]>(Multiply(
        proxy.get(),
        twice_bits_in_common.get(),
        has_value.get(),
        fields1.size,
        bit_precision_
      ));

     } else { // helper
      Compare(proxy.get(), nullptr, nullptr, fields1.size, bit_precision_);
      Multiplex(proxy.get(), nullptr, nullptr, nullptr, fields1.size, bit_precision_);
      Multiply(proxy.get(), nullptr, nullptr, fields1.size, bit_precision_);
      And(proxy.get(), nullptr, nullptr, boolean_array_byte_size_ * fields1.size);
      XorToArithmetic2(proxy.get(), nullptr, boolean_array_byte_size_*fields1.size*8,bit_precision_);
      Multiply(proxy.get(), nullptr, nullptr, fields1.size, bit_precision_);
    }
    return scores;
  }


  /**
   * @brief Compute the similarity score {0, 1} between two OptionalValue s.
   *
   * @param field1 p_field1: The first ExactField
   * @param field2 p_field2: The second ExactField
   * @return Score
   */
  Score ComputeBinarySimilarity(OptionalValue value1, OptionalValue value2) {
    Score score{};
    bool is_not_helper = proxy->GetPRole() != helper;
    OptionalValues values1(1, is_not_helper);
    if (is_not_helper) {
      OptionalValues values2(1);
      values1.has_value[0] = value1.has_value;
      values1.values[0] = value1.value;
      values2.has_value[0] = value2.has_value;
      values2.values[0] = value2.value;
      Scores scores = ComputeBinarySimilarities(values1,  values2);
      score.numerator = scores.numerators[0];
      score.denominator = scores.denominators[0];
    } else {
      ComputeBinarySimilarities(values1,  values1);
    }
    return score;
  }

  /**
  * @brief Compute the similarity score {0, 1} between two sets of OptionalValues.
  *
  * There are only as many comparisons as there are values in values1.
  * @param values1 p_values1: A set of OptionalValues
  * @param values2 p_values2: A set of OptionalValues
  * @return Scores
  */
  Scores ComputeBinarySimilarities(const OptionalValues &values1, const OptionalValues &values2) {
    bool is_not_helper = proxy->GetPRole() != helper;
    Scores scores(values1.size, is_not_helper);
    if (is_not_helper) {
      scores.denominators = std::unique_ptr<uint64_t[]>(Multiply(
        proxy.get(),
        values1.has_value.get(),
        values2.has_value.get(),
        values1.size,
        bit_precision_
      ));
      std::unique_ptr<uint64_t[]> is_equal(Equals(
        proxy.get(),
        values1.values.get(),
        values2.values.get(),
        values1.size,
        bit_precision_
      ));
      scores.numerators = std::unique_ptr<uint64_t[]>(Multiply(
        proxy.get(),
        is_equal.get(),
        scores.denominators.get(),
        values1.size,
        bit_precision_
      ));
    } else {
      Multiply(proxy.get(), nullptr, nullptr, values1.size, bit_precision_);
      Equals(proxy.get(), nullptr, nullptr, values1.size, bit_precision_);
      Multiply(proxy.get(), nullptr, nullptr, values1.size, bit_precision_);
    }
    return scores;
  }

  /**
   * @brief For each score, returns score1 >= score2
   * Instead of providing whole Score objects, the individual pointers are provided. This saves copying data in ComputeMaxScore.
   * @param numerator1 p_numerator1:...
   * @param numerator2 p_numerator2:...
   * @param denominator1 p_denominator1:...
   * @param denominator2 p_denominator2:...
   * @param count p_count: number of scores
   * @return uint64_t*
   */
  std::unique_ptr<uint64_t[]> CompareScores(
    const uint64_t *const numerator1,
    const uint64_t *const numerator2,
    const uint64_t *const denominator1,
    const uint64_t *const denominator2,
    size_t count
  ) {
    if (proxy->GetPRole() == helper) {
      Multiply(proxy.get(), nullptr, nullptr, count*2, bit_precision_);
      Compare(proxy.get(), nullptr, nullptr, count, bit_precision_);
      return nullptr;
    } else { // proxy1 and proxy2
      std::unique_ptr<uint64_t[]> compare1 = std::make_unique<uint64_t[]>(count*2);
      std::copy(numerator1,  numerator1+count,  compare1.get());
      std::copy(numerator2,  numerator2+count,  compare1.get()+count);
      std::unique_ptr<uint64_t[]> compare2 = std::make_unique<uint64_t[]>(count*2);
      std::copy(denominator2,  denominator2+count,  compare2.get());
      std::copy(denominator1,  denominator1+count,  compare2.get()+count);
      compare1 = std::unique_ptr<uint64_t[]>(Multiply(proxy.get(), compare1.get(), compare2.get(), count*2, bit_precision_));
      compare2 = nullptr;
      compare1 = std::unique_ptr<uint64_t[]>(Compare(proxy.get(), compare1.get(), compare1.get()+count, count, bit_precision_));
      return compare1;
    }
  }

  static void CopyExactFields(OptionalValues const &source, OptionalValues &destination, size_t source_index, size_t dest_index, size_t count) {
    std::copy(
      source.has_value.get() + source_index,
      source.has_value.get() + source_index + count,
      destination.has_value.get() + dest_index
    );
    std::copy(
      source.values.get() + source_index,
      source.values.get() + source_index + count,
      destination.values.get() + dest_index
    );
  }

  static void DuplicateExactFields(OptionalValues const &source, OptionalValues &destination, size_t source_index, size_t dest_index, size_t count, size_t duplication_count) {
    assertm(dest_index + count*duplication_count <= destination.size, "Field can't be duplicated because destination is too small");
    assertm(source_index + count <= source.size, "Index out of bound for source");
    // copy source record into destination
    if (duplication_count > 0) {
      CopyExactFields(source,  destination,  source_index,  dest_index, count);
    }
    // duplicate existing values until count is reached
    size_t copied_count = count;
    size_t length = count;
    size_t total_count_to_copy = count * duplication_count;
    while (copied_count < total_count_to_copy) {
      CopyExactFields(destination,  destination,  dest_index, dest_index + copied_count, length);
      copied_count += length;
      length = std::min(copied_count, total_count_to_copy - copied_count);
    }
  }

  /**
  * @brief Compute the maximum score for multiple blocks of scores.
  *
  * @param scores p_scores
  * @param scores_per_maximum p_scores_per_maximum: How large the blocks of scores are.
  * @param maxima_to_compute p_maxima_to_compute: How many blocks of scores should be processed.
  * @param start_index p_start_index: The start index of the first block of scores.
  * @return MaxScores
  */
  MaxScores ComputeMaxScores(const Scores &scores, size_t scores_per_maximum, size_t maxima_to_compute,  std::unique_ptr<uint64_t[]> &indices, size_t start_index = 0){
    /** MAIN IDEA:
     Repeatedly compare the first half of the remaining scores to the second
     half, using Multiplex to select the larger of the two scores.
     Have an array contain secret shared indices that are selected
     with the same selection vector to also get the argmax.
     */
    size_t old_length;
    size_t half_length;
    size_t length = scores_per_maximum;
    size_t total_length = scores_per_maximum*maxima_to_compute;
    size_t old_total_length;
    bool remainder_exists = false;
    bool is_not_helper = proxy->GetPRole() != helper;
    MaxScores max(maxima_to_compute, is_not_helper);
    if (is_not_helper) {
      // generate the initial data structures and values:
      size_t denominator_index = total_length;
      size_t old_denominator_index;
      size_t indices_index = total_length*2;
      size_t old_indices_index;
      size_t index_offset;
      std::unique_ptr<uint64_t[]> temp = std::make_unique<uint64_t[]>(total_length*3);
      std::unique_ptr<uint64_t[]> combined_selection = std::make_unique<uint64_t[]>(total_length*3);
      Scores remainder(maxima_to_compute);
      std::unique_ptr<uint64_t[]> remainder_index = std::make_unique<uint64_t[]>(maxima_to_compute);
      // combined_first contains the first half of the numerators,  denominators and indices
      std::unique_ptr<uint64_t[]> combined_first = std::make_unique<uint64_t[]>(total_length * 3);
      std::unique_ptr<uint64_t[]> combined_second = std::make_unique<uint64_t[]>(total_length * 3);
      //  we begin by copying the contents of scores into temp:
      std::copy(
        scores.numerators.get() + start_index,
        scores.numerators.get() + total_length + start_index,
        temp.get()
      );
      std::copy(
        scores.denominators.get() + start_index,
        scores.denominators.get() + total_length + start_index,
        temp.get()+total_length
      );
      std::copy(
        indices.get() + start_index,
        indices.get() + total_length + start_index,
        temp.get() + 2*total_length
      );
      //iteratively halve the array by comparing the first with the second half:
      while (length > 1) {
        // right now, temp has length maxima_to_compute*length*3 and contains the selected values (numerator, denominator, index)
        half_length = length/2;
        old_length = length;
        old_denominator_index = denominator_index;
        old_indices_index = indices_index;
        old_total_length = total_length;
        if (old_length % 2 == 1) {
          if (remainder_exists) {                           // the remainders are appended to complement combined_first
            length = half_length+1;
            index_offset = 0;
          } else {                                          //  the remainder is not copied over and stored for later
            length = half_length;
            index_offset = 1;
            }
        } else {
          length = half_length;
          index_offset = 0;
        }
        total_length = length*maxima_to_compute;
        denominator_index = total_length;
        indices_index = total_length*2;
        for (size_t i = 0; i < maxima_to_compute; i++) {
          if (old_length %  2 == 1) {
            if (remainder_exists) {                         // the remainders are added to the last position of combined_first. We later take care not to copy over that value.
              combined_first[(i+1)*length - 1] = remainder.numerators[i];
              combined_first[denominator_index + (i+1)*length - 1] = remainder.denominators[i];
              combined_first[indices_index + (i+1)*length - 1] = remainder_index[i];
            } else {                                        // the remainder is stored for later
              remainder.numerators[i] = temp[(i+1)*old_length - 1];
              remainder.denominators[i] = temp[old_total_length + (i+1)*old_length - 1];
              remainder_index[i] = temp[old_total_length*2 + (i+1)*old_length - 1];
            }
          }
          // copy numerators:
          std::copy(
            temp.get()+i*old_length,
            temp.get()+i*old_length+half_length, // if old_length is odd, this will select the smaller half (e.g. old_length=9 -> half_length=4, so only 4 values are copied over)
            combined_first.get()+i*length
           );
          std::copy(
            temp.get()+i*old_length+half_length,
            temp.get()+(i+1) * old_length - index_offset, // if we are storing the remainder, index_offset is 1, so we are not copying the remainder over
            combined_second.get()+i*length
           );
          //  copy denominators:
          std::copy(
            temp.get() + old_denominator_index + i*old_length,
            temp.get() + old_denominator_index + i*old_length + half_length,
            combined_first.get()+denominator_index + i*length
          );
          std::copy(
            temp.get() + old_denominator_index + i*old_length + half_length,
            temp.get() + old_denominator_index + (i+1)*old_length-index_offset,
            combined_second.get() + denominator_index + i*length
          );
          //  copy indices:
          std::copy(
            temp.get() + old_indices_index + i*old_length,
            temp.get() + old_indices_index + i*old_length + half_length,
            combined_first.get()+indices_index + i*length
          );
          std::copy(
            temp.get() + old_indices_index + i*old_length + half_length,
            temp.get() + old_indices_index + (i+1)*old_length-index_offset,
            combined_second.get() + indices_index + i*length
          );
        }
        if (old_length % 2 == 1) { // this can only be done now because we needed this information in the loop.
          remainder_exists = !remainder_exists;
        }

        temp = std::unique_ptr<uint64_t[]>(CompareScores(
          combined_first.get(),
          combined_second.get(),
          combined_first.get() + denominator_index,
          combined_second.get() + denominator_index,
          total_length
        ));
        // temp now contains the selection vector for multiplexing which is duplicated twice to select all values:
        std::copy(temp.get(),  temp.get() + total_length,  combined_selection.get());
        std::copy(temp.get(),  temp.get() + total_length,  combined_selection.get() + denominator_index);
        std::copy(temp.get(),  temp.get() + total_length,  combined_selection.get() + indices_index);
//        std::unique_ptr<double[]> numerators(ReconstructDouble(combined_selection.get(), total_length));
//        std::unique_ptr<double[]> denominators(ReconstructDouble(combined_selection.get()+denominator_index, total_length));
//        cout << "Block " << maxima_to_compute <<endl;
//        for (size_t i = 0; i < 10; i++) {
//            std::cout << i << "\t" << numerators[i] << "/" << denominators[i] << std::endl;
//        }
        temp = std::unique_ptr<uint64_t[]>(Multiplex(
          proxy.get(),
          combined_second.get(),
          combined_first.get(),
          combined_selection.get(),
          total_length*3,
          bit_precision_
        ));
        // temp now contains the selected values and has length maxima_to_compute*length*3
      }
      if (remainder_exists) {
        // do one last operation on the last score and the remainder:
        combined_first = std::move(temp);
        std::copy(remainder.numerators.get(),  remainder.numerators.get() + maxima_to_compute,  combined_second.get());
        std::copy(
          remainder.denominators.get(),
          remainder.denominators.get() + maxima_to_compute,
          combined_second.get() + denominator_index
         );
        std::copy(remainder_index.get(),  remainder_index.get() + maxima_to_compute,  combined_second.get() + indices_index);
        temp = std::unique_ptr<uint64_t[]>(CompareScores(
          combined_first.get(),
          combined_second.get(),
          combined_first.get() + denominator_index,
          combined_second.get() + denominator_index,
          total_length
         ));
        std::copy(temp.get(),  temp.get() + total_length,  combined_selection.get());
        std::copy(temp.get(),  temp.get() + total_length,  combined_selection.get() + denominator_index);
        std::copy(temp.get(),  temp.get() + total_length,  combined_selection.get() + indices_index);
        temp = std::unique_ptr<uint64_t[]>(Multiplex(
          proxy.get(),
          combined_second.get(),
          combined_first.get(),
          combined_selection.get(),
          total_length*3,
          bit_precision_
        ));
      }
      std::copy(temp.get(),  temp.get() + denominator_index,  max.numerators.get());
      std::copy(temp.get() + denominator_index,  temp.get() + indices_index,  max.denominators.get());
      std::copy(temp.get() + indices_index,  temp.get() + total_length*3,  max.indices.get());
      return max;
    } else { // helper
      while (length > 1) {
        // right now, temp has length maxima_to_compute*length*3 and contains the selected values
        half_length = length/2;
        old_length = length;
        if (old_length % 2 == 1) {
          if (remainder_exists) {                           // the remainders are appended to complement combined_first
            length = half_length+1;
          } else {                                          //  the remainder is not copied over and stored for later
            length = half_length;
          }
          remainder_exists = !remainder_exists;
        } else {
          length = half_length;
        }
        total_length = length*maxima_to_compute;
        CompareScores(
          nullptr,
          nullptr,
          nullptr,
          nullptr,
          total_length
        );
        Multiplex(
          proxy.get(),
          nullptr,
          nullptr,
          nullptr,
          total_length*3,
          bit_precision_
        );
      }
      if (remainder_exists) {
        // do one last operation on the last score and the remainder:
        CompareScores(
          nullptr,
          nullptr,
          nullptr,
          nullptr,
          total_length
         );
        Multiplex(
          proxy.get(),
          nullptr,
          nullptr,
          nullptr,
          total_length*3,
          bit_precision_
        );
      }
      return max;
    }
  }
};

#endif                                                      // PARTY_RL
