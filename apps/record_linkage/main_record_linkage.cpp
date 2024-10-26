#include <cstdlib>
#include <iostream>
#include <tuple>
#include "record_linkage.h"
#include "../../dependencies/fast-cpp-csv-parser/csv.h"

static const int REQUIRED_ARGS_COUNT = 7;

struct EpiLink{
  std::vector<std::string> ids;
  Records records;
};

static void PrintError() {
  std::clog
    << "Required arguments not provided. These are the required (and optional) arguments:\n"
    << "1: role (proxy1, proxy2 or helper)\n"
    << "2 & 3: helper's port and IP\n"
    << "4 & 5: proxy1's port and IP\n"
    << "6: record count\n"
    << "7: parallelism" << std::endl;
  exit(EXIT_FAILURE);
}

int SuppressStdOut() {
  fflush(stdout);
  int stdout_file_descriptor = dup(1);
  int null_file_descriptor = open("/dev/null", O_WRONLY);
  // check null_file_descriptor for error omitted
  dup2(null_file_descriptor, 1);
  close(null_file_descriptor);

  return stdout_file_descriptor;
}

void ResumeStdOut(int file_descriptor) {
  fflush(stdout);
  dup2(file_descriptor, 1);
  close(file_descriptor);
}

EpiLink ReadAndShareEpiLink(
  RecordLinkage &proxy,
  const std::string& file_name,
  size_t record_count,
  const std::string& empty_value = ""
) {
  if (proxy.GetPRole() == helper) {
    std::cerr << "Helper should not read data." << std::endl;
    exit(EXIT_FAILURE);
  }
  io::CSVReader<9> reader(file_name);
  reader.read_header(
    io::ignore_extra_column,
    "rec-id",
    "first-name",
    "surname",
    "birth-name",
    "city",
    "birth-day",
    "birth-month",
    "birth-year",
    "postcode"
  );
  size_t line = 0;
  std::string rec_id, first_name, surname, birth_name, birth_day, birth_month, birth_year, postcode, city;
  Records records = proxy.ConstructRecords(record_count);
  std::vector<std::string> ids(record_count);
  double fraction_of_values;
  BooleanArray boolean_array;
  std::string combined_string;
  size_t fuzzy_index = 0;
  size_t exact_index = 0;
  std::stringstream string_stream{};
  while (line < record_count and reader.read_row(
    rec_id,
    first_name,
    surname,
    birth_name,
    city,
    birth_day,
    birth_month,
    birth_year,
    postcode
  )) {
    ids[line] = rec_id;
    fraction_of_values = ((first_name != empty_value) + (surname != empty_value) + (birth_name != empty_value))/3.0;
    string_stream.str(first_name);
    string_stream<< first_name << " " << surname << " " << birth_name;
    boolean_array = RecordLinkage::ConvertStringToBooleanArray(string_stream.str());
    proxy.ShareFuzzyField(records.fuzzy_fields, fuzzy_index++, boolean_array, fraction_of_values);
    boolean_array = RecordLinkage::ConvertStringToBooleanArray(city);
    proxy.ShareFuzzyField(records.fuzzy_fields, fuzzy_index++, boolean_array, city != empty_value);
    proxy.ShareExactField(records.exact_fields, exact_index++, postcode, postcode != empty_value);
    proxy.ShareExactField(records.exact_fields, exact_index++, birth_day, birth_day != empty_value);
    proxy.ShareExactField(records.exact_fields, exact_index++, birth_month, birth_month != empty_value);
    proxy.ShareExactField(records.exact_fields, exact_index++, birth_year, birth_year != empty_value);
    line++;
  }
  if (line < record_count) {
    std::cerr << "CSV file does not have the requested number of lines. Returned " << line << " records instead" << std::endl;
    records.size = line;
    records.exact_fields.size = line * proxy.GetExactFieldsPerRecord();
    records.fuzzy_fields.size = line * proxy.GetExactFieldsPerRecord();
  }
  proxy.SendSize(line);
  return {ids, std::move(records)};
}

Records ReceiveEpiLink(RecordLinkage &proxy, Role sender) {
  if (proxy.GetPRole() == sender) {
    std::cerr << "Tried to receive data from itself" << std::endl;
    exit(EXIT_FAILURE);
  }
  size_t record_count = proxy.ReceiveSize(sender);
  Records records = proxy.ConstructRecords(record_count);
  if (proxy.GetPRole() != helper) {
    size_t fuzzy_index = 0;
    size_t exact_index = 0;
    for (size_t i = 0; i < record_count; i++) {
      proxy.ReceiveFuzzyField(records.fuzzy_fields, fuzzy_index++);
      proxy.ReceiveFuzzyField(records.fuzzy_fields, fuzzy_index++);
      proxy.ReceiveExactField(records.exact_fields, exact_index++);
      proxy.ReceiveExactField(records.exact_fields, exact_index++);
      proxy.ReceiveExactField(records.exact_fields, exact_index++);
      proxy.ReceiveExactField(records.exact_fields, exact_index++);
    }
  }
  return records;
}

int main(int argc, char* argv[]) {
  if (argc < REQUIRED_ARGS_COUNT) {
    PrintError();
  }
  //parse arguments:
  string role_string(argv[1]);
  Role proxy_role;
  if (role_string == "proxy1") {
    proxy_role = proxy1;
  } else if (role_string == "proxy2") {
    proxy_role = proxy2;
  } else if (role_string == "helper") {
    proxy_role = helper;
  } else {
    PrintError(); // this terminates the program
  }
  uint16_t helper_port = strtol(argv[2], nullptr, 10);
  string helper_ip(argv[3]);
  uint16_t p1_port = strtol(argv[4], nullptr, 10);
  string p1_ip(argv[5]);
  size_t record_count = strtol(argv[6], nullptr, 10);
  size_t parallelism = strtol(argv[7], nullptr, 10);
  std::vector<double> fuzzy_field_weights;
  fuzzy_field_weights.push_back(13.16); //combined name
  fuzzy_field_weights.push_back(5.44);  // city name (6.58)
  std::vector<double> exact_field_weights;
  exact_field_weights.push_back(5.47);      //postcode (6.58)
  exact_field_weights.push_back(3.6);       //birthday  (4.9)
  exact_field_weights.push_back(3.59);      //birth-month (3.58)
  exact_field_weights.push_back(3.52);      //birth-year (5.12)
  int file_descriptor = SuppressStdOut();
  RecordLinkage party(
          proxy_role,
          helper_port,
          helper_ip,
          p1_port,
          p1_ip,
          exact_field_weights,
          fuzzy_field_weights,
          0.8
          );
  ResumeStdOut(file_descriptor);
  cout << "Parties created" << endl;
  std::vector<std::string> ids;
  EpiLink epi_link{std::vector<std::string>(0), party.ConstructRecords(0, false)};
  Records shared_records = party.ConstructRecords(0, false);
  Records shared_records2 = party.ConstructRecords(0, false);
  std::cout << std::filesystem::current_path() << std::endl;

    switch(proxy_role) {
    case (proxy1):
      epi_link = ReadAndShareEpiLink(party, "../apps/record_linkage/dataset1.csv", record_count);
      ids = std::move(epi_link.ids);
      shared_records = std::move(epi_link.records);
      shared_records2 = ReceiveEpiLink(party, proxy2);
      break;
    case (proxy2):
      shared_records = ReceiveEpiLink(party, proxy1);
      epi_link = ReadAndShareEpiLink(party, "../apps/record_linkage/dataset2.csv", record_count);
      ids = std::move(epi_link.ids);
      shared_records2 = std::move(epi_link.records);
      break;
    case (helper): // the helper receives nothing but the sizes from proxy1 and proxy2.
      shared_records = ReceiveEpiLink(party, proxy1);
      shared_records2 = ReceiveEpiLink(party, proxy2);
      break;
  }
    cout << "Sharing completed, calling Record Linkage..." << endl;

    auto start = chrono::high_resolution_clock::now();
    OptionalValues matches = party.FindAllMatches(shared_records, shared_records2, parallelism);
    auto end = chrono::high_resolution_clock::now();
    double totaltime =
            chrono::duration_cast<chrono::nanoseconds>(end - start).count()*1e-9;
  if (proxy_role != helper) {
      cout << "Printing the results..." <<endl;
    std::unique_ptr<double[]> is_match(party.ReconstructDouble(matches.has_value.get(), matches.size));
    std::unique_ptr<double[]> match(party.ReconstructDouble(matches.values.get(), matches.size));
      cout << "Index" << "\t" << "isMatch?" << "\t" << "Matching" << endl;

      for (int i = 0; i < matches.size; ++i) {
          cout << i << "\t" << is_match[i] << "\t" << match[i] << endl;
      }
      cout<<"Total Time " << totaltime<<endl;
  }
}
