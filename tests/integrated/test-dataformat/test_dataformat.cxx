
#include <bout/output.hxx>
#include <dataformat.hxx>

int main() {
  const std::string izfilename = "sample.nc";

  // Create a file format handler
  auto izfile = data_format(izfilename.c_str());

  izfile->openr(izfilename);

  if (!izfile->is_valid()) {
    output << "\tERROR: Could not open file " << izfilename << endl;
  }

  izfile->close();

  return 0;
}
