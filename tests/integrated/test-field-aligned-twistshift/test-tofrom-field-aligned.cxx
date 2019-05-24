#include "bout.hxx"
#include "fft.hxx"
#include "field_factory.hxx"

std::pair<Field3D, Field3D> getModes(const Field3D& f) {
  const int nz = mesh->LocalNz;

  Array<dcomplex> buffer(nz/2 + 1);

  Field3D real{zeroFrom(f)};
  Field3D imag{zeroFrom(f)};

  BOUT_FOR_SERIAL(i, f.getRegion2D("RGN_ALL")) {
    rfft(&f(i, 0), nz, &buffer[0]);

    for (int k=0; k<nz/2+1; k++) {
      real(i, k) = buffer[k].real();
      imag(i, k) = buffer[k].imag();
    }
  }

  return std::make_pair(real, imag);
}

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  Field3D f = FieldFactory::get()->create3D("f");
  mesh->communicate(f);
  f = toFieldAligned(f);

  Field3D f_tofrom = copy(f);
  mesh->communicate(f_tofrom);

  auto modes = getModes(f);
  auto modes_tofrom = getModes(f_tofrom);

  dump.addOnce(f, "f");
  dump.addOnce(modes.first, "freal");
  dump.addOnce(modes.second, "fimag");
  dump.addOnce(f_tofrom, "f_tofrom");
  dump.addOnce(modes_tofrom.first, "freal_tofrom");
  dump.addOnce(modes_tofrom.second, "fimag_tofrom");

  dump.write();

  BoutFinalise();
  return 0;
}
