#include "../PortableMap.hpp"

int main(int argc, char * argv[])
{
  if ( argc != 2 ) {
    printf("Usage:\n%s portable_map_file\n", argv[0]);
    return 0;
  }

  PortableMap pnm;
  std::string file(argv[1]);
  pnm.load(file);
  pnm.print();
  pnm.writeBinary("binary_"+file);
  pnm.writeAscii("ascii_"+file);

  return 0;
}
