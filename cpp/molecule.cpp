#include "molecule.h"
#include <iostream>

void Molecule::print_geom()
{
  for(int i=0; i < natom; i++)
    std::cout << zvals[i] << " " << x_coords[i] << " " << y_coords[i] << " " << z_coords[i] << "\n";
}