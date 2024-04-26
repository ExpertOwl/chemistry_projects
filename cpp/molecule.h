#include <vector>

struct Molecule
{int natom;
int charge;
std::vector<int> zvals;
std::vector<int> x_coords;
std::vector<int> y_coords;
std::vector<int> z_coords;

void print_geom();
};