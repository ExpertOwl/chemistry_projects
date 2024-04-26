// open the filestream, read the first number which is the number of atoms, create list of tuples (vectors) with length = first number. Read in each line one by one and assign xyz's properly

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include "molecule.h"

using namespace std;

int read_natoms(string filename){
  int natom;
  ifstream input(filename);
  input >> natom;
  input.close();
  return natom;
}

int main(){
  cout << "hello world";
  Molecule mol;
  int charge, x, y, z;
  string filename = "acetaldehyde.dat";
  ifstream input(filename);
  input >> mol.natom;
  for (int i = 0; i < mol.natom; i++) {
    input >> charge >> x >> y >> z;
    input >> charge >> x >> y >> z;
    mol.zvals.push_back(charge);
    mol.x_coords.push_back(x);
    mol.y_coords.push_back(y);
    mol.z_coords.push_back(z);}
    mol.print_geom();
  return 0;
}

//  cout << "Number of atoms: " << natom <<endl;
  //cout << "input cartesian coordinates \n";
 // for (int i=0; i<natom; i++)
   // printf("%d %20.12f %20.12f %20.12f\n", (int) zval[i], x[i], y[i], z[i]);