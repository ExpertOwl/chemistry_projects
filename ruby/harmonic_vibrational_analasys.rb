require 'matrix'
require './molecule'

mol = Molecule.new
mol.read_geom('./input/h2o_geom.txt')
mol.read_hessian('./input/h2o_hessian.txt')
p mol.hessian
