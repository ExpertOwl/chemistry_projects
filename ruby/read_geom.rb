# frozen_string_literal: true
require 'csv'
require './molecule'

# filename = 'h2o_geom.txt'
filename = 'acetaldehyde.dat'

atomic_weights = []
CSV.foreach('../shared_data/atomic_weights.csv') { |row| atomic_weights << row[0].to_f }
ATOMIC_WEIGHTS = atomic_weights.freeze

mol = Molecule.new
mol.read_geom_from_file(filename)

# puts mol
# mol.print_atomic_distances
# mol.print_bond_angles
# mol.print_torsional_angles
# puts mol.center_of_mass
mol.center
mol.moments_of_inertia

# TODO:
# Calculate center of mass then continue from lesson 1
# Move geometeric functions to a molecule class
