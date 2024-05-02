# frozen_string_literal: true

require 'csv'
require 'matrix'
require './molecule'

# filename = 'h2o_geom.txt'
filename = 'acetaldehyde.dat'

atomic_weights = []
CSV.foreach('../shared_data/atomic_weights.csv') { |row| atomic_weights << row[0].to_f }
ATOMIC_WEIGHTS = atomic_weights.freeze

def close?(floata, floatb, tol = 1e-4)
  (floata - floatb).abs < tol
end

def classify_rotor_type(inertia_eigenvectors)
  i_a, i_b, i_c = inertia_eigenvectors

  if close?(i_a, 0)
    :linear
  elsif close?(i_a, i_c)
    :spherical_top
  elsif !(close?(i_a, i_b) && close?(i_b, i_c))
    :asymeteric_top
  else
    close?(i_a, i_b) ? :oblate_top : prolate_top
  end
end

mol = Molecule.new
mol.read_geom_from_file(filename)

# puts mol
# mol.print_atomic_distances
# mol.print_bond_angles
# mol.print_torsional_angles
# puts mol.center_of_mass

mol.center_cm_on_origin
p classify_rotor_type(mol.moments_of_inertia.eigen.eigenvalues)

# TODO: report units for moments of inertia output
# TODO: Calculate rotational constants
# TODO:# Move all printing (and maybe some other) functions from mol into somewhere else, possibly a GeomAnalysis class?
