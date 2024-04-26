# frozen_string_literal: true

require 'matrix'
require 'csv'

RAD_TO_DEG = 180 / Math::PI
ATOMIC_WEIGHTS = []
CSV.foreach('../shared_data/atomic_weights.csv') { |row| ATOMIC_WEIGHTS << row[0].to_f }
ATOMIC_WEIGHTS.freeze
# filename = 'h2o_geom.txt'
filename = 'acetaldehyde.dat'

def read_geom(filename)
  f = File.new(filename)
  atomic_number = []
  geom = []
  num_atoms = f.readline.to_i
  num_atoms.times do
    atom_xyz = f.readline.split.map!(&:to_f)
    atomic_number << atom_xyz[0]
    geom << atom_xyz[1..3]
  end
  geom = Matrix[*geom]
  [atomic_number, geom]
end

def distance(atom_i, atom_j)
  Math.sqrt((atom_i[0] - atom_j[0])**2 + (atom_i[1] - atom_j[1])**2 + (atom_i[2] - atom_j[2])**2)
end

def print_d_matrix(distance_matrix)
  p 'Interatomic Distances (bohr)'
  distance_matrix.each_with_index do |_d1, i|
    (0..i).each do |j|
      next if i == j

      pp [i, j, distance_matrix[i][j]]
    end
  end
end

def distance_matrix(geom)
  distance_matrix = []
  geom.each do |atom1|
    row = []
    geom do |atom2|
      row << distance(atom1, atom2)
    end
    distance_matrix << row
  end
  distance_matrix
end

def bond_angle(atom_i, atom_j, atom_k)
  cos_phi_ijk = unit_vector(atom_j, atom_i).inner_product(unit_vector(atom_j, atom_k))
  Math.acos(cos_phi_ijk)
end

def unit_vector(atom_i, atom_j)
(Vector[*atom_j] - Vector[*atom_i]).normalize
end

def out_of_plane_angle(atom_i, atom_j, atom_k, atom_l)
  e_kj = unit_vector(atom_k, atom_j)
  e_kl = unit_vector(atom_k, atom_l)
  e_ki = unit_vector(atom_k, atom_i)
  sin_phi_jkl = Math.sin(bond_angle(atom_j, atom_k, atom_l))
  sin_theta = e_kj.cross_product(e_kl) / sin_phi_jkl
  sin_theta = sin_theta.inner_product(e_ki)
  sin_theta.clamp(-1, 1)
  Math.asin(sin_theta)
end

def torsion_angle(atom_i, atom_j, atom_k, atom_l)
  e_ij = unit_vector(atom_i, atom_j)
  e_jk = unit_vector(atom_j, atom_k)
  e_kl = unit_vector(atom_k, atom_l)
  numerator = e_ij.cross_product(e_jk).inner_product(e_jk.cross_product(e_kl))
  denominator = Math.sin(bond_angle(atom_i, atom_j, atom_k)) * Math.sin(bond_angle(atom_j, atom_k, atom_l))
  cos_tau = numerator / denominator
  cos_tau.clamp(-1, 1)
  Math.acos(cos_tau)
end

_atomic_number, geom = read_geom(filename)


p (ATOMIC_WEIGHTS[1] - 1.008).abs < Float::EPSILON
p (distance(geom.row(0), geom.row(1)) - 2.845112131228).abs < Float::EPSILON
p (bond_angle(geom.row(0), geom.row(1), geom.row(2)) * RAD_TO_DEG - 124.26830826072018).abs < Float::EPSILON
p (torsion_angle(geom.row(6), geom.row(5), geom.row(4), geom.row(0)) * RAD_TO_DEG - 36.36679948827571).abs < Float::EPSILON

#TODO: Change distance matrix method to use an actual matrix
#Calculate center of mass then continue from lesson 1
#Create moleulce class to make getting atom attributes easier
