# frozen_string_literal: true

require 'csv'
require 'matrix'
require './molecule'

# Constants
PI = Math::PI
SPEED_OF_LIGHT = 2.99792458E10 # cm s-1
PLANCK = 6.6261E-27 # cm^2 g s^-1
ROTATIONAL_CONSTANT = PLANCK / (8 * (PI**2) * SPEED_OF_LIGHT) # g cm
# Converstions
RAD_TO_DEG = 180 / PI
BOHR_TO_ANGSTROM = 0.529177249
AMU_TO_G = 1.66054e-24
ANGSTROM_TO_CM = 1e-8
AMU_BOHR2_TO_G_CM = AMU_TO_G * BOHR_TO_ANGSTROM**2 * ANGSTROM_TO_CM**2
WAVENUMBER_TO_MHZ = 2.99724589E4
# Constants from files


filename = './input/acetaldehyde.dat'

def print_atomic_distances(mol)
  puts "Atomic Distances for #{mol} \n"
  (0..mol.num_atom - 1).to_a.combination(2).each do |(i, j)|
    puts "#{i} #{j}  #{distance(mol[i], mol[j])}" if i != j
  end
end

def distance(atom_i, atom_j)
  (atom_j - atom_i).magnitude
end

def print_bond_angles(mol)
  puts "Bond Angles for #{mol} \n"
  triplets = mol.combinations(3)
  triplets.each do |(i, j, k)|
    puts "#{i} #{j} #{k}  #{angle(mol[i], mol[j], mol[k]) * RAD_TO_DEG}"
  end
end

def angle(atom_i, atom_j, atom_k)
  cos_phi_ijk = unit_vector(atom_j, atom_i).inner_product(unit_vector(atom_j, atom_k))
  Math.acos(cos_phi_ijk)
end

def print_torsional_angles(mol)
  puts "Torsional Angles for #{mol} \n"
  quads = mol.combinations(4)
  quads.each do |(i, k, j, l)|
    puts "#{i} #{j} #{k} #{l}  #{torsion_angle(mol[i], mol[j], mol[k], mol[l]) * RAD_TO_DEG}"
  end
end

def torsion_angle(atom_i, atom_j, atom_k, atom_l)
  e_ij = unit_vector(atom_i, atom_j)
  e_jk = unit_vector(atom_j, atom_k)
  e_kl = unit_vector(atom_k, atom_l)
  # cos(t) = ((e_ij X e_jk) dot (e_jk X e_kl)) / (sin(theta_ijk) sin(theta_jkl))
  numerator = e_ij.cross_product(e_jk).inner_product(e_jk.cross_product(e_kl))
  denominator = Math.sin(angle(atom_i, atom_j, atom_k)) * Math.sin(angle(atom_j, atom_k, atom_l))
  cos_tau = numerator / denominator
  Math.acos(cos_tau.clamp(-1, 1))
end

def unit_vector(atom_i, atom_j)
  (atom_i - atom_j).normalize
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

def build_distance_matrix(mol)
  mol.distance_matrix = Matrix.build(mol.geom.row_count) { |i, j| distance(self[i], self[j]) }
end

def center_of_mass(mol)
  cm = Vector.zero(3)
  mol.atomic_weights.each_with_index do |weight, index|
    cm += mol[index] * weight
  end
  cm /= mol.atomic_weights.sum
end

def moments_of_inertia(mol)
  intertia_tensor = Matrix.zero(3)
  mol.geom.row_vectors.each_with_index do |r, i|
    # I = m(r.t r) I - r r.T
    intertia_tensor += ((r.inner_product(r) * Matrix.identity(3)) - (r * r.to_matrix.t)) * mol.atomic_weights[i]
  end
  intertia_tensor.eigensystem.eigenvalues
end

def print_rotor_type(inertia_eigenvalues)
  puts "Inertia eigenvalues (amu bohr^2) \n #{inertia_eigenvalues}"
  puts "Inertia eigenvalues (amu Å^2) \n #{convert(inertia_eigenvalues, BOHR_TO_ANGSTROM**2)}"
  puts "Inertia eigenvalues (g cm^2)\n #{convert(inertia_eigenvalues, AMU_BOHR2_TO_G_CM)}"
  puts "Molecule is #{rotor_type(*inertia_eigenvalues)}"
end

def convert(array_of_values, constant)
  array_of_values.map { |i| i * constant }
end

def rotor_type(i_a, i_b, i_c)
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

def close?(floata, floatb, tol = 1e-4)
  (floata - floatb).abs < tol
end

def print_rotational_constants(inertia_eigenvalues)
  rot_cons = inertia_eigenvalues.map { |val| ROTATIONAL_CONSTANT / val }
  rot_cons = convert(rot_cons, 1 / AMU_BOHR2_TO_G_CM)
  puts "Rotational constants (cm^-1) \n #{rot_cons}"
  rot_cons = convert(rot_cons, *WAVENUMBER_TO_MHZ)
  puts "Rotational constants (MHz) \n #{rot_cons}"
end

mol = Molecule.new
mol.read_geom(filename)
print_atomic_distances(mol)
print_bond_angles(mol)
print_torsional_angles(mol)
mol.translate(-center_of_mass(mol))
inertia_eigenvalues = moments_of_inertia(mol)
print_rotor_type(inertia_eigenvalues)
print_rotational_constants(inertia_eigenvalues)
