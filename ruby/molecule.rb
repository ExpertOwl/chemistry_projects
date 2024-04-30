# frozen_string_literal: true

require 'matrix'

RRAD_TO_DEG = 180 / Math::PI
# Class represemts geometry information for molecule as well as transformations of geometric coordinates
class Molecule
  # Constants
  # attr names
  attr_reader :num_atom, :list_of_atoms, :distance_matrix
  attr_accessor :geom, :atomic_weights

  def to_s
    @geom
  end

  def [](index)
    @geom.row(index)
  end

  def translate(xyz)
    @geom += Matrix.rows([xyz] * num_atom)
  end

  def add_atom(atom_xyz)
    @num_atom += 1
    @list_of_atoms << atom_xyz[0]
    @geom << atom_xyz[1..3]
  end

  def initialize
    @list_of_atoms = []
    @num_atom = 0
    @geom = []
    @atomic_weights = []
  end

  def read_geom_from_file(filename)
    f = File.new(filename)
    num_atom = f.readline.to_i
    num_atom.times do |_i|
      atom_xyz = f.readline.split.map!(&:to_f)
      add_atom(atom_xyz)
    end
    self.atomic_weights = Vector[*@list_of_atoms.map { |i| ATOMIC_WEIGHTS[i] }]
    self.geom = Matrix.rows(@geom)
  end

  def build_distance_matrix
    @distance_matrix = Matrix.build(@geom.row_count) { |i, j| distance(self[i], self[j]) }
  end

  def center_of_mass
    cm = Vector.zero(3)
    @atomic_weights.each_with_index do |weight, index|
      cm += self[index] * weight
    end
    cm /= atomic_weights.sum
  end

  def print_atomic_distances
    (0..@num_atom - 1).to_a.combination(2).each do |(i, j)|
      puts "#{i} #{j}  #{distance(self[i], self[j])}" if i != j
    end
  end

  def print_bond_angles
    triplets = generate_index_combinations(3)
    triplets.each do |(i, j, k)|
      puts "#{i} #{j} #{k}  #{bond_angle(self[i], self[j], self[k]) * RAD_TO_DEG}"
    end
  end

  def print_torsional_angles
    quads = generate_index_combinations(4)
    quads.each do |(i, k, j, l)|
      puts "#{i} #{j} #{k} #{l}  #{torsion_angle(self[i], self[j], self[k], self[l]) * RAD_TO_DEG}"
    end
  end

  def center
    translate(-center_of_mass)
  end

  def moments_of_inertia
    intertia_tensor = Matrix.zero(3)
    @geom.row_vectors.each_with_index do |r, i|
      # I = m(r.t r) I - r r.T
      intertia_tensor += ((r.inner_product(r) * Matrix.identity(3)) - (r * r.to_matrix.t)) * atomic_weights[i]
    end
    p intertia_tensor
  end

  private

  def unit_vector(atom_i, atom_j)
    (atom_i - atom_j).normalize
  end

  def distance(atom_i, atom_j)
    (atom_j - atom_i).magnitude
    # Math.sqrt((atom_i[0] - atom_j[0])**2 + (atom_i[1] - atom_j[1])**2 + (atom_i[2] - atom_j[2])**2)
  end

  def generate_index_combinations(tuple_length)
    (0..@num_atom - 1).to_a.combination(tuple_length)
  end

  def filter_index_combinations(combinations, max_distance)
    combinations.select do |tuple|
      tuple.uniq &&
        tuple.combination(2).map { |i, j| distance(self[i], self[j]) < max_distance }.all?
    end
  end

  def bond_angle(atom_i, atom_j, atom_k)
    cos_phi_ijk = unit_vector(atom_j, atom_i).inner_product(unit_vector(atom_j, atom_k))
    Math.acos(cos_phi_ijk)
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
    Math.acos(cos_tau.clamp(-1, 1))
  end
end

# TODO: refactor this class somehow
