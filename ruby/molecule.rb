# frozen_string_literal: true

require 'csv'
require 'matrix'
ATOMIC_WEIGHTS = CSV.read('../shared_data/atomic_weights.csv').map { |weight| weight[0].to_f }
# Class represemts geometry information for molecule as well as transformations of geometric coordinates
class Molecule
  attr_reader :num_atom, :list_of_atoms, :distance_matrix
  attr_accessor :geom, :atomic_weights, :hessian

  def initialize
    @num_atom = 0
    @list_of_atoms = []
    @geom = []
    @atomic_weights = []
  end

  def read_geom(filename)
    f = File.new(filename)
    num_atom = f.readline.to_i
    num_atom.times do |_i|
      atom_xyz = f.readline.split.map!(&:to_f)
      add_atom(atom_xyz)
    end
    self.atomic_weights = Vector[*@list_of_atoms.map { |i| ATOMIC_WEIGHTS[i] }]
    self.geom = Matrix.rows(@geom)
  end

  def read_hessian(filename)
    f = File.new(filename)
    puts 'Wrong number of atoms' unless @num_atom == f.readline.to_i
    input_data = lines_to_a(f)
    @hessian = Matrix.build(@num_atom * 3, @num_atom * 3) { |i, j| input_data[3 * i + j / 3][j % 3] }
  end

  def lines_to_a(filestream)
    input_data = filestream.read.split("\n")
    input_data.map!(&:split)
    input_data.map! { |row| row.map!(&:to_f) }
  end

  def to_s
    @geom
  end

  def [](index)
    @geom.row(index)
  end

  def add_atom(atom_xyz)
    @num_atom += 1
    @list_of_atoms << atom_xyz[0]
    @geom << atom_xyz[1..3]
  end

  def combinations(tuple_length, max_distance = 4)
    combinations = (0..@num_atom - 1).to_a.combination(tuple_length)
    combinations.select do |tuple|
      tuple.uniq &&
        tuple.combination(2).map { |i, j| distance(self[i], self[j]) < max_distance }.all?
    end
  end

  def center_cm_on_origin
    translate(-center_of_mass)
  end

  def translate(xyz)
    @geom += Matrix.rows([xyz] * num_atom)
  end
end
