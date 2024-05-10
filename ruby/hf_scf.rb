# frozen_string_literal: true

require 'matrix'
require './molecule'

mol = Molecule.new
mol.read_geom('./input/H2O/STO-3G/geom.dat')

def two_center_quantity_from_file(filename)
  f = File.new(filename)
  data = f.readlines.map do |line|
    line.split.map(&:to_f)
  end
  matr_size = data[-1][0]
  two_center_quantity = Matrix.zero(matr_size.to_i)
  data.each do |i, j, quant|
    two_center_quantity [i - 1, j - 1] = quant
    two_center_quantity [j - 1, i - 1] = quant
  end
  two_center_quantity
end

def read_2e_integrals(filename)
  f = File.new(filename)
  data = f.readlines.map do |line|
    line.split.map(&:to_f)
  end
  matr_size = data[-1][0]
  qq = [[[[0] * matr_size] * matr_size] * matr_size] * matr_size
  data.each do |u, v, l, o, quant|
    u -= 1
    v -= 1
    l -= 1
    o -= 1
    qq[u][v][l][o] = quant
    qq[v][u][l][o] = quant
    qq[u][v][o][l] = quant
    qq[v][u][o][l] = quant
    qq[l][o][u][v] = quant
    qq[o][l][u][v] = quant
    qq[l][o][v][u] = quant
    qq[o][l][v][u] = quant
  end
  qq
end

def compound(i, j)
  if i > j
    i * (i + 1) / 2 + j
  else
    j * (j + 1) / 2 + i
  end
end

mol.s = two_center_quantity_from_file('./input/H2O/STO-3G/s.dat')
mol.t = two_center_quantity_from_file('./input/H2O/STO-3G/t.dat')
mol.v = two_center_quantity_from_file('./input/H2O/STO-3G/v.dat')
mol.eri = read_2e_integrals('./input/H2O/STO-3G/eri.dat')
# p mol.eri.size
# p mol.eri
h_core = mol.t + mol.v
def printmat(matrix)
  matrix.row_vectors.each do |i|
    puts i
  end
end

v, eig, v_inv =  mol.s.eigensystem
s_half = v * (eig ** (-0.5)) * v_inv

fock = s_half.t * h_core * s_half

# c_0_prime, e, c_0_prime_inv = fock.eigensystem
# c_0 = s_half.t * c_0_prime.t
# p s_half * s_half.t
# p s_half
