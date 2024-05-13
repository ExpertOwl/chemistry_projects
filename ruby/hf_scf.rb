# frozen_string_literal: true

require 'matrix'
require './molecule'

def read_enuc_from_file(filename)
  f = File.new(filename)
  f.readline.to_f
end

def print_matrix(matrix)
  matrix.row_vectors.each do |row|
    row_str = []
    row.each do |item|
      row_str << if item.abs < 1E-12
                   format('% .10f', 0)
                 else
                   format('% .10f', item)
                 end
    end
    puts row_str.join('  ')
  end
  puts "\n"
end

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
  res = [[[[0] * matr_size] * matr_size] * matr_size] * matr_size
  data.each do |i, j, k, l, quant|
    ii = i - 1
    jj = j - 1
    kk = k - 1
    ll = l - 1
    res[ii][jj][kk][ll] = quant
    res[ii][jj][ll][kk] = quant
    res[jj][ii][kk][ll] = quant
    res[jj][ii][ll][kk] = quant
    res[kk][ll][ii][jj] = quant
    res[ll][kk][ii][jj] = quant
    res[kk][ll][jj][ii] = quant
    res[ll][kk][jj][ii] = quant
  end
  res
end

# def read_2e_integrals(filename)
#   f = File.new(filename)
#   data = f.readlines.map do |line|
#     line.split.map(&:to_f)
#   end
#   matr_size = data[-1][0]

#   qq = {}
#   data.each do |i, j, k, l, quant|
#     i -= 1
#     j -= 1
#     k -= 1
#     l -= 1
#     ij = compound(i,j)
#     kl = compound(k,l)
#     ijkl = compound(ij,kl)
#     qq[ijkl.to_i] = quant
#     qq[ijkl.to_i] = quant
#   end
#   qq
# end

def compound(i, j)
  i >= j ? ((i * (i + 1) / 2) + j) : ((j * (j + 1) / 2) + i)
end

mol = Molecule.new
mol.read_geom('./input/H2O/STO-3G/geom.dat')
mol.e_nuc = read_enuc_from_file('./input/H2O/STO-3G/enuc.dat')
mol.s = two_center_quantity_from_file('./input/H2O/STO-3G/s.dat')
mol.t = two_center_quantity_from_file('./input/H2O/STO-3G/t.dat')
mol.v = two_center_quantity_from_file('./input/H2O/STO-3G/v.dat')
eri = read_2e_integrals('./input/H2O/STO-3G/eri.dat')
# p mol.eri.size
# p mol.eri
h_core = mol.t + mol.v

v, eig, v_inv =  mol.s.eigensystem
s_half = v * (eig**-0.5) * v_inv
fock = s_half.t * h_core * s_half
c_0_prime, e, = fock.eigensystem
c = s_half * c_0_prime

a = e.max_by(e.row_count, &:-@)

a.map! do |i|
  e.find_index(i)[0]
end

c_ordered = Matrix.columns(a.map { |col_index| c.column_vectors[col_index] })

density_matrix = Matrix.build(c_ordered.row_count, c_ordered.column_count) do |i, j|
  num = 0
  5.times do |m| # Need to get orbital occupiation somehow instead of using 5!!!
    num += c_ordered[i, m] * c_ordered[j, m]
  end
  num
end
e_elec = 0

# initial guess
density_matrix.each_with_index do |d, i, j|
  e_elec += d * (h_core[i, j] + h_core[i, j])
end
puts e_elec

# build new fock
p compound(compound(0, 0), compound(0, 3))
nao = fock.row_count
fock_new = Matrix.zero(fock.row_count, fock.column_count)
2.times do |i|
  fock_new[i, 0] = h_core[i, 0]
  nao.times do |k|
    nao.times do |l|
      # fock_new[i,0] += (density_matrix[k,l] * eri[i][0][k][l] - 0.5 * eri[i][k][0][l])
      puts eri[i][0][k][l] == eri[i][k][0][l]
    end
  end
end
# number += value * (2*mol.eri[ijkl] - mol.eri[ikjl])
p fock_new

8.times do |q|
  p eri[q][q][q][q]
end
p

# TODO: ERI is being read in wrong. need to figure out why
get orbital occupiation somehow instead of using hard coded number
