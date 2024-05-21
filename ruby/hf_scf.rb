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
  res = {}
  data.each do |i, j, k, l, quant|
    i = i.to_i - 1
    j = j.to_i - 1
    k = k.to_i - 1
    l = l.to_i - 1
    res[[i,j,k,l]] = quant
    res[[i,j,l,k]] = quant
    res[[j,i,k,l]] = quant
    res[[j,i,l,k]] = quant
    res[[k,l,i,j]] = quant
    res[[l,k,i,j]] = quant
    res[[k,l,j,i]] = quant
    res[[l,k,j,i]] = quant
  end
  res
end

def read_one_electron_integrals(mol, folder)
  mol.read_geom("#{folder}/geom.dat")
  mol.e_nuc = read_enuc_from_file("#{folder}/enuc.dat")
  mol.s = two_center_quantity_from_file("#{folder}/s.dat")
  mol.t = two_center_quantity_from_file("#{folder}/t.dat")
  mol.v = two_center_quantity_from_file("#{folder}/v.dat")
end

def read_two_electron_integrals(mol,folder)
  mol.eri = read_2e_integrals("#{folder}/eri.dat")
end

def get_s_half(mol)
  v, eig, v_inv =  mol.s.eigensystem
  v * (eig**-0.5) * v_inv
end

def order_by_eigenvalues(matrix, eigenvalues)
 values = eigenvalues.max_by(eigenvalues.row_count, &:-@)
 values.map! do |i|
    eigenvalues.find_index(i)[0]
  end
  p values
  Matrix.columns(values.map { |col_index| matrix.column_vectors[col_index] })
end

def density_from_fock(mol, sort = false)
  mol.s_half = get_s_half(mol)
  fock_transformed = mol.s_half.t * mol.fock * mol.s_half

  c_0_prime, energy_eig, c_0_prime_inv = fock_transformed.eigensystem
  pp energy_eig.to_a
  c = mol.s_half * c_0_prime
  if sort
      c = order_by_eigenvalues(c, energy_eig) #Eigensystem does not guarantee order of eigenvalues, sort matrix
  end
  density_matrix = Matrix.build(c.row_count, c.column_count) do |i, j|

    num = 0
    mol.n_occ.times do |m|
      num += c[j, m] * c[i, m]
    end
    num
  end
  density_matrix
end

def electronic_energy(mol)
  e_elec = 0
    mol.density.each_with_index do |d, i, j|
      e_elec += d * (mol.fock[i, j] + mol.h_core[i, j])
    end
    e_elec
end

def fock_from_density(mol)
  nao = mol.nao
  eri = mol.eri
  h_core = mol.h_core
  density = mol.density
  fock_new = Matrix.zero(nao, nao)
  nao.times do |i|
    nao.times do |j|
      fock_new[i,j] = h_core[i,j]
      nao.times do |k|
        nao.times do |l|
          if mol.eri[[i,j,k,l]] && mol.eri[[i,k,j,l]]
            fock_new[i,j] += density[k,l] * (2 * eri[[i,j,k,l]] - eri[[i,k,j,l]])
          end
        end
      end
    end
  end
  fock_new
end

mol = Molecule.new
input_folder = "./input/H2O/STO-3G"
mol.n_occ = 5

read_one_electron_integrals(mol, input_folder)
read_two_electron_integrals(mol, input_folder)

mol.h_core = mol.t + mol.v
mol.nao = mol.v.row_count
mol.s_half = get_s_half(mol)

#Initial guess
mol.fock = mol.h_core
#SCF
puts " inital density \n"
mol.density = density_from_fock(mol, true)
mol.electronic_energy = electronic_energy(mol)
pp mol.electronic_energy
5.times do |i|
  puts "iter #{i}"
  mol.fock = fock_from_density(mol)
  mol.density = density_from_fock(mol, true)
  mol.electronic_energy = electronic_energy(mol)
  pp mol.electronic_energy

end





# TODO: Either density is being calcualated wrong or fock is wrong. Not sure. Check S&O.
