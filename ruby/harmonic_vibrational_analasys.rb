# frozen_string_literal: true

require 'matrix'
require './molecule'
FREQ_CONSTANT = 5140.484532

def print_hessian_eigen(mol)
  puts "Hessian eigenvalues (hartree/amu-bohr^2) for #{mol} \n"
  puts hessian_eigenval(mol)
  puts "\n"
end

def hessian_eigenval(mol)
  mol.hessian.eigensystem.eigenvalues
end

def print_freq(mol)
  puts "Vibrational frequencies for #{mol} (cm^-1) \n"
  hessian_eigenval(mol).each do |eigan|
    if  eigan.negative?
    puts Complex(0, (Math.sqrt(eigan.abs) * FREQ_CONSTANT).round(-3))
    else
    puts (Math.sqrt(eigan) * FREQ_CONSTANT).round(3)
    end
  end
  puts "\n"
end

mol = Molecule.new
mol.read_geom('./input/benzene_geom.txt')
mol.read_hessian('./input/benzene_hessian.txt')
print_hessian_eigen(mol)
print_freq(mol)
