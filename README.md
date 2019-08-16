# Doslib

a python2 code for plotting VASP DOSCAR output

## One-liner vasp_dos

installation: ln doslib/vasp_dos ~/bin/vasp_dos

usage: vasp_dos [-h] [-s S] [-p P] [-d D] [-f F] [--write_dos0] [--write_pdos]
                [--center_ef] [--peratom] [--perspecies] [-name NAME]
                [-path PATH] [-doscar DOSCAR] [-poscar POSCAR]
                [-zoomin ZOOMIN ZOOMIN] [-energyshift ENERGYSHIFT]
                [-atom_ntype [ATOM_NTYPE [ATOM_NTYPE ...]]]

A quick wrapper of the doslib python code. It only enables the most
common/basic functions.

optional arguments:

 *  -h, --help            show this help message and exit
 *  -s S                  analyze the s orbital from element n. ie: -s 1
 *  -p P                  analyze the s orbital from element n. ie: -p 1
 *  -d D                  analyze the s orbital from element n. ie: -d 1
 *  -f F                  analyze the s orbital from element n. ie: -f 1
 *  --write_dos0          write the matrix of total dos
 *  --write_pdos          write the matrix of pdos
 *  --center_ef           shift fermi level to zero
 *  --peratom             wirte_peratom information
 *  --perspecies          plot dos for each species
 *  -name NAME            name of the job
 *  -path PATH            pathway of the doscar, poscar
 *  -doscar DOSCAR        filename of the doscar
 *  -poscar POSCAR        filename of the poscar
 *  -zoomin ZOOMIN ZOOMIN only plot the dos from value1 to value2
 *  -energyshift ENERGYSHIFT shift the total energy by the value
 *  -atom_ntype [ATOM_NTYPE [ATOM_NTYPE ...]] declare the number of different atom types type1 type2 type3


details of the control setting can be found in the doslib/control.py. To use
the code as library, see dos-lib/example folder.

## call it as a library

more examples in the example directory

## about

This is a ridiculous code that I wrote to analyze doscar. I wanted it to be able to do some simple plotting and analysis. Originally, it is a very improvised code. Somehow, I was mad at the ugly presentation and decided to rewrite it. So it looks slightly neater. I'm not sure whether the current version is better to read, since I'm still learning pythony writing style. But at least, it is an attempt.
