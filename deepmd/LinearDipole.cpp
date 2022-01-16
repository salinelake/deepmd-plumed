#include "LinearDipole.h"

namespace PLMD {
namespace colvar {


PLUMED_REGISTER_ACTION(LinearDipole,"LinearDipole")
void
LinearDipole::
registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  keys.add("atoms","ATOMS","the group of atoms we are calculating the CV for (defaults to the whole system)");
  // keys.addFlag("COMPONENTS",false,"calculate the x, y and z components of the dipole separately and store them as label.x, label.y and label.z");
  keys.add("compulsory","BONDS","bonds.raw","the bonds topology");
  keys.add("compulsory","BORNS","borns.raw","the born tensor associated to each bond");
  for (unsigned k = 0; k < odim; ++k) {
    keys.addOutputComponent(cpnts[k],"COMPONENTS","the "+cpnts[k]+"-component of the dipole");
  }
}

void
LinearDipole::
load_bonds(std::vector<int> & bonds, const std::string & filename ) {
  std::ifstream fin(filename);
  int element;
  bonds.clear();
  while (fin >> element) {
    bonds.push_back(element);
  }
}

void
LinearDipole::
load_borns(std::vector<double> & borns, const std::string & filename ) {
  std::ifstream fin(filename);
  double element;
  borns.clear();
  while (fin >> element) {
    borns.push_back(element);
  }
}

LinearDipole::
LinearDipole(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  nopbc(false)
  components(true)
{
  // parseFlag("COMPONENTS",components);
  parseAtomList("ATOMS",atoms);
  parseFlag("NOPBC",nopbc);
  std::string bonds_file;
  parse("BONDS", bonds_file);
  std::string borns_file;
  parse("BORNS", borns_file);

  checkRead();
  // default use all atoms; otherwise warn the user
  unsigned nmdatoms = plumed.getAtoms().getNatoms();
  if (atoms.size() == 0) {
    atoms.resize(nmdatoms);
    for(unsigned i = 0; i < nmdatoms; ++i) {
      atoms[i].setIndex(i);
    }
    log.printf("  of all %u atoms in the system\n",static_cast<unsigned>(atoms.size()));
  } else {
    log.printf("  of %u atoms\n",static_cast<unsigned>(atoms.size()));
  }
  for(unsigned i = 0; i < atoms.size(); ++i) {
    log.printf("  %d", atoms[i].serial());
  }
  log.printf("  \n");
  if (atoms.size() != nmdatoms) { 
    log.printf("  # of atoms provided: %d != # of atoms in the system: %d\n", atoms.size(), nmdatoms); 
    log.printf("  Please make sure you know what you are doing!\n");
  }
  // print periodic boundary condition
  if(nopbc) { log.printf("  without periodic boundary conditions\n"); }
  else      { log.printf("  using periodic boundary conditions\n"); }

  // load bonds topology
  load_bonds(bonds, bonds_file);
  log.printf("  assign bonds to atoms\n");
  for(unsigned i = 0; i < bonds.size()/2; ++i) {
    log.printf(" (%d, %d) ", bonds[2*i], bonds[2*i+1]);
  }
  log.printf("  \n");

  // load bonds topology
  load_bonds(borns, borns_file);
  log.printf("  assign Born tensor to bonds\n");
  // for(unsigned i = 0; i < borns.size(); ++i) {
  //   log.printf("  %d", borns[i]);
  // }
  log.printf("  \n");

  // by default all components require derivarives
  addComponentWithDerivatives("x"); componentIsNotPeriodic("x");
  addComponentWithDerivatives("y"); componentIsNotPeriodic("y");
  addComponentWithDerivatives("z"); componentIsNotPeriodic("z");

  requestAtoms(atoms);
}

// calculator
void LinearDipole::calculate()
{
  unsigned natoms = getNumberOfAtoms();
  unsigned nbonds = bonds.size/2;
  unsigned adof = natoms * 3  // atomic degrees of freedom
  std::vector<double> _dipole( odim, 0);
  std::vector<double> _force ( odim * adof, 0);
  std::vector<double> _virial( odim * 9, 0);
  // iterate over all bonds, accumulate dipole, force and virial
  for(unsigned i = 0; i < nbonds; ++i) {
    int atom1 = bonds[2*i];
    int atom2 = bonds[2*i+1];
    std::vector<double> r1 = getPosition(atom1);
    std::vector<double> r2 = getPosition(atom2);
    std::vector<double> r12(3);
    if(!nopbc) {
        r12 = pbcDistance(r1,r2);
      } else{
        for(unsigned j = 0; j < 3; ++j) {
          r12[j] = r2[j] - r1[j];
          };
      }
    for(unsigned j = 0; j < odim; ++j) {
      std::vector<double> born_vector(borns.begin() + odim*3*i + 3*j, borns.begin() + odim*3*i + 3*(j+1)); // borns[idx_bond, j, *]
      // this bond's contribution to the force = - d(dipole)/dr
      for(unsigned k = 0; k < 3; ++k) {
        _dipole[j] += born_vector[k] * r12[k];  //  dipole[j] = born@(r2-r1)
        _force [j*adof + 3*atom2+k] -= born_vector[k]; //  F[j, atom2, k] += - d(dipole[j])/d(r2[k]) 
        _force [j*adof + 3*atom1+k] += born_vector[k]; //  F[j, atom1, k] += - d(dipole[j])/d(r1[k]) 
        }
      // this bond's contribution to the virial v = r x F
      for(unsigned k = 0; k < 3; ++k) {
        for(unsigned l = 0; l < 3; ++l) {
          _virial[j*9 + k*3 + l] -= r12[k] * born_vector[l];    // virial[idx_odim,k,l] += (r2-r1)[k] * borns[idx_bond,  idx_odim, l]
          }
      }
    }
  }
  // get global dipole
  for (unsigned j = 0; j < ODIM; ++j) {
    getPntrToComponent(j)->set(_dipole[j]);
  };
  // get atomic detivative = minus force
  for (unsigned j = 0; j < ODIM; ++j) {
    for(unsigned n = 0; n < natoms; ++n) {
      setAtomsDerivatives(
        getPntrToComponent(j), n,
        - Vector(_force[j*adof + 3*n + 0],
                 _force[j*adof + 3*n + 1], 
                 _force[j*adof + 3*n + 2])
      );
    };
  };
  // get box derivative = minus virial
  for (unsigned j = 0; j < ODIM; ++j) {
    
    setBoxDerivatives(
      getPntrToComponent(j), // what is plumed's convention?
      - Tensor(_virial[j*9 +0], _virial[j*9 +1], _virial[j*9 +2],
               _virial[j*9 +3], _virial[j*9 +4], _virial[j*9 +5],
               _virial[j*9 +6], _virial[j*9 +7], _virial[j*9 +8])
    );
  };
}


}
}
