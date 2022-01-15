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
load_borns(std::vector<FLOAT_PREC> & borns, const std::string & filename ) {
  std::ifstream fin(filename);
  FLOAT_PREC element;
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
  unsigned adof = natoms * 3  // atomic degree of freedom
  std::vector<FLOAT_PREC> _dipole( odim );
  std::vector<FLOAT_PREC> _force ( odim * adof);
  std::vector<FLOAT_PREC> _virial( odim * 9);
  // iterate over all bonds, accumulate dipole and force
  for(unsigned i = 0; i < nbonds; ++i) {
    int atom1 = bonds[2*i]
    int atom2 = bonds[2*i+1]
    std::vector<FLOAT_PREC> r1 = getPosition(atom1);
    std::vector<FLOAT_PREC> r2 = getPosition(atom2);
    std::vector<FLOAT_PREC> r12(3);
    if(!nopbc) {
        r12 = pbcDistance(r1,r2);
      } else{
        for(unsigned j = 0; j < 3; ++j) {
          r12[j] = r2[j] - r1[j];
          };
      }
    std::vector<FLOAT_PREC> born_tensor(borns.begin()+ odim*3*i, borns.begin()+odim*3*(i+1));
    for(unsigned j = 0; j < odim; ++j) {
      for(unsigned k = 0; k < 3; ++k) {
        _dipole[j] += born_tensor[3*j+k] * r12[k];  //  dipole[j] = born@(r2-r1)
        _force [j*adof + 3*atom2+k] -= born_tensor[3*j+k] //  F[j, atom2, k] += - d(dipole[j])/d(r2[k]) 
        _force [j*adof + 3*atom1+k] += born_tensor[3*j+k] //  F[j, atom1, k] += - d(dipole[j])/d(r1[k]) 
        for(unsigned l = 0; l < 3; ++l) {
          _virial[j, k, l] = 
          }
        }
      }
    }
  // get  dipole
  for (unsigned j = 0; j < ODIM; ++j) {
    getPntrToComponent(j)->set(_dipole[j]);
  };
  // get  force
  for (unsigned j = 0; j < ODIM; ++j) {
    for(unsigned n = 0; n < natoms; ++n) {
      setAtomsDerivatives(
        getPntrToComponent(j), 
        n,
        Vector(_force[j*adof + 3*n + 0],
               _force[j*adof + 3*n + 1], 
               _force[j*adof + 3*n + 2]) 
      );
    };
  };
  // get stress
  for (unsigned j = 0; j < ODIM; ++j) {
    setBoxDerivatives(
      getPntrToComponent(j),
      Tensor(_virial[j])
  };

}


}
}
