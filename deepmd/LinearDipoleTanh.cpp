#include "LinearDipoleTanh.h"

namespace PLMD {
namespace colvar {


PLUMED_REGISTER_ACTION(LinearDipoleTanh,"LINEARDIPOLETANH")
void
LinearDipoleTanh::
registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  keys.add("atoms","ATOMS","the group of atoms we are calculating the CV for (defaults to the whole system)");
  // keys.addFlag("COMPONENTS",false,"calculate the x, y and z components of the dipole separately and store them as label.x, label.y and label.z");
  keys.add("compulsory","BONDS","bonds.raw","the bonds topology");
  keys.add("compulsory","BORNS","borns.raw","the born tensor associated to each bond");
  keys.add("compulsory","RESCALE", "1.0", "rescaling coefficient for the dipole");
  for (unsigned k = 0; k < odim; ++k) {
    keys.addOutputComponent(cpnts[k],"COMPONENTS","the "+cpnts[k]+"-component of the dipole");
    keys.addOutputComponent("e"+cpnts[k],"COMPONENTS","the "+cpnts[k]+"-component of the effective (tanh) dipole");
    keys.addOutputComponent("n"+cpnts[k],"COMPONENTS","the number of unitcell dipole that has positive "+cpnts[k]+"-component");
  }
}

void
LinearDipoleTanh::
load_bonds(std::vector<int> & bonds, const std::string & filename ) {
  std::ifstream fin(filename);
  int element;
  bonds.clear();
  while (fin >> element) {
    bonds.push_back(element);
  }
}

void
LinearDipoleTanh::
load_borns(std::vector<double> & borns, const std::string & filename ) {
  std::ifstream fin(filename);
  double element;
  borns.clear();
  while (fin >> element) {
    borns.push_back(element);
  }
}


LinearDipoleTanh::
LinearDipoleTanh(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  nopbc(false),
  components(true)
{

  parseAtomList("ATOMS",atoms);
  parseFlag("NOPBC",nopbc);
  std::string bonds_file;
  parse("BONDS", bonds_file);
  std::string borns_file;
  parse("BORNS", borns_file);
  parse("RESCALE", rescale);

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
  // for(unsigned i = 0; i < atoms.size(); ++i) {
  //   log.printf("  %d", atoms[i].serial());
  // }
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
  nbonds = bonds.size()/2;
  // log.printf("  assign bonds to atoms\n");
  // for(unsigned i = 0; i < bonds.size()/2; ++i) {
  //   log.printf(" (%d, %d) ", bonds[2*i], bonds[2*i+1]);
  // }
  // log.printf("  \n");

  // load born tensor
  load_borns(borns, borns_file);
  if (nbonds != borns.size()/ odim /3){
    log.printf("#Bonds and #Borns tensor do not match");
  }
  log.printf("  assign Born tensor to bonds\n");
  for(unsigned i = 0; i < borns.size()/9; ++i) {
    log.printf(" (%d, %d) :", bonds[2*i], bonds[2*i+1]);
    for(unsigned j = 0; j < 9; ++j){ log.printf(" %f ", borns[9*i+j]); };
    log.printf("  \n ");
    if (i>10){
      log.printf(" ...... \n ");
      break;
    }
  }
  log.printf(" Rescaling coefficient for modified dipole = %f \n", rescale);

  addComponentWithDerivatives("x"); componentIsNotPeriodic("x");
  addComponentWithDerivatives("y"); componentIsNotPeriodic("y");
  addComponentWithDerivatives("z"); componentIsNotPeriodic("z");
  addComponentWithDerivatives("ex"); componentIsNotPeriodic("ex");
  addComponentWithDerivatives("ey"); componentIsNotPeriodic("ey");
  addComponentWithDerivatives("ez"); componentIsNotPeriodic("ez");
  addComponent("nx"); componentIsNotPeriodic("nx");
  addComponent("ny"); componentIsNotPeriodic("ny");
  addComponent("nz"); componentIsNotPeriodic("nz");
  requestAtoms(atoms);
  natoms = getNumberOfAtoms();
  adof = natoms * 3;  // atomic degrees of freedom
  unsigned nt=OpenMP::getNumThreads();
  log.printf(" Using %d threads \n", nt);
}
  
Vector 
LinearDipoleTanh::
_get_distance(const std::vector< Vector >  & positions, int & a1, int & a2)
{
  Vector r1 = positions[a1];
  Vector r2 = positions[a2];
  // Vector r1 = getPosition(a1) ;
  // Vector r2 = getPosition(a2) ;
  Vector r12;
  if(!nopbc) {
      r12 = pbcDistance(r1,r2);
    } else{
      r12 = r2 - r1;
    }
  return r12;
}
void
LinearDipoleTanh::
_feed_plumed(std::vector<double> & dipole, 
             std::vector<double> & force, 
             std::vector<double> & virial, 
             std::vector<double> & eff_dipole, 
             std::vector<double> & eff_force, 
             std::vector<double> & eff_virial, 
             std::vector<long> & dipole_upcount)
{
  // unsigned natoms = getNumberOfAtoms();
  for (unsigned j = 0; j < odim; ++j) {
    getPntrToComponent(cpnts[j])->set(dipole[j]);
    getPntrToComponent("e"+cpnts[j])->set(eff_dipole[j]);
    getPntrToComponent("n"+cpnts[j])->set(dipole_upcount[j]);
  };
  // get atomic detivative = minus force
  for (unsigned j = 0; j < odim; ++j) {
    for(unsigned n = 0; n < natoms; ++n) {
      setAtomsDerivatives(
        getPntrToComponent(cpnts[j]), n,
        - Vector(force[j*adof + 3*n + 0],
                 force[j*adof + 3*n + 1], 
                 force[j*adof + 3*n + 2])
                 );
      setAtomsDerivatives(
        getPntrToComponent("e"+cpnts[j]), n,
        - Vector(eff_force[j*adof + 3*n + 0],
                 eff_force[j*adof + 3*n + 1], 
                 eff_force[j*adof + 3*n + 2])
                 );
    };
  };
  // get box derivative = minus virial
  for (unsigned j = 0; j < odim; ++j) {
    setBoxDerivatives(
      getPntrToComponent(cpnts[j]), 
      - Tensor(virial[j*9 +0], virial[j*9 +1], virial[j*9 +2],
               virial[j*9 +3], virial[j*9 +4], virial[j*9 +5],
               virial[j*9 +6], virial[j*9 +7], virial[j*9 +8])
    );
    setBoxDerivatives(
      getPntrToComponent("e"+cpnts[j]), 
      - Tensor(eff_virial[j*9 +0], eff_virial[j*9 +1], eff_virial[j*9 +2],
               eff_virial[j*9 +3], eff_virial[j*9 +4], eff_virial[j*9 +5],
               eff_virial[j*9 +6], eff_virial[j*9 +7], eff_virial[j*9 +8])
    );
  };
}

void 
LinearDipoleTanh::
calculate()
{
  unsigned nt=OpenMP::getNumThreads();
  // log.printf(" Using %d threads \n", nt);
  const std::vector< Vector > positions = getPositions();
  // natoms = getNumberOfAtoms();
  std::vector<double> dipole( odim, 0);
  std::vector<double> force ( odim * adof, 0);
  std::vector<double> virial( odim * 9, 0);
  std::vector<double> unitcell_dipole(odim * natoms, 0);
  ///// compute the regular dipole and the derivatives
  if (nt==1){
    for(unsigned i = 0; i < nbonds; ++i) {
      _calculate_bond(i, positions, borns, dipole, force, virial, unitcell_dipole);
    }  
  }
  else{
    _calculate_omp(positions, borns, dipole, force, virial, unitcell_dipole );
  }
  ///// compute the modified dipole and the derivatives
  std::vector<double> eff_borns(borns.begin(),borns.end());
  std::vector<double> eff_dipole( odim, 0);
  std::vector<double> eff_force ( odim * adof, 0);
  std::vector<double> eff_virial( odim * 9, 0);
  std::vector<long> dipole_upcount(odim, 0);
  // get the modifiled dipole 
  for(unsigned i = 0; i < natoms; ++i) {
    for(unsigned j = 0; j < odim; ++j) {
      double x = unitcell_dipole[odim*i+j];
      if (x==0) {continue;}
      else{
        if (x > 0)  {dipole_upcount[j] += 1;}
        x = tanh(x * rescale);
        unitcell_dipole[odim*i+j] = x;
        eff_dipole[j] += x;
      }
    }
  }
  // get the modified born charges
  for(unsigned i = 0; i < nbonds; ++i) {
    for(unsigned j = 0; j < odim; ++j) {
      for(unsigned k = 0; k < 3; ++k) {
        int atom1 = bonds[2*i];
        double x = unitcell_dipole[odim*atom1+j];
        eff_borns[odim * 3 * i + 3*j+k] *= (1 - x*x )*rescale;
      }
    }
  }
  // get effective force and virial with the old routine
  std::vector<double>  placeholder( odim, 0);
  if (nt==1){
    for(unsigned i = 0; i < nbonds; ++i) {
      _calculate_bond(i, positions, eff_borns, placeholder, eff_force, eff_virial, unitcell_dipole);
    }
  }
  else{
    _calculate_omp(positions, eff_borns, placeholder, eff_force, eff_virial, unitcell_dipole );
  }
  ///// return everything
  _feed_plumed( dipole, force, virial, eff_dipole, eff_force, eff_virial, dipole_upcount);
}

void 
LinearDipoleTanh::
_calculate_omp(const std::vector< Vector >  & positions, 
               std::vector<double>          & _borns,
               std::vector<double>          & dipole, 
               std::vector<double>          & force, 
               std::vector<double>          & virial, 
               std::vector<double>          & unitcell_dipole)
{
  unsigned nt=OpenMP::getNumThreads();
  // natoms = getNumberOfAtoms();
  std::vector<std::vector<double>> omp_force_allrank(nt);
  double start = omp_get_wtime();
  #pragma omp parallel num_threads(nt)
  {
    unsigned rank = OpenMP::getThreadNum();
    std::vector<double> omp_dipole( odim, 0);
    std::vector<double> omp_force ( odim * adof, 0);
    std::vector<double> omp_virial( odim * 9, 0);  
    std::vector<double> omp_unitcell_dipole(odim * natoms, 0);
    // iterate over all bonds, accumulate dipole, force and virial
    #pragma omp for
    for(unsigned i = 0; i < nbonds; ++i) {
      _calculate_bond(i, positions, _borns, omp_dipole, omp_force, omp_virial, omp_unitcell_dipole);
    }
    if (rank==0) {log.printf(" -- OMP thread time = %f s  \n", omp_get_wtime() - start);}
    // add up dipole and virial in serial, leave force for parallel reduction
    #pragma omp critical
    {
      for(unsigned j = 0; j < odim; ++j) {
        dipole[j] += omp_dipole[j];
        for(unsigned k = 0; k < natoms; ++k) {
          unitcell_dipole[odim*k+j] += omp_unitcell_dipole[odim*k+j];
        }
        for(unsigned k = 0; k < 9; ++k) {
          virial[9*j+k] += omp_virial[9*j+k];
          };
      };
    }
    // initialize force buffer
    omp_force_allrank[rank].resize(odim * adof);
    for(unsigned n = 0; n < odim * adof; ++n) {
      omp_force_allrank[rank][n] = omp_force[n];
    }
    // converge and sum with different for-loop schedule 
    #pragma omp barrier
    #pragma omp for
    for(unsigned n = 0; n < odim * adof; ++n) {
      for (unsigned j = 0; j < nt; ++j) {
        force[n] += omp_force_allrank[j][n];
      }
    }
  }
  log.printf(" OMP total time = %f s  \n", omp_get_wtime() - start);
}

void 
LinearDipoleTanh::
_calculate_bond(unsigned                    & id_bond, 
                const std::vector< Vector > & positions, 
                std::vector<double>         & _borns,
                std::vector<double>         & _dipole, 
                std::vector<double>         & _force, 
                std::vector<double>         & _virial, 
                std::vector<double>         & _unitcell_dipole)
{
  int atom1 = bonds[2*id_bond];
  int atom2 = bonds[2*id_bond+1];
  Vector r12 = _get_distance( positions, atom1, atom2);
  std::vector<double> born_tensor(
    _borns.begin() + odim*3 * id_bond, 
    _borns.begin() + odim*3 * (id_bond+1)
    ); // borns[idx_bond, j, *]
  for(unsigned j = 0; j < odim; ++j) {
    ////////////////////   use plumed's Vector/Tensor class      ///////////////////////
    Vector born_vector(born_tensor[3*j],born_tensor[3*j+1], born_tensor[3*j+2]);
    Tensor bond_virial(-r12, born_vector);
    double bdr = dotProduct(born_vector, r12);
    _dipole[j] += bdr;
    _unitcell_dipole[odim * atom1+j] += bdr;
    // this bond's contribution to the virial v = r x F
    for(unsigned k = 0; k < 3; ++k) {
      _force [j*adof + 3*atom2+k] -= born_vector[k]; //  F[j, atom2, k] += - d(dipole[j])/d(r2[k]) 
      _force [j*adof + 3*atom1+k] += born_vector[k]; //  F[j, atom1, k] += - d(dipole[j])/d(r1[k]) 
      for(unsigned l = 0; l < 3; ++l) {
        _virial[9*j + 3*k + l] += bond_virial[k][l];    // virial[idx_odim,k,l] += (r2-r1)[k] * borns[idx_bond,  idx_odim, l]
        };
    }
  }
}


}
}
