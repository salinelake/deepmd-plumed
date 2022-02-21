#pragma once
#include "Common.h"
#include "colvar/Colvar.h"
#include "colvar/ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/OpenMP.h"
#include "tools/Vector.h"

#include <omp.h>

#include <string>
#include <cmath>

constexpr unsigned int odim = 3;
const std::array<std::string, odim> cpnts = {"x", "y", "z"};

namespace PLMD {
namespace colvar {
 //+PLUMEDOC COLVAR LinearDipole
/*
Calculate a collective variable for a group of atoms, linearly dependent on the displacement between atoms.

 \par Examples

The following tells plumed to print the global dipole linearly depending on static Born charge tensors specified by bonds.raw and borns.raw
\plumedfile
dipole: LinearDipole BONDS=bonds.raw BORNS=borns.raw
PRINT ARG=dipole STRIDE=100 FILE=COLVAR
\endplumedfile

 */
//+ENDPLUMEDOC

class LinearDipoleCount : public Colvar {
public:
  explicit LinearDipoleCount (const ActionOptions&);
  void calculate () override;
  static void registerKeywords (Keywords& keys);
protected:
  std::vector<AtomNumber> atoms; 
  bool nopbc;
  bool components;
private:
  unsigned nbonds;
  unsigned natoms;
  unsigned adof;
  std::vector<int> bonds;
  std::vector<double> borns;
  void load_bonds(std::vector<int> & bonds, 
           const std::string & filename );
  void load_borns(std::vector<double> & borns, 
           const std::string & filename );
  Vector _get_distance(const std::vector< Vector >  & positions, int & a1, int & a2);
  void _feed_plumed(std::vector<double> & dipole, std::vector<double> & force, std::vector<double> & virial, std::vector<long> & dipole_upcount);
  void _calculate_omp(const std::vector< Vector >  & positions, std::vector<double> & dipole, std::vector<double> & force, std::vector<double> & virial, std::vector<double> &  unitcell_dipole);
  void _calculate_bond( unsigned & id_bond, const std::vector< Vector >  & positions, std::vector<double> & _dipole, std::vector<double> & _force, std::vector<double> & _virial, std::vector<double> &  unitcell_dipole);

};


}
}
