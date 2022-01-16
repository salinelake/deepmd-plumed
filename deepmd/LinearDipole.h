#pragma once
#include "Common.h"
#include "colvar/Colvar.h"
#include "colvar/ActionRegister.h"
#include "core/PlumedMain.h"

#include <string>
#include <cmath>

constexpr unsigned int odim = 3;
const std::array<std::string, odim> cpnts = {"x", "y", "z"};

namespace PLMD {
namespace colvar {
 //+PLUMEDOC COLVAR LinearDipole
/*
Calculate a 3D collective variable for a group of atoms, linearly dependent on the displacement between atoms.

 \par Examples

The following tells plumed to calculate the dipole of the group of atoms containing
the atoms from 1-10 and print it every 5 steps
\plumedfile
d: LinearDipole
PRINT FILE=output STRIDE=5 ARG=d
\endplumedfile

 */
//+ENDPLUMEDOC

class LinearDipole : public Colvar {
public:
  explicit LinearDipole (const ActionOptions&);
  void calculate () override;
  static void registerKeywords (Keywords& keys);
protected:
  std::vector<AtomNumber> atoms; 
  bool components;
  bool nopbc;
  std::vector<int> model;
private:
  std::vector<int> bonds;
  std::vector<double> borns;
  void load_bonds(std::vector<int> & bonds, 
           const std::string & filename );
  void load_borns(std::vector<double> & borns, 
           const std::string & filename );
};


}
}
