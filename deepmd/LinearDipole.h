#include "Common.h"
#include "Colvar.h"
#include "ActionRegister.h"

#include <string>
#include <cmath>

namespace PLMD {
namespace colvar {
 
constexpr unsigned int ODIM = 3;
const std::array<std::string, ODIM> cpnts = {"x", "y", "z"};

 //+PLUMEDOC COLVAR DIPOLE
/*
Calculate a 3D collective variable for a group of atoms, linearly dependent on the displacement between atoms.

When running with periodic boundary conditions, the atoms should be
in the proper periodic image. This is done automatically since PLUMED 2.5,
by considering the ordered list of atoms and rebuilding the molecule with a procedure
that is equivalent to that done in \ref WHOLEMOLECULES . Notice that
rebuilding is local to this action. This is different from \ref WHOLEMOLECULES
which actually modifies the coordinates stored in PLUMED.

 \par Examples

The following tells plumed to calculate the dipole of the group of atoms containing
the atoms from 1-10 and print it every 5 steps
\plumedfile
d: DIPOLE GROUP=1-10
PRINT FILE=output STRIDE=5 ARG=d
\endplumedfile

 */
//+ENDPLUMEDOC

class LinearDipole : public Colvar {
  std::vector<AtomNumber> atoms; 
  bool components;
  bool nopbc;
  int  odim = ODIM;
  std::vector<int> bonds;
  std::vector<FLOAT_PREC> borns;
public:
  explicit LinearDipole (const ActionOptions&);
  void calculate () override;
  static void registerKeywords (Keywords& keys);
protected:
  std::vector<int> model;
private:
  void load_bonds(std::vector<int> & bonds, 
           const std::string & filename );
  void load_borns(std::vector<FLOAT_PREC> & borns, 
           const std::string & filename );

};


}
}
