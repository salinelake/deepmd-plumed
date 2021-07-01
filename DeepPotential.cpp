/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)
   See http://www.plumed.org for more information.
   This file is part of plumed, version 2.
   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.
   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "core/ActionAtomistic.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/Units.h"
#include "colvar/Colvar.h"
#include "colvar/ActionRegister.h"

#include "Common.h"
#include "deepmd/DeepPot.h"

namespace PLMD {
namespace dp_plmd { // to avoid conflicts with libdeepmd
//+PLUMEDOC COLVAR DEEPPOTENTIAL
/*
Calculate the potential energy for the system using a deep potential model.

A binary graph file of the model and a text file specify the type of each atom are needed.
The type file should be a list of integers with the i-th element correspond to the type 
in the deep tensor model of the i-th atom in the system (or specified by you).
The output is scaled by a factor given by UNIT_CVT, which defaults to 1.
By default will use periodic boundary conditions, which will be handled by 
the deep tensor model automatically. In case NOPBC flag is specified, the box
will be ignored and there will be no the pbc handling.

\par Examples

Here's a simple example showing how to use this CV (also the default values of the keywords):
\plumedfile
dp: DEEPPOTENTIAL MODEL=graph.pb ATYPE=type.raw UNIT_CVT=96.487
\endplumedfile
*/
//+ENDPLUMEDOC
class DeepPotential : public Colvar {
  std::vector<AtomNumber> atoms; 
  bool nopbc;
public:
  explicit DeepPotential(const ActionOptions&);
  void calculate() override;
  static void registerKeywords(Keywords& keys);
private:
  deepmd::DeepPot dp; 
  std::vector<int> atype;
  double energy_unit;
  double length_unit;
};

PLUMED_REGISTER_ACTION(DeepPotential,"DEEPPOTENTIAL")

void DeepPotential::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  keys.add("atoms","ATOMS","the group of atoms we are calculating the potential energy for (defaults to the whole system)");
  keys.add("compulsory","MODEL","graph.pb","the potential model binary graph file");
  keys.add("compulsory","ATYPE","type.raw" ,"the file specify the type (in the model) of each atom");
  keys.add("optional","UNIT_CVT","the unit conversion constant of output energy (will be multiplied to the graph output, default is converting to plumed unit)");
}

DeepPotential::DeepPotential(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  nopbc(false),
  energy_unit(global_energy_unit),
  length_unit(global_length_unit)
{
  // make sure the length unit passed to graph is Angstrom
  length_unit /= plumed.getAtoms().getUnits().getLength();
  // by default use the plumed unit. but user can change this
  energy_unit /= plumed.getAtoms().getUnits().getEnergy();

  parseAtomList("ATOMS",atoms);
  parseFlag("NOPBC",nopbc);
  std::string graph_file;
  parse("MODEL", graph_file);
  std::string type_file;
  parse("ATYPE", type_file);
  parse("UNIT_CVT", energy_unit); // will overwrite the default unit

  checkRead();
  
  addValueWithDerivatives(); setNotPeriodic();

  // default use all atoms; otherwise warn the user
  if (atoms.size() == 0) {
    atoms.resize(getTotAtoms());
    for(unsigned i = 0; i < getTotAtoms(); ++i) {
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
  if (atoms.size() != getTotAtoms()) { 
    log.printf("  # of atoms provided: %d != # of atoms in the system: %d\n", atoms.size(), getTotAtoms()); 
    log.printf("  Please make sure you know what you are doing!\n");
  }

  if(nopbc) { log.printf("  without periodic boundary conditions\n"); }
  else      { log.printf("  using periodic boundary conditions\n"); }

  // load deepmd model from graph file; check dimention is correct
  log.printf("  using graph file:  %s \n", graph_file.c_str());
  dp.init(graph_file);
  // dp.print_summary("  DeePMD: ");
  log.printf("  DeepPotential model initialized successfully\n");

  // load atom types that will be feed into the graph
  load_atype(atype, type_file);
  if (atype.size() != atoms.size()) { 
    throw std::runtime_error( "invalid atom type file! the size should be equal to the number of atoms" ); 
  }
  log.printf("  assign type to atoms\n");
  for(unsigned i = 0; i < atype.size(); ++i) {
    log.printf("  %d", atype[i]);
  }
  log.printf("  \n");

  // the output unit will be multiple to the graph output
  log.printf("  output unit conversion set to %f\n", energy_unit);

  requestAtoms(atoms); 
}

void DeepPotential::calculate()
{
  if (!nopbc) { makeWhole(); }
  unsigned N = getNumberOfAtoms();
  double _energy;
  std::vector<FLOAT_PREC> _force (N * 3);
  std::vector<FLOAT_PREC> _virial(9);
  std::vector<FLOAT_PREC> _coord (N * 3);
  std::vector<FLOAT_PREC> _box   (9); 
  IndexConverter ic(N, 3);
  
  // copy atom coords
  for (unsigned i = 0; i < N; ++i) {
    for (unsigned j = 0; j < 3; ++j) {
      _coord[ic.f(i,j)] = getPosition(i)[j] / length_unit;
    }
  }
  // copy box tensor
  Tensor box = getBox();
  if (nopbc) { 
    _box.clear(); // if size is not 9 then no pbc
  } else {
    for (unsigned i = 0; i < 3; ++i) {
      for (unsigned j = 0; j < 3; ++j) {
        _box[i * 3 + j] = box(i, j) / length_unit;
      }
    }
  }

  dp.compute(_energy, _force, _virial, _coord, atype, _box);

  // get back energy
  setValue(_energy * energy_unit);
  // get back force
  for(unsigned i = 0; i < N; i++) {
    setAtomsDerivatives(
      i,
      - Vector(_force[ic.f(i,0)],
               _force[ic.f(i,1)], 
               _force[ic.f(i,2)]) 
      * energy_unit / length_unit
    );
  }
  // get back virial
  // setBoxDerivativesNoPbc(); // ask plumed to calculate virial
  setBoxDerivatives(
    Tensor( // we manually transpose here
      _virial[0], _virial[3], _virial[6], 
      _virial[1], _virial[4], _virial[7], 
      _virial[2], _virial[5], _virial[8]) 
  );
  }

}
}