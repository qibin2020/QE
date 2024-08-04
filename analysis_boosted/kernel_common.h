#include <ROOT/RVec.hxx>
#include <TLorentzVector.h>
#include <map>

using Vec_i = const ROOT::VecOps::RVec<int>&;
using Vec_d = const ROOT::VecOps::RVec<double>&;
using Vec_f = const ROOT::VecOps::RVec<float>&;

using Four_Mom = TLorentzVector;
using Charged_Lepton = std::pair<int,Four_Mom>;
using Four_Mom_ref = const Four_Mom &;
using Three_Mom = TVector3;
using Three_Mom_ref = const TVector3 &;
using Frame = std::map<std::string, Three_Mom>;
using Frame_ref = const Frame &;

using Jet_indices = ROOT::VecOps::RVec<int>;
using Jet_indices_ref = const Jet_indices &;

#pragma link C++ class Charged_Lepton;
#pragma link C++ class Frame;
