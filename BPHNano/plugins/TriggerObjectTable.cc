// system include files
#include <memory>
#include <sstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

class TriggerObjectTableBParkProducer : public edm::stream::EDProducer<> {
public:
  explicit TriggerObjectTableBParkProducer(const edm::ParameterSet &iConfig) :
    name_(iConfig.getParameter<std::string>("name")),
    src_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("src"))),
    l1Muon_(consumes<l1t::MuonBxCollection>(iConfig.getParameter<edm::InputTag>("l1Muon")))
  {
    std::vector<edm::ParameterSet> selPSets = iConfig.getParameter<std::vector<edm::ParameterSet>>("selections");
    sels_.reserve(selPSets.size());
    std::stringstream idstr, qualitystr;
    idstr << "ID of the object: ";
    for (auto & pset : selPSets) {
      sels_.emplace_back(pset);
      idstr << sels_.back().id << " = " << sels_.back().name;
      if (sels_.size() < selPSets.size()) idstr << ", ";
      if (!sels_.back().qualityBitsDoc.empty()) {
        qualitystr << sels_.back().qualityBitsDoc << " for " << sels_.back().name << "; ";
      }
    }
    idDoc_ = idstr.str();
    bitsDoc_ = qualitystr.str();

    produces<nanoaod::FlatTable>();
  }

  ~TriggerObjectTableBParkProducer() override {}

private:
  void produce(edm::Event&, edm::EventSetup const&) override ;

  std::string name_;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> src_;
  std::string idDoc_, bitsDoc_;

  edm::EDGetTokenT<l1t::MuonBxCollection> l1Muon_;

  struct SelectedObject {
    std::string name;
    int id;
    StringCutObjectSelector<pat::TriggerObjectStandAlone> cut;
    StringCutObjectSelector<pat::TriggerObjectStandAlone> l1cut, l1cut_2, l2cut;
    float       l1DR2, l1DR2_2, l2DR2;
    StringObjectFunction<pat::TriggerObjectStandAlone> qualityBits;
    std::string qualityBitsDoc;

    SelectedObject(const edm::ParameterSet & pset) :
      name(pset.getParameter<std::string>("name")),
      id(pset.getParameter<int>("id")),
      cut(pset.getParameter<std::string>("sel")),
      l1cut(""), l1cut_2(""), l2cut(""),
      l1DR2(-1), l1DR2_2(-1), l2DR2(-1),
      qualityBits(pset.getParameter<std::string>("qualityBits")),
      qualityBitsDoc(pset.getParameter<std::string>("qualityBitsDoc"))
    {
      if (pset.existsAs<std::string>("l1seed")) {
        l1cut = StringCutObjectSelector<pat::TriggerObjectStandAlone>(pset.getParameter<std::string>("l1seed"));
        l1DR2 = std::pow(pset.getParameter<double>("l1deltaR"), 2);
      }
      if (pset.existsAs<std::string>("l1seed_2")) {
        l1cut_2 = StringCutObjectSelector<pat::TriggerObjectStandAlone>(pset.getParameter<std::string>("l1seed_2"));
        l1DR2_2 = std::pow(pset.getParameter<double>("l1deltaR_2"), 2);
      }
      if (pset.existsAs<std::string>("l2seed")) {
        l2cut = StringCutObjectSelector<pat::TriggerObjectStandAlone>(pset.getParameter<std::string>("l2seed"));
        l2DR2 = std::pow(pset.getParameter<double>("l2deltaR"), 2);
      }
    }

    bool match(const pat::TriggerObjectStandAlone & obj) const {
      return cut(obj);
    }
  };

  std::vector<SelectedObject> sels_;
};

// ------------ method called to produce the data  ------------
void
TriggerObjectTableBParkProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> src;
  iEvent.getByToken(src_, src);

  std::vector<std::pair<const pat::TriggerObjectStandAlone *, const SelectedObject *>> selected;
  for (const auto & obj : *src) {
    for (const auto & sel : sels_) {
      if (sel.match(obj)) {
        selected.emplace_back(&obj, &sel);
        break;
      }
    }
  }

  // Self-cleaning
  std::map<const pat::TriggerObjectStandAlone *, int> selected_bits;
  for (unsigned int i = 0; i < selected.size(); ++i) {
    const auto & obj = *selected[i].first;
    const auto & sel = *selected[i].second;
    selected_bits[&obj] = int(sel.qualityBits(obj));

    for (unsigned int j = 0; j < i; ++j) {
      const auto & obj2 = *selected[j].first;
      const auto & sel2 = *selected[j].second;
      if (sel.id == sel2.id && abs(obj.pt() - obj2.pt()) < 1e-6 && deltaR2(obj, obj2) < 1e-6) {
        selected_bits[&obj2] |= selected_bits[&obj]; //Keep filters from all the objects
        selected.erase(selected.begin() + i);
        i--;
      }
    }
  }

  edm::Handle<l1t::MuonBxCollection> l1Muon;
  iEvent.getByToken(l1Muon_, l1Muon);

  std::vector<pair<pat::TriggerObjectStandAlone, int>> l1Objects;


  for (l1t::MuonBxCollection::const_iterator it = l1Muon->begin(0); it != l1Muon->end(0); it++) {
    pat::TriggerObjectStandAlone l1obj(it->p4());
    l1obj.setCollection("L1Mu");
    l1obj.addTriggerObjectType(trigger::TriggerL1Mu);
    l1obj.setCharge(it->charge());
    l1Objects.emplace_back(l1obj, it->hwIso());
  }



  unsigned int nobj = selected.size();
  std::vector<float> pt(nobj, 0), eta(nobj, 0), phi(nobj, 0), l1pt(nobj, 0), l1pt_2(nobj, 0), l2pt(nobj, 0);
  std::vector<int>   id(nobj, 0), bits(nobj, 0), l1iso(nobj, 0), l1charge(nobj, 0);
  for (unsigned int i = 0; i < nobj; ++i) {
    const auto & obj = *selected[i].first;
    const auto & sel = *selected[i].second;
    pt[i] = obj.pt();
    eta[i] = obj.eta();
    phi[i] = obj.phi();
    id[i] = sel.id;
    bits[i] = selected_bits[&obj];
    if (sel.l1DR2 > 0) {
      float best = sel.l1DR2;
      for (const auto & l1obj : l1Objects) {
        const auto & seed = l1obj.first;
        float dr2 = deltaR2(seed, obj);
        if (dr2 < best && sel.l1cut(seed)) {
          l1pt[i] = seed.pt();
          l1iso[i] = l1obj.second;
          l1charge[i] = seed.charge();
        }
      }
    }
    if (sel.l1DR2_2 > 0) {
      float best = sel.l1DR2_2;
      for (const auto & l1obj : l1Objects) {
        const auto & seed = l1obj.first;
        float dr2 = deltaR2(seed, obj);
        if (dr2 < best && sel.l1cut_2(seed)) {
          l1pt_2[i] = seed.pt();
        }
      }
    }
    if (sel.l2DR2 > 0) {
      float best = sel.l2DR2;
      for (const auto & seed : *src) {
        float dr2 = deltaR2(seed, obj);
        if (dr2 < best && sel.l2cut(seed)) {
          l2pt[i] = seed.pt();
        }
      }
    }
  }

  auto tab  = std::make_unique<nanoaod::FlatTable>(nobj, name_, false, false);
  tab->addColumn<int>("id", id, idDoc_ );
  tab->addColumn<float>("pt", pt, "pt", 12);
  tab->addColumn<float>("eta", eta, "eta",  12);
  tab->addColumn<float>("phi", phi, "phi", 12);
  tab->addColumn<float>("l1pt", l1pt, "pt of associated L1 seed", 8);
  tab->addColumn<int>("l1iso", l1iso, "iso of associated L1 seed" );
  tab->addColumn<int>("l1charge", l1charge, "charge of associated L1 seed");
  tab->addColumn<float>("l1pt_2", l1pt_2, "pt of associated secondary L1 seed", 8);
  tab->addColumn<float>("l2pt", l2pt, "pt of associated 'L2' seed (i.e. HLT before tracking/PF)", 10);
  tab->addColumn<int>("filterBits", bits, "extra bits of associated information: " + bitsDoc_);
  iEvent.put(std::move(tab));
}


//define this as a plug-in
DEFINE_FWK_MODULE(TriggerObjectTableBParkProducer);
