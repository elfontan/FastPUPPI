// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"

#include "L1Trigger/Phase2L1ParticleFlow/interface/L1TPFUtils.h"

#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <unordered_set>

// ---------------------------------------------------------------
// Class: L1PFCandTableProducer
// ---------------------------------------------------------------
class L1PFCandTableProducer : public edm::global::EDProducer<>  {
    public:
        explicit L1PFCandTableProducer(const edm::ParameterSet&);
        ~L1PFCandTableProducer();

    private:
        virtual void produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const override;

        StringCutObjectSelector<reco::Candidate> sel_;

        struct ExtraVar {
            std::string name, expr;
            StringObjectFunction<reco::Candidate> func;
            ExtraVar(const std::string & n, const std::string & expr) : name(n), expr(expr), func(expr, true) {}
        };
        std::vector<ExtraVar> extraVars_;

        struct CandRecord {
            public:
                std::string coll;
                edm::EDGetTokenT<reco::CandidateView> src;
                StringCutObjectSelector<reco::Candidate> sel;
                
                CandRecord(const std::string & name, const edm::EDGetTokenT<reco::CandidateView> & tag, const edm::ParameterSet & pset) :
                    coll(name), src(tag), 
                    sel(pset.existsAs<std::string>(name+"_sel") ? pset.getParameter<std::string>(name+"_sel") : "", true) {}
        };
        std::vector<CandRecord> cands_;
        std::vector<CandRecord> gencands_;
};

L1PFCandTableProducer::L1PFCandTableProducer(const edm::ParameterSet& iConfig) :
    sel_(iConfig.getParameter<std::string>("commonSel"), true)
{
    edm::ParameterSet cands = iConfig.getParameter<edm::ParameterSet>("cands");
    auto candnames = cands.getParameterNamesForType<edm::InputTag>();
    for (const std::string & name : candnames) {
        cands_.emplace_back(name, consumes<reco::CandidateView>(cands.getParameter<edm::InputTag>(name)), cands);
        produces<nanoaod::FlatTable>(name+"Cands");
    }

    gencands_.emplace_back("GenCands", consumes<reco::CandidateView>(edm::InputTag("genParticlesForMETAllVisible")), iConfig);
    //gencands_.emplace_back("GenCands", consumes<reco::CandidateView>(edm::InputTag("genInAcceptance")), iConfig);

    if (iConfig.existsAs<edm::ParameterSet>("moreVariables")) {
        edm::ParameterSet vars = iConfig.getParameter<edm::ParameterSet>("moreVariables");
        auto morenames = vars.getParameterNamesForType<std::string>();
        for (const std::string & name : morenames) {
            extraVars_.emplace_back(name, vars.getParameter<std::string>(name));
        }
    }
 }

double calculate_deltaR(double eta1, double phi1, double eta2, double phi2) {
    return reco::deltaR(eta1, phi1, eta2, phi2);
}

L1PFCandTableProducer::~L1PFCandTableProducer() { }


// =====================================================================================
// PUPPI ML at L1: Target studies
// Studying gen ratios + multiplicities around *neutral RECO seeds only* (dR = 0.2)
//   - Three categories (both RECO and GEN):
//       1) charged-around-neutral-seed   (includes seed only if it is charged -> never, since seed is neutral)
//       2) neutral-around-neutral-seed   (includes seed)
//       3) all-around-neutral-seed       (includes seed)
//   - GEN thresholds kept: charged >2 GeV, neutral >1 GeV (status==1 for sums/counts)
//   - Updated matching logic (isGenMatched) using GEN status1 and passing thresholds, within dR < 0.1 and a minimum pt compatibility
// =====================================================================================

static inline double dR(double eta1,double phi1,double eta2,double phi2){
  return reco::deltaR(eta1,phi1,eta2,phi2);
}

static inline std::pair<double,double> caloEtaPhi(const reco::Candidate* c, double bz) {
  math::XYZTLorentzVector vtx(c->vx(), c->vy(), c->vz(), 0.);
  return l1tpf::propagateToCalo(c->p4(), vtx, c->charge(), bz);
}

// A "good enough" uniqueness signature w/o references/pointers.
// If you later have unique keys in ntuples, swap this for those.
static inline std::string recoKeyNoRef(const reco::Candidate* c, double bz) {
  auto calo = caloEtaPhi(c,bz);
  std::ostringstream ss;
  ss << c->pdgId() << "_"
     << c->charge() << "_"
     << std::fixed << std::setprecision(3)
     << c->pt() << "_"
     << calo.first << "_"
     << calo.second;
  return ss.str();
}

struct ConeMult {
  // RECO (seed excluded)
  int nRecoAll=0, nRecoCh=0, nRecoNe=0;
  int nRecoPho=0, nRecoNHad=0, nRecoChHad=0;

  // GEN (all in cone, seed doesn't exist on gen side)
  int nGenAll=0, nGenCh=0, nGenNe=0;
  int nGenPho=0, nGenNHad=0, nGenChHad=0;

  // GEN (status1 + thresholds) in cone  <-- THIS is what you requested
  int nGenS1ThrAll=0, nGenS1ThrCh=0, nGenS1ThrNe=0;
  int nGenS1ThrPho=0, nGenS1ThrNHad=0, nGenS1ThrChHad=0;
};

struct CatSums {
  // "Sum" means pT sum passing thresholds (GEN) and simply pT sum (RECO).
  // GEN uses: status==1 AND (charged>2 or neutral>1), and then category filter.
  double recoSum=0.0;
  double genSum=0.0;
  double ratio() const { return (recoSum>0 ? genSum/recoSum : -1.0); }
};

struct ConeResult {
  CatSums ch;   // (closure) charged-only in cone
  CatSums ne;   // (closure) neutral-only in cone
  CatSums all;  // (closure) all in cone
  ConeMult mult;

  // 3 categories around the neutral seed: (PU-correction denominators) seed + X
  double recoDen_seedPlusCh  = 0.0;
  double recoDen_seedPlusNe  = 0.0;
  double recoDen_seedPlusAll = 0.0;

  int isGenMatched = 0; 
  double seedCaloEta=0., seedCaloPhi=0.;

  // new: matched neutral GEN info
  float matchedGenNeutralPt = 0.f;   // 0 if none
  int   hasMatchedGenNeutral = 0;    // 1 if found
  double genChargedPlusMatchedNeutral = 0.0;
  int matchedGenNeutralPdgId = 0;
  float matchedGenNeutralDR  = -1.f;
};

// --- GEN classification (status1 + thresholds) ---
static inline bool genIsStatus1(const reco::GenParticle* gp){ return gp && gp->status()==1; }
static inline bool genIsNeutral(const reco::Candidate* c){ return c->charge()==0; }
static inline bool genPassThr(const reco::Candidate* c){
  // thresholds requested:
  //   - neutral >1 GeV
  //   - charged >2 GeV
  //if (c->charge()==0) return c->pt() > 1.0;
  //return c->pt() > 2.0;
  if (c->charge()==0) return c->pt() > 1.0;
  return c->pt() > 2.0;
}

static inline bool recoIsNeutral(const reco::Candidate* c){ return c->charge()==0; }

static inline void printConeDebugFull(
    const edm::EventID& eid,
    const reco::Candidate* seed,
    const std::vector<const reco::Candidate*>& reco_selected,
    const std::vector<const reco::Candidate*>& gen_selected,
    const ConeResult& out,
    double coneSize,
    double bz
){
  auto seedCalo = caloEtaPhi(seed,bz);

  const float r_ch  = (out.recoDen_seedPlusCh  > 0) ? (out.genChargedPlusMatchedNeutral / out.recoDen_seedPlusCh)  : -1.f;
  const float r_ne  = (out.recoDen_seedPlusNe  > 0) ? (out.ne.genSum / out.recoDen_seedPlusNe)  : -1.f;
  const float r_all = (out.recoDen_seedPlusAll > 0) ? (out.all.genSum / out.recoDen_seedPlusAll) : -1.f;

  std::cout << "\n============================================================\n";
  std::cout << "EVENT " << eid.run() << ":" << eid.luminosityBlock() << ":" << eid.event() << "\n";
  std::cout << "SEED  pdgId=" << seed->pdgId()
            << " q=" << seed->charge()
            << " pt=" << seed->pt()
            << " eta=" << seed->eta()
            << " phi=" << seed->phi()
            << " calo(eta,phi)=" << seedCalo.first << "," << seedCalo.second
            << "\n";

  std::cout << "-- MATCH --\n";
  std::cout << " hasMatchedGenNeutral=" << out.hasMatchedGenNeutral
            << " \n matchedGenNeutralPt=" << out.matchedGenNeutralPt
            << " matchedGenNeutralPdgId=" << out.matchedGenNeutralPdgId
            << " matchedGenNeutralDR=" << out.matchedGenNeutralDR
            << "\n";

  std::cout << "\n-- NUMERATORS (status1+thr) --\n";
  std::cout << " genChargedPtSum0p2=" << out.ch.genSum << "\n";
  std::cout << " genNeutralPtSum0p2=" << out.ne.genSum << "\n";
  std::cout << " genPtSum0p2=" << out.all.genSum << "\n";
  std::cout << " genChPlusMatchedNe0p2=" << out.genChargedPlusMatchedNeutral << "\n";

  std::cout << "\n-- DENOMINATORS (seed+RECO in cone) --\n";
  std::cout << " recoDen_seedPlusCh0p2=" << out.recoDen_seedPlusCh << "\n";
  std::cout << " recoDen_seedPlusNe0p2=" << out.recoDen_seedPlusNe << "\n";
  std::cout << " recoDen_seedPlusAll0p2=" << out.recoDen_seedPlusAll << "\n";

  std::cout << "\n-- RATIOS --\n";
  std::cout << " ratioPU_ch=" << r_ch << "  ratioPU_ne=" << r_ne << "  ratioPU_all=" << r_all << "\n";

  std::cout << "-- MULT (RECO, seed excluded) --\n";
  std::cout << " RECO all=" << out.mult.nRecoAll << " ch=" << out.mult.nRecoCh << " ne=" << out.mult.nRecoNe
            << " pho=" << out.mult.nRecoPho << " nHad=" << out.mult.nRecoNHad << " chHad=" << out.mult.nRecoChHad
            << "\n";

  std::cout << "-- MULT (GEN all in cone) --\n";
  std::cout << " GEN all=" << out.mult.nGenAll << " ch=" << out.mult.nGenCh << " ne=" << out.mult.nGenNe
            << " pho=" << out.mult.nGenPho << " nHad=" << out.mult.nGenNHad << " chHad=" << out.mult.nGenChHad
            << "\n";

  std::cout << "-- MULT (GEN status1+thr in cone) --\n";
  std::cout << " GEN_s1thr all=" << out.mult.nGenS1ThrAll << " ch=" << out.mult.nGenS1ThrCh << " ne=" << out.mult.nGenS1ThrNe
            << " \n\t pho=" << out.mult.nGenS1ThrPho << " nHad=" << out.mult.nGenS1ThrNHad << " chHad=" << out.mult.nGenS1ThrChHad
            << "\n";

  // ---- full lists ----
  std::cout << "\n-- RECO in cone (seed excluded) --\n";
  {
    std::unordered_set<std::string> seen;
    for (const auto* c : reco_selected) {
      if (c == seed) continue;
      const std::string k = recoKeyNoRef(c,bz);
      if (!seen.insert(k).second) continue;

      auto calo = caloEtaPhi(c,bz);
      const double dr = dR(seedCalo.first,seedCalo.second,calo.first,calo.second);
      if (dr >= coneSize) continue;

      std::cout << " reco: pt=" << std::setw(7) << c->pt()
                << " eta=" << std::setw(7) << c->eta()
                << " phi=" << std::setw(7) << c->phi()
                << " pdgId=" << std::setw(6) << c->pdgId()
                << " q=" << std::setw(2) << c->charge()
                << " dR=" << dr
                << "\n";
    }
  }

  std::cout << "\n-- GEN in cone (ALL) --\n";
  for (const auto* g : gen_selected) {
    const auto* gp = dynamic_cast<const reco::GenParticle*>(g);
    if (!gp) continue;

    auto calo = caloEtaPhi(g,bz);
    const double dr = dR(seedCalo.first,seedCalo.second,calo.first,calo.second);
    if (dr >= coneSize) continue;

    const bool thr = genPassThr(g);

    std::cout << " gen:  pt=" << std::setw(7) << g->pt()
              << " eta=" << std::setw(7) << g->eta()
              << " phi=" << std::setw(7) << g->phi()
              << " pdgId=" << std::setw(6) << gp->pdgId()
              << " q=" << std::setw(2) << g->charge()
              << " st=" << std::setw(2) << gp->status()
              << " passThr=" << (thr ? 1 : 0)
              << " dR=" << dr
              << "\n";
  }

  std::cout << "============================================================\n";
}

// ----------------------------------------------------------------------
// Main worker: neutral RECO seed only
// ----------------------------------------------------------------------
static inline ConeResult computeConeAroundNeutralRecoSeed(
    const reco::Candidate* seed,
    const std::vector<const reco::Candidate*>& reco_selected,
    const std::vector<const reco::Candidate*>& gen_selected,
    double coneSize, // use 0.2
    double bz
){
  ConeResult out;

  // enforce: seed must be neutral reco
  auto seedCalo = caloEtaPhi(seed,bz);
  out.seedCaloEta = seedCalo.first;
  out.seedCaloPhi = seedCalo.second;

  // --------------------------
  // RECO loop (unique-by-key)
  // --------------------------
  std::unordered_set<std::string> seenRecoKeys;
  seenRecoKeys.reserve(reco_selected.size()*2);

  out.recoDen_seedPlusCh  = seed->pt();
  out.recoDen_seedPlusNe  = seed->pt();
  out.recoDen_seedPlusAll = seed->pt();
  
  for (const auto* c : reco_selected) {
    if (c == seed) continue;
    
    // unique protection (no refs, no pointers)
    const std::string key = recoKeyNoRef(c,bz);
    if (!seenRecoKeys.insert(key).second) continue;

    auto calo = caloEtaPhi(c,bz);
    const double dr = dR(out.seedCaloEta,out.seedCaloPhi,calo.first,calo.second);
    if (dr >= coneSize) continue;

    // closure sums (NO seed here)
    out.all.recoSum += c->pt();
    if (c->charge()==0) out.ne.recoSum += c->pt();
    else                out.ch.recoSum += c->pt();

    // seed+X denominators (seed already included above)
    out.recoDen_seedPlusAll += c->pt();
    if (c->charge()==0) out.recoDen_seedPlusNe += c->pt();
    else                out.recoDen_seedPlusCh += c->pt();
    
    // multiplicities
    out.mult.nRecoAll++;
    if (c->charge()==0) out.mult.nRecoNe++; else out.mult.nRecoCh++;
    const int apdg = std::abs(c->pdgId());
    if (apdg==22)  out.mult.nRecoPho++;
    if (apdg==130) out.mult.nRecoNHad++;
    if (apdg==211) out.mult.nRecoChHad++;

  }

  // -------------------------------------------------------------------------
  // GEN loop (all-in-cone mult) + (status1+thr sums+mult) + improved matching
  // -------------------------------------------------------------------------
  const double matchDR  = 0.2;
  const double relPtTol = 0.5;

  double bestDR = 1e9;
  const reco::GenParticle* bestGenNeutral = nullptr;

  for (const auto* g : gen_selected) {
    const auto* gp = dynamic_cast<const reco::GenParticle*>(g);
    if (!gp) continue;

    auto calo = caloEtaPhi(g,bz);
    const double dr = dR(out.seedCaloEta,out.seedCaloPhi,calo.first,calo.second);
    if (dr >= coneSize) continue;

    // --- GEN multiplicities (ALL, no status/thr) ---
    out.mult.nGenAll++;
    if (g->charge()==0) out.mult.nGenNe++; else out.mult.nGenCh++;
    const int apdg = std::abs(gp->pdgId());
    if (apdg==22)  out.mult.nGenPho++;
    if (apdg==130) out.mult.nGenNHad++;
    if (apdg==211) out.mult.nGenChHad++;

    // --- status1 + thresholds gate ---
    if (!genIsStatus1(gp)) continue;
    if (!genPassThr(g)) continue;

    // --- GEN multiplicities (status1+thr) ---
    out.mult.nGenS1ThrAll++;
    if (g->charge()==0) out.mult.nGenS1ThrNe++; else out.mult.nGenS1ThrCh++;
    if (apdg==22)  out.mult.nGenS1ThrPho++;
    if (apdg==130) out.mult.nGenS1ThrNHad++;
    if (apdg==211) out.mult.nGenS1ThrChHad++;

    // --- sums (status1+thr only) ---
    out.all.genSum += g->pt();
    if (g->charge()==0) out.ne.genSum += g->pt();
    else                out.ch.genSum += g->pt();

    // --- improved matching: neutral only, pt-compatible, best DR ---
    if (g->charge() != 0) continue;

    //const double rel = (seed->pt()>0 ? std::abs(g->pt()-seed->pt())/seed->pt() : 999.);
    const double rel = (seed->pt()>0 ? (g->pt()-seed->pt())/seed->pt() : 999.);
    if (rel > relPtTol) continue;

    if (dr < bestDR) { bestDR = dr; bestGenNeutral = gp; }

    // --- Matching candidate: GEN neutral (status1+thr already enforced above) ---
    /* //GOOD
       if (g->charge() == 0) {
       if (dr < bestMatchDR) {
       bestMatchDR = dr;
       bestGenNeutral = g;
	}
    }
    */  

    // keep matching logic: "any passing status1 gen within cone" -> matched
    //out.isGenMatched = 1;

    // [1] bestMatchDR
    /*
    if (genIsStatus1(gp) && g->charge()==0 && genPassThr(g)) {
      //std::cout << "dr = " << dr << " and bestMatchDR = " << bestMatchDR << std:: endl;
      if (dr < bestMatchDR) {
	bestMatchDR = dr;
	bestGenMatch = g;
      }
    }
    */

  }

  // finalize match result
  out.hasMatchedGenNeutral = (bestGenNeutral && bestDR < matchDR) ? 1 : 0;
  out.isGenMatched         = out.hasMatchedGenNeutral;

  out.matchedGenNeutralPt    = out.hasMatchedGenNeutral ? bestGenNeutral->pt()    : 0.f;
  out.matchedGenNeutralPdgId = out.hasMatchedGenNeutral ? bestGenNeutral->pdgId() : 0;
  out.matchedGenNeutralDR    = out.hasMatchedGenNeutral ? bestDR                 : -1.f;

  // charged numerator + matched neutral (no double counting)
  out.genChargedPlusMatchedNeutral = out.ch.genSum + out.matchedGenNeutralPt;

  return out;  
}


void
L1PFCandTableProducer::produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const
{
    edm::Handle<reco::CandidateView> src;
    std::vector<const reco::Candidate *> selected;
    std::vector<const reco::Candidate *> gen_selected;
    std::vector<float> vals_pt, vals_eta, vals_phi, vals_mass;
    
    // ---- collect GEN candidates ----
    for (auto & gencands : gencands_) {
        iEvent.getByToken(gencands.src, src);
        for (const auto & j : *src) {
            if (sel_(j) && gencands.sel(j)) gen_selected.push_back(&j);
        }
    }

    // ---- loop over candidate collections ----
    for (auto & cands : cands_) {
        iEvent.getByToken(cands.src, src);
        for (const auto & j : *src) {
            if (sel_(j) && cands.sel(j)) selected.push_back(&j);
        }

        unsigned int ncands = selected.size();
	std::vector<float> vals_caloeta(ncands, 0.f);
	std::vector<float> vals_calophi(ncands, 0.f);
	//std::cout << "---------------------> Selected size=" << selected.size() << std::endl;
	
        auto out = std::make_unique<nanoaod::FlatTable>(ncands, cands.coll+"Cands", false);

        // ---- fill basic info ----
        vals_pt.resize(ncands);
        vals_eta.resize(ncands);
        vals_phi.resize(ncands);
        vals_mass.resize(ncands);
        for (unsigned int i = 0; i < ncands; ++i) {
            vals_pt[i]   = selected[i]->pt();
            vals_eta[i]  = selected[i]->eta();
            vals_phi[i]  = selected[i]->phi();
            vals_mass[i] = selected[i]->mass();
	    //std::cout << "Hello selected[i]->pt() = " << selected[i]->pt() << std::endl;
        }
        out->addColumn<float>("pt",   vals_pt,   "pt of cand");
        out->addColumn<float>("eta",  vals_eta,  "eta of cand");
        out->addColumn<float>("phi",  vals_phi,  "phi of cand");
        out->addColumn<float>("mass", vals_mass, "mass of cand");

        // ---- extra user-defined variables ----
        for (const auto & evar : extraVars_) {
            for (unsigned int i = 0; i < ncands; ++i) vals_pt[i] = evar.func(*selected[i]);
            out->addColumn<float>(evar.name, vals_pt, evar.expr);
        }

        // ---- allocate output vectors ----
        std::vector<int> vals_isGenMatched(ncands, 0);
	std::vector<float> vals_genPtSum0p2(ncands, -99.f);
	std::vector<float> vals_genNeutralPtSum0p2(ncands, -99.f);
	std::vector<float> vals_genChargedPtSum0p2(ncands, -99.f);
	
	std::vector<float> vals_recoPtSum0p2(ncands, -99.f);
	std::vector<float> vals_recoNeutralPtSum0p2(ncands, -99.f);
	std::vector<float> vals_recoChargedPtSum0p2(ncands, -99.f);

	std::vector<float> vals_recoDen_seedPlusCh0p2(ncands, -99.f);
	std::vector<float> vals_recoDen_seedPlusNe0p2(ncands, -99.f);
	std::vector<float> vals_recoDen_seedPlusAll0p2(ncands, -99.f);
	
	std::vector<float> vals_genChPlusMatchedNe0p2(ncands, -99.f);
	
	std::vector<float> vals_ratioPU_ch0p2(ncands, -99.f);
	std::vector<float> vals_ratioPU_ne0p2(ncands, -99.f);
	std::vector<float> vals_ratioPU_all0p2(ncands, -99.f);
	
	std::vector<int>   vals_hasMatchedGenNeutral(ncands, 0);
	std::vector<float> vals_matchedGenNeutralPt(ncands, 0.f);
	
	// ---- multiplicities (0p2) ----
	std::vector<int> nRecoInCone0p2(ncands, -1);
	std::vector<int> nRecoChInCone0p2(ncands, -1);
	std::vector<int> nRecoNeInCone0p2(ncands, -1);
	std::vector<int> nRecoPhoInCone0p2(ncands, -1);
	std::vector<int> nRecoNHadInCone0p2(ncands, -1);
	std::vector<int> nRecoChHadInCone0p2(ncands, -1);
	
	std::vector<int> nGenInCone0p2(ncands, -1);
	std::vector<int> nGenChInCone0p2(ncands, -1);
	std::vector<int> nGenNeInCone0p2(ncands, -1);
	std::vector<int> nGenPhoInCone0p2(ncands, -1);
	std::vector<int> nGenNHadInCone0p2(ncands, -1);
	std::vector<int> nGenChHadInCone0p2(ncands, -1);
	

        // ---- main candidate loop ----
        const float bz = 3.8112;
	const double cone = 0.2;

        for (unsigned int i = 0; i < ncands; ++i) {
	  const auto* cand = selected[i];
	  
	  // always fill calo for completeness
	  math::XYZTLorentzVector vertex(cand->vx(), cand->vy(), cand->vz(), 0.);
	  auto caloetaphi = l1tpf::propagateToCalo(cand->p4(), vertex, cand->charge(), bz);
	  vals_caloeta[i] = caloetaphi.first;
	  vals_calophi[i] = caloetaphi.second;
	  
	  // Only neutral reco seeds
	  if (cand->charge() != 0) {
	    // leave sentinels
	    vals_isGenMatched[i] = 0;
	    continue;
	  }
	  
	  // ---Compute cone around neutral reco seed
	  auto coneRes = computeConeAroundNeutralRecoSeed(cand, selected, gen_selected, cone, bz);

	  const float r_ch  = (coneRes.recoDen_seedPlusCh  > 0) ? (coneRes.genChargedPlusMatchedNeutral / coneRes.recoDen_seedPlusCh) : -1.f;
	  const float r_ne  = (coneRes.recoDen_seedPlusNe  > 0) ? (coneRes.ne.genSum / coneRes.recoDen_seedPlusNe) : -1.f;
	  const float r_all = (coneRes.recoDen_seedPlusAll > 0) ? (coneRes.all.genSum / coneRes.recoDen_seedPlusAll) : -1.f;
	  
	  // ---Debug printout
	  //if (coneRes.hasMatchedGenNeutral ||
	  //    (r_ch  > 0.8 && r_ch  < 1.2) ||
	  //    (r_ne  > 0.8 && r_ne  < 1.2) ||
	  //    (r_all > 0.8 && r_all < 1.2)) {
	  //  printConeDebugFull(iEvent.id(), cand, selected, gen_selected, coneRes, cone, bz);
	  //}
	  
	  // ---- sums (same names where possible) ----
	  vals_genPtSum0p2[i]         = coneRes.all.genSum;
	  vals_genNeutralPtSum0p2[i]  = coneRes.ne.genSum;
	  vals_genChargedPtSum0p2[i]  = coneRes.ch.genSum;
	  
	  vals_recoPtSum0p2[i]        = coneRes.all.recoSum;
	  vals_recoNeutralPtSum0p2[i] = coneRes.ne.recoSum;
	  vals_recoChargedPtSum0p2[i] = coneRes.ch.recoSum;

	  vals_recoDen_seedPlusCh0p2[i]  = coneRes.recoDen_seedPlusCh;
	  vals_recoDen_seedPlusNe0p2[i]  = coneRes.recoDen_seedPlusNe;
	  vals_recoDen_seedPlusAll0p2[i] = coneRes.recoDen_seedPlusAll;
	  
	  vals_genChPlusMatchedNe0p2[i]  = coneRes.genChargedPlusMatchedNeutral;
	  
	  vals_ratioPU_ch0p2[i]  = (coneRes.recoDen_seedPlusCh  > 0) ? (coneRes.genChargedPlusMatchedNeutral / coneRes.recoDen_seedPlusCh)  : -1.f;
	  vals_ratioPU_ne0p2[i]  = (coneRes.recoDen_seedPlusNe  > 0) ? (coneRes.ne.genSum         / coneRes.recoDen_seedPlusNe)  : -1.f;
	  vals_ratioPU_all0p2[i] = (coneRes.recoDen_seedPlusAll > 0) ? (coneRes.all.genSum        / coneRes.recoDen_seedPlusAll) : -1.f;
	  
	  vals_hasMatchedGenNeutral[i] = coneRes.hasMatchedGenNeutral;
	  vals_matchedGenNeutralPt[i]  = coneRes.matchedGenNeutralPt;
 
	  // ---- matching ----
	  vals_isGenMatched[i] = coneRes.isGenMatched;
	  
	  // ---- multiplicities (0p2) ----
	  nRecoInCone0p2[i]     = coneRes.mult.nRecoAll;
	  nRecoChInCone0p2[i]   = coneRes.mult.nRecoCh;
	  nRecoNeInCone0p2[i]   = coneRes.mult.nRecoNe;
	  nRecoPhoInCone0p2[i]  = coneRes.mult.nRecoPho;
	  nRecoNHadInCone0p2[i] = coneRes.mult.nRecoNHad;
	  nRecoChHadInCone0p2[i]= coneRes.mult.nRecoChHad;
	  
	  nGenInCone0p2[i]      = coneRes.mult.nGenAll;
	  nGenChInCone0p2[i]    = coneRes.mult.nGenCh;
	  nGenNeInCone0p2[i]    = coneRes.mult.nGenNe;
	  nGenPhoInCone0p2[i]   = coneRes.mult.nGenPho;
	  nGenNHadInCone0p2[i]  = coneRes.mult.nGenNHad;
	  nGenChHadInCone0p2[i] = coneRes.mult.nGenChHad;
	}
	
	// ---- Add columns for target studies (0p2 cone) ----
	// ---------------------------------------------------
	// ---- sums  ----
	out->addColumn<float>("genPtSum0p2",         vals_genPtSum0p2,         "GEN all (status1+thr) pT sum in dR=0.2 around neutral RECO seed");
	out->addColumn<float>("genNeutralPtSum0p2",  vals_genNeutralPtSum0p2,  "GEN neutral (status1+thr) pT sum in dR=0.2 around neutral RECO seed");
	out->addColumn<float>("genChargedPtSum0p2",  vals_genChargedPtSum0p2,  "GEN charged (status1+thr) pT sum in dR=0.2 around neutral RECO seed");
	
	out->addColumn<float>("recoPtSum0p2",        vals_recoPtSum0p2,        "RECO all pT sum in dR=0.2 around neutral RECO seed");
	out->addColumn<float>("recoNeutralPtSum0p2", vals_recoNeutralPtSum0p2, "RECO neutral pT sum in dR=0.2 around neutral RECO seed");
	out->addColumn<float>("recoChargedPtSum0p2", vals_recoChargedPtSum0p2, "RECO charged pT sum in dR=0.2 around neutral RECO seed");

	out->addColumn<float>("recoDen_seedPlusCh0p2",  vals_recoDen_seedPlusCh0p2,  "Denominator: seed + reco charged in cone");
	out->addColumn<float>("recoDen_seedPlusNe0p2",  vals_recoDen_seedPlusNe0p2,  "Denominator: seed + reco neutral in cone");
	out->addColumn<float>("recoDen_seedPlusAll0p2", vals_recoDen_seedPlusAll0p2, "Denominator: seed + reco all in cone");
	
	out->addColumn<float>("genChPlusMatchedNe0p2",  vals_genChPlusMatchedNe0p2,  "Numerator: gen charged + matched neutral (if any)");
	
	out->addColumn<float>("ratioPU_ch0p2",  vals_ratioPU_ch0p2,  "PU-style ratio: (gen charged + matched neutral) / (seed + reco charged)");
	out->addColumn<float>("ratioPU_ne0p2",  vals_ratioPU_ne0p2,  "PU-style ratio: gen neutral / (seed + reco neutral)");
	out->addColumn<float>("ratioPU_all0p2", vals_ratioPU_all0p2, "PU-style ratio: gen all / (seed + reco all)");
	
	out->addColumn<int>  ("hasMatchedGenNeutral", vals_hasMatchedGenNeutral, "1 if a matched GEN neutral was found");
	out->addColumn<float>("matchedGenNeutralPt",  vals_matchedGenNeutralPt,  "pT of matched GEN neutral (0 if none)");
 
	// ---- matching + calo ----
	out->addColumn<int>("isGenMatched", vals_isGenMatched, "1 if any GEN status1 passing threshold is within dR<0.2 of neutral RECO seed");
	out->addColumn<float>("caloeta", vals_caloeta, "");
	out->addColumn<float>("calophi", vals_calophi, "");
	
	// ---- multiplicities (0p2) ----
	out->addColumn<int>("nRecoInCone0p2",     nRecoInCone0p2,     "");
	out->addColumn<int>("nRecoChInCone0p2",   nRecoChInCone0p2,   "");
	out->addColumn<int>("nRecoNeInCone0p2",   nRecoNeInCone0p2,   "");
	out->addColumn<int>("nRecoPhoInCone0p2",  nRecoPhoInCone0p2,  "");
	out->addColumn<int>("nRecoNHadInCone0p2", nRecoNHadInCone0p2, "");
	out->addColumn<int>("nRecoChHadInCone0p2",nRecoChHadInCone0p2,"");
	
	out->addColumn<int>("nGenInCone0p2",      nGenInCone0p2,      "");
	out->addColumn<int>("nGenChInCone0p2",    nGenChInCone0p2,    "");
	out->addColumn<int>("nGenNeInCone0p2",    nGenNeInCone0p2,    "");
	out->addColumn<int>("nGenPhoInCone0p2",   nGenPhoInCone0p2,   "");
	out->addColumn<int>("nGenNHadInCone0p2",  nGenNHadInCone0p2,  "");
	out->addColumn<int>("nGenChHadInCone0p2", nGenChHadInCone0p2, "");
	
	// ---- Save to event ----
        iEvent.put(std::move(out), cands.coll+"Cands");
        selected.clear();
    }
}

// define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1PFCandTableProducer);

