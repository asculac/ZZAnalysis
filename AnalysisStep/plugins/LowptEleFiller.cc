/** \class LowptEleFiller
 *
 *  No description available.
 *
 *  $Date: 2013/05/24 15:42:42 $
 *  $Revision: 1.28 $
 *  \author N. Amapane (Torino)
 *  \author S. Bolognesi (JHU)
 *  \author C. Botta (CERN)
 *  \author S. Casasso (Torino)
 */

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/ESHandle.h>

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include <ZZAnalysis/AnalysisStep/interface/CutSet.h>
#include <ZZAnalysis/AnalysisStep/interface/LeptonIsoHelper.h>

#include <vector>
#include <string>
#include "TRandom3.h"

using namespace edm;
using namespace std;
using namespace reco;



class LowptEleFiller : public edm::EDProducer {
 public:
  /// Constructor
  explicit LowptEleFiller(const edm::ParameterSet&);

  /// Destructor
  ~LowptEleFiller();

 private:
  virtual void beginJob(){};
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};

  edm::EDGetTokenT<pat::ElectronRefVector> electronToken;
  edm::EDGetTokenT<pat::ElectronRefVector> pf_electronToken;
  edm::EDGetTokenT<pat::ElectronRefVector> electronToken_bis;
  edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;
  int sampleType;
  int setup;
  const StringCutObjectSelector<pat::Electron, true> cut;
  const CutSet<pat::Electron> flags;
  edm::EDGetTokenT<double> rhoToken;
  edm::EDGetTokenT<vector<Vertex> > vtxToken;
  TRandom3 rgen_;
};


LowptEleFiller::LowptEleFiller(const edm::ParameterSet& iConfig) :
  electronToken(consumes<pat::ElectronRefVector>(iConfig.getParameter<edm::InputTag>("src"))),
  pf_electronToken(consumes<pat::ElectronRefVector>(iConfig.getParameter<edm::InputTag>("src_pf"))),
  vertexSrc_(consumes<reco::VertexCollection> ( iConfig.getParameter<edm::InputTag>("src_vertex") )),
  sampleType(iConfig.getParameter<int>("sampleType")),
  setup(iConfig.getParameter<int>("setup")),
  cut(iConfig.getParameter<std::string>("cut")),
  flags(iConfig.getParameter<ParameterSet>("flags")),
  rgen_(0)
{
  rhoToken = consumes<double>(LeptonIsoHelper::getEleRhoTag(sampleType, setup));
  vtxToken = consumes<vector<Vertex> >(edm::InputTag("goodPrimaryVertices"));
  produces<pat::ElectronCollection>();

}
LowptEleFiller::~LowptEleFiller(){
}


void
LowptEleFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // Get leptons and rho
  edm::Handle<pat::ElectronRefVector> electronHandle;
  iEvent.getByToken(electronToken, electronHandle);

  edm::Handle<pat::ElectronRefVector> pf_electronHandle;
  iEvent.getByToken(pf_electronToken, pf_electronHandle);

  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoToken, rhoHandle);
  double rho = *rhoHandle;

  edm::Handle<vector<Vertex> > vertices;
  iEvent.getByToken(vtxToken,vertices);

  edm::Handle<reco::VertexCollection> vertexHandle;
  iEvent.getByToken(vertexSrc_, vertexHandle);
  const reco::Vertex & PV = vertexHandle->front();

  edm::ESHandle<TransientTrackBuilder> theB ;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  // Output collection
  auto result = std::make_unique<pat::ElectronCollection>();

  for (unsigned int i = 0; i< electronHandle->size(); ++i){
    //---Clone the pat::Electron
    pat::Electron l(*((*electronHandle)[i].get()));


    //--- PF ISO -- not existing for lowpt collection!
    // for cone size R=0.3 :
    float PFChargedHadIso   = l.pfIsolationVariables().sumChargedHadronPt;
    float PFNeutralHadIso   = l.pfIsolationVariables().sumNeutralHadronEt;
    float PFPhotonIso       = l.pfIsolationVariables().sumPhotonEt;
    // cout<<"PFChargedHadIso "<<PFChargedHadIso<<endl;
    // cout<<"PFNeutralHadIso "<<PFNeutralHadIso<<endl;
    // cout<<"PFPhotonIso "<<PFPhotonIso<<endl;
    // cout<<"rho "<<rho<<endl;

    float SCeta = l.superCluster()->eta();
    float fSCeta = fabs(SCeta);

    float combRelIsoPF = LeptonIsoHelper::combRelIsoPF(sampleType, setup, rho, l);

    //--- SIP, dxy, dz
    // float IP      =  fabs(l.dB(pat::Electron::PV3D));
    // float IPError = llow.edB(pat::Electron::PV3D);

    // cout <<"IP: "<<IP<<endl;
    // cout<<"IPError: "<<IPError<<endl;
    // float SIP     = IP/IPError; //FIXME! why is IP inf in low pt ele?
    // // cout<<"SIP: "<< SIP<<endl;
     

    //compute IP for electrons: need transient track
    const reco::GsfTrackRef gsfTrk = l.gsfTrack();
    const reco::TransientTrack eleTT =(*theB).build( l.gsfTrack() );
    // PV3D
    std::pair<bool, Measurement1D>res = IPTools::signedImpactParameter3D(eleTT, GlobalVector(gsfTrk->px(), gsfTrk->py(), gsfTrk->pz()), PV);
    double d0_corr = res.second.value();
    double d0_err = PV.isValid() ? res.second.error() : -1.0;
    l.setDB(d0_corr, d0_err, pat::Electron::PV3D);
    //cout<<"IP: "<<endl<<d0_corr<<endl;
    //cout<<"IP_err: "<<endl<<d0_err<<endl;
    //cout<<"SIP: "<<endl<<d0_corr/d0_err<<endl;
  
    l.addUserFloat("SIP", abs(l.dB(pat::Electron::PV3D)/l.edB(pat::Electron::PV3D)));
    //cout<<"SIP: "<<abs(l.dB(pat::Electron::PV3D)/l.edB(pat::Electron::PV3D))<<endl<<"IP: "<< l.dB(pat::Electron::PV3D)<<"IP_err: "<< l.edB(pat::Electron::PV3D)<<endl;

    float dxy = 999.;
    float dz  = 999.;
    const Vertex* vertex = 0;
    if (vertices->size()>0) {
      vertex = &(vertices->front());
      dxy = fabs(l.gsfTrack()->dxy(vertex->position()));
      dz  = fabs(l.gsfTrack()->dz(vertex->position()));
    }

    float BDT = l.electronID("ID");
 
    l.setP4(reco::Particle::PolarLorentzVector(l.pt(), l.gsfTrack()->etaMode(), l.gsfTrack()->phiMode(), l.mass()));
      //take modes
  //  if (use_regression_for_p4_) {
  //    // pt from regression, eta and phi from gsf track mode
  //    l.setP4(reco::Particle::PolarLorentzVector(l.pt(), l.gsfTrack()->etaMode(), l.gsfTrack()->phiMode(), l.mass()));

  //  }else if(use_gsf_mode_for_p4_) {
  //    l.setP4(reco::Particle::PolarLorentzVector(l.gsfTrack()->ptMode(), l.gsfTrack()->etaMode(), l.gsfTrack()->phiMode(), l.mass()));
  //  } else {
  //    // Fix the mass to the proper one
  //    l.setP4(reco::Particle::PolarLorentzVector(l.pt(), l.eta(), l.phi(), l.mass()));    
  //  }
    float pt = l.pt();
    float electronID= l.electronID("ID");

    float trackIso = l.dr03TkSumPt(); //we are using trackIso instead of combRelIsoPF for low pt electrons

    bool isBDT = true;  // selectedLowElectrons already pass ID 
    bool isLOOSE = true; //flag for lowpt electrons
    bool isPFoverlap =false; //flag for overlaps with PF electrons
    
    //PF overlap cleaning    
    bool clean_out = false;
      for(unsigned int iEle=0; iEle<pf_electronHandle->size(); ++iEle) {
        pat::Electron lpf(*((*pf_electronHandle)[iEle].get()));
        clean_out |= (
                fabs(lpf.vz() - l.vz()) < 0.5 &&
                      reco::deltaR(l.eta(), l.phi(), lpf.eta(), lpf.phi() ) < 0.03   );

    }
    if(clean_out) 
    continue; // we are skipping the low pt electrons which overlap with PF electrons
    //if(clean_out) isPFoverlap =true; // we leave the low pt electron but with a flag on, this is problematic later when building a final candidate!

    // ISO cut 
    if(trackIso/pt>0.35) 
    continue; //  iso cut applied in a poor way

    //-- Missing hit
	  int missingHit;
	  missingHit = l.gsfTrack()->hitPattern().numberOfAllHits(HitPattern::MISSING_INNER_HITS);

    //-- Flag for crack electrons (which use different efficiency SFs)
    bool isCrack = l.isGap();

    //--- Trigger matching
    int HLTMatch = 0; //FIXME


    //-- Scale and smearing corrections are now stored in the miniAOD https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2#Energy_Scale_and_Smearing
    //-- Unchanged in UL implementation TWiki accessed on 27/04
    float uncorrected_pt = l.pt();
    float corr_factor = 1;//l.userFloat("ecalTrkEnergyPostCorr") / l.energy();//get scale/smear correction factor directly from miniAOD
    //scale and smsear electron


    //get all scale uncertainties and their breakdown
    float scale_total_up =1; //l.userFloat("energyScaleUp") / l.energy();
    float scale_stat_up = 1; //l.userFloat("energyScaleStatUp") / l.energy();
    float scale_syst_up = 1; //l.userFloat("energyScaleSystUp") / l.energy();
    float scale_gain_up = 1; //l.userFloat("energyScaleGainUp") / l.energy();
    float scale_total_dn =1; //l.userFloat("energyScaleDown") / l.energy();
    float scale_stat_dn = 1; //l.userFloat("energyScaleStatDown") / l.energy();
    float scale_syst_dn = 1; //l.userFloat("energyScaleSystDown") / l.energy();
    float scale_gain_dn = 1; //l.userFloat("energyScaleGainDown") / l.energy();
//1;
    //get all smearing unc1;e//rtainties and their breakdown
    float sigma_total_up =1; //l.userFloat("energySigmaUp") / l.energy();
    float sigma_rho_up =  1; //l.userFloat("energySigmaRhoUp") / l.energy();
    float sigma_phi_up =  1; //l.userFloat("energySigmaPhiUp") / l.energy();
    float sigma_total_dn =1; //l.userFloat("energySigmaDown") / l.energy();
    float sigma_rho_dn =  1; //l.userFloat("energySigmaRhoDown") / l.energy();
    float sigma_phi_dn =  1; //l.userFloat("energySigmaPhiDown") / l.energy();
//1;

    //--- Embed user variables
    l.addUserFloat("trackIso",trackIso);
    l.addUserFloat("PFChargedHadIso",PFChargedHadIso);
    l.addUserFloat("PFNeutralHadIso",PFNeutralHadIso);
    l.addUserFloat("PFPhotonIso",PFPhotonIso);
    l.addUserFloat("combRelIsoPF",combRelIsoPF);
    l.addUserFloat("SCeta",SCeta);
    l.addUserFloat("rho",rho);
    //l.addUserFloat("SIP",SIP);
    l.addUserFloat("dxy",dxy);
    l.addUserFloat("dz",dz);
    l.addUserFloat("BDT",BDT);
    l.addUserFloat("isBDT",isBDT);
    l.addUserFloat("isLOOSE",isLOOSE);
    l.addUserFloat("isPFoverlap",isPFoverlap);
    l.addUserFloat("isCrack",isCrack);
    l.addUserFloat("electronID",electronID);
    l.addUserFloat("HLTMatch", HLTMatch);
    l.addUserFloat("missingHit", missingHit);
    l.addUserFloat("uncorrected_pt",uncorrected_pt);
    l.addUserFloat("scale_total_up",scale_total_up);
    l.addUserFloat("scale_stat_up",scale_stat_up);
    l.addUserFloat("scale_syst_up",scale_syst_up);
    l.addUserFloat("scale_gain_up",scale_gain_up);
    l.addUserFloat("scale_total_dn",scale_total_dn);
    l.addUserFloat("scale_stat_dn",scale_stat_dn);
    l.addUserFloat("scale_syst_dn",scale_syst_dn);
    l.addUserFloat("scale_gain_dn",scale_gain_dn);
    l.addUserFloat("sigma_total_up",sigma_total_up);
    l.addUserFloat("sigma_total_dn",sigma_total_dn);
    l.addUserFloat("sigma_rho_up",sigma_rho_up);
    l.addUserFloat("sigma_rho_dn",sigma_rho_dn);
    l.addUserFloat("sigma_phi_up",sigma_phi_up);
    l.addUserFloat("sigma_phi_dn",sigma_phi_dn);

    //--- MC parent code
//     MCHistoryTools mch(iEvent);
//     if (mch.isMC()) {
//       int MCParentCode = 0;
//       //      int MCParentCode = mch.getParentCode(&l); //FIXME: does not work on cmg
//       l.addUserFloat("MCParentCode",MCParentCode);
//     }

    //--- Check selection cut. Being done here, flags are not available; but this way we
    //    avoid wasting time on rejected leptons.
    if (!cut(l)) continue;

    //--- Embed flags (ie flags specified in the "flags" pset)
    for(CutSet<pat::Electron>::const_iterator flag = flags.begin(); flag != flags.end(); ++flag) {
      l.addUserFloat(flag->first,int((*(flag->second))(l)));
    }

    result->push_back(l);
  }
  iEvent.put(std::move(result));
}


#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(LowptEleFiller);
