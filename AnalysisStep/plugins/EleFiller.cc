/** \class EleFiller
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

#include <ZZAnalysis/AnalysisStep/interface/CutSet.h>
#include <ZZAnalysis/AnalysisStep/interface/LeptonIsoHelper.h>

#include <ZZAnalysis/AnalysisStep/interface/MCHistoryTools.h>

#include <vector>
#include <string>
#include "TRandom3.h"

using namespace edm;
using namespace std;
using namespace reco;



class EleFiller : public edm::EDProducer {
 public:
  /// Constructor
  explicit EleFiller(const edm::ParameterSet&);

  /// Destructor
  ~EleFiller();

 private:
  virtual void beginJob(){};
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};

  void findFirstNonElectronMother(const reco::Candidate *particle,
						 int &ancestorPID, int &ancestorStatus);

  edm::EDGetTokenT<pat::ElectronRefVector> electronToken;
  edm::EDGetTokenT<pat::ElectronRefVector> electronToken_bis;
  edm::EDGetTokenT<reco::GenParticleCollection> genToken;
  int sampleType;
  int setup;
  const StringCutObjectSelector<pat::Electron, true> cut;
  const CutSet<pat::Electron> flags;
  edm::EDGetTokenT<double> rhoToken;
  edm::EDGetTokenT<vector<Vertex> > vtxToken;
  TRandom3 rgen_;
};


EleFiller::EleFiller(const edm::ParameterSet& iConfig) :
  electronToken(consumes<pat::ElectronRefVector>(iConfig.getParameter<edm::InputTag>("src"))),
  genToken(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("src_gen"))),
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
EleFiller::~EleFiller(){
}


void
EleFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // Get leptons and rho
  edm::Handle<pat::ElectronRefVector> electronHandle;
  iEvent.getByToken(electronToken, electronHandle);

  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoToken, rhoHandle);
  double rho = *rhoHandle;

  edm::Handle<vector<Vertex> > vertices;
  iEvent.getByToken(vtxToken,vertices);

  edm::Handle<reco::GenParticleCollection>  genParticles;
  iEvent.getByToken(genToken, genParticles);

  // Output collection
  auto result = std::make_unique<pat::ElectronCollection>();

  for (unsigned int i = 0; i< electronHandle->size(); ++i){

    //---Clone the pat::Electron
    pat::Electron l(*((*electronHandle)[i].get()));

    //--- PF ISO
    // for cone size R=0.3 :
    float PFChargedHadIso   = l.pfIsolationVariables().sumChargedHadronPt;
    float PFNeutralHadIso   = l.pfIsolationVariables().sumNeutralHadronEt;
    float PFPhotonIso       = l.pfIsolationVariables().sumPhotonEt;

    float SCeta = l.superCluster()->eta();
    float fSCeta = fabs(SCeta);

    float combRelIsoPF = LeptonIsoHelper::combRelIsoPF(sampleType, setup, rho, l);

    //--- SIP, dxy, dz
    float IP      = fabs(l.dB(pat::Electron::PV3D));
    float IPError = l.edB(pat::Electron::PV3D);
    float SIP     = IP/IPError;

    float dxy = 999.;
    float dz  = 999.;
    const Vertex* vertex = 0;
    if (vertices->size()>0) {
      vertex = &(vertices->front());
      dxy = fabs(l.gsfTrack()->dxy(vertex->position()));
      dz  = fabs(l.gsfTrack()->dz(vertex->position()));
    }


    // Load correct RunII BDT ID+iso
    float BDT = -99;
    if      ( setup == 2016 ) BDT = l.userFloat("ElectronMVAEstimatorRun2Summer16ULIdIsoValues");
    else if ( setup == 2017 ) BDT = l.userFloat("ElectronMVAEstimatorRun2Summer17ULIdIsoValues");
    else if ( setup == 2018 ) BDT = l.userFloat("ElectronMVAEstimatorRun2Summer18ULIdIsoValues");
    // cout << "BDT = " << BDT << endl;

    float pt = l.pt();


    bool isBDT = false;

    if ( setup==2016 )
    {
       //WP taken from https://github.com/asculac/cmssw/blob/Electron_XGBoost_MVA_16UL_17UL/RecoEgamma/ElectronIdentification/python/Identification/mvaElectronID_Summer16UL_ID_ISO_cff.py#L27-L34 and transfered with https://github.com/cms-sw/cmssw/blob/CMSSW_9_4_X/RecoEgamma/EgammaTools/interface/MVAValueMapProducer.h#L145 so that they are between -1 and 1
       isBDT         = (pt<=10 && ((fSCeta<0.8                  && BDT >  0.9557993256) ||
                                   (fSCeta>=0.8 && fSCeta<1.479 && BDT >  0.9475406570) ||
                                   (fSCeta>=1.479               && BDT >  0.9285158721)))
                    || (pt>10  && ((fSCeta<0.8                  && BDT >  0.3272075608) ||
                                   (fSCeta>=0.8 && fSCeta<1.479 && BDT >  0.2468345995) ||
                                   (fSCeta>=1.479               && BDT >  -0.5955762814)));
    }
	 else if (setup==2017)
	 {
	   //WP taken from https://github.com/asculac/cmssw/blob/Electron_XGBoost_MVA_16UL_17UL/RecoEgamma/ElectronIdentification/python/Identification/mvaElectronID_Summer17UL_ID_ISO_cff.py#L27-L34 and transfered with https://github.com/cms-sw/cmssw/blob/CMSSW_9_4_X/RecoEgamma/EgammaTools/interface/MVAValueMapProducer.h#L145 so that they are between -1 and 1
       isBDT         = (pt<=10 && ((fSCeta<0.8                  && BDT >  0.9128577458) ||
                                   (fSCeta>=0.8 && fSCeta<1.479 && BDT >  0.9056792368) ||
                                   (fSCeta>=1.479               && BDT >  0.9439440575)))
                    || (pt>10  && ((fSCeta<0.8                  && BDT >  0.1559788054) ||
                                   (fSCeta>=0.8 && fSCeta<1.479 && BDT >  0.0273863727) ||
                                   (fSCeta>=1.479               && BDT >  -0.5532483665)));
	 }
    else if ( setup==2018 )
    {
       //WP taken from https://github.com/asculac/cmssw/blob/e379e4dd45bb70374cae460a44caff784da92e94/RecoEgamma/ElectronIdentification/python/Identification/mvaElectronID_Summer18UL_ID_ISO_cff.py#L27-L34 and transfered with https://github.com/cms-sw/cmssw/blob/CMSSW_9_4_X/RecoEgamma/EgammaTools/interface/MVAValueMapProducer.h#L145 so that they are between -1 and 1
       isBDT         = (pt<=10 && ((fSCeta<0.8                  && BDT >  0.9044286167) ||
                                   (fSCeta>=0.8 && fSCeta<1.479 && BDT >  0.9094166886) ||
                                   (fSCeta>=1.479               && BDT >  0.9443653660)))
                    || (pt>10  && ((fSCeta<0.8                  && BDT >  0.1968600840) ||
                                   (fSCeta>=0.8 && fSCeta<1.479 && BDT >  0.0759172100) ||
                                   (fSCeta>=1.479               && BDT >  -0.5169136775)));
     }
	 else
	 {
	  	 std::cerr << "[ERROR] EleFiller: no BDT setup for: " << setup << " year!" << std::endl;
	 }

    //-- Missing hit
	 int missingHit;
	 missingHit = l.gsfTrack()->hitPattern().numberOfAllHits(HitPattern::MISSING_INNER_HITS);

    //-- Flag for crack electrons (which use different efficiency SFs)
    bool isCrack = l.isGap();

    //--- Trigger matching
    int HLTMatch = 0; //FIXME

    //dodano za usporedbu iso
    float trackIso = l.dr03TkSumPt();
    //cout<<"track iso"<<endl<<trackIso;


    //-- Scale and smearing corrections are now stored in the miniAOD https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2#Energy_Scale_and_Smearing
    //-- Unchanged in UL implementation TWiki accessed on 27/04
    float uncorrected_pt = l.pt();
    float corr_factor = l.userFloat("ecalTrkEnergyPostCorr") / l.energy();//get scale/smear correction factor directly from miniAOD
    //scale and smsear electron
    l.setP4(reco::Particle::PolarLorentzVector(l.pt()*corr_factor, l.eta(), l.phi(), l.mass()*corr_factor));
    
    int genMatch = 99;
    double dR = 999;
    double charge_reco=l.charge();
    double charge_gen=-99;
    reco::GenParticle const* closestElectron = nullptr;
    for (auto&& particle : *(genParticles.product())) {
      // Drop everything that is not electron or not status 1
      if (abs(particle.pdgId()) != 11 || particle.status() != 1)
        continue;
      //
      double dRtmp = reco::deltaR(l.p4(), particle.p4());
      if (dRtmp < dR) {
        dR = dRtmp;
        closestElectron = &particle;
         charge_gen=particle.charge();
      }
    }

    int ancestorPID = -999; 
    int ancestorStatus = -999;
    findFirstNonElectronMother(closestElectron, ancestorPID, ancestorStatus);
    // See if the closest electron is close enough. If not, no match found.
    if (closestElectron == nullptr || dR >= 0.1){ //PROVJERIT!! bilo dR_ umisto 0.1
      cout<<" UNMATCHED"<<endl;
      genMatch=0;
      }
    
    else if( ancestorPID == -999 && ancestorStatus == -999 ){
      // No non-electron parent??? This should never happen.
      // Complain.
      printf("ElectronNtupler: ERROR! Electron does not apper to have a non-electron parent\n");
      cout<<" UNMATCHED"<<endl;//return UNMATCHED;
      genMatch=0;
    }
    else if( abs(ancestorPID) > 50 ){// && ancestorStatus == 2 )
      cout<<" TRUE_NON_PROMPT_ELECTRON"<<endl;//return TRUE_NON_PROMPT_ELECTRON;
      genMatch=3;
    }
    else  if( abs(ancestorPID) == 23 && ancestorStatus == 22 ){
      cout<<" TRUE_ELECTRON_FROM_Z"<<endl;//return TRUE_ELECTRON_FROM_Z;
      genMatch=1;
    }
    // What remains is true prompt electrons
    else {
      cout<<" TRUE_PROMPT_ELECTRON"<<endl;
      genMatch=2;//return TRUE_PROMPT_ELECTRON;
    }

    //get all scale uncertainties and their breakdown
    float scale_total_up = l.userFloat("energyScaleUp") / l.energy();
    float scale_stat_up = l.userFloat("energyScaleStatUp") / l.energy();
    float scale_syst_up = l.userFloat("energyScaleSystUp") / l.energy();
    float scale_gain_up = l.userFloat("energyScaleGainUp") / l.energy();
    float scale_total_dn = l.userFloat("energyScaleDown") / l.energy();
    float scale_stat_dn = l.userFloat("energyScaleStatDown") / l.energy();
    float scale_syst_dn = l.userFloat("energyScaleSystDown") / l.energy();
    float scale_gain_dn = l.userFloat("energyScaleGainDown") / l.energy();

    //get all smearing uncertainties and their breakdown
    float sigma_total_up = l.userFloat("energySigmaUp") / l.energy();
    float sigma_rho_up = l.userFloat("energySigmaRhoUp") / l.energy();
    float sigma_phi_up = l.userFloat("energySigmaPhiUp") / l.energy();
    float sigma_total_dn = l.userFloat("energySigmaDown") / l.energy();
    float sigma_rho_dn = l.userFloat("energySigmaRhoDown") / l.energy();
    float sigma_phi_dn = l.userFloat("energySigmaPhiDown") / l.energy();



    //--- Embed user variables
    l.addUserFloat("genCharge",charge_gen);
    l.addUserFloat("recoCharge",charge_reco);
    l.addUserFloat("genMatch",genMatch);
    l.addUserFloat("trackIso",trackIso);
    l.addUserFloat("PFChargedHadIso",PFChargedHadIso);
    l.addUserFloat("PFNeutralHadIso",PFNeutralHadIso);
    l.addUserFloat("PFPhotonIso",PFPhotonIso);
    l.addUserFloat("combRelIsoPF",combRelIsoPF);
    l.addUserFloat("SCeta",SCeta);
    l.addUserFloat("rho",rho);
    l.addUserFloat("SIP",SIP);
    l.addUserFloat("dxy",dxy);
    l.addUserFloat("dz",dz);
    l.addUserFloat("BDT",BDT);
    l.addUserFloat("isBDT",isBDT);
    l.addUserFloat("isCrack",isCrack);
    l.addUserFloat("HLTMatch", HLTMatch);
    //l.addUserCand("MCMatch",genMatch); // FIXME
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
    // MCHistoryTools mch(iEvent);
    // if (mch.isMC()) {
    //   int MCParentCode = 0;
    //   int MCParentCode = mch.getParentCode(&l); //FIXME: does not work on cmg
    //   cout<<"Parent code"<<MCParentCode<<endl;
    //   //l.addUserFloat("MCParentCode",MCParentCode);
    // }

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


void EleFiller::findFirstNonElectronMother(const reco::Candidate *particle,
						 int &ancestorPID, int &ancestorStatus){

  if( particle == 0 ){
    printf("ElectronNtupler: ERROR! null candidate pointer, this should never happen\n");
    return;
  }

  // Is this the first non-electron parent? If yes, return, otherwise go deeper into recursion
  if( abs(particle->pdgId()) == 11 || (abs(particle->pdgId()) == 23 && (particle->status())!=22)){
    //std::cout << "Electron Ancestor " << particle->mother(0)->pdgId() << " Status " << particle->mother(0)->status() <<std::endl; 
    findFirstNonElectronMother(particle->mother(0), ancestorPID, ancestorStatus);
  }
  else{
    ancestorPID = particle->pdgId();
    ancestorStatus = particle->status();
    //std::cout << "Particle Id " << ancestorPID << " Status " << ancestorStatus <<std::endl;
  }
  return;
}
#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(EleFiller);
