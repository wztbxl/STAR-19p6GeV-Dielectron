#ifndef STMINITREEMAKER_HH
#define STMINITREEMAKER_HH

/***************************************************************************
 *
 * $Id: StMiniTreeMaker.h 2015/04/09  Exp $ 
 * StMiniTreeMaker - class to produce miniTree for mtd related analysis
 * Author: Shuai Yang
 *--------------------------------------------------------------------------
 *
 ***************************************************************************/

#include "StMaker.h"
#include "StTreeStructure.h"

#include "StThreeVectorF.hh"
#include "TLorentzVector.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"

#include <vector>
#include <map>
#ifndef ST_NO_NAMESPACES
using std::vector;
#endif

class TH1D;
class TH2D;
class TH3D;
class TString;
class TTree;
class TFile;

class StPicoDstMaker;
class StPicoDst;
class StPicoTrack;
class StRefMultCorr;

#if !defined(ST_NO_TEMPLATE_DEF_ARGS) || defined(__CINT__)
typedef vector<Int_t> IntVec;
typedef vector<Double_t> DoubleVec;
#else
typedef vector<Int_t, allocator<Int_t>> IntVec;
typedef vector<Double_t, allocator<Double_t>> DoubleVec;
#endif

class StMiniTreeMaker : public StMaker {
	public:
		StMiniTreeMaker(const Char_t *name = "StMiniTreeMaker");
		~StMiniTreeMaker();

		Int_t    Init();
		Int_t    InitRun(const Int_t runNumber);
		Int_t    Make();
		Int_t    Finish();

		void	 setTriggerIDs(const IntVec triggerids);
		void     setUseDefaultVtx(const Bool_t flag); 
		void     setMaxVtxR(const Double_t max);
		void     setMaxVtxZ(const Double_t max);
		void     setMaxVzDiff(const Double_t max);
		void     setMinTrackPt(const Double_t min);
		void     setMaxTrackEta(const Double_t max);
		void     setMinNHitsFit(const Int_t min);
		void     setMinNHitsFitRatio(const Double_t min);
		void     setMinNHitsDedx(const Int_t min);
		void     setMaxDca(const Double_t max);
		void     setMaxnSigmaE(const Double_t max);
		void     setMaxBeta2TOF(const Double_t max);
		void     setFillHisto(const Bool_t fill);
		void     setFillTree(const Bool_t fill);
		void     setOutFileName(const TString name);
		void     setStreamName(const TString name);
		void     setPrintMemory(const Bool_t pMem);
		void     setPrintCpu(const Bool_t pCpu);
		void     setPrintConfig(const Bool_t print);

	protected:
		void     printConfig();
		void     bookTree();
		void     bookHistos();
		Bool_t   processPicoEvent();
		void     calQxQy(StPicoTrack *pTrack, TVector3 vtexPos) const;
		void     fillEventPlane();
		Bool_t   isValidTrack(StPicoTrack *pTrack, TVector3 vtxPos) const;


	private:
		StPicoDstMaker *mPicoDstMaker;
		StPicoDst      *mPicoDst;
		StRefMultCorr  *refMultCorr; //decide centrality

		Bool_t         mPrintMemory;         // Flag to print out memory usage
		Bool_t         mPrintCpu;            // Flag to print out CPU usage
		Bool_t         mPrintConfig;         // Flag to print out task configuration
		TString        mStreamName;          // Data stream name
		Bool_t         mDefaultVtx;          // Use Default Vertex
		Bool_t         mSelectVtxRank;       // Vertex ranking > 0
		Double_t       mMaxVtxR;             // Maximum vertex r
		Double_t       mMaxVtxZ;             // Maximum vertex z
		Double_t       mMaxVzDiff;           // Maximum VpdVz-TpcVz 
		Double_t       mMinTrkPt;            // Minimum track pt
		Double_t       mMaxTrkEta;           // Maximum track eta
		Int_t          mMinNHitsFit;         // Minimum number of hits used for track fit
		Double_t       mMinNHitsFitRatio;    // Minimum ratio of hits used for track fit
		Int_t          mMinNHitsDedx;        // Minimum number of hits used for de/dx
		Double_t       mMaxDca;              // Maximum track dca
		Double_t       mMaxnSigmaE;          // Maximum nSigmaE cut
		Double_t       mMaxBeta2TOF;         // Maximum |1-1./beta| for TpcE

		Bool_t         mFillHisto;           // Flag of fill the histogram 
		Bool_t         mFillTree;            // Flag of fill the event tree
		TFile          *fOutFile;            // Output file
		TString        mOutFileName;         // Name of the output file 
		StEvtData      mEvtData;
		TTree          *mEvtTree;            // Pointer to the event tree

		IntVec         mTriggerIDs;

		DoubleVec      vEtaPlusCosPart;
		DoubleVec      vEtaPlusSinPart;
		DoubleVec      vEtaPlusPtWeight;
		DoubleVec      vEtaMinusCosPart;
		DoubleVec      vEtaMinusSinPart;
		DoubleVec      vEtaMinusPtWeight;

		//for pile up rejection
		TF1* PileupUplimit;
		TF1* PileupLowlimit;
		TF1* PileupLimit;


		//define histograms ongoing...
		TH1D           *hEvent;
		TH2D           *hVtxYvsVtxX;
		TH2D           *hVPDVzvsTPCVz;
		TH1D           *hVzDiff;
		TH2D           *hGRefMultvsGRefMultCorr;
		TH1D           *hCentrality;
		TH1D           *hVz;
		TH1D           *hVr;

		TH1D           *hRawQx;
		TH1D           *hRawQy;
		TH1D           *hRawEventPlane;

		TH2D           *hdEdxvsP;
		TH2D           *hdNdxvsP;
		TH2D           *hnSigEvsP;
		TH2D           *hBetavsP;
		//histograms for inclusive QA
		TH1D           *hNHitsFit;
		TH1D           *hNHitsPoss;
		TH1D           *hNHitsdEdx;
		TH1D           *hDca;
		TH2D           *hDcavsPt;
		// TH2D           *hBetavsP;
		TH2D           *hBetavsEta;
		TH2D           *hBetavsPhi;
		TH2D           *dEdxvsEta;
		TH2D           *dEdxvsPhi;
		TH2D           *hNSigmaEvsP;
		TH2D           *hNSigmaEvsEta;
		TH2D           *hNSigmaEvsPhi;
		TH2D           *hMSquarevsP;
		TH2D           *hEtavsPhi;
		TH2D           *hEtavsPt;
		TH2D           *hPhivsPt;
		TH2D           *hEEtavsPhi;
		TH2D           *hEEtavsPt;
		TH2D           *hEPhivsPt;
		TH3D           *hEVxvsVyvsVz;
		TH2D           *hnTOFMatchvsRefmult;
		// TH2D           *hEVzvsVx;
		// TH2D           *hEVzvsVy;
		//checking the electron origin  

		ClassDef(StMiniTreeMaker, 1)
};

inline void	StMiniTreeMaker::setTriggerIDs(const IntVec triggerids) { mTriggerIDs = triggerids; }
inline void StMiniTreeMaker::setUseDefaultVtx(const Bool_t flag) { mDefaultVtx = flag; }
inline void StMiniTreeMaker::setMaxVtxR(const Double_t max) { mMaxVtxR = max; }
inline void StMiniTreeMaker::setMaxVtxZ(const Double_t max) { mMaxVtxZ = max; }
inline void StMiniTreeMaker::setMaxVzDiff(const Double_t max) { mMaxVzDiff = max; }
inline void StMiniTreeMaker::setMinTrackPt(const Double_t min){ mMinTrkPt = min;}
inline void StMiniTreeMaker::setMaxTrackEta(const Double_t max){ mMaxTrkEta = max; }
inline void StMiniTreeMaker::setMinNHitsFit(const Int_t min) { mMinNHitsFit = min; }
inline void StMiniTreeMaker::setMinNHitsFitRatio(const Double_t min) { mMinNHitsFitRatio = min; }
inline void StMiniTreeMaker::setMinNHitsDedx(const Int_t min) { mMinNHitsDedx = min; }
inline void StMiniTreeMaker::setMaxDca(const Double_t max) { mMaxDca = max; }
inline void StMiniTreeMaker::setMaxnSigmaE(const Double_t max) { mMaxnSigmaE = max; }
inline void StMiniTreeMaker::setMaxBeta2TOF(const Double_t max) { mMaxBeta2TOF = max; }
inline void StMiniTreeMaker::setFillHisto(const Bool_t fill) { mFillHisto = fill; }
inline void StMiniTreeMaker::setFillTree(const Bool_t fill) { mFillTree = fill; }
inline void StMiniTreeMaker::setOutFileName(const TString name) { mOutFileName = name; }
inline void StMiniTreeMaker::setStreamName(const TString name) { mStreamName = name; }
inline void StMiniTreeMaker::setPrintMemory(const Bool_t pMem) { mPrintMemory = pMem; }
inline void StMiniTreeMaker::setPrintCpu(const Bool_t pCpu) { mPrintCpu = pCpu; }
inline void StMiniTreeMaker::setPrintConfig(const Bool_t print) { mPrintConfig = print; }
#endif
