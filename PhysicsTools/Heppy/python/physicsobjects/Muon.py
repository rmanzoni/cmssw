from PhysicsTools.Heppy.physicsobjects.Lepton import Lepton
from PhysicsTools.HeppyCore.utils.deltar import deltaR

class Muon( Lepton ):

    def __init__(self, *args, **kwargs):
        super(Muon, self).__init__(*args, **kwargs)
        self._trackForDxyDz = "muonBestTrack"

    def setTrackForDxyDz(self,what):
        if not hasattr(self,what):
            raise RuntimeError("I don't have a track called "+what)
        self._trackForDxyDz = what

    def looseId( self ):
        '''Loose ID as recommended by mu POG.'''
        return self.physObj.isLooseMuon()

    def tightId( self ):
        '''Tight ID as recommended by mu POG 
        (unless redefined in the lepton analyzer).

        If not using the LeptonAnalyzer, make sure to set self.associatedVertex, 
        that is necessary for tight muon identification. 
        '''
        return getattr(self,"tightIdResult",self.muonID("POG_ID_Tight"))

    def muonID(self, name, vertex=None):
        if name == "" or name is None: 
            return True
        if name.startswith("POG_"):
            if name == "POG_ID_Loose": return self.physObj.isLooseMuon()
            if vertex is None:
                vertex = getattr(self, 'associatedVertex', None)
            if name == "POG_ID_Tight":  return self.physObj.isTightMuon(vertex)
            if name == "POG_ID_HighPt": return self.physObj.isHighPtMuon(vertex)
            if name == "POG_ID_Soft":   return self.physObj.isSoftMuon(vertex)
            if name == "POG_ID_Soft_ICHEP":
                if not self.physObj.muonID("TMOneStationTight"): return False
                if not self.physObj.innerTrack().hitPattern().trackerLayersWithMeasurement() > 5: return False
                if not self.physObj.innerTrack().hitPattern().pixelLayersWithMeasurement() > 0: return False
                if not abs(self.physObj.innerTrack().dxy(vertex.position())) < 0.3 : return False
                if not abs(self.physObj.innerTrack().dz(vertex.position())) < 20 : return False
                return True
            if name == "POG_ID_TightNoVtx":  return self.looseId() and \
                                                 self.isGlobalMuon() and \
                                                 self.globalTrack().normalizedChi2() < 10 and \
                                                 self.globalTrack().hitPattern().numberOfValidMuonHits() > 0 and \
                                                 self.numberOfMatchedStations()>1 and \
                                                 self.innerTrack().hitPattern().numberOfValidPixelHits()>0 and \
                                                 self.innerTrack().hitPattern().trackerLayersWithMeasurement() > 5
            if name == "POG_ID_Medium":
                if not self.looseId(): return False
                goodGlb = self.physObj.isGlobalMuon() and self.physObj.globalTrack().normalizedChi2() < 3 and self.physObj.combinedQuality().chi2LocalPosition < 12 and self.physObj.combinedQuality().trkKink < 20;
                return self.physObj.innerTrack().validFraction() > 0.8 and self.physObj.segmentCompatibility() >= (0.303 if goodGlb else 0.451)
            if name == "POG_ID_Medium_ICHEP":
                #validFraction() > 0.49 changed from 0.8
                if not self.looseId(): return False
                goodGlb = self.physObj.isGlobalMuon() and self.physObj.globalTrack().normalizedChi2() < 3 and self.physObj.combinedQuality().chi2LocalPosition < 12 and self.physObj.combinedQuality().trkKink < 20;
                return self.physObj.innerTrack().validFraction() > 0.49 and self.physObj.segmentCompatibility() >= (0.303 if goodGlb else 0.451)
            if name == "POG_ID_Medium_Moriond":
                # Allow the run-dependent treatment from the muonPOG recommendation
                # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#Short_Term_Instructions_for_Mori
                # muon should have the muon.event attribute to determine the run'''
                run = self.event.eventAuxiliary().id().run()
                if run > 273016 and run < 278820:
                    return self.muonID("POG_ID_Medium_ICHEP")
                else:
                    return self.muonID("POG_ID_Medium")
                
            if name == "POG_Global_OR_TMArbitrated":
                return self.physObj.isGlobalMuon() or (self.physObj.isTrackerMuon() and self.physObj.numberOfMatchedStations() > 0)
        elif name.startswith("HZZ_"):
            if name == "HZZ_ID_TkHighPt":
                primaryVertex = vertex if vertex != None else getattr(self, 'associatedVertex', None) 
                return ( self.physObj.numberOfMatchedStations() > 1 
                         and (self.physObj.muonBestTrack().ptError()/self.physObj.muonBestTrack().pt()) < 0.3 
                         and abs(self.physObj.muonBestTrack().dxy(primaryVertex.position())) < 0.2 
                         and abs(self.physObj.muonBestTrack().dz(primaryVertex.position())) < 0.5 
                         and self.physObj.innerTrack().hitPattern().numberOfValidPixelHits() > 0 
                         and self.physObj.innerTrack().hitPattern().trackerLayersWithMeasurement() > 5 )
            if name == "HZZ_ID_LooseOrTkHighPt":
                if self.physObj.isLooseMuon(): return True
                return self.physObj.pt() > 200 and self.muonID("HZZ_ID_TkHighPt")
        return self.physObj.muonID(name)

    def mvaId(self):
        '''For a transparent treatment of electrons and muons. Returns -99'''
        return -99
    

    def dxy(self, vertex=None):
        '''either pass the vertex, or set associatedVertex before calling the function.
        note: the function does not work with standalone muons as innerTrack
        is not available.
        '''
        if vertex is None:
            vertex = self.associatedVertex
        return getattr(self,self._trackForDxyDz)().dxy( vertex.position() )

    def edxy(self):
        '''returns the uncertainty on dxy (from gsf track)'''
        return getattr(self,self._trackForDxyDz)().dxyError()
 

    def dz(self, vertex=None):
        '''either pass the vertex, or set associatedVertex before calling the function.
        note: the function does not work with standalone muons as innerTrack
        is not available.
        '''
        if vertex is None:
            vertex = self.associatedVertex
        return getattr(self,self._trackForDxyDz)().dz( vertex.position() )

    def edz(self):
        '''returns the uncertainty on dxz (from gsf track)'''
        return getattr(self,self._trackForDxyDz)().dzError()

    def chargedHadronIso(self,R=0.4):
        if   R == 0.3: return self.physObj.pfIsolationR03().sumChargedHadronPt 
        elif R == 0.4: return self.physObj.pfIsolationR04().sumChargedHadronPt 
        raise RuntimeError("Muon chargedHadronIso missing for R=%s" % R)

    def neutralHadronIso(self,R=0.4):
        if   R == 0.3: return self.physObj.pfIsolationR03().sumNeutralHadronEt 
        elif R == 0.4: return self.physObj.pfIsolationR04().sumNeutralHadronEt 
        raise RuntimeError("Muon neutralHadronIso missing for R=%s" % R)

    def photonIso(self,R=0.4):
        if   R == 0.3: return self.physObj.pfIsolationR03().sumPhotonEt 
        elif R == 0.4: return self.physObj.pfIsolationR04().sumPhotonEt 
        raise RuntimeError("Muon photonIso missing for R=%s" % R)

    def chargedAllIso(self,R=0.4):
        if   R == 0.3: return self.physObj.pfIsolationR03().sumChargedParticlePt 
        elif R == 0.4: return self.physObj.pfIsolationR04().sumChargedParticlePt 
        raise RuntimeError("Muon chargedAllIso missing for R=%s" % R)

    def puChargedHadronIso(self,R=0.4):
        if   R == 0.3: return self.physObj.pfIsolationR03().sumPUPt 
        elif R == 0.4: return self.physObj.pfIsolationR04().sumPUPt 
        raise RuntimeError("Muon chargedHadronIso missing for R=%s" % R)


    def absIsoWithFSR(self, R=0.4, puCorr="deltaBeta", dBetaFactor=0.5):
        '''
        Calculate Isolation, subtract FSR, apply specific PU corrections" 
        '''
        photonIso = self.photonIso(R)
        if hasattr(self,'fsrPhotons'):
            for gamma in self.fsrPhotons:
                dr = deltaR(gamma.eta(), gamma.phi(), self.physObj.eta(), self.physObj.phi())
                if dr > 0.01 and dr < R:
                    photonIso = max(photonIso-gamma.pt(),0.0)                
        if puCorr == "deltaBeta":
            offset = dBetaFactor * self.puChargedHadronIso(R)
        elif puCorr == "rhoArea":
            offset = self.rho*getattr(self,"EffectiveArea"+(str(R).replace(".","")))
        elif puCorr in ["none","None",None]:
            offset = 0
        else:
             raise RuntimeError("Unsupported PU correction scheme %s" % puCorr)
        return self.chargedHadronIso(R)+max(0.,photonIso+self.neutralHadronIso(R)-offset)            

    def ptErr(self):
        if "_ptErr" in self.__dict__: return self.__dict__['_ptErr']
        return self.bestTrack().ptError()

    ###################################
               ## HNL IDs ##
    ###################################

    def Medium(self):
        isGoodGlobal = self.isGlobalMuon()                           and \
                       self.combinedQuality().chi2LocalPosition < 12 and \
                       self.combinedQuality().trkKink < 20
        
        isMedium = 0
        if self.isLooseMuon() and self.segmentCompatibility() > (0.303 if isGoodGlobal else 0.451):
            isMedium = 1
        
        return isMedium

    ''' 
    def MediumID(self):

        Medium = 0
        bool goodGlob = False
        goodGlob = _lGlobalMuon[i] & _lCQChi2Position[i] < 12 & _lCQTrackKink[i] < 20;
        Medium = _lPOGLoose[i] && _muonSegComp[i] > (goodGlob ? 0.303 : 0.451);


        # time veto
        _passTimingVeto[i] = true;
        bool cmbok =( _lMuTimenDof[i] >7 );
        bool rpcok =( _lMuRPCTimenDof[i] >1 & _lMuRPCTimeErr[i]==0 );
            if (rpcok) {
                if ( (fabs(_lMuRPCTime[i])>10) && !(cmbok && fabs(_lMuTime[i])<10) )
                    _passTimingVeto[i]=false;
          } else { if  (cmbok && ( _lMuTime[i]>20 || _lMuTime[i]<-45) ) 
                    _passTimingVeto[i]=false; }
        
        return Medium
    '''
