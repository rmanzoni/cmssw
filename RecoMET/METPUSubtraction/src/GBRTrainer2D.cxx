#include "../interface/GBRTrainer2D.h"
#include "../interface/GBREvent2D.h"
#include "CondFormats/EgammaObjects/interface/GBRForest2D.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include <assert.h>


//_______________________________________________________________________
GBRTrainer2D::GBRTrainer2D() : 
  fTree(0),
  fMinEvents(2000),
  fShrinkage(0.1),
  fNQuantiles(std::numeric_limits<unsigned short>::max()+1),
  fNBinsMax(128),
  fTransitionQuantile(0.7),
  sepgains(0),
  ws(0)
{

}

//_______________________________________________________________________
GBRTrainer2D::~GBRTrainer2D() 
{
  fTree = 0;

  //clear arrays
  if (sepgains) {
    delete[] sepgains;
    delete[] cutvals;
    delete[] nlefts;
    delete[] nrights;
    delete[] sumwlefts;
    delete[] sumwrights;
    delete[] bestbins;  
  }
  
  if (ws) {

    for (unsigned int ivar=0; ivar<fInputVars.size(); ++ivar) {
      delete[] ws[ivar];
      delete[] ns[ivar];
      delete[] tgtxs[ivar];
      delete[] tgtys[ivar];
      delete[] tgt2s[ivar];
      
      delete[] sumws[ivar];
      delete[] sumns[ivar];
      delete[] sumtgtxs[ivar];
      delete[] sumtgtys[ivar];
      delete[] sumtgt2s[ivar];
      delete[] varvals[ivar];   
      delete[] bsepgains[ivar];
      
      delete[] quants[ivar];
      delete[] bins[ivar];
    }
    
    delete[] ws;
    delete[] ns;
    delete[] tgtxs;
    delete[] tgtys;
    delete[] tgt2s;
    
    delete[] sumws;
    delete[] sumns;
    delete[] sumtgtxs;
    delete[] sumtgtys;
    delete[] sumtgt2s;
    delete[] varvals;  
    delete[] bsepgains;
    
    delete[] quants;
    delete[] bins;
  }
  
}

//_______________________________________________________________________
const GBRForest2D *GBRTrainer2D::TrainForest(int ntrees)
{
  
  const int nvars = fInputVars.size();
  
  
  //initialize TTreeFormulas to read variables from TTree
  std::vector<TTreeFormula*> inputforms;
  for (std::vector<std::string>::const_iterator it = fInputVars.begin(); 
       it != fInputVars.end(); ++it) {
    inputforms.push_back(new TTreeFormula(it->c_str(),it->c_str(),fTree));
  }
  
  TTreeFormula targetxform(fTargetXVar.c_str(),fTargetXVar.c_str(),fTree);
  TTreeFormula targetyform(fTargetYVar.c_str(),fTargetYVar.c_str(),fTree);
  TTreeFormula cutform(fTrainingCut.c_str(),fTrainingCut.c_str(),fTree);
  
  
  Long64_t nev = 0;  
  
  //loop over tree to count training events with non-zero weight;
  for (Long64_t iev=0; iev<fTree->GetEntries(); ++iev) {
    fTree->LoadTree(iev);
    if (cutform.EvalInstance()!=0.) {
      ++nev;
    }
  }
  
  printf("nev = %i, nvar = %i\n",int(nev),nvars);

  //initialize arrays
  
  sepgains = new float[nvars];
  cutvals = new float[nvars];
  nlefts = new int[nvars];
  nrights = new int[nvars];
  sumwlefts = new float[nvars];
  sumwrights = new float[nvars];
  bestbins = new int[nvars];

  ws = new float*[nvars];
  ns = new int*[nvars];
  tgtxs = new float*[nvars];
  tgtys = new float*[nvars];
  tgt2s = new float*[nvars];  
  sumws = new float*[nvars];
  sumns = new int*[nvars];
  sumtgtxs = new float*[nvars];
  sumtgtys = new float*[nvars];
  sumtgt2s = new float*[nvars];
  varvals = new float*[nvars];    
  bsepgains = new float*[nvars];
  
  quants = new int*[nvars];
  bins = new int*[nvars];
  
  for (int ivar=0; ivar<nvars; ++ivar) {
    ws[ivar] = new float[fNBinsMax];
    ns[ivar] = new int[fNBinsMax];
    tgtxs[ivar] = new float[fNBinsMax];
    tgtys[ivar] = new float[fNBinsMax];
    tgt2s[ivar] = new float[fNBinsMax];  
    sumws[ivar] = new float[fNBinsMax];
    sumns[ivar] = new int[fNBinsMax];
    sumtgtxs[ivar] = new float[fNBinsMax];
    sumtgtys[ivar] = new float[fNBinsMax];
    sumtgt2s[ivar] = new float[fNBinsMax];
    varvals[ivar] = new float[fNBinsMax];  
    bsepgains[ivar] = new float[fNBinsMax];      
    
    quants[ivar] = new int[nev];
    bins[ivar] = new int[nev];
  }
    
  
  std::vector<GBREvent2D*> evts;
  evts.reserve(nev);
  
  double sumw = 0.;
  
  //loop over tree to fill arrays and event vector
  for (Long64_t iev=0; iev<fTree->GetEntries(); ++iev) {
    fTree->LoadTree(iev);
    
    float weight = cutform.EvalInstance();
    
    if (weight==0.) continue; //skip events with 0 weight
    
    sumw += weight;
    
    evts.push_back(new GBREvent2D(nvars));
    GBREvent2D *evt = evts.back();
    evt->SetWeight(weight);
    evt->SetTarget(targetxform.EvalInstance(), targetyform.EvalInstance());
    
    //printf("target = %5f\n",targetform.EvalInstance());
    
    for (int i=0; i<nvars; ++i) {
      evt->SetVar(i,inputforms[i]->EvalInstance());
    }

  }
  
  //map of input variable quantiles to values
  fQuantileMaps.resize(nvars, std::vector<float>(fNQuantiles));
  
  //parallelize building of quantiles for each input variable
  //(sorting of event pointer vector is cpu-intensive)
#pragma omp parallel for
  for (int ivar=0; ivar<nvars; ++ivar) {
    printf("sorting var %i\n",ivar);
        
    std::map<int,float,std::greater<float> > tmpmap;    
    std::vector<GBREvent2D*> evtsvarsort(evts.begin(),evts.end());
    
    std::sort(evtsvarsort.begin(),evtsvarsort.end(),GBRVarCMP(ivar));
    
    double sumwq = 0;
    for (unsigned int iev=0; iev<evtsvarsort.size(); ++iev) {
      sumwq += evtsvarsort[iev]->Weight();
      int quant = int((sumwq/sumw)*(fNQuantiles-1));
      float val = evtsvarsort[iev]->Var(ivar);
    
      //ensure that events with numerically identical values receive the same quantile
      if (iev>0 && val==evtsvarsort[iev-1]->Var(ivar)) quant = evtsvarsort[iev-1]->Quantile(ivar);
    
      evtsvarsort[iev]->SetQuantile(ivar,quant);
    
      tmpmap[quant] = val;
    
    }
    

    for (int i=0; i<fNQuantiles; ++i) {
      std::map<int,float,std::greater<float> >::const_iterator mit = tmpmap.lower_bound(i);
      
      float val;
      if (mit!=tmpmap.end()) val = mit->second;
      else val = -std::numeric_limits<float>::max();
      
      fQuantileMaps[ivar][i] = val;
      
      
    }
    
    
    
  }
    
  //sort events by target and compute median
  //std::sort(evts.begin(),evts.end(),GBRTargetCMP());
  
  //compute vector mean
  //float meansumw = 0.;
  float sumx = 0.;
  float sumy = 0.;
  for (std::vector<GBREvent2D*>::iterator it=evts.begin(); it!=evts.end(); ++it) {
    sumx += (*it)->Weight()*(*it)->TargetX();
    sumy += (*it)->Weight()*(*it)->TargetY();
  }
  float meanx = sumx/sumw;
  float meany = sumy/sumw;
  
  //set initial response and recompute targets
  GBRForest2D *forest = new GBRForest2D;
  forest->SetInitialResponse(meanx,meany);

  for (std::vector<GBREvent2D*>::iterator it=evts.begin(); it!=evts.end(); ++it) {
    (*it)->SetTarget( (*it)->TargetX()-meanx, (*it)->TargetY()-meany );
  }  
  
  //sort by vector mangnitude of the recomputed target and compute transformed target
  //according to huber loss function derivative (cutoff of outliers)
  //Definition for vectors is a bit ad hoc (direction of outliers is kept, but magnitude is truncated
  
  std::sort(evts.begin(),evts.end(),GBRAbsTargetCMP());
  double transumw = 0.;
  float transition = 0.;
  std::vector<GBREvent2D*>::const_iterator transit=evts.begin();
  while(transumw<(fTransitionQuantile*sumw) && transit!=evts.end()) {
    transumw += (*transit)->Weight();
    transition = (*transit)->TargetMag();  
    ++transit;
  } 
  
  for (std::vector<GBREvent2D*>::iterator it=evts.begin(); it!=evts.end(); ++it) {
    float tgtmag = (*it)->TargetMag();
    float tgtx = (*it)->TargetX();
    float tgty = (*it)->TargetY();
    if (tgtmag<transition) {
      (*it)->SetTransTarget(tgtx,tgty);
    }
    else {
      (*it)->SetTransTarget(tgtx*transition/tgtmag, tgty*transition/tgtmag);
    } 
  }

  
  //printf("sumw = %5f, median = %5f, transition = %5f\n",sumw, median,transition);
  
  printf("sumw = %5f, meanx = %5f, meany = %5f, transition = %5f\n",sumw, meanx, meany,transition);
  
  //loop over requested number of trees
  for (int itree=0; itree<ntrees; ++itree) {
    printf("tree %i\n",itree);

    //sort events by recomputed target, which is expected/required for correct computation
    //of median for each terminal mode
    //std::sort(evts.begin(),evts.end(),GBRTargetCMP());
      
    forest->Trees().push_back(GBRTree2D());
    GBRTree2D &tree = forest->Trees().back();

    //train a single tree recursively from the root node
    TrainTree(evts,sumw,tree,nvars,transition);
    
    //recompute transition point and transformed target
    std::sort(evts.begin(),evts.end(),GBRAbsTargetCMP());
    double transumw = 0.;
    float transition = 0.;
    std::vector<GBREvent2D*>::const_iterator transit=evts.begin();
    while(transumw<(fTransitionQuantile*sumw) && transit!=evts.end()) {
      transumw += (*transit)->Weight();
      transition = (*transit)->TargetMag();  
      ++transit;
    } 
    
    for (std::vector<GBREvent2D*>::iterator it=evts.begin(); it!=evts.end(); ++it) {
      float tgtmag = (*it)->TargetMag();
      float tgtx = (*it)->TargetX();
      float tgty = (*it)->TargetY();
      if (tgtmag<transition) {
        (*it)->SetTransTarget(tgtx, tgty);
      }
      else {
        (*it)->SetTransTarget(tgtx*transition/tgtmag, tgty*transition/tgtmag);
      } 
    }
    
  }
  
  //return fully trained GBRForest2D
  return forest;
  
}

//_______________________________________________________________________
void GBRTrainer2D::TrainTree(const std::vector<GBREvent2D*> &evts, double sumwtotal, GBRTree2D &tree, int nvars, double transition) {
  
  //index of current intermediate node
  int thisidx = tree.CutIndices().size();    
  
  //number of events input to node
  int nev = evts.size();
  
  //index of best cut variable
  int bestvar = 0;

  //trivial open-mp based multithreading of loop over input variables
  //The loop is thread safe since each iteration writes into its own
  //elements of the 2-d arrays
#pragma omp parallel for
  for (int ivar=0; ivar<nvars; ++ivar) {

    //fill temporary array of quantiles (to allow auto-vectorization of later loops)
    for (int iev = 0; iev<nev; ++iev) {
      quants[ivar][iev] = evts[iev]->Quantile(ivar);
    }
    
    int minquant = std::numeric_limits<int>::max();
    int maxquant = 0;
    
    //find max and min quantiles in the input events
    //(this loop should be vectorized by gcc with reasonable optimization options)
    for (int iev = 0; iev<nev; ++iev) {
      if (quants[ivar][iev]<minquant) minquant = quants[ivar][iev];
      if (quants[ivar][iev]>maxquant) maxquant = quants[ivar][iev];
    }    
    
    //calculate offset and scaling (powers of 2) to reduce the total number of quantiles
    //to the fNBinsMax for the search for the best split value
    int offset = minquant;
    unsigned int bincount = maxquant-minquant+1;
    unsigned int pscale = 0;
    while (bincount>fNBinsMax) {
      ++pscale;
      bincount >>= 1;
    }    
//    int scale = 1<<pscale;
    
    //final number of bins (guaranteed to be <= fNBinsMax) for best split search
    const unsigned int nbins = ((maxquant-offset)>>pscale) + 1;
    

    //zero arrays where necessary and compute map between bin numbers
    //and variable cut values
    //This loop should auto-vectorize in appropriate compiler/settings
    for (unsigned int ibin=0; ibin<nbins; ++ibin) {
      ws[ivar][ibin] = 0.;
      ns[ivar][ibin] = 0;
      tgtxs[ivar][ibin] = 0.;
      tgtys[ivar][ibin] = 0.;
      tgt2s[ivar][ibin] = 0.;
      
      int quant = ((1+ibin)<<pscale) + offset - 1;
      
      varvals[ivar][ibin] = fQuantileMaps[ivar][quant];

    }
    
    //compute reduced bin value for each event using bit-shift operations
    //This loop should auto-vectorize in appropriate compiler/settings
    for (int iev=0;iev<nev;++iev) {
      bins[ivar][iev] = (quants[ivar][iev]-offset)>>pscale;
    }

     
    //compute summed quantities differential in each bin
    //(filling 'histograms')
    //This loop is one of the most expensive in the algorithm for large training samples
    //This loop can unfortunately not be vectorized because the memory addressed 
    //are computed within the loop iteration
    //JOSH: Is this already fundamentally making vectorization impossible because the addresses to be incremented are
    //scattered, or is it just that the compiler can't resolve the dependencies?  If the latter, can we force gcc to vectorize
    //this loop)
    
    for (int iev=0;iev<nev;++iev) {
      int ibin = bins[ivar][iev];
      
      ws[ivar][ibin] += evts[iev]->Weight();
      ++ns[ivar][ibin];
      tgtxs[ivar][ibin] += evts[iev]->WeightedTransTargetX();
      tgtys[ivar][ibin] += evts[iev]->WeightedTransTargetY();
      tgt2s[ivar][ibin] += evts[iev]->WeightedTransTarget2();

    } 
 
    //convert differential arrays to cumulative arrays by summing over
    //each element
    //loop cannot be vectorized because this is an iterative calculation
    sumws[ivar][0] = ws[ivar][0];
    sumns[ivar][0] = ns[ivar][0];
    sumtgtxs[ivar][0] = tgtxs[ivar][0];
    sumtgtys[ivar][0] = tgtys[ivar][0];
    sumtgt2s[ivar][0] = tgt2s[ivar][0];    
    
    for (unsigned int ibin=1; ibin<nbins; ++ibin) {      
      sumws[ivar][ibin] = sumws[ivar][ibin-1] + ws[ivar][ibin];
      sumns[ivar][ibin] = sumns[ivar][ibin-1] + ns[ivar][ibin];
      sumtgtxs[ivar][ibin] = sumtgtxs[ivar][ibin-1] + tgtxs[ivar][ibin];
      sumtgtys[ivar][ibin] = sumtgtys[ivar][ibin-1] + tgtys[ivar][ibin];
      sumtgt2s[ivar][ibin] = sumtgt2s[ivar][ibin-1] + tgt2s[ivar][ibin];  
    }
    
    //int n = sumns[ivar][nbins-1];
    float sumw = sumws[ivar][nbins-1];
    float sumtgtx = sumtgtxs[ivar][nbins-1];
    float sumtgty = sumtgtys[ivar][nbins-1];
    float sumtgt2 = sumtgt2s[ivar][nbins-1];      
    
    //weighted variance of target in full dataset
    float fullvariance = sumtgt2 - (sumtgtx*sumtgtx + sumtgty*sumtgty)/sumw;
    
   // printf("fullrms = %5f, sumtgt2 = %5f, sumtgt = %5f, sumw = %5f\n",fullrms,sumtgt2,sumtgt,sumw);
    
    //printf("short loop\n");
    float maxsepgain = -std::numeric_limits<float>::max();
    float cutval = 0.;
    int nleft= 0;
    int nright = 0;
    float sumwleft=0.;
    float sumwright=0.;
    int bestbin=0;
    
    //loop over all bins and compute improvement in weighted variance of target for each split
    //This loop is relatively expensive and should auto-vectorize in the appropriate compiler/settings
    for (unsigned int ibin=0; ibin<nbins; ++ibin) {
      float leftvariance = sumtgt2s[ivar][ibin] - (sumtgtxs[ivar][ibin]*sumtgtxs[ivar][ibin] + sumtgtys[ivar][ibin]*sumtgtys[ivar][ibin])/sumws[ivar][ibin];
      float rightsumw = sumw - sumws[ivar][ibin];
      float righttgtxsum = sumtgtx - sumtgtxs[ivar][ibin];
      float righttgtysum = sumtgty - sumtgtys[ivar][ibin];
      float righttgt2sum = sumtgt2 - sumtgt2s[ivar][ibin];
      float rightvariance = righttgt2sum - (righttgtxsum*righttgtxsum + righttgtysum*righttgtysum)/rightsumw;
      
      //weighted improvement in variance from this split
      bsepgains[ivar][ibin] = fullvariance - rightvariance - leftvariance;
    }
    
    //loop over computed variance improvements and select best split, respecting also minimum number of events per node
    //This loop cannot auto-vectorize, at least in gcc 4.6x due to the mixed type conditional, but it's relatively fast
    //in any case
    for (unsigned int ibin=0; ibin<nbins; ++ibin) {   
      if (sumns[ivar][ibin]>=fMinEvents && (nev-sumns[ivar][ibin])>=fMinEvents && bsepgains[ivar][ibin]>maxsepgain) {
	maxsepgain = bsepgains[ivar][ibin];
        bestbin = ibin;
      }
    }
     
    cutval = varvals[ivar][bestbin];
    nleft = sumns[ivar][bestbin];
    nright = nev - nleft;
    sumwleft = sumws[ivar][bestbin];
    sumwright = sumw - sumwleft;        
    
    sepgains[ivar] = maxsepgain;
    cutvals[ivar] = cutval;
    nlefts[ivar] = nleft;
    nrights[ivar] = nright;
    sumwlefts[ivar] = sumwleft;
    sumwrights[ivar] = sumwright;
    bestbins[ivar] = bestbin;
        
  }
  

  
  float globalsepgain = -std::numeric_limits<float>::max();
  for (int ivar=0; ivar<nvars; ++ivar) {
    if (sepgains[ivar]>globalsepgain) {
      globalsepgain = sepgains[ivar];
      bestvar = ivar;
    }
  }    
  
  //if no appropriate split found, make this node terminal
  if (globalsepgain<=0.) {
    //no valid split found, making this node a leaf node
    //printf("globalsepgain = %5f, no valid split\n",globalsepgain);
    BuildLeaf(evts,sumwtotal,tree,transition);
    return;
  }
  
  //fill vectors of event pointers for left and right nodes below this one
  std::vector<GBREvent2D*> leftevts;
  std::vector<GBREvent2D*> rightevts;
  
  leftevts.reserve(nev);
  rightevts.reserve(nev);
  
  int nleft = 0;
  int nright = 0;
  double sumwleft = 0.;
  double sumwright = 0.;
  
  for (std::vector<GBREvent2D*>::const_iterator it = evts.begin(); it!=evts.end(); ++it) {
    if ((*it)->Var(bestvar)>cutvals[bestvar]) {
      ++nright;
      sumwright += (*it)->Weight();
      rightevts.push_back(*it);
    }
    else {
      ++nleft;
      sumwleft += (*it)->Weight();
      leftevts.push_back(*it);
    }    
  }
 
  //printf("thisidx = %i, bestvar = %i, cutval = %5f, n = %i, nleft = %i, nright = %i\n",thisidx,bestvar,cutvals[bestvar],nev,nlefts[bestvar],nrights[bestvar]);
  
  assert(nlefts[bestvar]==nleft);
  assert(nrights[bestvar]==nright);
  
  //fill intermediate node
  tree.CutIndices().push_back(bestvar);
  tree.CutVals().push_back(cutvals[bestvar]);
  tree.LeftIndices().push_back(0);   
  tree.RightIndices().push_back(0);  
  
  //check if left node is terminal
  bool termleft = nleft<=(2*fMinEvents);
  if (termleft) tree.LeftIndices()[thisidx] = -tree.ResponsesX().size();
  else tree.LeftIndices()[thisidx] = tree.CutIndices().size();
  
  //printf("this idx = %i, termleft = %i, nleft = %i, fMinEvents = %i\n",thisidx,  termleft,nleft,fMinEvents);  
  
  //build left node as appropriate
  if (termleft) {  
    BuildLeaf(leftevts,sumwleft,tree,transition);
  }
  else {  
    TrainTree(leftevts,sumwleft,tree,nvars,transition);  
  }
  
  //check if right node is terminal
  bool termright = nright<=(2*fMinEvents);
  if (termright) tree.RightIndices()[thisidx] = -tree.ResponsesX().size();
  else tree.RightIndices()[thisidx] = tree.CutIndices().size();
    
  //printf("this idx = %i, termright = %i, nright = %i, fMinEvents = %i\n",thisidx,  termright,nright,fMinEvents);    
  
  //build right node as appropriate
  if (termright) {  
    BuildLeaf(rightevts,sumwright,tree,transition);
  }
  else {  
    TrainTree(rightevts,sumwright,tree,nvars,transition);  
  }
  
}

  
  


//_______________________________________________________________________
void GBRTrainer2D::BuildLeaf(const std::vector<GBREvent2D*> &evts, double sumw, GBRTree2D &tree, double transition) {

  //printf("building leaf\n");
  
  //int thisidx = -tree.ResponsesX().size();
  //printf("thisidx = %i\n",thisidx);
  
 
  float sumx = 0.;
  float sumy = 0.;
  for (std::vector<GBREvent2D*>::const_iterator it=evts.begin(); it!=evts.end(); ++it) {
    sumx += (*it)->Weight()*(*it)->TargetX();
    sumy += (*it)->Weight()*(*it)->TargetY();
  }
  float meanx = sumx/sumw;
  float meany = sumy/sumw;
  
  float shiftx = 0.;
  float shifty = 0.;
  const float invsumw = 1.0/sumw;
  for (std::vector<GBREvent2D*>::const_iterator it = evts.begin(); it!=evts.end(); ++it) {
    float weight = (*it)->Weight();
    float diffx = (*it)->TargetX() - meanx;
    float diffy = (*it)->TargetY() - meany;
    float diffmag = sqrt(diffx*diffx + diffy*diffy);
    
    
    if (diffmag > transition) {
      diffx *= transition/diffmag;
      diffy *= transition/diffmag;
    }
    
    shiftx += weight*invsumw*diffx; 
    shifty += weight*invsumw*diffy; 
  
  }
  
  float responsex = fShrinkage*(meanx+shiftx);
  float responsey = fShrinkage*(meany+shifty);
  
  tree.ResponsesX().push_back(responsex);
  tree.ResponsesY().push_back(responsex);
  
  for (std::vector<GBREvent2D*>::const_iterator it = evts.begin(); it!=evts.end(); ++it) {
    (*it)->SetTarget((*it)->TargetX()-responsex, (*it)->TargetY()-responsey);
  }
  
  //printf("thisidx = %i, n = %i, responsex = %5f, responsey = %5f\n", thisidx, int(evts.size()) ,responsex, responsey);
  
}


