
#ifndef EGAMMAOBJECTS_GBRTrainer2D
#define EGAMMAOBJECTS_GBRTrainer2D

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// GBRTrainer2D                                                           //
//                                                                      //
// A fast minimal implementation of Gradient-Boosted Regression Trees   //
// which has been especially optimized for size on disk and in memory.  //                                                                  
//                                                                      //
// This class implements a cpu-time optimized training routine          //
// designed especially for regression training with a very large,       //
//  number of input variables and large training samples                //
//                                                                      //
// Expensive loops are optimized to ensure vectorization by the         //
// compiler, and the pre-sorting, as well as  main split-finding        //
//  search is multithreaded using OpenMP at the level of the            //
// loop over input variables                                            //
//                                                                      //
// Split-finding algorithm is a reworked version of the fixed binning   //
// method used in TMVA.  Input variables are binned in terms of         //
// (weighted) quantiles of the training sample.  In this way the binning //
// is robust against outliers in the training sample, and much fewer    //
// bins are needed for the split search on each node.                   //
//                                                                      //
//  Josh Bendavid - MIT                                                 //
//////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>
#include <stdio.h>

  class GBRForest2D;
  class GBRTree2D;
  class TTree;
  class TCut;
  class GBREvent2D;
  
  class GBRTrainer2D {

    public:

       GBRTrainer2D();
       ~GBRTrainer2D();
       
       void AddInputVar(std::string var)    { fInputVars.push_back(var); }
       void SetTargetXVar(std::string var)  { fTargetXVar = var;         }
       void SetTargetYVar(std::string var)  { fTargetYVar = var;         }
       void SetTree(TTree *tree)            { fTree = tree;              }
       void SetTrainingCut(std::string cut) { fTrainingCut = cut;        }
       void SetMinEvents(int n)             { fMinEvents = n;            }
       void SetShrinkage(float x)           { fShrinkage = x;            }
       
       const GBRForest2D *TrainForest(int ntrees);
       
    protected:
      
       //float WeightedMedian(std::vector<GBREvent*> &evts);
      
      void TrainTree(const std::vector<GBREvent2D*> &evts, double sumwtotal, GBRTree2D &tree, int nvars, double transition);      
      void BuildLeaf(const std::vector<GBREvent2D*> &evts, double sumw, GBRTree2D &tree, double transition);
      
      TTree                    *fTree;
      std::string               fTrainingCut;
      std::vector<std::string>  fInputVars;  
      std::string               fTargetXVar;
      std::string               fTargetYVar;
      int                       fMinEvents;
      float                     fShrinkage;
      std::vector<std::vector<float> > fQuantileMaps;
      int                       fNQuantiles;
      unsigned int              fNBinsMax;
      float                     fTransitionQuantile;
      
      float *sepgains;
      float *cutvals;  
      int *nlefts;
      int *nrights;
      float *sumwlefts;
      float *sumwrights;  
      int   *bestbins;
      
      
      float **ws;
      int **ns;
      float **tgtxs;
      float **tgtys;
      float **tgt2s;
      float **sumws;
      int **sumns;
      float **sumtgtxs;
      float **sumtgtys;
      float **sumtgt2s;
      float **varvals;
      float **bsepgains;
      
      int **quants;
      int **bins;
      
      //float *targets;
  
  };
#endif
