//===----------------- ModCriticalPath.cpp ----------------------===//
// This file implements the Scaffold Pass of counting the number 
//  of critical timesteps and gate parallelism in program
//  in callgraph post-order.
//
//        This file was created by Scaffold Compiler Working Group
//
//===----------------------------------------------------------------------===//

#define DEBUG_TYPE "ModCriticalPath"
#include <vector>
#include <limits>
#include "llvm/Pass.h"
#include "llvm/Function.h"
#include "llvm/Module.h"
#include "llvm/BasicBlock.h"
#include "llvm/Instruction.h"
#include "llvm/Instructions.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/ADT/Statistic.h"
#include "llvm/Support/InstIterator.h"
#include "llvm/PassAnalysisSupport.h"
#include "llvm/Analysis/CallGraph.h"
#include "llvm/Support/CFG.h"
#include "llvm/ADT/SCCIterator.h"
#include "llvm/Argument.h"
#include "llvm/ADT/ilist.h"
#include "llvm/Constants.h"
#include "llvm/IntrinsicInst.h"


using namespace llvm;
using namespace std;

#define MAX_GATE_ARGS 30
#define MAX_BT_COUNT 15 //max backtrace allowed - to avoid infinite recursive loops
#define NUM_QGATES 17
#define _CNOT 0
#define _H 1
#define _S 2
#define _T 3
#define _X 4
#define _Y 5
#define _Z 6
#define _MeasX 7
#define _MeasZ 8
#define _PrepX 9
#define _PrepZ 10
#define _Tdag 11
#define _Sdag 12
#define _Rz 13
#define _Toffoli 14
#define _Fredkin 15
#define _All 16

bool debugModCriticalPath = false;

namespace {
  
  struct qGateArg{ //arguments to qgate calls
    Value* argPtr;
    int argNum;
    bool isQbit;
    bool isCbit;
    bool isUndef;
    bool isPtr;
    int valOrIndex; //Value if not Qbit, Index if Qbit & not a Ptr
    double angle;
    qGateArg(): argPtr(NULL), argNum(-1), isQbit(false), isCbit(false), isUndef(false), isPtr(false), valOrIndex(-1), angle(0.0){ }
  };
  
struct qArgInfo{
  string name;
  int index;
  qArgInfo(): name("none"), index(-1){ }
};

struct qGate{
  Function* qFunc;
  int numArgs;
  qArgInfo args[MAX_GATE_ARGS];
  double angle;
  qGate():qFunc(NULL), numArgs(0), angle(0.0) { }
};

  struct ArrParGates{
    uint64_t parallel_gates[NUM_QGATES];
    vector<double> RzAngles;
  };

  struct allTSParallelism{
    uint64_t timesteps;
    vector<ArrParGates> gates;
    allTSParallelism(): timesteps(0){ }
  };  

  struct MaxInfo{ //TimeStepInfo
    uint64_t timesteps;
    uint64_t parallel_gates[NUM_QGATES];
    MaxInfo(): timesteps(0){ }
  };

  struct ModCriticalPath : public ModulePass {
    static char ID; // Pass identification
    
    string gate_name[NUM_QGATES];
    vector<qGateArg> tmpDepQbit;
    vector<Value*> vectQbit;
    
    int btCount; //backtrace count

    vector<qArgInfo> currTimeStep; //contains set of arguments operated on currently
    vector<string> currParallelFunc;
    MaxInfo maxParallelFactor; //overall max parallel factor
    map<string, int> gate_index;    
    vector<uint64_t> curr_parallel_ts; //vector of current timesteps that are parallel; used for comparing functions
    map<string, allTSParallelism > funcParallelFactor; //string is function name
    map<string, MaxInfo> funcMaxParallelFactor;

    map<string, map<int,uint64_t> > funcQbits; //qbits in current function
    map<Function*, map<unsigned int, map<int,uint64_t> > > tableFuncQbits;
    map<string, unsigned int> funcArgs;

    map<Function*, uint64_t> funcCritPath;
    allTSParallelism currTS;

    bool isFirstMeas;
       
       
    ModCriticalPath() : ModulePass(ID) {}
    
    bool backtraceOperand(Value* opd, int opOrIndex);
    void analyzeAllocInst(Function* F,Instruction* pinst);
    void analyzeCallInst(Function* F,Instruction* pinst);
    void getFunctionArguments(Function *F);
    
    void saveTableFuncQbits(Function* F);
    void print_tableFuncQbits();

    void init_gate_names(){
        gate_name[_CNOT] = "CNOT";
        gate_name[_H] = "H";
        gate_name[_S] = "S";
        gate_name[_T] = "T";
        gate_name[_Toffoli] = "Toffoli";
        gate_name[_X] = "X";
        gate_name[_Y] = "Y";
        gate_name[_Z] = "Z";
        gate_name[_MeasX] = "MeasX";
        gate_name[_MeasZ] = "MeasZ";
        gate_name[_PrepX] = "PrepX";
        gate_name[_PrepZ] = "PrepZ";
        gate_name[_Sdag] = "Sdag";
        gate_name[_Tdag] = "Tdag";
        gate_name[_Fredkin] = "Fredkin";
        gate_name[_Rz] = "Rz";
        gate_name[_All] = "All";                    
        
        gate_index["CNOT"] = _CNOT;        
        gate_index["H"] = _H;
        gate_index["S"] = _S;
        gate_index["T"] = _T;
        gate_index["Toffoli"] = _Toffoli;
        gate_index["X"] = _X;
        gate_index["Y"] = _Y;
        gate_index["Z"] = _Z;
        gate_index["Sdag"] = _Sdag;
        gate_index["Tdag"] = _Tdag;
        gate_index["MeasX"] = _MeasX;
        gate_index["MeasZ"] = _MeasZ;
        gate_index["PrepX"] = _PrepX;
        gate_index["PrepZ"] = _PrepZ;
        gate_index["Fredkin"] = _Fredkin;
        gate_index["Rz"] = _Rz;
        gate_index["All"] = _All;                    
        
        
        }
        

    void init_gates_as_functions();    
    void init_critical_path_algo(Function* F);
    void calc_critical_time(Function* F, qGate qg);        
    void print_funcQbits();
    void print_qgate(qGate qg);
    void print_critical_info(string func);
    void calc_max_parallelism_statistic();

    void update_critical_info(string currFunc, uint64_t ts, string fname, double angle);

    void print_scheduled_gate(qGate qg, uint64_t ts);

    uint64_t find_max_funcQbits();
    void memset_funcQbits(uint64_t val);

    void print_qgateArg(qGateArg qg)
    {
      errs()<< "Printing QGate Argument:\n";
      if(qg.argPtr) errs() << "  Name: "<<qg.argPtr->getName()<<"\n";
      errs() << "  Arg Num: "<<qg.argNum<<"\n"
	     << "  isUndef: "<<qg.isUndef
	     << "  isQbit: "<<qg.isQbit
	     << "  isCbit: "<<qg.isCbit
	     << "  isPtr: "<<qg.isPtr << "\n"
	     << "  Value or Index: "<<qg.valOrIndex<<"\n";
    }                    
    
    void CountCriticalFunctionResources (Function *F);
    
    bool runOnModule (Module &M);
    
    
    virtual void getAnalysisUsage(AnalysisUsage &AU) const {
      AU.setPreservesAll();  
      AU.addRequired<CallGraph>();    
    }
    
  }; // End of struct ModCriticalPath
} // End of anonymous namespace



char ModCriticalPath::ID = 0;
static RegisterPass<ModCriticalPath> X("ModCriticalPath", "Get Critical Path with strict modular info");

void ModCriticalPath::getFunctionArguments(Function* F)
{
  for(Function::arg_iterator ait=F->arg_begin();ait!=F->arg_end();++ait)
    {    
      //if(ait) errs() << "Argument: "<<ait->getName()<< " ";

      string argName = (ait->getName()).str();
      Type* argType = ait->getType();
      unsigned int argNum=ait->getArgNo();         

      qGateArg tmpQArg;
      tmpQArg.argPtr = ait;
      tmpQArg.argNum = argNum;

      if(argType->isPointerTy()){
	tmpQArg.isPtr = true;

	Type *elementType = argType->getPointerElementType();
	if (elementType->isIntegerTy(16)){ //qbit*
	  tmpQArg.isQbit = true;
	  vectQbit.push_back(ait);
	  
	  map<int,uint64_t> tmpMap;
	  tmpMap[-1] = 0; //add entry for entire array
	  tmpMap[-2] = 0; //add entry for max
	  funcQbits[argName]=tmpMap;

	  funcArgs[argName] = argNum;

	}
	else if (elementType->isIntegerTy(1)){ //cbit*
	  tmpQArg.isCbit = true;
	  //vectQbit.push_back(ait);
	  //funcArgs[argName] = argNum;
	}
      }
      else if (argType->isIntegerTy(16)){ //qbit
	tmpQArg.isQbit = true;
	vectQbit.push_back(ait);

	  map<int,uint64_t> tmpMap;
	  tmpMap[-1] = 0; //add entry for entire array
	  tmpMap[-2] = 0; //add entry for max
	  funcQbits[argName]=tmpMap;
	  funcArgs[argName] = argNum;
      }
      else if (argType->isIntegerTy(1)){ //cbit
	tmpQArg.isCbit = true;
	//vectQbit.push_back(ait);
	//funcArgs[argName] = argNum;
      }
      
    }
}

bool ModCriticalPath::backtraceOperand(Value* opd, int opOrIndex)
{
  if(opOrIndex == 0) //backtrace for operand
    {
      //search for opd in qbit/cbit vector
      vector<Value*>::iterator vIter=find(vectQbit.begin(),vectQbit.end(),opd);
      if(vIter != vectQbit.end()){
	tmpDepQbit[0].argPtr = opd;
	
	return true;
      }
      
      if(btCount>MAX_BT_COUNT)
	return false;
      
      if(GetElementPtrInst *GEPI = dyn_cast<GetElementPtrInst>(opd))
	{

	  if(GEPI->hasAllConstantIndices()){
	    Instruction* pInst = dyn_cast<Instruction>(opd);
	    unsigned numOps = pInst->getNumOperands();

	    backtraceOperand(pInst->getOperand(0),0);
	    
	    //NOTE: getelemptr instruction can have multiple indices. Currently considering last operand as desired index for qubit. Check this reasoning. 
	    if(ConstantInt *CI = dyn_cast<ConstantInt>(pInst->getOperand(numOps-1))){
	      if(tmpDepQbit.size()==1){
		tmpDepQbit[0].valOrIndex = CI->getZExtValue();
	      }
	    }
	  }
	  
	  else if(GEPI->hasIndices()){
	    
	    Instruction* pInst = dyn_cast<Instruction>(opd);
	    unsigned numOps = pInst->getNumOperands();
	    backtraceOperand(pInst->getOperand(0),0);

	    if(tmpDepQbit[0].isQbit && !(tmpDepQbit[0].isPtr)){     
	      //NOTE: getelemptr instruction can have multiple indices. consider last operand as desired index for qubit. Check if this is true for all.
	      backtraceOperand(pInst->getOperand(numOps-1),1);
	      
	    }
	  }
	  else{	    
	    Instruction* pInst = dyn_cast<Instruction>(opd);
	    unsigned numOps = pInst->getNumOperands();
	    for(unsigned iop=0;iop<numOps;iop++){
	      backtraceOperand(pInst->getOperand(iop),0);
	    }
	  }
	  return true;
	}
      
      if(Instruction* pInst = dyn_cast<Instruction>(opd)){
	unsigned numOps = pInst->getNumOperands();
	for(unsigned iop=0;iop<numOps;iop++){
	  btCount++;
	  backtraceOperand(pInst->getOperand(iop),0);
	  btCount--;
	}
	return true;
      }
      else{
	return true;
      }
    }
  else if(opOrIndex == 0){ //opOrIndex == 1; i.e. Backtracing for Index    
    if(btCount>MAX_BT_COUNT) //prevent infinite backtracing
      return true;

    if(ConstantInt *CI = dyn_cast<ConstantInt>(opd)){
      tmpDepQbit[0].valOrIndex = CI->getZExtValue();
      return true;
    }      

    if(Instruction* pInst = dyn_cast<Instruction>(opd)){
      unsigned numOps = pInst->getNumOperands();
      for(unsigned iop=0;iop<numOps;iop++){
	btCount++;
	backtraceOperand(pInst->getOperand(iop),1);
	btCount--;
      }
    }

  }
  else{ //opOrIndex == 2: backtracing to call inst MeasZ
    if(CallInst *endCI = dyn_cast<CallInst>(opd)){
      if(endCI->getCalledFunction()->getName().find("llvm.Meas") != string::npos){
	tmpDepQbit[0].argPtr = opd;

	return true;
      }
      else{
	if(Instruction* pInst = dyn_cast<Instruction>(opd)){
	  unsigned numOps = pInst->getNumOperands();
	  bool foundOne=false;
	  for(unsigned iop=0;(iop<numOps && !foundOne);iop++){
	    btCount++;
	    foundOne = foundOne || backtraceOperand(pInst->getOperand(iop),2);
	    btCount--;
	  }
	  return foundOne;
	}
      }
    }
    else{
      if(Instruction* pInst = dyn_cast<Instruction>(opd)){
	unsigned numOps = pInst->getNumOperands();
	bool foundOne=false;
	for(unsigned iop=0;(iop<numOps && !foundOne);iop++){
	  btCount++;
	  foundOne = foundOne || backtraceOperand(pInst->getOperand(iop),2);
	  btCount--;
	}
	return foundOne;
      }
    }
  }
  return false;
}


void ModCriticalPath::analyzeAllocInst(Function* F, Instruction* pInst){
  if (AllocaInst *AI = dyn_cast<AllocaInst>(pInst)) {
    Type *allocatedType = AI->getAllocatedType();
    
    if(ArrayType *arrayType = dyn_cast<ArrayType>(allocatedType)) {      
      qGateArg tmpQArg;
      
      Type *elementType = arrayType->getElementType();
      uint64_t arraySize = arrayType->getNumElements();
      if (elementType->isIntegerTy(16)){
	vectQbit.push_back(AI);
	tmpQArg.isQbit = true;
	tmpQArg.argPtr = AI;
	tmpQArg.valOrIndex = arraySize;

	map<int,uint64_t> tmpMap; //add qbit to funcQbits
	tmpMap[-1] = 0; //entry for entire array
	tmpMap[-2] = 0; //entry for max
	funcQbits[AI->getName()]=tmpMap;

      }
      
      if (elementType->isIntegerTy(1)){
	vectQbit.push_back(AI); //Cbit added here
	tmpQArg.isCbit = true;
	tmpQArg.argPtr = AI;
	tmpQArg.valOrIndex = arraySize;
      }
    }
  }
}


void ModCriticalPath::init_critical_path_algo(Function* F){

  MaxInfo structMaxInfo;
  allTSParallelism vectGateInfo; //initialize an entry in the function map
  //ArrParGates tmpParGates;
  

  isFirstMeas = true;

  currTimeStep.clear(); //initialize critical time steps   

  curr_parallel_ts.clear();

  maxParallelFactor.timesteps = 0;   //initialize critical time steps

  for(int i=0; i< NUM_QGATES ; i++){
    maxParallelFactor.parallel_gates[i] = 0;
    structMaxInfo.parallel_gates[i] = 0;
    //tmpParGates.parallel_gates[i] = 0;
  }
  //vectGateInfo.gates.push_back(tmpParGates); //dummy first entry
  funcParallelFactor[F->getName().str()] = vectGateInfo;
  funcMaxParallelFactor[F->getName().str()] = structMaxInfo;
  currParallelFunc.clear();
}

void ModCriticalPath::print_funcQbits(){
  for(map<string, map<int,uint64_t> >::iterator mIter = funcQbits.begin(); mIter!=funcQbits.end(); ++mIter){
    errs() << "Var "<< (*mIter).first << " ---> ";
    for(map<int,uint64_t>::iterator indexIter  = (*mIter).second.begin(); indexIter!=(*mIter).second.end(); ++indexIter){
      errs() << (*indexIter).first << ":"<<(*indexIter).second<< "  ";
    }
    errs() << "\n";
  }
}


void ModCriticalPath::print_qgate(qGate qg){
  errs() << qg.qFunc->getName() << " : ";
  for(int i=0;i<qg.numArgs;i++){
    errs() << qg.args[i].name << qg.args[i].index << ", "  ;
  }
  errs() << "\n";
}


void ModCriticalPath::print_critical_info(string func){
    map<string, allTSParallelism>::iterator fitr = funcParallelFactor.find(func);
    assert(fitr!=funcParallelFactor.end() && "Func not found in funcParFac");
    errs() << "Timesteps = " << (*fitr).second.timesteps << "\n";
    for(unsigned int i = 0; i<(*fitr).second.gates.size(); i++){
        errs() << i << " :";
         for(int k=0;k<NUM_QGATES;k++)
                errs() << " " << (*fitr).second.gates[i].parallel_gates[k];
        errs() << "\n";
        }
        }

void ModCriticalPath::update_critical_info(string currFunc, uint64_t ts, string fname, double angle){
    map<string, allTSParallelism>::iterator fitr = funcParallelFactor.find(fname);
    map<string, allTSParallelism>::iterator citr = funcParallelFactor.find(currFunc);
    assert(fitr!=funcParallelFactor.end() && "parallel func not found in funcParallelFactor");

    assert(citr!=funcParallelFactor.end() && "curr func not found in funcParallelFactor");

    unsigned currTSsize = (*citr).second.gates.size();
    unsigned newTSsize = (*fitr).second.gates.size();
    
    //errs() << "ts = " << ts << " currsize = " << currTSsize << " newsize = "<<newTSsize <<"\n";
    
    if(ts+newTSsize > currTSsize)
        (*citr).second.timesteps = ts+newTSsize;
    
    if(ts + newTSsize <= currTSsize){
        for(unsigned i = 0; i<newTSsize; i++){
	  for(int k=0; k<NUM_QGATES;k++){
                (*citr).second.gates[ts+i].parallel_gates[k] += (*fitr).second.gates[i].parallel_gates[k];            
		/*
		if(angle!=0.0){
		  //TBD: search the vector first
		  vector<double>::iterator ait = (*citr).second.gates[ts+i].angles.find(angle);
		  if(ait==(*citr).second.gates[ts+i].angles.end())
		    (*citr).second.gates[ts+i].angles.push_back(angle);
		    }*/
	  }
	}
    }
    else{
      for(unsigned i = 0; i<currTSsize-ts; i++){
	//errs() << "i = " << i << " ts+i=" << ts+i << "\n";
	for(int k=0; k<NUM_QGATES;k++)
                (*citr).second.gates[ts+i].parallel_gates[k] += (*fitr).second.gates[i].parallel_gates[k];            
            }
            
      for(unsigned i = currTSsize-ts; i<newTSsize; i++){
	//errs() << "i = " << i << "\n";
            ArrParGates tmpParGates;
            for(int k=0; k<NUM_QGATES;k++){
                tmpParGates.parallel_gates[k] = (*fitr).second.gates[i].parallel_gates[k];                              
            }
            (*citr).second.gates.push_back(tmpParGates); 
        }
        
        }

    //print_critical_info(currFunc);
}


uint64_t ModCriticalPath::find_max_funcQbits(){
  uint64_t max_timesteps = 0;
  for(map<string, map<int,uint64_t> >::iterator mIter = funcQbits.begin(); mIter!=funcQbits.end(); ++mIter){
    map<int,uint64_t>::iterator arrIter = (*mIter).second.find(-2);
    if((*arrIter).second > max_timesteps)
      max_timesteps = (*arrIter).second;
  }

  return max_timesteps;

}

void ModCriticalPath::memset_funcQbits(uint64_t val){
  for(map<string, map<int,uint64_t> >::iterator mIter = funcQbits.begin(); mIter!=funcQbits.end(); ++mIter){
    for(map<int,uint64_t>::iterator arrIter = (*mIter).second.begin(); arrIter!=(*mIter).second.end();++arrIter)
      (*arrIter).second = val;
  }
}

void ModCriticalPath::print_scheduled_gate(qGate qg, uint64_t ts){
  string tmpGateName = qg.qFunc->getName();
  if(tmpGateName.find("llvm.")!=string::npos)
    tmpGateName = tmpGateName.substr(5);
  errs() << ts << " : " << tmpGateName;
  for(int i = 0; i<qg.numArgs; i++){
    //if(qg.args[i].index != -1)
      errs() << " " << qg.args[i].name << qg.args[i].index;
  }

  if(tmpGateName == "PrepX" || tmpGateName == "PrepZ"){
    if(qg.angle > 0)
      errs() << ", 1";
    else
      errs() << ", 0";
  }
  else if(tmpGateName == "Rz" || tmpGateName == "Ry" || tmpGateName == "Rx")
    errs() << ", "<<qg.angle;

  errs() << "\n";
}

void ModCriticalPath::print_tableFuncQbits(){
  for(map<Function*, map<unsigned int, map<int, uint64_t> > >::iterator m1 = tableFuncQbits.begin(); m1!=tableFuncQbits.end(); ++m1){
    errs() << "Function " << (*m1).first->getName() << " \n  ";
    for(map<unsigned int, map<int, uint64_t> >::iterator m2 = (*m1).second.begin(); m2!=(*m1).second.end(); ++m2){
      errs() << "\tArg# "<< (*m2).first << " -- ";
      for(map<int, uint64_t>::iterator m3 = (*m2).second.begin(); m3!=(*m2).second.end(); ++m3){
	errs() << " ; " << (*m3).first << " : " << (*m3).second;
      }
      errs() << "\n";
    }
  }
}

void ModCriticalPath::calc_max_parallelism_statistic()
{

    map<string, allTSParallelism>::iterator fitr = funcParallelFactor.find("main");
    assert(fitr!=funcParallelFactor.end() && "Func not found in funcParFac");
    errs() << "Timesteps = " << (*fitr).second.timesteps << "\n";

    uint64_t max_par = 0;
    uint64_t ts_par = 0;

    for(unsigned int i = 0; i<(*fitr).second.gates.size(); i++){

      uint64_t nonZero = 0;

      for(int k=0;k<NUM_QGATES-1;k++){ //skip the 'All' entry
	//count non-zero entries
	if((*fitr).second.gates[i].parallel_gates[k] > 0)
	  nonZero++;
      }

      if(nonZero>max_par){
	max_par = nonZero;
	ts_par = i;
      }
        
    }

    errs() << "Max parallelism in types of gates = " << max_par << " in TS: " << ts_par << "\n";

}


void ModCriticalPath::calc_critical_time(Function* F, qGate qg){
  string fname = qg.qFunc->getName();

  //print_qgate(qg);

  if(isFirstMeas && (fname == "llvm.MeasX" || fname == "llvm.MeasZ")){
    uint64_t maxFQ = find_max_funcQbits();
    memset_funcQbits(maxFQ);
    
    //--print_scheduled_gate(qg,maxFQ+1);


    map<string, map<int,uint64_t> >::iterator mIter = funcQbits.find(qg.args[0].name);
	
    int argIndex = qg.args[0].index;
    
    //update the timestep number for that argument
    map<int,uint64_t>::iterator indexIter = (*mIter).second.find(argIndex);
    (*indexIter).second =  maxFQ + 1;
    
    //update -2 entry for the array, i.e. max ts over all indices
    indexIter = (*mIter).second.find(-2);
    (*indexIter).second = maxFQ + 1;

    //update_critical_info(F->getName().str(), maxFQ, qg.qFunc->getName(), qg.angle);   
    isFirstMeas = false;
  }
  else{ //not a Meas gate
    uint64_t max_ts_of_all_args = 0;
    
    //find last timestep for all arguments of qgate
    for(int i=0;i<qg.numArgs; i++){
      map<string, map<int,uint64_t> >::iterator mIter = funcQbits.find(qg.args[i].name);
      assert(mIter!=funcQbits.end()); //should already have an entry for the name of the qbit
      
      int argIndex = qg.args[i].index;
      
      //find the index of argument in the map<int,int>
      if(argIndex == -1) //operation on entire array
	{
	  //find max for the array
	  map<int,uint64_t>::iterator indexIter = (*mIter).second.find(-2);
	  if((*indexIter).second > max_ts_of_all_args)
	    max_ts_of_all_args = (*indexIter).second;	  
	}
      else
	{
	  map<int,uint64_t>::iterator indexIter = (*mIter).second.find(argIndex);
	  if(indexIter!=(*mIter).second.end()){
	    if((*indexIter).second > max_ts_of_all_args)
	      max_ts_of_all_args = (*indexIter).second;
	  }
	  else{
	    //find the value for entire array
	    map<int,uint64_t>::iterator fullArrayIndexIter = (*mIter).second.find(-1);
	    
	    ((*mIter).second)[argIndex] = (*fullArrayIndexIter).second;
	    if((*fullArrayIndexIter).second > max_ts_of_all_args)
	      max_ts_of_all_args = (*fullArrayIndexIter).second;
	    //((*mIter).second)[argIndex] = 0;
	  }
	}
    }
    
    if(debugModCriticalPath){
      errs() << "Before Scheduling: \n";
      print_funcQbits();
      }

    //errs() << "Max timestep for all args = " << max_ts_of_all_args << "\n";
    
    
    if(fname.find("llvm.")!=string::npos){ //is intrinsic
      
      //schedule gate in max_ts_of_all_args + 1th timestep
      //--print_scheduled_gate(qg,max_ts_of_all_args+1);
      
      //find last timestep for all arguments of qgate
      for(int i=0;i<qg.numArgs; i++){
	map<string, map<int,uint64_t> >::iterator mIter = funcQbits.find(qg.args[i].name);
	
	int argIndex = qg.args[i].index;
	
	if(argIndex == -1){
	  for(map<int,uint64_t>::iterator entryIter = (*mIter).second.begin(); entryIter!=(*mIter).second.end();++entryIter){
	    (*entryIter).second = max_ts_of_all_args + 1;
	  }
	  
	}
	else{
	  //update the timestep number for that argument
	  map<int,uint64_t>::iterator indexIter = (*mIter).second.find(argIndex);
	  (*indexIter).second =  max_ts_of_all_args + 1;
	  
	  //update -2 entry for the array, i.e. max ts over all indices
	  indexIter = (*mIter).second.find(-2);
	  if((*indexIter).second < max_ts_of_all_args + 1)
	    (*indexIter).second = max_ts_of_all_args + 1;
	}  
      }
    } //intrinsic func
    
    else{
      //not an intrinsic function
      //errs() << "Non intrinsic \n";

      map<Function*, uint64_t>::iterator cpit = funcCritPath.find(qg.qFunc);
      assert(cpit != funcCritPath.end() && "func not found in funcCritPath");
      uint64_t lenCritPath = (*cpit).second;
      
      for(int i=0;i<qg.numArgs; i++){
	map<string, map<int,uint64_t> >::iterator mIter = funcQbits.find(qg.args[i].name);
	assert(mIter!=funcQbits.end()); //should already have an entry for the name of the qbit
	
	//differentiate for qbit and qbit*
	
	if(qg.args[i].index == -1){ //qbit*
	  for(map<int,uint64_t>::iterator entryIter = (*mIter).second.begin(); entryIter!=(*mIter).second.end();++entryIter){
	    (*entryIter).second = max_ts_of_all_args + lenCritPath;
	  }
	}
      
	else{ //qbit was passed
	  map<int, uint64_t>::iterator currEntryIt = (*mIter).second.find(qg.args[i].index);
	  if(currEntryIt!=(*mIter).second.end())
	    (*currEntryIt).second = max_ts_of_all_args + lenCritPath;
	  else
	    (*mIter).second[qg.args[i].index] = max_ts_of_all_args + lenCritPath;
	  
	  //update -2 entry for the array, i.e. max ts over all indices
	  currEntryIt = (*mIter).second.find(-2);
	  if((*currEntryIt).second < max_ts_of_all_args + lenCritPath)
	    (*currEntryIt).second = max_ts_of_all_args + lenCritPath;
	  
	}
      }
    } //not intrinsic Gate
  } // not a MeasX gate
  
  if(debugModCriticalPath)
  {   
    errs() << "\nAfter Scheduling: \n";
    print_funcQbits();
    errs() << "\n";
  }
  
}

void ModCriticalPath::analyzeCallInst(Function* F, Instruction* pInst){
  if(CallInst *CI = dyn_cast<CallInst>(pInst))
    {      
      if(debugModCriticalPath)
	errs() << "Call inst: " << CI->getCalledFunction()->getName() << "\n";

      if(CI->getCalledFunction()->getName() == "store_cbit"){	//trace return values
	return;
      }

      vector<qGateArg> allDepQbit;                                  
      
      bool tracked_all_operands = true;

      int myPrepState = -1;
      double myRotationAngle = 0.0;
      
      for(unsigned iop=0;iop<CI->getNumArgOperands();iop++){
	tmpDepQbit.clear();
	
	qGateArg tmpQGateArg;
	btCount=0;
	
	tmpQGateArg.argNum = iop;
	
	
	if(isa<UndefValue>(CI->getArgOperand(iop))){
	  errs() << "WARNING: LLVM IR code has UNDEF values. \n";
	  tmpQGateArg.isUndef = true;	
	  //exit(1);
	}
	
	Type* argType = CI->getArgOperand(iop)->getType();
	if(argType->isPointerTy()){
	  tmpQGateArg.isPtr = true;
	  Type *argElemType = argType->getPointerElementType();
	  if(argElemType->isIntegerTy(16))
	    tmpQGateArg.isQbit = true;
	  //if(argElemType->isIntegerTy(1))
	  //tmpQGateArg.isCbit = true;
	}
	else if(argType->isIntegerTy(16)){
	  tmpQGateArg.isQbit = true;
	  tmpQGateArg.valOrIndex = 0;	 
	}	  	
	//else if(argType->isIntegerTy(1)){
	//tmpQGateArg.isCbit = true;
	//tmpQGateArg.valOrIndex = 0;	 
	//}

	
        //if(tmpQGateArg.isQbit || tmpQGateArg.isCbit){
	if(tmpQGateArg.isQbit){
            tmpDepQbit.push_back(tmpQGateArg);	
            tracked_all_operands &= backtraceOperand(CI->getArgOperand(iop),0);
	}

        if(tmpDepQbit.size()>0){	  
	  allDepQbit.push_back(tmpDepQbit[0]);
	  assert(tmpDepQbit.size() == 1 && "tmpDepQbit SIZE GT 1");
	  tmpDepQbit.clear();
	}
	
      }
      
      if(allDepQbit.size() > 0){
	if(debugModCriticalPath)
	{
	    errs() << "\nCall inst: " << CI->getCalledFunction()->getName();	    
	    errs() << ": Found all arguments: ";       
	    for(unsigned int vb=0; vb<allDepQbit.size(); vb++){
	      if(allDepQbit[vb].argPtr)
		errs() << allDepQbit[vb].argPtr->getName() <<" Index: ";
                                
	      //else
		errs() << allDepQbit[vb].valOrIndex <<" ";
	    }
	    errs()<<"\n";
	    
	}
          
       string fname =  CI->getCalledFunction()->getName();  
       qGate thisGate;
       thisGate.qFunc =  CI->getCalledFunction();

       if(myPrepState!=-1) thisGate.angle = (float)myPrepState;
       if(myRotationAngle!=0.0) thisGate.angle = myRotationAngle;

       for(unsigned int vb=0; vb<allDepQbit.size(); vb++){
            if(allDepQbit[vb].argPtr){
                qGateArg param =  allDepQbit[vb];       
                thisGate.args[thisGate.numArgs].name = param.argPtr->getName();
		if(!param.isPtr)
		  thisGate.args[thisGate.numArgs].index = param.valOrIndex;
                thisGate.numArgs++;
	    }
       }
       
       calc_critical_time(F,thisGate);       

      }    
      allDepQbit.erase(allDepQbit.begin(),allDepQbit.end());
    }
}


void ModCriticalPath::saveTableFuncQbits(Function* F){
  map<unsigned int, map<int, uint64_t> > tmpFuncQbitsMap;

  if(F->getName() == "PARSENODEROOT")
    print_funcQbits();

  for(map<string, map<int, uint64_t> >::iterator mapIt = funcQbits.begin(); mapIt!=funcQbits.end(); ++mapIt){
    map<string, unsigned int>::iterator argIt = funcArgs.find((*mapIt).first);
    if(argIt!=funcArgs.end()){
      unsigned int argNum = (*argIt).second;
      tmpFuncQbitsMap[argNum] = (*mapIt).second;
    }
  }
  tableFuncQbits[F] = tmpFuncQbitsMap;
}

void ModCriticalPath::CountCriticalFunctionResources (Function *F) {
      // Traverse instruction by instruction
  init_critical_path_algo(F);
  
  
  for (inst_iterator I = inst_begin(*F), E = inst_end(*F); I != E; ++I) {
    Instruction *Inst = &*I;                            // Grab pointer to instruction reference
    analyzeAllocInst(F,Inst);          
    analyzeCallInst(F,Inst);	
  }

  //saveTableFuncQbits(F);
  //print_tableFuncQbits();
  //funcParallelFactor[F] = currTS;
}


void ModCriticalPath::init_gates_as_functions(){
  for(int  i =0; i< NUM_QGATES ; i++){
    string gName = gate_name[i];
    string fName = "llvm.";
    fName.append(gName);

    allTSParallelism tmp_info;
    MaxInfo tmp_max_info;

    tmp_info.timesteps = 1;
    tmp_max_info.timesteps = 1;

    ArrParGates tmp_gate_info;
    for(int  k=0; k< NUM_QGATES ; k++){
      tmp_gate_info.parallel_gates[k] = 0;
      tmp_max_info.parallel_gates[k] = 0;
    }
    tmp_gate_info.parallel_gates[i] = 1;
    tmp_gate_info.parallel_gates[_All] = 1;

    tmp_max_info.parallel_gates[i] = 1;
    tmp_max_info.parallel_gates[_All] = 1;
    
    tmp_info.gates.push_back(tmp_gate_info);

    funcParallelFactor[fName] = tmp_info;
    funcMaxParallelFactor[fName] = tmp_max_info;
    
  }
}


bool ModCriticalPath::runOnModule (Module &M) {
  init_gate_names();
  init_gates_as_functions();
  
  // iterate over all functions, and over all instructions in those functions
  CallGraphNode* rootNode = getAnalysis<CallGraph>().getRoot();
  
  //Post-order
  for (scc_iterator<CallGraphNode*> sccIb = scc_begin(rootNode), E = scc_end(rootNode); sccIb != E; ++sccIb) {
    const std::vector<CallGraphNode*> &nextSCC = *sccIb;
    for (std::vector<CallGraphNode*>::const_iterator nsccI = nextSCC.begin(), E = nextSCC.end(); nsccI != E; ++nsccI) {
      Function *F = (*nsccI)->getFunction();	  
            
      if(F && !F->isDeclaration()){
	errs() << "\nFunction: " << F->getName() << "\n";      

	if(F->getName()=="measure_z")
	  debugModCriticalPath = true;

	funcQbits.clear();
	//errs() << "Pt2 \n";
	funcArgs.clear();

	//errs() << "pt3 \n";
	getFunctionArguments(F);

	// count the critical resources for this function
	CountCriticalFunctionResources(F);

	//if(F->getName() == "main")
	//errs() << F->getName() << ": " << "Critical Path Length : " << find_max_funcQbits() << "\n";
	errs() << F->getName() << " " << find_max_funcQbits() << "\n";

	funcCritPath[F] = find_max_funcQbits();
	errs() << "Pt1 \n";
      }
      else{
	    if(debugModCriticalPath)
	      errs() << "WARNING: Ignoring external node or dummy function.\n";
	  }
    }
  }
  //print_tableFuncQbits();
  //print_critical_info("main");

  //calc_max_parallelism_statistic();
  
  return false;
} // End runOnModule
