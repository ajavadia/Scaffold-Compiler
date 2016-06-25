//===----------------------- RollupLoops.cpp --------------------------===///
// This file implements the Scaffold pass of traversing basic blocks, finding
// purely quantum loops, and transforming them so they get called only once.


//===----------------------------------------------------------------------===//

#define DEBUG_TYPE "RollupLoops"
#include "llvm/Module.h"
#include "llvm/Pass.h"
#include "llvm/Function.h"
#include "llvm/BasicBlock.h"
#include "llvm/Instruction.h"
#include "llvm/Instructions.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/ADT/Statistic.h"
#include "llvm/Support/InstIterator.h"
#include "llvm/PassAnalysisSupport.h"
#include "llvm/Analysis/LoopInfo.h"
#include "llvm/Analysis/ScalarEvolution.h"
#include "llvm/Analysis/ScalarEvolutionExpressions.h"
#include "llvm/Transforms/Utils/Local.h"
#include "llvm/Intrinsics.h"
#include <sstream>
#include <climits>

using namespace llvm;
using namespace std;


bool debugRollupLoops = false;


STATISTIC(NumLoopsRolled, "Number of loops rolled up");
STATISTIC(NumTotalLoops, "Number of total loops");

// An anonymous namespace for the pass. Things declared inside it are
// only visible to the current file.
namespace {

    // Derived from FunctionPass to count qbits in functions
  struct RollupLoops : public ModulePass {
    static char ID; // Pass identification
    
    map<BasicBlock*, bool> blockCollapse; //roll this loop
    map<BasicBlock*, string> blockRepCond; //repeat condition
    map<BasicBlock*, BasicBlock*> blockBodyInc; //map for.body and for.inc
    vector<BasicBlock*> blockForAll; //list of forall basicblocks
    map<BasicBlock*, int> bbTripCount; //non zero trip count of loops
    
    RollupLoops() : ModulePass(ID) {}
    //AnalysisUsage AU;
    
    virtual void getAnalysisUsage(AnalysisUsage &AU) const {
      AU.addRequired<LoopInfo>();
      AU.addPreserved<LoopInfo>();
      AU.addRequired<ScalarEvolution>();
      AU.addPreserved<ScalarEvolution>();            
    }
    

    /*void print_blockForAll(){
      errs() << "Printing blockForAll: \n";
      for(vector<BasicBlock*>::iterator vit = blockForAll.begin(); vit!=blockForAll.end();++vit)
	errs() << (*vit) << ":" << (*vit)->getName() << " ";
      errs() << "\n";
     
      }*/
    
    bool checkIfQuantumType(Type* t) {
      bool isQtmTy = true;
      
      if(t->isIntegerTy(16) || t->isIntegerTy(1))
	isQtmTy &= true;
      else if(ArrayType *arrayType = dyn_cast<ArrayType>(t))
	{
	  Type *elementType = arrayType->getElementType();
	  if(elementType->isIntegerTy(16) || elementType->isIntegerTy(1))
	    isQtmTy &= true;
	  else isQtmTy &= false;
	}	      
      else if(t->isPointerTy()){
	Type* elementType = t->getPointerElementType();
	if(elementType->isIntegerTy(16) || elementType->isIntegerTy(1))
	  isQtmTy &= true;
      }
      else{
	isQtmTy = false;
	
      }
      
      return isQtmTy;
    }
    
    
    bool checkIfQuantum(BasicBlock* BB) {

      //Does not allow arithmetic/compare instructions
      //Does not allow function calls to have non-qbit/non-cbit arguments
      //all gepi instructions muct depend directly on indvar

      bool isQtm = true;

      string indVarStr = "";

      /*if(BB->size() <= 2){ //phi inst and branch inst 
	errs() << "Num of insts <= 2 \n";
	return false;
	}*/

      //errs() << "BB name: " << BB->getName() << "\n";

      bool hasOtherThanPhiAndBrInst = false;      
      //bool bodyHasIndVarComputation = false;
      
      for(BasicBlock::iterator bbi = BB->begin(); bbi != BB->end(); ++bbi)
	{
	  if(debugRollupLoops)
	    errs() << *bbi << "\n";
	  
	  if(GetElementPtrInst *GEPI = dyn_cast<GetElementPtrInst>(&*bbi)){
	    if(debugRollupLoops)
	      errs() << "Is a GEPI \n";

	    hasOtherThanPhiAndBrInst = true;
	    
	    //check type of GEPI inst first
	    Type* gepiType = GEPI->getType();
	    
	    if(!checkIfQuantumType(gepiType)){
	      isQtm = false;
	      if(debugRollupLoops)
		errs() << "Unrecognized Type of GEPI Inst. \n";
	      break;
	    }
	    
	    //check indices of GEPI inst next: should be constant or indvars
	    if(isQtm && GEPI->hasIndices())
	      {
		unsigned numOps = GEPI->getNumOperands();
		for(unsigned opIter=1; opIter < numOps; opIter++)
		  {
		    Value* opInst = GEPI->getOperand(opIter);

		    if(opInst->getName().str()==indVarStr){
		      //bodyHasIndVarComputation = true;
		      //errs() << "IndVar is present in block \n";
		      blockForAll.push_back(BB);
		      //print_blockForAll();
		    }      

		    if(isa<ConstantInt>(opInst)){
		      isQtm &= true;
		    }
		    //else if(opInst->getName().find("indvars.iv")!=string::npos)
		    else if(opInst->getName().str()==indVarStr)
		      isQtm &= true;
		    else{
		      isQtm = false;
		      if(debugRollupLoops)
			errs() << "GEPI indices do not satisfy qtm block criteria \n";
		      return isQtm;
		      //break;
		    }
		    
		  }
	      }
	  }
	  
	  else if(isa<LoadInst>(*bbi)){
	    if(debugRollupLoops)
	      errs() << "Is a Load \n";
	  }
	  

	  else if(isa<TruncInst>(*bbi)){
	    if(debugRollupLoops)
	      errs() << "Is Trunc \n";
	  }

	  else if(isa<PHINode>(*bbi)){
	    if(debugRollupLoops)
	      errs() << "Is a PHI Node \n";
	      //errs() << *bbi << "\n";
	    indVarStr = (*bbi).getName();
	      //errs() << "Phi Node Var = " << (*bbi).getName() << "\n";
	  }
	  
	  else if(CallInst *CI = dyn_cast<CallInst>(&*bbi))
	    {
	      hasOtherThanPhiAndBrInst = true;
	      /*
	      if((CI->getCalledFunction()->getName() == "llvm.PrepZ") || (CI->getCalledFunction()->getName() == "llvm.PrepX")){
		isQtm &= true;
		//errs() << "Found PrepZ or PrepX\n";		
	      }
	      else{
	      */

		//check if all operands of call inst are i16 or i16*
		//if not i16 or i16*, they should be constant integers i64 or i32
		for(unsigned iop=0;iop<CI->getNumArgOperands();iop++){
		  Type* argType = CI->getArgOperand(iop)->getType();
		  
		  //errs() << "Called Func = " << CI->getCalledFunction()->getName() << "\n";
		  //errs() << "Call inst = " << *CI << "\n";

		  if(argType->isIntegerTy(32) || argType->isIntegerTy(64)){
		    //must be constant int
		    if(!isa<ConstantInt>(CI->getArgOperand(iop))){
		      if(debugRollupLoops)
		      errs() << "Not a constant int \n";
		      return false;
		    }
		  }

		  else if(!checkIfQuantumType(argType)){
		    if(debugRollupLoops)
		      errs() << "Operand Type of Call inst is non-quantum \n";
		    //isQtm = false;
		    return false;
		    
		  }
		  //else errs() << "Is qtm \n";
		}
		//}
	    }
	  
	  else if(isa<BranchInst>(*bbi)){
	    if(debugRollupLoops)
	      errs() << "Is a Branch Inst \n";
	  }
	  
	  else{
	    if(debugRollupLoops)
	      errs() << "Unrecognized in quantum module. \n";
	    //isQtm = false;
	    return false;
	    //break;
	  }
	}

      if(hasOtherThanPhiAndBrInst)
	return isQtm;
      else{
	//errs() << "Only contains Phi and Br \n";
	return false;

      }
    }
    
    
    
    void getIncAndCmpConditions(Loop *L, BasicBlock *BB, BasicBlock* Header){
      if(debugRollupLoops)
      errs() << "Basic Block: " << BB->getName() << "\n";

      map<BasicBlock*, bool>::iterator mapIter = blockCollapse.find(Header);
      if(mapIter!=blockCollapse.end())
	if((*mapIter).second == true)
	  {

	    if(debugRollupLoops)
	      errs() << "Get Inc and Cmp conditions \n";
	    
	    stringstream loopRepStr;
	    stringstream incVal;
	    stringstream cmpVal;
	    stringstream initVal;

	    //PHINode* PN;
	    //unsigned IncomingEdge;
	    //unsigned BackEdge;


	    //get initVal str
	    for(BasicBlock::iterator hIter = Header->begin(); hIter!=Header->end(); ++hIter)
	      {
		//find phinode
		if(PHINode *PN = dyn_cast<PHINode>(&*hIter)){ //assuming we will get one PHI node at the beginning of the block
		  //errs() << "Found PHI Node \n";
		  //errs() << *PN << "\n";

		  unsigned IncomingEdge = L->contains(PN->getIncomingBlock(0));	       
		  unsigned BackEdge     = IncomingEdge^1;
		  if(ConstantInt *CInt = dyn_cast<ConstantInt>(PN->getIncomingValue(IncomingEdge))){
		    int val = CInt->getZExtValue();
		    if(debugRollupLoops)
		      errs() << "Incoming Value = "<< val << "\n";
		    initVal << "i="<<val;		
		    
		  }

		  if(BinaryOperator *Incr = dyn_cast<BinaryOperator>(PN->getIncomingValue(BackEdge))){
		    //Instruction* bInst = &*bbi;
		    unsigned instOp = Incr->getOpcode();
		    
		    //for(BasicBlock::iterator bbi = BB->begin(); bbi != BB->end(); ++bbi)
		    //{
		    if(debugRollupLoops)
		    errs() << *Incr << "\n";
		
		    //Instruction* bInst = &*bbi;
		    //unsigned instOp = bInst->getOpcode();
		    if(instOp == Instruction::Add || instOp == Instruction::Sub){
		      
		      unsigned numOps = Incr->getNumOperands();
		      for(unsigned opIter=0; opIter < numOps; opIter++)
			{
			  Value* opInst = Incr->getOperand(opIter);
			  
			  if(ConstantInt *CI = dyn_cast<ConstantInt>(opInst)){
			    int val = CI->getZExtValue();
			    if(debugRollupLoops)
			      errs() << "Arith value = "<< val << "\n";
			    
			    if(instOp == Instruction::Add){
			      if(val == 1) incVal << "i++";
			      else if(val == -1) incVal << "i--";
			      else if(val>0) incVal << "i+=" << val;
			      else if(val<0) incVal << "i-=" << val;		  
			    }
			    
			    if(instOp == Instruction::Sub){
			      if(val == 1) incVal << "i--";
			      else if(val == -1) incVal << "i++";
			      else if(val>0) incVal << "i-=" << val;
			      else if(val<0) incVal << "i+=" << val;		  
			    }
			    
			  }
			  //else if(Incr->getName().find("indvars.iv")==string::npos)
			    //errs() << "WARNING: Not an inc instruction. \n";
			}
		      
		    } //end of inst::Add or Inst::Sub		    
		  } //end of Incr
		}
	      }
		    		    
	    for(BasicBlock::iterator bbi = BB->begin(); bbi != BB->end(); ++bbi)
	      {
		if(debugRollupLoops)
		  errs() << "Inst : " << *bbi << "\n";
		
		Instruction* bInst = &*bbi;
		
		if(ICmpInst *IC = dyn_cast<ICmpInst>(bInst))
		  {
		    if(debugRollupLoops)
		      errs() << "Compare Inst = "<<*IC<<"\n";
		    int val = INT_MAX;
		    
		    unsigned numOps = IC->getNumOperands();
		    for(unsigned opIter=0; opIter < numOps; opIter++)
		      {
			Value* opInst = IC->getOperand(opIter);
			
			if(ConstantInt *CI = dyn_cast<ConstantInt>(opInst)){
			  val = CI->getZExtValue();
			  if(debugRollupLoops)
			    errs() << "Arith value = "<< val << "\n";
			  
			}
			//else if(opInst->getName().find("indvars.iv")==string::npos)
			//errs() << "WARNING: Not an inc instruction. \n";
		      }
		    
		    if(val != INT_MAX)
		      {
			unsigned icmp_cond = IC->getUnsignedPredicate();
			switch(icmp_cond){
			case CmpInst::ICMP_EQ:
			  cmpVal << "i=" << val;
			  //errs() << "equal \n";
			  break;
			case CmpInst::ICMP_NE:
			  cmpVal << "i!=" << val;
			  //errs() << "not equal \n";
			  break;
			case CmpInst::ICMP_UGT:
			case CmpInst::ICMP_SGT:
			  cmpVal << "i>" << val;
			  //errs() << "greater \n";
			  break;
			case CmpInst::ICMP_UGE:
			case CmpInst::ICMP_SGE:
			  cmpVal << "i>=" << val;
			  //errs() << "greater or equal \n";
			  break;
			case CmpInst::ICMP_ULT:
			case CmpInst::ICMP_SLT:
			  cmpVal << "i<" << val;
			  //errs() << "less \n";
			  break;
			case CmpInst::ICMP_ULE:
			case CmpInst::ICMP_SLE:
			  cmpVal << "i<=" << val;
			  //errs() << "less equal \n";
			  break;		  
			} //end of switch
		      }	       	      
		    
		  } //end of Icmp Inst	
	      } //end of BB iterator       		      
	    	    
	    if(initVal.str()!="" && cmpVal.str()!="" && incVal.str()!=""){
	      //errs() << "Searching blockForAll for " << Header << ":" << Header->getName() << "\n";
	      vector<BasicBlock*>::iterator isForAll = find(blockForAll.begin(), blockForAll.end(),Header);
	      if(isForAll!=blockForAll.end())		
		loopRepStr << "forall "<< initVal.str() << "; " << cmpVal.str() << "; " << incVal.str();
	     
	      else
		loopRepStr << "repeat "<< initVal.str() << "; " << cmpVal.str() << "; " << incVal.str();
		
	    }
	    if(debugRollupLoops)
	      errs() << "Rep String1 = " << loopRepStr.str() << "\n";
	    
	    blockRepCond[Header] = loopRepStr.str();
	    blockBodyInc[Header] = BB;
	  }
      
    }

    bool checkIfQuantumHeaderLatch(Loop *L, BasicBlock *BB){
      //One arithmetic/compare instruction on indvar
      //Does not allow function calls to have non-qbit/non-cbit arguments
      //all gepi instructions muct depend directly on indvar

      //distinguishes between repeat loops and forall loops
      //repeat loops do not have indvar in its body
      
      //errs() << "In checkIfQuantumHeaderLatch \n";

      bool isQtm = true;

      bool bodyHasIndVarComputation = false;

      string indVarStr = "";
      stringstream loopRepStr;
      stringstream incVal;
      stringstream cmpVal;
      stringstream initVal;
      
      bool hasOtherThanPhiAndBrInst = false;      

      Instruction* incrInst = NULL;
      
      for(BasicBlock::iterator bbi = BB->begin(); bbi != BB->end(); ++bbi)
	{
	  if(debugRollupLoops)
	    errs() << *bbi << "\n";
	  
	  if(GetElementPtrInst *GEPI = dyn_cast<GetElementPtrInst>(&*bbi)){
	    if(debugRollupLoops)
	      errs() << "Is a GEPI \n";

	    hasOtherThanPhiAndBrInst = true;
	    
	    //check type of GEPI inst first
	    Type* gepiType = GEPI->getType();
	    
	    if(!checkIfQuantumType(gepiType)){
	      isQtm = false;
	      if(debugRollupLoops)
		errs() << "Unrecognized Type of GEPI Inst. \n";
	      break;
	    }
	    
	    //check indices of GEPI inst next: should be constant or indvars
	    if(isQtm && GEPI->hasIndices())
	      {
		unsigned numOps = GEPI->getNumOperands();
		for(unsigned opIter=1; opIter < numOps; opIter++)
		  {
		    Value* opInst = GEPI->getOperand(opIter);

		    if(opInst->getName().str()==indVarStr){
		      //errs() << "2: Body has IndVar \n";
		      bodyHasIndVarComputation = true;
		    }    

		    if(isa<ConstantInt>(opInst)){
		      isQtm &= true;
		    }
		    //else if(opInst->getName().find("indvars.iv")!=string::npos)
		    else if(opInst->getName().str()==indVarStr)
		      isQtm &= true;
		    else{
		      //isQtm = false;
		      if(debugRollupLoops)
			errs() << "GEPI indices do not satisfy qtm block criteria \n";
		      return false;
		      //break;
		    }
		    
		  }
	      }
	  }
	  
	  else if(isa<LoadInst>(*bbi)){
	    if(debugRollupLoops)
	      errs() << "Is a Load \n";
	  }
	  

	  else if(isa<TruncInst>(*bbi)){
	    if(debugRollupLoops)
	      errs() << "Is Trunc \n";
	  }


	  else if(PHINode *PN = dyn_cast<PHINode>(&*bbi)){ //assuming we will get one PHI node at the beginning of the block
	    
	    if(debugRollupLoops)
	      errs() << "Is a PHI Node \n";
	    //errs() << *bbi << "\n";
	    indVarStr = (*bbi).getName();
	    //errs() << "Phi Node Var = " << (*bbi).getName() << "\n";
	    
	    
	    //errs() << "Found PHI Node \n";
	    //errs() << *PN << "\n";
	    
	    unsigned IncomingEdge = L->contains(PN->getIncomingBlock(0));	       
	    unsigned BackEdge     = IncomingEdge^1;
	    if(ConstantInt *CInt = dyn_cast<ConstantInt>(PN->getIncomingValue(IncomingEdge))){
	      int val = CInt->getZExtValue();
	      if(debugRollupLoops)
		errs() << "Incoming Value = "<< val << "\n";
	      initVal << "i="<<val;		
	      
	    }
	    
	    if(BinaryOperator *Incr = dyn_cast<BinaryOperator>(PN->getIncomingValue(BackEdge))){
	      //Instruction* bInst = &*bbi;
	      unsigned instOp = Incr->getOpcode();
	      incrInst = Incr; //save incr inst
	      
	      if(debugRollupLoops)
		errs() << *Incr << "\n";
	      
	      if(instOp == Instruction::Add || instOp == Instruction::Sub){
		
		unsigned numOps = Incr->getNumOperands();
		for(unsigned opIter=0; opIter < numOps; opIter++)
		  {
		    Value* opInst = Incr->getOperand(opIter);
		    
		    if(ConstantInt *CI = dyn_cast<ConstantInt>(opInst)){
		      int val = CI->getZExtValue();
		      if(debugRollupLoops)
			errs() << "Arith value = "<< val << "\n";
		      
		      if(instOp == Instruction::Add){
			if(val == 1) incVal << "i++";
			else if(val == -1) incVal << "i--";
			else if(val>0) incVal << "i+=" << val;
			else if(val<0) incVal << "i-=" << val;		  
		      }
		      
		      if(instOp == Instruction::Sub){
			if(val == 1) incVal << "i--";
			else if(val == -1) incVal << "i++";
			else if(val>0) incVal << "i-=" << val;
			else if(val<0) incVal << "i+=" << val;		  
		      }
		      
		    }
		    else if(Incr->getName().find("indvars.iv")==string::npos)
		      errs() << "WARNING: Not an inc instruction. \n";
		  }
		
	      } //end of inst::Add or Inst::Sub		    
	    } //end of Incr
	  }	    	  
	  
	  else if(CallInst *CI = dyn_cast<CallInst>(&*bbi))
	    {
	      hasOtherThanPhiAndBrInst = true;

	      /*if((CI->getCalledFunction()->getName() == "llvm.PrepZ") || (CI->getCalledFunction()->getName() == "llvm.PrepX")){
		isQtm &= true;
		//errs() << "Found PrepZ or PrepX\n";		
	      }
	      else{
	      */
		//check if all operands of call inst are i16 or i16*
		//if not i16 or i16*, they should be constant integers i64 or i32
		for(unsigned iop=0;iop<CI->getNumArgOperands();iop++){
		  //if(debugRollupLoops)
		  //errs() << "Processing arg num = " << iop << "\n";

		  Type* argType = CI->getArgOperand(iop)->getType();
		  
		  //errs() << "Called Func = " << CI->getCalledFunction()->getName() << "\n";
		  //errs() << "Call inst = " << *CI << "\n";
		  if(argType->isIntegerTy(32) || argType->isIntegerTy(64)){
		    //errs() << "Arg type integer \n";
		    //must be constant int
		    if(!isa<ConstantInt>(CI->getArgOperand(iop))){
		      if(debugRollupLoops)
			errs() << "Not a constant int \n";
		      return false;
		    }
		  }

		  else if(!checkIfQuantumType(argType)){
		    if(debugRollupLoops)
		      errs() << "Operand Type of Call inst is non-quantum \n";
		    //isQtm = false;
		    return false;
		    
		  }
		  //else errs() << "Call inst is qtm \n";
		}
		//}
	    }
	  

	  
	  else if(isa<BranchInst>(*bbi)){
	    if(debugRollupLoops)
	      errs() << "Is a Branch Inst \n";
	  }
	  

	  
	  else if(ICmpInst *IC = dyn_cast<ICmpInst>(&*bbi))
	    {
	      if(debugRollupLoops)
		errs() << "Compare Inst = "<<*IC<<"\n";
	      int val = INT_MAX;
	      
	      unsigned numOps = IC->getNumOperands();
	      for(unsigned opIter=0; opIter < numOps; opIter++)
		{
		  Value* opInst = IC->getOperand(opIter);
		  
		  if(ConstantInt *CI = dyn_cast<ConstantInt>(opInst)){
		    val = CI->getZExtValue();
		    if(debugRollupLoops)
		      errs() << "Arith value = "<< val << "\n";
		    
		  }
		  //else if(opInst->getName().find("indvars.iv")==string::npos)
		  //errs() << "WARNING: Not an inc instruction. \n";
		}
	      
	      if(val != INT_MAX)
		{
		  unsigned icmp_cond = IC->getUnsignedPredicate();
		  switch(icmp_cond){
		  case CmpInst::ICMP_EQ:
		    cmpVal << "i=" << val;
		    //errs() << "equal \n";
		    break;
		  case CmpInst::ICMP_NE:
		    cmpVal << "i!=" << val;
		    //errs() << "not equal \n";
		    break;
		  case CmpInst::ICMP_UGT:
		  case CmpInst::ICMP_SGT:
		    cmpVal << "i>" << val;
		    //errs() << "greater \n";
		    break;
		  case CmpInst::ICMP_UGE:
		  case CmpInst::ICMP_SGE:
		    cmpVal << "i>=" << val;
		    //errs() << "greater or equal \n";
		    break;
		  case CmpInst::ICMP_ULT:
		  case CmpInst::ICMP_SLT:
		    cmpVal << "i<" << val;
		    //errs() << "less \n";
		    break;
		  case CmpInst::ICMP_ULE:
		  case CmpInst::ICMP_SLE:
		    cmpVal << "i<=" << val;
		    //errs() << "less equal \n";
		    break;		  
		  } //end of switch
		}	       	      
	      
	    } //end of Icmp Inst	

	  else if(BinaryOperator *BIncr = dyn_cast<BinaryOperator>(&*bbi)){
	    if(BIncr != incrInst){
	      if(debugRollupLoops)
		errs() << "unrecognized incr instr \n";
	      //debugRollupLoops = false;
	      return false;

	    }
	  }
	  else{
	    if(debugRollupLoops)
	      errs() << "Unrecognized in quantum module. \n";
	    //isQtm = false;
	    //debugRollupLoops = false;
	    return false;
	    //break;
	  }
	}      
      
      if(hasOtherThanPhiAndBrInst && isQtm){
	if(initVal.str()!="" && cmpVal.str()!="" && incVal.str()!=""){
	  if(bodyHasIndVarComputation)
	    loopRepStr << "forall "<< initVal.str() << "; " << cmpVal.str() << "; " << incVal.str();
	  else
	    loopRepStr << "repeat "<< initVal.str() << "; " << cmpVal.str() << "; " << incVal.str();
	}

	if(debugRollupLoops)
	  errs() << "Rep String2 = " << loopRepStr.str() << "\n";
	
	blockRepCond[BB] = loopRepStr.str();
	blockBodyInc[BB] = BB;
	//debugRollupLoops = false;
	return isQtm;
      }
      else{
	//errs() << "Only contains Phi and Br \n";
	//debugRollupLoops = false;
	return false;

      }
    }    
    
    void getLoopInfo(Function& F) {
      LoopInfo *LI = &getAnalysis<LoopInfo> ( F );
      ScalarEvolution *SE = &getAnalysis<ScalarEvolution>( F );

      for (Function::iterator BB = F.begin(), E = F.end(); BB != E; ++BB)
	{           
	  
	  string BBname = BB->getName();
	  //errs() << "Basic Block Name = "<<BBname << "\n";	    
	  
	  
	  Loop* L = LI->getLoopFor(&*BB);
	  if (L != NULL)
	    {    
	      
	      if(BBname.find("for.body")!=string::npos){
		++NumTotalLoops;
		if(debugRollupLoops)
		  errs() << "Counting this: "<<BBname<<"\n";

		BasicBlock* Latch = L->getLoopLatch();
                if(Latch){

		  int tripCount = SE->getSmallConstantTripCount(L, Latch);
		  if(debugRollupLoops)
		    errs() << "Trip Count of " << BBname << " is " << tripCount <<"\n";

		  if(tripCount != 0)
		    bbTripCount[BB] = tripCount;

		  if(Latch!=&*BB){
		    bool isCollapsable = checkIfQuantum(&*BB);
		    if(debugRollupLoops)
		      errs() << "isCollapsable = " << isCollapsable << "\n";

		    if(isCollapsable){
		      blockCollapse[&*BB] = isCollapsable;		    
		      
		      //errs() << "Latch is: " << Latch->getName() << "\n";
		      
		      getIncAndCmpConditions(L, Latch, &*BB);
		    }
		  }
		  else{ //Latch and BB are the same basic blocks
		    bool isCollapsable = checkIfQuantumHeaderLatch(L,&*BB);
		    if(debugRollupLoops)
		      errs() << "isCollapsable = " << isCollapsable << "\n";
		    
		    blockCollapse[&*BB] = isCollapsable;
		    

		    //getIncAndCmpConditions(L, Latch, &*BB);


		    }

		}

	      }
	      /*
	      else if(BBname.find("for.inc") != string::npos){	      
		BasicBlock *Header = L->getHeader();                        
		errs() << "Parent: F[" << Header->getParent()->getName()
		       << "] Loop %" << Header->getName() << "\n";
		
		getIncAndCmpConditions(L, &*BB, Header);
		
		}*/
	    }
	  
	  //errs() << "\n";
	}
    }
    
   
    void addDummyStart(Module &M, BasicBlock* BB, string repStr){
      
      BasicBlock::iterator BBiter = BB->getFirstNonPHI();

      while(isa<AllocaInst>(BBiter))
	++BBiter;

      int tripCount = 0;

      map<BasicBlock*,int>::iterator tcIter = bbTripCount.find(BB);
      if(tcIter!=bbTripCount.end())
	tripCount = (*tcIter).second;

      //errs() << "Adding TripCount of " << BB->getName() << " is "<<tripCount<<"\n";

      string dummyFnName = "qasmRepLoopStart";
      dummyFnName.append(repStr);
      if(debugRollupLoops)
	errs() << "Dummy Fn Name = " << dummyFnName << "\n";

      /*
      //addDummyFunc
      Function* dummyStartLoop; //dummyFunc to mark start of loop
      dummyStartLoop = cast<Function>(M.getOrInsertFunction(dummyFnName, Type::getVoidTy(BB->getContext()), (Type*)0));
            
      CallInst::Create(dummyStartLoop, "",(Instruction*)&(*BBiter));      
      */

      
	//add dummy function with trip count as argument
	//will need to change function cloning pass to ignore the argument
      Function* dummyStartLoop; //dummyFunc to mark start of loop
      dummyStartLoop = cast<Function>(M.getOrInsertFunction(dummyFnName, Type::getVoidTy(BB->getContext()), Type::getInt32Ty(BB->getContext()), (Type*)0));
            
      Value* intArg = ConstantInt::get(Type::getInt32Ty(BB->getContext()),tripCount);

      CallInst::Create(dummyStartLoop, intArg, "",(Instruction*)&(*BBiter));      
      
      
    }


    void addDummyEnd(Module &M, BasicBlock* BB, string repStr){

     TerminatorInst *BBTerm = BB->getTerminator();

     //while(isa<AllocaInst>(BBiter))
     //++BBiter;

      string dummyFnName = "qasmRepLoopEnd";
      dummyFnName.append(repStr);
      if(debugRollupLoops)
	errs() << "Dummy Fn Name = " << dummyFnName << "\n";

      //addDummyFunc
      Function* dummyEndLoop; //dummyFunc to mark start of loop
      dummyEndLoop = cast<Function>(M.getOrInsertFunction(dummyFnName, Type::getVoidTy(BB->getContext()), (Type*)0));
            
      CallInst::Create(dummyEndLoop, "",(Instruction*)BBTerm);      




      return;

    }
    void modifyIndVars(Loop* L, BasicBlock* BB){
      
      string indVarStr = "";

      for(BasicBlock::iterator bbi = BB->begin(); bbi!=BB->end(); ++bbi)
	{
	  if(debugRollupLoops)
	    errs()<<"\t" << *bbi<<"\n";

	  if(PHINode *PN = dyn_cast<PHINode>(&*bbi)){
            if(debugRollupLoops)
	      errs() << "Is a PHI Node: " << *bbi << "\n";                                               
            indVarStr = (*bbi).getName();
	    //errs() << "Phi Node Var = " << (*bbi).getName() << "\n";    	 

	    //rewrite Init Value
	    unsigned IncomingEdge = L->contains(PN->getIncomingBlock(0));
	    //unsigned BackEdge     = IncomingEdge^1;

	    Value* currVal = PN->getIncomingValue(IncomingEdge);
	    	    
	    PN->setIncomingValue(IncomingEdge,Constant::getNullValue(currVal->getType()));	    

	  }
	  else if(GetElementPtrInst *GEPI = dyn_cast<GetElementPtrInst>(&*bbi)){
	    //check indices of GEPI inst next: should be constant or indvars          
            if(GEPI->hasIndices())
              {
		unsigned numOps = GEPI->getNumOperands();
                for(unsigned opIter=1; opIter < numOps; opIter++)
                  {
                    Value* opInst = GEPI->getOperand(opIter);               
                    if(opInst->getName().str()==indVarStr){
		      //GEPI->setOperand(opIter, Constant::getNullValue(Type::getInt16Ty(BB->getContext())));
		      GEPI->setOperand(opIter, ConstantInt::get(Type::getInt16Ty(BB->getContext()),-2));
		    }           
                  }
              }
	  }
	}           
    }

    void transformInc(Loop* L, BasicBlock* BB){
      //errs() << "In Transform Inc\n";
            

      Instruction* IncrInst;
      
      for(BasicBlock::iterator bbi = BB->begin(); bbi != BB->end(); ++bbi)
	{
	  if(debugRollupLoops)
	    errs() << *bbi << "\n";
	    
	  Instruction* bInst = &*bbi;	  

	  if(BinaryOperator *Incr = dyn_cast<BinaryOperator>(bInst)){

	  unsigned instOp = Incr->getOpcode();
	  if(instOp == Instruction::Add) // || instOp == Instruction::Sub){
	    {

	      IncrInst = Incr;

	      //errs() << "Inc inst name = "<<bInst->getName() << "\n";
	      //assert(Incr->getName().find("indvars.iv")!=string::npos && "Not an inc instr in for.inc");

	      unsigned numOps = Incr->getNumOperands();
	      for(unsigned opIter=0; opIter < numOps; opIter++)
		{
		  Value* opInst = Incr->getOperand(opIter);

		  if(isa<ConstantInt>(opInst)){
		    Incr->setOperand(opIter, ConstantInt::get(opInst->getType(),1));
		    break;		    
		  }

		}
	    } //end of Inst::Add
	  else if(instOp == Instruction::Sub){ //should not contain a sub inst
	    unsigned numOps = Incr->getNumOperands();
	    for(unsigned opIter=0; opIter < numOps; opIter++)
	      {
		Value* opInst = Incr->getOperand(opIter);

		if(isa<ConstantInt>(opInst)){
		  Incr->setOperand(opIter, ConstantInt::get(opInst->getType(),-1));
		  break;
		}

	      }
	  } //end of Inst::Sub     
	  } //end of if Binary Operator
	  
      //Rewrite Exit Value and Condition
	  if(ICmpInst *IC = dyn_cast<ICmpInst>(bInst))
	    {
	      
	      CmpInst::Predicate NewPred = CmpInst::ICMP_SLT;	     
	      BranchInst *Br = cast<BranchInst>(IC->use_back());
	      assert(Br->isConditional() && "Did not find a branch");

	      /*
	      if(ConstantInt *Val1 = dyn_cast<ConstantInt>(IC->getOperand(1))){
		errs() << "Opd1: Val: "<< Val1->getZExtValue() << "\n";		
	      }
	      else{

		errs() << "Opd1: Name: "<< IC->getOperand(1)->getName() << "\n";		

	      }

	      if(ConstantInt *Val0 = dyn_cast<ConstantInt>(IC->getOperand(0))){
		errs() << "Opd0: Val: "<< Val0->getZExtValue() << "\n";		
	      }
	      else{

		errs() << "Opd0: Name: "<< IC->getOperand(0)->getName() << "\n";		

		}	*/      

	      Value* ICop = IC->getOperand(0);

	      ICmpInst *NewCmp = new ICmpInst(Br,NewPred,ICop,ConstantInt::get(ICop->getType(),1),IC->getName());

	      NewCmp->takeName(IC);
	      IC->replaceAllUsesWith(NewCmp);
	      RecursivelyDeleteTriviallyDeadInstructions(IC);
	      break;
	    }

	} //end of BB iterator      

      return;
    }


    void getLoopInfoAndModify(Function &F,BasicBlock* thisBB,BasicBlock* thisInc){


    //errs() << "BB = " << thisBB->getName() << " Func=" << F->getName() <<"\n";                   
              LoopInfo *LI = &getAnalysis<LoopInfo> ( F );                                                                Loop* L1 = LI->getLoopFor(&*thisBB);                                                           
              //assert(L1 != NULL); //must be a loop by this point                                           
              modifyIndVars(L1, thisBB);                                                                                                                                                                                
              Loop* L2 = LI->getLoopFor(&*thisInc);                                                          
              transformInc(L2, thisInc);                                                                   
    }

    void processLoopsToRoll(Module &M)
    {
      for(map<BasicBlock*, bool>::iterator bbc = blockCollapse.begin(); bbc!=blockCollapse.end(); ++bbc){
	if((*bbc).second == true){

	  map<BasicBlock*, string>::iterator bbrep = blockRepCond.find((*bbc).first);
	  if(bbrep!=blockRepCond.end())
	    if((*bbrep).second != ""){
	      

	      BasicBlock* thisBB = (*bbc).first;

	      ++NumLoopsRolled;

	      
	      map<BasicBlock*, BasicBlock*>::iterator bbInc = blockBodyInc.find(thisBB);
	      assert(bbInc!=blockBodyInc.end());
	      BasicBlock* thisInc = (*bbInc).second;
	      if(debugRollupLoops){
		errs() << "Found basicblock to collapse : " << thisBB->getName() <<"\n";
		for(BasicBlock::iterator bbi = thisBB->begin(); bbi!=thisBB->end(); ++bbi)
		  {
		    errs()<<"\t" << *bbi<<"\n";		  		
		  }
	      }
  
	      addDummyStart(M,thisBB,(*bbrep).second);
	      addDummyEnd(M,thisBB,(*bbrep).second);


	      Function* F = thisBB->getParent();
	      getLoopInfoAndModify(*F,thisBB,thisInc);

	    } //end of found basic block to collapse

	}



      } //blockCollapse iterator



    }



    void processFunctions(Function &F) {
      if(F.isDeclaration())
	errs() << "Fn is declaration\n";      
      
      //blockCollapse.clear();
      //blockRepCond.clear();
      
      //errs() << "Function: " << F.getName() << '\n';
      
      getLoopInfo(F);
      //errs() << "\n";
                 
    } // End runOnFunction
    

    bool runOnModule(Module &M){
      for(Module::iterator mIter = M.begin(); mIter != M.end(); ++mIter) {
	Function* F = &(*mIter);

	if(F && !F->isDeclaration()){
	  processFunctions((*F));
	}	
      }

      //add global variable to replace globalRepValue=0
      //Value* initGV = ConstantInt::get(Type::getInt64Ty(M.getContext()),0);
      //GlobalVariable *GV = new GlobalVariable(M,Type::getInt64Ty(M.getContext()), true, GlobalValue::ExternalLinkage, (Constant*) initGV, "globalRepValue");


      //process loops identified as candidates
      processLoopsToRoll(M);


      return true;
      }
   

  }; // End of struct FunctionRollupLoops
} // End of anonymous namespace

char RollupLoops::ID = 0;
static RegisterPass<RollupLoops> X("rollup-loops", "Do not unroll loops if possible.");
