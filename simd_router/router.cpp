//===----------------------------- router.cpp -----------------------------===//
// This file implements the fine-mesh router for the 
// Multi-SIMD tiled quantum architecture.
//
//                              Ali JavadiAbhari
//                            Princeton University
//                                January 2015
//
//===----------------------------------------------------------------------===//

// Usage:
// $ ./router </path/to/benchmark_name/without/.lpfs/.cg/suffix>
// $ ./router grovers   OR    $ ./router ../grovers

// 1/29:  Read LPFS
//        Get logical qubit moves/src/dest
//        Set up qbit data structure

// 1/30:  Pseudocode

// 2/2:   Finished Logical/Physical Inst structures

// 2/3:   Annotate Globabl memory partitions
//        Inject EPR
         
// 2/4:   Generate full physical instructions
//        Event Driven Simulator

// 2/5:   Group talk

// 3/8:   One-to-One marking of Log/Phys instructions

// 3/11:  Finished DAG-based ready queue

// 3/16:  Smoothing

// 4/27:  get_dag 20% speedup

// 5/4:   boost serialization (incomplete: segfault)

// 5/5:   Parse Coarse Grain Instructions

// 5/17:  ~1000X speedup by removing qbit->Instruction map

// 5/20:  Backward/Backforth Smoothing

// 6/11:  Window Smoothing

// 6/23:  Local Memory

// 6/25:  Call Graph

#include <iostream>   //std::cout, std::cin
#include <fstream>    //std::ifstream
#include <vector>
#include <algorithm>  //std::find_if, std::remove_if, std::sort
#include <stdlib.h>   //system
#include <cmath>      //sqrt
#include <time.h>     //clock
#include <sstream>    //std::stringstream
#include <map>
#include <unordered_map>
#include <tuple>
#include <stack>
#include <cstring>    //strcmp
#include <iterator>
#include <memory>     //std::shared_ptr, std::make_unique
#include <limits>     //std::numeric_limits
#include <boost/graph/graphviz.hpp>
#include <boost/archive/text_oarchive.hpp> //native binary output archive used for saving.
#include <boost/archive/text_iarchive.hpp> //native binary input archive used for loading.
#include <boost/serialization/export.hpp>  
#include <boost/serialization/vector.hpp>  //serialize vectors
#include <boost/graph/adj_list_serialize.hpp> //serialize adjacency list
#include <boost/serialization/unordered_map.hpp>  //serialize unordered_map
#include <boost/serialization/base_object.hpp> // serialize in polymorphism
#include <boost/serialization/shared_ptr.hpp> // serialize shared_ptr
#include <boost/graph/depth_first_search.hpp> // DFS for graph traversal

#define _DEBUG_PROGRESS
//#define _DEBUG_SERIALIZATION
//#define _DEBUG_PERFORMANCE
//#define _DEBUG_QBIT_MANAGEMENT
//#define _DEBUG_MESH
//#define _DEBUG_LOGICAL_INSTS
//#define _DEBUG_CG_INSTS
//#define _DEBUG_ANNOTATE_GLOBAL_MEM
//#define _DEBUG_EPR_INJECTION
//#define _DEBUG_PHYSICAL_INSTS
//#define _DEBUG_DAG
//#define _DEBUG_CALL_GRAPH
//#define _DEBUG_READY_QUEUE       
//#define _DEBUG_SIMULATOR

bool forward = false;
bool backward = false;

bool cap = false;
bool window = false;
unsigned int mov_cap = std::numeric_limits<unsigned int>::max();
unsigned int window_size = std::numeric_limits<unsigned int>::max();

bool report_usage = false;
bool report_ages = false;
bool report_storage = false;

//TODO: make these input arguments later
#define distribute 1
#define LEAF_SIMULATION_MAX 3

#define P_th 4              // Steane code threshold = 10^-4   
#define epsilon 0.5         // total desired logical error
unsigned long long total_logical_gates;    // KQ parameter, needed for calculating L_error_rate
double L_error_rate;           // desired logical error rate = 10^-(L_error_rate)
int P_error_rate;           // device error rate = 10^-(P_error_rate)
int concatenation_level;    // concatenation_level of Steane code


unsigned int SIMD_K=0;
unsigned int SIMD_D=0;
unsigned int num_magic_factories = 1;
unsigned int num_zero_factories = 1;
unsigned int num_epr_factories = 1;
unsigned int SIMD_rows=0;
unsigned int SIMD_cols=0;

std::vector<unsigned int> magic_factories;
std::vector<unsigned int> zero_factories;
std::vector<unsigned int> epr_factories;

const std::unordered_map<std::string, unsigned int> op_delays 
      = {{"PrepZ",1}, {"X",1}, {"Z",1}, {"H",1}, {"CNOT",10}, {"T",1}, {"Tdag",1}, {"S",1}, {"Sdag",1}, {"MeasZ",10}};

// For a given concatenation_level, what is the latency and area of creating one output from the ancilla factories?
unsigned int zero_delay = 0;
unsigned int epr_delay = 0;
unsigned int magic_delay = 0;
unsigned int zero_size = 0;
unsigned int epr_size = 0;
unsigned int magic_size = 0;

// ------------------------- qbit class --------------------------
enum sub_loc_t {G, T, TU_G, TU_T, L};
enum state_t {IDLE, IN_OP, IN_MOV};
enum kind_t {DATA, EPR1, EPR2, ZERO, MAGIC};

const char* sub_loc_to_string(sub_loc_t sub_loc) {
  switch (sub_loc)
  {
    case G:     return "G";
    case T:     return "T";
    case TU_G:  return "TU_G";
    case TU_T:  return "TU_T";
    case L:     return "L";
  }
}

const char* state_to_string(state_t state) {
  switch (state)
  {
    case IDLE:    return "idle";
    case IN_OP:   return "in_op";
    case IN_MOV:  return "in_mov";
  }
}

const char* kind_to_string(kind_t kind) {
  switch (kind)
  {
    case DATA:    return "data";
    case EPR1:    return "epr1";
    case EPR2:    return "epr2";
    case ZERO:    return "zero";
    case MAGIC:   return "magic";                  
  }
}

class qbit
{
  public:  
    std::string q_id;   // "a1","a2","a3",...
    kind_t kind;        // 0(data),1(epr1),2(epr2),3(zero),4(magic)
    unsigned int age;   // 0,1,2,3,...
    unsigned int loc;   // 1,2,3,...
    sub_loc_t sub_loc;  // 0(G), 1(T), 2(TU_G), 3(TU_T), 4(L)
    state_t state;      // 0(IDLE), 1(IN_OP), 2(IN_MOV)

    unsigned int op_time_remaining;
    unsigned int random;
    unsigned int destination;
    sub_loc_t destination_sub;
    
  public:
    qbit(){}
    qbit(std::string q_id, kind_t kind, int age, int loc, sub_loc_t sub_loc, state_t state)
      :q_id(q_id), kind(kind), age(age), loc(loc), sub_loc(sub_loc), state(state){
        #ifdef _DEBUG_QBIT_MANAGEMENT
          std::cout<<"creating qubit "<<q_id<<" in "<<loc<<"|"<<sub_loc_to_string(sub_loc)<<std::endl;
        #endif  
        insert_q_map();       // insert into q_map upon construction
      }

    ~qbit(){
      #ifdef _DEBUG_QBIT_MANAGEMENT
        std::cout<<"deleting qubit "<<q_id<<" in "<<loc<<"|"<<sub_loc_to_string(sub_loc)<<std::endl;
      #endif  
      erase_q_map();  // delete from q_map upon destruction
      record_age();
    }
    
    void status();

    void insert_q_map();

    void erase_q_map();

    void record_age();
};

// --------------------------- Instruction Classes -------------------------
// Instruction: parent for MOVinst/OPinst (logical) and BMOVinst (physical)
class Instruction
{
  private:
    friend class boost::serialization::access; // In order to make Instruction serializable
    template<class Archive>
      void serialize(Archive & ar, const unsigned int version) {
        // choose what class fields you want to serialize
        ar & seq;
        ar & ts;
        ar & qid;
        ar & is_executing;
        ar & is_complete;
      }

  public:
    unsigned int seq;
    unsigned int ts; 
    std::vector<std::string> qid;  

    bool is_executing;  
    bool is_complete;
    bool no_child;

    Instruction(){} // Default ctor
    Instruction(unsigned int timestep, std::vector<std::string> qbit_id) // Ctor with timestep and qbit args
      : seq(0), ts(timestep), qid(qbit_id), is_executing(false), is_complete(false), no_child(false) {}

    virtual ~Instruction() {} // virtual destructor to allow dynamic_cast (b/c only allowed in polymorphism)

    virtual void print() {};
};

// CGinst: Coarse-Grain instruction (non-leaf)
class CGinst : public Instruction
{
  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version){
      // invoke base class serialization
      ar & boost::serialization::base_object<Instruction>(*this);      
      // choose what class fields you want to serialize
      ar & module_name;
      ar & is_leaf;  
      ar & ordered_leaf_names;
    }

  public:
    std::string module_name;
    bool is_leaf;
    std::vector<std::string> ordered_leaf_names;  // all leaves that will execute if this were the call graph root

    CGinst(){} // Default ctor
    CGinst(unsigned int timestep, std::string module_name, std::vector<std::string> qreg, bool is_leaf)
      : Instruction(timestep, qreg), module_name(module_name), is_leaf(is_leaf) {}

    void print();
};
BOOST_CLASS_EXPORT_GUID(CGinst, "CGinst")

// Operation -- Logical or Physical
class OPinst : public Instruction
{
  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version){
      // invoke base class serialization
      ar & boost::serialization::base_object<Instruction>(*this);
      // choose what class fields you want to serialize
      ar & zone;
      ar & op_type;
    }

  public:
    unsigned int zone;
    std::string op_type;

    OPinst(){} // Default ctor
    OPinst(unsigned int timestep, unsigned int zone, 
           std::string operation_type, std::vector<std::string> qbit_id) // Ctor accepting all data
      : Instruction(timestep, qbit_id), zone(zone), op_type(operation_type) {}

    ~OPinst() {}

    void print();
};
BOOST_CLASS_EXPORT_GUID(OPinst, "OPinst")


// Ballistic Move -- Physical
class BMOVinst : public Instruction
{
  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version){
      // invoke base class serialization
      ar & boost::serialization::base_object<Instruction>(*this);
      // choose what class fields you want to serialize
      ar & src;
      ar & src_sub;  
      ar & dest;
      ar & dest_sub;
    }
  
  public:
    unsigned int src;
    sub_loc_t src_sub;
    unsigned int dest;
    sub_loc_t dest_sub;

    BMOVinst(){} // Default ctor
    BMOVinst(unsigned int timestep, 
              unsigned int source, sub_loc_t source_sub, 
              unsigned int destination, sub_loc_t destination_sub, 
              std::vector<std::string> qbit_id) // Ctor accepting all data
      : Instruction(timestep, qbit_id), src(source), src_sub(source_sub), dest(destination), dest_sub(destination_sub){}

    ~BMOVinst() {}

    void print();
};
BOOST_CLASS_EXPORT_GUID(BMOVinst, "BMOVinst")


// Teleportation Move -- Logical
class MOVinst : public Instruction
{
  public:
  unsigned int src, dest;

  MOVinst(){} // Default ctor
  MOVinst(unsigned int timestep, unsigned int source, unsigned int destination, std::vector<std::string> qbit_id) // Ctor accepting all data
    :  Instruction(timestep, qbit_id), src(source), dest(destination) {}

  ~MOVinst() {}

  void print();
};

void CGinst::print() {
  std::cout<<seq<<' ';  
  std::cout<<module_name;
  for (auto &r : qid)
    std::cout<<' '<<r;
  std::cout<<std::endl;
}

void OPinst::print() {
  std::cout<<seq<<' ';
  std::cout<<op_type<<' '<<zone<<' '<<qid[0];
  if(qid.size()>1)
    std::cout<<' '<<qid[1]<<std::endl;
  else
    std::cout<<std::endl;   
  return;
}

void MOVinst::print() {
  std::cout<<seq<<' ';  
  std::cout<<"MOV "
           <<src<<' '
           <<dest<<' '
           <<qid[0]
           <<std::endl;   
  return;
}

void BMOVinst::print() {
  std::cout<<seq<<' ';  
  std::cout<<"BMOV "
           <<src<<'|'<<sub_loc_to_string(src_sub)
           <<' '<<dest<<'|'<<sub_loc_to_string(dest_sub)
           <<' '<<qid[0]
           <<std::endl;   
  return;
}

// -------------------------------- Typedefs ---------------------------------

typedef std::vector<std::shared_ptr<Instruction>> InstVecTy;
typedef std::unordered_map<std::string, InstVecTy> InstTableTy;
typedef std::tuple<unsigned int, unsigned int, unsigned int> SeqTupleTy;
typedef std::vector<SeqTupleTy> SeqTupleVecTy;
typedef std::unordered_map<std::string, SeqTupleVecTy> SeqTupleTableTy;
typedef std::tuple<std::shared_ptr<Instruction>,std::shared_ptr<Instruction>,std::shared_ptr<Instruction>> InstTupleTy;
typedef std::vector<InstTupleTy> InstTupleVecTy;
typedef std::unordered_map<std::string, InstTupleVecTy> InstTupleTableTy;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS> DAGTy;
typedef boost::adjacency_list<boost::setS, boost::vecS, boost::bidirectionalS> CallGraphTy;
typedef std::unordered_map<std::string, DAGTy> DAGTableTy;
typedef std::unordered_map<std::string, CallGraphTy> CallGraphTableTy;
typedef DAGTy::vertex_descriptor VertexID;
typedef DAGTy::edge_descriptor EdgeID;


// ------------------------- Global Data Structures ---------------------------

// structures to keep track of live qbits
std::map<std::string, qbit*> q_map;
// Result vectors for experiments
unsigned long long cycle = 0;
unsigned long long current_leaf_cycle = 0;
std::map<std::string, std::vector<unsigned long long>> leaf_cycles; // a leaf may execute a few times, we average its length.
std::map<unsigned int, unsigned int> qbits_per_cycle;
std::map<unsigned int, SeqTupleTy> ancillas_per_cycle;
std::map<unsigned int, unsigned int> zeros_per_cycle;
std::map<unsigned int, unsigned int> eprs_per_cycle;
std::map<unsigned int, unsigned int> magics_per_cycle;
std::map<std::string, unsigned int>  qbit_ages;
std::map<unsigned int, std::vector<unsigned int>> storage_per_cycle;
// structures to keep track of leaf/non-leaf modules
std::vector<std::string> leaves;
std::vector<std::string> non_leaves;
std::vector<std::string> all_ordered_leaf_names; // all leaves in order of execution
std::vector<std::string> sub_ordered_leaf_names; // small subset of leaves, respecting same order of execution
std::stack<std::string> call_stack;             // call stack to keep track of module calls in call graph traversal  


void qbit::insert_q_map() {
  if (q_map.find(q_id)!=q_map.end()) {
    std::cerr<<"Error: Attempting to create an existing qubit: "<<q_id<<std::endl;
    return;
  }
  else {
    #ifdef _DEBUG_QBIT_MANAGEMENT
      std::cout<<"Inserted in q_map: "<<q_id<<std::endl;
    #endif
      q_map[q_id] = this;    
  }
  return;
}

void qbit::erase_q_map() {
  q_map.erase(q_id);
  return;
}

void qbit::record_age() {
  qbit_ages[q_id] = age;
  return;
}

void qbit::status() {
  std::cout<<"qubit "<<q_id<<":"
           <<"\n\tkind: "<<kind_to_string(kind)
           <<", age: "<<age  
           <<", loc: "<<loc
           <<", subloc: "<<sub_loc_to_string(sub_loc)
           <<", state: "<<state_to_string(state)
           <<std::endl;
  return;
}

// ------------------------- Print Helper Functions ---------------------------

// which Instruction in PhysicalInst_v has the sequence tag of seq?
std::shared_ptr<Instruction> which_Instruction(InstVecTy &PhysicalInst_v, unsigned int seq) {
  auto I = std::find_if(PhysicalInst_v.begin(), PhysicalInst_v.end(), [&seq](const std::shared_ptr<Instruction> &i){return i->seq == seq;});
  if (I == PhysicalInst_v.end()) {
    std::cerr<<"seq #"<<seq<<" not found in PhysicalInst_v."<<std::endl;
    exit(1);
  }
  return *I;  
}

// print Instruction* vector contents
void print_logical_insts(const InstTableTy &mapOfLogicalInst_v) {
  std::cout<<"---------------- Printing Logical Instructions ---------------"<<std::endl;
  for (auto &F : mapOfLogicalInst_v) {
    std::cout<<F.first<<" (Leaf)"<<std::endl;
    InstVecTy LogicalInst_v;
    LogicalInst_v = F.second;
    for (auto &I : LogicalInst_v) {
      if (std::shared_ptr<MOVinst> l_mov_inst= std::dynamic_pointer_cast<MOVinst>(I))
        l_mov_inst->print();
      else if (std::shared_ptr<BMOVinst> l_bmov_inst= std::dynamic_pointer_cast<BMOVinst>(I))  //local mem moves
        l_bmov_inst->print();
      else if (std::shared_ptr<OPinst> l_op_inst= std::dynamic_pointer_cast<OPinst>(I))
        l_op_inst->print();   
    }
    std::cout<<std::endl;
  }
}

// print CGinst* vector contents
void print_cg_insts (const InstTableTy &mapOfCGInst_v) {
  std::cout<<"---------------- Printing Coarse-Grain Instructions ---------------"<<std::endl;
  for (auto &F : mapOfCGInst_v) {
    std::cout<<F.first<<" (Nonleaf)"<<std::endl;
    InstVecTy CGInst_v;
    CGInst_v = F.second;
    for (auto &I : CGInst_v) {
      if (std::shared_ptr<CGinst> cg_inst= std::dynamic_pointer_cast<CGinst>(I))
        cg_inst->print();
    }
    std::cout<<std::endl;
  }
}

void print_physical_insts(const InstTableTy &mapOfPhysicalInst_v) {
  std::cout<<"---------------- Printing Physical Instructions ---------------"<<std::endl;
  for (auto &F : mapOfPhysicalInst_v) {    
    std::cout<<F.first<<" (Leaf)"<<std::endl;
    InstVecTy PhysicalInst_v;
    PhysicalInst_v = F.second;    
    for (auto &I : PhysicalInst_v) {
      if (std::shared_ptr<BMOVinst> ph_mov_inst= std::dynamic_pointer_cast<BMOVinst>(I))
        ph_mov_inst->print();  
      else if (std::shared_ptr<OPinst> ph_op_inst= std::dynamic_pointer_cast<OPinst>(I))
        ph_op_inst->print();     
    }
    std::cout<<std::endl;
  }
}

void print_dag(const DAGTableTy &mapOfDAG) {
  for (auto &F : mapOfDAG) {
    std::cout<<"------------------------ Printing DAG -------------------------"<<std::endl;      
    std::cout<<F.first<<" (Leaf)"<<std::endl;
    DAGTy g;   
    g = F.second;        
    boost::write_graphviz(std::cout, g);
    std::cout<<std::endl;
  }
}

void print_call_graph(const CallGraphTableTy &mapOfCallGraph, InstTableTy &mapOfCGInst_v) {
  std::cout<<"--------------------- Printing Call Graph ----------------------"<<std::endl;
  for (auto &F : mapOfCallGraph) {
    std::cout<<F.first<<" (Nonleaf)"<<std::endl;
    InstVecTy CGInst_v = mapOfCGInst_v[F.first];
    CallGraphTy g;   
    g = F.second;     
    for (int seq=0; seq<boost::num_vertices(g); seq++) {
      std::shared_ptr<Instruction> I = which_Instruction(CGInst_v, seq);
      std::shared_ptr<CGinst> cg_inst = std::dynamic_pointer_cast<CGinst>(I);
      cg_inst->print();
    }
    std::cout<<std::endl;
  }
}

void print_mesh() {
  std::cout<<"----------------------- Printing Mesh --------------------------"<<std::endl;
  bool done_printing = false;
  for (unsigned int r=0; r<2*SIMD_rows; r++) {
    if (done_printing) break;    
    for (unsigned int c=0; c <SIMD_cols; c++) {
      if (done_printing) break;      
      unsigned int zone_num = r*(SIMD_cols)+c+1; 
      if (zone_num >= SIMD_K) done_printing = true;     
      if (c%2==0) std::cout<<"T"<<zone_num;
      else std::cout<<"G"<<zone_num;
      if (std::find(magic_factories.begin(),magic_factories.end(),zone_num) != magic_factories.end()) std::cout<<" (magic)";
      else if (std::find(epr_factories.begin(),epr_factories.end(),zone_num) != epr_factories.end()) std::cout<<" (epr)";
      else if (std::find(zero_factories.begin(),zero_factories.end(),zone_num) != zero_factories.end()) std::cout<<" (zero)";
      std::cout<<"\t";      
    }
    done_printing = false;
    std::cout<<std::endl;
    for (unsigned int c=0; c <SIMD_cols; c++) {   
      if (done_printing) break;            
      unsigned int zone_num = r*(SIMD_cols)+c+1;
      if (zone_num >= SIMD_K) done_printing = true;           
      if (c%2==0) std::cout<<"G"<<zone_num;
      else std::cout<<"T"<<zone_num;      
      if (std::find(magic_factories.begin(),magic_factories.end(),zone_num) != magic_factories.end()) std::cout<<" (magic)";
      else if (std::find(epr_factories.begin(),epr_factories.end(),zone_num) != epr_factories.end()) std::cout<<" (epr)";
      else if (std::find(zero_factories.begin(),zero_factories.end(),zone_num) != zero_factories.end()) std::cout<<" (zero)";
      std::cout<<"\t";       
    }
    std::cout<<std::endl;
  }
}


// ------------------------- Parser and Parser Helpers ------------------------

// is there any of several words in a given string?
template<typename T, size_t N>
T * end(T (&ra)[N]) {
    return ra + N;
}

std::string::size_type is_there(std::vector<std::string> needles, std::string haystack) {
  std::vector<std::string>::iterator needle;
  std::string::size_type pos;
  for(needle=needles.begin(); needle!=needles.end(); ++needle){
    pos = haystack.find(*needle);
    if(pos != std::string::npos){
      return pos;
    }
  }  
  return std::string::npos;
}

// tokenize string
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
      if (!item.empty())             
        elems.push_back(item);
    }
    return elems;
}

// Set the latency and size of ancilla factories
void create_factories(unsigned int concatenation_level) {
  unsigned int magic_cnots= 0;
  for (int i=0;i<concatenation_level;i++)
    magic_cnots += ( (unsigned int)pow(7.0,(double)(i)) - 1 );

  zero_delay = op_delays.find("PrepZ")->second + concatenation_level * ( 3*op_delays.find("H")->second+9*op_delays.find("CNOT")->second+1*op_delays.find("MeasZ")->second );
  epr_delay  = op_delays.find("PrepZ")->second + (concatenation_level+1) * ( 1*op_delays.find("H")->second+1*op_delays.find("CNOT")->second);;
  magic_delay= op_delays.find("PrepZ")->second + (concatenation_level+1) * ( 
                2*op_delays.find("H")->second+4*op_delays.find("CNOT")->second+1*op_delays.find("T")->second+1*op_delays.find("MeasZ")->second) + magic_cnots*op_delays.find("CNOT")->second;

  zero_size  = (concatenation_level==0) ? 1 : 8*(unsigned int)pow(7.0,(double)(concatenation_level-1));
  epr_size  = (concatenation_level==0) ? 2 : 2*8*(unsigned int)pow(7.0,(double)(concatenation_level-1));
  magic_size = (concatenation_level==0) ? 2 : 2*8*(unsigned int)pow(7.0,(double)(concatenation_level-1));

  std::cerr<<"zero_delay:"<<zero_delay<<" epr_delay:"<<epr_delay<<" magic_delay:"<<magic_delay<<" zero_size:"<<zero_size<<" epr_size:"<<epr_size<<" magic_size:"<<magic_size<<std::endl;
}

// Create mesh based on number of compute regions and number of zero/magic/epr regions
// optimizations of layout help latency
InstTableTy create_mesh(const InstTableTy &mapOfLogicalInst_v) {
  #ifdef _DEBUG_PROGRESS
    std::cerr<<"creating mesh..."<<std::endl;
  #endif

  unsigned int num_compute_regions = SIMD_K - (num_zero_factories + num_epr_factories + num_magic_factories);  

  // naive layout: just append ancilla factories as extra regions
  for (unsigned int i=0; i<num_zero_factories; i++)
    zero_factories.push_back(num_compute_regions+i+1);
  for (unsigned int i=0; i<num_epr_factories; i++)
    epr_factories.push_back(num_compute_regions+num_zero_factories+i+1);
  for (unsigned int i=0; i<num_magic_factories; i++)
    magic_factories.push_back(num_compute_regions+num_zero_factories+num_epr_factories+i+1);

  InstTableTy mapOfLogicalInst_v2;                               
  std::string leaf_func;   
  for (auto &F : mapOfLogicalInst_v) {
    leaf_func = F.first;
    InstVecTy LogicalInst_v;
    LogicalInst_v = F.second;
    for (auto &I : LogicalInst_v) {
      if (std::shared_ptr<MOVinst> l_mov_inst = std::dynamic_pointer_cast<MOVinst>(I)) {
        mapOfLogicalInst_v2[leaf_func].push_back(I);
      }
      if (std::shared_ptr<BMOVinst> ph_mov_inst = std::dynamic_pointer_cast<BMOVinst>(I)) {
        mapOfLogicalInst_v2[leaf_func].push_back(I);
      }  
      if (std::shared_ptr<OPinst> l_op_inst = std::dynamic_pointer_cast<OPinst>(I)) {
        mapOfLogicalInst_v2[leaf_func].push_back(I);
      }  
    }
  }
  return mapOfLogicalInst_v2;   
}

// Read in LPFS schedule, parse: SIMD_K, MOVinst(qbits, src, dest), OPinst(optype)
InstTableTy parse_LPFS_file (const std::string file_path) {
  #ifdef _DEBUG_PROGRESS
    std::cerr<<"parsing LPFS..."<<std::endl;
  #endif
  
  std::ifstream LPFSfile (file_path);
  std::string line;
  InstTableTy mapOfLogicalInst_v;
  std::string leaf_func;    
  
  if (LPFSfile.is_open())
  {    
    // op_strings -- can be expanded to other 1-qubit / 2-qubit gates
    const char* all_ops[] = {"PrepZ ", "X ", "Z ", "H ", "CNOT ", "T ", "Tdag ", "S ", "Sdag ", "MeasZ "};
    std::vector<std::string> op_strings(all_ops, end(all_ops));
     
    while ( std::getline (LPFSfile,line) ) {
      // FunctionHeaders
      if (line.find("Function") != std::string::npos) {
        std::vector<std::string> elems;
        split(line, ' ', elems);
        leaf_func = elems[1];
        if (SIMD_K == 0) {
          SIMD_K = atoi(elems[5].c_str()) + num_zero_factories + num_magic_factories + num_epr_factories;
          SIMD_D = atoi(elems[7].c_str());    
          SIMD_rows = (int)(ceil(sqrt(2.0*SIMD_K)));
          SIMD_cols = (int)(ceil(2.0*SIMD_K/SIMD_rows));
          std::cerr<<"Topology : SIMD("<<SIMD_K<<","<<SIMD_D<<") : "<<SIMD_rows<<"*"<<SIMD_cols<<std::endl;
        }        
      }

      // MOVinst : Teleport
      else if (line.find("TMOV") != std::string::npos) {
        std::vector<std::string> elems;
        split(line, ' ', elems);
        unsigned int timestep = atoi(elems[0].substr(0,elems[0].find(',')).c_str());
        unsigned int source = atoi(elems[3].c_str());
        unsigned int destination = atoi(elems[2].c_str());
        std::vector<std::string> qbit_id;        
        std::string qbit_id1 = elems[4];
        qbit_id.push_back(qbit_id1);
        std::shared_ptr<Instruction> logical_inst_cur (std::make_shared <MOVinst> (timestep, source, destination, qbit_id));        
        mapOfLogicalInst_v[leaf_func].push_back(logical_inst_cur);
      }      

      // BMOVinst: Local Memory Move
      else if (line.find("BMOV") != std::string::npos) {
        std::vector<std::string> elems;
        split(line, ' ', elems);
        unsigned int timestep = atoi(elems[0].substr(0,elems[0].find(',')).c_str());
        unsigned int source = atoi(elems[3].c_str());
        unsigned int destination = atoi(elems[2].c_str());
        sub_loc_t source_sub, destination_sub;
        if (source % 10 == 0) {
          source = (unsigned int)(source / 10);
          source_sub = L;
          destination_sub = T;
        }
        else if (destination % 10 == 0) {
          destination = (unsigned int)(destination / 10);
          destination_sub = L;
          source_sub = T;
        }
        else
          std::cerr<<"Error: incorrect Local Memory move."<<std::endl;
        std::vector<std::string> qbit_id;        
        std::string qbit_id1 = elems[4];
        qbit_id.push_back(qbit_id1);
        std::shared_ptr<Instruction> logical_inst_cur (std::make_shared <BMOVinst> (timestep, source, source_sub, destination, destination_sub, qbit_id));        
        mapOfLogicalInst_v[leaf_func].push_back(logical_inst_cur);
      } 

      // OPinsts
      else if (is_there(op_strings, line) != std::string::npos) {
        std::vector<std::string> elems;
        std::vector<std::string> ts_zone;        
        split(line, ' ', elems);
        split(elems[0], ',', ts_zone);
        unsigned int timestep = atoi(ts_zone[0].c_str());
        unsigned int zone = atoi(ts_zone[1].c_str());
        std::string operation_type = elems[1];
        std::vector<std::string> qbit_id;
        std::string qbit_id1 = elems[2]; 
        qbit_id.push_back(qbit_id1);        
        if (elems.size() == 4) {
          std::string qbit_id2 = elems[3];          
          qbit_id.push_back(qbit_id2);
        }        
        std::shared_ptr<Instruction> logical_inst_cur (std::make_shared <OPinst> (timestep, zone, operation_type, qbit_id));
        if (operation_type != "X" && operation_type != "Z" && operation_type != "T" && operation_type != "Tdag")
          mapOfLogicalInst_v[leaf_func].push_back(logical_inst_cur);        
      }
    }
    LPFSfile.close();
  }
  else {
    std::cerr<<"Error: Unable to open file."<<std::endl;
    exit(1);
  }

  return mapOfLogicalInst_v; 
}

// parse profile of module frequencies
std::unordered_map<std::string, unsigned long long> parse_freq (const std::string file_path) {
  #ifdef _DEBUG_PROGRESS
    std::cerr<<"parsing Freq..."<<std::endl;
  #endif  
  std::ifstream profile_freq_file (file_path);
  std::string line;
  std::string module_name = "";  
  std::unordered_map<std::string, unsigned long long> mapOfFreq;
  mapOfFreq.clear();
  unsigned long long freq = 0;
  if (profile_freq_file.is_open()) {
    while ( getline (profile_freq_file,line) ) {
      std::vector<std::string> elems; 
      elems.clear();     
      split(line, ' ', elems);              
      module_name = elems[0];
      freq = std::stoull(elems[9]);
      mapOfFreq[module_name] = freq;
    }
  }
  else {
    std::cerr << "Unable to open .freq file" << std::endl;
    exit(1);
  }
  return mapOfFreq;
}

// Read in CG schedule, parse: CGinst (name, width, length, is_leaf)
// and add instruction sequence tags
InstTableTy parse_CG_file (const std::string file_path) {
  #ifdef _DEBUG_PROGRESS
    std::cerr<<"parsing CG..."<<std::endl;
  #endif
  
  std::ifstream CGfile (file_path);
  std::string line;
  InstTableTy mapOfCGInst_v;
  InstVecTy bodyInst_v;
  unsigned int seq = 0;          // tag Instruction sequence numbers as we build bodyInst_v


  if (CGfile.is_open())
  {
    while ( std::getline (CGfile,line) ) {
      if(line.empty())
        continue;
      // Function Bodies: CGInsts
      else if (line.find("#") == std::string::npos && line.find("SIMD") == std::string::npos) {
        
        std::vector<std::string> elems;
        split(line, ' ', elems);
        std::string callee_name = elems[0];
        if (callee_name.find("llvm.")!=std::string::npos)
          continue;
        std::vector<std::string> args;
        for (int i = 2; i<elems.size(); i++)
          args.push_back( elems[i].substr(0, elems[i].find_first_of('(') ) );
        bool is_callee_leaf;
        if (std::find(leaves.begin(), leaves.end(), callee_name) != leaves.end())
          is_callee_leaf = false;
        else
          is_callee_leaf = true;
        std::shared_ptr<Instruction> cg_inst_cur (std::make_shared <CGinst> (0, callee_name, args, is_callee_leaf));
        cg_inst_cur->seq = seq++;
        bodyInst_v.push_back(cg_inst_cur);
      }
      // Function Summary
      else if (line.find("SIMD k") != std::string::npos) {
        std::vector<std::string> elems;
        split(line, ' ', elems);
        std::string module_name = elems[3];
        std::string leaf_status = elems[5];
        bool is_leaf = atoi(leaf_status.c_str());
        if (is_leaf == 0) {
          if (bodyInst_v.size() == 0) {
            std::cerr<<"Error: Bad .cg file format."<<std::endl;
            exit(1);
          }
          mapOfCGInst_v[module_name] = bodyInst_v;
          non_leaves.push_back(module_name);  // add module to non-leaves
        }
        else 
          leaves.push_back(module_name);      // add module to leaves
        bodyInst_v.clear();                    // all of this non-leaf's body is read
        seq = 0;
      }      
    }

    CGfile.close();
  }
  else {
    std::cerr<<"Error: Unable to open file."<<std::endl;
    exit(1);
  }
  return mapOfCGInst_v;
}

// --------------------- LogicalInst 2 PhysicalInst -----------------------

// Determine and track location in Distributed Global Memory
// and modify all instructions that have "0" for global memory 
InstTableTy annotate_global_memory (const InstTableTy &mapOfLogicalInst_v) {

  #ifdef _DEBUG_PROGRESS
    std::cerr<<"annotating global memory..."<<std::endl;
  #endif
    
  InstTableTy mapOfPhysicalInst_v;                               
  std::string leaf_func;  
  std::unordered_map<std::string, unsigned int> G_map; // maps q_id to location. last location determines which G partition   

  for (auto &F : mapOfLogicalInst_v) {
    leaf_func = F.first;
    InstVecTy LogicalInst_v;
    LogicalInst_v = F.second;
    G_map.clear();
    for (auto &I : LogicalInst_v) {
      if (std::shared_ptr<MOVinst> l_mov_inst = std::dynamic_pointer_cast<MOVinst>(I)) {
        unsigned int timestep = l_mov_inst->ts;
        unsigned int source = l_mov_inst->src;
        unsigned int destination = l_mov_inst->dest;
        std::vector<std::string> q_id = l_mov_inst->qid;
        sub_loc_t source_sub, destination_sub;
        #ifdef _DEBUG_ANNOTATE_GLOBAL_MEM
          std::cout<<"MOV "<<source<<' '<<destination<<' '<<q_id<<std::endl;              
        #endif
        if (source==0) {
          source_sub = G; // GlobalMem-to-Tile
          destination_sub = T;
          if(G_map.find(q_id[0])!=G_map.end()) {
            unsigned int G_cur = G_map.find(q_id[0])->second;
            source = G_cur; // source is not 0 anymore, but rather the current qbit.loc
            G_map[q_id[0]] = destination;
            if (G_cur==destination) {
              #ifdef _DEBUG_ANNOTATE_GLOBAL_MEM
                std::cout<<q_id[0]<<" staying in "<<G_cur<<std::endl; 
              #endif
            }
            else {
              #ifdef _DEBUG_ANNOTATE_GLOBAL_MEM
                std::cout<<q_id[0]<<" going from G"<<G_cur<<" to G"<<destination<<std::endl; 
              #endif
            }
          }
          else {
            G_map[q_id[0]] = destination;  // for fresh qubits map to where they are used first (first-touch)          
            source = destination; // source is not 0 anymore, but rather the destination loc
            #ifdef _DEBUG_ANNOTATE_GLOBAL_MEM 
              std::cout<<q_id[0]<<" initially mapped to G"<<destination<<std::endl; 
            #endif
          } 
        }
        else if (destination==0) {
          source_sub = T; // Tile-to-GlobalMem
          destination_sub = G;
          if(G_map.find(q_id[0])!=G_map.end()) {
            int G_cur = G_map.find(q_id[0])->second;
            destination = G_cur;  // destination is not 0 anymore, but rather qbit.loc
            #ifdef _DEBUG_ANNOTATE_GLOBAL_MEM
              std::cout<<q_id[0]<<" already in G"<<G_cur<<std::endl; 
            #endif
          }
          else {
            std::cerr<<"Error: "<<q_id[0]<<" NOT TRACKED!"<<std::endl; 
          }
        }
        else {
          source_sub = T; // Tile-to-Tile
          destination_sub = T;
          if(G_map.find(q_id[0])!=G_map.end()) {
            unsigned int G_cur = G_map.find(q_id[0])->second;
            if (G_cur!=source)  // if tracked right, it should be at "source" now.
              std::cerr<<"Error: "<<q_id[0]<<" MAL-TRACKED!"<<std::endl; 
            G_map[q_id[0]] = destination;
            #ifdef _DEBUG_ANNOTATE_GLOBAL_MEM 
              std::cout<<q_id[0]<<" moved from G"<<G_cur<<" to G"<<destination<<std::endl; 
            #endif
          }
          else {
            std::cerr<<"Error: "<<q_id[0]<<" NOT TRACKED!"<<std::endl; 
          }
        }
        std::shared_ptr<Instruction> physical_inst_cur (new BMOVinst(timestep, source, source_sub, destination, destination_sub, q_id));
        mapOfPhysicalInst_v[leaf_func].push_back(physical_inst_cur);
      }
      if (std::shared_ptr<BMOVinst> ph_mov_inst = std::dynamic_pointer_cast<BMOVinst>(I)) {
        mapOfPhysicalInst_v[leaf_func].push_back(I);
      }  
      if (std::shared_ptr<OPinst> l_op_inst = std::dynamic_pointer_cast<OPinst>(I)) {
        mapOfPhysicalInst_v[leaf_func].push_back(I);
      }  
    }
  }
  return mapOfPhysicalInst_v;   
}

// Add QEC redundancies for level-l Steane error correction
// 7^l physical qubits for each logical qubit (increases channel widths, increases mesh size -- actual qubits not added)
// 7^l physical gates of transversal PrepZ/H/CNOT/S/Sdag/MeasZ gates for each logical gate (can all be done in parallel so actual gates not added. X/Z: software)
// one encoded magic state for each T gate, plus BMOV to bring it onto data
// three encoded zero state for each data region
InstTableTy add_qec (const InstTableTy &mapOfPhysicalInst_v) {

  #ifdef _DEBUG_PROGRESS
  std::cerr<<"adding QEC..."<<std::endl;
  #endif

  InstTableTy mapOfPhysicalInst_v2;                                      
  std::string leaf_func;   
  unsigned long int zero_count;  
  zero_count = 0;  
  for (auto &F : mapOfPhysicalInst_v) {  
    leaf_func = F.first;
    InstVecTy PhysicalInst_v = F.second;
    InstVecTy PhysicalInst_v2;      // construct this new vector, to avoid lots of inserts
    PhysicalInst_v2.clear();

    unsigned int seq = 0;          // tag Instruction sequence numbers as we build PhysicalInst_v2
    unsigned int timestep;
    unsigned int qec_loc;
    sub_loc_t qec_sub_loc;
    std::vector<std::string> qec_id; 
    bool qec_necessary = false;  
    for (auto &I : PhysicalInst_v) {
      if (std::shared_ptr<BMOVinst> ph_mov_inst = std::dynamic_pointer_cast<BMOVinst>(I)) {            
        timestep = ph_mov_inst->ts;        
        qec_loc = ph_mov_inst->dest;
        qec_sub_loc = ph_mov_inst->dest_sub;   
        qec_id = ph_mov_inst->qid;     
        unsigned int source_sub = ph_mov_inst->src_sub;   
        PhysicalInst_v2.push_back(I);
        I->seq = seq++; 
        // if this is a local memory move, no qec injection is necessary
        if (source_sub == L || qec_sub_loc == L) {
          qec_necessary = false;          
        }           
        else 
          qec_necessary = true;           
      }      
      else if (std::shared_ptr<OPinst> ph_op_inst = std::dynamic_pointer_cast<OPinst>(I)) {
        timestep = ph_op_inst->ts;        
        qec_loc = ph_op_inst->zone;
        qec_sub_loc = G;
        qec_id = ph_op_inst->qid;        
        std::string op_type = ph_op_inst->op_type;
        PhysicalInst_v2.push_back(I);
        I->seq = seq++;                
        qec_necessary = true;
      } 
      if (qec_necessary) {
        for (int i=0; i<qec_id.size(); i++) {
          // Zero naming scheme: 
          // In order to not create false dependency between Zero qubits
          // add a unique "zero_count" number to the end of the name 
          zero_count++;
          int factory_loc = zero_count % zero_factories.size();
          std::string zero_id1 = qec_id[0]+"_zero1_"+std::to_string(zero_count);
          std::string zero_id2 = qec_id[0]+"_zero2_"+std::to_string(zero_count);
          //std::string zero_id3 = qec_id[0]+"_zero3_"+std::to_string(zero_count);
          //std::string zero_id4 = qec_id[0]+"_zero4_"+std::to_string(zero_count);

          // placeholder vector for q_id arguments
          std::vector<std::string> qbit_id (2); 

          // BMOV to QEC location
          qbit_id = {zero_id1};                
          std::shared_ptr<Instruction> physical_inst_1 (std::make_shared <BMOVinst> (timestep,zero_factories[factory_loc],  T, qec_loc,  (sub_loc_t)(qec_sub_loc),     qbit_id));
          qbit_id = {zero_id2};                       
          std::shared_ptr<Instruction> physical_inst_2 (std::make_shared <BMOVinst> (timestep,zero_factories[factory_loc],  T, qec_loc,  (sub_loc_t)(qec_sub_loc),     qbit_id));
          /* uncomment for conservative qec
          qbit_id = {zero_id3};  
          std::shared_ptr<Instruction> physical_inst_3 (std::make_shared <BMOVinst> (timestep,zero_factories[factory_loc],  T, qec_loc,  (sub_loc_t)(qec_sub_loc),     qbit_id));          
          qbit_id = {zero_id4};
          std::shared_ptr<Instruction> physical_inst_4 (std::make_shared <BMOVinst> (timestep,zero_factories[factory_loc],  T, qec_loc,  (sub_loc_t)(qec_sub_loc),     qbit_id));          
          */

          // QEC Operations
          qbit_id = {zero_id1};        
          std::shared_ptr<Instruction> physical_inst_5 (std::make_shared <OPinst> (timestep,qec_loc,     "H",    qbit_id ));          
          qbit_id = {qec_id[i],zero_id1};
          std::shared_ptr<Instruction> physical_inst_6 (std::make_shared <OPinst> (timestep,qec_loc,     "CNOT", qbit_id ));
          qbit_id = {zero_id1};                
          std::shared_ptr<Instruction> physical_inst_7 (std::make_shared <OPinst> (timestep,qec_loc,     "MeasZ",qbit_id ));
          qbit_id = {zero_id2};        
          std::shared_ptr<Instruction> physical_inst_8 (std::make_shared <OPinst> (timestep,qec_loc,     "CNOT", qbit_id ));          
          qbit_id = {zero_id2,qec_id[i]};
          std::shared_ptr<Instruction> physical_inst_9 (std::make_shared <OPinst> (timestep,qec_loc,     "H",    qbit_id ));
          qbit_id = {zero_id2};                
          std::shared_ptr<Instruction> physical_inst_10 (std::make_shared <OPinst>(timestep,qec_loc,     "MeasZ",qbit_id ));
          /*  uncomment for conservative qec
          qbit_id = {zero_id3};        
          std::shared_ptr<Instruction> physical_inst_11 (std::make_shared <OPinst> (timestep,qec_loc,    "H",    qbit_id ));          
          qbit_id = {qec_id[i],zero_id3};
          std::shared_ptr<Instruction> physical_inst_12 (std::make_shared <OPinst> (timestep,qec_loc,    "CNOT", qbit_id ));
          qbit_id = {zero_id3};                
          std::shared_ptr<Instruction> physical_inst_13 (std::make_shared <OPinst> (timestep,qec_loc,    "MeasZ",qbit_id ));
          qbit_id = {zero_id4};        
          std::shared_ptr<Instruction> physical_inst_14 (std::make_shared <OPinst> (timestep,qec_loc,    "CNOT", qbit_id ));          
          qbit_id = {zero_id4,qec_id[i]};
          std::shared_ptr<Instruction> physical_inst_15 (std::make_shared <OPinst> (timestep,qec_loc,    "H",    qbit_id ));
          qbit_id = {zero_id4};                
          std::shared_ptr<Instruction> physical_inst_16 (std::make_shared <OPinst>(timestep,qec_loc,     "MeasZ",qbit_id ));
          */
          /* X and Z corrections done in software
          qbit_id = {qec_id[i]};        
          std::shared_ptr<Instruction> physical_inst_17 (std::make_shared <OPinst> (timestep,qec_loc,     "X",    qbit_id ));  
          qbit_id = {qec_id[i]};        
          std::shared_ptr<Instruction> physical_inst_18 (std::make_shared <OPinst> (timestep,qec_loc,     "Z",    qbit_id ));            
          */
          
          // Recycling BMOVs
          qbit_id = {zero_id1};                
          std::shared_ptr<Instruction> physical_inst_19 (std::make_shared <BMOVinst> (timestep,qec_loc,  (sub_loc_t)(qec_sub_loc),zero_factories[factory_loc],  G,     qbit_id));
          qbit_id = {zero_id2};                       
          std::shared_ptr<Instruction> physical_inst_20 (std::make_shared <BMOVinst> (timestep,qec_loc,  (sub_loc_t)(qec_sub_loc),zero_factories[factory_loc],  G,     qbit_id));
          /* uncomment for conservative qec
          qbit_id = {zero_id3};  
          std::shared_ptr<Instruction> physical_inst_21 (std::make_shared <BMOVinst> (timestep,qec_loc,  (sub_loc_t)(qec_sub_loc),zero_factories[factory_loc],  G,     qbit_id));          
          qbit_id = {zero_id4};
          std::shared_ptr<Instruction> physical_inst_22 (std::make_shared <BMOVinst> (timestep,qec_loc,  (sub_loc_t)(qec_sub_loc),zero_factories[factory_loc],  G,     qbit_id));          
          */

          // tag sequence numbers, and increment
          physical_inst_1->seq = seq++;
          physical_inst_2->seq = seq++;
          //physical_inst_3->seq = seq++;
          //physical_inst_4->seq = seq++;
          physical_inst_5->seq = seq++;
          physical_inst_6->seq = seq++;
          physical_inst_7->seq = seq++;
          physical_inst_8->seq = seq++;
          physical_inst_9->seq = seq++;
          physical_inst_10->seq = seq++;
          //physical_inst_11->seq = seq++;
          //physical_inst_12->seq = seq++;
          //physical_inst_13->seq = seq++;
          //physical_inst_14->seq = seq++;
          //physical_inst_15->seq = seq++;
          //physical_inst_16->seq = seq++;
          //physical_inst_17->seq = seq++;
          //physical_inst_18->seq = seq++;
          physical_inst_19->seq = seq++;
          physical_inst_20->seq = seq++;
          //physical_inst_21->seq = seq++;
          //physical_inst_22->seq = seq++;

          // add to the final PhysicalInst_v2 vector
          PhysicalInst_v2.push_back(physical_inst_1);
          PhysicalInst_v2.push_back(physical_inst_2);
          //PhysicalInst_v2.push_back(physical_inst_3);
          //PhysicalInst_v2.push_back(physical_inst_4);
          PhysicalInst_v2.push_back(physical_inst_5);
          PhysicalInst_v2.push_back(physical_inst_6);
          PhysicalInst_v2.push_back(physical_inst_7);
          PhysicalInst_v2.push_back(physical_inst_8);
          PhysicalInst_v2.push_back(physical_inst_9);
          PhysicalInst_v2.push_back(physical_inst_10);
          //PhysicalInst_v2.push_back(physical_inst_11);
          //PhysicalInst_v2.push_back(physical_inst_12);  
          //PhysicalInst_v2.push_back(physical_inst_13);
          //PhysicalInst_v2.push_back(physical_inst_14);
          //PhysicalInst_v2.push_back(physical_inst_15); 
          //PhysicalInst_v2.push_back(physical_inst_16);
          //PhysicalInst_v2.push_back(physical_inst_17);
          //PhysicalInst_v2.push_back(physical_inst_18); 
          PhysicalInst_v2.push_back(physical_inst_19);
          PhysicalInst_v2.push_back(physical_inst_20);
          //PhysicalInst_v2.push_back(physical_inst_21); 
          //PhysicalInst_v2.push_back(physical_inst_22);
        }
      }  

    }
    mapOfPhysicalInst_v2[leaf_func] = PhysicalInst_v2;
  }  

  return mapOfPhysicalInst_v2;
}

// Translate each logical teleport to 3 physical moves
// and add teleport operations (CNOT/H/X/Z)
// and add instruction sequence tags
SeqTupleTableTy inject_epr (InstTableTy &mapOfPhysicalInst_v) {

  #ifdef _DEBUG_PROGRESS
    std::cerr<<"injecting EPR..."<<std::endl;
  #endif
  
  // -- Possible Teleportations, and Their Translation --
    // Tel_Mov = G1 ---> T1   /   T1 ----> G1
    // Tel_Mov = T1 ---> T2
    // Tel_Mov = G1 ---> T3
  
  // -- Translation of a Teleportation Move --
  
    // -- Ballistic Moves to Tel. Units --
      // Bal_Mov = src          ---> TU_src  (data)
      // Bal_Mov = T_epr_factory  ---> TU_src  (epr1)
      // Bal_Mov = T_epr_factory  ---> TU_dest (epr2)

    // -- Followed by These Operations --
      // CNOT  (data) (epr1)
      // H     (data)
      // MeasZ (data)
      // MeasZ (epr1)
      // Z     (epr2)         <runtime dependant>
      // X     (epr2)         <runtime dependant>

    // -- Ballistic Moves to Final Locations --
      // Bal_Mov = TU_T1 ---> T1   (epr2)
      // Bal_Mov = TU_G1 ---> G_epr_factory* (epr1)  (*in fact: cooling EEU)
      // Bal_Mov = TU_G1 ---> G_epr_factory* (data)  (*in fact: cooling EEU)

    // -- At This Point --
      // recycle epr1, data 
      // epr2 is data   
  
  SeqTupleTableTy mapOfTeleport_seq_tuples;  
  std::string leaf_func;   
  unsigned long int epr_count;  
  epr_count = 0;    
  for (auto &F : mapOfPhysicalInst_v) {  
    leaf_func = F.first;
    InstVecTy PhysicalInst_v;    
    PhysicalInst_v = F.second;
    InstVecTy PhysicalInst_v2;     // construct this new vector, to avoid lots of inserts
    SeqTupleVecTy teleport_seq_tuples;
    unsigned int seq = 0;          // tag Instruction sequence numbers as we build PhysicalInst_v2
    
    for (auto &I : PhysicalInst_v) {
      if (std::shared_ptr<BMOVinst> ph_mov_inst = std::dynamic_pointer_cast<BMOVinst>(I)) {
        unsigned int timestep = ph_mov_inst->ts;
        unsigned int source = ph_mov_inst->src;
        unsigned int destination = ph_mov_inst->dest;
        std::vector<std::string> q_id = ph_mov_inst->qid;
        sub_loc_t source_sub = ph_mov_inst->src_sub;
        sub_loc_t destination_sub = ph_mov_inst->dest_sub;   

        // if this is a local memory move, or move to/from zero/magic factories, no epr injection is necessary
        if (source_sub == L || destination_sub == L || 
            std::find(zero_factories.begin(),zero_factories.end(),source)!=zero_factories.end() || 
            std::find(zero_factories.begin(),zero_factories.end(),destination)!=zero_factories.end() ||
            std::find(magic_factories.begin(),magic_factories.end(),source)!=magic_factories.end() || 
            std::find(magic_factories.begin(),magic_factories.end(),destination)!=magic_factories.end()) {
          PhysicalInst_v2.push_back(I);
          I->seq = seq++;          
        }     

        // otherwise, check for teleport cases 1 & 2 & 3, and if matched, translate them by injecting eprs       
        // Teleport Cases 1 & 2 & 3
        else if ((source_sub==G && destination_sub==T && source==destination) ||  
                (source_sub==T && destination_sub==G && source==destination) ||
                (source_sub==T && destination_sub==T) ||       
                (source_sub==G && destination_sub==T && source!=destination)) {

          #ifdef _DEBUG_EPR_INJECTION
            std::cout<<"\nInjecting EPRs for..."<<std::endl;
            ph_mov_inst->print();
          #endif            

          // EPR naming scheme: 
          // In order to not create false dependency between epr qubits
          // add a unique "epr_count" number to the end of the name 
          epr_count++;
          int factory_loc = epr_count % epr_factories.size();
          std::string epr_id1 = q_id[0]+"_epr1_"+std::to_string(epr_count);
          std::string epr_id2 = q_id[0]+"_epr2_"+std::to_string(epr_count);  

          // placeholder vector for q_id arguments
          std::vector<std::string> qbit_id (2);        

          // BMOV to Tel. Units
          qbit_id = {q_id[0]};                
          std::shared_ptr<Instruction> physical_inst_1 (std::make_shared <BMOVinst> (timestep,source,                      source_sub,source,     (sub_loc_t)(source_sub+2),     qbit_id));
          qbit_id = {epr_id1};                       
          std::shared_ptr<Instruction> physical_inst_2 (std::make_shared <BMOVinst> (timestep,epr_factories[factory_loc],  T,         source,     (sub_loc_t)(source_sub+2),     qbit_id));
          qbit_id = {epr_id2};        
          std::shared_ptr<Instruction> physical_inst_3 (std::make_shared <BMOVinst> (timestep,epr_factories[factory_loc],  T,         destination,(sub_loc_t)(destination_sub+2),qbit_id));
          
          // Tel. Operations
          qbit_id = {q_id[0],epr_id1};
          std::shared_ptr<Instruction> physical_inst_4 (std::make_shared <OPinst> (timestep,source,     "CNOT", qbit_id ));
          qbit_id = {q_id[0]};        
          std::shared_ptr<Instruction> physical_inst_5 (std::make_shared <OPinst> (timestep,source,     "H",    qbit_id ));
          qbit_id = {q_id[0]};                
          std::shared_ptr<Instruction> physical_inst_6 (std::make_shared <OPinst> (timestep,source,     "MeasZ",qbit_id ));
          qbit_id = {epr_id1};                
          std::shared_ptr<Instruction> physical_inst_7 (std::make_shared <OPinst> (timestep,source,     "MeasZ",qbit_id ));
          /* X and Z corrections done in software
          qbit_id = {epr_id2};                
          std::shared_ptr<Instruction> physical_inst_8 (std::make_shared <OPinst> (timestep,destination,"X",    qbit_id ));
          qbit_id = {epr_id2};                
          std::shared_ptr<Instruction> physical_inst_9 (std::make_shared <OPinst> (timestep,destination,"Z",    qbit_id ));        
          */
        
          // BMOV to final destination
          qbit_id = {epr_id2};                
          std::shared_ptr<Instruction> physical_inst_10 (std::make_shared <BMOVinst> (timestep,destination,(sub_loc_t)(destination_sub+2),destination,           destination_sub,qbit_id)); 
          qbit_id = {epr_id1};                
          std::shared_ptr<Instruction> physical_inst_11 (std::make_shared <BMOVinst> (timestep,source,     (sub_loc_t)(source_sub+2),     epr_factories[factory_loc], G,         qbit_id));
          qbit_id = {q_id[0]};                
          std::shared_ptr<Instruction> physical_inst_12 (std::make_shared <BMOVinst> (timestep,source,     (sub_loc_t)(source_sub+2),     epr_factories[factory_loc], G,         qbit_id));
         
          // tag sequence numbers, and increment
          physical_inst_1->seq = seq++;
          physical_inst_2->seq = seq++;
          physical_inst_3->seq = seq++;
          physical_inst_4->seq = seq++;
          physical_inst_5->seq = seq++;
          physical_inst_6->seq = seq++;
          physical_inst_7->seq = seq++;
          //physical_inst_8->seq = seq++;
          //physical_inst_9->seq = seq++;
          physical_inst_10->seq = seq++;
          physical_inst_11->seq = seq++;
          physical_inst_12->seq = seq++;

          // add to the final PhysicalInst_v2 vector
          PhysicalInst_v2.push_back(physical_inst_1);
          PhysicalInst_v2.push_back(physical_inst_2);
          PhysicalInst_v2.push_back(physical_inst_3);
          PhysicalInst_v2.push_back(physical_inst_4);
          PhysicalInst_v2.push_back(physical_inst_5);
          PhysicalInst_v2.push_back(physical_inst_6);
          PhysicalInst_v2.push_back(physical_inst_7);
          //PhysicalInst_v2.push_back(physical_inst_8);
          //PhysicalInst_v2.push_back(physical_inst_9);
          PhysicalInst_v2.push_back(physical_inst_10);
          PhysicalInst_v2.push_back(physical_inst_11);
          PhysicalInst_v2.push_back(physical_inst_12);    

          // important to know these in order to correctly switch data/epr2 qubits during simulation
          // tuples of <data_mov, epr1_mov, epr2_mov>
          SeqTupleTy T = std::make_tuple(physical_inst_12->seq, physical_inst_11->seq, physical_inst_10->seq);        
          teleport_seq_tuples.push_back(T);
        
        }
        else {
          std::cerr<<"Error: Illegal Teleportation Move!!"<<std::endl;
          exit(1);
        }
      }
      else if (std::shared_ptr<OPinst> ph_op_inst = std::dynamic_pointer_cast<OPinst>(I)) {
        PhysicalInst_v2.push_back(I);
        I->seq = seq++;        
      }   
    } 
    mapOfPhysicalInst_v[leaf_func] = PhysicalInst_v2;    
    mapOfTeleport_seq_tuples[leaf_func] = teleport_seq_tuples;              
  }

  return mapOfTeleport_seq_tuples;
}

// --------------------- Simulator Helper Functions ---------------------

// input: vector of physical instruction pointers
// action: for each leaf, iterate over all its instructions.
// for each qbit, iterate over other insts in order to find an out edge for that inst.
// 1-qbit inst has 1 in edge and 1 out edge. 2-qbit inst has 2 in edges and 2 out edges.
DAGTableTy get_dag(InstTableTy &mapOfPhysicalInst_v) {

  #ifdef _DEBUG_PROGRESS
    std::cerr<<"computing DAG..."<<std::endl;
  #endif
  
  std::vector<std::string> q_id;
  std::vector<std::string> q_id_next;
  std::vector<std::string> q_id_after_teleport;  
  DAGTableTy mapOfDAG;  // return this

  for (auto &F : mapOfPhysicalInst_v) {
    std::string leaf_name = F.first;

    // create a DAG with the size of PhysicalInst_v and populate it
    DAGTy g(F.second.size());   
    for (InstVecTy::const_iterator I = (F.second).begin(); I != (F.second).end(); ++I) {
      bool no_child = true;
      q_id = (*I)->qid;
      
      // for each argument of this Instruction
      for (int i=0; i<q_id.size(); i++) {    
        std::string q_id_i = q_id[i];        
        // iterate through later instructions until 
        // it either finds an inst with the same argument, or runs out of insts      
        auto I_next = std::next(I);
        bool found_next = 0;      
        while (I_next != (F.second).end()) {
          q_id_next = (*I_next)->qid;
          // for each argument of the later instruction
          for (int j=0; j<q_id_next.size(); j++) {
            if (q_id_next[j] == q_id_i) {
              std::pair < EdgeID, bool > p = boost::edge((*I)->seq, (*I_next)->seq, g);
              if ( !p.second ) {  // don't draw duplicate edge 
                boost::add_edge((*I)->seq, (*I_next)->seq, g);
              }
              found_next = 1;   // set this flag. while loop will end.              
            }
            if(found_next) 
              break;            
          }       
          if(found_next) {
            no_child = false;
            break;
          }
          std::advance(I_next,1);
        }          
      }
      (*I)->no_child = no_child;
    }

    // for teleportations, some special dependencies exist. iterate one more time to add/remove these special edges.
    for (InstVecTy::const_iterator I = (F.second).begin(); I != (F.second).end(); ++I) {
      if(std::shared_ptr<OPinst> ph_op_inst = std::dynamic_pointer_cast<OPinst>(*I)) {        
        std::string q_id_i = (ph_op_inst->qid)[0];
        if (ph_op_inst->op_type == "MeasZ" && q_id_i.find("epr1")!=std::string::npos) {
          // the last three of teleport ops (I+1/I+2/I+3) is BMOV of epr2/epr1/data to final locations.
          // they MUST complete first before the teleportation is considered complete, since data and epr switch. So add dep to I_after_teleport.
          std::string q_teleport_data = ((*(I+3))->qid)[0];  // this is the "BMOV data" instruction.
          auto I_after_teleport = I+4;                 // do a similar while loop as above to add edges out of the last 3 ops.
          bool found_after_teleport = 0;            
          while (I_after_teleport != (F.second).end()) {
            q_id_after_teleport = (*I_after_teleport)->qid;
            for (int k=0; k<q_id_after_teleport.size(); k++) {
              if (q_id_after_teleport[k] == q_teleport_data) {
                boost::add_edge((*(I+1))->seq, (*I_after_teleport)->seq, g);
                boost::add_edge((*(I+2))->seq, (*I_after_teleport)->seq, g);                                
                //boost::remove_edge((*(I+3))->seq, (*I_after_teleport)->seq, g);
                // these BMOV instructions now have a child. The epr1/epr2 qubits will be deleted together after the completion of the entire teleportation
                (*(I+1))->no_child = false;
                (*(I+2))->no_child = false;                
                found_after_teleport = 1;
              }
            }
            if(found_after_teleport)
              break;
            I_after_teleport++;
          }
        }
      }
    } 
    mapOfDAG[leaf_name] = g;
  } 
  return mapOfDAG;
}


// input: vector of CG instruction pointers
// action: for each non-leaf, iterate over all its CGinsts.
// for each arg, iterate over other insts in order to find an out edge for that inst.
CallGraphTableTy get_call_graph(const InstTableTy &mapOfCGInst_v) {

  #ifdef _DEBUG_PROGRESS
    std::cerr<<"computing call graph..."<<std::endl;
  #endif

  std::vector<std::string> q_reg;
  std::vector<std::string> q_reg_next;
  CallGraphTableTy mapOfCallGraph;  // return this    
    
  for (auto &F : mapOfCGInst_v) {
    std::string nonleaf_name = F.first;

    // create a DAG with the size of bodyInst_v and populate it
    CallGraphTy g(F.second.size());
    for (InstVecTy::const_iterator I = (F.second).begin(); I != (F.second).end(); ++I) {
      q_reg = (*I)->qid;
      
      // for each argument of this CGInst
      for (int i=0; i<q_reg.size(); i++) {
        std::string q_reg_i = q_reg[i];
        
        // iterate through later CG instructions until
        // it either finds an inst with the same arguments, or runs out of insts
        auto I_next = std::next(I);
        bool found_next = 0;
        while (I_next != (F.second).end()) {
          q_reg_next = (*I_next)->qid;
          // for each argument fo the later instruction
          for (int j=0; j<q_reg_next.size(); j++) {
            if (q_reg_next[j] == q_reg_i) { 
              boost::add_edge((*I)->seq, (*I_next)->seq, g);  // parallel edges won't be added since EdgeList is chosen as setS in boost::adjacency_list
              found_next = 1;   // set this flag, while loop will end
            }
          }
          if (found_next)
            break;
          std::advance(I_next, 1);
        }
      }
    }
    mapOfCallGraph[nonleaf_name] = g;
  }
  return mapOfCallGraph;
}

// Get all leaves that will execute as a result of root (traverse tree with root as root node)
void traverse_call_graph (std::string root_name, const InstTableTy &mapOfCGInst_v) {
  if (mapOfCGInst_v.find(root_name) == mapOfCGInst_v.end()) {  // if root is a leaf module  
    all_ordered_leaf_names.push_back(root_name);
  }
  else {                                                        // if root is a nonleaf
    InstVecTy active_module_body = mapOfCGInst_v.find(root_name)->second;
    call_stack.push(root_name);
    for (auto &I : active_module_body)
      if (std::shared_ptr<CGinst> cg_inst= std::dynamic_pointer_cast<CGinst>(I))
        traverse_call_graph(cg_inst->module_name, mapOfCGInst_v);  
  }
  return;
}

/*
void serialize_dag(InstTableTy &mapOfCGInst_v, InstTableTy &mapOfPhysicalInst_v, SeqTupleTableTy &mapOfTeleport_seq_tuples, DAGTableTy &mapOfDAG, std::string filename) {

  #ifdef _DEBUG_PROGRESS
    std::cerr<<"serializing..."<<std::endl;
  #endif

  // save serialized data
  {
    std::ofstream ofs(filename);
    boost::archive::text_oarchive oa(ofs);
    oa & mapOfCGInst_v;
    oa & mapOfPhysicalInst_v;
    //oa & mapOfTeleport_seq_tuples;
    oa & mapOfDAG;
  }
}

void deserialize_dag(InstTableTy &mapOfCGInst_v, InstTableTy &mapOfPhysicalInst_v, SeqTupleTableTy &mapOfTeleport_seq_tuples, DAGTableTy &mapOfDAG, std::string filename) {

  #ifdef _DEBUG_PROGRESS
    std::cerr<<"deserializing..."<<std::endl;
  #endif

  mapOfCGInst_v.clear();
  mapOfPhysicalInst_v.clear();
  mapOfTeleport_seq_tuples.clear();    
  mapOfDAG.clear();
  
  // load serialized vector
  {
    std::ifstream ifs(filename);
    boost::archive::text_iarchive ia(ifs);
    ia & mapOfCGInst_v;
    ia & mapOfPhysicalInst_v;
    //ia & mapOfTeleport_seq_tuples;    
    ia & mapOfDAG;
  }

  // printout v2 values
  #ifdef _DEBUG_SERIALIZATION   
    print_physical_insts(mapOfPhysicalInst_v);
  #endif   
}*/

// creates IDLE qbits for MOVinsts if they haven't been created before
bool create_instruction_qbits(std::shared_ptr<Instruction> I) {
  std::vector<std::string> q_id = (I)->qid;
  unsigned int loc;
  sub_loc_t sub_loc;
  if (std::shared_ptr<BMOVinst> ph_mov_inst = std::dynamic_pointer_cast<BMOVinst>(I)) {
    loc = ph_mov_inst->src;
    sub_loc = ph_mov_inst->src_sub;       
  } 
  else if (std::shared_ptr<OPinst> ph_op_inst = std::dynamic_pointer_cast<OPinst>(I)) {
    for (int i=0; i<q_id.size(); i++){
      if (q_map.find(q_id[i]) == q_map.end()) {     // all op insts MUST have been created
        std::cerr<<"Error: OPinst qubits not yet created."<<std::endl;     
        return 1;
      }
      else      // already created
        continue;
    }
  }  
  for (int i=0; i<q_id.size(); i++) {
    if (q_map.find(q_id[i]) == q_map.end()) {
      kind_t q_kind;
      if(q_id[i].find("epr1") != std::string::npos) 
        q_kind = EPR1;
      else if (q_id[i].find("epr2") != std::string::npos)
        q_kind = EPR2;
      else if (q_id[i].find("zero") != std::string::npos)
        q_kind = ZERO;
      else if (q_id[i].find("magic") != std::string::npos)
        q_kind = MAGIC;
      else
        q_kind = DATA;  
      #ifdef _DEBUG_READY_QUEUE        
      #ifdef _DEBUG_QUBIT_MANAGEMENT
        std::cout<<"\t\t";
      #endif
      #endif  
      //TODO: this should always be created in EEU. Insert (MOV EEU-->src) instructions in PhysInsts.                     
      qbit* q = new qbit(q_id[i],q_kind,0,loc,sub_loc,IDLE);    
    }
    else      // already created
      continue;    
  }     
  return 0;
}

// copy the data structures of a module. used right before the execution of that module.
void copy_module_data(InstVecTy &PhysicalInst_v_original, InstVecTy &PhysicalInst_v,
                      SeqTupleVecTy &teleport_seq_tuples_original, InstTupleVecTy &teleport_inst_tuples,
                      DAGTy &g_original, DAGTy &g) {
    // create a working set of PhysicalInst_v, by copying the Instruction objects (get new pointers to all)
    std::shared_ptr<Instruction> I_copy;
    for (auto &I_original : PhysicalInst_v_original) {  
      if (std::shared_ptr<BMOVinst> ph_mov_inst_original = std::dynamic_pointer_cast<BMOVinst>(I_original))
        I_copy = std::make_shared <BMOVinst> (*ph_mov_inst_original);
      else if (std::shared_ptr<OPinst> ph_op_inst_original = std::dynamic_pointer_cast<OPinst>(I_original))
        I_copy = std::make_shared <OPinst> (*ph_op_inst_original);
      PhysicalInst_v.push_back(I_copy);
    }
    I_copy = NULL;  // release the temporary pointer    

    // create a working teleport_inst_tuples, from the teleport_seq_tuples.
    // this should contain the same pointers as PhysicalInst_v instructions. construct it from the tuples of seq numbers found before.
    for (auto &t : teleport_seq_tuples_original) {  
      unsigned int one = std::get<0>(t);
      unsigned int two = std::get<1>(t);
      unsigned int three = std::get<2>(t);      
      InstTupleTy T = std::make_tuple(PhysicalInst_v[one], PhysicalInst_v[two], PhysicalInst_v[three]);
      teleport_inst_tuples.push_back(T);
    }

    // just copy because DAG contains pure unsigned int, not vectors
    g = g_original; 
}

// ------------------------ Simulator Functions --------------------------
// input: constant ref. to vector of all physical instructions and vector of all currently executing instructions 
// action: identify the physical instructions whose qbit(s) have no prior dependency (first in chain)
// no action: does NOT delete any inst from ready_queue. Those will be deleted UPON COMPLETION by the simulator.
// important: this must preserve order of instructions in the ready_queue to be similar to that of PhysicalInst_v, otherwise smoothing won't happen.
// bool is_next: signifies whether this is for a current module or prefetching the next module. If next module, only EPR movs should be considered.
InstVecTy initialize_ready_queue(InstVecTy &PhysicalInst_v, DAGTy &g, bool is_next) {
  
  InstVecTy ready_queue;
  for (auto const &I : PhysicalInst_v){
    
    #ifdef _DEBUG_READY_QUEUE   
    //  std::cout<<"\nevaluating..."<<std::endl<<"\t";
    //  I->print();
    #endif

    unsigned int seq = I->seq;  

    // no in edge means not waiting for a previous instruction: means "ready"      
    if (boost::in_degree(seq, g)==0) {   

      // BMOVinst: error catching
      if (std::shared_ptr<BMOVinst> ph_mov_inst = std::dynamic_pointer_cast<BMOVinst>(I)) {
        std::string q = (I->qid)[0];        
        if (is_next && 
            (q.find("zero") == std::string::npos) && 
            (q.find("epr") == std::string::npos) && 
            (q.find("magic") == std::string::npos)) // not interested if this is not an ancilla mov being prefetched from the next module
          continue;
        if (q_map.find(q) != q_map.end()) {   // only check this for the qbits already created, others are automatically ready.
          // FIXME: between modules, a qubit may persist. There should be an instruction taking it from end-location of 1stMod to start-location of 2ndMod.
          
          //int loc = ph_mov_inst->src;
          //sub_loc_t sub_loc = ph_mov_inst->src_sub;   
          //int q_loc = q_map[q]->loc;    
          //sub_loc_t q_sub_loc = q_map[q]->sub_loc;
          //if ((q_loc != loc) || (q_sub_loc != sub_loc)) { //FIXME:RESTORE THIS
          //  std::cerr<<"Error: BMOV inst qbit "<<q<<" not in current loc and subloc."<<std::endl;
          //  exit(1);
          //}
          q_map[q]->loc = ph_mov_inst->src; // FIXME:REMOVE THIS
          q_map[q]->sub_loc = ph_mov_inst->src_sub; // FIXME:REMOVE THIS          
        }
      } 
            
      // OPinst: error catching
      else if (std::shared_ptr<OPinst> ph_op_inst = std::dynamic_pointer_cast<OPinst>(I)) {
        if (is_next)  // not interested if this is an OP being prefetched from the next module
          continue;
        bool q_not_ready = 0;
        for (int i=0; i<(I->qid).size(); i++) {
          std::string q = (I->qid)[i];
          if (q_map.find(q) == q_map.end()) {
            std::cerr<<"Error: OPinst qbit "<<q<<" not yet created, but all edge dependencies met."<<std::endl;
            exit(1);
          }
          int loc = ph_op_inst->zone;  
          int q_loc = q_map[q]->loc;    
          sub_loc_t q_sub_loc = q_map[q]->sub_loc;
          if ((q_loc != loc)) {   // TODO: op sub_loc can right now be T, as well as TU_T and TU_G. change?          
            q_not_ready = 1;
          }
        }
        if (q_not_ready) {
          std::cerr<<"Error: BMOV inst qbit not in current loc and subloc."<<std::endl;
          exit(1);
        }
      }     
      
      ready_queue.push_back(I);        
      #ifdef _DEBUG_READY_QUEUE    
      //  std::cout<<"\t\tinstruction is ready!"<<std::endl;
      #endif  
      
      // qbits will be created just before execution, not here.
    }
  }

  return ready_queue;
}

// input: ref. to ready_queue to be filled in place, and the seq number of the Instruction to be insrerted
// action: insert the Instruction with ->seq=seq in the correct space in ready_queue (ready_queue order must always preserve PhysicalInst_v order)
void insert_ready_queue(InstVecTy &PhysicalInst_v, InstVecTy &ready_queue, unsigned int seq) {
  // use "seq" field of Instruction object to recognize the place of insertion.  
  // can't find the instruction by doing PhysicalInst_v[seq] because Instructions are continuously removed from PhysicalInst_v.
  std::shared_ptr<Instruction> I = which_Instruction(PhysicalInst_v, seq);
  auto place = std::find_if(ready_queue.begin(), ready_queue.end(), [&seq](const std::shared_ptr<Instruction> I_after){return (I_after->seq > (unsigned int)seq);});
  ready_queue.insert(place, I);
}

// input: pointer to a recently completed instruction ce_I, and the dag of operations g.
// action: update the dag by removing dependency edges that no longer exist.
void update_completed_instruction (InstVecTy &PhysicalInst_v, InstVecTy &ready_queue, DAGTy &g, std::shared_ptr<Instruction> &ce_I) {
  ce_I->is_complete = true;  
  ce_I->is_executing = false;  
  
  // for the instructions connected by out edges, find and delete the in edge from this instruction
  unsigned int seq = ce_I->seq;
  #ifdef _DEBUG_SIMULATOR
    std::cout<<"#"<<seq<<" completed!"<<std::endl;
  #endif
  DAGTy::adjacency_iterator neighborIt, neighborEnd;
  boost::tie(neighborIt, neighborEnd) = adjacent_vertices(seq, g);
  for (; neighborIt != neighborEnd; ++neighborIt) {
    VertexID seq_out = *neighborIt;
    #ifdef _DEBUG_SIMULATOR
      std::cout<<"removing dep. to #"<<seq_out<<std::endl;        
    #endif
    remove_edge(seq, seq_out, g);// FIXME: this is WRONG! invalidates both neighborIt and neighborEnd. put it out of loop.
    if(boost::in_degree(seq_out, g) == 0)
      insert_ready_queue(PhysicalInst_v, ready_queue, (unsigned int)seq_out);  
  }   
  clear_out_edges(seq, g);             // no more out edges for this instruction either
}

// input: exec_queue and exec_queue_n, containing instructions executing from current and next modules.
// action: update the status of qbits pertainig to each instruction
// 1. one-hop change in qbit position for BMOVinst 2. decrement qbit op remaining time for OPinst
// action: if an instruction completes, update edges of dag and add next ones to the ready queue
void do_one_cycle(InstVecTy &PhysicalInst_v, InstVecTy &ready_queue, InstVecTy &exec_queue, DAGTy &g, 
                  InstVecTy &exec_queue_n, std::vector<unsigned int> &previously_completed) {

  // exec_queue of the current module 
  for (auto &ce_I : exec_queue) {
    std::vector<std::string> q_id = ce_I->qid;              
    for (int i=0; i<q_id.size(); i++) {  
      qbit *q = (q_map.find(q_id[i]))->second;
      // update qbit state based on BMOVinst content or OPinst content
      if (q->state == IN_MOV) {
        // TODO: Route qbit one hop ahead
        q->random--;
        if (q->random == 0) {
          q->loc = q->destination;  //TODO: remove after fixing this part of code
          q->sub_loc = q->destination_sub;
        }
        
        if (q->loc==q->destination && q->sub_loc==q->destination_sub)
          q->state = IDLE;     

        #ifdef _DEBUG_SIMULATOR
          std::cout<<"Simulating move (from current module) on: "<<std::endl;
          q->status();
        #endif        
      }
      else if (q->state == IN_OP) {        
        q->op_time_remaining--; 
        if (q->op_time_remaining == 0)
          q->state = IDLE;
        
        #ifdef _DEBUG_SIMULATOR
          std::cout<<"Simulated operation on: "<<std::endl;
          q->status();
        #endif        
      } 
      else {
        std::cerr<<"Error: Executing qubit (current) "<<q_id[i]<<" is neither IN_MOV nor IN_OP."<<std::endl;
        exit(1);
      }

      // if the qbit is IDLE, then the BMOVinst/OPinst has finished
      // this flow holds for two-qbit gates, since they are acted on simultaneously
      if (q->state == IDLE)
        update_completed_instruction (PhysicalInst_v, ready_queue, g, ce_I);    
    }
  }
  
  // exec_queue of the next module
  for (auto &ce_I : exec_queue_n) {
    std::vector<std::string> q_id = ce_I->qid;              
    for (int i=0; i<q_id.size(); i++) {  
      qbit *q = (q_map.find(q_id[i]))->second;
      // update qbit state based on BMOVinst content or OPinst content
      if (q->state == IN_MOV) {
        // TODO: Route qbit one hop ahead
        q->random--;
        if (q->random == 0) {
          q->loc = q->destination;  //TODO: remove after fixing this part of code
          q->sub_loc = q->destination_sub;
        }
        
        if (q->loc==q->destination && q->sub_loc==q->destination_sub)
          q->state = IDLE; 

        #ifdef _DEBUG_SIMULATOR
          std::cout<<"Simulated move (from next module) on: "<<std::endl;      
          q->status();
        #endif        
      }
      else if (q->state == IN_OP) {    
        std::cerr<<"Error: OP instruction has been prefetched from the next module."<<std::endl;
        exit(1);
      }  
      else {
        std::cerr<<"Error: Executing qubit (next) "<<q_id[i]<<" is neither IN_MOV nor IN_OP."<<std::endl;
        exit(1);        
      }  

      //if the qbit is IDLE, then the BMOVinst has finished. insert its seq so next module can read it.
      if (q->state == IDLE) {
        previously_completed.push_back(ce_I->seq);
        ce_I->is_complete = true;
        ce_I->is_executing = false;
      }
    }
  }

  // increment the age of all existing qbits
  for (auto &it : q_map)
    (it.second)->age++;
}

// ------------------------------- main -----------------------------------
int main (int argc, char *argv[]) {
  
  for (int i = 0; i<argc; i++) {
    if (strcmp(argv[i],"--p")==0) {
      if (argc > (i+1)) {
        P_error_rate = atoi(argv[i+1]);
      }
      else {
        std::cout<<"Usage: $ router <benchmark> --p <P_error_rate>";
        return 1;
      }
    }        
    if (strcmp(argv[i],"--cap")==0) {
      if (argc > (i+1)) {
        if (strcmp(argv[i+1],"inf")!=0) {
          cap = true;
          forward = true;
          mov_cap = atoi(argv[i+1]);
        }
      }
      else {
        std::cout<<"Usage: $ router <benchmark> --cap <move_cap>";
        return 1;
      }
    }
    if (strcmp(argv[i],"--window")==0) {
      if (argc > (i+1)) {
        if (strcmp(argv[i+1],"inf")!=0) {
          window = true;
          forward = true;
          window_size = atoi(argv[i+1]);
        }
      }
      else {
        std::cout<<"Usage: $ router <benchmark> --window <window_size>";
        return 1;
      }
    }    
    if (strcmp(argv[i],"--forward")==0) {
      backward = false;
      forward = true;
    }      
    if (strcmp(argv[i],"--backward")==0) {
      backward = true;
      forward = false;
    }       
    if (strcmp(argv[i],"--backforth")==0) {
      backward = true;
      forward = true;
    }         
    if (strcmp(argv[i],"--usage")==0) {
      report_usage = true;
    }
    if (strcmp(argv[i],"--ages")==0) {
      report_ages = true;
    } 
    if (strcmp(argv[i],"--storage")==0) {
      report_storage = true;
    }        
  }

  // -------- Start processing benchmark
  // -------- Input: Logical LPFS schedule for leaves
  // -------- Input: Logical coarse-grain list schedule for modules
  std::string benchmark_path(argv[1]); 
  std::string benchmark_dir = benchmark_path.substr(0, benchmark_path.find_last_of('/'));
  std::string benchmark_name = benchmark_path.substr(benchmark_path.find_last_of('/')+1, benchmark_path.length());
  std::string LPFSPath = benchmark_path+".lpfs";
  std::string FreqPath = benchmark_path+".freq";
  std::string CGPath = benchmark_path+".cg";  
  #ifdef _DEBUG_PROGRESS
    std::cerr<<"---------------- Processing ---------------"<<std::endl
             << benchmark_name<<std::endl
             <<"-------------------------------------------"<<std::endl;
  #endif

  // -------- Major data-structures
  InstTableTy mapOfLogicalInst_v;
  std::unordered_map<std::string, unsigned long long> mapOfFreq;
  InstTableTy mapOfCGInst_v;
  InstTableTy mapOfPhysicalInst_v;  
  SeqTupleTableTy mapOfTeleport_seq_tuples;
  DAGTableTy mapOfDAG;  
  CallGraphTableTy mapOfCallGraph;  
    
  // -------- Using boost::serialization library, major data-structures will be serialized to file
  // -------- If that file exists already read from it. If not, write the file for the first time.
  std::string serial_file_name = benchmark_path+".serial";  
  std::ifstream serial_file(serial_file_name);
  //if (serial_file.good())
  //  deserialize_dag(mapOfCGInst_v, mapOfPhysicalInst_v, mapOfTeleport_seq_tuples, mapOfDAG, serial_file_name);
    
  //else {
    // -------- From: LPFS Leaves schedule
    // -------- Parse: SIMD-K, MOVinst(ts, qbit, src, dest), OPinst(ts, optype, opzone, qbits)
    // -------- Return: Map between leaf function name and vector of logical Instruction* for that leaf
    mapOfLogicalInst_v = parse_LPFS_file(LPFSPath);
    #ifdef _DEBUG_LOGICAL_INSTS
      print_logical_insts(mapOfLogicalInst_v);
    #endif

    // -------- From: Runtime frequency estimation
    // -------- Parse: leaf modules
    // -------- Return: Map between leaf function name and aggregate number of occurances of that leaf    
    mapOfFreq = parse_freq(FreqPath);

    // -------- From: Coarse-Grain schedule
    // -------- Parse: modules
    // -------- Return: vector of Instruction* (all of them CGinsts*)
    mapOfCGInst_v = parse_CG_file(CGPath);
    #ifdef _DEBUG_CG_INSTS
      print_cg_insts(mapOfCGInst_v);
    #endif

    // -------- Calculate concatenation level
    // -------- Find number of Logical OPinst in each leaf module
    // -------- Multiply by the frequency of that leaf module (mapOfFreq)
    if (P_error_rate < P_th) {
      std::cerr << "Physical error rate is higher than the code threshold. Terminating..\n";
      exit(1);
    }  
    
    total_logical_gates = 0;
    for (auto &m_it : mapOfLogicalInst_v) {
      std::string leaf_name = m_it.first;         
      unsigned int leaf_size = 0;
      for (auto &v_it : m_it.second)
        if (std::shared_ptr<OPinst> op_inst = std::dynamic_pointer_cast<OPinst>(v_it)) 
          leaf_size++;
      unsigned long long leaf_freq = mapOfFreq[leaf_name];
    #ifdef _DEBUG_PROGRESS
      std::cerr<<"leaf: "<<leaf_name<<" - size: "<<leaf_size<<" - freq: "<<leaf_freq<<std::endl;
    #endif    
      total_logical_gates += (unsigned long long)leaf_size * leaf_freq;
    }

    std::cerr << "total logical gates: " << total_logical_gates << std::endl;
    
    L_error_rate = (double)epsilon/(double)total_logical_gates; 
    double c = 1.0 / pow(10.0,(-1.0*(double)P_th));
    double two_to_l = ( log(c * L_error_rate) / log(c * pow(10.0,(-1.0*(double)P_error_rate))) );
    concatenation_level = (int) ceil ( log2(two_to_l) );
    if ( L_error_rate > pow(10.0,(-1.0*(double)P_error_rate)) ) concatenation_level = 0; // very small circuit (large L_error_rate), means no concatenation
    #ifdef _DEBUG_PROGRESS
      std::cerr<<"Physical error rate (p): " << P_error_rate << std::endl;  
      std::cerr<<"Logical error rate (p_L): " << L_error_rate << std::endl;  
      std::cerr<<"Concatenation level (l): "<<concatenation_level<<std::endl;
    #endif
      
    if (concatenation_level < 0) {
      std::cerr << "Inoperable parameters for concatenated code operation. Try changing physical error rate.\n";
      exit(1);
    }  

    // -------- Set the latency and size of ancilla factories
    create_factories(concatenation_level);

    // -------- Create mesh based on number of comuptation regions and number of magic/zero/epr factories
    // -------- Naive model: put ancilla factories on the sides (K+1, K+2, ...)
    mapOfLogicalInst_v = create_mesh(mapOfLogicalInst_v);
    #ifdef _DEBUG_MESH 
      print_mesh();
    #endif  

    // -------- Create Physical Instructions from Logical Instructions
    // -------- Determine Global memory map, add QEC redundancies, add EPRs
    // -------- Return: Map between leaf function name and vector of physical Instruction* for that leaf
    mapOfPhysicalInst_v = annotate_global_memory(mapOfLogicalInst_v);
    mapOfPhysicalInst_v = add_qec(mapOfPhysicalInst_v);
    mapOfTeleport_seq_tuples = inject_epr(mapOfPhysicalInst_v);

    #ifdef _DEBUG_PHYSICAL_INSTS  
      print_physical_insts(mapOfPhysicalInst_v);
    #endif  
      
    // -------- Create DAG dependency graph for all leaf modules (basic blocks)
    // -------- Construct a boost::adjacancy_list graph for each module
    mapOfDAG = get_dag(mapOfPhysicalInst_v);
    #ifdef _DEBUG_DAG
      print_dag(mapOfDAG);
    #endif

    // -------- Create DAG dependency graph for the coarse-grain modules, rooted at main
    // -------- Construct a boost::adjacancy_list graph for each module
    mapOfCallGraph = get_call_graph(mapOfCGInst_v);    
    #ifdef _DEBUG_CALL_GRAPH
      print_call_graph(mapOfCallGraph, mapOfCGInst_v);
    #endif

  //  serialize_dag(mapOfCGInst_v, mapOfPhysicalInst_v, mapOfTeleport_seq_tuples, mapOfDAG, serial_file_name);
  //}


  // -------- Event-driven simulator
  // -------- keep a ready queue  
  // -------- for all instructions that are ready keep simulating cycle-by-cycle
  // -------- keep an eye on completions; must re-evaluate ready queue
  #ifdef _DEBUG_PROGRESS
    std::cerr<<"starting simulation..."<<std::endl;
  #endif

  // Traverse the leaves in pre-order depth-first, starting at main.
  // Get list of all leaves in order of execution. Knowing this statically helps in backward smoothing later on.    
  std::string root_module_name = "main";
  traverse_call_graph(root_module_name, mapOfCGInst_v); 


  // Essential data-structures of a simulating leaf module
  InstVecTy ready_queue;    // which ones have all dependencies met and are ready to go? 
  InstVecTy ready_queue_n;  // which ones are ready to go from the next module? (just calculated once in current module, not updated) 
  InstVecTy current_exec;   // which ones are currently being executed? 
  InstVecTy current_exec_n; // which one are currently being executed from the next module?
  std::vector<unsigned int> previously_completed; // which ones belonging to the next module are completed in this module?
 
  // The original, non-executed data-structures of a particular leaf module: Its PhysicaInstruction_v, teleport_seq_tuples, and DAG
  InstVecTy PhysicalInst_v_original;
  InstVecTy PhysicalInst_v_original_n;  
  SeqTupleVecTy teleport_seq_tuples_original;
  SeqTupleVecTy teleport_seq_tuples_original_n;  
  DAGTy g_original;
  DAGTy g_original_n;
  
  #ifdef _DEBUG_PERFORMANCE
    unsigned int base_time = clock();
    unsigned int agg_prev_time = 0;  
  #endif

  // find a small subset of all_ordered_leaf_names, to speed up simulation
  std::map<std::string, int> sub_ordered_leaf_count;
  sub_ordered_leaf_count.clear();  
  for (auto &i : all_ordered_leaf_names) {
    if (sub_ordered_leaf_count.find(i)==sub_ordered_leaf_count.end()) {
      sub_ordered_leaf_count[i] = 1;
      sub_ordered_leaf_names.push_back(i);
    }
    else if (sub_ordered_leaf_count[i] < LEAF_SIMULATION_MAX) {
      sub_ordered_leaf_count[i]++;
      sub_ordered_leaf_names.push_back(i);
    }      
    else
      continue;
  }

  // -------- Traverse tree of coarse-grain instructions
  // -------- Find leaf modules to simulate 
  for (std::vector<std::string>::const_iterator F = sub_ordered_leaf_names.begin(); F != sub_ordered_leaf_names.end(); F++) {
    std::string leaf_name = (*F);    
 
    #ifdef _DEBUG_PROGRESS
      std::cerr<<"---------- Leaf Module: "<<leaf_name<<" ----------"<<std::endl;
    #endif

    // copy the originals into a working set, so the original stays intact 
    // because PhysicalInsts, teleport_inst_tuples and DAG edges get deleted in simulation
    PhysicalInst_v_original = mapOfPhysicalInst_v.find(leaf_name)->second;
    teleport_seq_tuples_original = mapOfTeleport_seq_tuples.find(leaf_name)->second;
    g_original = mapOfDAG.find(leaf_name)->second;

    InstVecTy PhysicalInst_v; 
    InstTupleVecTy teleport_inst_tuples;
    DAGTy g;

    copy_module_data (PhysicalInst_v_original, PhysicalInst_v, teleport_seq_tuples_original, teleport_inst_tuples, g_original, g);

    // maybe some of this module's Instructions were completed during the simulation of the last module?
    // If so read them and update g. Delete them from PhysicalInst_v
    for (auto &s : previously_completed) {
      std::shared_ptr<Instruction> I = which_Instruction(PhysicalInst_v, s);
      #ifdef _DEBUG_BACK_SMOOTH
        std::cout<<"previously completed: "<<I->seq<<std::endl;      
      #endif      
      update_completed_instruction (PhysicalInst_v, ready_queue, g, I);          
      PhysicalInst_v.erase(std::remove_if(PhysicalInst_v.begin(),PhysicalInst_v.end(),
                          [&s](const std::shared_ptr<Instruction> &I) { return (I->seq == s); }),
                          PhysicalInst_v.end() );      
    }
    previously_completed.clear();   // used this information. reset it to record the prefetches of next module.

    // mark the first set of ready instructions (those with no in_edge dependence)
    ready_queue = initialize_ready_queue(PhysicalInst_v, g, false);

    // get ready queue of the next module. pre-fetch ancilla moves from this ready_queue_n in case of backward smoothing.
    if (backward && *F!=sub_ordered_leaf_names.back()) {
      std::string leaf_name_n = (*(F+1));   
      if (leaf_name != leaf_name_n) {      // don't cross smooth two similar modules: qubits can go wrong.

      PhysicalInst_v_original_n = mapOfPhysicalInst_v.find(leaf_name_n)->second; 
      teleport_seq_tuples_original_n = mapOfTeleport_seq_tuples.find(leaf_name_n)->second;
      g_original_n = mapOfDAG.find(leaf_name_n)->second;

      InstVecTy PhysicalInst_v_n; 
      InstTupleVecTy teleport_inst_tuples_n;
      DAGTy g_n;
      
      copy_module_data (PhysicalInst_v_original_n, PhysicalInst_v_n, teleport_seq_tuples_original_n, teleport_inst_tuples_n, g_original_n, g_n);
            
      // for backward smoothing: mark the ancilla movs of the next module.
      ready_queue_n = initialize_ready_queue(PhysicalInst_v_n, g_n, true); 
      }     
    }
    
    // simulate until all PhysicalInsts are executed. 
    // don't stop if some prefetching of ancilla from the next module is in the midst of execution.
    while (!PhysicalInst_v.empty() || !current_exec_n.empty()) { 
      #ifdef _DEBUG_SIMULATOR    
        std::cout<<"\ncycle: "<<cycle<<std::endl;    
      #endif

      #ifdef _DEBUG_PERFORMANCE
        if(cycle%10000==0) {
          unsigned long long agg_cycle_time = (clock()-base_time)/(CLOCKS_PER_SEC);
          std::cerr<<cycle<<" "<<agg_cycle_time-agg_prev_time<<std::endl;
          agg_prev_time = agg_cycle_time;
        }
      #endif
       
  
      // Greedy Policy: make all "ready" instructions "executing"
      // Cap Policy: no more than mov_cap instructions enter "executing"
      // Window Policy: only those whose out_edge is within window_size of head of PhysicalInst_v enter "executing"
      unsigned int mov_count = 0;
      bool mov_cap_reached = false; 
      unsigned int back_mov_count = 0;      
      bool back_mov_cap_reached = false;

      // -------- Issue to executing queue + (Forward Smoothing)
      for (InstVecTy::iterator rq_I=ready_queue.begin(); rq_I!=ready_queue.end(); ++rq_I) {
        // BMOVinst
        std::string q = ((*rq_I)->qid)[0];                
        if (std::shared_ptr<BMOVinst> ph_mov_inst = std::dynamic_pointer_cast<BMOVinst>(*rq_I)) {
          // in case of forward smoothing, think about when to issue those ancilla movs whose qbits haven't been created.
          // reason: we want to route those JIT. With ancillas that are created, we want to be greedy and finish them soon.
          if(forward && q_map.find(q)==q_map.end() && (q.find("zero")!=std::string::npos || q.find("epr")!=std::string::npos || q.find("epr")!=std::string::npos ) 
              ) {
            if(cap) {
              if (mov_cap_reached)   // skip this one. too many BMOVinsts currently happening.
                continue;
              else 
                mov_count++;
              if (mov_count >= mov_cap)
                mov_cap_reached = true;
            }
          
            if (window) {
              bool upcoming_out_edge = false;
              unsigned int checked_insts = 0;
              unsigned int seq = (*rq_I)->seq;

              if (boost::in_degree(seq, g) != 0 || boost::out_degree(seq, g) != 1)
                std::cerr<<"Error: Window method checking EPR move that either doesn't contain 0 in_edge or 1 out_edge."<<std::endl;

              DAGTy::adjacency_iterator neighborIt, neighborEnd;
              boost::tie(neighborIt, neighborEnd) = adjacent_vertices(seq, g);
              VertexID seq_out = *neighborIt;

              for (InstVecTy::iterator I=PhysicalInst_v.begin(); I!=PhysicalInst_v.end(); ++I) {
                checked_insts++;
                if (checked_insts >= window_size)
                  break;
                if ((*I)->seq == (unsigned int)seq_out) {
                  upcoming_out_edge = true;
                  break;
                }
              }
              if (!upcoming_out_edge)   // skip this one. too early to create and send this ancilla qubit.
                continue;
            }
          }

          // BMOVinst is issued: make sure the qbit is in fact currently in the source
          if (q_map.find(q) != q_map.end()) {   // only check for qbits already created, some may be created with execution
            int loc = ph_mov_inst->src;
            sub_loc_t sub_loc = ph_mov_inst->src_sub;   
            int q_loc = q_map[q]->loc;    
            sub_loc_t q_sub_loc = q_map[q]->sub_loc;
            if ((q_loc != loc) || (q_sub_loc != sub_loc)) {
              std::cerr<<"Error: MOVinst qbit "<<q<<" not in correct loc."<<std::endl;
              return 1;
            }
          }
        } 
              
        // OPinst is always issued: make sure the qbit is in fact currently in the op location
        else if (std::shared_ptr<OPinst> ph_op_inst = std::dynamic_pointer_cast<OPinst>(*rq_I)) {
          for (int i=0; i<((*rq_I)->qid).size(); i++) {
            std::string q = ((*rq_I)->qid)[i];
            if (q_map.find(q) == q_map.end()) {
              std::cerr<<"Error: OPinst qbit not yet created, but marked as ready."<<std::endl;
              return 1;
            }
            int loc = ph_op_inst->zone;  
            int q_loc = q_map[q]->loc;    
            if (q_loc != loc) {   // TODO: sub_loc can right now be T, as well as TU_T and TU_G. change?          
              std::cerr<<"Error: OPinst qbit "<<q<<" not in correct loc."<<std::endl;
              return 1;
            }
          }
        }     
        
        std::shared_ptr<Instruction> ce_I = *rq_I;
        current_exec.push_back(ce_I);
        ce_I->is_executing = true;

        // create qbits just before execution
        // set attributes in all qbits in all instructions that will execute          
        create_instruction_qbits(ce_I);      
        std::vector<std::string> q_id = ce_I->qid;              
        if (std::shared_ptr<BMOVinst> ph_mov_inst = std::dynamic_pointer_cast<BMOVinst>(ce_I)) {
          unsigned int source = ph_mov_inst->src;
          sub_loc_t source_sub = ph_mov_inst->src_sub;           
          unsigned int destination = ph_mov_inst->dest;
          sub_loc_t destination_sub = ph_mov_inst->dest_sub;      
          for (int i=0; i<q_id.size(); i++) {  
            qbit *q = (q_map.find(q_id[i]))->second;
            q->state = IN_MOV;
            q->destination = destination; 
            q->destination_sub = destination_sub;
            if ( source_sub==L || destination_sub==L || 
                (source==destination && (source_sub==TU_G||source_sub==TU_T||destination_sub==TU_G||destination_sub==TU_T)) )  // these BMOVs don't use the mesh
              q->random = 1;          
            else {                                                                                                             // these BMOVs do use the mesh
              int src_row = 2 * (int)((source-1) / (SIMD_cols)) + 1;
              int src_col = (source-1) % (SIMD_cols) + 1;
              int dest_row = 2 * (int)((destination-1) / (SIMD_cols)) + 1;
              int dest_col = (destination-1) % (SIMD_cols) + 1;         
              if ( (src_col%2==0&&(source_sub==T||source_sub==TU_T)) || (src_col%2==1&&(source_sub==G||source_sub==TU_G)) ) src_row++;                          // adjust for subloc
              if ( (dest_col%2==0&&(destination_sub==T||destination_sub==TU_T)) || (dest_col%2==1&&(destination_sub==G||destination_sub==TU_G)) ) dest_row++;   // adjust for subloc 
              int route_distance = abs(src_row - dest_row) + abs(src_col - dest_col);
              q->random = rand() % ( route_distance * (int)pow(7.0,(double)concatenation_level) ) + 1;                  
            }
          }
        }
        else if (std::shared_ptr<OPinst> ph_op_inst= std::dynamic_pointer_cast<OPinst>(ce_I)) {
          std::string op_type = ph_op_inst->op_type;      
          for (int i=0; i<q_id.size(); i++) {
            qbit *q = (q_map.find(q_id[i]))->second;
            q->state = IN_OP;
            q->op_time_remaining = op_delays.find(op_type)->second;   
          }
        }      
      }

      // erase "is_executing" insts from "ready_queue": use erase-remove idiom
      ready_queue.erase( std::remove_if(ready_queue.begin(), ready_queue.end(), 
                            [](const std::shared_ptr<Instruction> I) { return I->is_executing; }), 
                            ready_queue.end() );
      if (!forward && !ready_queue.empty()) {
        std::cerr<<"Error: greedy policy failed to empty the ready_queue."<<std::endl;
        return 1;
      }

      // -------- (Backward Smoothing)
      if (backward) {
        // only look at BMOVInst for ancilla qubits
        for (InstVecTy::iterator rq_I=ready_queue_n.begin(); rq_I!=ready_queue_n.end(); ++rq_I) {
          // BMOVinst: make sure the qbit is in fact currently in the source
          if (std::shared_ptr<BMOVinst> ph_mov_inst = std::dynamic_pointer_cast<BMOVinst>(*rq_I)) {
            if (back_mov_cap_reached)   // skip this one. too many BMOVinsts currently happening.
              break;                    // no need to continue because there aren't any OPinsts here. All are MOVinsts.
            else 
              back_mov_count++;
            if (back_mov_count >= mov_cap)
              back_mov_cap_reached = true;      
            
            std::shared_ptr<Instruction> ce_I = *rq_I;
            current_exec_n.push_back(ce_I);
            ce_I->is_executing = true;

            // create qbits just before execution
            create_instruction_qbits(ce_I);      
            std::vector<std::string> q_id = ce_I->qid;              
            if (std::shared_ptr<BMOVinst> ph_mov_inst = std::dynamic_pointer_cast<BMOVinst>(ce_I)) {
              unsigned int source = ph_mov_inst->src;
              sub_loc_t source_sub = ph_mov_inst->src_sub;                  
              unsigned int destination = ph_mov_inst->dest;
              sub_loc_t destination_sub = ph_mov_inst->dest_sub;      
              for (int i=0; i<q_id.size(); i++) {  
                qbit *q = (q_map.find(q_id[i]))->second;
                q->state = IN_MOV;
                q->destination = destination; 
                q->destination_sub = destination_sub; 
                if ( source_sub==L || destination_sub==L || 
                   (source==destination && (source_sub==TU_G||source_sub==TU_T||destination_sub==TU_G||destination_sub==TU_T)) )  // these BMOVs don't use the mesh
                  q->random = 1;          
                else {                                                                                                            // these BMOVs do use the mesh
                  int src_row = 2 * (int)((source-1) / (SIMD_cols)) + 1;
                  int src_col = (source-1) % (SIMD_cols) + 1;
                  int dest_row = 2 * (int)((destination-1) / (SIMD_cols)) + 1;
                  int dest_col = (destination-1) % (SIMD_cols) + 1;         
                  if ( (src_col%2==0&&(source_sub==T||source_sub==TU_T)) || (src_col%2==1&&(source_sub==G||source_sub==TU_G)) ) src_row++;                          // adjust for subloc
                  if ( (dest_col%2==0&&(destination_sub==T||destination_sub==TU_T)) || (dest_col%2==1&&(destination_sub==G||destination_sub==TU_G)) ) dest_row++;   // adjust for subloc     
                  int route_distance = abs(src_row - dest_row) + abs(src_col - dest_col);
                  q->random = rand() % ( route_distance * (int)pow(7.0,(double)concatenation_level) ) + 1;                  
                }
              }
            }          
          } 
          else {
            std::cout<<"Error: Trying to execute a non-BMOVinst from the next module."<<std::endl;
            exit(1);
          }
        } 

        // erase "is_executing" insts from "ready_queue_n": use erase-remove idiom
        ready_queue_n.erase( std::remove_if(ready_queue_n.begin(), ready_queue_n.end(), 
                              [](const std::shared_ptr<Instruction> I) { return I->is_executing; }), 
                              ready_queue_n.end() );
      }



      // advance simulation cycle by one
      do_one_cycle(PhysicalInst_v, ready_queue, current_exec, g, current_exec_n, previously_completed);
      cycle++;
      current_leaf_cycle++;
      // record usage of all qubits, and eprs/zeros/magics specifically.
      qbits_per_cycle[cycle] = q_map.size();  
      unsigned int live_zero_count = 0;      
      unsigned int live_epr_count = 0;
      unsigned int live_magic_count = 0;
      std::vector<unsigned int> cycle_storage(SIMD_K, 0);
      for (auto &it : q_map) {
        if (it.first.find("zero") != std::string::npos)
          live_zero_count++;        
        if (it.first.find("epr") != std::string::npos)
          live_epr_count++;
        if (it.first.find("magic") != std::string::npos)
          live_magic_count++;        
        unsigned int loc = (it.second)->loc;
        cycle_storage[(loc-1)]++;
      }
      SeqTupleTy T = std::make_tuple(live_zero_count, live_epr_count, live_magic_count);
      ancillas_per_cycle[cycle] = T;
      zeros_per_cycle[cycle] = live_zero_count;      
      eprs_per_cycle[cycle] = live_epr_count;
      magics_per_cycle[cycle] = live_magic_count;      
      storage_per_cycle[cycle] = cycle_storage;

      #ifdef _DEBUG_SIMULATOR
        std::cout<<"*************************************"<<std::endl;      
        std::cout<<"current_exec.size() = "<<current_exec.size()<<std::endl;
        std::cout<<"current_exec_n.size() = "<<current_exec_n.size()<<std::endl;        
        std::cout<<"PhysicalInst_v.size() = "<<PhysicalInst_v.size()<<std::endl;
        std::cout<<"teleport_inst_tuples.size(): "<<teleport_inst_tuples.size()<<std::endl;      
        std::cout<<"*************************************"<<std::endl;        
      #endif

      // One or more currently executing insts. may have completed in this cycle
      bool completion_flag = 0;
      unsigned int completed_count = 0;    
      unsigned int completed_count_n = 0;    
      
      for (InstVecTy::iterator ce_I = current_exec.begin(); ce_I != current_exec.end(); ++ce_I) {
        if ((*ce_I)->is_complete) {
          completed_count++;
          completion_flag = 1;
          if ((*ce_I)->no_child) {
            if( q_map.find(((*ce_I)->qid)[0]) == q_map.end() ) {
              std::cerr<<"Error: Childless instruction complete, but no qbit to delete."<<std::endl;    
              return 1;
            }
            qbit* free_qbit = q_map[((*ce_I)->qid)[0]];  // need to only check first qubit. only BMOVs can be no_child instructions
            delete free_qbit;
          } 
        }
      }
      for (InstVecTy::iterator ce_I = current_exec_n.begin(); ce_I != current_exec_n.end(); ++ce_I) {
        if ((*ce_I)->is_complete) {
          completed_count_n++;
          completion_flag = 1;
          if ((*ce_I)->no_child) {          
            if( q_map.find(((*ce_I)->qid)[0]) == q_map.end() ) {
              std::cerr<<"Error: Childless instruction complete, but no qbit to delete."<<std::endl;    
              return 1;
            }          
            qbit* free_qbit = q_map[((*ce_I)->qid)[0]];  // need to only check first qubit. only BMOVs can be no_child instructions
            delete free_qbit;            
          }
        }
      }
      
      #ifdef _DEBUG_SIMULATOR            
        std::cout<<"Completed in this cycle: "<<completed_count<<" (current), "<<completed_count_n<<" (next)"<<std::endl;
      #endif
      
      // The use of flag has ensured that the following book-keeping is invoked only upon the completion of a current inst
      if (completion_flag){
        // which teleportations finished? tuples of teleport_inst_tuples that became is_complete.
        // 1. switch data & epr2; 2. delete epr1 and epr2 qbits 
        for (auto t = teleport_inst_tuples.begin(); t != teleport_inst_tuples.end(); ++t) {
          // can't find the instruction by doing PhysicalInst_v[seq] because Instructions are continuously removed from PhysicalInst_v.

          std::shared_ptr<Instruction> data_mov = std::get<0>(*t);
          std::shared_ptr<Instruction> epr1_mov = std::get<1>(*t);
          std::shared_ptr<Instruction> epr2_mov = std::get<2>(*t);

          if( (q_map.find(((data_mov)->qid)[0]) == q_map.end()) ||
              (q_map.find(((epr1_mov)->qid)[0]) == q_map.end()) ||
              (q_map.find(((epr2_mov)->qid)[0]) == q_map.end()) ) {
            continue;
          }

          qbit* data = q_map[((data_mov)->qid)[0]];
          qbit* epr1 = q_map[((epr1_mov)->qid)[0]];
          qbit* epr2 = q_map[((epr2_mov)->qid)[0]];       

          if (((data_mov)->qid).size()!=1 || data->kind != DATA) {
            std::cerr<<"Error: In tuple, "<<data->q_id<<" must be DATA."<<std::endl;    
            return 1;          
          }
          if (((epr1_mov)->qid).size()!=1 || epr1->kind != EPR1) {
            std::cerr<<"Error: In tuple, "<<epr1->q_id<<" must be EPR1."<<std::endl;                  
            return 1;
          }
          if (((epr2_mov)->qid).size()!=1 || epr2->kind != EPR2) {
            std::cerr<<"Error: In tuple, "<<epr2->q_id<<" must be EPR2."<<std::endl;                  
            return 1;
          }

          if ((data_mov)->is_complete && (epr1_mov)->is_complete && (epr2_mov)->is_complete) {
            #ifdef _DEBUG_QBIT_MANAGEMENT
              std::cout<<"Teleportation complete: "<<"<"<<data->q_id<<","<<epr1->q_id<<","<<epr2->q_id<<">"<<std::endl;
            #endif
            // 1. everything that points to data should now point to epr2 (switch qbit objects, not pointers!)          
            #ifdef _DEBUG_QBIT_MANAGEMENT
              std::cout<<"Exchanging: "<<data->q_id<<" & "<<epr2->q_id<<std::endl;            
              data->status();
              epr2->status();
            #endif            
            *data = *epr2;
            data->q_id = epr2->q_id.substr(0, epr2->q_id.find_first_of("_epr"));
            data->kind = DATA;

            // 2. the original data and epr1 are now obsolete qbits and will be deleted (recycled)          
            delete epr1; 
            delete epr2;    
 
          }
          else {
            // teleportation not finished yet and thus not ready for erasure
            continue;
          }
          
        }
        // finally, delete from teleport_inst_tuples itself so next time we iterate over less
        teleport_inst_tuples.erase(std::remove_if(teleport_inst_tuples.begin(), teleport_inst_tuples.end(), 
              [](const std::tuple<std::shared_ptr<Instruction>,std::shared_ptr<Instruction>,std::shared_ptr<Instruction>> T) 
                                      { return ((std::get<0>(T))->is_complete &&
                                                (std::get<1>(T))->is_complete &&
                                                (std::get<2>(T))->is_complete); }), teleport_inst_tuples.end() ); 

        // delete the completed instruction from the following vectors: use c++11 erase-remove idiom
        current_exec.erase( std::remove_if(current_exec.begin(), current_exec.end(), 
                            [](const std::shared_ptr<Instruction> &I) { return I->is_complete; }), 
                            current_exec.end() );      
        current_exec_n.erase( std::remove_if(current_exec_n.begin(), current_exec_n.end(), 
                            [](const std::shared_ptr<Instruction> &I) { return I->is_complete; }), 
                            current_exec_n.end() );          
        PhysicalInst_v.erase(std::remove_if(PhysicalInst_v.begin(),PhysicalInst_v.end(),
                            [](const std::shared_ptr<Instruction> &I) { return I->is_complete; }),
                            PhysicalInst_v.end() );

      }   
    }
    leaf_cycles[leaf_name].push_back(current_leaf_cycle);
    current_leaf_cycle = 0;
  }  

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Report overall benchmark results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  std::string output_dir = benchmark_dir+"/simd_simulation/";
  std::string mkdir_command = "mkdir -p "+output_dir;
  system(mkdir_command.c_str());

  // Direction of Smoothing:
  std::string direction;
  //  Forward: Each move in its own module. Put a cap on the number of qubits released for move at each cycle.
  //           cap.inf.window.inf.forward: greedily send all movements of CURRENT module in its first cycle.
  if (!backward && forward)
    direction = std::string(".forward");
  //  Backward: Release 'cap' moves from the next module. No limit on moves from this module.
  //            cap.inf.window.inf.backward: greedily send all the possible qubits of CURRENT AND EPRs of NEXT module in the first cycle of this module. 
  else if (backward && !forward)
    direction = std::string(".backward");
  // Backforth: Release 'cap' moves from this module and 'cap' moves from the next module.
  //            cap.inf.window.inf.backforth: greedily send all the possible qubits of CURRENT AND EPRs of NEXT module in the first cycle of this module. 
  else if (backward && forward)
    direction = std::string(".backforth");
 
  // KQ: total number of logical gates
  // k: total number of physical timesteps
  // q: total number of physical qubits
  unsigned long long total_cycles = 0;
  for (auto &m_it : leaf_cycles) {
    unsigned long long leaf_avg_cycles = 0;
    unsigned long long leaf_freq = 0;
    std::string leaf_name = m_it.first; 
    std::cout << leaf_name << "\t\t";
    for (auto &v_it : m_it.second) {
      std::cout << v_it << "\t";
      leaf_avg_cycles += v_it;
    }
    std::cout << std::endl;  
    leaf_avg_cycles /= m_it.second.size();
    leaf_freq = (unsigned long long) mapOfFreq[leaf_name];    
    total_cycles += leaf_avg_cycles * leaf_freq;    
  }

  unsigned int max_q_count = 0;
  for (auto &it : qbits_per_cycle) {
    if (it.second > max_q_count) max_q_count = it.second;
  }
  unsigned long long num_physical_qbits = (unsigned long long) max_q_count * (unsigned long long)pow(7.0,(double)concatenation_level);
  

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Report results for plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  std::string kq_file_path;  
  std::string usage_file_path;
  std::string ages_file_path;
  std::string storage_file_path;
  std::ofstream kq_file;    
  std::ofstream usage_file;
  std::ofstream ages_file;
  std::ofstream storage_file;
 
  kq_file_path = output_dir+benchmark_name
                    +".p."+(std::to_string(P_error_rate))
                    +".cap."+(cap ? (std::to_string(mov_cap)) : std::string("inf"))
                    +".window."+(window ? (std::to_string(window_size)) : std::string("inf"))
                    +direction                             
                    +".kq";
  kq_file.open(kq_file_path);

  kq_file << "error rate: " << "10^-" << P_error_rate << std::endl;
  kq_file << "concatenation level: " << concatenation_level << std::endl;
  kq_file << "total cycles: " << total_cycles << std::endl;
  kq_file << "max qubits: " << num_physical_qbits << std::endl;
  kq_file << "logical KQ: " << total_logical_gates << std::endl;
  kq_file << "physical kq: " << (total_cycles * num_physical_qbits) << std::endl;  
  kq_file.close();  
  std::cerr << "kq report written to:\t" << kq_file_path << std::endl; 

  if (report_usage) {
    usage_file_path = output_dir+benchmark_name
                        +".p."+(std::to_string(P_error_rate))      
                        +".cap."+(cap ? (std::to_string(mov_cap)) : std::string("inf"))
                        +".window."+(window ? (std::to_string(window_size)) : std::string("inf"))
                        +direction                              
                        +".usage";
    usage_file.open(usage_file_path);
    for (auto &it : ancillas_per_cycle) {
      usage_file << it.first << "\t" << std::get<0>(it.second) + std::get<1>(it.second) + std::get<2>(it.second) << std::endl; 
      //usage_file << it.first << "\t" << std::get<0>(it.second) << "\t" << std::get<1>(it.second) << "\t" << std::get<2>(it.second) << std::endl; 
    }
    usage_file.close();
    std::cerr << "usage report written to:\t" << usage_file_path << std::endl;      
  }

  if (report_ages) {
    ages_file_path = output_dir+benchmark_name
                        +".p."+(std::to_string(P_error_rate))      
                        +".cap."+(cap ? (std::to_string(mov_cap)) : std::string("inf"))
                        +".window."+(window ? (std::to_string(window_size)) : std::string("inf"))
                        +direction                             
                        +".ages";
    ages_file.open(ages_file_path);
    double avg_ancilla_age=0;
    int peak_ancilla_age=0;
    int ancilla_count=0;
    for (auto &it : qbit_ages) {
      if( (it.first.find("zero") != std::string::npos) ||
          (it.first.find("epr") != std::string::npos) ||          
          (it.first.find("magic") != std::string::npos) ) {
        avg_ancilla_age += it.second;
        ancilla_count++;
        if (it.second > peak_ancilla_age)
          peak_ancilla_age = it.second;
      }
    } 
    ages_file<<avg_ancilla_age/ancilla_count;
    ages_file<<"\t"<<peak_ancilla_age<<"\n";    
    ages_file.close();  
    std::cerr << "ages report written to:\t" << ages_file_path << std::endl;  
  }

  if (report_storage) {
    storage_file_path = output_dir+benchmark_name
                        +".p."+(std::to_string(P_error_rate))      
                        +".cap."+(cap ? (std::to_string(mov_cap)) : std::string("inf"))
                        +".window."+(window ? (std::to_string(window_size)) : std::string("inf"))
                        +direction
                        +".storage";        
    storage_file.open(storage_file_path);
    std::vector<double> peak_storage (SIMD_K, 0);
    for (auto &it : storage_per_cycle)
      for (int k = 0; k<SIMD_K; k++)
        if (peak_storage[k] < (it.second)[k])
          peak_storage[k] = (it.second)[k];
    for (int k = 0; k<SIMD_K; k++)
      storage_file<<peak_storage[k]<<"\t";
    storage_file.close();
    std::cerr << "storage report written to:\t" << storage_file_path << std::endl;      
  }

  return 0;  
}
