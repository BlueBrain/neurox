#include <cstring>

#include "neurox/neurox.h"

#include "coreneuron/nrnoc/membfunc.h"

using namespace std;
using namespace neurox;
using namespace neurox::interpolators;

Mechanism::Mechanism(const int type, const short int data_size,
                     const short int pdata_size, const char is_artificial,
                     char pnt_map, const char is_ion,
                     const short int sym_length, const char *sym,
                     Memb_func &memb_func, const short int dependencies_count,
                     const int *dependencies, const short int successors_count,
                     const int *successors)
    : type_(type),
      data_size_(data_size),
      pdata_size_(pdata_size),
      vdata_size_(0),
      successors_count_(successors_count),
      dependencies_count_(dependencies_count),
      sym_length_(sym_length),
      pnt_map_(pnt_map),
      is_artificial_(is_artificial),
      is_ion_(is_ion),
      dependencies_(nullptr),
      successors_(nullptr),
      state_vars_(nullptr),
      pnt_receive_(nullptr),
      ode_matsol_(nullptr),
      ode_spec_(nullptr),
      div_capacity_(nullptr),
      mul_capacity_(nullptr),
      current_func_runtime_(-1),
      state_func_runtime_(-1) {
  // to be set by neuronx::UpdateMechanismsDependencies
  this->dependency_ion_index_ = Mechanism::IonTypes::kNoIon;

  // set function pointers
  memcpy(&this->memb_func_, &memb_func, sizeof(Memb_func));
  assert(sym_length > 0);

  // set name (overwrites existing, pointing to CN data structures)
  this->memb_func_.sym = new char[sym_length + 1];
  std::memcpy(this->memb_func_.sym, sym, sym_length);
  this->memb_func_.sym[sym_length] = '\0';

  // copy dparam_semantics (overwrites existing, pointing to CN data structures)
  this->memb_func_.dparam_semantics = new int[pdata_size];
  memcpy(this->memb_func_.dparam_semantics, memb_func.dparam_semantics,
         sizeof(int) * pdata_size);

  // non-coreneuron functions
  if (this->type_ != CAP && !this->is_ion_) {
    this->memb_func_.current = get_cur_function(this->memb_func_.sym);
    this->pnt_receive_ = get_net_receive_function(this->memb_func_.sym);
    this->memb_func_.jacob = get_jacob_function(this->memb_func_.sym);
  } else if (this->type_ == CAP)  // capacitance: capac.c
  {
    this->memb_func_.current = nrn_cur_capacitance;
    this->memb_func_.jacob = nrn_jacob_capacitance;
    this->div_capacity_ = nrn_div_capacity;  // CVODE-specific
    this->mul_capacity_ = nrn_mul_capacity;  // CVODE-specific
  } else if (this->is_ion_)                  // ion: eion.c
  {
    this->memb_func_.current = nrn_cur_ion;
  }

  vdata_size_ = input::DataLoader::HardCodedVdataCount(this->type_, pnt_map);

  // CVODES-specific
  if (input_params_->interpolator_ !=
      interpolators::InterpolatorIds::kBackwardEuler) {
    this->state_vars_ = new StateVars();
    if (!this->is_ion_ && this->type_ != MechanismTypes::kCapacitance) {
      // get state variables count, values and offsets
      state_vars_f_t stf = get_ode_state_vars_function(this->memb_func_.sym);
      // if (stf != NULL)
      stf(&this->state_vars_->count_, &this->state_vars_->var_offsets_,
          &this->state_vars_->dv_offsets_);

      // state variables diagonal at given point
      this->ode_matsol_ = get_ode_matsol_function(this->memb_func_.sym);

      // derivative description
      this->ode_spec_ = get_ode_spec_function(this->memb_func_.sym);
    }
  }

  this->memb_func_.is_point = pnt_map > 0 ? 1 : 0;

  if (dependencies != nullptr) {
    assert(dependencies_count > 0);
    this->dependencies_ = new int[dependencies_count];
    std::memcpy(this->dependencies_, dependencies,
                dependencies_count * sizeof(int));
  }

  if (successors != nullptr) {
    assert(successors_count > 0);
    this->successors_ = new int[successors_count];
    std::memcpy(this->successors_, successors, successors_count * sizeof(int));
  }
};

Mechanism::StateVars::StateVars()
    : count_(0), var_offsets_(nullptr), dv_offsets_(nullptr) {}

Mechanism::StateVars::StateVars(int count, int *offsets, int *dv_offsets)
    : count_(count), var_offsets_(offsets), dv_offsets_(dv_offsets) {}

Mechanism::StateVars::~StateVars() {
  delete[] var_offsets_;
  delete[] dv_offsets_;
}

Mechanism::IonTypes Mechanism::GetIonIndex() {
  assert(this->memb_func_.sym);
  if (strcmp("na_ion", this->memb_func_.sym) == 0)
    return Mechanism::IonTypes::kNa;
  if (strcmp("k_ion", this->memb_func_.sym) == 0)
    return Mechanism::IonTypes::kK;
  if (strcmp("ttx_ion", this->memb_func_.sym) == 0)
    return Mechanism::IonTypes::kTTX;
  if (strcmp("ca_ion", this->memb_func_.sym) == 0)
    return Mechanism::IonTypes::kCa;
  return Mechanism::IonTypes::kNoIon;
}

void Mechanism::RegisterBeforeAfterFunctions() {
  // Copy Before-After functions
  // register_mech.c::hoc_reg_ba()
  for (int i = 0; i < BEFORE_AFTER_SIZE; i++)
    this->before_after_functions_[i] =
        get_BA_function(this->memb_func_.sym, i);  // NULL at the moment
  // this->BAfunctions[i] = nrn_threads[0].tbl[i]->bam->f;
}

Mechanism::~Mechanism() {
  delete[] memb_func_.sym;
  delete[] successors_;
  delete[] dependencies_;
}

int Mechanism::ModFunctionStateThread(int i, void *args_ptr) {
  MembListThreadArgs *args = (MembListThreadArgs *)args_ptr;
  args->memb_func->state(args->nt, &args->ml_state[i], args->mech_type);
  return 0;
}

int Mechanism::ModFunctionCurrentThread(int i, void *args_ptr) {
  MembListThreadArgs *args = (MembListThreadArgs *)args_ptr;
  if (args->requires_shadow_vectors)
    args->memb_func->current(args->nt, &args->ml_current[i], args->mech_type,
                             args->acc_rhs_d, args->acc_di_dv, args->acc_args);
  else
    args->memb_func->current(args->nt, &args->ml_current[i], args->mech_type,
                             NULL, NULL, NULL);
  return 0;
}

void Mechanism::CallModFunction(
    const void *branch_ptr, const Mechanism::ModFunctions function_id,
    Memb_list *other_ml,    // other Memb_list (if any)
    const NetconX *netcon,  // for net_receive only
    const floble_t tt       // for net_receive only
) {
  const Branch *branch = (Branch *)branch_ptr;
  assert(branch);
  NrnThread *nt = branch->nt_;
  assert(nt);

  // TODO: support parallel mechs for CVODE
  // TODO in Vectorizer we need to support 'other_ml' for threaded execution
  // TODO instance parallelism for mech_receive?
  if (other_ml) {
    assert(input_params_->interpolator_ == InterpolatorIds::kBackwardEuler);
  }

  /* if net receive, it's a special function type: process now */
  if (function_id == Mechanism::ModFunctions::kNetReceive ||
      function_id == Mechanism::ModFunctions::kNetReceiveInit) {
    assert(function_id != Mechanism::ModFunctions::kNetReceiveInit);  // N/A yet
    assert(this->pnt_receive_);

    Memb_list *memb_list =
        other_ml
            ? &other_ml[mechanisms_map_[netcon->mech_type_]]
            : &branch->mechs_instances_[mechanisms_map_[netcon->mech_type_]];
    assert(memb_list);
    int iml = netcon->mech_instance_;
    int weight_index = netcon->weight_index_;
    nt->_t = tt;  // as seen in netcvode.cpp (NetCon::deliver)
    this->pnt_receive_(nt, this->type_, iml, weight_index, 0,
                       Branch::CoreneuronNetSend, Branch::CoreneuronNetEvent);
    return;
  }

  /* memb_list to be used is the one with all instances or
   * an user-provided one (eg CVODE provides non-capacitors memb_list */
  const int mech_offset = neurox::mechanisms_map_[this->type_];
  Memb_list *memb_list = other_ml ? &other_ml[mech_offset]
                                  : &branch->mechs_instances_[mech_offset];
  assert(memb_list);

  /* create argument for thread instances of current and state functions */
  MembListThreadArgs *threads_args = nullptr;

  if (branch->mechs_instances_parallel_)
    if (function_id == ModFunctions::kState ||
        function_id == ModFunctions::kCurrent ||
        function_id == ModFunctions::kCurrentCapacitance)
      threads_args = &branch->mechs_instances_parallel_[mech_offset];

  if (memb_list->nodecount > 0) {
    switch (function_id) {
      case Mechanism::ModFunctions::kBeforeInitialize:
      case Mechanism::ModFunctions::kAfterInitialize:
      case Mechanism::ModFunctions::kBeforeBreakpoint:
      case Mechanism::ModFunctions::kAfterSolve:
      case Mechanism::ModFunctions::kBeforeStep:
        if (before_after_functions_[(int)function_id])
          before_after_functions_[(int)function_id](nt, memb_list, type_);
        break;
      case Mechanism::ModFunctions::kAlloc:
        if (memb_func_.alloc)
          memb_func_.alloc(memb_list->data, memb_list->pdata, type_);
        break;
      case Mechanism::ModFunctions::kCurrentCapacitance: {
        assert(type_ == MechanismTypes::kCapacitance);
        assert(memb_func_.current != NULL);
        if (threads_args && threads_args->ml_current_count > 1)
          hpx_par_for_sync(Mechanism::ModFunctionCurrentThread, 0,
                           threads_args->ml_current_count, threads_args);
        else
          memb_func_.current(nt, memb_list, type_, NULL, NULL, NULL);
        break;
      }
      case Mechanism::ModFunctions::kCurrent: {
        assert(type_ != MechanismTypes::kCapacitance);
        if (memb_func_.current) {
          if (threads_args && threads_args->ml_current_count > 1) {
            hpx_par_for_sync(Mechanism::ModFunctionCurrentThread, 0,
                             threads_args->ml_current_count, threads_args);
          } else {
            bool requires_shadow_vectors =
                /* graph-parallelism */
                branch->mechs_graph_
                /* not CaDynamics_E2 (no updates in cur function) */
                && this->type_ != MechanismTypes::kCaDynamics_E2
                /* not ion (updates in nrn_cur_ion function) */
                && (!this->is_ion_);

            if (requires_shadow_vectors) {
              mod_acc_f_t acc_rhs_d =
                  Branch::MechanismsGraph::AccumulateRHSandD;
              mod_acc_f_t acc_di_dv =
                  this->dependencies_count_
                      ? Branch::MechanismsGraph::AccumulateIandDIDV
                      : nullptr;
              void *args = branch->mechs_graph_;
              memb_func_.current(nt, memb_list, type_, acc_rhs_d, acc_di_dv,
                                 args);
            } else {
              memb_func_.current(nt, memb_list, type_, NULL, NULL, NULL);
            }
          }
        }
      } break;
      case Mechanism::ModFunctions::kState: {
        if (memb_func_.state) {
          if (threads_args && threads_args->ml_state_count > 1) {
            hpx_par_for_sync(Mechanism::ModFunctionStateThread, 0,
                             threads_args->ml_state_count, threads_args);
          } else {
            memb_func_.state(nt, memb_list, type_);
          }
        }
      } break;
      case Mechanism::ModFunctions::kJacobCapacitance:
        assert(type_ == MechanismTypes::kCapacitance);
        assert(memb_func_.jacob != NULL);
        nrn_jacob_capacitance(nt, memb_list, type_);
        break;
      case Mechanism::ModFunctions::kJacob:
        assert(type_ != MechanismTypes::kCapacitance);
        if (memb_func_.jacob) {
          assert(0);  // No jacob function pointers yet
                      // (get_jacob_function(xxx))
          memb_func_.jacob(nt, memb_list, type_);
        }
        break;
      case Mechanism::ModFunctions::kInitialize:
        if (memb_func_.initialize) memb_func_.initialize(nt, memb_list, type_);
        break;
      case Mechanism::ModFunctions::kDestructor:
        if (memb_func_.destructor) memb_func_.destructor();
        break;
      case Mechanism::ModFunctions::kThreadMemInit:
        assert(0);  // should be called only by constructor Branch(...)
        if (memb_func_.thread_mem_init_)
          memb_func_.thread_mem_init_(memb_list->_thread);
        break;
      case Mechanism::ModFunctions::kThreadTableCheck:
        if (memb_func_.thread_table_check_)
          memb_func_.thread_table_check_(0, memb_list->nodecount,
                                         memb_list->data, memb_list->pdata,
                                         memb_list->_thread, nt, type_);
        break;
      case Mechanism::ModFunctions::kThreadCleanup:
        assert(0);  // should only be called by destructor ~Branch(...)
        if (memb_func_.thread_cleanup_)
          memb_func_.thread_cleanup_(memb_list->_thread);
        break;
      case Mechanism::ModFunctions::kODEMatsol:  // CVODE-specific
        if (this->ode_matsol_ && this->state_vars_->count_ > 0)
          tools::Vectorizer::CallVecFunction(this->ode_matsol_, nt, memb_list,
                                             type_);
        break;
      case Mechanism::ModFunctions::kODESpec:  // CVODE-specific
        if (this->ode_spec_ && this->state_vars_->count_ > 0) {
          for (int i = 0; i < memb_list->nodecount; i++) {
            int nd = memb_list->nodeindices[i];
            // fprintf(stderr, "== odespec: mech %d, node %d\n", this->type_,
            // nd);
          }
          tools::Vectorizer::CallVecFunction(this->ode_spec_, nt, memb_list,
                                             type_);
        }
        break;
      case Mechanism::ModFunctions::kDivCapacity:  // CVODE-specific
        if (this->div_capacity_) assert(type_ == MechanismTypes::kCapacitance);
        assert(this->div_capacity_ != NULL);
        nrn_div_capacity(nt, memb_list, type_);
        break;
      case Mechanism::ModFunctions::kMulCapacity:  // CVODE-specific
        if (this->mul_capacity_) assert(type_ == MechanismTypes::kCapacitance);
        assert(this->div_capacity_ != NULL);
        nrn_mul_capacity(nt, memb_list, type_);
        break;
      default:
        printf("ERROR: Unknown ModFunction with id %d.\n", function_id);
        exit(1);
    }
  }
}
