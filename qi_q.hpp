
// Code generated by stanc 0016615
#include <stan/model/model_header.hpp>
namespace qi_q_model_namespace {

using stan::io::dump;
using stan::model::assign;
using stan::model::index_uni;
using stan::model::index_max;
using stan::model::index_min;
using stan::model::index_min_max;
using stan::model::index_multi;
using stan::model::index_omni;
using stan::model::model_base_crtp;
using stan::model::rvalue;
using namespace stan::math;


stan::math::profile_map profiles__;
static constexpr std::array<const char*, 23> locations_array__ = 
{" (found before start of program)",
 " (in '/Users/sean00002/Dropbox/PhD/Andrew_Rotation/qi_q.stan', line 28, column 3 to column 39)",
 " (in '/Users/sean00002/Dropbox/PhD/Andrew_Rotation/qi_q.stan', line 29, column 3 to column 20)",
 " (in '/Users/sean00002/Dropbox/PhD/Andrew_Rotation/qi_q.stan', line 35, column 3 to column 26)",
 " (in '/Users/sean00002/Dropbox/PhD/Andrew_Rotation/qi_q.stan', line 37, column 6 to column 26)",
 " (in '/Users/sean00002/Dropbox/PhD/Andrew_Rotation/qi_q.stan', line 39, column 9 to column 41)",
 " (in '/Users/sean00002/Dropbox/PhD/Andrew_Rotation/qi_q.stan', line 38, column 6 to line 39, column 41)",
 " (in '/Users/sean00002/Dropbox/PhD/Andrew_Rotation/qi_q.stan', line 36, column 26 to line 40, column 4)",
 " (in '/Users/sean00002/Dropbox/PhD/Andrew_Rotation/qi_q.stan', line 36, column 3 to line 40, column 4)",
 " (in '/Users/sean00002/Dropbox/PhD/Andrew_Rotation/qi_q.stan', line 22, column 3 to column 22)",
 " (in '/Users/sean00002/Dropbox/PhD/Andrew_Rotation/qi_q.stan', line 23, column 3 to column 27)",
 " (in '/Users/sean00002/Dropbox/PhD/Andrew_Rotation/qi_q.stan', line 24, column 28 to column 38)",
 " (in '/Users/sean00002/Dropbox/PhD/Andrew_Rotation/qi_q.stan', line 24, column 39 to column 44)",
 " (in '/Users/sean00002/Dropbox/PhD/Andrew_Rotation/qi_q.stan', line 24, column 3 to column 46)",
 " (in '/Users/sean00002/Dropbox/PhD/Andrew_Rotation/qi_q.stan', line 28, column 27 to column 37)",
 " (in '/Users/sean00002/Dropbox/PhD/Andrew_Rotation/qi_q.stan', line 5, column 6 to column 43)",
 " (in '/Users/sean00002/Dropbox/PhD/Andrew_Rotation/qi_q.stan', line 6, column 6 to column 19)",
 " (in '/Users/sean00002/Dropbox/PhD/Andrew_Rotation/qi_q.stan', line 7, column 6 to column 34)",
 " (in '/Users/sean00002/Dropbox/PhD/Andrew_Rotation/qi_q.stan', line 4, column 51 to line 8, column 4)",
 " (in '/Users/sean00002/Dropbox/PhD/Andrew_Rotation/qi_q.stan', line 12, column 6 to column 54)",
 " (in '/Users/sean00002/Dropbox/PhD/Andrew_Rotation/qi_q.stan', line 11, column 51 to line 13, column 4)",
 " (in '/Users/sean00002/Dropbox/PhD/Andrew_Rotation/qi_q.stan', line 16, column 6 to column 50)",
 " (in '/Users/sean00002/Dropbox/PhD/Andrew_Rotation/qi_q.stan', line 15, column 51 to line 17, column 4)"};


template <bool propto__, typename T0__, typename T1__, typename T2__>
stan::promote_args_t<T0__, T1__,
T2__>
gammaModeSD_lpdf(const T0__& parm, const T1__& m, const T2__& sd,
                 std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<T0__, T1__, T2__>;
  int current_statement__ = 0; 
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    local_scalar_t__ r;
    r = DUMMY_VAR__;
    
    current_statement__ = 15;
    r = ((m + stan::math::sqrt((pow(m, 2) + (4 * pow(sd, 2))))) /
          (2 * pow(sd, 2)));
    local_scalar_t__ s;
    s = DUMMY_VAR__;
    
    current_statement__ = 16;
    s = (1 + (m * r));
    current_statement__ = 17;
    return gamma_lpdf<false>(parm, r, s);
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
  }
  
}

struct gammaModeSD_lpdf_functor__ {
template <bool propto__, typename T0__, typename T1__, typename T2__>
stan::promote_args_t<T0__, T1__,
T2__>
operator()(const T0__& parm, const T1__& m, const T2__& sd,
           std::ostream* pstream__)  const 
{
return gammaModeSD_lpdf<propto__>(parm, m, sd, pstream__);
}
};

template <bool propto__, typename T0__, typename T1__, typename T2__>
stan::promote_args_t<T0__, T1__,
T2__>
betaModeConc_lpdf(const T0__& parm, const T1__& m, const T2__& c,
                  std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<T0__, T1__, T2__>;
  int current_statement__ = 0; 
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    current_statement__ = 19;
    return beta_lpdf<false>(parm, ((m * (c - 2)) + 1),
             (((1 - m) * (c - 2)) + 1));
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
  }
  
}

struct betaModeConc_lpdf_functor__ {
template <bool propto__, typename T0__, typename T1__, typename T2__>
stan::promote_args_t<T0__, T1__,
T2__>
operator()(const T0__& parm, const T1__& m, const T2__& c,
           std::ostream* pstream__)  const 
{
return betaModeConc_lpdf<propto__>(parm, m, c, pstream__);
}
};

template <bool propto__, typename T0__, typename T1__, typename T2__>
stan::promote_args_t<T0__, T1__,
T2__>
betaMeanConc_lpdf(const T0__& parm, const T1__& m, const T2__& c,
                  std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<T0__, T1__, T2__>;
  int current_statement__ = 0; 
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  try {
    current_statement__ = 21;
    return beta_lpdf<false>(parm, (m * (c - 2)), ((1 - m) * (c - 2)));
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
  }
  
}

struct betaMeanConc_lpdf_functor__ {
template <bool propto__, typename T0__, typename T1__, typename T2__>
stan::promote_args_t<T0__, T1__,
T2__>
operator()(const T0__& parm, const T1__& m, const T2__& c,
           std::ostream* pstream__)  const 
{
return betaMeanConc_lpdf<propto__>(parm, m, c, pstream__);
}
};

class qi_q_model final : public model_base_crtp<qi_q_model> {

 private:
  int N_RNA;
  int N_VARIANTS;
  std::vector<std::vector<double>> qi; 
  
 
 public:
  ~qi_q_model() { }
  
  inline std::string model_name() const final { return "qi_q_model"; }

  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = stanc3 0016615", "stancflags = "};
  }
  
  
  qi_q_model(stan::io::var_context& context__,
             unsigned int random_seed__ = 0,
             std::ostream* pstream__ = nullptr) : model_base_crtp(0) {
    int current_statement__ = 0;
    using local_scalar_t__ = double ;
    boost::ecuyer1988 base_rng__ = 
        stan::services::util::create_rng(random_seed__, 0);
    (void) base_rng__;  // suppress unused var warning
    static constexpr const char* function__ = "qi_q_model_namespace::qi_q_model";
    (void) function__;  // suppress unused var warning
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    try {
      int pos__;
      pos__ = std::numeric_limits<int>::min();
      
      pos__ = 1;
      current_statement__ = 9;
      context__.validate_dims("data initialization","N_RNA","int",
           std::vector<size_t>{});
      N_RNA = std::numeric_limits<int>::min();
      
      current_statement__ = 9;
      N_RNA = context__.vals_i("N_RNA")[(1 - 1)];
      current_statement__ = 9;
      check_greater_or_equal(function__, "N_RNA", N_RNA, 1);
      current_statement__ = 10;
      context__.validate_dims("data initialization","N_VARIANTS","int",
           std::vector<size_t>{});
      N_VARIANTS = std::numeric_limits<int>::min();
      
      current_statement__ = 10;
      N_VARIANTS = context__.vals_i("N_VARIANTS")[(1 - 1)];
      current_statement__ = 10;
      check_greater_or_equal(function__, "N_VARIANTS", N_VARIANTS, 1);
      current_statement__ = 11;
      validate_non_negative_index("qi", "N_VARIANTS", N_VARIANTS);
      current_statement__ = 12;
      validate_non_negative_index("qi", "N_RNA", N_RNA);
      current_statement__ = 13;
      context__.validate_dims("data initialization","qi","double",
           std::vector<size_t>{static_cast<size_t>(N_VARIANTS),
            static_cast<size_t>(N_RNA)});
      qi = std::vector<std::vector<double>>(N_VARIANTS, std::vector<double>(N_RNA, std::numeric_limits<double>::quiet_NaN()));
      
      
      {
        std::vector<local_scalar_t__> qi_flat__;
        current_statement__ = 13;
        qi_flat__ = context__.vals_r("qi");
        current_statement__ = 13;
        pos__ = 1;
        current_statement__ = 13;
        for (int sym1__ = 1; sym1__ <= N_RNA; ++sym1__) {
          current_statement__ = 13;
          for (int sym2__ = 1; sym2__ <= N_VARIANTS; ++sym2__) {
            current_statement__ = 13;
            assign(qi, qi_flat__[(pos__ - 1)],
              "assigning variable qi", index_uni(sym2__), index_uni(sym1__));
            current_statement__ = 13;
            pos__ = (pos__ + 1);
          }
        }
      }
      current_statement__ = 13;
      for (int sym1__ = 1; sym1__ <= N_VARIANTS; ++sym1__) {
        current_statement__ = 13;
        for (int sym2__ = 1; sym2__ <= N_RNA; ++sym2__) {
          current_statement__ = 13;
          check_greater_or_equal(function__, "qi[sym1__, sym2__]",
                                 qi[(sym1__ - 1)][(sym2__ - 1)], 0);
        }
      }
      current_statement__ = 13;
      for (int sym1__ = 1; sym1__ <= N_VARIANTS; ++sym1__) {
        current_statement__ = 13;
        for (int sym2__ = 1; sym2__ <= N_RNA; ++sym2__) {
          current_statement__ = 13;
          check_less_or_equal(function__, "qi[sym1__, sym2__]",
                              qi[(sym1__ - 1)][(sym2__ - 1)], 1);
        }
      }
      current_statement__ = 14;
      validate_non_negative_index("q", "N_VARIANTS", N_VARIANTS);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    num_params_r__ = N_VARIANTS + 1;
    
  }
  
  template <bool propto__, bool jacobian__ , typename VecR, typename VecI, 
  stan::require_vector_like_t<VecR>* = nullptr, 
  stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr> 
  inline stan::scalar_type_t<VecR> log_prob_impl(VecR& params_r__,
                                                 VecI& params_i__,
                                                 std::ostream* pstream__ = nullptr) const {
    using T__ = stan::scalar_type_t<VecR>;
    using local_scalar_t__ = T__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    static constexpr const char* function__ = "qi_q_model_namespace::log_prob";
    (void) function__;  // suppress unused var warning
    
    try {
      std::vector<local_scalar_t__> q;
      q = std::vector<local_scalar_t__>(N_VARIANTS, DUMMY_VAR__);
      
      current_statement__ = 1;
      q = in__.template read_constrain_lub<std::vector<local_scalar_t__>, jacobian__>(
            0, 1, lp__, N_VARIANTS);
      local_scalar_t__ c2;
      c2 = DUMMY_VAR__;
      
      current_statement__ = 2;
      c2 = in__.template read_constrain_lb<local_scalar_t__, jacobian__>(2,
             lp__);
      {
        current_statement__ = 3;
        lp_accum__.add(gamma_lpdf<propto__>(c2, 1.1, 0.005));
        current_statement__ = 8;
        for (int j = 1; j <= N_VARIANTS; ++j) {
          current_statement__ = 4;
          lp_accum__.add(
            uniform_lpdf<propto__>(rvalue(q, "q", index_uni(j)), 0, 1));
          current_statement__ = 6;
          for (int i = 1; i <= N_RNA; ++i) {
            current_statement__ = 5;
            lp_accum__.add(
              betaModeConc_lpdf<propto__>(
                rvalue(qi, "qi", index_uni(j), index_uni(i)),
                rvalue(q, "q", index_uni(j)), c2, pstream__));
          }
        }
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
    } // log_prob_impl() 
    
  template <typename RNG, typename VecR, typename VecI, typename VecVar, 
  stan::require_vector_like_vt<std::is_floating_point, VecR>* = nullptr, 
  stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr, 
  stan::require_std_vector_vt<std::is_floating_point, VecVar>* = nullptr> 
  inline void write_array_impl(RNG& base_rng__, VecR& params_r__,
                               VecI& params_i__, VecVar& vars__,
                               const bool emit_transformed_parameters__ = true,
                               const bool emit_generated_quantities__ = true,
                               std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    vars__.resize(0);
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    static constexpr bool propto__ = true;
    (void) propto__;
    double lp__ = 0.0;
    (void) lp__;  // dummy to suppress unused var warning
    int current_statement__ = 0; 
    stan::math::accumulator<double> lp_accum__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    constexpr bool jacobian__ = false;
    (void) DUMMY_VAR__;  // suppress unused var warning
    static constexpr const char* function__ = "qi_q_model_namespace::write_array";
    (void) function__;  // suppress unused var warning
    
    try {
      std::vector<double> q;
      q = std::vector<double>(N_VARIANTS, std::numeric_limits<double>::quiet_NaN());
      
      
      current_statement__ = 1;
      q = in__.template read_constrain_lub<std::vector<local_scalar_t__>, jacobian__>(
            0, 1, lp__, N_VARIANTS);
      double c2;
      c2 = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 2;
      c2 = in__.template read_constrain_lb<local_scalar_t__, jacobian__>(2,
             lp__);
      for (int sym1__ = 1; sym1__ <= N_VARIANTS; ++sym1__) {
        vars__.emplace_back(q[(sym1__ - 1)]);
      }
      vars__.emplace_back(c2);
      if (logical_negation((primitive_value(emit_transformed_parameters__) ||
            primitive_value(emit_generated_quantities__)))) {
        return ;
      } 
      if (logical_negation(emit_generated_quantities__)) {
        return ;
      } 
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // write_array_impl() 
    
  template <typename VecVar, typename VecI, 
  stan::require_std_vector_t<VecVar>* = nullptr, 
  stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr> 
  inline void transform_inits_impl(const stan::io::var_context& context__,
                                   VecI& params_i__, VecVar& vars__,
                                   std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    vars__.clear();
    vars__.reserve(num_params_r__);
    int current_statement__ = 0; 
    
    try {
      int pos__;
      pos__ = std::numeric_limits<int>::min();
      
      pos__ = 1;
      std::vector<double> q;
      q = std::vector<double>(N_VARIANTS, std::numeric_limits<double>::quiet_NaN());
      
      
      current_statement__ = 1;
      q = context__.vals_r("q");
      std::vector<double> q_free__;
      q_free__ = std::vector<double>(N_VARIANTS, std::numeric_limits<double>::quiet_NaN());
      
      
      current_statement__ = 1;
      for (int sym1__ = 1; sym1__ <= N_VARIANTS; ++sym1__) {
        current_statement__ = 1;
        assign(q_free__, stan::math::lub_free(q[(sym1__ - 1)], 0, 1),
          "assigning variable q_free__", index_uni(sym1__));
      }
      double c2;
      c2 = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 2;
      c2 = context__.vals_r("c2")[(1 - 1)];
      double c2_free__;
      c2_free__ = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 2;
      c2_free__ = stan::math::lb_free(c2, 2);
      for (int sym1__ = 1; sym1__ <= N_VARIANTS; ++sym1__) {
        vars__.emplace_back(q_free__[(sym1__ - 1)]);
      }
      vars__.emplace_back(c2_free__);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // transform_inits_impl() 
    
  inline void get_param_names(std::vector<std::string>& names__) const {
    
    names__ = std::vector<std::string>{"q", "c2"};
    
    } // get_param_names() 
    
  inline void get_dims(std::vector<std::vector<size_t>>& dimss__) const {
    
    dimss__ = std::vector<std::vector<size_t>>{std::vector<size_t>{
                                                                   static_cast<size_t>(N_VARIANTS)
                                                                   },
      std::vector<size_t>{}};
    
    } // get_dims() 
    
  inline void constrained_param_names(
                                      std::vector<std::string>& param_names__,
                                      bool emit_transformed_parameters__ = true,
                                      bool emit_generated_quantities__ = true) const
    final {
    
    for (int sym1__ = 1; sym1__ <= N_VARIANTS; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "q" + '.' + std::to_string(sym1__));
      } 
    }
    param_names__.emplace_back(std::string() + "c2");
    if (emit_transformed_parameters__) {
      
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // constrained_param_names() 
    
  inline void unconstrained_param_names(
                                        std::vector<std::string>& param_names__,
                                        bool emit_transformed_parameters__ = true,
                                        bool emit_generated_quantities__ = true) const
    final {
    
    for (int sym1__ = 1; sym1__ <= N_VARIANTS; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "q" + '.' + std::to_string(sym1__));
      } 
    }
    param_names__.emplace_back(std::string() + "c2");
    if (emit_transformed_parameters__) {
      
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // unconstrained_param_names() 
    
  inline std::string get_constrained_sizedtypes() const {
    
    return std::string("[{\"name\":\"q\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(N_VARIANTS) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"parameters\"},{\"name\":\"c2\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"}]");
    
    } // get_constrained_sizedtypes() 
    
  inline std::string get_unconstrained_sizedtypes() const {
    
    return std::string("[{\"name\":\"q\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(N_VARIANTS) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"parameters\"},{\"name\":\"c2\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"}]");
    
    } // get_unconstrained_sizedtypes() 
    
  
    // Begin method overload boilerplate
    template <typename RNG>
    inline void write_array(RNG& base_rng,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                            const bool emit_transformed_parameters = true,
                            const bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      std::vector<double> vars_vec;
      vars_vec.reserve(vars.size());
      std::vector<int> params_i;
      write_array_impl(base_rng, params_r, params_i, vars_vec,
          emit_transformed_parameters, emit_generated_quantities, pstream);
      vars = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>>(
        vars_vec.data(), vars_vec.size());
    }

    template <typename RNG>
    inline void write_array(RNG& base_rng, std::vector<double>& params_r,
                            std::vector<int>& params_i,
                            std::vector<double>& vars,
                            bool emit_transformed_parameters = true,
                            bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      write_array_impl(base_rng, params_r, params_i, vars,
       emit_transformed_parameters, emit_generated_quantities, pstream);
    }

    template <bool propto__, bool jacobian__, typename T_>
    inline T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
                       std::ostream* pstream = nullptr) const {
      Eigen::Matrix<int, -1, 1> params_i;
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }

    template <bool propto__, bool jacobian__, typename T__>
    inline T__ log_prob(std::vector<T__>& params_r,
                        std::vector<int>& params_i,
                        std::ostream* pstream = nullptr) const {
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }


    inline void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream = nullptr) const final {
      std::vector<double> params_r_vec;
      params_r_vec.reserve(params_r.size());
      std::vector<int> params_i;
      transform_inits_impl(context, params_i, params_r_vec, pstream);
      params_r = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>>(
        params_r_vec.data(), params_r_vec.size());
    }
    inline void transform_inits(const stan::io::var_context& context,
                                std::vector<int>& params_i,
                                std::vector<double>& vars,
                                std::ostream* pstream = nullptr) const final {
      transform_inits_impl(context, params_i, vars, pstream);
    }

};
}
using stan_model = qi_q_model_namespace::qi_q_model;

#ifndef USING_R

// Boilerplate
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}

stan::math::profile_map& get_stan_profile_data() {
  return qi_q_model_namespace::profiles__;
}

#endif

