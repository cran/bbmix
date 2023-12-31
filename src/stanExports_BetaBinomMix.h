// Generated by rstantools.  Do not edit by hand.

/*
    bbmix is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    bbmix is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with bbmix.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.18.0
#include <stan/model/model_header.hpp>
namespace model_BetaBinomMix_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_BetaBinomMix");
    reader.add_event(44, 42, "end", "model_BetaBinomMix");
    return reader;
}
#include <stan_meta_header.hpp>
class model_BetaBinomMix : public prob_grad {
private:
    int N;
    int K;
    vector<int> n;
    vector<int> m;
    vector_d alpha_p;
    vector_d beta_p;
public:
    model_BetaBinomMix(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_BetaBinomMix(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_BetaBinomMix_namespace::model_BetaBinomMix";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        // initialize member variables
        try {
            current_statement_begin__ = 5;
            context__.validate_dims("data initialization", "N", "int", context__.to_vec());
            N = int(0);
            vals_i__ = context__.vals_i("N");
            pos__ = 0;
            N = vals_i__[pos__++];
            current_statement_begin__ = 6;
            context__.validate_dims("data initialization", "K", "int", context__.to_vec());
            K = int(0);
            vals_i__ = context__.vals_i("K");
            pos__ = 0;
            K = vals_i__[pos__++];
            current_statement_begin__ = 7;
            validate_non_negative_index("n", "N", N);
            context__.validate_dims("data initialization", "n", "int", context__.to_vec(N));
            validate_non_negative_index("n", "N", N);
            n = std::vector<int>(N,int(0));
            vals_i__ = context__.vals_i("n");
            pos__ = 0;
            size_t n_limit_0__ = N;
            for (size_t i_0__ = 0; i_0__ < n_limit_0__; ++i_0__) {
                n[i_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 8;
            validate_non_negative_index("m", "N", N);
            context__.validate_dims("data initialization", "m", "int", context__.to_vec(N));
            validate_non_negative_index("m", "N", N);
            m = std::vector<int>(N,int(0));
            vals_i__ = context__.vals_i("m");
            pos__ = 0;
            size_t m_limit_0__ = N;
            for (size_t i_0__ = 0; i_0__ < m_limit_0__; ++i_0__) {
                m[i_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 9;
            validate_non_negative_index("alpha_p", "K", K);
            context__.validate_dims("data initialization", "alpha_p", "vector_d", context__.to_vec(K));
            validate_non_negative_index("alpha_p", "K", K);
            alpha_p = vector_d(static_cast<Eigen::VectorXd::Index>(K));
            vals_r__ = context__.vals_r("alpha_p");
            pos__ = 0;
            size_t alpha_p_i_vec_lim__ = K;
            for (size_t i_vec__ = 0; i_vec__ < alpha_p_i_vec_lim__; ++i_vec__) {
                alpha_p[i_vec__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 10;
            validate_non_negative_index("beta_p", "K", K);
            context__.validate_dims("data initialization", "beta_p", "vector_d", context__.to_vec(K));
            validate_non_negative_index("beta_p", "K", K);
            beta_p = vector_d(static_cast<Eigen::VectorXd::Index>(K));
            vals_r__ = context__.vals_r("beta_p");
            pos__ = 0;
            size_t beta_p_i_vec_lim__ = K;
            for (size_t i_vec__ = 0; i_vec__ < beta_p_i_vec_lim__; ++i_vec__) {
                beta_p[i_vec__] = vals_r__[pos__++];
            }
            // validate, data variables
            current_statement_begin__ = 5;
            check_greater_or_equal(function__,"N",N,0);
            current_statement_begin__ = 6;
            check_greater_or_equal(function__,"K",K,0);
            current_statement_begin__ = 7;
            current_statement_begin__ = 8;
            current_statement_begin__ = 9;
            current_statement_begin__ = 10;
            // initialize data variables
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 14;
            validate_non_negative_index("theta", "K", K);
            num_params_r__ += (K - 1);
            current_statement_begin__ = 15;
            validate_non_negative_index("mu", "K", K);
            num_params_r__ += K;
            current_statement_begin__ = 16;
            validate_non_negative_index("lambda", "K", K);
            num_params_r__ += K;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_BetaBinomMix() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        stan::io::writer<double> writer__(params_r__,params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        if (!(context__.contains_r("theta")))
            throw std::runtime_error("variable theta missing");
        vals_r__ = context__.vals_r("theta");
        pos__ = 0U;
        validate_non_negative_index("theta", "K", K);
        context__.validate_dims("initialization", "theta", "vector_d", context__.to_vec(K));
        vector_d theta(static_cast<Eigen::VectorXd::Index>(K));
        for (int j1__ = 0U; j1__ < K; ++j1__)
            theta(j1__) = vals_r__[pos__++];
        try {
            writer__.simplex_unconstrain(theta);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable theta: ") + e.what());
        }
        if (!(context__.contains_r("mu")))
            throw std::runtime_error("variable mu missing");
        vals_r__ = context__.vals_r("mu");
        pos__ = 0U;
        validate_non_negative_index("mu", "K", K);
        context__.validate_dims("initialization", "mu", "vector_d", context__.to_vec(K));
        vector_d mu(static_cast<Eigen::VectorXd::Index>(K));
        for (int j1__ = 0U; j1__ < K; ++j1__)
            mu(j1__) = vals_r__[pos__++];
        try {
            writer__.ordered_unconstrain(mu);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable mu: ") + e.what());
        }
        if (!(context__.contains_r("lambda")))
            throw std::runtime_error("variable lambda missing");
        vals_r__ = context__.vals_r("lambda");
        pos__ = 0U;
        validate_non_negative_index("lambda", "K", K);
        context__.validate_dims("initialization", "lambda", "vector_d", context__.to_vec(K));
        vector_d lambda(static_cast<Eigen::VectorXd::Index>(K));
        for (int j1__ = 0U; j1__ < K; ++j1__)
            lambda(j1__) = vals_r__[pos__++];
        try {
            writer__.vector_lb_unconstrain(0.10000000000000001,lambda);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable lambda: ") + e.what());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(vector<T__>& params_r__,
                 vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            // model parameters
            stan::io::reader<local_scalar_t__> in__(params_r__,params_i__);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  theta;
            (void) theta;  // dummy to suppress unused var warning
            if (jacobian__)
                theta = in__.simplex_constrain(K,lp__);
            else
                theta = in__.simplex_constrain(K);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  mu;
            (void) mu;  // dummy to suppress unused var warning
            if (jacobian__)
                mu = in__.ordered_constrain(K,lp__);
            else
                mu = in__.ordered_constrain(K);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  lambda;
            (void) lambda;  // dummy to suppress unused var warning
            if (jacobian__)
                lambda = in__.vector_lb_constrain(0.10000000000000001,K,lp__);
            else
                lambda = in__.vector_lb_constrain(0.10000000000000001,K);
            // transformed parameters
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            // model body
            {
            current_statement_begin__ = 23;
            validate_non_negative_index("log_theta", "K", K);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  log_theta(static_cast<Eigen::VectorXd::Index>(K));
            (void) log_theta;  // dummy to suppress unused var warning
            stan::math::initialize(log_theta, DUMMY_VAR__);
            stan::math::fill(log_theta,DUMMY_VAR__);
            stan::math::assign(log_theta,stan::math::log(theta));
            current_statement_begin__ = 26;
            validate_non_negative_index("lps", "K", K);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  lps(static_cast<Eigen::VectorXd::Index>(K));
            (void) lps;  // dummy to suppress unused var warning
            stan::math::initialize(lps, DUMMY_VAR__);
            stan::math::fill(lps,DUMMY_VAR__);
            current_statement_begin__ = 28;
            for (int k = 1; k <= K; ++k) {
                current_statement_begin__ = 29;
                lp_accum__.add(beta_log<propto__>(get_base1(mu,k,"mu",1), get_base1(alpha_p,k,"alpha_p",1), get_base1(beta_p,k,"beta_p",1)));
                current_statement_begin__ = 30;
                lp_accum__.add(gamma_log<propto__>(get_base1(lambda,k,"lambda",1), (get_base1(alpha_p,k,"alpha_p",1) + get_base1(beta_p,k,"beta_p",1)), 1));
            }
            current_statement_begin__ = 33;
            lp_accum__.add(dirichlet_log<propto__>(theta, rep_vector(1.0,K)));
            current_statement_begin__ = 35;
            for (int i = 1; i <= N; ++i) {
                current_statement_begin__ = 36;
                stan::math::assign(lps, log_theta);
                current_statement_begin__ = 37;
                for (int k = 1; k <= K; ++k) {
                    current_statement_begin__ = 38;
                    stan::model::assign(lps, 
                                stan::model::cons_list(stan::model::index_uni(k), stan::model::nil_index_list()), 
                                (stan::model::rvalue(lps, stan::model::cons_list(stan::model::index_uni(k), stan::model::nil_index_list()), "lps") + beta_binomial_log(get_base1(n,i,"n",1),get_base1(m,i,"m",1),(get_base1(lambda,k,"lambda",1) * get_base1(mu,k,"mu",1)),(get_base1(lambda,k,"lambda",1) * (1 - get_base1(mu,k,"mu",1))))), 
                                "assigning variable lps");
                }
                current_statement_begin__ = 40;
                lp_accum__.add(log_sum_exp(lps));
            }
            }
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("theta");
        names__.push_back("mu");
        names__.push_back("lambda");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(K);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(K);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(K);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__,params_i__);
        static const char* function__ = "model_BetaBinomMix_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        vector_d theta = in__.simplex_constrain(K);
        vector_d mu = in__.ordered_constrain(K);
        vector_d lambda = in__.vector_lb_constrain(0.10000000000000001,K);
            for (int k_0__ = 0; k_0__ < K; ++k_0__) {
            vars__.push_back(theta[k_0__]);
            }
            for (int k_0__ = 0; k_0__ < K; ++k_0__) {
            vars__.push_back(mu[k_0__]);
            }
            for (int k_0__ = 0; k_0__ < K; ++k_0__) {
            vars__.push_back(lambda[k_0__]);
            }
        // declare and define transformed parameters
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // validate transformed parameters
            // write transformed parameters
            if (include_tparams__) {
            }
            if (!include_gqs__) return;
            // declare and define generated quantities
            // validate generated quantities
            // write generated quantities
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng,params_r_vec,params_i_vec,vars_vec,include_tparams,include_gqs,pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    static std::string model_name() {
        return "model_BetaBinomMix";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "theta" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "mu" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "lambda" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        for (int k_0__ = 1; k_0__ <= (K - 1); ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "theta" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "mu" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "lambda" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
    }
}; // model
}
typedef model_BetaBinomMix_namespace::model_BetaBinomMix stan_model;
#endif
