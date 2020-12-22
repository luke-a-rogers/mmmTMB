#include <TMB.hpp>

// Methods called by the model -------------------------------------------------

template <class Type>
array<Type> create_movement_rates(array<Type> tp, matrix<int> tz) {
  // Initialize values
  int na = tz.rows(); // Matrix tz has rows = na, and cols = na
  int npt = tp.dim(1); // Array tp has dim = c(np, npt, ng)
  int ng = tp.dim(2);
  int cp_num; // Current parameter in the numerator
  int cp_den; // Current parameter in the denominator
  Type rates_sum; // Movement rate sum
  Type denom_sum; // Denominator sum
  // Initialize array
  array<Type> r(na, na, npt, ng);
  r.setZero();
  // Compute movement probabilities from movement parameters
  for (int mg = 0; mg < ng; mg++) {
    for (int cpt = 0; cpt < npt; cpt++) {
      // Set the current parameter indexes to zero
      cp_num = 0;
      cp_den = 0;
      for (int pa = 0; pa < na; pa++) {
        // Set the denominator sum to zero for this row
        denom_sum = 0;
        for (int ca = 0; ca < na; ca++) {
          if (tz(ca, pa) != 0) {
            denom_sum += exp(tp(cp_den, cpt, mg));
            cp_den +=1;
          }
        }
        // Set the movement rate sum to zero
        rates_sum = Type(0);
        // Compute the movement rate for ca != pa
        for (int ca = 0; ca < na; ca++) {
          if (tz(ca, pa) != 0) {
            r(pa, ca, cpt, mg) = exp(tp(cp_num, cpt, mg))/(Type(1) + denom_sum);
            rates_sum += r(pa, ca, cpt, mg);
            cp_num += 1;
          }
        }
        // Compute the movement rate for pa == ca
        r(pa, pa, cpt, mg) = Type(1) - rates_sum;
      }
    }
  }
  return r;
}

// Enumerate valid constants ---------------------------------------------------

enum valid_family {
  poisson_family = 0,
  nbinom1_family = 1,
  nbinom2_family = 2
};

enum valid_process {
  no_process = 0,
  rw_process = 1
};

// Model definition ------------------------------------------------------------

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data ----------------------------------------------------------------------

  DATA_IMATRIX(tx); // Tag release integer matrix (transpose)
  DATA_IMATRIX(ty); // Tag recover integer matrix (transpose)
  DATA_IMATRIX(tz); // Area index integer matrix (transpose)
  DATA_ARRAY(tl); // Tag reporting rate matrix (transpose)
  DATA_SCALAR(h); // Instantaneous tag loss rate
  DATA_SCALAR(d); // Tags attached after release (proportion)
  // Settings
  DATA_INTEGER(error_family); // Error family
  DATA_INTEGER(max_liberty); // Max steps at liberty
  // Index limits
  DATA_INTEGER(np); // Number of movement parameters per step and group
  DATA_INTEGER(nt); // Number of time steps
  DATA_INTEGER(na); // Number of areas
  DATA_INTEGER(ng); // Number of groups
  DATA_INTEGER(npt); // Number of parameter time steps
  DATA_INTEGER(nft); // Number of fishing rate time steps
  DATA_INTEGER(nfa); // Number of fishing rate areas
  DATA_INTEGER(nfg); // Number of fishing rate areas
  DATA_INTEGER(nlt); // Number of tag reporting rate time steps
  DATA_INTEGER(nla); // Number of tag reporting rate areas
  DATA_INTEGER(nlg); // Number of tag reporting rate areas
  // Index vectors
  DATA_IVECTOR(vpt); // Time step for movement parameters
  DATA_IVECTOR(vft); // Time step for fishing rate
  DATA_IVECTOR(vfa); // Area for fishing rate
  DATA_IVECTOR(vfg); // Area for fishing rate
  DATA_IVECTOR(vlt); // Time step for tag reporting rate
  DATA_IVECTOR(vla); // Area for tag reporting rate
  DATA_IVECTOR(vlg); // Area for tag reporting rate

  // Parameters ----------------------------------------------------------------

  PARAMETER_ARRAY(tp); // Movement parameters dim = c(np, npt, ng)
  PARAMETER_ARRAY(logit_exp_neg_tf); // Fishing mortality rate
  PARAMETER(logit_exp_neg_m); // Natural mortality
  PARAMETER(log_k); // Negative binomial dispersion
  // PARAMETER(log_f_sd);

  // Initialize the nll --------------------------------------------------------

  parallel_accumulator<Type> jnll(this);
  // Type jnll = 0;

  // Instantiate arrays --------------------------------------------------------

  // Matrices and arrays use column-major indexing in TMB (first index fastest)
  array<Type> Y(na, max_liberty, ng, na, nt); // Tag recoveries
  array<Type> N(na, max_liberty, ng, na, nt); // Tag abundance
  array<Type> R(na, na, npt, ng); // Movement rates
  array<Type> S(na, nt, ng); // Survival rates
  array<Type> Yhat(na, max_liberty, ng, na, nt); // Tag recoveries (predicted)
  array<Type> tF(nfa, nft, nfg); // Fishing mortality rates
  array<Type> tL(nla, nlt, nlg); // Tag reporting rates

  // Initialize arrays ---------------------------------------------------------

  Y.setZero(); // Initialize to zero
  N.setZero(); // Initialize to zero
  R.setZero(); // Initialize to zero
  S.setZero(); // Initialize to zero
  Yhat.setZero(); // Initialize to zero
  tL = tl;

  // Inverse transform parameters ----------------------------------------------

  // Fishing mortality rates
  for (int mg = 0; mg < ng; mg++) {
    for (int ct = 0; ct < nt; ct++) {
      for (int ca = 0; ca < na; ca++) {
        tF(vfa(ca), vft(ct), vfg(mg)) =
          -log(invlogit(logit_exp_neg_tf(vfa(ca), vft(ct), vfg(mg))));
          // TODO incorp F_dev in here at parameter level
      }
    }
  }
  // Natural mortality
  Type M = -log(invlogit(logit_exp_neg_m));
  // Fishing mortality bias standard deviation
  // Type F_sd = exp(log_f_sd);

  // Compute constants ---------------------------------------------------------

  int txcols = tx.cols();
  int tycols = ty.cols();

  // Populate tag abundances with releases -------------------------------------

  for (int i = 0; i < txcols; i++) {
    N(tx(1, i), 0, tx(2, i), tx(1, i), tx(0, i)) = d * tx(3, i);
  }

  // Populate tag recoveries ---------------------------------------------------

  for (int i = 0; i < tycols; i++) {
    Y(ty(4, i), ty(3, i) - ty(0, i), ty(2, i), ty(1, i), ty(0, i)) = ty(5, i);
  }

  // Populate movement rates ---------------------------------------------------

  R = create_movement_rates(tp, tz);

  // Compute the random effects nll --------------------------------------------

  // vector<Type> F_mu(nfa);
  // F_mu.setZero();
  // for (int ca = 0; ca < na; ca++) {
  //   F_mu(vfa(ca)) = tF.row(vfa(ca)).mean();
  // }
  // for (int ct = 0; ct < nt; ct++) {
  //   for (int ca = 0; ca < na; ca++) {
  //     jnll += dnorm(tF(vfa(ca), vft(ct)), F_mu(vfa(ca)), F_sd, TRUE);
  //   }
  // }

  // Populate the survival array -----------------------------------------------

  for (int mg = 0; mg < ng; mg++) {
    for (int ct = 0; ct < nt; ct++) {
      for (int ca = 0; ca < na; ca++) {
        S(ca, ct, mg) = exp(-tF(vfa(ca), vft(ct), vfg(mg)) - M - h);
      }
    }
  }

  // -------------- Predict recoveries -----------------------------------------

  // Uses ct = mt + dt
  // ct: current time step
  // mt: release time step
  // dt: time step difference
  // Initialize dt_max
  Type dt_max = Type(0);
  for (int mt = 0; mt < nt; mt++) {
    for (int ma = 0; ma < na; ma++) {
      for (int mg = 0; mg < ng; mg++) {
        // Were tags released in this stratum?
        if (N(ma, 0, mg, ma, mt) > 0) {
          // Define dt_max
          dt_max = max_liberty;
          if (mt + dt_max > nt) dt_max = nt - mt;
          // Project the abundance array forward
          for (int dt = 1; dt < dt_max; dt++) {
            for (int ca = 0; ca < na; ca++) {
              for (int pa = 0; pa < na; pa++) {
                N(ca, dt, mg, ma, mt) += N(pa, dt - 1, mg, ma, mt) *
                  S(pa, mt + dt - 1, mg) * R(pa, ca, vpt(mt + dt), mg);
              }
            }
          }
          // Compute predicted recoveries and jnll
          for (int dt = 1; dt < dt_max; dt++) {
            for (int ra = 0; ra < na; ra++) {
              // Populate the predicted recovery array
              Yhat(ra, dt, mg, ma, mt) = N(ra, dt, mg, ma, mt) *
                (Type(1) - exp(-tF(vfa(ra), vft(mt + dt), vfg(mg)))) *
                tL(vla(ra), vlt(mt + dt), vlg(mg));
              // Shelter error distribution from mean zero
              if (Yhat(ra, dt, mg, ma, mt) > 0) {
                // Calculate the joint negative log likelihood
                switch (error_family) {
                case poisson_family:
                  jnll -= dpois(
                    Y(ra, dt, mg, ma, mt),
                    Yhat(ra, dt, mg, ma, mt),
                    true);
                  break;
                case nbinom1_family:
                  jnll -= dnbinom_robust(
                    Y(ra, dt, mg, ma, mt),
                    log(Yhat(ra, dt, mg, ma, mt)),
                    log(Yhat(ra, dt, mg, ma, mt)) - log_k,
                    true);
                  break;
                case nbinom2_family:
                  jnll -= dnbinom_robust(
                    Y(ra, dt, mg, ma, mt),
                    log(Yhat(ra, dt, mg, ma, mt)),
                    log(Yhat(ra, dt, mg, ma, mt)) * Type(2) - log_k,
                    true);
                  break;
                default:
                  error("error_family not implemented");
                }
              }
            }
          }
        }
      }
    }
  }

  // Return the nll ------------------------------------------------------------

  return jnll;
}
