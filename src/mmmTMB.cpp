#include <TMB.hpp>

// ---------------- Methods called by the model --------------------------------

template <class Type>
array<Type> create_aK(array<Type> ta, matrix<int> tm) {
  // Initialize values
  int na = tm.rows(); // Matrix tm has rows = na, and cols = na
  int npt = ta.dim(1); //  Array ta has dim = c(np, npt, ng)
  int ng = ta.dim(2);
  int cp_num; // Current parameter in the numerator (?)
  int cp_den; // Current parameter in the denominator (?)
  Type k_sum; // Probability sum
  Type d_sum; // Denominator sum
  // Initialize array
  array<Type> k(na, na, npt, ng);
  k.setZero();
  // Compute movement probabilities from movement parameters
  for (int mg = 0; mg < ng; mg++) {
    for (int cpt = 0; cpt < npt; cpt++) {
      // Set the current parameter indexes to zero
      cp_num = 0;
      cp_den = 0;
      for (int pa = 0; pa < na; pa++) {
        // Set the denominator sum to zero for this row
        d_sum = 0;
        for (int ca = 0; ca < na; ca++) {
          if (tm(ca, pa) != 0) {
            d_sum += exp(ta(cp_den, cpt, mg));
            cp_den +=1;
          }
        }
        // Set the probability sum to zero
        k_sum = Type(0);
        // Compute the probability for ca != pa
        for (int ca = 0; ca < na; ca++) {
          if (tm(ca, pa) != 0) {
            k(pa, ca, cpt, mg) = exp(ta(cp_num, cpt, mg)) / (Type(1) + d_sum);
            k_sum += k(pa, ca, cpt, mg);
            cp_num += 1;
          }
        }
        // Compute the probability for pa == ca
        k(pa, pa, cpt, mg) = Type(1) - k_sum;
      }
    }
  }
  return k;
}

// ---------------- Enumerate valid constants ----------------------------------

enum valid_family {
  poisson_family = 0,
  nbinom1_family = 1
};

enum valid_process {
  no_process = 0,
  rw_process = 1
};

// ---------------- Model definition -------------------------------------------

template<class Type>
Type objective_function<Type>::operator() ()
{
  // -------------- Data -------------------------------------------------------

  DATA_IMATRIX(tmT); // Tag release integer matrix (transpose)
  DATA_IMATRIX(tmR); // Tag recover integer matrix (transpose)
  DATA_IMATRIX(tmI); // Area index integer matrix (transpose)
  DATA_MATRIX(tmL); // Fishing mortality rate bias matrix (transpose)
  DATA_MATRIX(tmW); // Fishing mortality rate weights matrix (transpose)
  // Settings
  DATA_INTEGER(error_family); // Error family
  DATA_IVECTOR(span_liberty); // Min and max steps at liberty
  DATA_INTEGER(time_varying); // Time varying?
  DATA_INTEGER(time_process); // None or random walk
  // Index limits
  DATA_INTEGER(np); // Number of movement parameters per step and group
  DATA_INTEGER(nt); // Number of time steps
  DATA_INTEGER(na); // Number of areas
  DATA_INTEGER(ng); // Number of groups
  DATA_INTEGER(npt); // Number of parameter time steps
  DATA_INTEGER(nft); // Number of fishing rate time steps
  DATA_INTEGER(nfa); // Number of fishing rate areas
  DATA_INTEGER(nlt); // Number of tag reporting rate time steps
  DATA_INTEGER(nla); // Number of tag reporting rate areas
  DATA_INTEGER(nwt); // Number of fishing rate weight time steps
  DATA_INTEGER(nwa); // Number of fishing rate weight areas
  // Index vectors
  DATA_IVECTOR(vpt); // Time step for movement parameters
  DATA_IVECTOR(vft); // Time step for fishing rate
  DATA_IVECTOR(vfa); // Area for fishing rate
  DATA_IVECTOR(vlt); // Time step for tag reporting rate
  DATA_IVECTOR(vla); // Area for tag reporting rate
  DATA_IVECTOR(vwt); // Time step for fishing rate weight
  DATA_IVECTOR(vwa); // Area for fishing rate weight

  // -------------- Parameters -------------------------------------------------

  PARAMETER_ARRAY(taP); // Movement parameters dim = c(np, npt, ng)
  PARAMETER_ARRAY(logit_exp_neg_tmF); // Fishing mortality
  PARAMETER(logit_exp_neg_sM); // Natural mortality
  PARAMETER(logit_exp_neg_sH); // Tag loss
  PARAMETER(logit_sA); // 1 - Initial tag loss
  PARAMETER_VECTOR(log_vB); // Fishing bias by group
  PARAMETER(log_sD); // Negative binomial dispersion

  // -------------- Initialize the nll -----------------------------------------

  parallel_accumulator<Type> jnll(this);
  // Type jnll = 0;
  Type re_nll = 0;

  // -------------- Instantiate arrays -----------------------------------------

  // Arrays (3D+) use column-major indexing in TMB (first index fastest)
  array<Type> aR(na, nt, ng, na, nt); // Tag recoveries
  array<Type> aN(na, nt, ng, na, nt); // Tag abundance (predicted)
  array<Type> aK(na, na, npt, ng); // Movement rates
  array<Type> aS(na, nt, ng); // Survival rates
  array<Type> aRhat(na, nt, ng, na, nt); // Tag recoveries (predicted)
  matrix<Type> tmF(nfa, nft); // Fishing mortality rates
  vector<Type> vB(ng); // Fishing mortality bias

  // -------------- Initialize arrays ------------------------------------------

  aR.setZero(); // Initialize to zero
  aN.setZero(); // Initialize to zero
  aK.setZero(); // Initialize to zero
  aS.setZero(); // Initialize to zero
  aRhat.setZero(); // Initialize to zero

  // -------------- Inverse transform parameters -------------------------------

  // Fishing mortality
  for (int ct = 0; ct < nt; ct++) {
    for (int ca = 0; ca < na; ca++) {
      tmF(vfa(ca), vft(ct)) = -log(invlogit(logit_exp_neg_tmF(vfa(ca), vft(ct))));
    }
  }
  Type sM = -log(invlogit(logit_exp_neg_sM)); // Natural mortality
  Type sH = -log(invlogit(logit_exp_neg_sH)); // Tag loss rate
  Type sA = invlogit(logit_sA); // Tags attached after release (proportion)
  vB = exp(log_vB); // Fishing bias
  Type sD = exp(log_sD); // Negative binomial dispersion

  // -------------- Compute constants ------------------------------------------

  int tmTcols = tmT.cols();
  int tmRcols = tmR.cols();

  // -------------- Populate tag abundances with releases ----------------------

  for (int i = 0; i < tmTcols; i++) {
    aN(tmT(1, i), tmT(0, i), tmT(2, i), tmT(1, i), tmT(0, i)) = sA * tmT(3, i);
  }

  // -------------- Populate tag recoveries ------------------------------------

  for (int i = 0; i < tmRcols; i++) {
    aR(tmR(4, i), tmR(3, i), tmR(2, i), tmR(1, i), tmR(0, i)) = tmR(5, i);
  }

  // -------------- Populate movement rates ------------------------------------

  aK = create_aK(taP, tmI);

  // -------------- Compute the random effects nll -----------------------------



  // -------------- Populate the survival array --------------------------------

  for (int mg = 0; mg < ng; mg++) {
    for (int ct = 0; ct < nt; ct++) {
      for (int ca = 0; ca < na; ca++) {
        aS(ca, ct, mg) = exp(-(vB(mg) * tmF(vfa(ca), vft(ct)) *
          tmW(vwa(ca), vwt(ct)) + sM + sH));
      }
    }
  }

  // -------------- Predict recoveries -----------------------------------------

  // TODO: Add simulation blocks
  for (int mt = 0; mt < nt; mt++) {
    for (int ma = 0; ma < na; ma++) {
      for (int mg = 0; mg < ng; mg++) {
        // Were tags released in this stratum?
        if (aN(ma, mt, mg, ma, mt) > 0) {
          // Project the abundance array forward
          for (int ct = (mt + 1); ct < nt; ct++) {
            for (int ca = 0; ca < na; ca++) {
              for (int pa = 0; pa < na; pa++) {
                aN(ca, ct, mg, ma, mt) += aN(pa, ct - 1, mg, ma, mt) *
                  aS(pa, ct - 1, mg) * aK(pa, ca, vpt(ct), mg);
              }
            }
          }
          // Compute predicted recoveries and jnll
          for (int rt = mt; rt < nt; rt++) {
            // Was the recapture time within the span at liberty?
            if (rt > (mt + span_liberty(1) - Type(1))) {
              if (rt < (mt + span_liberty(2) + Type(1))) {
                for (int ra = 0; ra < na; ra++) {
                  // Populate the predicted recovery array
                  aRhat(ra, rt, mg, ma, mt) = aN(ra, rt, mg, ma, mt) *
                    (Type(1) - exp(-vB(mg) * tmF(vfa(ra), vft(rt)) *
                    tmW(vwa(ra), vwt(rt)))) *
                    tmL(vla(ra), vlt(rt));
                  // Shelter error distribution from mean zero
                  if (aRhat(ra, rt, mg, ma, mt) > 0) {
                    // Calculate the joint negative log likelihood
                    switch (error_family) {
                    case poisson_family:
                      jnll += -dpois(
                        aR(ra, rt, mg, ma, mt),
                        aRhat(ra, rt, mg, ma, mt),
                        true);
                      // Simulate {aR(ra, rt, mg, ma, mt) =
                      // rpois(aRhat(ra, rt, mg, ma, mt));}
                      break;
                    case nbinom1_family:
                      jnll += -dnbinom_robust(
                        aR(ra, rt, mg, ma, mt),
                        log(aRhat(ra, rt, mg, ma, mt)),
                        log(aRhat(ra, rt, mg, ma, mt)) + log_sD,
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
    }
  }

  // ------------- AD Report ---------------------------------------------------

  ADREPORT(sD); // Negative binomial dispersion
  ADREPORT(vB); // Fishing mortality bias

  // ------------- Report ------------------------------------------------------

  // REPORT(taP);
  // REPORT(aK);

  // ------------- Report simulation -------------------------------------------

  // SIMULATE {REPORT(aR);}

  // ------------- Return the nll ----------------------------------------------

  return jnll;
}
