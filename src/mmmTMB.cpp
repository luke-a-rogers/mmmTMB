#include <TMB.hpp>

// ---------------- Methods called by the model --------------------------------

template <class Type>
array<Type> create_aK(array<Type> a, matrix<int> m) {
  // Initialize values
  int na = a.dim(0);
  int npt = a.dim(1);
  int ng = a.dim(2);
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
          if (m(ca, pa) != 0) {
            d_sum += exp(a(cp_den, cpt, mg));
            cp_den +=1;
          }
        }
        // Set the probability sum to zero
        k_sum = Type(0);
        // Compute the probability for ca != pa
        for (int ca = 0; ca < na; ca++) {
          if (m(ca, pa) != 0) {
            k(ca, pa, cpt, mg) = exp(a(cp_num, cpt, mg)) / (Type(1) + d_sum);
            k_sum += k(ca, pa, cpt, mg);
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

  DATA_IMATRIX(mT); // Tag release integer matrix
  DATA_IMATRIX(mR); // Tag recover integer matrix
  DATA_IMATRIX(mI); // Area index integer matrix
  DATA_MATRIX(mL); // Fishing mortality rate bias matrix
  DATA_MATRIX(mW); // Fishing mortality rate weights matrix
  DATA_SCALAR(sA); // Proportion of tags still attached following release
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

  PARAMETER_ARRAY(aP); // Movement parameters dim = c(npt, np, ng)
  PARAMETER_ARRAY(lmF); // Log fishing mortality rates dim = c(nft, nfa)
  PARAMETER_VECTOR(lvB); // Log fishing bias length = ng
  PARAMETER(lsM); // Log natural mortality
  PARAMETER(lsH); // Log tag loss rate
  PARAMETER(lsC); // Log initial tag loss rate
  PARAMETER(lsD); // Log negative binomial dispersion

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
  matrix<Type> mF(nfa, nft); // Fishing mortality rates
  vector<Type> vB(ng); // Fishing mortality bias

  // -------------- Initialize arrays ------------------------------------------

  aR.setZero(); // Initialize to zero
  aN.setZero(); // Initialize to zero
  aK.setZero(); // Initialize to zero
  aS.setZero(); // Initialize to zero
  aRhat.setZero(); // Initialize to zero

  // -------------- Inverse transform parameters -------------------------------

  mF = exp(lmF); // Fishing mortality rates
  vB = exp(lvB); // Fishing bias
  Type sM = exp(lsM); // Natural mortality
  Type sH = exp(lsH); // Tag loss rate
  Type sD = exp(lsD); // Negative binomial dispersion

  // -------------- Compute constants ------------------------------------------

  int mTrows = mT.rows();
  int mRrows = mR.rows();

  // -------------- Populate tag abundances with releases ----------------------

  for (int i = 0; i < mTrows; i++) {
    aN(mT(1, i), mT(2, i), mT(0, i), mT(1, i), mT(2, i)) = sA * mT(3, i);
  }

  // -------------- Populate tag recoveries ------------------------------------

  for (int i = 0; i < mRrows; i++) {
    aR(mR(0, i), mR(1, i), mR(2, i), mR(3, i), mR(4, i)) = mR(5, i);
  }

  // -------------- Populate movement rates ------------------------------------

  aK = create_aK(aP, mI);

  // -------------- Compute the random effects nll -----------------------------



  // -------------- Populate the survival array --------------------------------

  for (int mg = 0; mg < ng; mg++) {
    for (int ct = 0; ct < nt; ct++) {
      for (int ca = 0; ca < na; ca++) {
        aS(ca, ct, mg) = exp(-(vB(mg) * mF(vfa(ca), vft(ct)) *
          mW(vwa(ca), vwt(ct)) + sM + sH));
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
                    (Type(1) - exp(-vB(mg) * mF(vfa(ra), vft(rt)) *
                    mW(vwa(ra), vwt(rt)))) *
                    mL(vla(ra), vlt(rt));
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
                        log(aRhat(ra, rt, mg, ma, mt)) + lsD,
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

  // REPORT(aP);
  // REPORT(aK);

  // ------------- Report simulation -------------------------------------------

  // SIMULATE {REPORT(aR);}

  // ------------- Return the nll ----------------------------------------------

  return jnll;
}
