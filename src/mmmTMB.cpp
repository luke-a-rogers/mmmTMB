#include <TMB.hpp>


// ---------------- Enumerate valid constants ----------------------------------

enum valid_family {
  poisson_family = 0,
  nbinom1_family = 1,
  nbinom2_family = 2
};

enum valid_process {
  no_process = 0,
  rw_process = 1,
  ar_process = 2
};

template<class Type>
Type objective_function<Type>::operator() ()
{
  // -------------- Empirical counts -------------------------------------------

  DATA_ARRAY(released_3d);
  DATA_ARRAY(recovered_5d);

  // -------------- Estimates by area and time  --------------------------------

  DATA_MATRIX(capture_rate_2d);
  DATA_MATRIX(report_ratio_2d);

  // -------------- Global estimates--------------------------------------------

  DATA_SCALAR(tag_loss_rate);
  DATA_SCALAR(imm_loss_ratio);

  // -------------- Movement template ------------------------------------------

  DATA_IMATRIX(template_2d);

  // -------------- Model structure --------------------------------------------

  DATA_INTEGER(recapture_delay);
  DATA_INTEGER(error_family);
  DATA_INTEGER(result_units);
  DATA_INTEGER(time_process);
  // DATA_INTEGER(normalize);
  DATA_IVECTOR(movement_index);

  // -------------- Parameters -------------------------------------------------

  PARAMETER_ARRAY(movement_parameters_3d);
  PARAMETER_ARRAY(log_capture_bias_2d);
  PARAMETER(log_natural_mortality);
  PARAMETER(log_dispersion);

  // -------------- Initialize the nll -----------------------------------------

  parallel_accumulator<Type> jnll(this);
  // Type jnll = 0;
  Type re_nll = 0;

  // -------------- Initialize array dimensions --------------------------------

  int np = template_2d.sum(); // Number of parameters for one group at one time
  int nt = released_3d.dim(0); // Number of time steps in release and recovery data
  int na = released_3d.dim(1); // Number of release and recovery areas
  int ng = released_3d.dim(2); // Number of release and recovery groups
  int nv = movement_parameters_3d.dim(0); // Number of movement parameter time steps

  // -------------- Initialize arrays ------------------------------------------

  array<Type> capture_bias_2d(na, ng);
  array<Type> survival_3d(nt, na, ng);
  array<Type> movement_probability_4d(na, na, nv, ng);
  array<Type> predicted_abundance_5d(nt, nt, na, na, ng);
  array<Type> predicted_occupancy_5d(nt, nt, na, na, ng);
  array<Type> predicted_recovered_5d(nt, nt, na, na, ng);
  // Initialize to zero
  predicted_abundance_5d.setZero();
  predicted_occupancy_5d.setZero();
  predicted_recovered_5d.setZero();

  // -------------- Inverse transform: Scalar parameters -----------------------

  Type natural_mortality = exp(log_natural_mortality);
  Type dispersion = exp(log_dispersion);

  // -------------- Inverse transform: Bias parameters -------------------------

  for (int ca = 0; ca < na; ca++) {
    for (int mg = 0; mg < ng; mg++) {
      capture_bias_2d(ca, mg) = exp(log_capture_bias_2d(ca, mg));
    }
  }

  // -------------- Populate the movement probability array --------------------

  // Initialize  values
  int par_numer;
  int par_denom;
  Type denom_sum;
  Type prob_sum;
  // Compute probabilities
  for (int mg = 0; mg < ng; mg++) {
    for (int cv = 0; cv < nv; cv++) {
      par_numer = 0;
      par_denom = 0;
      for (int pa = 0; pa < na; pa++) {
        // Set the denom sum to 0 for this row
        denom_sum = 0;
        // Compute the denom sum
        for (int ca = 0; ca < na; ca++) {
          if (template_2d(pa, ca) > 0) {
            denom_sum += exp(movement_parameters_3d(cv, par_denom, mg));
            par_denom += 1;
          }
        }
        // Set the prob sum to 0 for this row
        prob_sum = Type(0);
        // Compute the probability for pa != ca
        for (int ca = 0; ca < na; ca++) {
          if (template_2d(pa, ca) > 0) {
            movement_probability_4d(pa, ca, cv, mg) =
              exp(movement_parameters_3d(cv, par_numer, mg)) / (Type(1) + denom_sum);
            prob_sum += movement_probability_4d(pa, ca, cv, mg);
            par_numer += 1;
          }
        }
        // Compute the probability for pa == ca
        movement_probability_4d(pa, pa, cv, mg) = Type(1) - prob_sum;
      }
    }
  }

  // -------------- Compute the random effects nll -----------------------------



  // -------------- Populate the survival array --------------------------------

  for (int ct = 0; ct < nt; ct++) {
    for (int ca = 0; ca < na; ca++) {
      for (int mg = 0; mg < ng; mg++) {
        survival_3d(ct, ca, mg) =
          exp(-(capture_bias_2d(ca, mg) * capture_rate_2d(ct, ca) +
          natural_mortality + tag_loss_rate));
      }
    }
  }

  // -------------- Predict recoveries -----------------------------------------
  // TODO: Add simulation blocks
  for (int mt = 0; mt < nt; mt++) {
    for (int ma = 0; ma < na; ma++) {
      for (int mg = 0; mg < ng; mg++) {
        // Were tags released in this stratum?
        if (released_3d(mt, ma, mg) > 0) {
          // Populate the abundance array with releases
          predicted_abundance_5d(mt, mt, ma, ma, mg) = (Type(1) - imm_loss_ratio) *
            released_3d(mt, ma, mg);
          // Project the abundance array forward
          for (int ct = (mt + 1); ct < nt; ct++) {
            for (int ca = 0; ca < na; ca++) {
              for (int pa = 0; pa < na; pa++) {
                predicted_abundance_5d(mt, ct, ma, ca, mg) +=
                  predicted_abundance_5d(mt, ct - 1, ma, pa, mg) *
                  survival_3d(ct - 1, pa, mg) *
                  movement_probability_4d(pa, ca, movement_index(ct), mg);
              }
              // Record predicted occupancy
              if (predicted_abundance_5d(mt, ct, ma, ca, mg) > 0) {
                predicted_occupancy_5d(mt, ct, ma, ca, mg) = Type(1);
              }
            }
          }
          // Compute predicted recoveries and jnll
          for (int rt = mt; rt < nt; rt++) {
            for (int ra = 0; ra < na; ra++) {
              // Populate the predicted recovery array
              predicted_recovered_5d(mt, rt, ma, ra, mg) =
                predicted_abundance_5d(mt, rt, ma, ra, mg) *
                (1 - exp(-capture_bias_2d(ra, mg) *
                capture_rate_2d(rt, ra))) * report_ratio_2d(rt, ra);
              // Was occupancy predicted?
              // Quick fix: replace predicted_occupancy_5d by predicted_recovered_5d
              if (predicted_recovered_5d(mt, rt, ma, ra, mg) > 0) {
                // Was the recapture time later than the delay?
                if (rt > (mt + recapture_delay - Type(1))) {
                  // Calculate the negative log likelihood
                  switch (error_family) {
                  case poisson_family:
                    jnll += -dpois(
                      recovered_5d(mt, rt, ma, ra, mg),
                      predicted_recovered_5d(mt, rt, ma, ra, mg),
                      true);
                    // SIMULATE {recovered_5d(mt, rt, ma, ra, mg) =
                    // 	rpois(predicted_recovered_5d(mt, rt, ma, ra, mg));}
                    break;
                  case nbinom1_family:
                    jnll += -dnbinom_robust(
                      recovered_5d(mt, rt, ma, ra, mg),
                      log(predicted_recovered_5d(mt, rt, ma, ra, mg)),
                      log(predicted_recovered_5d(mt, rt, ma, ra, mg)) - log_dispersion,
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

  // ------------- Transform to results units ----------------------------------

  Type natural_mortality_results = natural_mortality * result_units;

  // ------------- AD Report ---------------------------------------------------

  ADREPORT(natural_mortality_results);
  ADREPORT(dispersion);
  ADREPORT(capture_bias_2d);

  // ------------- Report ------------------------------------------------------

  // REPORT(movement_parameters_3d);
  // REPORT(movement_probability_4d);

  // ------------- Report simulation -------------------------------------------

  // SIMULATE {REPORT(recovered_5d);}

  // ------------- Return the nll ----------------------------------------------

  // if (normalize) {return re_nll;}
  return jnll;
}
