data {

  // Flag for sampling from the prior alone
  int only_prior;

  // Dataset related data
  int<lower=0> N_datasets;

  // Antigen related data
  int<lower=0> N_ags;

  // Sera related data
  int<lower=0> N_srs;
  int<lower=0> N_sr_groups;
  int sr_groups[N_srs];

  // Titer-measurement related data
  real logtiters[N_ags, N_srs, N_datasets];
  int titertypes[N_ags, N_srs, N_datasets];

}

parameters {

  // Sera group related parameters
  matrix[N_ags, N_sr_groups] sr_group_gmts;

  // Sera related parameters
  vector[N_srs] sr_effects;

  // dataset related parameters
  vector[N_datasets] dataset_effects;
  array[N_ags, N_sr_groups, N_datasets] real sr_group_gmts_deviation;
  vector<lower=0.01>[N_datasets] sr_group_gmts_deviation_sd;

  // Allow different noise related parameters for each dataset
  vector[N_datasets] logtiter_error_sigma;

}

model {

  // Work out your prior likelihoods for titer noise
  for (dataset in 1:N_datasets) {
    logtiter_error_sigma[dataset] ~ inv_gamma(3, 1.5);
  }

  // Set prior for the sr_group_gmts
  for (ag in 1:N_ags) {
    sr_group_gmts[ag] ~ normal(7.0, 20.0);
  }

  // Set prior for the sr_effects
  sr_effects ~ normal(0.0, 6.0);

  // Set prior for the dataset_effects
  dataset_effects ~ normal(0.0, 6.0);

  // Set prior for the sr_group_gmts_deviation_sd
  for (dataset in 1:N_datasets) {
    sr_group_gmts_deviation_sd[dataset] ~ inv_gamma(2, 5);
  }

  // Inform the sr_group_gmts_deviation_sd
  for (dataset in 1:N_datasets) {
    for (ag in 1:N_ags) {
      for (sr_group in 1:N_sr_groups) {
        target += normal_lpdf(
          sr_group_gmts_deviation[ag, sr_group, dataset] | 0, sr_group_gmts_deviation_sd[dataset]
        );
      }
    }
  }

  // Cycle through each antigen and sera to get the predicted titer
  // and likelihood of the measurements given that prediction
  if (only_prior == 0) {
    for (dataset in 1:N_datasets) {
      for (sr in 1:N_srs) {
        for (ag in 1:N_ags) {

          // Only calculate predicted titer if there is a measured titer
          if (titertypes[ag, sr, dataset] > 0) {

            // Work out predicted log titer
            int sr_group = sr_groups[sr];
            real predicted_logtiter = sr_group_gmts[ag, sr_group] + sr_group_gmts_deviation[ag, sr_group, dataset] + dataset_effects[dataset] + sr_effects[sr];

            // Work out likelihood of the predicted titer
            if (titertypes[ag, sr, dataset] == 1) {
              target += normal_lpdf(
                logtiters[ag, sr, dataset] | predicted_logtiter, logtiter_error_sigma[dataset]
              );
            } else if (titertypes[ag, sr, dataset] == 2) {
              target += normal_lcdf(
                logtiters[ag, sr, dataset] | predicted_logtiter, logtiter_error_sigma[dataset]
              );
            }
          }
        }
      }
    }
  }
}
