# jonashaslbeck@protonmail; Feb 15, 2023

# ------------------------------------------------------------
# -------- Function to Process mlVAR Outputs -----------------
# ------------------------------------------------------------


Process_mlVAR <- function(object1,
                          object2) {

  # a) Between network
  btw_1 <- object1$results$Gamma_Omega_mu$mean
  btw_1 <- (btw_1+t(btw_1)) / 2 # Apply AND-rule
  btw_2 <- object2$results$Gamma_Omega_mu$mean
  btw_2 <- (btw_2+t(btw_2)) / 2 # Apply AND-rule
  btw_diff <- btw_1 - btw_2


  # b.1) VAR: fixed effects
  phi_fix_1 <- object1$results$Beta$mean
  phi_fix_2 <- object2$results$Beta$mean
  phi_fix_diff <- phi_fix_1 - phi_fix_2


  # b.2) VAR: RE sds
  phi_RE_sd_1 <- object1$results$Beta$SD
  phi_RE_sd_2 <- object2$results$Beta$SD
  phi_RE_sd_diff <- phi_RE_sd_1 - phi_RE_sd_2


  # c.1) Contemp: fixed effects
  Gam_fix_1 <- object1$results$Gamma_Theta$mean
  Gam_fix_1 <- (Gam_fix_1 + t(Gam_fix_1)) / 2 # Apply AND-rule
  Gam_fix_2 <- object2$results$Gamma_Theta$mean
  Gam_fix_2 <- (Gam_fix_2 + t(Gam_fix_2)) / 2 # Apply AND-rule
  Gam_fix_diff <- Gam_fix_1 - Gam_fix_2


  # c.2) Contemp: RE sds
  Gam_RE_sd_1 <- object1$results$Gamma_Theta$SD
  Gam_RE_sd_1 <- (Gam_RE_sd_1 + t(Gam_RE_sd_1)) / 2 # Apply AND-rule
  Gam_RE_sd_2 <- object2$results$Gamma_Theta$SD
  Gam_RE_sd_2 <- (Gam_RE_sd_2 + t(Gam_RE_sd_2)) / 2 # Apply AND-rule
  Gam_RE_sd_diff <- Gam_RE_sd_1 - Gam_RE_sd_2


  outlist <- list("diff_between" = btw_diff,
                  "diff_phi_fix" = phi_fix_diff,
                  "diff_phi_RE_sd" = phi_RE_sd_diff,
                  "diff_gam_fix" = Gam_fix_diff,
                  "diff_gam_RE_sd" = Gam_RE_sd_diff)

  return(outlist)

} # eoF


# ------------------------------------------------------------
# -------- Function for Permutation Test ---------------------
# ------------------------------------------------------------

# Inputs:
# - two nested datasets
# - number of samples in permutation test

# Output:
# - sampling distribution for each parameter
# - test statistic based on the mlVAR models estimated on the two datasets

# vars <- c("V1", "V2", "V3")
# idvar <- "id"
# nP <- 5
# saveModels = TRUE
# verbose = TRUE


mlVAR_GC <- function(data1, # dataset of group 1
                     data2, # dataset of group 1
                     vars, # variables to be included in mlVAR (same in both data sets)
                     idvar, # variable indicating the nesting/subject id (same in both data sets)
                     dayvar = NULL,
                     beepvar = NULL,
                     estimator,
                     contemporaneous,
                     temporal,
                     nCores = 1,
                     nP = 500, # number of samples in permutation test
                     saveModels = FALSE, # if TRUE, all models are saved; defaults to FALSE to save memory
                     verbose = TRUE # if TRUE, progress bar is mapped on permutations
) {


  # ------ Input Checks -----

  # TODO:
  # Making sure that ids are unique across both datasets
  # ... do this by adding unique char strings in both groups

  # TODO:
  # All the usual input checks ...

  # ------ Collect passed down arguments -----
  if(missing(estimator)) estimator <- "default"
  if(missing(contemporaneous)) contemporaneous <- "orthogonal"
  if(missing(temporal)) temporal <- "orthogonal"
  if(missing(nCores)) nCores <- 1


  # ------ Get Basic Info -----

  # number of variables
  p <- length(vars)

  # Combine in list/single matrix for shorter code below
  l_data <- list(data1, data2)
  m_data_cmb <- rbind(data1, data2)

  # Get IDs
  ids1 <- sapply(data1[, idvar], as.character)
  ids2 <- sapply(data2[, idvar], as.character)
  v_ids <- c(ids1, ids2)
  v_u_ids <- unique(v_ids)
  u_ids1 <- unique(ids1)
  u_ids2 <- unique(ids2)
  # Get Number of subjects
  v_Ns <- c(length(u_ids1), length(u_ids2))
  totalN <- sum(v_Ns)


  # ------ Loop Over Permutations -----

  # Storage
  if(saveModels) l_out_mods <- list(vector("list", length = nP),
                                    vector("list", length = nP))

  # Setup cores, if multi-core
  if(nCores > 1) {
    cl <- makeCluster(nCores, outfile = "")
    registerDoParallel(cl)
  }

  # Progress bar
  # if(nCores==1)
  # if(verbose == TRUE) pb <- txtProgressBar(min = 0, max = nP + 1, initial = 0, char = "-", style = 3) # plus 1, because we also estimate on the true group split below

  out_P <- foreach(b = 1:nP,
                   .packages = c("mlVAR", "mixnet"),
                   .export = c("m_data_cmb", "vars", "idvar", "estimator",
                               "contemporaneous", "temporal"),
                   .verbose = TRUE) %dopar% {

                     # --- Make permutation ---
                     # This is done in a way that keeps the size in each group exactly the same as in the real groups
                     v_ids_rnd <- v_u_ids[sample(1:totalN, size=totalN, replace=FALSE)]
                     v_ids_1 <- v_ids_rnd[1:v_Ns[1]]
                     v_ids_2 <- v_ids_rnd[(v_Ns[1]+1):totalN]

                     # Split data based on permutations
                     data_h0_1 <- m_data_cmb[v_ids %in% v_ids_1, ]
                     data_h0_2 <- m_data_cmb[v_ids %in% v_ids_2, ]
                     l_data_h0 <- list(data_h0_1, data_h0_2)

                     # --- Fit mlVAR models ---

                     # Output list
                     l_pair_b <- list()
                     l_models <- list()

                     for(j in 1:2) {

                       # browser()

                       # TODO: make this variable specification of dayvar/beepvar less hacky
                       if(is.null(dayvar)) {
                         l_pair_b[[j]] <- mlVAR(data = l_data_h0[[j]],
                                                vars = vars,
                                                idvar = idvar,
                                                estimator = estimator,
                                                contemporaneous = contemporaneous,
                                                temporal = temporal,
                                                nCores = 1, # we use parallelization across resamples (not nodes in nodewise estimation here)
                                                verbose = FALSE,
                                                lags = 1) # TODO: later allow also higher order lags (see also below)
                       } else {
                         l_pair_b[[j]] <- mlVAR(data = l_data_h0[[j]],
                                                vars = vars,
                                                idvar = idvar,
                                                estimator = estimator,
                                                contemporaneous = contemporaneous,
                                                temporal = temporal,
                                                nCores = 1,
                                                dayvar = dayvar,
                                                beepvar = beepvar,
                                                verbose = FALSE,
                                                lags = 1) # TODO: later allow also higher order lags (see also below)
                       } # end if: dayvar specified

                       if(saveModels) l_models[[j]] <- l_pair_b[[j]]


                     } # end loop: J=2 groups

                     # All differences are: Group 1 - Group 2
                     diffs_b <- Process_mlVAR(object1 = l_pair_b[[1]],
                                              object2 = l_pair_b[[2]])

                     outlist_b <- list("diff_between" = diffs_b$diff_between,
                                       "diff_phi_fix" = diffs_b$diff_phi_fix,
                                       "diff_phi_RE_sd" = diffs_b$diff_phi_RE_sd,
                                       "diff_gam_fix" = diffs_b$diff_gam_fix,
                                       "diff_gam_RE_sd" = diffs_b$diff_between,
                                       "Models"=l_models)

                     return(outlist_b)

                     # Update progress bar
                     # if(nCores==1)
                     # if(verbose == TRUE) setTxtProgressBar(pb, b) # no idea what this will do in parallel processing

                   } # end foreach: over permutations

  # Close down cores, if multi-core
  if(nCores>1) stopCluster(cl)

  # ------ Loop results into objects for Sampling Distribution -----

  # Collect sampling distributions (for now) for:
  # a) between-person partial correlations
  # b.1) VAR/phi fixed effects
  # b.2) VAR/phi random effects sds
  # c.1) Contemp./Gamma fixed effects
  # c.2) Contemp./Gamma random effects sds
  # - I guess mean estimates make no sense, since we centered within-person
  # TODO: Later: Also output sampling distributions for correlations between REs, if specified
  # Create Storage
  a_between <- array(NA, dim=c(p, p, nP))
  a_phi_fixed <- array(NA, dim=c(p, p, nP))
  a_phi_RE_sd <- array(NA, dim=c(p, p, nP))
  a_gam_fixed <- array(NA, dim=c(p, p, nP))
  a_gam_RE_sd <- array(NA, dim=c(p, p, nP))

  # browser()

  for(b in 1:nP) {

    # browser()

    # Fill into arrays
    a_between[, , b] <- out_P[[b]]$diff_between
    a_phi_fixed[, , b] <- out_P[[b]]$diff_phi_fix[, , 1] # TODO: adapt also to higher order lags
    a_phi_RE_sd[, , b] <- out_P[[b]]$diff_phi_RE_sd[, , 1] # TODO: adapt also to higher order lags
    a_gam_fixed[, , b] <- out_P[[b]]$diff_gam_fix
    a_gam_RE_sd[, , b] <- out_P[[b]]$diff_gam_RE_sd

    if(saveModels) l_out_mods[[b]] <- diffs_b$Models

  } # end for: loop in permutations


  # ------ Create Test Statistics -----

  l_out_emp <- list()

  for(j in 1:2) {

    # TODO: make this variable specification of dayvar/beepvar less hacky
    if(is.null(dayvar)) {
      l_out_emp[[j]] <-  mlVAR(data = l_data[[j]],
                               vars = vars,
                               idvar = idvar,
                               estimator = estimator,
                               contemporaneous = contemporaneous,
                               temporal = temporal,
                               nCores = 1,
                               verbose = FALSE,
                               lags = 1)
    } else {
      l_out_emp[[j]] <-  mlVAR(data = l_data[[j]],
                               vars = vars,
                               idvar = idvar,
                               estimator = estimator,
                               contemporaneous = contemporaneous,
                               temporal = temporal,
                               nCores = 1,
                               dayvar = dayvar,
                               beepvar = beepvar,
                               verbose = FALSE,
                               lags = 1)
    } # end if: dayvar specified

    # if(nCores==1)
    # if(verbose == TRUE) setTxtProgressBar(pb, b+1)

  } # Loop: 2 groups


  # --- Matrices with True differences ---

  diffs_true <- Process_mlVAR(object1 = l_out_emp[[1]],
                              object2 = l_out_emp[[2]])


  # ------ Compute p-values -----

  # a) between
  m_pval_btw <- matrix(NA, p, p)
  for(i in 2:p) for(j in 1:(i-1)) m_pval_btw[i,j] <- mean(abs(a_between[i,j,])>abs(diffs_true$diff_between[i,j]))

  # b.1) VAR: fixed effects
  m_pval_phi_fix <- matrix(NA, p, p)
  for(i in 1:p) for(j in 1:p) m_pval_phi_fix[i,j] <- mean(abs(a_phi_fixed[i,j,])>abs(diffs_true$diff_phi_fix[i,j,]))

  # b.2) VAR: RE sds
  m_pval_phi_RE_sd <- matrix(NA, p, p)
  for(i in 1:p) for(j in 1:p) m_pval_phi_RE_sd[i,j] <- mean(abs(a_phi_RE_sd[i,j,])>abs(diffs_true$diff_phi_RE_sd[i,j,]))

  # c.1) Contemp: fixed effects
  m_pval_gam_fixed <- matrix(NA, p, p)
  for(i in 2:p) for(j in 1:(i-1)) m_pval_gam_fixed[i,j] <- mean(abs(a_gam_fixed[i,j,])>abs(diffs_true$diff_gam_fix[i,j]))

  # c.2) Contemp: RE sds
  m_pval_gam_RE_sd <- matrix(NA, p, p)
  for(i in 2:p) for(j in 1:(i-1)) m_pval_gam_RE_sd[i,j] <- mean(abs(a_gam_RE_sd[i,j,])>abs(diffs_true$diff_gam_RE_sd[i,j]))


  # ------ Create Output List -----

  if(saveModels) l_out_ret <- l_out else l_out_ret <- NULL

  outlist <- list("TrueDiffs" = list("Between" = diffs_true$diff_between,
                                     "Phi_mean" = diffs_true$diff_phi_fix,
                                     "Phi_sd" = diffs_true$diff_phi_RE_sd,
                                     "Gam_mean" = diffs_true$diff_gam_fix,
                                     "Gam_sd" = diffs_true$diff_gam_RE_sd),
                  "Pval" = list("Between" = m_pval_btw,
                                "Phi_mean" = m_pval_phi_fix,
                                "Phi_sd" = m_pval_phi_RE_sd,
                                "Gam_mean" = m_pval_gam_fixed,
                                "Gam_sd" = m_pval_gam_RE_sd),
                  "SampDist" = list("Between" = a_between,
                                    "Phi_mean" = a_phi_fixed,
                                    "Phi_sd" = a_phi_RE_sd,
                                    "Gam_mean" = a_gam_fixed,
                                    "Gam_sd" = a_gam_RE_sd),
                  "Models" = l_out_ret)

  # ------ Return Output -----

  return(outlist)


} # eoF





