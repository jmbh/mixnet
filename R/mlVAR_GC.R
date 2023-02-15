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
# nB <- 5
# saveModels = TRUE
# verbose = TRUE

mlVAR_GC <- function(data1, # dataset of group 1
                     data2, # dataset of group 1
                     vars, # variables to be included in mlVAR (same in both data sets)
                     idvar, # variable indicating the nesting/subject id (same in both data sets)
                     nB = 500, # number of samples in permutation test
                     saveModels = TRUE, # if TRUE, all models are saved
                     verbose = TRUE, # if TRUE, progress bar is mapped on permutations
                     ... # arguments passed down to mlVAR()
) {


  # ------ Input Checks -----

  # TODO:
  # Making sure that ids are unique across both datasets
  # ... do this by adding unique char strings in both groups

  # TODO:
  # All the usual input checks ...


  # ------ Get Basic Info -----

  # number of variables
  p <- length(vars)

  # Combine in list/single matrix for shorter code below
  l_data <- list(data1, data2)
  m_data_cmb <- rbind(data1, data2)

  # Get IDs
  l_ids <- lapply(l_data, function(x) unique(x$id))
  v_ids <- unlist(l_ids)
  # Get Number of subjects
  N_ids <- lapply(l_ids, length)
  v_Ns <- unlist(N_ids)
  totalN <- sum(v_Ns)


  # ------ Loop Over Permutations -----

  # Storage
  l_out <- list(vector("list", length = nB),
                vector("list", length = nB))

  # Progress bar
  if(verbose == TRUE) pb <- txtProgressBar(min = 0, max=nB, initial = 0, char="-", style = 3)

  for(b in 1:nB) {

    # --- Make permutation ---
    # This is done in a way that keeps the size in each group exactly the same as in the real groups
    v_ids_rnd <- v_ids[sample(1:totalN, size=totalN, replace=FALSE)]
    v_ids_1 <- v_ids_rnd[1:v_Ns[1]]
    v_ids_2 <- v_ids_rnd[(v_Ns[1]+1):totalN]

    # Split data based on permutations
    data_h0_1 <- m_data_cmb[m_data_cmb$id %in% v_ids_1, ]
    data_h0_2 <- m_data_cmb[m_data_cmb$id %in% v_ids_2, ]
    l_data_h0 <- list(data_h0_1, data_h0_2)

    # --- Fit mlVAR models ---

    for(j in 1:2) {

      l_out[[j]][[b]] <- mlVAR(data = l_data_h0[[j]],
                               vars = vars,
                               idvar = idvar,
                               verbose = FALSE,
                               lags = 1) # TODO: later allow also higher order lags (see also below)

    } # end loop: J=2 groups

    # Update progress bar
    if(verbose==TRUE) setTxtProgressBar(pb, b)

  } # end loop: nB permutations


  # ------ Collect Sampling Distribution -----

  # Collect sampling distributions (for now) for:
  # a) between-person partial correlations
  # b.1) VAR/phi fixed effects
  # b.2) VAR/phi random effects sds
  # c.1) Contemp./Gamma fixed effects
  # c.2) Contemp./Gamma random effects sds
  # - I guess mean estimates make no sense, since we centered within-person
  # TODO: Later: Also output sampling distributions for correlations between REs, if specified

  # Create Storage
  a_between <- array(NA, dim=c(p, p, nB))
  a_phi_fixed <- array(NA, dim=c(p, p, nB))
  a_phi_RE_sd <- array(NA, dim=c(p, p, nB))
  a_gam_fixed <- array(NA, dim=c(p, p, nB))
  a_gam_RE_sd <- array(NA, dim=c(p, p, nB))

  for(b in 1:nB) {

    # All differences are: Group 1 - Group 2
    diffs_b <- Process_mlVAR(object1 = l_out[[1]][[b]],
                             object2 = l_out[[2]][[b]])

    # Fill into arrays
    a_between[, , b] <- diffs_b$diff_between
    a_phi_fixed[, , b] <- diffs_b$diff_phi_fix[, , 1] # TODO: adapt also to higher order lags
    a_phi_RE_sd[, , b] <- diffs_b$diff_phi_RE_sd[, , 1] # TODO: adapt also to higher order lags
    a_gam_fixed[, , b] <- diffs_b$diff_gam_fix
    a_gam_RE_sd[, , b] <- diffs_b$diff_gam_RE_sd

  } # for: permutations


  # ------ Create Test Statistics -----
  # Group 1
  out_emp_1 <- mlVAR(data = data1,
                     vars = vars,
                     idvar = idvar,
                     verbose = FALSE,
                     lags = 1)
  # Group 2
  out_emp_2 <- mlVAR(data = data2,
                     vars = vars,
                     idvar = idvar,
                     verbose = FALSE,
                     lags = 1)


  # --- Matrices with True differences ---

  diffs_true <- Process_mlVAR(object1 = out_emp_1,
                              object2 = out_emp_2)



  # ------ Compute p-values -----

  # a) between
  m_pval_btw <- matrix(NA, p, p)
  for(i in 2:p) for(j in 1:(i-1)) m_pval_btw[i,j] <- mean(abs(a_between[i,j,])>abs(diffs_true$diff_between[i,j]))

  # b.1) VAR: fixed effects
  m_pval_phi_fix <- matrix(NA, p, p)
  for(i in 1:p) for(j in 1:p) m_pval_phi_fix[i,j] <- mean(abs(a_phi_fixed[i,j,])>abs(diffs_true$diff_phi_fix[i,j,]))
  # DEVV:
  # a_phi_fixed[1,3,]
  # hist(a_phi_fixed[1,3,], xlim = c(-.5, .5))
  # abline(v = diffs_true$diff_phi_fix[1,3,], col="red")

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





