if (!exists("fms")) {
  ########################################################
  ## test objects
  fms <- list(Dyestuff = lmer(Yield ~ (1 | Batch), Dyestuff),
              Dyestuff2 = lmer(Yield ~ (1 | Batch), Dyestuff2)
              , Penicillin = lmer(diameter ~ (1 | plate) + (1 | sample), Penicillin)
              , sleepstudy = lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
              ## , InstEval = lmer(y ~ studage + lectage + service + dept + (1 | s) + (1 | d), InstEval)
              ## FIXME: make tests pass with those and write new tests comparing classical case
              # , DyestuffOffs = lmer(Yield ~ (1 | Batch), Dyestuff,
              #                         weights = 1:nrow(Dyestuff), offset = rep(1, nrow(Dyestuff)))
              # , Dyestuff2Off = lmer(Yield ~ (1 | Batch), Dyestuff2,
              #                         offset = rep(-1, nrow(Dyestuff2)))
              # , PenicillinOff = lmer(diameter ~ (1 | plate) + (1 | sample), Penicillin,
              #                          offset = as.numeric(Penicillin$plate))
              # , sleepstudyOff = lmer(Reaction ~ Days + (Days|Subject), sleepstudy,
              #                          offset = as.numeric(sleepstudy$Subject))
              # , DyestuffWt = lmer(Yield ~ (1 | Batch), Dyestuff,
              #                     weights = 1:nrow(Dyestuff))
              # , Dyestuff2Wt = lmer(Yield ~ (1 | Batch), Dyestuff2,
              #                      weights = as.numeric(Dyestuff2$Batch) / 6)
              # , PenicillinWt = lmer(diameter ~ (1 | plate) + (1 | sample), Penicillin,
              #                       weights = with(Penicillin, as.numeric(plate) * as.numeric(sample) / 100))
              # , sleepstudyWt = lmer(Reaction ~ Days + (Days|Subject), sleepstudy,
              #                       weights = (sleepstudy$Days + 1) / 10)
              # , DyestuffWtOffs = lmer(Yield ~ (1 | Batch), Dyestuff,
              #                     weights = 1:nrow(Dyestuff), offset = rep(1, nrow(Dyestuff)))
              # , Dyestuff2WtOff = lmer(Yield ~ (1 | Batch), Dyestuff2,
              #                      weights = as.numeric(Dyestuff2$Batch) / 6,
              #                      offset = rep(-1, nrow(Dyestuff2)))
              # , PenicillinWtOff = lmer(diameter ~ (1 | plate) + (1 | sample), Penicillin,
              #                       weights = with(Penicillin, as.numeric(plate) * as.numeric(sample) / 100),
              #                       offset = as.numeric(Penicillin$plate))
              # , sleepstudyWtOff = lmer(Reaction ~ Days + (Days|Subject), sleepstudy,
              #                       weights = (sleepstudy$Days + 1) / 10,
              #                       offset = as.numeric(sleepstudy$Subject))
  )

  cat("Created source lme4 objects\n")

  ## FIXME should use the same functions as in the package
  rPDs <- lapply(fms, convToRlmerPredD, smoothPsi)
  rPDTests <- lapply(rPDs, convToRlmerPredD_test)
  rPDASs <- lapply(fms, convToRlmerPredD_DAS, smoothPsi)
  rPDASTests <- lapply(rPDs, convToRlmerPredD_DAS_test)
  rPDAStaus <- lapply(fms, convToRlmerPredD_DAStau, smoothPsi)
  resps <- lapply(fms, convToRlmerResp)
  objs <- sapply(seq_along(resps), conv2rlmerObj_test)
  names(objs) <- names(fms)

  cat("Created test objects\n")
}
