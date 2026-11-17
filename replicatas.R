# =====================================================================
# Replicatas Monte Carlo do MTBD.SIR.SIRS
# Baseada na sua mc_MTBD_SIR_SIRS, com quatro acrescimos:
#   - aceita model = "BD" (o simulador ja aceitava, o wrapper nao)
#   - repassa aceitar_extincao
#   - limite de tempo por replicata (evita travar o app num BD explosivo)
#   - callback ao_progredir(r, re) para a barra de progresso do Shiny
# Continua funcionando fora do Shiny, igual a versao original.
# =====================================================================

# Avalia 'expr' com limite de tempo decorrido (segundos). Inf = sem limite.
com_limite <- function(expr, segundos = Inf) {
  if (is.finite(segundos) && segundos > 0) {
    setTimeLimit(elapsed = segundos, transient = TRUE)
    on.exit(setTimeLimit(elapsed = Inf, transient = FALSE), add = TRUE)
  }
  expr
}

mc_MTBD_SIR_SIRS <- function(re = 40, model = c("SIR", "SIRS", "BD"),
                             seed = NULL, verbose = TRUE,
                             n, lambdasky_1, lambdasky_2,
                             deathsky_1, deathsky_2,
                             timesky, sampprobsky, N,
                             rho = 0, timestop = 0,
                             omegasky = NULL, timeIntro2 = 0,
                             aceitar_extincao = FALSE,
                             limite_seg = Inf,
                             ao_progredir = NULL) {

  model <- match.arg(model)
  if (re < 1) stop("re deve ser >= 1")
  if (!is.null(seed) && !is.na(seed)) set.seed(seed)
  if (is.null(omegasky)) omegasky <- rep(0, length(timesky))
  if (model == "BD" && missing(N)) N <- 0

  resultados <- vector("list", re)
  tempos     <- numeric(re)
  status     <- character(re)
  logs       <- vector("list", re)

  tInicial <- Sys.time()

  for (r in seq_len(re)) {
    saida <- NULL
    log   <- character()
    erro  <- NULL

    t0 <- proc.time()[["elapsed"]]
    tryCatch({
      # capture.output guarda os print() do simulador (ex.: "extinct tree
      # at time ...") em vez de inundar o console
      log <- utils::capture.output(
        saida <- com_limite(
          MTBD.SIR.SIRS(
            n           = n,
            lambdasky_1 = lambdasky_1, lambdasky_2 = lambdasky_2,
            deathsky_1  = deathsky_1,  deathsky_2  = deathsky_2,
            timesky     = timesky,     sampprobsky = sampprobsky,
            rho         = rho,         timestop    = timestop,
            model       = model,       N           = N,
            omegasky    = omegasky,    timeIntro2  = timeIntro2,
            trackinfecteds   = TRUE,
            aceitar_extincao = aceitar_extincao
          ),
          limite_seg
        )
      )
    }, error = function(e) erro <<- conditionMessage(e))
    tempos[r] <- proc.time()[["elapsed"]] - t0

    # classifica o que aconteceu nesta replicata
    status[r] <-
      if (!is.null(erro) && is.finite(limite_seg) && tempos[r] >= limite_seg * 0.98) {
        "tempo esgotado"
      } else if (!is.null(erro)) {
        paste("erro:", erro)
      } else if (is.list(saida)) {
        "ok"
      } else if (any(grepl("nenhuma sub-arvore", log))) {
        "< 2 amostras"
      } else if (identical(saida, 0) || any(grepl("extinct tree", log))) {
        "extinta antes de timestop"
      } else {
        "parametros invalidos"
      }

    resultados[[r]] <- if (is.list(saida)) saida else NULL
    logs[[r]]       <- log

    if (verbose)
      cat(sprintf("[%s] rep %2d/%-2d  %5.2fs  %s\n", model, r, re, tempos[r], status[r]))
    if (is.function(ao_progredir)) ao_progredir(r, re)
  }

  elapsed_total <- as.numeric(difftime(Sys.time(), tInicial, units = "secs"))
  nFalhou <- sum(status != "ok")

  if (verbose) {
    cat(sprintf("\n===== %s | %d replicatas =====\n", model, re))
    cat(sprintf("tempo total:   %.2f s  (media %.3f s/rep)\n",
                elapsed_total, mean(tempos)))
    cat(sprintf("falhas:        %d\n\n", nFalhou))
  }

  invisible(list(
    model         = model,
    re            = re,
    elapsed_total = elapsed_total,
    tempos        = tempos,
    n_falhou      = nFalhou,
    resultados    = resultados,
    status        = status,
    logs          = logs,
    parametros    = list(n = n, N = N, rho = rho, timestop = timestop,
                         timeIntro2 = timeIntro2, timesky = timesky,
                         lambdasky_1 = lambdasky_1, lambdasky_2 = lambdasky_2,
                         deathsky_1 = deathsky_1, deathsky_2 = deathsky_2,
                         sampprobsky = sampprobsky, omegasky = omegasky,
                         aceitar_extincao = aceitar_extincao, seed = seed)
  ))
}
