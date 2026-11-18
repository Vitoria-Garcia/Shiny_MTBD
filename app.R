# =====================================================================
# App Shiny: simulacao Monte Carlo MTBD (duas linhagens) SIR / SIRS / BD
#
# ARQUIVO UNICO: tudo (simulador, replicatas, graficos e app) esta aqui,
# entao o app roda de qualquer jeito:
#   - botao "Run App" do RStudio
#   - shiny::runApp("caminho/para/app.R")
#   - selecionar tudo (Ctrl+A) e rodar (Ctrl+Enter)
#
# Secoes:
#   1. Simulador MTBD.SIR.SIRS (sua funcao + teto de infectados, marcado [ALTERADO: teto])
#   2. Replicatas mc_MTBD_SIR_SIRS
#   3. Graficos e tabela-resumo
#   4. Auxiliares do app
#   5. Interface
#   6. Servidor
#
# Pacotes: shiny, bslib, ape, phytools  (+ showtext, opcional, p/ PDF em Montserrat)
# =====================================================================

library(shiny)
library(bslib)
library(ape)
library(phytools)


# #####################################################################
# 1. SIMULADOR (sua funcao; unica mudanca: teto de infectados,
#    procure por [ALTERADO: teto])
# #####################################################################

MTBD.SIR.SIRS <- function(n, lambdasky_1, lambdasky_2, deathsky_1, deathsky_2,
                          timesky, sampprobsky, rho = 0, timestop = 0,
                          model = "SIR", N = 0, trackinfecteds = FALSE,
                          omegasky = rep(0, length(timesky)),
                          timeIntro2 = 0,
                          aceitar_extincao = FALSE,     # [ALTERADO] novo argumento
                          max_infectados = Inf) {       # [ALTERADO: teto]
  # [ALTERADO] aceitar_extincao = TRUE: se os infectados zerarem antes de
  #   timestop (no SIR isso e o FIM NATURAL da epidemia, nao so extincao
  #   precoce), a simulacao termina e as arvores sao devolvidas. Quem chama
  #   decide o que descartar (ex.: < k amostras). Com FALSE, comportamento
  #   original: a replicata inteira e jogada fora (devolve 0).
  # [ALTERADO: teto] max_infectados: teto de infectados vivos (as duas
  #   linhagens somadas). Ao atingir o teto a simulacao para NAQUELE instante,
  #   como um timestop antecipado: as amostras feitas ate ali formam as
  #   arvores e, com rho > 0, os vivos sao amostrados nesse instante.
  #   Evita a explosao do BD. Inf = sem teto (comportamento original).
  
  detail <- 0
  
  # Define o modelo
  # BD: birth-death puro. Taxa de nascimento = lambda CONSTANTE (sem susc/N).
  #     Sem suscetiveis, sem recuperados, sem omega. Crescimento exponencial.
  #     As duas linhagens crescem de forma INDEPENDENTE (nao competem).
  if (model == "SIR") {
    SIRS <- FALSE; SIR <- TRUE;  BD <- FALSE
  } else if (model == "SIRS") {
    SIR <- FALSE;  SIRS <- TRUE; BD <- FALSE
  } else if (model == "BD") {
    SIR <- FALSE;  SIRS <- FALSE; BD <- TRUE
  } else {
    stop("Apenas os modelos 'SIR', 'SIRS' e 'BD' sao aceitos")
  }
  
  # Remocao dividida em morte (mu) e amostragem (psi), por linhagem
  # amostragem IGUAL para as duas
  psampsky <- sampprobsky
  
  # Linhagem A
  musky_1  <- deathsky_1 * (1 - psampsky)
  psisky_1 <- deathsky_1 * psampsky
  
  # Linhagem B
  musky_2  <- deathsky_2 * (1 - psampsky)
  psisky_2 <- deathsky_2 * psampsky
  
  extant <- 0
  
  # Condicoes de parada
  if (rho > 0 && timestop == 0) { print("rho > 0 and timestop = 0 is not a valid stopping condition")
    return()
  }
  
  if (n == 0 && timestop == 0) { print("n = 0 and timestop = 0 is not a valid stopping condition")
    return()
  }
  
  if (n > 0 && timestop > 0)   { print("n > 0 and timestop > 0 are two stopping condition. Condition n is ignored.")
  }
  
  # n > N so faz sentido no SIR (populacao finita); BD nao tem N
  if (SIR && n > N)            { print("n > N cannot be ran")
    return()
  }
  
  # omega so existe no SIRS
  if (SIR || BD) {
    for (i in seq_along(omegasky)) {
      if (omegasky[i] != 0) { print("omega is not a parameter for SIR/BD models")
        return() }
    }
  }
  
  if (timeIntro2 > 0 && timestop > 0 && timeIntro2 >= timestop) {
    print("timeIntro2 >= timestop: a Linhagem B nao chega a ser introduzida")
  }
  
  # AVISO BD: se lambda >> mu + psi e a parada e por timestop, a populacao
  # cresce exponencialmente e a simulacao pode ficar lenta / gigante.
  
  timesky <- c(timesky, 10^100)
  timesky <- timesky[-1]
  extincttree <- 1
  
  # Função que registra o estado atual
  record <- function() {
    timeevent  <<- c(timeevent, time)
    totalinfecteds <<- c(totalinfecteds, length(leaves))
    totalsusceptibles <<- c(totalsusceptibles, susc)
    totalrecovereds <<- c(totalrecovereds, recovered)
    totalsampleds <<- c(totalsampleds, length(sampled))
    cumulativeinfected <<- c(cumulativeinfected, cuminf)
    n1_hist <<- c(n1_hist, sum(lineages == 1))
    n2_hist <<- c(n2_hist, sum(lineages == 2))
    leaves_hist   <<- c(leaves_hist,   list(leaves))
    lineages_hist <<- c(lineages_hist, list(lineages))
  }
  
  
  #####################################################
  ########## Começa simulação #########################
  #####################################################
  
  while (extincttree == 1) {
    edge <- rbind(c(-1, -2), c(-3, -4))
    
    timeCreation4 <- max(timeIntro2, 0)
    timecreation <- c(0, 0, 0, timeCreation4)
    sampled <- vector()
    sampledtype <- vector()
    extinct <- vector()
    time <- 0
    maxspecies <- -4
    edge.length <- c(0, 0)
    extincttree <- 0
    stop <- 0
    recovered <- 0
    totalrecovereds <- 0
    timeevent <- 0
    leaves_hist   <- list()
    lineages_hist <- list()
    pl_tipo <- integer(); pl_p1 <- numeric()   # [ALTERADO] substitui prop_log <- data.frame()
    parou_teto <- FALSE                          # [ALTERADO: teto]
    
    # Quando a Linhagem B entra
    if (timeIntro2 <= 0) {
      leaves <- c(-2, -4)
      lineages <- c(1, 2)
      susc <- N - 2
      totalsusceptibles <- N - 2
      totalinfecteds <- 2
      cumulativeinfected <- 2
      cuminf <- 2
      n1_hist <- 1
      n2_hist <- 1
      totalsampleds <- 0
      introduzida2 <- TRUE
    } else {
      leaves <- c(-2)
      lineages <- c(1)
      susc <- N - 1
      totalsusceptibles <- N - 1
      totalinfecteds <- 1
      cumulativeinfected <- 1
      cuminf <- 1
      n1_hist <- 1
      n2_hist <- 0
      totalsampleds <- 0
      introduzida2 <- FALSE
    }
    
    # No BD nao existe pool de suscetiveis: susc e irrelevante (NA para
    # deixar explicito que nao deve ser interpretado como suscetiveis).
    if (BD) { susc <- NA; totalsusceptibles <- NA }
    
    lambda_1 <- lambdasky_1[1]
    mu_1 <- musky_1[1]
    psi_1 <- psisky_1[1]
    lambda_2 <- lambdasky_2[1]
    mu_2 <- musky_2[1]
    psi_2 <- psisky_2[1]
    omega  <- omegasky[1]
    timecut <- timesky[1]
    shift <- 2
    
    while (stop == 0) {
      if (length(leaves) == 0) {
        phy2 <- 0
        extincttree <- 1
        print(paste("extinct tree at time ", time, sep = ""))
        if (timestop > 0 && aceitar_extincao) {   # [ALTERADO]
          extincttree <- 0                          # [ALTERADO] mantem a arvore
        } else if (timestop > 0) { phy2 <- 0
        return(phy2) }
        stop <- 1
        
      } else {
        lam_i <- ifelse(lineages == 1, lambda_1, lambda_2)
        mu_i  <- ifelse(lineages == 1, mu_1, mu_2)
        psi_i <- ifelse(lineages == 1, psi_1, psi_2)
        
        # >>> taxa de nascimento por folha <<<
        #  BD : lambda CONSTANTE
        #  SIR/SIRS : lambda/N * susc (depende de suscetiveis)
        birth_i <- if (BD) lam_i else lam_i / N * susc
        rate_i  <- birth_i + mu_i + psi_i
        
        taxa_folhas  <- sum(rate_i)
        
        # perda de imunidade so no SIRS; no SIR/BD omega = 0 -> 0
        taxa_SIRS <- omega * recovered
        
        totalrate <- taxa_folhas + taxa_SIRS
        
        timestep <- rexp(1, totalrate)
        time <- time + timestep
        
        prox_intro <- if (!introduzida2) timeIntro2 else Inf
        
        if (timestop > 0 && time > timestop) {
          stop <- 1
        } else if (time >= prox_intro && prox_intro <= timecut) {
          # Introdução da Linhagem B
          time <- timeIntro2
          leaves <- c(leaves, -4)
          lineages <- c(lineages, 2)
          if (!BD) susc <- susc - 1     # BD nao consome suscetivel
          cuminf <- cuminf + 1
          introduzida2 <- TRUE
          record()
        } else if (time < timecut) {
          
          if (runif(1, 0, totalrate) < taxa_SIRS) {
            # Perda de imunidade (so SIRS)
            if (SIRS) { recovered <- recovered - 1
            susc <- susc + 1
            record() }
            
          } else {
            
            del  <- sample.int(length(leaves), 1, prob = rate_i)
            
            # [ALTERADO] antes: prop_log <<- rbind(...). O '<<-' gravava no
            # ambiente global (o prop_log devolvido saia sempre vazio) e o
            # rbind de data.frame a cada evento era o gargalo da simulacao.
            pl_tipo <- c(pl_tipo, lineages[del])
            pl_p1   <- c(pl_p1, sum(rate_i[lineages == 1]) / sum(rate_i))
            
            species <- leaves[del]
            tp <- lineages[del]
            lam <- lam_i[del]
            mu <- mu_i[del]
            psi <- psi_i[del]
            
            # >>> taxa da folha escolhida (mesma logica birth) <<<
            birth <- if (BD) lam else lam / N * susc
            taxa  <- birth + mu + psi
            edgespecevent <- which(edge == species) - length(edge.length)
            specevent <- runif(1, 0, 1)
            
            # Nascimento / Transmissao
            if (birth / taxa > specevent) {
              edge.length[edgespecevent] <- time - timecreation[-species]
              edge <- rbind(edge, c(species, maxspecies - 1))
              edge <- rbind(edge, c(species, maxspecies - 2))
              edge.length <- c(edge.length, 0, 0)
              leaves   <- c(leaves, maxspecies - 1, maxspecies - 2)
              lineages <- c(lineages, tp, tp)
              maxspecies <- maxspecies - 2
              leaves   <- leaves[-del]
              lineages <- lineages[-del]
              timecreation <- c(timecreation, time, time)
              if (!BD) susc <- susc - 1     # BD nao consome suscetivel
              cuminf <- cuminf + 1
              record()
              
              # Amostragem
            } else if ((birth + psi) / taxa > specevent) {
              sampled <- c(sampled, species)
              sampledtype <- c(sampledtype, tp)
              leaves <- leaves[-del]
              lineages <- lineages[-del]
              edge.length[edgespecevent] <- time - timecreation[-species]
              if (length(sampled) == n && timestop == 0) stop <- 1
              if (SIRS) recovered <- recovered + 1
              record()
              
              # Morte / Extincao sem amostragem
            } else {
              extinct <- c(extinct, species)
              leaves <- leaves[-del]
              lineages <- lineages[-del]
              edge.length[edgespecevent] <- time - timecreation[-species]
              if (SIRS) recovered <- recovered + 1
              record()
            }
          }
        } else {
          time <- timecut
          timecut <- timesky[shift]
          lambda_1 <- lambdasky_1[shift]
          mu_1 <- musky_1[shift]
          psi_1 <- psisky_1[shift]
          lambda_2 <- lambdasky_2[shift]
          mu_2 <- musky_2[shift]
          psi_2 <- psisky_2[shift]
          omega <- omegasky[shift]
          shift <- shift + 1
        }
        
        # [ALTERADO: teto] para quando os infectados vivos atingem o teto
        if (stop == 0 && length(leaves) >= max_infectados) {
          stop <- 1
          parou_teto <- TRUE
          print(paste("teto de ", max_infectados, " infectados atingido no tempo ", time, sep = ""))
        }
      }
    } # fim while (stop == 0)
  }   # fim while (extincttree == 1)
  
  # individuos vivos no tempo final
  if (rho == 0) {
    extant <- length(leaves)
    while (length(leaves) > 0) {
      del <- 1
      extinct <- c(extinct, leaves[del])
      k <- which(edge == leaves[del]) - length(edge.length)
      edge.length[k] <- time - timecreation[-leaves[del]]
      leaves <- leaves[-del]
      lineages <- lineages[-del]
    }
  } else {
    if (!parou_teto) time <- timestop   # [ALTERADO: teto] se parou no teto, amostra no instante da parada
    todiesample <- runif(length(leaves), 0, 1)
    todie <- which(todiesample > rho)
    extinct <- c(extinct, leaves[todie])
    if (length(todie) == length(leaves)) print("no sampled extant individuals")
    for (j in 1:length(leaves)) {
      k <- which(edge == leaves[j]) - length(edge.length)
      edge.length[k] <- time - timecreation[-leaves[j]]
    }
    if (length(todie) > 0) { leaves <- leaves[-todie]
    lineages <- lineages[-todie] }
    sampled     <- c(sampled, leaves)
    sampledtype <- c(sampledtype, lineages)
    recovered   <- recovered + length(leaves)
    leaves <- vector()
    lineages <- vector()
    record()
  }
  
  prop_log <- data.frame(tipo = pl_tipo, p1 = pl_p1)   # [ALTERADO]
  
  # =====================================================================
  # RECONSTRUCAO multiPhylo
  # =====================================================================
  montar_phylo <- function(edge, edge.length, raiz, sampled, extinct) {
    if (length(sampled) < 2) return(NULL)
    if (is.null(dim(edge))) edge <- matrix(edge, ncol = 2, byrow = TRUE)
    
    ext_sub <- extinct[extinct %in% edge[, 2]]
    if (length(ext_sub) > 0) {
      for (j in seq_along(ext_sub)) {
        del      <- which(edge == ext_sub[j]) - length(edge.length)
        surpress <- edge[del, 1]
        edge.length <- edge.length[-del]
        edge  <- edge[-del, , drop = FALSE]
        del2   <- which(edge[, 1] == surpress)
        modify <- which(edge[, 2] == surpress)
        edge[modify, 2] <- edge[del2, 2]
        edge.length[modify] <- edge.length[modify] + edge.length[del2]
        edge.length <- edge.length[-del2]
        edge <- edge[-del2, , drop = FALSE]
      }
    }
    
    S <- length(sampled)
    leaf <- 1
    interior <- S + 1
    temp <- sort(unique(c(edge)), decreasing = TRUE)
    for (j in temp) {
      if (sum(match(sampled, j, 0)) == 0 || j == raiz) {
        posvalue <- interior
        interior <- interior + 1
      } else {
        posvalue <- leaf
        leaf <- leaf + 1
      }
      replacel <- which(edge == j)
      for (k in seq_along(replacel)) {
        if (replacel[k] <= length(edge.length)) edge[replacel[k], 1] <- posvalue
        else edge[replacel[k] - length(edge.length), 2] <- posvalue
      }
    }
    
    phy <- list(edge = edge)
    phy$tip.label   <- paste0("t", sample(S))
    phy$edge.length <- edge.length
    phy$Nnode       <- S
    class(phy) <- "phylo"
    # [ALTERADO] o ramo-tronco (raiz -> MRCA das amostras) vira root.edge
    # direto do collapse.singles. Da o mesmo valor que o getx dava, mas
    # sem depender do TreeSim e sem o risco do sort() misturar tempos com
    # a coluna ttype (0/1) do getx em arvores com altura < 1.
    # Tempo absoluto do MRCA = root.time + root.edge.
    collapse.singles(phy, root.edge = TRUE)
  }
  
  descendentes <- function(raiz, edge) {
    fila <- raiz
    visitou <- raiz
    while (length(fila) > 0) {
      atual <- fila[1]
      fila <- fila[-1]
      filhos <- edge[edge[, 1] == atual, 2]
      novos  <- setdiff(filhos, visitou)
      visitou <- c(visitou, novos)
      fila <- c(fila, novos)
    }
    visitou
  }
  
  edgeM <- if (is.null(dim(edge))) matrix(edge, ncol = 2, byrow = TRUE) else edge
  
  arvores <- list()
  tipos_ok <- integer()
  roots   <- c(-1, -3)
  for (r in seq_along(roots)) {
    raiz <- roots[r]
    t0_sub <- if (raiz == -3) timeCreation4 else 0
    nos_sub <- descendentes(raiz, edgeM)
    linhas <- which(edgeM[, 1] %in% nos_sub)
    edge_sub <- edgeM[linhas, , drop = FALSE]
    el_sub <- edge.length[linhas]
    samp_sub <- sampled[sampled %in% nos_sub]
    phy <- montar_phylo(edge_sub, el_sub, raiz, samp_sub, extinct)
    if (!is.null(phy)) {
      phy$root.time <- t0_sub
      arvores[[length(arvores) + 1]] <- phy
      tipos_ok <- c(tipos_ok, r)
    }
  }
  
  if (length(arvores) == 0) {
    print("nenhuma sub-arvore com amostras suficientes")
    return(0)
  }
  class(arvores) <- "multiPhylo"
  names(arvores) <- paste0("tipo", tipos_ok)
  
  if (trackinfecteds) {
    list(arvores, list(
      eventtimes          = timeevent,
      infecteds           = totalinfecteds,
      cumulativeinfecteds = cumulativeinfected,
      cumulativesampleds  = totalsampleds,
      susceptibles        = totalsusceptibles,
      recovereds          = totalrecovereds,
      N_hist              = totalsusceptibles,
      time_hist           = timeevent,
      n1                  = n1_hist,
      n2                  = n2_hist,
      sampledtype         = sampledtype,
      parou_teto          = parou_teto,       # [ALTERADO: teto]
      prop_log = prop_log ))}
  else {
    arvores
  }
}

# #####################################################################
# 2. REPLICATAS
# #####################################################################

# =====================================================================
# Replicatas Monte Carlo do MTBD.SIR.SIRS
# Baseada na sua mc_MTBD_SIR_SIRS, com quatro acrescimos:
#   - aceita model = "BD" (o simulador ja aceitava, o wrapper nao)
#   - repassa aceitar_extincao
#   - limite de tempo por replicata (evita travar o app num BD explosivo)
#   - callback ao_progredir(r, re) para a barra de progresso do Shiny
#   - teto de infectados vivos (max_infectados): a replicata que atinge o
#     teto e truncada naquele instante ou, com descartar_teto = TRUE, descartada
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
                             max_infectados = Inf,
                             descartar_teto = FALSE,
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
            aceitar_extincao = aceitar_extincao,
            max_infectados   = max_infectados
          ),
          limite_seg
        )
      )
    }, error = function(e) erro <<- conditionMessage(e))
    tempos[r] <- proc.time()[["elapsed"]] - t0
    
    # classifica o que aconteceu nesta replicata
    teto <- is.list(saida) && isTRUE(saida[[2]]$parou_teto)
    status[r] <-
      if (!is.null(erro) && is.finite(limite_seg) && tempos[r] >= limite_seg * 0.98) {
        "tempo esgotado"
      } else if (!is.null(erro)) {
        paste("erro:", erro)
      } else if (teto && descartar_teto) {
        "descartada (teto)"
      } else if (teto) {
        "ok · parou no teto"
      } else if (is.list(saida)) {
        "ok"
      } else if (any(grepl("nenhuma sub-arvore", log))) {
        "< 2 amostras"
      } else if (identical(saida, 0) || any(grepl("extinct tree", log))) {
        "extinta antes de timestop"
      } else {
        "parametros invalidos"
      }
    
    # resultados[r] <- list(...) (e nao [[r]] <- NULL, que APAGA o elemento e
    # encurta a lista quando as ultimas replicatas falham)
    resultados[r] <- list(if (is.list(saida) && !(teto && descartar_teto)) saida else NULL)
    logs[[r]]       <- log
    
    if (verbose)
      cat(sprintf("[%s] rep %2d/%-2d  %5.2fs  %s\n", model, r, re, tempos[r], status[r]))
    if (is.function(ao_progredir)) ao_progredir(r, re)
  }
  
  elapsed_total <- as.numeric(difftime(Sys.time(), tInicial, units = "secs"))
  nFalhou <- sum(!vapply(resultados, is.list, logical(1)))
  
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
                         aceitar_extincao = aceitar_extincao,
                         max_infectados = max_infectados,
                         descartar_teto = descartar_teto, seed = seed)
  ))
}

# #####################################################################
# 3. GRAFICOS E TABELA-RESUMO
# #####################################################################

# =====================================================================
# Extracao dos resultados e graficos (base R), usados pelo app e pelos PDFs
# =====================================================================

COR_TRAJ <- c("#935999", "#073A70")   # trajetorias (cores do seu PDF)
COR_FENO <- c("#935999", "#073A70")   # arvores e figura da aba Sobre

# ---- extracao ---------------------------------------------------------

# Trajetorias (2o elemento de cada replicata bem-sucedida)
trajetorias_de <- function(resultados) {
  Filter(Negate(is.null),
         lapply(resultados, function(r) if (is.list(r)) r[[2]] else NULL))
}

# Tempo (absoluto) em que a arvore comeca a ser desenhada: o MRCA das
# amostras. A raiz da phylo e o MRCA; o ramo-tronco desde a origem da
# linhagem (root.time) ficou em root.edge (collapse.singles). Por isso o
# deslocamento e root.time + root.edge, e nao so root.time.
inicio_mrca <- function(tr) {
  t0 <- if (!is.null(tr$root.time)) tr$root.time else 0
  re <- if (!is.null(tr$root.edge)) tr$root.edge else 0
  t0 + re
}

# Todas as arvores validas, com a linhagem e a replicata de origem
extrai_arvores <- function(resultados, min_tips = 3) {
  arvores <- list(); tipos <- integer(); reps <- integer()
  for (r in seq_along(resultados)) {
    res <- resultados[[r]]
    if (is.null(res) || !is.list(res)) next
    mp <- res[[1]]
    if (!inherits(mp, "multiPhylo")) next
    nm <- names(mp)
    for (i in seq_along(mp)) {
      tr <- mp[[i]]
      if (!inherits(tr, "phylo") || ape::Ntip(tr) < min_tips) next
      arvores[[length(arvores) + 1]] <- tr
      tipos <- c(tipos, if (!is.null(nm) && nm[i] == "tipo2") 2L else 1L)
      reps  <- c(reps, r)
    }
  }
  list(arvores = arvores, tipos = tipos, reps = reps,
       inicio = vapply(arvores, inicio_mrca, numeric(1)))
}

# Introducao da Linhagem B a partir dos dados (como no seu script)
tempo_intro2 <- function(trajs) {
  if (!length(trajs)) return(NA_real_)
  v <- suppressWarnings(min(sapply(trajs, function(z) {
    w <- z$eventtimes[z$n2 > 0]; if (length(w)) min(w) else Inf
  })))
  if (is.finite(v)) v else NA_real_
}

# Ultimo instante com dados (trajetorias e pontas das arvores)
tempo_maximo <- function(trajs, info = NULL) {
  t1 <- if (length(trajs)) max(sapply(trajs, function(z) max(z$eventtimes))) else 0
  t2 <- 0
  if (!is.null(info) && length(info$arvores))
    t2 <- max(mapply(function(tr, o) o + max(phytools::nodeHeights(tr)),
                     info$arvores, info$inicio))
  max(t1, t2, 1e-8)
}

# Serie em grade regular (funcao-degrau); NA depois do fim da replicata
serie_em_grade <- function(z, campo, grade) {
  idx <- findInterval(grade, z$eventtimes)
  y <- rep(NA_real_, length(grade))
  ok <- idx > 0 & grade <= max(z$eventtimes)
  y[ok] <- z[[campo]][idx[ok]]
  y
}

mediana_em_grade <- function(trajs, campo, grade) {
  M <- sapply(trajs, serie_em_grade, campo = campo, grade = grade)
  if (is.null(dim(M))) M <- matrix(M, ncol = 1)
  apply(M, 1, function(l) if (all(is.na(l))) NA else stats::median(l, na.rm = TRUE))
}

# ---- graficos ---------------------------------------------------------

estilo_eixos <- list(col.axis = "grey45", col.lab = "grey45", fg = "grey60")

grafico_vazio <- function(msg, familia = "") {
  op <- par(mar = c(0, 0, 0, 0), family = familia); on.exit(par(op))
  plot.new(); text(0.5, 0.5, msg, col = "grey50", cex = 1.1)
}

linha_intro2 <- function(intro2) {
  if (!is.na(intro2) && intro2 > 0) abline(v = intro2, lty = 2, col = "grey65")
}

VARIAVEIS_TRAJ <- c("Infectados por linhagem" = "inf",
                    "Infectados (total)"      = "infecteds",
                    "Suscetíveis"             = "susceptibles",
                    "Recuperados (SIRS)"      = "recovereds",
                    "Amostrados acumulados"   = "cumulativesampleds")

plot_trajetorias <- function(trajs, variavel = "inf", cores = COR_TRAJ,
                             alpha = 0.4, intro2 = NA, xmax = NULL,
                             mediana = FALSE, familia = "") {
  if (!length(trajs)) return(grafico_vazio("Nenhuma trajetória disponível", familia))
  op <- par(mar = c(4.5, 4.5, 1.5, 1.5), family = familia,
            font.axis = 2, font.lab = 2)
  on.exit(par(op))
  if (is.null(xmax)) xmax <- tempo_maximo(trajs)
  grade <- seq(0, xmax, length.out = 600)
  
  if (variavel == "inf") {
    ymax <- max(sapply(trajs, function(z) max(c(z$n1, z$n2))))
    do.call(plot, c(list(NA, xlim = c(0, xmax), ylim = c(0, ymax), bty = "n",
                         xlab = "Tempo", ylab = "Infectados vivos"), estilo_eixos))
    linha_intro2(intro2)
    c1 <- adjustcolor(cores[1], alpha.f = alpha)
    c2 <- adjustcolor(cores[2], alpha.f = alpha)
    for (z in trajs) {
      lines(z$eventtimes, z$n1, type = "s", col = c1)
      i2 <- which(z$n2 > 0)
      if (length(i2)) {                       # da entrada ate o fim (inclui a queda a 0)
        sel <- seq(min(i2), length(z$n2))
        lines(z$eventtimes[sel], z$n2[sel], type = "s", col = c2)
      }
    }
    if (mediana) {
      for (k in 1:2) {
        m <- mediana_em_grade(trajs, c("n1", "n2")[k], grade)
        lines(grade, m, type = "s", lwd = 5, col = "white")
        lines(grade, m, type = "s", lwd = 2.5, col = cores[k])
      }
    }
    legend("topright", c("Linhagem A", "Linhagem B"), col = cores,
           text.col = cores, lwd = 2, bty = "n", text.font = 2)
  } else {
    vals <- unlist(lapply(trajs, `[[`, variavel))
    if (!length(vals) || all(is.na(vals)))
      return(grafico_vazio("Esta variável não existe neste modelo (ex.: S no BD, R no SIR)", familia))
    rotulo <- names(VARIAVEIS_TRAJ)[VARIAVEIS_TRAJ == variavel]
    do.call(plot, c(list(NA, xlim = c(0, xmax), ylim = c(0, max(vals, na.rm = TRUE)),
                         bty = "n", xlab = "Tempo", ylab = rotulo), estilo_eixos))
    linha_intro2(intro2)
    cc <- adjustcolor(cores[2], alpha.f = alpha)
    for (z in trajs) lines(z$eventtimes, z[[variavel]], type = "s", col = cc)
    if (mediana) {
      m <- mediana_em_grade(trajs, variavel, grade)
      lines(grade, m, type = "s", lwd = 5, col = "white")
      lines(grade, m, type = "s", lwd = 2.5, col = cores[2])
    }
  }
  invisible()
}

# Arvores de uma replicata, no eixo de tempo absoluto
plot_arvores_rep <- function(mp, cores = COR_FENO, xmax, intro2 = NA,
                             familia = "") {
  if (!inherits(mp, "multiPhylo") || !length(mp))
    return(grafico_vazio("Replicata sem árvores", familia))
  n <- length(mp)
  op <- par(mfrow = c(n, 1), mar = c(2, 1, 2.2, 1), oma = c(2.5, 0, 0, 0),
            family = familia, font.main = 2)
  on.exit(par(op))
  for (i in seq_len(n)) {
    tr   <- mp[[i]]
    tipo <- if (identical(names(mp)[i], "tipo2")) 2L else 1L
    t0   <- if (!is.null(tr$root.time)) tr$root.time else 0
    ape::plot.phylo(tr, root.edge = TRUE, show.tip.label = FALSE,
                    edge.color = cores[tipo], edge.width = 1.2,
                    x.lim = c(-t0, xmax - t0), no.margin = FALSE)
    ape::tiplabels(pch = 16, cex = 0.55, col = cores[tipo])
    title(sprintf("Linhagem %s  ·  %d amostras  ·  MRCA em t = %.2f",
                  c("A", "B")[tipo], ape::Ntip(tr), inicio_mrca(tr)),
          col.main = "grey30", cex.main = 1)
    if (!is.na(intro2) && intro2 > 0)
      abline(v = intro2 - t0, lty = 2, col = "grey65")
    at <- pretty(c(0, xmax))
    axis(1, at = at - t0, labels = at, col.axis = "grey45", col = "grey60", font = 2)
  }
  mtext("Tempo", side = 1, outer = TRUE, line = 1.2, col = "grey45", font = 2)
  invisible()
}

# ---- tabela-resumo ----------------------------------------------------

tabela_resumo <- function(dados) {
  res <- dados$resultados
  st  <- if (!is.null(dados$status)) dados$status else
    ifelse(vapply(res, is.list, logical(1)), "ok", "falhou")
  tempos <- if (!is.null(dados$tempos)) dados$tempos else rep(NA_real_, length(res))
  n_rep <- max(length(res), length(st))   # .rds antigos podem ter lista encurtada
  linhas <- lapply(seq_len(n_rep), function(r) {
    x <- if (r <= length(res)) res[[r]] else NULL
    base <- data.frame(Replicata = r, Status = st[r],
                       `Tempo (s)` = round(tempos[r], 2), check.names = FALSE)
    if (!is.list(x)) {
      return(data.frame(base, `t final` = NA_real_,
                        `Amostras L1` = NA_integer_, `Amostras L2` = NA_integer_,
                        `Pico L1` = NA_integer_, `Pico L2` = NA_integer_,
                        `Pontas árv. L1` = NA_integer_, `Pontas árv. L2` = NA_integer_,
                        check.names = FALSE))
    }
    z <- x[[2]]; mp <- x[[1]]
    pontas <- function(nm) if (nm %in% names(mp)) ape::Ntip(mp[[nm]]) else 0L
    data.frame(base,
               `t final`        = round(max(z$eventtimes), 2),
               `Amostras L1`    = sum(z$sampledtype == 1),
               `Amostras L2`    = sum(z$sampledtype == 2),
               `Pico L1`        = as.integer(max(z$n1)),
               `Pico L2`        = as.integer(max(z$n2)),
               `Pontas árv. L1` = pontas("tipo1"),
               `Pontas árv. L2` = pontas("tipo2"),
               check.names = FALSE)
  })
  do.call(rbind, linhas)
}

# #####################################################################
# 4. AUXILIARES DO APP
# #####################################################################

`%||%` <- function(a, b) if (is.null(a)) b else a

# Replicatas com resultado (inclui as que pararam no teto)
n_validas <- function(d) sum(vapply(d$resultados, is.list, logical(1)))

# ---- auxiliares -------------------------------------------------------

# Le "1.5" ou "1.5; 0.8" ou "1,5 0,8" -> numeric. Separador: ";" ou espaco.
# Decimal: "." ou ",".
le_vetor <- function(txt, nome) {
  if (is.null(txt)) txt <- ""
  tok <- strsplit(trimws(txt), "[;[:space:]]+")[[1]]
  tok <- sub(",$", "", tok[nzchar(tok)])
  v <- suppressWarnings(as.numeric(gsub(",", ".", tok)))
  if (anyNA(v)) stop(sprintf("Não entendi o valor de \"%s\": \"%s\"", nome, txt), call. = FALSE)
  v
}

numero <- function(x, nome, min = -Inf) {
  if (is.null(x) || length(x) != 1 || is.na(x))
    stop(sprintf("Preencha \"%s\".", nome), call. = FALSE)
  if (x < min) stop(sprintf("\"%s\" deve ser ≥ %s.", nome, min), call. = FALSE)
  x
}

# Junta e confere os parametros da interface; para com mensagem clara
monta_parametros <- function(input) {
  modelo <- input$model
  mud <- le_vetor(input$timesky, "tempos de mudança")
  if (length(mud) && (any(mud <= 0) || is.unsorted(mud, strictly = TRUE)))
    stop("Os tempos de mudança devem ser > 0 e crescentes.", call. = FALSE)
  timesky <- c(0, mud)
  k <- length(timesky)
  
  vetor_sky <- function(id, nome, min = 0, max = Inf) {
    v <- le_vetor(input[[id]], nome)
    if (!length(v)) stop(sprintf("Preencha \"%s\".", nome), call. = FALSE)
    if (length(v) == 1) v <- rep(v, k)
    if (length(v) != k)
      stop(sprintf("\"%s\" tem %d valores, mas há %d intervalo(s) de tempo (use 1 valor ou %d).",
                   nome, length(v), k, k), call. = FALSE)
    if (any(v < min | v > max))
      stop(sprintf("\"%s\" deve estar entre %s e %s.", nome, min, max), call. = FALSE)
    v
  }
  
  p <- list(
    re          = numero(input$re, "replicatas", 1),
    model       = modelo,
    seed        = input$seed,
    n           = numero(input$n, "n (amostras)", 0),
    timestop    = numero(input$timestop, "timestop", 0),
    rho         = numero(input$rho, "rho", 0),
    timeIntro2  = numero(input$timeIntro2, "introdução da Linhagem B", 0),
    N           = if (modelo == "BD") 0 else numero(input$N, "N (população)", 2),
    timesky     = timesky,
    lambdasky_1 = vetor_sky("lambda1", "λ₁"),
    lambdasky_2 = vetor_sky("lambda2", "λ₂"),
    deathsky_1  = vetor_sky("death1",  "μ₁"),
    deathsky_2  = vetor_sky("death2",  "μ₂"),
    sampprobsky = vetor_sky("samp",    "prob. de amostragem", 0, 1),
    omegasky    = if (modelo == "SIRS") vetor_sky("omega", "ω") else rep(0, k),
    aceitar_extincao = isTRUE(input$aceitar_extincao),
    max_infectados = if (modelo == "BD" && !is.na(input$max_inf) && input$max_inf > 0)
      input$max_inf else Inf,
    descartar_teto = modelo == "BD" && isTRUE(input$descartar_teto),
    limite_seg  = if (is.na(input$limite) || input$limite <= 0) Inf else input$limite
  )
  
  if (p$rho > 1) stop("rho é uma probabilidade (0 a 1).", call. = FALSE)
  if (is.finite(p$max_infectados) && p$max_infectados < 3)
    stop("O teto de infectados deve ser ≥ 3 (ou 0 para não usar teto).", call. = FALSE)
  if (p$timestop == 0 && p$n == 0)
    stop("Defina uma condição de parada: timestop > 0 ou n > 0.", call. = FALSE)
  if (p$rho > 0 && p$timestop == 0)
    stop("rho > 0 exige timestop > 0.", call. = FALSE)
  if (modelo == "SIR" && p$timestop == 0 && p$n > p$N)
    stop("No SIR, n não pode ser maior que N.", call. = FALSE)
  if (p$timestop > 0 && p$n > 0)
    showNotification("timestop > 0: a condição n é ignorada pelo simulador.",
                     type = "warning", duration = 6)
  if (p$timestop > 0 && p$timeIntro2 >= p$timestop)
    showNotification("timeIntro2 ≥ timestop: a Linhagem B não chega a entrar.",
                     type = "warning", duration = 6)
  p
}

# Fonte dos graficos (tela e PDF): Montserrat. Tenta, nesta ordem:
#   1) showtext + Google Fonts, como no seu script original (precisa de
#      internet ao abrir o app). O Shiny ja ajusta a resolucao do showtext
#      sozinho no renderPlot.
#   2) Montserrat instalada no computador: o Shiny desenha com o 'ragg',
#      que acha as fontes do sistema (funciona offline).
#   3) fonte padrao do R.
.fonte <- new.env()
fonte_graficos <- function() {
  if (!is.null(.fonte$familia)) return(.fonte$familia)
  .fonte$showtext <- FALSE
  fam <- ""
  if (requireNamespace("showtext", quietly = TRUE) &&
      requireNamespace("sysfonts", quietly = TRUE)) {
    ok <- tryCatch({
      sysfonts::font_add_google("Montserrat", "Montserrat"); TRUE
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (ok) {
      showtext::showtext_auto()
      .fonte$showtext <- TRUE
      fam <- "Montserrat"
    }
  }
  if (!nzchar(fam) && requireNamespace("ragg", quietly = TRUE) &&
      requireNamespace("systemfonts", quietly = TRUE) &&
      !isFALSE(getOption("shiny.useragg")) &&
      "Montserrat" %in% systemfonts::system_fonts()$family) {
    fam <- "Montserrat"
  }
  .fonte$familia <- fam
  fam
}

salva_pdf <- function(arquivo, largura, altura, desenha) {
  fam <- fonte_graficos()
  if (nzchar(fam) && !.fonte$showtext && capabilities("cairo")) {
    # fonte do sistema: cairo_pdf embute a Montserrat instalada
    grDevices::cairo_pdf(arquivo, width = largura, height = altura,
                         bg = "transparent", family = fam)
  } else {
    # com showtext ativo (showtext_auto), o pdf() ja sai em Montserrat
    grDevices::pdf(arquivo, width = largura, height = altura, bg = "transparent")
  }
  on.exit(grDevices::dev.off())
  desenha(fam)
}

FONTE <- fonte_graficos()   # carrega a fonte uma vez, ao abrir o app


# #####################################################################
# 5. INTERFACE
# #####################################################################

# ---- interface --------------------------------------------------------

ajuda <- function(...) tags$small(class = "text-muted d-block mb-2", ...)

barra_lateral <- sidebar(
  width = 340,
  radioButtons("model", "Modelo", c("SIR", "SIRS", "BD"),
               selected = "SIRS", inline = TRUE),
  conditionalPanel(
    "input.model == 'BD'",
    numericInput("max_inf", "Teto de infectados vivos (BD)", 1000, min = 3, step = 100),
    checkboxInput("descartar_teto", "Descartar replicatas que atingem o teto", FALSE),
    ajuda("Ao atingir o teto, a simulação para naquele instante (como um timestop ",
          "antecipado). Marcado: a replicata é descartada em vez de truncada. 0 = sem teto.")
  ),
  layout_columns(
    numericInput("re", "Replicatas", 20, min = 1, step = 1),
    numericInput("seed", "Semente", 1, step = 1)
  ),
  ajuda("Semente vazia = aleatória."),
  
  accordion(
    open = c("Taxas"),
    accordion_panel(
      "Parada",
      layout_columns(
        numericInput("timestop", "timestop", 20, min = 0),
        numericInput("n", "n (amostras)", 0, min = 0, step = 1)
      ),
      ajuda("Com timestop > 0, n é ignorado. Para parar por n, use timestop = 0."),
      numericInput("rho", "rho (amostragem final)", 0, min = 0, max = 1, step = 0.05),
      checkboxInput("aceitar_extincao", "Aceitar extinção antes de timestop", TRUE),
      ajuda("Desmarcado: a replicata que se extingue é descartada (conta como falha).")
    ),
    accordion_panel(
      "População e introdução",
      conditionalPanel("input.model != 'BD'",
                       numericInput("N", "N (tamanho da população)", 500, min = 2, step = 1)),
      numericInput("timeIntro2", "Introdução da Linhagem B (tempo)", 5, min = 0)
    ),
    accordion_panel(
      "Taxas",
      ajuda("Um valor, ou um por intervalo separados por ponto-e-vírgula (ex.: 1.5; 0.8)."),
      layout_columns(
        textInput("lambda1", "λ₁", "1.5"),
        textInput("lambda2", "λ₂", "2")
      ),
      layout_columns(
        textInput("death1", "μ₁", "1"),
        textInput("death2", "μ₂", "1")
      ),
      ajuda("μ = taxa de morte (deathsky). Cada morte vira amostra com a prob. de amostragem."),
      layout_columns(
        textInput("samp", "Prob. amostragem", "0.3"),
        conditionalPanel("input.model == 'SIRS'", textInput("omega", "ω (perda imunidade)", "0.1"))
      ),
      textInput("timesky", "Tempos de mudança (skyline)", "", placeholder = "vazio = taxas constantes")
    ),
    accordion_panel(
      "Execução",
      numericInput("limite", "Tempo máx. por replicata (s)", 60, min = 0),
      ajuda("Interrompe replicatas que explodem (comum no BD). 0 = sem limite.")
    )
  ),
  actionButton("rodar", "Rodar simulação", class = "btn-primary w-100", icon = icon("play")),
  hr(class = "my-1"),
  fileInput("arquivo", "…ou carregue um resultado salvo (.rds)", accept = ".rds",
            buttonLabel = "Abrir", placeholder = "ex.: mc_sirs_1.rds")
)

aba <- function(titulo, controles, saida) {
  nav_panel(titulo, layout_sidebar(
    sidebar = sidebar(position = "right", width = 240, open = "desktop", controles),
    saida
  ))
}

# ---- aba "Sobre o modelo" ---------------------------------------------
# Texto de introducao + duas figuras em SVG (desenhadas na propria pagina,
# na fonte do tema e nas cores das linhagens definidas em COR_TRAJ/COR_FENO).

# Figura 1: compartimentos S, I1, I2, R e as taxas entre eles
svg_compartimentos <- function(cores = COR_TRAJ, cinza = "#8A8F98") {
  caixa <- function(x, y, cor, letra, sub, w = 130, h = 64) sprintf(paste0(
    '<rect x="%.0f" y="%.0f" width="%.0f" height="%.0f" rx="10" fill="%s"/>',
    '<text x="%.0f" y="%.0f" text-anchor="middle" class="fig-letra">%s</text>',
    '<text x="%.0f" y="%.0f" text-anchor="middle" class="fig-sub">%s</text>'),
    x, y, w, h, cor, x + w / 2, y + 31, letra, x + w / 2, y + 51, sub)
  seta <- function(x1, y1, x2, y2, tracejado = FALSE) sprintf(
    '<line x1="%.0f" y1="%.0f" x2="%.0f" y2="%.0f" stroke="#6B7280" stroke-width="2.2" marker-end="url(#seta-c)"%s/>',
    x1, y1, x2, y2, if (tracejado) ' stroke-dasharray="6 5"' else "")
  rotulo <- function(x, y, txt, classe = "fig-taxa", ancora = "middle") sprintf(
    '<text x="%.0f" y="%.0f" text-anchor="%s" class="%s">%s</text>', x, y, ancora, classe, txt)
  
  paste0(
    '<svg viewBox="0 0 760 274" xmlns="http://www.w3.org/2000/svg" role="img" ',
    'aria-label="Diagrama de compartimentos S, I1, I2 e R">',
    '<defs><marker id="seta-c" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" ',
    'markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="#6B7280"/></marker></defs>',
    caixa(20, 95, cinza, "S", "suscetíveis", w = 120),
    caixa(320, 20, cores[1], "I₁", "Linhagem A"),
    caixa(320, 170, cores[2], "I₂", "Linhagem B"),
    caixa(620, 95, cinza, "R", "removidos", w = 120),
    seta(142, 112, 316, 60),  rotulo(212, 72, "λ₁"),
    seta(142, 142, 316, 194), rotulo(212, 194, "λ₂"),
    seta(452, 60, 616, 112),  rotulo(548, 74, "μ₁"),
    seta(452, 194, 616, 142), rotulo(548, 192, "μ₂"),
    '<path d="M680,161 V250 H80 V163" fill="none" stroke="#6B7280" stroke-width="2.2" ',
    'stroke-dasharray="6 5" marker-end="url(#seta-c)"/>',
    rotulo(380, 268, "ω — perda de imunidade (só SIRS)", "fig-nota"),
    '</svg>')
}

# Figura 2: exemplo de arvore filogenetica com amostragem ao longo do tempo.
# Arvore "em pe": origem em cima, presente embaixo. Cada simbolo e um tipo
# de evento do modelo; a lista ao lado (legenda_arvore) explica cada um.
COR_EVENTO <- c(origem = "#1F2937", transmissao = "#4472C4", morte = "#E0352B",
                passado = "#F28C28", presente = "#7B3294", nao_amostrado = "#6AAE5A")

# desenha o simbolo de um evento centrado em (x, y); o triangulo fica com a
# base em y (apoiado na linha do presente)
simbolo <- function(tipo, x, y, r = 6) {
  cor <- COR_EVENTO[[tipo]]
  switch(tipo,
         origem = , transmissao = , nao_amostrado =
           sprintf('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="%s"/>', x, y, r, cor),
         morte = sprintf(paste0('<path d="M%.1f %.1f l%.1f %.1f m0 %.1f l%.1f %.1f" stroke="%s" ',
                                'stroke-width="2.6" stroke-linecap="round" fill="none"/>'),
                         x - r, y - r, 2 * r, 2 * r, -2 * r, -2 * r, 2 * r, cor),
         passado = sprintf('<path d="M%.1f %.1f l%.1f %.1f l%.1f %.1f l%.1f %.1f z" fill="%s"/>',
                           x, y - 1.35 * r, 1.35 * r, 1.35 * r, -1.35 * r, 1.35 * r,
                           -1.35 * r, -1.35 * r, cor),
         presente = sprintf('<path d="M%.1f %.1f L%.1f %.1f L%.1f %.1f z" fill="%s"/>',
                            x, y - 1.9 * r, x + 1.15 * r, y, x - 1.15 * r, y, cor))
}

# icone de um simbolo, para usar na legenda e no meio do texto
icone <- function(tipo)
  HTML(sprintf('<svg class="ic" viewBox="0 0 20 20" aria-hidden="true">%s</svg>',
               simbolo(tipo, 10, if (tipo == "presente") 16 else 10)))

svg_arvore_exemplo <- function() {
  P <- 318                                         # altura da linha do presente
  linha <- function(x1, y1, x2, y2) sprintf(
    '<line x1="%.0f" y1="%.0f" x2="%.0f" y2="%.0f" stroke="#1F2937" stroke-width="2.4" stroke-linecap="square"/>',
    x1, y1, x2, y2)
  v <- function(x, y1, y2) linha(x, y1, x, y2)     # ramo (vertical)
  h <- function(y, x1, x2) linha(x1, y, x2, y)     # ligacao entre os filhos de um no
  txt <- function(x, y, s, classe = "fig-nota", ancora = "start") sprintf(
    '<text x="%.0f" y="%.0f" text-anchor="%s" class="%s">%s</text>', x, y, ancora, classe, s)
  
  paste0(
    '<svg class="arvore" viewBox="0 0 488 348" xmlns="http://www.w3.org/2000/svg" role="img" ',
    'aria-label="Exemplo de árvore filogenética com amostragem ao longo do tempo">',
    '<defs><marker id="seta-t" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" ',
    'markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="#6B7280"/></marker></defs>',
    # seta do tempo (esquerda) e linha do presente
    '<line x1="20" y1="40" x2="20" y2="300" stroke="#6B7280" stroke-width="2" marker-end="url(#seta-t)"/>',
    '<text transform="translate(13 170) rotate(-90)" text-anchor="middle" class="fig-nota">tempo</text>',
    '<line x1="40" y1="318" x2="390" y2="318" stroke="#6B7280" stroke-width="1.6" stroke-dasharray="5 4"/>',
    # ramos
    v(205, 34, 78), h(78, 110, 300),                        # origem -> raiz (MRCA)
    v(110, 78, 232), h(232, 75, 145), v(75, 232, P), v(145, 232, P),
    v(300, 78, 150), h(150, 250, 355),
    v(355, 150, 258),                                        # termina amostrado no passado
    v(250, 150, 188), h(188, 215, 290),
    v(215, 188, 252),                                        # termina sem amostra
    v(290, 188, 282), h(282, 268, 312), v(268, 282, P), v(312, 282, P),
    # eventos (desenhados por cima dos ramos)
    simbolo("origem", 205, 34),
    simbolo("transmissao", 205, 78),  simbolo("transmissao", 110, 232),
    simbolo("transmissao", 300, 150), simbolo("transmissao", 250, 188),
    simbolo("transmissao", 290, 282),
    simbolo("morte", 215, 252),
    simbolo("passado", 355, 258),
    simbolo("nao_amostrado", 75, P),
    simbolo("presente", 145, P), simbolo("presente", 268, P), simbolo("presente", 312, P),
    # rotulos
    txt(217, 34, "Origem", "fig-rot"), txt(217, 50, "t = 0"),
    txt(196, 70, "MRCA (raiz)", ancora = "end"),
    txt(398, 322, "Presente", "fig-rot"), txt(398, 338, "t = timestop"),
    '</svg>')
}

legenda_arvore <- function() {
  item <- function(tipo, titulo, texto)
    tags$li(icone(tipo), tags$div(tags$b(HTML(titulo)), tags$br(), HTML(texto)))
  tags$ul(class = "simbolos",
          item("origem", "Origem",
               "O primeiro infectado da linhagem, em t = 0 (na Linhagem B, em t = <code>timeIntro2</code>)."),
          item("transmissao", "Nascimento (λ) → transmissão",
               "Um infectado transmite o patógeno para outra pessoa e o ramo se divide em dois."),
          item("morte", "Morte (μ) → fim do período infeccioso",
               "O infectado deixa de transmitir (recupera, é isolado ou morre) <b>sem</b> ser amostrado,
         com taxa μ · (1 − δ). O ramo termina e não aparece nos dados."),
          item("passado", "Amostragem no passado (ψ)",
               "O infectado é removido e o patógeno é sequenciado antes do presente, com taxa
         ψ = μ · δ. Vira uma ponta da árvore."),
          item("presente", "Amostragem no presente (ρ)",
               "No fim da simulação (t = <code>timestop</code>), cada infectado ainda vivo é
         sequenciado com probabilidade ρ."),
          item("nao_amostrado", "Não amostrado no presente (1 − ρ)",
               "Infectado ainda vivo em <code>timestop</code> que não foi sequenciado. Também não
         aparece nos dados.")
  )
}

# Uma linha da tabela de parametros
linha_par <- function(nome, codigo, descricao)
  tags$tr(tags$td(HTML(nome)), tags$td(tags$code(codigo)), tags$td(HTML(descricao)))
grupo_par <- function(titulo)
  tags$tr(class = "grupo", tags$td(colspan = 3, titulo))

aba_sobre <- nav_panel(
  "Sobre o modelo",
  div(class = "sobre",
      
      tags$h3("Sobre o modelo"),
      p(HTML("Este app simula uma epidemia em que <b>duas linhagens</b> (por exemplo, duas variantes
      de um vírus) circulam ao mesmo tempo. Cada simulação produz a <b>trajetória</b> do número de
      infectados ao longo do tempo e a
      <b>árvore filogenética</b> das amostras de cada linhagem.")),
      
      p(HTML("O modelo é um <b>processo de nascimento e morte multitipo</b> (<i>multi-type
      birth-death</i>, MTBD). Cada infectado pode “nascer” (transmitir, gerando um novo infectado)
      ou “morrer” (deixar de ser infeccioso). O <b>tipo</b> é a linhagem do infectado (A ou B):
      quem ele infecta herda o mesmo tipo, e cada tipo tem suas próprias taxas de nascimento (λ) e
      de morte (μ).")),
      
      # ---------------------------------------------------------------
      tags$h4("1. O modelo epidemiológico"),
      p(HTML("Os eventos acontecem um de cada
      vez, em instantes sorteados ao acaso (simulação estocástica em tempo contínuo).
      Há três tipos de evento:")),
      tags$ul(
        tags$li(HTML("<b>Transmissão.</b> Um infectado da linhagem <i>i</i> infecta um suscetível,
        com taxa <b>λ<sub><i>i</i></sub></b> por infectado (no Birth-Death puro, com taxa λ<sub><i>i</i></sub>
        constante). Na árvore, isso vira um nó: a linhagem se divide em duas.")),
        tags$li(HTML("<b>Remoção.</b> O infectado deixa de transmitir (recupera, é isolado ou morre),
        com taxa <b>μ<sub><i>i</i></sub></b>. Em cada remoção, com probabilidade <b>δ</b>, o
        patógeno dessa pessoa é amostrado (sequenciado) e vira uma ponta da árvore. Assim, a taxa
        de amostragem é ψ = μ · δ e a de remoção sem amostra é μ · (1 − δ). Ser amostrado também
        remove o infectado.")),
        tags$li(HTML("<b>Perda de imunidade</b> (só no SIRS). Cada recuperado volta a ser suscetível
        com taxa <b>ω</b>."))
      ),
      div(class = "figura", HTML(svg_compartimentos()),
          p(class = "legenda-fig", "Fluxo entre os compartimentos no SIR e no SIRS. Nesses modelos, a
          taxa efetiva de transmissão é λ · S/N, que cai à medida que os suscetíveis acabam. No BD
          não existem S nem R: cada infectado transmite com taxa λ constante.")),
      
      p(HTML("<b>Os três modelos</b>")),
      tags$ul(
        tags$li(HTML("<b>SIR:</b> população fixa de N pessoas. Cada transmissão gasta um suscetível;
        quando eles ficam escassos, a transmissão cai e a epidemia termina. Os recuperados ficam
        imunes para sempre.")),
        tags$li(HTML("<b>SIRS:</b> como o SIR, mas a imunidade dura em média 1/ω. Os recuperados
        voltam a ser suscetíveis, o que permite novas ondas e circulação endêmica.")),
        tags$li(HTML("<b>BD (nascimento e morte):</b> sem população finita, cada infectado transmite
        com taxa λ constante. Se λ &gt; μ, o número de infectados cresce exponencialmente sem
        limite; por isso existe o teto de infectados."))
      ),
      p(HTML("<b>As duas linhagens.</b> A Linhagem A começa com um infectado em t = 0; a Linhagem B
      é introduzida com um infectado em t = <code>timeIntro2</code>. Cada linhagem tem seus próprios
      λ e μ. No SIR e no SIRS elas competem pelos mesmos
      suscetíveis; no BD crescem de forma independente. Se a Linhagem A se extinguir antes de
      <code>timeIntro2</code>, a simulação termina e a Linhagem B não chega a entrar.")),
      
      # ---------------------------------------------------------------
      tags$h4("2. O que é uma árvore filogenética"),
      p(HTML("Uma árvore filogenética mostra a história de transmissão entre os patógenos: quem
      descende de quem e há quanto tempo eles compartilham um ancestral. No exemplo abaixo, o tempo
      corre de cima (origem) para baixo (presente), e cada símbolo marca um tipo de evento do
      modelo.")),
      div(class = "figura figura-arvore",
          div(class = "arvore-flex", HTML(svg_arvore_exemplo()), legenda_arvore()),
          p(class = "legenda-fig", "Exemplo de árvore filogenética com amostragem ao longo do tempo.")),
      p(HTML("<b>A árvore que vemos nos dados.</b> Só conhecemos os patógenos que foram sequenciados
      ("), icone("passado"), " e ", icone("presente"), HTML("). Por isso, a árvore usada na análise
      (e desenhada na aba <b>Árvores</b>) é a <b>árvore reconstruída</b>: os ramos que terminam em "),
        icone("morte"), " ou ", icone("nao_amostrado"), HTML(" são podados, e um nó que fica com um só
      descendente amostrado deixa de ser nó. No exemplo, a árvore reconstruída tem 4 pontas e 3 nós.")),
      tags$ul(
        tags$li(HTML("<b>Pontas (folhas):</b> Cada ponta é um infectado amostrado, no instante da
        amostragem (no passado ou no presente).")),
        tags$li(HTML("<b>Nós internos:</b> Cada nó é um evento de transmissão que deu origem a duas
        cadeias com descendentes amostrados.")),
        tags$li(HTML("<b>Ramos:</b> O comprimento de um ramo é o tempo entre dois eventos.")),
        tags$li(HTML("<b>MRCA:</b> o MRCA é o ancestral comum mais recente de todas as
        amostras de uma linhagem, ou seja, a raiz da árvore reconstruída."))
      ),
      p(HTML("O app monta uma árvore separada para cada linhagem, desde que ela tenha pelo menos 2
      amostras. Na aba <b>Árvores</b>, elas aparecem deitadas, com o tempo no eixo horizontal. A forma
      da árvore guarda a marca da dinâmica da epidemia: quando a transmissão é rápida, muitos nós se
      concentram em pouco tempo; quando ela diminui, os nós ficam mais espaçados.")),
      
      # ---------------------------------------------------------------
      tags$h4("3. Parâmetros"),
      p(HTML("Os parâmetros do modelo são:")),
      tags$table(class = "table table-sm align-middle",
                 tags$thead(tags$tr(tags$th("Parâmetro"), tags$th("No código"), tags$th("O que faz"))),
                 tags$tbody(
                   grupo_par("Modelo e replicatas"),
                   linha_par("Modelo", "model", "SIR, SIRS ou BD (veja acima)."),
                   linha_par("Replicatas", "re", "Quantas epidemias simular com os mesmos parâmetros."),
                   linha_par("Semente", "seed", "Fixa o sorteio aleatório: a mesma semente repete os mesmos
          resultados. Vazio = aleatória."),
                   linha_par("Teto de infectados vivos <small>(BD)</small>", "max_infectados", "Para a
          simulação quando o total de infectados vivos (duas linhagens) atinge esse valor.
          0 = sem teto."),
                   linha_par("Descartar no teto <small>(BD)</small>", "descartar_teto", "Marcado: a replicata
          que atinge o teto é descartada. Desmarcado: ela é mantida até o instante em que parou."),
                   
                   
                   grupo_par("População e introdução"),
                   linha_par("N", "N", "Tamanho da população (SIR e SIRS). No início, todos são suscetíveis,
          menos os infectados iniciais."),
                   linha_par("Introdução da Linhagem B", "timeIntro2", "Instante em que entra o primeiro
          infectado da Linhagem B. timeIntro2 = 0: As duas começam juntas."),
                   
                   grupo_par("Taxas"),
                   linha_par("λ<sub>1</sub>, λ<sub>2</sub>", "lambdasky_1, lambdasky_2", "Taxa de transmissão
          de cada linhagem, por infectado e por unidade de tempo."),
                   linha_par("μ<sub>1</sub>, μ<sub>2</sub>", "deathsky_1, deathsky_2", "Taxa de remoção de cada linhagem."),
                   linha_par("Prob. amostragem (δ)", "sampprobsky", "Probabilidade de cada remoção virar
          amostra. É a mesma para as duas linhagens."),
                   linha_par("ω", "omegasky", "Taxa de perda de imunidade (só SIRS)."),
                   linha_par("Tempos de mudança", "timesky", "Instantes em que as taxas mudam (skyline).")
                 )
      ),
      
      # ---------------------------------------------------------------
      tags$h4("4. O que cada aba mostra"),
      tags$ul(
        tags$li(HTML("<b>Trajetórias:</b> Número de infectados vivos de cada linhagem ao longo do
        tempo; o tracejado marca a entrada da Linhagem B. Dá para trocar
        por suscetíveis, recuperados ou amostras acumuladas e mostrar a mediana.")),
        tags$li(HTML("<b>Árvores:</b> As árvores reconstruídas de uma replicata, no eixo de tempo
        absoluto.")),
        tags$li(HTML("<b>Resumo:</b> Uma linha por replicata, com status, número de amostras, pico de
        infectados e tamanho das árvores.")),
        tags$li(HTML("<b>Log:</b> as mensagens do simulador em cada replicata."))
      )
  )
)

# Cabecalhos que abrem/fecham (Parada, Populacao, Taxas, Execucao) na cor
# primaria do tema (#073A70), abertos ou fechados, com texto e seta brancos.
# O conteudo da secao aberta ganha um fundo azul bem claro.
# var(--bs-primary) segue o 'primary' do bs_theme: mudou la, muda aqui.
css_secoes <- "
.sidebar .accordion .accordion-button,
.sidebar .accordion .accordion-button:not(.collapsed) {
  background-color: var(--bs-primary);
  color: #fff;
  font-weight: 600;
  box-shadow: none;
}
.sidebar .accordion .accordion-button:hover {
  background-color: rgba(var(--bs-primary-rgb), .85);
}
.sidebar .accordion .accordion-button:focus {
  box-shadow: 0 0 0 .2rem rgba(var(--bs-primary-rgb), .35);
}
.sidebar .accordion .accordion-button::after {
  filter: brightness(0) invert(1);              /* seta branca */
}
.sidebar .accordion .accordion-item {
  border-color: #fff;                           /* separa as faixas azuis */
}
.sidebar .accordion .accordion-body {
  background-color: rgba(var(--bs-primary-rgb), .05);
}
/* aba Sobre o modelo */
.sobre { max-width: 920px; margin: 0 auto; padding: .25rem 1rem 2rem; line-height: 1.65; }
.sobre h3 { color: var(--bs-primary); font-weight: 700; margin-bottom: .75rem; }
.sobre h4 { color: var(--bs-primary); font-weight: 700; margin-top: 2.2rem;
            padding-bottom: .35rem; border-bottom: 2px solid rgba(var(--bs-primary-rgb), .15); }
.sobre li { margin-bottom: .35rem; }
.sobre code { color: var(--bs-primary); }
.sobre .figura { margin: 1.25rem 0 1.5rem; text-align: center; }
.sobre .figura svg { width: 100%; max-width: 760px; height: auto; }
.sobre .legenda-fig { font-size: .85rem; color: var(--bs-secondary-color); margin-top: .4rem; }
.sobre table td:first-child { white-space: nowrap; }
.sobre tr.grupo td { background: rgba(var(--bs-primary-rgb), .07); color: var(--bs-primary);
                     font-weight: 700; }
.fig-letra { font-size: 24px; font-weight: 700; fill: #fff; }
.fig-sub   { font-size: 12px; fill: #fff; }
.fig-taxa  { font-size: 15px; font-weight: 600; fill: #374151; }
.fig-nota  { font-size: 12.5px; fill: #6B7280; }
.fig-tit   { font-size: 15px; font-weight: 700; fill: var(--bs-primary); }
.fig-rot   { font-size: 13.5px; font-weight: 700; fill: #1F2937; }
/* figura da arvore: arvore e legenda lado a lado; em tela estreita a legenda
   desce para baixo da arvore (nao precisa arrastar para o lado) */
.sobre .figura-arvore .arvore-flex { display: flex; flex-wrap: wrap; align-items: center;
  justify-content: center; gap: 1rem 2rem; }
.sobre .figura.figura-arvore svg.arvore { flex: 1 1 340px; width: 100%; max-width: 440px;
  min-width: 0; }
.sobre .simbolos { flex: 1 1 260px; max-width: 400px; list-style: none; padding: 0; margin: 0;
  text-align: left; font-size: .88rem; line-height: 1.4; }
.sobre .simbolos li { display: flex; gap: .65rem; align-items: flex-start; margin-bottom: .75rem; }
.sobre svg.ic { width: 1em !important; height: 1em !important; min-width: 0 !important;
  max-width: none !important; vertical-align: -.15em; }
.sobre .simbolos svg.ic { flex: 0 0 auto; width: 18px !important; height: 18px !important;
  margin-top: .1rem; }
.sobre .figura.figura-arvore .legenda-fig::before { content: none; }
.irs, .selectize-input, .selectize-dropdown, .shiny-notification,
.form-control, .btn, table {
  font-family: inherit;                         /* tudo na fonte do tema */
}
/* barra do titulo: azul com letras brancas */
.bslib-page-title,
.bslib-page-sidebar > .navbar {
  background-color: var(--bs-primary) !important;
  color: #fff !important;
  border-bottom: none !important;
  font-size: 1.2rem !important;                 /* titulo longo: cabe numa linha */
}
.bslib-page-sidebar > .navbar .navbar-brand { color: #fff !important; font-size: 1.2rem !important; }

/* titulo longo quebra linha em qualquer versao do bslib. Nas versoes novas
   o titulo e um 'navbar-brand' do Bootstrap, que por padrao NUNCA quebra
   linha: o titulo ficava mais largo que a janela e empurrava a pagina toda
   para o lado (barra lateral cortada, faixa azul terminando antes da borda) */
.bslib-page-title,
.bslib-page-sidebar > .navbar .navbar-brand {
  white-space: normal !important;
  min-width: 0;
  max-width: 100%;
}
.bslib-page-sidebar > .navbar .container-fluid { min-width: 0; }

/* janelas estreitas (computador com janela pequena, tablet): a tabela de
   parametros da aba Sobre quebra linha em vez de sair do cartao */
@media (max-width: 991.98px) {
  .sobre table td:first-child { white-space: normal; }
  .sobre table code { word-break: break-word; }
}

/* aviso que aparece so quando o app esta encolhido num celular */
#aviso-celular { display: none; position: fixed; top: 0; left: 0; right: 0; z-index: 3000;
  background: #FFF3CD; color: #5C4400; padding: .7em 2.6em .7em 1em; line-height: 1.35;
  box-shadow: 0 2px 6px rgba(0,0,0,.2); }
#aviso-celular a { color: var(--bs-primary); font-weight: 700; }
#aviso-celular .fechar { position: absolute; right: .6em; top: .35em; font-size: 1.4em;
  cursor: pointer; background: none; border: 0; color: #5C4400; }

/* celular: a aba Sobre ocupa so a largura da tela */
.sobre { width: 100%; min-width: 0; }
#tabela, #log { width: 100%; min-width: 0; overflow-x: auto; }   /* Resumo e Log rolam para o lado */
@media (max-width: 575.98px) {
  .bslib-page-title,
  .bslib-page-sidebar > .navbar .navbar-brand { font-size: 1rem !important; line-height: 1.3; }
  .sobre { padding: 0 .25rem 1.5rem; }
  .sobre .figura { overflow-x: auto; }                 /* figura: arrasta para o lado */
  .sobre .figura svg { min-width: 620px; }
  .sobre .legenda-fig::before { content: '↔ arraste a figura para o lado. '; font-weight: 600; }
  .sobre table td:first-child { white-space: normal; }
  .sobre table th:nth-child(2),
  .sobre table td:nth-child(2) { display: none; }      /* esconde a coluna 'No codigo' */
}
"

# Montserrat em toda a interface: o navegador baixa do Google Fonts; sem
# internet, usa a Montserrat instalada no computador (se houver) ou a Segoe UI.
tema <- bs_theme(version = 5, primary = "#073A70",
                 base_font = font_collection(
                   font_google("Montserrat", local = FALSE, wght = c(400, 500, 600, 700)),
                   "Segoe UI", "Roboto", "Helvetica Neue", "Arial", "sans-serif")) |>
  bs_add_rules(css_secoes)

ui <- page_sidebar(
  title = "Modelos epidemiológicos de Nascimento e Morte dependentes da densidade para árvores filogenéticas · SIR / SIRS / BD",
  theme = tema,
  fillable = FALSE,   # a pagina rola em vez de espremer os graficos (celular deitado)
  sidebar = barra_lateral,
  navset_card_tab(
    id = "abas",
    aba_sobre,
    aba("Trajetórias",
        tagList(
          selectInput("variavel", "Variável", VARIAVEIS_TRAJ),
          sliderInput("alpha_traj", "Opacidade", 0.05, 1, 0.4, step = 0.05),
          checkboxInput("mediana", "Mostrar mediana", FALSE),
          downloadButton("pdf_traj", "PDF", class = "btn-sm w-100")
        ),
        plotOutput("g_traj", height = "520px")),
    aba("Árvores",
        tagList(
          selectInput("rep_arvore", "Replicata", choices = NULL),
          ajuda("Eixo em tempo absoluto; o tracejado marca a entrada da Linhagem B."),
          downloadButton("nexus", "Todas as árvores (.nex)", class = "btn-sm w-100")
        ),
        plotOutput("g_arvores", height = "680px")),
    nav_panel("Resumo",
              uiOutput("caixas"),
              div(class = "d-flex justify-content-end mb-2",
                  downloadButton("rds", "Baixar resultados (.rds)", class = "btn-sm")),
              tableOutput("tabela")),
    nav_panel("Log", verbatimTextOutput("log"))
  ),
  # telas pequenas (celular, tablet): as barras laterais se recolhem para
  # sobrar espaco para os graficos. Os controles de cada aba comecam
  # fechados, e a barra de parametros fecha ao clicar em "Rodar simulacao".
  # No computador (tela >= 992 px) nada muda. As setas reabrem as barras.
  div(id = "aviso-celular"),
  tags$script(HTML("
    // Celular mostrando o app 'encolhido': dentro de outra pagina (ex.: painel
    // do Connect Cloud) ou com o navegador em 'Site para computador'.
    // Nesses casos aparece um aviso no topo explicando o que fazer.
    $(function() {
      var menorLado = Math.min(screen.width, screen.height);
      if (menorLado >= 600) return;                          // nao e celular
      var paisagem = (typeof window.orientation === 'number' && Math.abs(window.orientation) === 90) ||
                     screen.width > screen.height;
      var larguraTela = paisagem ? Math.max(screen.width, screen.height) : menorLado;
      var fator = window.innerWidth / larguraTela;           // quanto a pagina foi encolhida
      var emOutraPagina = window.self !== window.top;
      if (!emOutraPagina && fator < 1.5) return;             // tudo certo
      var msg = emOutraPagina
        ? 'O app está aberto dentro de outra página e fica pequeno no celular. ' +
          '<a href=\"' + location.href + '\" target=\"_blank\" rel=\"noopener\">Toque aqui para abrir em tela cheia</a>.'
        : 'Seu navegador está mostrando a versão para computador, por isso o app fica pequeno. ' +
          'No Chrome, toque em ⋮ e desmarque <b>Site para computador</b>.';
      var aviso = document.getElementById('aviso-celular');
      aviso.innerHTML = msg + '<button class=\"fechar\" aria-label=\"Fechar\">×</button>';
      aviso.style.fontSize = Math.round(15 * Math.max(fator, 1)) + 'px';
      aviso.style.display = 'block';
      aviso.querySelector('.fechar').onclick = function() { aviso.style.display = 'none'; };
    });

    function fecharPaineis(incluirPrincipal) {
      if (window.innerWidth >= 992) return;
      document.querySelectorAll('.bslib-sidebar-layout > .collapse-toggle').forEach(function(t) {
        var layout = t.parentElement;
        var principal = !layout.classList.contains('sidebar-right');
        var aberta = !layout.classList.contains('sidebar-collapsed');
        if ((incluirPrincipal || !principal) && aberta) t.click();
      });
    }
    $(document).on('shiny:connected', function() { setTimeout(function() { fecharPaineis(false); }, 300); });
    $(document).on('click', '#rodar', function() { fecharPaineis(true); });
  "))
)


# #####################################################################
# 6. SERVIDOR
# #####################################################################

# ---- servidor ---------------------------------------------------------

server <- function(input, output, session) {
  
  dados <- reactiveVal(NULL)
  
  # ---- rodar ----
  observeEvent(input$rodar, {
    p <- tryCatch(monta_parametros(input), error = function(e) {
      showNotification(conditionMessage(e), type = "error", duration = 8); NULL
    })
    if (is.null(p)) return()
    
    res <- withProgress(message = "Simulando", value = 0, {
      do.call(mc_MTBD_SIR_SIRS, c(p, list(
        verbose = FALSE,
        ao_progredir = function(r, re)
          incProgress(1 / re, detail = sprintf("replicata %d de %d", r, re))
      )))
    })
    dados(res)
    if (identical(input$abas, "Sobre o modelo")) nav_select("abas", "Trajetórias")
    
    n_ok   <- n_validas(res)
    n_teto <- sum(grepl("teto", res$status))
    showNotification(
      paste0(sprintf("%d de %d replicatas válidas", n_ok, res$re),
             if (n_teto) sprintf(" · %d atingiram o teto", n_teto) else "",
             sprintf(" · %.1f s", res$elapsed_total)),
      type = if (n_ok == res$re) "message" else "warning", duration = 6)
  })
  
  # ---- carregar .rds ----
  observeEvent(input$arquivo, {
    obj <- tryCatch(readRDS(input$arquivo$datapath), error = function(e) NULL)
    if (is.list(obj) && is.null(obj$resultados) && length(obj) && is.list(obj[[1]]))
      obj <- list(resultados = obj)        # aceita tambem a lista crua de replicatas
    if (!is.list(obj) || is.null(obj$resultados)) {
      showNotification("Arquivo não reconhecido: esperava a saída de mc_MTBD_SIR_SIRS.",
                       type = "error", duration = 8)
      return()
    }
    if (is.null(obj$re)) obj$re <- length(obj$resultados)
    if (is.null(obj$status))
      obj$status <- ifelse(vapply(obj$resultados, is.list, logical(1)), "ok", "falhou")
    dados(obj)
    if (identical(input$abas, "Sobre o modelo")) nav_select("abas", "Trajetórias")
    showNotification(sprintf("Carregado: %d replicatas (%d ok).",
                             length(obj$resultados), n_validas(obj)),
                     type = "message")
  })
  
  # ---- derivados ----
  pronto <- reactive({
    validate(need(dados(), "Clique em “Rodar simulação” ou carregue um .rds salvo."))
    dados()
  })
  trajs   <- reactive(trajetorias_de(pronto()$resultados))
  info    <- reactive(extrai_arvores(pronto()$resultados, min_tips = 3))
  intro2  <- reactive(tempo_intro2(trajs()))
  xmax    <- reactive(tempo_maximo(trajs(), info()))
  
  observe({
    d <- dados(); req(d)
    ok <- which(vapply(d$resultados, is.list, logical(1)))
    updateSelectInput(session, "rep_arvore",
                      choices = setNames(ok, sprintf("Replicata %d", ok)))
  })
  
  # ---- graficos ----
  desenha_traj <- function(fam = "") plot_trajetorias(
    trajs(), input$variavel, COR_TRAJ, input$alpha_traj, intro2(), xmax(),
    input$mediana, fam)
  
  output$g_traj <- renderPlot(desenha_traj(FONTE), res = 96)
  output$g_arvores <- renderPlot({
    req(input$rep_arvore)
    r <- as.integer(input$rep_arvore)
    x <- pronto()$resultados[[r]]
    validate(need(is.list(x), "Replicata sem resultado."))
    plot_arvores_rep(x[[1]], COR_FENO, xmax(), intro2(), FONTE)
  }, res = 96)
  
  # ---- resumo ----
  output$caixas <- renderUI({
    d <- pronto()
    i <- info()
    tempo <- if (!is.null(d$elapsed_total)) sprintf("%.1f s", d$elapsed_total) else "—"
    layout_columns(
      fill = FALSE,
      value_box("Replicatas válidas",
                sprintf("%d / %d", n_validas(d), d$re %||% length(d$resultados)),
                if (any(grepl("teto", d$status)))
                  tags$small(sprintf("%d truncadas no teto · %d descartadas",
                                     sum(d$status == "ok · parou no teto"),
                                     sum(d$status == "descartada (teto)")))),
      value_box("Tempo total", tempo)
    )
  })
  output$tabela <- renderTable(tabela_resumo(pronto()), na = "—",
                               striped = TRUE, hover = TRUE, spacing = "s")
  
  output$log <- renderText({
    d <- pronto()
    if (is.null(d$logs)) return("Sem log (resultado carregado de arquivo).")
    linhas <- vapply(seq_along(d$logs), function(r) {
      l <- sub('^\\[1\\] "?(.*?)"?$', "\\1", d$logs[[r]], perl = TRUE)
      cab <- sprintf("rep %02d  %-28s %6.2f s", r, d$status[r], d$tempos[r])
      if (!length(l)) return(cab)
      extra <- if (length(l) > 4) sprintf("\n    … (+%d linhas)", length(l) - 4) else ""
      paste0(cab, "\n", paste0("    ", head(l, 4), collapse = "\n"), extra)
    }, character(1))
    paste(linhas, collapse = "\n")
  })
  
  # ---- downloads ----
  output$pdf_traj <- downloadHandler(
    filename = function() sprintf("trajetorias_%s.pdf", pronto()$model %||% "MC"),
    content  = function(file) salva_pdf(file, 8, 5, desenha_traj)
  )
  output$rds <- downloadHandler(
    filename = function() sprintf("mc_%s.rds", tolower(pronto()$model %||% "resultado")),
    content  = function(file) saveRDS(pronto(), file)
  )
  output$nexus <- downloadHandler(
    filename = function() sprintf("arvores_%s.nex", pronto()$model %||% "MC"),
    content  = function(file) {
      res <- pronto()$resultados
      todas <- list()
      for (r in seq_along(res)) {
        if (!is.list(res[[r]])) next
        mp <- res[[r]][[1]]
        for (i in seq_along(mp)) todas[[sprintf("rep%02d_%s", r, names(mp)[i])]] <- mp[[i]]
      }
      class(todas) <- "multiPhylo"
      ape::write.nexus(todas, file = file, translate = FALSE)
    }
  )
}

shinyApp(ui, server)