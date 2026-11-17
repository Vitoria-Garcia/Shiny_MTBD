# =====================================================================
# Simulador MTBD (duas linhagens) para SIR / SIRS / BD
# Funcao original do usuario, SEM alteracoes. Se voce atualizar a funcao,
# basta substituir este arquivo: o app le tudo que estiver na pasta R/.
# =====================================================================

MTBD.SIR.SIRS <- function(n, lambdasky_1, lambdasky_2, deathsky_1, deathsky_2,
                          timesky, sampprobsky, rho = 0, timestop = 0,
                          model = "SIR", N = 0, trackinfecteds = FALSE,
                          omegasky = rep(0, length(timesky)),
                          timeIntro2 = 0,
                          aceitar_extincao = FALSE) {   # [ALTERADO] novo argumento
  # [ALTERADO] aceitar_extincao = TRUE: se os infectados zerarem antes de
  #   timestop (no SIR isso e o FIM NATURAL da epidemia, nao so extincao
  #   precoce), a simulacao termina e as arvores sao devolvidas. Quem chama
  #   decide o que descartar (ex.: < k amostras). Com FALSE, comportamento
  #   original: a replicata inteira e jogada fora (devolve 0).
  
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
  
  # Linhagem 1
  musky_1  <- deathsky_1 * (1 - psampsky)
  psisky_1 <- deathsky_1 * psampsky
  
  # Linhagem 2
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
    print("timeIntro2 >= timestop: a linhagem 2 nao chega a ser introduzida")
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
    
    # Quando a linhagem 2 entra
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
          # Introdução da linhagem 2
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
    time <- timestop
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
      prop_log = prop_log ))}
  else {
    arvores
  }
}