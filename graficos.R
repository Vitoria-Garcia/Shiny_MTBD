# =====================================================================
# Extracao dos resultados e graficos (base R), usados pelo app e pelos PDFs
# =====================================================================

COR_TRAJ <- c("#935999", "#073A70")   # trajetorias (cores do seu PDF)
COR_FENO <- c("#935999", "#073A70")   # fenogramas e arvores

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

# Caráter BM em cada arvore. Devolve, para cada arvore, um vetor indexado
# pelo numero do no (1..Ntip = pontas, Ntip+1.. = nos internos).
#   internos = "fastAnc"   -> pontas simuladas + nos reconstruidos (seu script)
#   internos = "simulados" -> valores verdadeiros simulados em todos os nos
simula_caracteres <- function(arvores, semente = 1, sig2 = 1,
                              internos = c("fastAnc", "simulados")) {
  internos <- match.arg(internos)
  if (!is.null(semente) && !is.na(semente)) set.seed(semente)
  lapply(arvores, function(tr) {
    tryCatch({
      nt <- ape::Ntip(tr)
      vv <- numeric(nt + tr$Nnode)
      if (internos == "simulados") {
        v <- phytools::fastBM(tr, sig2 = sig2, internal = TRUE)
        vv[match(names(v)[seq_len(nt)], tr$tip.label)] <- v[seq_len(nt)]
        nos <- v[-seq_len(nt)]
        vv[as.integer(names(nos))] <- nos
      } else {
        x   <- phytools::fastBM(tr, sig2 = sig2)
        anc <- phytools::fastAnc(tr, x)
        vv[seq_len(nt)] <- x[tr$tip.label]
        vv[as.integer(names(anc))] <- anc
      }
      vv
    }, error = function(e) NULL)
  })
}

# Introducao da linhagem 2 a partir dos dados (como no seu script)
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

grafico_vazio <- function(msg) {
  op <- par(mar = c(0, 0, 0, 0)); on.exit(par(op))
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
  if (!length(trajs)) return(grafico_vazio("Nenhuma trajetória disponível"))
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
        lines(grade, m, type = "s", lwd = 2.5, col = if (k == 1) "grey30" else cores[2])
      }
    }
    legend("topright", c("Linhagem 1", "Linhagem 2"), col = cores,
           text.col = c("grey45", cores[2]), lwd = 2, bty = "n", text.font = 2)
  } else {
    vals <- unlist(lapply(trajs, `[[`, variavel))
    if (!length(vals) || all(is.na(vals)))
      return(grafico_vazio("Esta variável não existe neste modelo (ex.: S no BD, R no SIR)"))
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

desenha_fenograma <- function(tr, vv, inicio = 0, col = "#0000ff20", lwd = 1) {
  H <- phytools::nodeHeights(tr)
  segments(H[, 1] + inicio, vv[tr$edge[, 1]],
           H[, 2] + inicio, vv[tr$edge[, 2]], col = col, lwd = lwd)
}

plot_fenogramas <- function(info, valores, cores = COR_FENO, alpha = 0.4,
                            intro2 = NA, xmax = NULL, familia = "") {
  ok <- !vapply(valores, is.null, logical(1))
  if (!any(ok)) return(grafico_vazio("Nenhuma árvore com >= 3 pontas"))
  if (is.null(xmax)) xmax <- tempo_maximo(list(), info)
  ylim <- range(unlist(valores[ok]))               # mesmo eixo Y nos 3 paineis

  op <- par(mfrow = c(3, 1), mar = c(0.5, 4.5, 2.0, 1.5), oma = c(4, 0, 0, 0),
            family = familia, font.axis = 2, font.lab = 2, font.main = 2)
  on.exit(par(op))

  painel <- function(idx, titulo, eixo_x) {
    do.call(plot, c(list(NA, xlim = c(0, xmax), ylim = ylim, bty = "n",
                         xaxt = if (eixo_x) "s" else "n", xlab = "",
                         ylab = "caráter (BM)", main = titulo, col.main = "grey30"),
                    estilo_eixos))
    linha_intro2(intro2)
    if (!length(idx))
      text(xmax / 2, mean(ylim), "sem árvores com >= 3 pontas", col = "grey60")
    for (i in idx)
      desenha_fenograma(info$arvores[[i]], valores[[i]], info$inicio[i],
                        col = adjustcolor(cores[info$tipos[i]], alpha.f = alpha))
  }
  painel(which(ok & info$tipos == 1L), "Linhagem 1", FALSE)
  painel(which(ok & info$tipos == 2L), "Linhagem 2", FALSE)
  painel(which(ok), "Sobrepostas", TRUE)
  legend("topleft", c("Linhagem 1", "Linhagem 2"), col = cores, lwd = 3,
         bty = "n", text.col = "grey45", text.font = 2)
  mtext("Tempo", side = 1, outer = TRUE, line = 2.5, col = "grey45", font = 2)
  invisible()
}

# Arvores de uma replicata, no eixo de tempo absoluto
plot_arvores_rep <- function(mp, cores = COR_FENO, xmax, intro2 = NA,
                             familia = "") {
  if (!inherits(mp, "multiPhylo") || !length(mp))
    return(grafico_vazio("Replicata sem árvores"))
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
    title(sprintf("Linhagem %d  ·  %d amostras  ·  MRCA em t = %.2f",
                  tipo, ape::Ntip(tr), inicio_mrca(tr)),
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
  linhas <- lapply(seq_along(res), function(r) {
    x <- res[[r]]
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
