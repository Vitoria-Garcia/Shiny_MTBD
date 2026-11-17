# =====================================================================
# App Shiny: simulacao Monte Carlo MTBD (duas linhagens) SIR / SIRS / BD
#
# Estrutura da pasta:
#   app.R                 <- este arquivo (interface + servidor)
#   R/MTBD_SIR_SIRS.R     <- seu simulador (sem alteracoes)
#   R/replicatas.R        <- mc_MTBD_SIR_SIRS (replicatas + progresso)
#   R/graficos.R          <- trajetorias, fenogramas, arvores, resumo
# O Shiny carrega sozinho todos os arquivos da pasta R/.
#
# Para rodar:  shiny::runApp("caminho/para/app_mtbd")
# Pacotes:     shiny, bslib, ape, phytools  (+ showtext, opcional, p/ PDF em Montserrat)
# =====================================================================

library(shiny)
library(bslib)
library(ape)
library(phytools)

`%||%` <- function(a, b) if (is.null(a)) b else a

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
    timeIntro2  = numero(input$timeIntro2, "introdução da linhagem 2", 0),
    N           = if (modelo == "BD") 0 else numero(input$N, "N (população)", 2),
    timesky     = timesky,
    lambdasky_1 = vetor_sky("lambda1", "λ linhagem 1"),
    lambdasky_2 = vetor_sky("lambda2", "λ linhagem 2"),
    deathsky_1  = vetor_sky("death1",  "δ linhagem 1"),
    deathsky_2  = vetor_sky("death2",  "δ linhagem 2"),
    sampprobsky = vetor_sky("samp",    "prob. de amostragem", 0, 1),
    omegasky    = if (modelo == "SIRS") vetor_sky("omega", "ω") else rep(0, k),
    aceitar_extincao = isTRUE(input$aceitar_extincao),
    limite_seg  = if (is.na(input$limite) || input$limite <= 0) Inf else input$limite
  )

  if (p$rho > 1) stop("rho é uma probabilidade (0 a 1).", call. = FALSE)
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
    showNotification("timeIntro2 ≥ timestop: a linhagem 2 não chega a entrar.",
                     type = "warning", duration = 6)
  p
}

# Fonte do PDF: Montserrat via showtext (baixa do Google 1x); senao, padrao
.fonte <- new.env()
fonte_pdf <- function() {
  if (is.null(.fonte$familia)) {
    ok <- requireNamespace("showtext", quietly = TRUE) &&
          requireNamespace("sysfonts", quietly = TRUE)
    if (ok) ok <- tryCatch({
      sysfonts::font_add_google("Montserrat", "Montserrat"); TRUE
    }, error = function(e) FALSE)
    .fonte$familia <- if (ok) "Montserrat" else ""
  }
  .fonte$familia
}

salva_pdf <- function(arquivo, largura, altura, desenha) {
  fam <- fonte_pdf()
  grDevices::pdf(arquivo, width = largura, height = altura, bg = "transparent")
  if (nzchar(fam)) showtext::showtext_begin()
  on.exit({
    if (nzchar(fam)) showtext::showtext_end()
    grDevices::dev.off()
  })
  desenha(fam)
}

# ---- interface --------------------------------------------------------

ajuda <- function(...) tags$small(class = "text-muted d-block mb-2", ...)

barra_lateral <- sidebar(
  width = 340,
  radioButtons("model", "Modelo", c("SIR", "SIRS", "BD"),
               selected = "SIRS", inline = TRUE),
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
      numericInput("timeIntro2", "Introdução da linhagem 2 (tempo)", 5, min = 0)
    ),
    accordion_panel(
      "Taxas",
      ajuda("Um valor, ou um por intervalo separados por ponto-e-vírgula (ex.: 1.5; 0.8)."),
      layout_columns(
        textInput("lambda1", "λ linhagem 1", "1.5"),
        textInput("lambda2", "λ linhagem 2", "2")
      ),
      layout_columns(
        textInput("death1", "δ linhagem 1", "1"),
        textInput("death2", "δ linhagem 2", "1")
      ),
      ajuda("δ = remoção total (morte + amostragem)."),
      layout_columns(
        textInput("samp", "Prob. amostragem", "0.3"),
        conditionalPanel("input.model == 'SIRS'", textInput("omega", "ω (perda imunidade)", "0.1"))
      ),
      textInput("timesky", "Tempos de mudança (skyline)", "", placeholder = "vazio = taxas constantes"),
      uiOutput("r0")
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

ui <- page_sidebar(
  title = "MTBD de duas linhagens · SIR / SIRS / BD",
  theme = bs_theme(version = 5, primary = "#073A70",
                   base_font = font_collection("Montserrat", "Segoe UI", "Roboto",
                                               "Helvetica Neue", "Arial", "sans-serif")),
  sidebar = barra_lateral,
  navset_card_tab(
    id = "abas",
    aba("Trajetórias",
        tagList(
          selectInput("variavel", "Variável", VARIAVEIS_TRAJ),
          sliderInput("alpha_traj", "Opacidade", 0.05, 1, 0.4, step = 0.05),
          checkboxInput("mediana", "Mostrar mediana", FALSE),
          downloadButton("pdf_traj", "PDF", class = "btn-sm w-100")
        ),
        plotOutput("g_traj", height = "520px")),
    aba("Fenogramas",
        tagList(
          numericInput("bm_seed", "Semente do caráter", 1, step = 1),
          numericInput("bm_sig2", "σ² (taxa BM)", 1, min = 0, step = 0.1),
          radioButtons("bm_internos", "Nós internos",
                       c("Reconstruídos (fastAnc)" = "fastAnc",
                         "Simulados (verdadeiros)" = "simulados")),
          sliderInput("alpha_feno", "Opacidade", 0.05, 1, 0.4, step = 0.05),
          downloadButton("pdf_feno", "PDF", class = "btn-sm w-100")
        ),
        plotOutput("g_feno", height = "820px")),
    aba("Árvores",
        tagList(
          selectInput("rep_arvore", "Replicata", choices = NULL),
          ajuda("Eixo em tempo absoluto; o tracejado marca a entrada da linhagem 2."),
          downloadButton("nexus", "Todas as árvores (.nex)", class = "btn-sm w-100")
        ),
        plotOutput("g_arvores", height = "680px")),
    nav_panel("Resumo",
              uiOutput("caixas"),
              div(class = "d-flex justify-content-end mb-2",
                  downloadButton("rds", "Baixar resultados (.rds)", class = "btn-sm")),
              tableOutput("tabela")),
    nav_panel("Log", verbatimTextOutput("log"))
  )
)

# ---- servidor ---------------------------------------------------------

server <- function(input, output, session) {

  dados <- reactiveVal(NULL)

  # R0 do primeiro intervalo (SIR/SIRS: λ/δ; BD: λ/δ = razao nascimento/remocao)
  output$r0 <- renderUI({
    l1 <- tryCatch(le_vetor(input$lambda1, "")[1], error = function(e) NA)
    l2 <- tryCatch(le_vetor(input$lambda2, "")[1], error = function(e) NA)
    d1 <- tryCatch(le_vetor(input$death1,  "")[1], error = function(e) NA)
    d2 <- tryCatch(le_vetor(input$death2,  "")[1], error = function(e) NA)
    if (anyNA(c(l1, l2, d1, d2)) || d1 <= 0 || d2 <= 0) return(NULL)
    rotulo <- if (input$model == "BD") "λ/δ" else "R₀"
    ajuda(sprintf("%s (1º intervalo): L1 = %.2f · L2 = %.2f", rotulo, l1 / d1, l2 / d2))
  })

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

    n_ok <- sum(res$status == "ok")
    showNotification(
      sprintf("%d de %d replicatas ok · %.1f s", n_ok, res$re, res$elapsed_total),
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
    showNotification(sprintf("Carregado: %d replicatas (%d ok).",
                             length(obj$resultados), sum(obj$status == "ok")),
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
  valores <- reactive({
    sig2 <- if (is.na(input$bm_sig2) || input$bm_sig2 <= 0) 1 else input$bm_sig2
    simula_caracteres(info()$arvores, input$bm_seed, sig2, input$bm_internos)
  })

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
  desenha_feno <- function(fam = "") plot_fenogramas(
    info(), valores(), COR_FENO, input$alpha_feno, intro2(), xmax(), fam)

  output$g_traj <- renderPlot(desenha_traj(), res = 96)
  output$g_feno <- renderPlot(desenha_feno(), res = 96)
  output$g_arvores <- renderPlot({
    req(input$rep_arvore)
    r <- as.integer(input$rep_arvore)
    x <- pronto()$resultados[[r]]
    validate(need(is.list(x), "Replicata sem resultado."))
    plot_arvores_rep(x[[1]], COR_FENO, xmax(), intro2())
  }, res = 96)

  # ---- resumo ----
  output$caixas <- renderUI({
    d <- pronto()
    i <- info()
    tempo <- if (!is.null(d$elapsed_total)) sprintf("%.1f s", d$elapsed_total) else "—"
    layout_columns(
      fill = FALSE,
      value_box("Replicatas ok", sprintf("%d / %d", sum(d$status == "ok"), length(d$status))),
      value_box("Árvores (≥ 3 pontas)",
                sprintf("%d L1 · %d L2", sum(i$tipos == 1L), sum(i$tipos == 2L))),
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
  output$pdf_feno <- downloadHandler(
    filename = function() sprintf("fenogramas_%s.pdf", pronto()$model %||% "MC"),
    content  = function(file) salva_pdf(file, 8, 9, desenha_feno)
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
