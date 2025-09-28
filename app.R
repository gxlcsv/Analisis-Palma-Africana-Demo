# app.R
# ============================================================
# Palmatica — Shiny (sin raster pesado)
# - Capa principal: Palmatica_Promedio_Plantas.gpkg (POINT recomendado)
# - Overlays: Sub_Lote_Palmatica_Plantas.gpkg (opcional), Palmatica_Promedio.gpkg (opcional)
# - Fondo: Esri.WorldImagery
# - Clasificación: BAJO / MEDIO / ALTO
# - Extras:
#   * Zoom auto + botones "Zoom a capa" y "Zoom al sublote"
#   * Filtro: seleccionado=VERDE, no seleccionado=ROJO (CSS)
#   * Tabla: PLANT_ID solo 1 vez (primera) + encabezados MAYÚSCULAS; sin 'CLASS'
#   * Histogramas con cortes + valores (horizontal)
#   * Promedio de Lote bajo KPIs (renombrado)
#   * Umbrales manuales con fallback si faltan t1/t2
#   * SIN avisos jsonlite: evitamos vectores nombrados en fitBounds()
#   * KPI nuevos: ÁREA DEL LOTE (ha) y DENSIDAD (plantas/ha)
# ============================================================

# ---------- Paquetes ----------
req_pkgs <- c(
  "shiny","shinyWidgets","sf","dplyr","stringr","janitor",
  "leaflet","leaflet.extras","DT","ggplot2","tidyr","scales","readr","rsconnect"
)
to_install <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)

library(shiny); library(shinyWidgets)
library(sf); library(dplyr); library(stringr); library(janitor)
library(leaflet); library(leaflet.extras); library(DT)
library(ggplot2); library(tidyr); library(scales); library(readr)

options(shiny.maxRequestSize = 500*1024^2)
`%||%` <- function(a,b) if(!is.null(a)) a else b

# ---------- RUTA DE DATOS ----------
data_dir <- "C:/Users/Usuario/Documents/PalmaticaApp/data"
path_main <- file.path(data_dir, "Nombre_Promedio_Plantas.gpkg")
path_subl <- file.path(data_dir, "Sub_Lote_Nombre_Plantas.gpkg")
path_prom <- file.path(data_dir, "Nombre_Promedio.gpkg")

# Servir archivos del directorio de datos si hace falta (para el logo)
if (dir.exists(data_dir)) shiny::addResourcePath("pdata", data_dir)

# Helper para pintar el logo (intenta primero en ./www, luego en /data)
gxl_logo_tag <- function() {
  src <- NULL
  if (file.exists("www/Logo_GXL.jpg")) {
    src <- "Logo_GXL.jpg"              # servido automáticamente desde ./www
  } else if (file.exists(file.path(data_dir, "Logo_GXL.jpg"))) {
    src <- "pdata/Logo_GXL.jpg"        # servido desde /data
  }
  if (is.null(src)) return(NULL)
  tags$div(class = "gxl-logo",
           tags$img(src = src, alt = "Green Xpo Lab", title = "Green Xpo Lab",
                    style = "height:80px;"))
}

setwd("C:/Users/Usuario/Documents/PalmaticaApp")  # carpeta donde están PalmaApp.R y data/
# Si tu archivo principal es PalmaApp.R:
rsconnect::writeManifest(appPrimaryDoc = "app.R")

# (Alternativa) Si usas un wrapper app.R:
# rsconnect::writeManifest(appPrimaryDoc = "app.R")

# ---------- Utilitarios ----------
ensure_crs_4326 <- function(x) {
  if (is.na(sf::st_crs(x))) sf::st_crs(x) <- 4326
  sf::st_transform(x, 4326)
}
to_points <- function(sfobj) {
  gtypes <- unique(as.character(sf::st_geometry_type(sfobj)))
  if (length(gtypes) == 1 && gtypes %in% c("POINT","MULTIPOINT")) return(ensure_crs_4326(sfobj))
  tmp3857 <- sf::st_transform(sfobj, 3857)
  pts <- sf::st_point_on_surface(sf::st_geometry(tmp3857))
  sf_pts <- sfobj
  sf::st_geometry(sf_pts) <- sf::st_transform(sf::st_sfc(pts, crs = 3857), 4326)
  ensure_crs_4326(sf_pts)
}
guess_id_col <- function(df) {
  cand <- names(df); cand <- setdiff(cand, attr(df, "sf_column"))
  target <- cand[stringr::str_detect(cand, "(^plant_id$|^id$|id_|_id$|plant|planta|tree|codigo|cod|name)")]
  if (length(target) > 0) target[which.max(vapply(target, function(x) dplyr::n_distinct(df[[x]]), integer(1)))] else NULL
}
standardize_param_cols <- function(df) {
  nm <- names(df)
  map <- c(
    biomasa   = nm[stringr::str_detect(nm, "(?i)bioma|biomasa")][1]              %||% NA,
    clorofila = nm[stringr::str_detect(nm, "(?i)clor|chlor|chl")][1]             %||% NA,
    estres    = nm[stringr::str_detect(nm, "(?i)estr[eé]s|stress|psii|pri")][1]  %||% NA,
    nitrogeno = nm[stringr::str_detect(nm, "(?i)nitro|nitr[oó]geno|n_?index")][1]%||% NA
  )
  for (std in names(map)) if (!is.na(map[[std]]) && map[[std]] != std) {
    df <- dplyr::rename(df, !!std := dplyr::all_of(map[[std]]))
  }
  df
}
suggested_cutpoints_01 <- function(param) {
  switch(param,
         biomasa   = c(0.45, 0.65),
         clorofila = c(0.35, 0.55),
         nitrogeno = c(0.40, 0.60),
         estres    = c(0.40, 0.60),
         c(1/3, 2/3)
  )
}
backtransform_cuts <- function(x, t1_01, t2_01) {
  q <- quantile(x, probs = c(0.02, 0.98), na.rm = TRUE, names = FALSE)
  c(q[1] + t1_01*(q[2]-q[1]), q[1] + t2_01*(q[2]-q[1]))
}
classify_bma <- function(x, method = c("cuantiles","sugeridos","umbrales"), param = NULL, t1 = NA, t2 = NA) {
  method <- match.arg(method); labs <- c("BAJO","MEDIO","ALTO"); brks <- c(-Inf, NA, NA, Inf)
  if (method == "cuantiles") {
    v <- x[is.finite(x)]; if (length(v) < 3) return(factor(rep(NA_character_, length(x)), levels=labs, ordered=TRUE))
    q <- stats::quantile(v, probs = c(1/3, 2/3), na.rm = TRUE, names = FALSE, type = 7); brks[2:3] <- q
  } else if (method == "sugeridos") {
    cut01 <- suggested_cutpoints_01(param %||% ""); brks[2:3] <- backtransform_cuts(x, cut01[1], cut01[2])
  } else {
    if (!is.finite(t1) || !is.finite(t2)) {
      v <- x[is.finite(x)]
      if (length(v) >= 3) {
        q <- stats::quantile(v, probs = c(1/3, 2/3), na.rm = TRUE, names = FALSE, type = 7); brks[2:3] <- q
      } else {
        return(factor(rep(NA_character_, length(x)), levels=labs, ordered=TRUE))
      }
    } else brks[2:3] <- sort(c(t1, t2))
  }
  cut(x, breaks = brks, labels = labs, right = TRUE, include.lowest = TRUE, ordered_result = TRUE)
}
class_palette <- c("BAJO"="#d73027","MEDIO"="#fee08b","ALTO"="#1a9850")

pick_main_layer <- function(gpkg_path) {
  ls <- try(sf::st_layers(gpkg_path)$name, silent = TRUE)
  if (inherits(ls, "try-error") || length(ls)==0) stop("No se pudieron listar capas de: ", basename(gpkg_path))
  best <- ls[1]; best_score <- -1L
  for (ly in ls) {
    sf0 <- suppressWarnings(try(sf::st_read(gpkg_path, layer = ly, quiet = TRUE), silent = TRUE))
    if (inherits(sf0, "try-error")) next
    nm <- names(sf0); sc <- 0L
    sc <- sc + as.integer(any(stringr::str_detect(nm, "(?i)bioma|biomasa")))
    sc <- sc + as.integer(any(stringr::str_detect(nm, "(?i)clor|chlor|chl")))
    sc <- sc + as.integer(any(stringr::str_detect(nm, "(?i)estr[eé]s|stress|psii|pri")))
    sc <- sc + as.integer(any(stringr::str_detect(nm, "(?i)nitro|nitr[oó]geno|n_?index")))
    if (sc > best_score) { best <- ly; best_score <- sc }
  }
  best
}
pick_sublote_layer <- function(gpkg_path) {
  ls <- try(sf::st_layers(gpkg_path)$name, silent = TRUE)
  if (inherits(ls, "try-error") || length(ls)==0) return(NULL)
  pick <- NULL; poly_first <- NULL
  for (ly in ls) {
    sf0 <- suppressWarnings(try(sf::st_read(gpkg_path, layer = ly, quiet = TRUE), silent = TRUE))
    if (inherits(sf0, "try-error")) next
    gtypes <- unique(as.character(sf::st_geometry_type(sf0)))
    if (any(gtypes %in% c("POLYGON","MULTIPOLYGON"))) { poly_first <- ly; break }
    if (is.null(pick)) pick <- ly
  }
  poly_first %||% pick
}

# ---------- Carga de datos ----------
status_msgs <- c()
if (!file.exists(path_main)) stop("No se encontró 'Palmatica_Promedio_Plantas.gpkg' en: ", data_dir)

main_layer_name <- pick_main_layer(path_main)
plantas_raw <- sf::st_read(path_main, layer = main_layer_name, quiet = TRUE) |> janitor::clean_names()
plantas_raw <- ensure_crs_4326(plantas_raw)
plantas_pts <- to_points(plantas_raw)

id_col <- guess_id_col(plantas_pts)
if (is.null(id_col)) plantas_pts <- plantas_pts |> dplyr::mutate(plant_id = as.character(dplyr::row_number())) else
  if (id_col != "plant_id") plantas_pts <- plantas_pts |> dplyr::mutate(plant_id = as.character(.data[[id_col]]))

plantas_pts <- standardize_param_cols(plantas_pts)
to_num <- function(v) suppressWarnings(as.numeric(v))
for (c in c("biomasa","clorofila","estres","nitrogeno")) if (c %in% names(plantas_pts) && !is.numeric(plantas_pts[[c]])) plantas_pts[[c]] <- to_num(plantas_pts[[c]])

plantas_pts <- plantas_pts |> dplyr::arrange(plant_id) |> dplyr::group_by(plant_id) |> dplyr::slice(1) |> dplyr::ungroup()
status_msgs <- c(status_msgs, sprintf("Capa principal: %s (registros=%s)", main_layer_name, nrow(plantas_pts)))

sublote <- NULL; sub_layer_name <- NULL
if (file.exists(path_subl)) {
  sub_layer_name <- pick_sublote_layer(path_subl)
  if (!is.null(sub_layer_name)) {
    sublote <- ensure_crs_4326(sf::st_read(path_subl, layer = sub_layer_name, quiet = TRUE) |> janitor::clean_names())
    status_msgs <- c(status_msgs, sprintf("Sub Lote: %s (registros=%s)", sub_layer_name, nrow(sublote)))
  } else status_msgs <- c(status_msgs, "Sub Lote: no se halló capa utilizable.")
} else status_msgs <- c(status_msgs, "Sub Lote: archivo no encontrado (opcional).")

promedio_tbl <- NULL
if (file.exists(path_prom)) {
  pr_layers <- sf::st_layers(path_prom)$name
  if (length(pr_layers) > 0) {
    pr <- suppressWarnings(sf::st_read(path_prom, layer = pr_layers[1], quiet = TRUE))
    if (inherits(pr, "sf")) pr <- sf::st_drop_geometry(pr)
    promedio_tbl <- janitor::clean_names(pr)
    status_msgs <- c(status_msgs, sprintf("Promedio: %s (registros=%s)", pr_layers[1], nrow(promedio_tbl)))
  } else status_msgs <- c(status_msgs, "Promedio: archivo sin capas.")
} else status_msgs <- c(status_msgs, "Promedio: archivo no encontrado (opcional).")

# bbox inicial (como NUMÉRICOS sin nombres para evitar jsonlite warnings)
get_initial_bbox <- function() {
  if (nrow(plantas_pts) > 0) return(sf::st_bbox(ensure_crs_4326(plantas_pts)))
  if (!is.null(sublote) && nrow(sublote) > 0) return(sf::st_bbox(ensure_crs_4326(sublote)))
  c(xmin = -80, ymin = -10, xmax = -70, ymax = 5)
}
initial_bbox <- get_initial_bbox()
bbox_vals <- function(bb) {
  c(
    as.numeric(bb["xmin"]),
    as.numeric(bb["ymin"]),
    as.numeric(bb["xmax"]),
    as.numeric(bb["ymax"])
  )
}

# ---------- Métricas de Lote: Área (ha) y Densidad (plantas/ha) ----------
# - Área basada en sublote (geodésica con s2) → m2 → ha
# - Densidad: plantas dentro del sublote / área_ha
lote_stats <- local({
  if (!is.null(sublote) && nrow(sublote) > 0) {
    sub_union <- sf::st_union(sublote)
    area_m2 <- as.numeric(sf::st_area(sub_union))
    area_ha <- area_m2 / 10000
    if (is.finite(area_ha) && area_ha > 0) {
      inters <- sf::st_intersects(plantas_pts, sub_union, sparse = TRUE)
      n_pl_in <- sum(lengths(inters) > 0)
      dens <- n_pl_in / area_ha
      list(area_ha = area_ha, n_pl_sublote = n_pl_in, dens_ha = dens)
    } else {
      list(area_ha = NA_real_, n_pl_sublote = NA_integer_, dens_ha = NA_real_)
    }
  } else {
    list(area_ha = NA_real_, n_pl_sublote = NA_integer_, dens_ha = NA_real_)
  }
})

# ---------- UI ----------
ui <- fluidPage(
  tags$head(tags$style(HTML("
    .leaflet-container { background: #f8f9fb; }
    .value-box { padding: 6px 10px; border-radius: 12px; background: #f5f7fa; margin-bottom: 8px; }
    .small { font-size: 12px; color: #555; }

    /* Filtro: no seleccionado ROJO, seleccionado VERDE */
    #filtro_cat .btn { background-color:#e74c3c; border-color:#e74c3c; color:white; }
    #filtro_cat .btn.active { background-color:#1a9850; border-color:#1a9850; color:white; }

    /* Logo fijo esquina superior derecha */
    .gxl-logo {
      position: fixed;
      top: 10px;
      right: 12px;
      z-index: 1000;
      background: rgba(255,255,255,.85);
      padding: 6px 8px;
      border-radius: 8px;
      box-shadow: 0 1px 4px rgba(0,0,0,.2);
    }
    .gxl-logo img { display:block; height:40px; }
    
    /* Chips de conteo por clase */
    .chips { display:flex; gap:8px; flex-wrap:wrap; margin-top:6px; }
    .chip  { padding:4px 10px; border-radius:999px; font-weight:700; font-size:12px; }
    .chip.verde    { background:#1a9850; color:#fff; }
    .chip.amarillo { background:#fee08b; color:#333; }
    .chip.rojo     { background:#d73027; color:#fff; }

  "))),
  # Logo (usa gxl_logo_tag() definido en el paso 1)
  gxl_logo_tag(),
  
  # Título actualizado
  titlePanel("Análisis de Palma Africana - Demostración (Visualización y Clasificación)"),
  
  sidebarLayout(
    sidebarPanel(width = 4,
                 # (eliminado: Ruta de datos / Estado de carga)
                 
                 awesomeRadio(
                   "metodo", "Método de clasificación",
                   choices = list(
                     "Umbrales sugeridos (Palma Africana)" = "sugeridos",
                     "Cuantiles (tercios)" = "cuantiles",
                     "Umbrales manuales" = "umbrales"
                   ),
                   selected = "sugeridos"
                 ),
                 conditionalPanel("input.metodo == 'umbrales'",
                                  helpText("Defina t1 (BAJO→MEDIO) y t2 (MEDIO→ALTO). Si faltan, se usarán cuantiles temporalmente."),
                                  numericInput("t_bio_1","Biomasa t1", value = NA),
                                  numericInput("t_bio_2","Biomasa t2", value = NA),
                                  numericInput("t_clo_1","Clorofila t1", value = NA),
                                  numericInput("t_clo_2","Clorofila t2", value = NA),
                                  numericInput("t_est_1","Estrés t1", value = NA),
                                  numericInput("t_est_2","Estrés t2", value = NA),
                                  numericInput("t_nit_1","Nitrógeno t1", value = NA),
                                  numericInput("t_nit_2","Nitrógeno t2", value = NA)
                 ),
                 pickerInput("param_color", "Parámetro a colorear",
                             choices = list("Índice compuesto","Biomasa","Clorofila","Estrés","Nitrógeno"),
                             selected = "Índice compuesto"),
                 checkboxGroupButtons(
                   "filtro_cat", "Filtrar clases",
                   choices = list("ALTO","MEDIO","BAJO"),
                   selected = c("ALTO","MEDIO","BAJO"),
                   status = "default", size = "xs", justified = TRUE
                 ),
                 fluidRow(
                   column(6, actionButton("zoom_plants", "Zoom a capa", class = "btn btn-primary btn-sm", width = "100%")),
                   column(6, actionButton("zoom_sublote","Zoom al sublote", class = "btn btn-secondary btn-sm", width = "100%"))
                 ),
                 hr(),
                 downloadButton("dl_tabla","Descargar tabla (.csv)", class = "btn-success")
    ),
    mainPanel(width = 8,
              tabsetPanel(id = "tabs",
                          tabPanel("Mapa",
                                   leafletOutput("mapa", height = 650),
                                   br(),
                                   fluidRow(
                                     column(4, div(class="value-box",
                                                   strong("Total Plantas: "), textOutput("n_total", inline = TRUE))),
                                     column(4, div(class="value-box",
                                                   strong("Área del Lote (ha): "), textOutput("area_ha_txt", inline = TRUE))),
                                     column(4, div(class="value-box",
                                                   strong("Densidad (Plantas/ha): "), textOutput("dens_ha_txt", inline = TRUE)))
                                   ),
                                   fluidRow(
                                     column(12, div(class="value-box",
                                                    strong("Conteo por clase: "),
                                                    uiOutput("conteo_cats")))
                                   ),
                                   plotOutput("bar100_clases", height = 140),
                                   br(),
                                   uiOutput("titulo_prom_top"),
                                   tableOutput("prom_resumen")
                          ),
                          tabPanel("Tabla", DTOutput("tabla")),
                          tabPanel("Histogramas",
                                   helpText("Histogramas del parámetro seleccionado con cortes (se muestran los valores de corte)."),
                                   plotOutput("hist_param", height = 380),
                                   br(),
                                   helpText("Distribución del índice compuesto (promedio de clases mapeadas a 1/2/3)."),
                                   plotOutput("hist_comp", height = 340)
                          )
              )
    )
  )
)

# ---------- SERVER ----------
server <- function(input, output, session) {
  
  # ---- Clasificación ----
  r_classified <- reactive({
    df <- plantas_pts
    clasif <- function(v, p) {
      if (!(p %in% names(df))) return(factor(NA, levels=c("BAJO","MEDIO","ALTO"), ordered=TRUE))
      if (input$metodo == "cuantiles") classify_bma(v, "cuantiles")
      else if (input$metodo == "sugeridos") classify_bma(v, "sugeridos", param = p)
      else {
        t1 <- switch(p, biomasa=input$t_bio_1, clorofila=input$t_clo_1, estres=input$t_est_1, nitrogeno=input$t_nit_1)
        t2 <- switch(p, biomasa=input$t_bio_2, clorofila=input$t_clo_2, estres=input$t_est_2, nitrogeno=input$t_nit_2)
        classify_bma(v, "umbrales", t1=t1, t2=t2)
      }
    }
    out <- df |>
      dplyr::mutate(
        clase_biomasa   = clasif(biomasa,   "biomasa"),
        clase_clorofila = clasif(clorofila, "clorofila"),
        clase_estres    = clasif(estres,    "estres"),
        clase_nitrogeno = clasif(nitrogeno, "nitrogeno")
      )
    map_score <- c(BAJO=1,MEDIO=2,ALTO=3)
    present <- c("clase_biomasa","clase_clorofila","clase_estres","clase_nitrogeno")
    mat <- do.call(cbind, lapply(present[present %in% names(out)], function(col) map_score[as.character(out[[col]])]))
    avg <- rowMeans(mat, na.rm = TRUE)
    comp <- cut(avg, breaks = c(1-1e-6,1.6667,2.3333,3+1e-6),
                labels = c("BAJO","MEDIO","ALTO"), include.lowest = TRUE, right = TRUE, ordered_result = TRUE)
    out$compuesto <- comp; out$compuesto_score <- avg
    out
  })
  
  r_param_sel <- reactive({
    sel <- input$param_color
    d <- r_classified()
    cat_col <- switch(sel,
                      "Biomasa"   = "clase_biomasa",
                      "Clorofila" = "clase_clorofila",
                      "Estrés"    = "clase_estres",
                      "Nitrógeno" = "clase_nitrogeno",
                      "Índice compuesto" = "compuesto"
    )
    d |> dplyr::mutate(.clase_mapa = .data[[cat_col]])
  })
  
  r_filtrado <- reactive({
    d <- r_param_sel()
    keep <- input$filtro_cat
    d |> dplyr::filter(is.na(.clase_mapa) | as.character(.clase_mapa) %in% keep)
  })
  
  # ---- KPIs ----
  output$n_total    <- renderText({ format(nrow(r_param_sel()), big.mark = ",") })
  output$area_ha_txt <- renderText({
    if (is.null(sublote) || nrow(sublote) == 0) return("—")
    number(lote_stats$area_ha, accuracy = 0.01, big.mark = ",")
  })
  output$dens_ha_txt <- renderText({
    if (is.null(sublote) || nrow(sublote) == 0) return("—")
    number(lote_stats$dens_ha, accuracy = 0.01, big.mark = ",")
  })
  
  # Chips de colores con conteos y porcentajes
  output$conteo_cats <- renderUI({
    d <- r_param_sel()
    if (nrow(d) == 0) return(tags$div(class="chips"))  # vacío
    
    tb <- table(as.character(d$.clase_mapa))
    n_alto  <- if ("ALTO"  %in% names(tb)) as.integer(tb[["ALTO"]])  else 0L
    n_medio <- if ("MEDIO" %in% names(tb)) as.integer(tb[["MEDIO"]]) else 0L
    n_bajo  <- if ("BAJO"  %in% names(tb)) as.integer(tb[["BAJO"]])  else 0L
    tot     <- n_alto + n_medio + n_bajo
    
    p <- function(n) if (tot > 0) sprintf("%.1f%%", 100*n/tot) else "0.0%"
    
    tags$div(class="chips",
             tags$span(class="chip verde",    sprintf("ALTO: %s (%s)",  format(n_alto,  big.mark=","), p(n_alto))),
             tags$span(class="chip amarillo", sprintf("MEDIO: %s (%s)", format(n_medio, big.mark=","), p(n_medio))),
             tags$span(class="chip rojo",     sprintf("BAJO: %s (%s)",  format(n_bajo,  big.mark=","), p(n_bajo)))
    )
  })
  
  # Barra 100% apilada (en vez del pastel)
  output$bar100_clases <- renderPlot({
    d <- r_param_sel()
    if (nrow(d) == 0) return(NULL)
    df <- data.frame(clase = as.character(d$.clase_mapa))
    df <- df[!is.na(df$clase), , drop = FALSE]
    if (!nrow(df)) return(NULL)
    
    df$clase <- factor(df$clase, levels = c("BAJO","MEDIO","ALTO"))
    df <- df |> dplyr::count(clase) |> dplyr::arrange(clase)
    
    ggplot(df, aes(x = 1, y = n, fill = clase)) +
      geom_col(width = 0.6, position = "fill") +
      scale_fill_manual(values = class_palette) +
      coord_flip() +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      geom_text(aes(label = scales::percent(n/sum(n), accuracy = 0.1)),
                position = position_fill(vjust = 0.5), size = 4) +
      labs(x = NULL, y = NULL, title = "Distribución por clase (100%)") +
      theme_minimal(base_size = 12) +
      theme(axis.text.y = element_blank(),
            axis.ticks  = element_blank(),
            panel.grid  = element_blank(),
            legend.position = "none")
  })
  
  
  # ---- Promedio de Lote (en MAPA) ----
  prettify_prom <- function(tbl) {
    if (is.null(tbl) || nrow(tbl) == 0) return(NULL)
    out <- tbl
    if ("class" %in% names(out)) out <- dplyr::select(out, -class)
    rn_map <- c(
      ndre_prommean = "Clorofila",
      gndvi_prommean = "Nitrogeno",
      ndvi_prommean = "Estres",
      osavi_prommean = "Biomasa"
    )
    for (k in names(rn_map)) if (k %in% names(out)) out <- dplyr::rename(out, !!rn_map[[k]] := !!sym(k))
    out
  }
  prom_pretty <- reactive({ prettify_prom(promedio_tbl) })
  output$titulo_prom_top <- renderUI({
    if (is.null(prom_pretty())) return(NULL)
    tags$div(tags$b("Promedio de Lote"))
  })
  output$prom_resumen <- renderTable({
    tbl <- prom_pretty(); if (is.null(tbl)) return(NULL)
    cols_pref <- c("Biomasa","Clorofila","Estres","Nitrogeno")
    cols_show <- intersect(cols_pref, names(tbl))
    if (length(cols_show) == 0) return(tbl[1, , drop = FALSE])
    as.data.frame(tbl[, cols_show, drop = FALSE])
  }, striped = TRUE, spacing = "xs")
  
  # ---- Tabla principal ----
  output$tabla <- renderDT({
    d <- r_classified() |> sf::st_drop_geometry()
    if ("class" %in% names(d)) d <- d |> dplyr::select(-class)
    if (!"plant_id" %in% names(d)) d$plant_id <- seq_len(nrow(d))
    d <- d |> dplyr::relocate(plant_id, .before = 1) |>
      dplyr::rename(
        `Clase Biomasa`   = clase_biomasa,
        `Clase Clorofila` = clase_clorofila,
        `Clase Estrés`    = clase_estres,
        `Clase Nitrógeno` = clase_nitrogeno,
        `Clase compuesta` = compuesto,
        `Score compuesto` = compuesto_score
      )
    dup_idx <- duplicated(toupper(names(d))); d <- d[, !dup_idx, drop = FALSE]
    names(d) <- toupper(names(d))
    if ("PLANT_ID" %in% names(d)) {
      d <- d |> dplyr::relocate(PLANT_ID, .before = 1)
      dup_idx2 <- duplicated(names(d)); d <- d[, !dup_idx2, drop = FALSE]
    }
    DT::datatable(d, rownames = FALSE, filter = "top",
                  options = list(scrollX = TRUE, pageLength = 15))
  })
  output$dl_tabla <- downloadHandler(
    filename = function() sprintf("palmatica_clasificacion_%s.csv", format(Sys.Date(), "%Y%m%d")),
    content  = function(file) {
      d <- r_classified() |> sf::st_drop_geometry()
      if ("class" %in% names(d)) d <- d |> dplyr::select(-class)
      readr::write_csv(d, file)
    }
  )
  
  # ---- Histogramas ----
  # --- REEMPLAZA output$hist_param POR ESTE ---
  output$hist_param <- renderPlot({
    d <- r_param_sel(); sel <- input$param_color
    vcol <- switch(sel,
                   "Biomasa"="biomasa","Clorofila"="clorofila",
                   "Estrés"="estres","Nitrógeno"="nitrogeno",
                   "Índice compuesto"="compuesto_score")
    if (!vcol %in% names(d)) return(NULL)
    v_draw <- d[[vcol]]; v_draw <- v_draw[is.finite(v_draw)]
    if (length(v_draw) < 3) return(NULL)
    
    # Cortes sobre todos los datos del parámetro (estable)
    all_df <- r_classified()
    v_all <- all_df[[vcol]]; v_all <- v_all[is.finite(v_all)]
    if (length(v_all) < 3 && sel != "Índice compuesto") return(NULL)
    
    # Determinar t1/t2 y etiquetas
    if (sel == "Índice compuesto") {
      cuts <- c(1.6667, 2.3333); lab1 <- "BAJO→MEDIO"; lab2 <- "MEDIO→ALTO"
    } else if (input$metodo == "cuantiles") {
      cuts <- as.numeric(quantile(v_all, probs = c(1/3, 2/3), na.rm = TRUE))
      lab1 <- "BAJO→MEDIO"; lab2 <- "MEDIO→ALTO"
    } else if (input$metodo == "sugeridos") {
      param <- switch(sel,"Biomasa"="biomasa","Clorofila"="clorofila","Estrés"="estres","Nitrógeno"="nitrogeno")
      cut01 <- suggested_cutpoints_01(param); cuts <- backtransform_cuts(v_all, cut01[1], cut01[2])
      lab1 <- "BAJO→MEDIO (sugerido)"; lab2 <- "MEDIO→ALTO (sugerido)"
    } else {
      t1 <- switch(sel, Biomasa=input$t_bio_1, Clorofila=input$t_clo_1, `Estrés`=input$t_est_1, `Nitrógeno`=input$t_nit_1)
      t2 <- switch(sel, Biomasa=input$t_bio_2, Clorofila=input$t_clo_2, `Estrés`=input$t_est_2, `Nitrógeno`=input$t_nit_2)
      if (!is.finite(t1) || !is.finite(t2)) {
        cuts <- as.numeric(quantile(v_all, probs = c(1/3, 2/3), na.rm = TRUE))
        lab1 <- "BAJO→MEDIO (cuantiles por falta de t1/t2)"; lab2 <- "MEDIO→ALTO (cuantiles)"
      } else {
        cuts <- sort(c(t1,t2)); lab1 <- "BAJO→MEDIO (manual)"; lab2 <- "MEDIO→ALTO (manual)"
      }
    }
    
    # Binéo manual para colorear por umbrales
    hinfo <- hist(v_draw, breaks = 30, plot = FALSE)
    df <- data.frame(
      xmin  = head(hinfo$breaks, -1),
      xmax  = tail(hinfo$breaks, -1),
      count = hinfo$counts
    )
    df$xmid <- (df$xmin + df$xmax)/2
    df$clase_bin <- cut(df$xmid,
                        breaks = c(-Inf, cuts[1], cuts[2], Inf),
                        labels = c("BAJO","MEDIO","ALTO"),
                        right = TRUE, ordered_result = TRUE)
    
    ymax <- max(df$count); if (!is.finite(ymax)) ymax <- 1
    fmt  <- scales::label_number(accuracy = 0.01, big.mark = ",")
    
    ggplot(df) +
      geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = count, fill = clase_bin), color = NA) +
      scale_fill_manual(values = class_palette) +
      labs(x = sel, y = "Frecuencia", title = paste("Histograma —", sel)) +
      theme_minimal(base_size = 12) +
      guides(fill = "none") +
      geom_vline(xintercept = cuts, linetype = "dashed") +
      annotate("text", x = cuts[1], y = ymax*0.92,
               label = paste0(lab1, ": ", fmt(cuts[1])),
               angle = 0, vjust = -0.5, size = 3.7) +
      annotate("text", x = cuts[2], y = ymax*0.84,
               label = paste0(lab2, ": ", fmt(cuts[2])),
               angle = 0, vjust = -0.5, size = 3.7)
  })
  
  # --- REEMPLAZA output$hist_comp POR ESTE ---
  output$hist_comp <- renderPlot({
    v <- r_classified()$compuesto_score; v <- v[is.finite(v)]
    if (length(v) < 3) return(NULL)
    cuts <- c(1.6667, 2.3333)
    
    hinfo <- hist(v, breaks = 30, plot = FALSE)
    df <- data.frame(
      xmin  = head(hinfo$breaks, -1),
      xmax  = tail(hinfo$breaks, -1),
      count = hinfo$counts
    )
    df$xmid <- (df$xmin + df$xmax)/2
    df$clase_bin <- cut(df$xmid,
                        breaks = c(-Inf, cuts[1], cuts[2], Inf),
                        labels = c("BAJO","MEDIO","ALTO"),
                        right = TRUE, ordered_result = TRUE)
    
    ymax <- max(df$count); if (!is.finite(ymax)) ymax <- 1
    
    ggplot(df) +
      geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = count, fill = clase_bin), color = NA) +
      scale_fill_manual(values = class_palette) +
      labs(x = "Score compuesto (1=BAJO, 2=MEDIO, 3=ALTO)", y = "Frecuencia", title = "Histograma — Índice compuesto") +
      theme_minimal(base_size = 12) +
      guides(fill = "none") +
      geom_vline(xintercept = cuts, linetype = "dashed") +
      annotate("text", x = 1.6667, y = ymax*0.92, label = "BAJO→MEDIO: 1.6667", angle = 0, vjust = -0.5, size = 3.7) +
      annotate("text", x = 2.3333, y = ymax*0.84, label = "MEDIO→ALTO: 2.3333", angle = 0, vjust = -0.5, size = 3.7)
  })
  
  # ---- Mapa ----
  output$mapa <- renderLeaflet({
    bounds <- bbox_vals(initial_bbox)  # NUMÉRICOS sin nombres
    leaflet() |>
      addProviderTiles(providers$Esri.WorldImagery, group = "Satélite") |>
      addProviderTiles(providers$OpenStreetMap, group = "OSM") |>
      addScaleBar(position = "bottomleft") |>
      addLayersControl(baseGroups = c("Satélite","OSM"),
                       overlayGroups = c("Sub Lote","Plantas"),
                       options = layersControlOptions(collapsed = FALSE)) |>
      fitBounds(bounds[1], bounds[2], bounds[3], bounds[4])
  })
  
  # Sub Lote
  observe({
    leafletProxy("mapa") |> clearGroup("Sub Lote")
    if (is.null(sublote) || nrow(sublote)==0) return()
    sub <- ensure_crs_4326(sublote)
    gtypes <- unique(as.character(sf::st_geometry_type(sub)))
    if (any(gtypes %in% c("POLYGON","MULTIPOLYGON"))) {
      leafletProxy("mapa") |> addPolygons(data = sub, group = "Sub Lote", fill = FALSE, weight = 2, color = "#2b8cbe")
    } else if (any(gtypes %in% c("LINESTRING","MULTILINESTRING"))) {
      leafletProxy("mapa") |> addPolylines(data = sub, group = "Sub Lote", weight = 2, color = "#2b8cbe")
    } else {
      leafletProxy("mapa") |> addCircleMarkers(data = sub, group = "Sub Lote", radius = 3, color = "#2b8cbe",
                                               fill = TRUE, fillOpacity = 0.6)
    }
  })
  
  # Plantas (popup precalculado; sin label hover)
  observe({
    d <- r_filtrado(); req(nrow(d) > 0)
    add_col <- function(nm, lab) if (nm %in% names(d)) paste0("<b>", lab, ": </b>", ifelse(is.na(d[[nm]]),"—", scales::number(d[[nm]])), "<br/>") else ""
    d$popup_html <- paste0(
      "<b>PLANTA: </b>", d$plant_id, "<br/>",
      add_col("biomasa","BIOMASA"),
      add_col("clorofila","CLOROFILA"),
      add_col("estres","ESTRÉS"),
      add_col("nitrogeno","NITRÓGENO"),
      "<b>CLASE COMPUESTA: </b>", as.character(d$compuesto),
      " (", sprintf("%.2f", d$compuesto_score %||% NA_real_), ")"
    )
    pal <- colorFactor(class_palette, levels = c("BAJO","MEDIO","ALTO"))
    leafletProxy("mapa") |>
      clearGroup("Plantas") |>
      addCircleMarkers(
        data = to_points(d) |> sf::st_transform(4326),
        group = "Plantas",
        radius = 5, stroke = TRUE, weight = 1, color = "#333333",
        fillOpacity = 0.85, fillColor = ~pal(as.character(.clase_mapa)),
        popup = ~popup_html
      )
  })
  
  # Botón: Zoom a capa
  observeEvent(input$zoom_plants, {
    proxy <- leafletProxy("mapa")
    d <- r_filtrado()
    if (nrow(d) > 0) {
      bb <- sf::st_bbox(ensure_crs_4326(d)); vals <- bbox_vals(bb)
      proxy %>% fitBounds(vals[1], vals[2], vals[3], vals[4])
    } else if (nrow(plantas_pts) > 0) {
      bb <- sf::st_bbox(ensure_crs_4326(plantas_pts)); vals <- bbox_vals(bb)
      proxy %>% fitBounds(vals[1], vals[2], vals[3], vals[4])
    } else if (!is.null(sublote) && nrow(sublote) > 0) {
      bb <- sf::st_bbox(ensure_crs_4326(sublote)); vals <- bbox_vals(bb)
      proxy %>% fitBounds(vals[1], vals[2], vals[3], vals[4])
    } else {
      proxy %>% setView(lng = -74, lat = 4, zoom = 5)
    }
  })
  
  # Botón: Zoom al sublote
  observeEvent(input$zoom_sublote, {
    proxy <- leafletProxy("mapa")
    if (!is.null(sublote) && nrow(sublote) > 0) {
      bb <- sf::st_bbox(ensure_crs_4326(sublote)); vals <- bbox_vals(bb)
      proxy %>% fitBounds(vals[1], vals[2], vals[3], vals[4])
    } else {
      showNotification("No hay capa de Sub Lote cargada.", type = "warning")
    }
  })
}

shinyApp(ui, server)
