##################################################################
######### 2° Workshop: Experimentação na agricultura 4.0 #########
##################################################################

# Instruções:
# Baixar o repositório do github
# Instalar a versão de desenvolvimento do pliman com o comando abaixo
# Descomente para rodar

#devtools::install_github("TiagoOlivoto/pliman")

############################### Pacotes e diretório #####################
library(pliman)
library(tidyverse)
set_wd_here()


################### Severidade de doenças ##################
## Interativa
set_wd_here()
img <- image_import("sev_orange.jpg", plot = TRUE)
sev <- measure_disease_iter(img, viewer = "mapview")
sev$results


## Segmentar e contar lesões
leaves <- image_import("ideal.jpg", plot = TRUE) |> image_vertical()
sev <- measure_disease_iter(leaves,
                            show_features = TRUE,
                            lesion_size = "small",
                            show_segmentation = TRUE,
                            viewer = "mapview",
                            watershed = TRUE,
                            save_image = TRUE)


## Utilizando índices de cores
# Severidade de doenças

set_wd_here("sevsoja")
img <- image_pliman("sev_leaf.jpg")
image_index(img, "B")
image_index(img, "NGRDI")

image_segment_iter(img,
                   index = c("B", "NGRDI"),
                   ncol = 3)

## Processamento em lote
sev_lote <- 
  measure_disease(pattern = "soy",
                  index_lb = "B",
                  index_dh = "NGRDI",
                  threshold = c("Otsu", -0.03),
                  plot =  FALSE,
                  save_image = TRUE,
                  dir_processed = "proc",
                  show_contour = FALSE,
                  parallel = TRUE,
                  col_lesions = "brown")

sev_lote$severity |> 
  ggplot(aes(x = symptomatic)) +
  geom_histogram(bins = 8)






# Bactérias
set_wd_here()
img <- image_import("bacterias.jpg")
plot(img)
image_index(img, "B")

res <-
  analyze_objects(img,
                  threshold = 0.42,
                  # object_size = "small",
                  tolerance = 0.1,
                  extension = 1,
                  # filter = 2,
                  marker = "point",
                  show_contour = FALSE,
                  index = "B")
res$statistics




# Waterpaper
set_wd_here("watterpaper")
img <- image_import("IMG_5600.jpeg")
cropped <- 
  image_crop(img,
             width = 415:1325,
             height = 355:1687,
             plot = TRUE)
# dpi 462
res <- 
  analyze_objects(cropped, 
                  threshold = -0.3,
                  lower_noise = 0,
                  index = "B-R",
                  marker = "point",
                  marker_col = "red",
                  show_contour = FALSE,
                  invert = TRUE)
meas <- get_measures(res, dpi = 46.2)
dfplot <- 
  meas |> 
  as_tibble() |> 
  select(diam_mean, perimeter, area) |> 
  pivot_longer(everything())



ggplot(dfplot, aes(x = value)) +
  geom_histogram(bins = 10,
                 fill = "#921614") +
  facet_wrap(~name, scales = "free") +
  theme_minimal(base_size = 18) +
  labs(x = "Valor observado (mm)",
       y = "Número de gotas")





################# Fenotipagem de alto rendimento   ##############
mosaic <- mosaic_input("lettuce.tif")
# indexes to be computed
indexes <- c("NGRDI", "GLI", "SCI", "BI", "VARI", "MGVRI")
# index used to segment plants
mosaic_index(mosaic, "GLI",
             r = 1, g = 2, b = 3)
# draw four blocks of 12 plots
# each block will be analyzed separately
an <- mosaic_analyze(mosaic,
                     r = 1, g = 2, b = 3,
                     nrow = 18,
                     ncol = 4,
                     segment_individuals = TRUE,
                     segment_index = "GLI",
                     plot_index = indexes)


# plot the results at plot and individual level
mosaic_view(mosaic,
            r = 1, g = 2, b = 3,
            shapefile = an$result_indiv,
            attribute = "mean.NGRDI")


res <- readRDS("res_lettuce.rds")
plotid <- df <- rio::import("plotid.xlsx")

b1 <- c(16, 11, 6, 3, 13, 17, 9, 15, 10, 5, 7, 2, 18, 4, 12, 14, 1, 8)
b2 <- c(13, 7, 18, 11, 4, 15, 1, 14, 9, 17, 5, 12, 2, 16, 8, 3, 10, 6)
b3 <- c(12, 6, 13, 2, 8, 1, 17, 4, 3, 10, 14, 9, 16, 11, 15, 5, 7, 18)
b4 <- c(2, 7, 14, 11, 4, 17, 10, 8, 16, 6, 13, 3, 9, 5, 18, 1, 12, 15)
res$plot_id <- leading_zeros(c(b1, b2, b3, b4), 4)

library(tidyverse)
library(AgroR)
res <-
  res |>
  left_join(plotid)

diametro <- with(res, FAT2DBC(inoculante, p, block, diam_mean))
p1 <-
  diametro$graph1 +
  labs(x = "Inoculação",
       y = "Diâmetro (cm)")
p2 <-
  diametro$graph2 +
  labs(x = "Fósforo",
       y = "Diâmetro (cm)")
p1 / p2




# epidemiologia temporal
library(epifitter)
library(tidyverse)
set.seed(1)
ep1 <-
  sim_logistic(N = 100,
               y0 = 0.001,
               dt = 5,
               r = 0.12,
               alpha = 0.2,
               n = 3)
ep2 <-
  sim_logistic(N = 100,
               y0 = 0.001,
               dt = 5,
               r = 0.23,
               alpha = 0.2,
               n = 3)


data <-
  data.frame(time =  ep1[,2],
             Cult1 = ep1[,4],
             Cult2 = ep2[,4]) |>
  mutate(time = as.numeric(time)) |>
  pivot_longer(-time,
               names_to = "Cultivar",
               values_to = "Severidade") |>
  group_by(time, Cultivar)


curva <-
  data %>%
  ggplot(aes(time, Severidade, color = Cultivar))+
  geom_jitter(width = 0.2) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Dias",
       y = "Severidade")

# Usando o ggplot2
formula <- as.formula("y ~ K * (1 + ((K - y0)/y0) * exp(-r * x))^-1")

# fit the model
# Package epifitter
modfit <-
  fit_multi(
    data = data,
    time_col = "time",
    intensity_col = "Severidade",
    strata_cols = "Cultivar",
    starting_par = list(y0 = 0.01, r = 0.03, K = 1),
    nlin = TRUE,
    estimate_K = TRUE,
  )$Parameters |>
  filter(model == "Logistic")
modfit


plot_curv <-
  ggplot(data, aes(time, Severidade, color = Cultivar))+
  geom_jitter(width = 0.2) +
  geom_smooth(method = "nls",
              data = data |> subset(Cultivar == "Cult1"),
              size = 1,
              method.args = list(formula = formula,
                                 start = c(y0 = 0.001,
                                           r = 0.11,
                                           K = 1)),
              se = FALSE) +
  geom_smooth(method = "nls",
              data = data |> subset(Cultivar == "Cult2"),
              size = 1,
              method.args = list(formula = formula,
                                 start = c(y0 = 0.001,
                                           r = 0.23,
                                           K = 1)),
              se = FALSE) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Dias da epidemia",
       y = "Severidade (%)")


D(expression(K * (1 + ((K - y0)/y0) * exp(-r * x))^-1), "x")
dy <- function(x,y0,K,r){
  K * ((1 + ((K - y0)/y0) * exp(-r * x))^-(1 + 1) * (((K - y0)/y0) *
                                                       (exp(-r * x) * r)))
}



plot_deriv <-
  ggplot() +
  stat_function(fun = dy,
                n = 400,
                aes(color = "Cult1"),
                size = 1,
                xlim = c(0, 100),
                args = c(y0 = modfit[1, 3],
                         r = modfit[1, 5],
                         K = modfit[1, 7])) +
  stat_function(fun = dy,
                n = 400,
                aes(color = "Cult2"),
                size = 1,
                xlim = c(0, 100),
                args = c(y0 = modfit[2, 3],
                         r = modfit[2, 5],
                         K = modfit[2, 7])) +
  labs(x = "Dias da epidemia",
       y = "Taxa de aumento da severidade (% por dia)",
       color = "Cultivar") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

library(patchwork)
plot_curv + plot_deriv

optimise(dy,
         c(0, 100),
         y0 = modfit[1, 3],
         r = modfit[1, 5],
         K = modfit[1, 7],
         maximum = TRUE)
optimise(dy,
         c(0, 100),
         y0 = modfit[2, 3],
         r = modfit[2, 5],
         K = modfit[2, 7],
         maximum = TRUE)






# DOSE-RESPOSTA
library(drda)
set_wd_here("dose-resposta")

# Computar a severidade por folha
sev <-
  measure_disease_byl(pattern = "img",
                      index = "B",
                      index_dh = "NGRDI",
                      parallel = TRUE,
                      opening = c(25, 0))
sev_ <-
  sev$severity |>
  separate(img, into = c("img", "produto", "dose"), sep = "_") |>
  mutate(dose = as.numeric(str_replace_all(dose, ",", ".")),
         symptomatic = symptomatic / 100)
rio::export(sev_, "drcdata.xlsx")

sev_ <- rio::import("drcdata.xlsx")




models <- 
  sev_ |> 
  group_by(produto) |> 
  nest() |> 
  mutate(models = map(data, 
                      ~drda(symptomatic ~ dose,
                            data = .,
                            mean_function = "ll4"))) |> # definir o modelo aqui
  dplyr::select(-data)

# função para obter os coeficientes
get_results <- function(model,
                        resplevel = 0.5,
                        type = "relative"){
  coefs <- coef(model) |> t()
  ed <- effective_dose(model, y = resplevel) |> as.data.frame()
  integ <- data.frame(nauc = nauc(model, range(model$model[[2]])))
  cbind(coefs, ed, integ)
}

# Obter os coeficientes
# alpha:  the value of the function at x = 0
# delta: height of the curve
# eta: the steepness (growth rate) of the curve
# phi: the x value at which the curve is equal to its mid-point

coefs <- 
  models |> 
  mutate(coefs = map_dfr(
    .x = models,
    .f = ~get_results(., resplevel = 0.5)) # DL50
  ) |> 
  dplyr::select(-models) |> 
  unnest(coefs) |> 
  ungroup() |> 
  as.data.frame()
coefs


plot(models$models[[1]], models$models[[2]],
     level = 0,
     base = "10",
     ylim = c(0, 0.5),
     xlim = c(0, 100),
     legend = c("P1", "P2"),
     xlab = "Dose (ppm)",
     ylab = "Severidade da doença",
     col = metan::ggplot_color(2),
     cex = 2)

models

# derivada em relação a dose do modelo
D(expression(alpha + delta * x^eta / (x^eta + phi^eta)), "x")

dy <- function(x,alpha,  delta,   eta,   phi){
  delta * (x^(eta - 1) * eta)/(x^eta + phi^eta) - delta * x^eta * 
    (x^(eta - 1) * eta)/(x^eta + phi^eta)^2
}

# Primeira derivada
ggplot(data.frame(x = c(0, 5)), aes(x = x)) +
  pmap(coefs |> select(produto:phi), function(produto, alpha, delta, eta, phi) {
    stat_function(fun = function(x) dy(x, alpha, delta, eta, phi),
                  aes(color = produto),
                  linewidth = 1)
  }) + 
  geom_vline(aes(xintercept = phi,
                 color = produto),
             data = coefs,
             linetype = 2) +
  labs(x = "Dose (ppm)",
       y = "Taxa de redução da severidade (% por ppm)",
       color = "Produto") +
  ggthemes::theme_base()



# Integração da área reduzida até o ponto médio * 2
# Igual ao parâmetro delta
integrate(dy,
          lower = 0,
          upper = 0.4335226,
          alpha = coefs$alpha[2],
          delta = coefs$delta[2],
          eta = coefs$eta[2],
          phi = coefs$phi[2])



# Espacial
## Análise FOCI
pak::pkg_install("emdelponte/r4pde")
foci <- tibble::tribble(
  ~x, ~`1`, ~`2`, ~`3`, ~`4`, ~`5`, ~`6`, ~`7`, ~`8`, ~`9`,
  1,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  2,   1,   1,   1,   0,   0,   0,   0,   1,   0,
  3,   1,   1,   1,   0,   0,   0,   1,   1,   1,
  4,   0,   1,   1,   0,   0,   0,   0,   1,   0,
  5,   0,   1,   1,   0,   0,   0,   0,   0,   0,
  6,   0,   0,   0,   1,   0,   0,   0,   0,   0,
  7,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  8,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  9,   0,   0,   0,   0,   0,   1,   0,   1,   0,
  10,   0,   0,   0,   0,   0,   0,   1,   0,   0,
  11,   0,   1,   0,   0,   0,   1,   0,   1,   0,
  12,   0,   0,   0,   0,   0,   0,   0,   0,   0
)

library(tidyr)

foci2 <- foci |> 
  pivot_longer(2:10, names_to = "y", values_to = "i")
foci2

library(ggplot2)
foci2 |> 
  ggplot(aes(x, y, fill = factor(i)))+
  geom_tile(color = "black")+
  scale_fill_manual(values = c("grey70", "darkred"))+
  theme_void()+
  coord_fixed()+
  theme(legend.position = "none")


foci2$y <- as.integer(foci2$y) # transform to numeric

library(r4pde)
result_foci <- AFSD(foci2)

plot_AFSD(result_foci[[3]])






# Autocorrelação espacial
# Necessário um objeto sf

library(spdep)
map <- sf::st_as_sf(sf::st_read(system.file("ex/lux.shp", package="terra")))
nb <- poly2nb(map, queen = TRUE) # queen shares point or border
nbw <- nb2listw(nb, style = "W", type = "queen", torus = FALSE)

# Global Moran's I
gmoran <- moran.test(map$POP, nbw)
lmoran <- localmoran(map$POP, nbw, alternative = "greater")
map$lmI <- lmoran[, "Ii"]
shapefile_view(map, attribute = "lmI")





# Índices de seleção
library(metan)

# simulate a data set
# 10 genotypes
# 5 replications
# 4 traits
df <-
  g_simula(ngen = 10,
           nrep = 5,
           nvars = 4,
           gen_eff = 35,
           seed = c(1, 2, 3, 4))

# run a mixed-effect model (genotype as random effect)
mod <-
  gamem(df,
        gen = GEN,
        rep = REP,
        resp = everything())
# BLUPs for genotypes
gmd(mod, "blupg")
mod <- mgidi(mod)
