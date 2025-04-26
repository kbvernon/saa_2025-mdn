library(arcgislayers)
library(colorspace)
library(magick)
library(patchwork)
library(sf)
library(terra)
library(tidyterra)
library(tidyverse)

# default theme ----------------------------------------------------------

theme_set(theme_bw(12))

theme_update(
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.grid.major = element_line(linewidth = 0.3, color = "#E9ECEF"),
  panel.grid.minor = element_blank(),
  text = element_text(family = "Public Sans")
)

coolors <- list(
  "darkblue" = "#2A324B",
  "yellow"   = "#FFD23F",
  "orange"   = "#CD5334",
  "purple"   = "#533747",
  "blue"     = "#28666E",
  "gray"     = "#829399"
)

coolors_lt <- lighten(coolors, 0.5)
coolors_dk <- darken(coolors, 0.5)

public_sans <- file.path(
  "https://fonts.googleapis.com",
  "css2?family=Public+Sans:ital,wght@0,100..900;1,100..900&display=swap"
)

# ggplot2 3.5.2.9000 now supports absolute panel sizing!
ph <- unit(3, "in")

get_dims <- function(x, panel_width, panel_height) {
  # get dimensions around plot panel in millimeters
  d <- patchwork::get_dim(x)

  # convert to inches and sum
  d <- lapply(unclass(d), \(z) ggplot2::unit(0.0393701 * sum(z), "in"))

  w <- panel_width + d[["l"]] + d[["r"]]
  h <- panel_height + d[["t"]] + d[["b"]]

  c(w, h)
}

# data -------------------------------------------------------------------

tree_rings <- read_csv("_data/tree-rings.csv")
room_counts <- read_csv("_data/room-counts.csv")
ceramics <- read_csv("_data/ceramic-presence.csv")
lookup <- read_csv("_data/ceramic-lookup.csv")

# four corners
furl <- file.path(
  "https://tigerweb.geo.census.gov",
  "arcgis/rest/services",
  "TIGERweb",
  "State_County",
  "MapServer/0"
) |> arc_open()

four_corners <- furl |>
  arc_select(
    fields = "NAME",
    where = "NAME IN ('Utah', 'Colorado', 'New Mexico', 'Arizona')"
  ) |>
  rename_with(tolower)

# probability distributions
region <- read_csv("_data/region.csv")
sites <- read_csv("_data/important-sites.csv")

# random sites -----------------------------------------------------------

N <- 15

set.seed(1701)

sites <- tibble(
  x = runif(N, min = 0, max = 10),
  y = runif(N, min = 0, max = 10)
)

s <- sites |>
  filter(x < 5, y > 8) |>
  slice(1)

p1 <- ggplot(sites, aes(x, y)) +
  geom_point(
    size = 3,
    shape = 21,
    fill = coolors_lt[["purple"]],
    color = coolors[["purple"]]
  ) +
  annotate(
    "text",
    x = s[["x"]] + 0.35,
    y = s[["y"]],
    hjust = 0,
    label = "s\u1D62",
    size = 18 / .pt
  ) +
  labs(
    x = "Easting",
    y = "Northing"
  ) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 10)) +
  theme(
    panel.heights = ph,
    panel.widths = ph,
    plot.margin = margin(r = 12)
  )

# sample gmm
rgmm <- function(n, theta) {
  # rmultinom returns a length(prob) x n matrix
  y <- apply(
    rmultinom(n, size = 1, prob = theta[["pi"]]),
    MARGIN = 2,
    FUN = \(.k) which(.k == 1)
  )

  rnorm(n, mean = theta[["mu"]][y], sd = theta[["sigma"]][y])
}

dgmm <- function(x, theta) {
  dN <- sapply(
    seq_along(theta[["pi"]]),
    FUN = \(.k) dnorm(x, mean = theta[["mu"]][.k], sd = theta[["sigma"]][.k])
  )
  colSums(theta[["pi"]] * t(dN))
}

chrono <- vector(mode = "list", length = N)

for (i in 1:N) {
  K <- rpois(1, 9)
  w <- rexp(K)

  n_mu <- floor(K / 3)

  mu <- c(
    rpois(n_mu, 50),
    rpois(n_mu, 100),
    rpois(K - (2 * n_mu), 150)
  )

  theta <- list(
    pi = w / sum(w),
    mu = mu,
    sigma = rgamma(K, 30, 3)
  )

  chrono[[i]] <- tibble(
    s = i,
    x = 1:200,
    y = dgmm(x, theta)
  )

  rm(i, K, w, n_mu, mu, theta)
}

chrono <- bind_rows(chrono)

chrono <- bind_rows(
  chrono |> mutate(type = "Site Chronology"),
  chrono |> group_by(x) |> summarize(
    type = "Regional Chronology",
    s = 999,
    y = mean(y)
  )
) |> mutate(type = fct_relevel(factor(type), "Site Chronology"))

p2 <- ggplot() +
  geom_line(
    data = chrono,
    aes(x, y, group = s, color = type, linewidth = type),
    lineend = "round"
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Site Chronology" = alpha(coolors_lt[["purple"]], 0.35),
      "Regional Chronology" = coolors[["orange"]]
    )
  ) +
  scale_linewidth_manual(
    name = NULL,
    values = c(
      "Site Chronology" = 0.3,
      "Regional Chronology" = 2
    )
  ) +
  labs(
    x = "Time",
    y = "Probability Density"
  ) +
  theme(
    legend.position = "none",
    panel.heights = ph,
    panel.widths = 1.3 * ph,
    plot.margin = margin(l = 12)
  )

free(p1) + free(p2)

dims <- get_dims(p1, ph, ph) + get_dims(p2, 1.3 * ph, ph)
dims[2] <- dims[2] / 2

fn <- "figures/goal.svg"

ggsave(
  filename = fn,
  width = dims[1],
  height = dims[2],
  web_fonts = list(public_sans)
)

remove(N, sites, s, rgmm, dgmm, chrono, p1, p2, dims, fn)

# mesa verde -------------------------------------------------------------

bb8 <- tree_rings |>
  st_as_sf(coords = c("x", "y"), crs = 26912) |>
  st_buffer(20000) |>
  st_transform(4269) |>
  st_bbox()

dy <- bb8[["ymax"]] - bb8[["ymin"]]
dx <- dy * 3

mx <- mean(bb8[c("xmax", "xmin")])

bb8[["xmin"]] <- mx - (dx / 2)
bb8[["xmax"]] <- mx + (dx / 2)

vrt <- file.path(
  "https://prd-tnm.s3.amazonaws.com/StagedProducts",
  "Elevation/13/TIFF",
  "USGS_Seamless_DEM_13.vrt"
)

dem <- vrt |>
  rast(win = ext(bb8)) |>
  aggregate(fact = 4)

slope <- terrain(dem, "slope", unit = "radians")
aspect <- terrain(dem, "aspect", unit = "radians")

hill <- shade(slope, aspect, 30, 45)
hill <- setValues(hill, scales::rescale(values(hill), to = c(1, 1000)))
hill <- round(hill)

bb9 <- tree_rings |>
  st_as_sf(coords = c("x", "y"), crs = 26912) |>
  st_buffer(1000) |>
  st_transform(4269) |>
  st_bbox() |>
  st_as_sfc()

cover <- st_sym_difference(bb9, st_as_sfc(bb8))

p1 <- ggplot() +
  geom_spatraster(
    data = hill,
    fill = hcl.colors(1000, "Grays")[values(hill)],
    maxcell = Inf
  ) +
  geom_spatraster(data = dem, alpha = 0.5) +
  scale_fill_hypso_c("utah_1") +
  geom_sf(
    data = cover,
    fill = alpha("white", 0.5),
    color = "transparent"
  ) +
  geom_sf(
    data = four_corners,
    fill = "transparent",
    color = coolors[["gray"]],
    linewidth = 0.2
  ) +
  geom_sf(
    data = bb9,
    fill = "transparent",
    color = coolors[["purple"]],
    linewidth = 0.6
  ) +
  annotate(
    "text",
    x = c(rep(bb8[["xmin"]] + 0.02, 2), rep(bb8[["xmax"]] - 0.02, 2)),
    y = 36.999 + c(0.02, -0.02, 0.02, -0.02),
    label = c("Utah", "Arizona", "Colorado", "New Mexico"),
    size = 11 / .pt,
    hjust = c(0, 0, 1, 1),
    vjust = c(0, 1, 0, 1),
    color = coolors[["gray"]]
  ) +
  coord_sf(
    crs = 4269,
    datum = 4269,
    xlim = bb8[c("xmin", "xmax")],
    ylim = bb8[c("ymin", "ymax")],
    expand = FALSE
  ) +
  theme(
    axis.title = element_blank(),
    legend.position = "none",
    panel.heights = ph,
    panel.widths = 3 * ph,
    plot.margin = margin()
  )

p1

dims <- get_dims(p1, 3 * ph, ph)

fn <- "figures/mesa-verde.png"

ggsave(
  filename = fn,
  width = dims[1],
  height = dims[2],
  bg = "white"
)

image_read(fn) |>
  image_negate() |>
  image_modulate(hue = 200) |>
  image_write("figures/mesa-verde-invert.png")

remove(
  bb8, bb9, dx, dy, mx, furl,
  hill, vrt, slope, aspect, dem, cover,
  p1, dims, fn
)

# mesa verde data --------------------------------------------------------

p1 <- ggplot(ceramics, aes(x, y)) +
  geom_point(
    shape = 21,
    fill = coolors[["purple"]],
    color = coolors_dk[["purple"]],
    alpha = 0.33,
    size = 2
  ) +
  coord_cartesian(
    xlim = range(ceramics[["x"]]),
    ylim = range(ceramics[["y"]])
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Ceramics XY"
  ) +
  theme(
    axis.title = element_text(size = rel(0.8)),
    panel.heights = ph,
    panel.widths = ph,
    plot.margin = margin(r = 4),
    plot.title = element_text(color = coolors[["purple"]])
  )

p2 <- ggplot(tree_rings, aes(x, y)) +
  geom_point(
    shape = 21,
    fill = coolors[["orange"]],
    color = coolors_dk[["orange"]],
    alpha = 0.33,
    size = 2
  ) +
  coord_cartesian(
    xlim = range(ceramics[["x"]]),
    ylim = range(ceramics[["y"]])
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Tree Rings XY"
  ) +
  theme(
    axis.title = element_text(size = rel(0.8)),
    panel.heights = ph,
    panel.widths = ph,
    plot.margin = margin(l = 4, r = 4),
    plot.title = element_text(color = coolors[["orange"]])
  )

p3 <- ggplot(room_counts) +
  geom_point(
    aes(x, y),
    shape = 21,
    fill = coolors[["blue"]],
    color = coolors_dk[["blue"]],
    alpha = 0.33,
    size = 2
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Rooms XY"
  ) +
  coord_cartesian(
    xlim = range(ceramics[["x"]]),
    ylim = range(ceramics[["y"]])
  ) +
  theme(
    axis.title = element_text(size = rel(0.8)),
    legend.position = "none",
    panel.heights = ph,
    panel.widths = ph,
    plot.margin = margin(l = 4),
    plot.title = element_text(color = coolors[["blue"]])
  )

p1 + p2 + p3

dims <- get_dims(p1, ph, ph) + get_dims(p2, ph, ph) + get_dims(p3, ph, ph)
dims[2] <- get_dims(p3, ph, ph)[2]

fn <- "figures/mesa-verde-sites.svg"

ggsave(
  filename = fn,
  width = dims[1],
  height = dims[2],
  bg = "white",
  web_fonts = list(public_sans)
)

remove(p1, p2, p3, dims, fn)

# spatial aggregates -----------------------------------------------------

bluff_br <- ceramics |>
  filter(c23 == 1) |>
  select(x, y)

p1 <- ggplot() +
  geom_density2d_filled(
    data = bluff_br,
    aes(x, y),
    bins = 9,
    color = "white",
    linewidth = 0.1
  ) +
  scale_fill_brewer(palette = "YlOrBr") +
  geom_point(
    data = tree_rings,
    aes(x, y),
    color = alpha(coolors[["blue"]], 0.1)
  ) +
  coord_cartesian(
    xlim = range(ceramics[["x"]]),
    ylim = range(ceramics[["y"]]),
    expand = FALSE
  ) +
  labs(
    x = "Easting",
    y = "Northing",
    title = "Bluff Black-on-red"
  ) +
  theme(
    legend.position = "none",
    panel.heights = ph,
    panel.widths = ph,
    panel.background = element_rect(fill = "#FFFFE5"),
    plot.margin = margin(r = 6)
  )

d <- MASS::kde2d(bluff_br[["x"]], bluff_br[["y"]])
d <- rast(d)

names(d) <- "density"

p2 <- d |>
  terra::extract(as.matrix(tree_rings[c("x", "y")])) |>
  as_tibble() |>
  mutate(date = tree_rings[["date"]]) |>
  filter(!is.na(density)) |>
  ggplot(aes(density, date)) +
  geom_point(alpha = 0.1, color = coolors[["blue"]]) +
  labs(
    x = "Density of Bluff B-r Sites",
    y = "Tree Ring Date"
  ) +
  theme(
    panel.heights = ph,
    panel.widths = ph,
    plot.margin = margin(l = 6)
  )

p1 + p2

dims <- get_dims(p1, ph, ph) + get_dims(p2, ph, ph)
dims[2] <- get_dims(p1, ph, ph)[2]

fn <- "figures/bluff-br.svg"

ggsave(
  filename = fn,
  width = dims[1],
  height = dims[2],
  bg = "white",
  web_fonts = list(public_sans)
)

remove(p1, d, p2, dims, fn)

# regional chronology ----------------------------------------------------

trd <- density(tree_rings[["date"]])
trd <- trd[c("x", "y")] |>
  as_tibble() |>
  rename("year" = x, "density" = y) |>
  mutate(density = density / max(density)) |>
  filter(year >= 575, year <= 1400)

p1 <- ggplot(trd, aes(year)) +
  geom_ribbon(
    aes(ymin = 0, ymax = density),
    color = "transparent",
    fill = alpha(coolors[["orange"]], 0.5)
  ) +
  geom_line(
    aes(y = density),
    color = coolors_dk[["orange"]]
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Kernel Density"
  ) +
  scale_x_continuous(
    limits = c(575, 1400),
    breaks = c(800, 1000, 1200)
  ) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  theme(
    axis.text.x = element_text(),
    panel.heights = ph,
    panel.widths = ph,
    plot.margin = margin(r = 3)
  )

rd <- region |> mutate(density = density / max(density))

p2 <- ggplot(rd, aes(year)) +
  geom_ribbon(
    aes(ymin = 0, ymax = density),
    color = "transparent",
    fill = alpha(coolors[["blue"]], 0.5)
  ) +
  geom_line(
    aes(y = density),
    color = coolors_dk[["blue"]]
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Mixture Density"
  ) +
  scale_x_continuous(
    limits = c(575, 1400),
    breaks = c(800, 1000, 1200)
  ) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  theme(
    axis.text.x = element_text(),
    panel.heights = ph,
    panel.widths = ph,
    plot.margin = margin(l = 3)
  )

p1 + p2

dims <- get_dims(p1, ph, ph) + get_dims(p2, ph, ph)
dims[2] <- get_dims(p1, ph, ph)[2]

fn <- "figures/regional-density.svg"

ggsave(
  filename = fn,
  width = dims[1],
  height = dims[2],
  bg = "white",
  web_fonts = list(public_sans)
)

remove(trd, p1, rd, p2, dims, fn)

# site chronologies ------------------------------------------------------

sites1 <- sites |>
  filter(
    site_name %in% c(
      "Castle Rock Pueblo", "Grass Mesa Village", "Yellow Jacket Pueblo"
    )
  )

tr1 <- tree_rings |>
  filter(trinomial %in% unique(sites1[["site_trinomial"]])) |>
  select(trinomial, date) |>
  rename(site_trinomial = trinomial) |>
  left_join(
    sites1 |> select(site_trinomial, site_name) |> distinct(),
    by = "site_trinomial"
  )

p1 <- ggplot(sites1, aes(year)) +
  geom_ribbon(
    aes(ymin = 0, ymax = density),
    color = "transparent",
    fill = alpha(coolors[["blue"]], 0.5)
  ) +
  geom_line(
    aes(y = density),
    color = coolors_dk[["blue"]]
  ) +
  geom_point(
    data = tr1,
    aes(x = date, y = 0.05),
    size = 1.5,
    alpha = 0.4,
    color = coolors[["yellow"]]
  ) +
  facet_wrap(vars(site_name)) +
  labs(
    x = NULL,
    y = NULL
  ) +
  scale_x_continuous(
    limits = c(575, 1400),
    breaks = c(800, 1000, 1200)
  ) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  theme(
    axis.text.x = element_text(),
    panel.heights = rep(ph, 3),
    panel.widths = rep(ph, 3),
    plot.margin = margin(),
    strip.background = element_blank(),
    strip.text = element_text(size = rel(1.1), hjust = 0)
  )

p1

dims <- get_dims(p1, 3.1 * ph, ph)

fn <- "figures/important-sites-1.svg"

ggsave(
  filename = fn,
  width = dims[1],
  height = dims[2],
  bg = "white",
  web_fonts = list(public_sans)
)

remove(sites1, tr1, p1, dims, fn)

sites2 <- sites |>
  filter(
    !site_name %in% c(
      "Castle Rock Pueblo", "Grass Mesa Village", "Yellow Jacket Pueblo"
    )
  )

tr2 <- tree_rings |>
  filter(trinomial %in% unique(sites2[["site_trinomial"]])) |>
  select(trinomial, date) |>
  rename(site_trinomial = trinomial) |>
  left_join(
    sites2 |> select(site_trinomial, site_name) |> distinct(),
    by = "site_trinomial"
  )

p2 <- ggplot(sites2, aes(year)) +
  geom_ribbon(
    aes(ymin = 0, ymax = density),
    color = "transparent",
    fill = alpha(coolors[["blue"]], 0.5)
  ) +
  geom_line(
    aes(y = density),
    color = coolors_dk[["blue"]]
  ) +
  geom_point(
    data = tr2,
    aes(x = date, y = 0.05),
    size = 1.5,
    alpha = 0.4,
    color = coolors[["yellow"]]
  ) +
  facet_wrap(vars(site_name)) +
  labs(
    x = NULL,
    y = NULL
  ) +
  scale_x_continuous(
    limits = c(575, 1400),
    breaks = c(800, 1000, 1200)
  ) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  theme(
    axis.text.x = element_text(),
    panel.heights = rep(ph, 2),
    panel.widths = rep(ph, 2),
    plot.margin = margin(),
    strip.background = element_blank(),
    strip.text = element_text(size = rel(1.1), hjust = 0)
  )

p2

dims <- get_dims(p2, 2.1 * ph, ph)

fn <- "figures/important-sites-2.svg"

ggsave(
  filename = fn,
  width = dims[1],
  height = dims[2],
  bg = "white",
  web_fonts = list(public_sans)
)

remove(sites2, tr2, p2, dims, fn)
