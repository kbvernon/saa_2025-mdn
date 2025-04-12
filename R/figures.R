library(colorspace)
library(magick)
library(patchwork)
library(sf)
library(terra)
library(tidyterra)
library(tidyverse)
library(tigris)

# default theme ----------------------------------------------------------

theme_set(theme_bw(14))

theme_update(
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.grid.major = element_line(linewidth = 0.3, color = "#E9ECEF"),
  panel.grid.minor = element_blank(),
  text = element_text(family = "Public Sans")
)

# data -------------------------------------------------------------------

tree_rings <- read_csv("_data/tree-rings.csv")

ceramics <- read_csv("_data/ceramic-presence.csv")
lookup <- read_csv("_data/ceramic-lookup.csv")

four_corners <- states() |>
  rename_with(tolower) |>
  filter(name %in% c("Utah", "Colorado", "New Mexico", "Arizona")) |>
  select(name)

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
    fill = lighten("#533747", 0.5),
    color = "#533747"
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
  theme(plot.margin = margin(r = 6))

# sample gmm -------------------------------------------------------------

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

  rm(i, K, w, theta)
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
      "Site Chronology" = adjust_transparency(lighten("#533747", 0.5), 0.35),
      "Regional Chronology" = "#CD5334"
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
    plot.margin = margin(l = 6)
  )

p1 + p2 + plot_layout(widths = c(1, 2), heights = 1)

ggsave(
  "figures/goal.svg",
  width = 7.5,
  height = 2.8
)

# mesa verde -------------------------------------------------------------

bb8 <- tree_rings |>
  st_as_sf(coords = c("x", "y"), crs = 26912) |>
  st_buffer(20000) |>
  st_transform(4269) |>
  st_bbox()

dx <- bb8[["xmax"]] - bb8[["xmin"]]
dy <- bb8[["ymax"]] - bb8[["ymin"]]

dem <- file.path(
  "https://prd-tnm.s3.amazonaws.com/StagedProducts",
  "Elevation/13/TIFF",
  "USGS_Seamless_DEM_13.vrt"
) |>
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

ggplot() +
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
    color = "#829399",
    linewidth = 0.2
  ) +
  annotate(
    "text",
    x = c(rep(bb8[["xmin"]] + 0.02, 2), rep(bb8[["xmax"]] - 0.02, 2)),
    y = 36.999 + c(0.02, -0.02, 0.02, -0.02),
    label = c("UT", "AZ", "CO", "NM"),
    size = 11 / .pt,
    hjust = c(0, 0, 1, 1),
    vjust = c(0, 1, 0, 1),
    color = "#829399"
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
    legend.position = "none"
  )

ggsave(
  "figures/mesa-verde.png",
  width = 6,
  height = 6 * (dy / dx),
  bg = "white"
)

image_read("figures/mesa-verde.png") |>
  image_trim() |>
  image_write("figures/mesa-verde.png")

image_read("figures/mesa-verde.png") |>
  image_negate() |>
  image_modulate(hue = 200) |>
  image_write("figures/mesa-verde-invert.png")

# mesa verde data --------------------------------------------------------

p1 <- ggplot(ceramics, aes(x, y)) +
  geom_point(
    shape = 21,
    fill = "#533747",
    color = darken("#533747", 0.6),
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
    plot.margin = margin(r = 2),
    plot.title = element_text(color = "#533747")
  )

p2 <- ggplot(tree_rings, aes(x, y)) +
  geom_point(
    shape = 21,
    fill = "#CD5334",
    color = darken("#CD5334", 0.6),
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
    plot.margin = margin(l = 2, r = 2),
    plot.title = element_text(color = "#CD5334")
  )

p3 <- ggplot(tree_rings) +
  geom_density(
    aes(date),
    color = darken("#CD5334", 0.6),
    fill = alpha("#CD5334", 0.5)
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Cutting Date Density"
  ) +
  scale_x_continuous(breaks = c(800, 1000, 1200)) +
  theme(
    axis.text.x = element_text(),
    legend.direction = "vertical",
    legend.position = "bottom",
    plot.margin = margin(l = 2)
  )

p1 + p2 + p3

ggsave(
  filename = "figures/cmv-sites.svg",
  width = 9,
  height = 3.5
)

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
    color = alpha("#28666E", 0.1)
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
    panel.background = element_rect(fill = "#FFFFE5")
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
  geom_point(alpha = 0.1, color = "#28666E") +
  labs(
    x = "Density of Bluff B-r Sites",
    y = "Tree Ring Date"
  )

p1 + p2

ggsave(
  "figures/bluff-br.svg",
  width = 6.8,
  height = 3.5
)

# models -----------------------------------------------------------------



tree_rings <- density(tree_rings$date)
tree_rings <- tibble(year = tree_rings$x, density = tree_rings$y)
tree_rings <- tree_rings |>
  mutate(density = density / max(density)) |>
  filter(year >= 575, year <= 1400)

p1 <- ggplot(tree_rings, aes(year)) +
  geom_ribbon(
    aes(ymin = 0, ymax = density),
    color = "transparent",
    fill = alpha("#CD5334", 0.5)
  ) +
  geom_line(
    aes(y = density),
    color = darken("#CD5334", 0.6)
  ) +
  labs(
    x = NULL,
    y = "Probability Density",
    title = "Kernel Density"
  ) +
  scale_x_continuous(
    limits = c(575, 1400),
    breaks = c(600, 800, 1000, 1200, 1400)
  ) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  theme(
    axis.text.x = element_text(),
    plot.margin = margin(r = 2)
  )

region <- region |> mutate(density = density / max(density))

p2 <- ggplot(region, aes(year)) +
  geom_ribbon(
    aes(ymin = 0, ymax = density),
    color = "transparent",
    fill = alpha("#28666E", 0.5)
  ) +
  geom_line(
    aes(y = density),
    color = darken("#28666E", 0.6)
  ) +
  labs(
    x = NULL,
    y = "Probability Density",
    title = "Mixture Density"
  ) +
  scale_x_continuous(
    limits = c(575, 1400),
    breaks = c(600, 800, 1000, 1200, 1400)
  ) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  theme(
    axis.text.x = element_text(),
    plot.margin = margin(l = 2)
  )

p3 + p4 + plot_layout(axis_titles = "collect")

ggsave(
  "figures/regional-density.svg",
  width = 7,
  height = 3.5
)

tree_rings <- read_csv("_data/tree-rings.csv")


sites1 <- sites |>
  filter(
    site_name %in% c(
      "Castle Rock Pueblo", "Grass Mesa Village", "Yellow Jacket Pueblo"
    )
  )

tr1 <- tree_rings |>
  filter(trinomial %in% unique(sites1$site_trinomial)) |>
  select(trinomial, date) |>
  rename(site_trinomial = trinomial) |>
  left_join(
    sites1 |> select(site_trinomial, site_name) |> distinct(),
    by = "site_trinomial"
  )


ggplot(sites1, aes(year)) +
  geom_ribbon(
    aes(ymin = 0, ymax = density),
    color = "transparent",
    fill = alpha("#28666E", 0.5)
  ) +
  geom_line(
    aes(y = density),
    color = darken("#28666E", 0.6)
  ) +
  geom_point(
    data = tr1,
    aes(x = date, y = 0.05),
    size = 1.5,
    alpha = 0.4,
    color = "#FFD23F"
  ) +
  facet_wrap(vars(site_name)) +
  labs(
    x = NULL,
    y = "Probability Density"
  ) +
  scale_x_continuous(
    limits = c(575, 1400),
    breaks = c(600, 800, 1000, 1200, 1400),
    labels = c("", 800, 1000, 1200, "")
  ) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  theme(
    axis.text.x = element_text(),
    plot.margin = margin(l = 2),
    strip.background = element_blank(),
    strip.text = element_text(size = rel(1.1), hjust = 0)
  )

ggsave(
  filename = "figures/important-sites-1.svg",
  width = 9,
  height = 3.5
)

sites2 <- sites |>
  filter(
    !site_name %in% c(
      "Castle Rock Pueblo", "Grass Mesa Village", "Yellow Jacket Pueblo"
    )
  )

tr2 <- tree_rings |>
  filter(trinomial %in% unique(sites2$site_trinomial)) |>
  select(trinomial, date) |>
  rename(site_trinomial = trinomial) |>
  left_join(
    sites2 |> select(site_trinomial, site_name) |> distinct(),
    by = "site_trinomial"
  )

ggplot(sites2, aes(year)) +
  geom_ribbon(
    aes(ymin = 0, ymax = density),
    color = "transparent",
    fill = alpha("#28666E", 0.5)
  ) +
  geom_line(
    aes(y = density),
    color = darken("#28666E", 0.6)
  ) +
  geom_point(
    data = tr2,
    aes(x = date, y = 0.05),
    size = 1.5,
    alpha = 0.4,
    color = "#FFD23F"
  ) +
  facet_wrap(vars(site_name)) +
  labs(
    x = NULL,
    y = "Probability Density"
  ) +
  scale_x_continuous(
    limits = c(575, 1400),
    breaks = c(600, 800, 1000, 1200, 1400),
    labels = c("", 800, 1000, 1200, "")
  ) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  theme(
    axis.text.x = element_text(),
    plot.margin = margin(l = 2),
    strip.background = element_blank(),
    strip.text = element_text(size = rel(1.1), hjust = 0)
  )

ggsave(
  filename = "figures/important-sites-2.svg",
  width = 6.3,
  height = 3.5
)
