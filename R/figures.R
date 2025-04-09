library(colorspace)
library(patchwork)
library(tidyverse)

theme_set(theme_bw(14))

theme_update(
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.grid.major = element_line(linewidth = 0.3, color = "#E9ECEF"),
  panel.grid.minor = element_blank(),
  text = element_text(family = "Public Sans")
)

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

tree_rings <- read_csv("_data/tree-rings.csv")

ceramics <- read_csv("_data/ceramic-presence.csv")

p1 <- ggplot(ceramics, aes(x, y)) +
  geom_point(
    shape = 21,
    fill = "#533747",
    color = darken("#533747", 0.6),
    alpha = 0.33,
    size = 2
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

region <- read_csv("_data/region.csv")

tree_rings <- density(tree_rings$date)
tree_rings <- tibble(year = tree_rings$x, density = tree_rings$y)
tree_rings <- tree_rings |>
  mutate(density = density / max(density)) |>
  filter(year >= 575, year <= 1400)

p3 <- ggplot(tree_rings, aes(year)) +
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

p4 <- ggplot(region, aes(year)) +
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
sites <- read_csv("_data/important-sites.csv")

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
