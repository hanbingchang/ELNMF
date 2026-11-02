# Run with Rscript or source(); locate the repository independently of the working directory.
repo_root <- local({
  candidates <- c(Sys.getenv("ELNMF_ROOT", unset = ""), getwd())
  for (frame in sys.frames()) {
    if (!is.null(frame$ofile)) candidates <- c(dirname(frame$ofile), candidates)
  }
  script_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(script_arg)) candidates <- c(dirname(sub("^--file=", "", script_arg[[1]])), candidates)
  found <- NULL
  for (candidate in candidates[nzchar(candidates)]) {
    current <- normalizePath(candidate, winslash = "/", mustWork = FALSE)
    repeat {
      if (dir.exists(file.path(current, "main", "shared")) && dir.exists(file.path(current, "data"))) {
        found <- current; break
      }
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
    if (!is.null(found)) break
  }
  if (is.null(found)) stop("Cannot locate ELNMF repository. Run from its root or set ELNMF_ROOT.")
  found
})
setwd(repo_root)
# Outputs are created only when this script is run, outside main/.

args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(args) >= 1) args[[1]] else "results/HBC_ELNMF_scaling"

suppressPackageStartupMessages(library(ggplot2))

benchmark_path <- file.path(output_root, "HBC_ELNMF_scaling_benchmark.csv")
failure_path <- file.path(output_root, "HBC_ELNMF_memory_test.csv")

dat <- read.csv(benchmark_path, check.names = FALSE)
failed <- dat[dat$Status != "Completed", , drop = FALSE]
if (file.exists(failure_path)) {
  historical_failure <- read.csv(failure_path, check.names = FALSE)
  common <- intersect(names(failed), names(historical_failure))
  failed <- rbind(failed[!failed$N_spots %in% historical_failure$N_spots, common, drop = FALSE],
                  historical_failure[, common, drop = FALSE])
}

dat$N_spots <- as.integer(dat$N_spots)
dat$Runtime_min <- as.numeric(dat$Runtime_min)
dat$Peak_RAM_GB <- as.numeric(dat$Peak_RAM_GB)
completed <- dat[dat$Status == "Completed", , drop = FALSE]
completed <- completed[order(completed$N_spots), ]

failed$N_spots <- as.integer(failed$N_spots)
failed$Runtime_min <- as.numeric(failed$Runtime_min)
failed$Peak_RAM_GB <- as.numeric(failed$Peak_RAM_GB)
failed$Peak_RAM_percent <- as.numeric(failed$Peak_RAM_percent)

base_theme <- theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(size = 14, colour = "#202124"),
    axis.text = element_text(size = 11.5, colour = "#3C4043"),
    axis.text.x = element_text(angle = 25, hjust = 1),
    axis.line = element_line(linewidth = 0.65, colour = "#3C4043"),
    plot.margin = margin(15, 30, 15, 15)
  )

x_breaks <- sort(unique(c(completed$N_spots, failed$N_spots)))

p_memory <- ggplot(completed, aes(x = N_spots, y = Peak_RAM_GB)) +
  geom_line(linewidth = 1.15, colour = "#2A9D8F") +
  geom_point(
    size = 3.2,
    shape = 21,
    stroke = 0.85,
    fill = "#A8DADC",
    colour = "#176B60"
  ) +
  geom_point(
    data = failed,
    aes(x = N_spots, y = Peak_RAM_GB),
    inherit.aes = FALSE,
    shape = 4,
    size = 5.2,
    stroke = 1.6,
    colour = "#C44E52"
  ) +
  geom_text(
    data = failed,
    aes(x = N_spots, y = Peak_RAM_GB,
        label = sprintf("Not completed\nobserved peak %.1f GiB", Peak_RAM_GB)),
    inherit.aes = FALSE, hjust = 1, vjust = -0.7,
    size = 4, colour = "#C44E52"
  ) +
  scale_x_continuous(
    breaks = x_breaks,
    labels = scales::comma,
    expand = expansion(mult = c(0.04, 0.13))
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.25))
  ) +
  labs(
    x = "Number of spatial locations",
    y = "Peak process RAM (GiB)"
  ) +
  base_theme

p_runtime <- ggplot(completed, aes(x = N_spots, y = Runtime_min)) +
  geom_line(linewidth = 1.15, colour = "#2878B5") +
  geom_point(
    size = 3.2,
    shape = 21,
    stroke = 0.85,
    fill = "#9AC9DB",
    colour = "#165A84"
  ) +
  geom_point(
    data = failed,
    aes(x = N_spots, y = Runtime_min),
    inherit.aes = FALSE,
    shape = 4,
    size = 5.2,
    stroke = 1.6,
    colour = "#C44E52"
  ) +
  geom_text(
    data = failed,
    aes(x = N_spots, y = Runtime_min,
        label = sprintf("Not completed\nelapsed %.1f min", Runtime_min)),
    inherit.aes = FALSE, hjust = 1, vjust = -0.7,
    size = 4, colour = "#C44E52"
  ) +
  scale_x_continuous(
    breaks = x_breaks,
    labels = scales::comma,
    expand = expansion(mult = c(0.04, 0.13))
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.25))
  ) +
  labs(
    x = "Number of spatial locations",
    y = "Elapsed runtime (min)"
  ) +
  base_theme

ggsave(
  file.path(output_root, "HBC_ELNMF_peak_memory_scaling.png"),
  p_memory,
  width = 8.2,
  height = 5.6,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(output_root, "HBC_ELNMF_runtime_scaling.png"),
  p_runtime,
  width = 8.2,
  height = 5.6,
  dpi = 600,
  bg = "white"
)
