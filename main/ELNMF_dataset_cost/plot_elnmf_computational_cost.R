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

suppressPackageStartupMessages({
  library(ggplot2)
})

output_dir <- "results/ELNMF_computational_cost"
input_csv <- file.path(output_dir, "ELNMF_runtime_peak_memory.csv")

if (!file.exists(input_csv)) {
  stop("Benchmark CSV not found: ", input_csv)
}

dat <- read.csv(input_csv, check.names = FALSE, stringsAsFactors = FALSE)
dat <- dat[dat$Status == "Completed", , drop = FALSE]
if (nrow(dat) == 0L) stop("No completed benchmark runs are available.")

dataset_order <- c("MVC", "HBC", "DLPFC 151672")
dat$Dataset <- factor(dat$Dataset, levels = dataset_order)

palette <- c(
  "MVC" = "#D97757",
  "HBC" = "#2A9D8F",
  "DLPFC 151672" = "#7562A8"
)

base_theme <- theme_classic(base_size = 12, base_family = "Arial") +
  theme(
    axis.title = element_text(size = 11.5, face = "bold"),
    axis.text = element_text(size = 10.5, colour = "#1F2937"),
    axis.line = element_line(linewidth = 0.45, colour = "#374151"),
    axis.ticks = element_line(linewidth = 0.4, colour = "#374151"),
    legend.position = "none",
    plot.margin = margin(12, 16, 12, 12)
  )

runtime_max <- max(dat$Runtime_minutes, na.rm = TRUE)
memory_max <- max(dat$Peak_RAM_percent, na.rm = TRUE)

p_runtime <- ggplot(dat, aes(x = Dataset, y = Runtime_minutes, fill = Dataset)) +
  geom_col(width = 0.62, colour = "#27313D", linewidth = 0.35) +
  geom_text(
    aes(label = sprintf("%.1f min", Runtime_minutes)),
    vjust = -0.45,
    size = 3.7,
    family = "Arial",
    colour = "#1F2937"
  ) +
  scale_fill_manual(values = palette) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.17)),
    limits = c(0, runtime_max * 1.18)
  ) +
  labs(
    x = "Dataset",
    y = "Runtime (minutes)"
  ) +
  base_theme

p_memory <- ggplot(dat, aes(x = Dataset, y = Peak_RAM_percent, fill = Dataset)) +
  geom_col(width = 0.62, colour = "#27313D", linewidth = 0.35) +
  geom_text(
    aes(label = sprintf("%.1f%%\n(%.2f GB)", Peak_RAM_percent, Peak_RAM_GB)),
    vjust = -0.25,
    lineheight = 0.95,
    size = 3.6,
    family = "Arial",
    colour = "#1F2937"
  ) +
  scale_fill_manual(values = palette) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.22)),
    limits = c(0, memory_max * 1.24)
  ) +
  labs(
    x = "Dataset",
    y = "Peak memory usage (%)"
  ) +
  base_theme

ggsave(
  file.path(output_dir, "ELNMF_runtime.png"),
  p_runtime,
  width = 6.2,
  height = 5.2,
  dpi = 400,
  bg = "white"
)
ggsave(
  file.path(output_dir, "ELNMF_peak_memory.png"),
  p_memory,
  width = 6.2,
  height = 5.2,
  dpi = 400,
  bg = "white"
)

message("Saved two PNG computational-cost panels to: ", output_dir)
