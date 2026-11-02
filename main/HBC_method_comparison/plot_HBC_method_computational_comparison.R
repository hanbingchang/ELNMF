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

result_root <- "results"
baseline_csv <- file.path(
  result_root,
  "HBC_baseline_computational_benchmark",
  "HBC_computational_benchmark.csv"
)
elnmf_csv <- file.path(
  result_root,
  "ELNMF_computational_cost",
  "ELNMF_runtime_peak_memory.csv"
)

baseline <- read.csv(baseline_csv, stringsAsFactors = FALSE, check.names = FALSE)
elnmf_all <- read.csv(elnmf_csv, stringsAsFactors = FALSE, check.names = FALSE)
elnmf <- elnmf_all[elnmf_all$Dataset == "HBC" & elnmf_all$Status == "Completed", , drop = FALSE]
if (nrow(elnmf) != 1L) stop("Expected exactly one completed ELNMF HBC benchmark row.")

total_ram_gb <- elnmf$Total_physical_RAM_GB[[1]]
baseline$Runtime_min <- baseline$Runtime_sec / 60
baseline$Peak_RAM_percent <- 100 * baseline$Peak_RAM_GB / total_ram_gb

plot_data <- rbind(
  data.frame(
    Method = "ELNMF",
    Runtime_sec = elnmf$Runtime_seconds,
    Runtime_min = elnmf$Runtime_minutes,
    Peak_RAM_GB = elnmf$Peak_RAM_GB,
    Peak_RAM_percent = elnmf$Peak_RAM_percent,
    Device = "CPU",
    stringsAsFactors = FALSE
  ),
  baseline[, c(
    "Method", "Runtime_sec", "Runtime_min", "Peak_RAM_GB",
    "Peak_RAM_percent", "Device"
  )]
)

method_order <- c("ELNMF", "CGNMF", "stHGC", "MNMST", "Seurat")
plot_data$Method <- factor(plot_data$Method, levels = method_order)
plot_data <- plot_data[order(plot_data$Method), , drop = FALSE]

# Color-blind-friendly palette; the proposed method is visually emphasized.
method_colors <- c(
  "ELNMF" = "#0072B2",
  "CGNMF" = "#D55E00",
  "stHGC" = "#009E73",
  "MNMST" = "#CC79A7",
  "Seurat" = "#E69F00"
)

base_theme <- theme_classic(base_size = 12, base_family = "sans") +
  theme(
    axis.title = element_text(size = 12, face = "bold", colour = "#1F2937"),
    axis.text = element_text(size = 10.5, colour = "#1F2937"),
    axis.line = element_line(linewidth = 0.5, colour = "#374151"),
    axis.ticks = element_line(linewidth = 0.45, colour = "#374151"),
    legend.position = "none",
    plot.margin = margin(14, 16, 10, 14)
  )

runtime_max <- max(plot_data$Runtime_min)
memory_max <- max(plot_data$Peak_RAM_percent)

p_runtime <- ggplot(plot_data, aes(Method, Runtime_min, fill = Method)) +
  geom_col(width = 0.68, colour = "#27313D", linewidth = 0.4) +
  geom_text(
    aes(label = sprintf("%.1f", Runtime_min)),
    vjust = -0.45,
    size = 3.5,
    colour = "#1F2937"
  ) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(
    limits = c(0, runtime_max * 1.16),
    expand = expansion(mult = c(0, 0.04))
  ) +
  labs(x = "Method", y = "Runtime (minutes)") +
  base_theme

p_memory <- ggplot(plot_data, aes(Method, Peak_RAM_percent, fill = Method)) +
  geom_col(width = 0.68, colour = "#27313D", linewidth = 0.4) +
  geom_text(
    aes(label = sprintf("%.1f%%\n(%.2f GB)", Peak_RAM_percent, Peak_RAM_GB)),
    vjust = -0.25,
    lineheight = 0.95,
    size = 3.3,
    colour = "#1F2937"
  ) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(0, memory_max * 1.19),
    expand = expansion(mult = c(0, 0.04))
  ) +
  labs(x = "Method", y = "Peak memory usage (%)") +
  base_theme

write.csv(
  transform(
    plot_data,
    Method = as.character(Method),
    Runtime_sec = round(Runtime_sec, 3),
    Runtime_min = round(Runtime_min, 3),
    Peak_RAM_GB = round(Peak_RAM_GB, 3),
    Peak_RAM_percent = round(Peak_RAM_percent, 2)
  ),
  file.path(result_root, "HBC_ELNMF_baseline_computational_comparison.csv"),
  row.names = FALSE
)

ggsave(
  file.path(result_root, "HBC_method_runtime_comparison.png"),
  p_runtime,
  width = 6.4,
  height = 5.4,
  dpi = 500,
  bg = "white"
)

ggsave(
  file.path(result_root, "HBC_method_peak_memory_comparison.png"),
  p_memory,
  width = 6.4,
  height = 5.4,
  dpi = 500,
  bg = "white"
)

message("Saved separate HBC runtime and peak-memory PNG figures.")
