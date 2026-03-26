# =============================================================================
# 00_load_data.R
#
# Bypassing the Truncation Problem of Truncated Levy Flights
#
# PURPOSE:
#   Load and combine all USD/CHF 1-minute CSV files into a single data frame.
#   All other scripts source this file at the top.
#
# DATA SOURCE:
#   TickData Inc. — 1-minute bid/ask quotes for USD/CHF, 2000-01 to 2015-06.
#   See data/README.md for the expected file names and column layout.
# =============================================================================

# -- Packages ------------------------------------------------------------------
library(ggplot2)
library(extrafont)
library(grid)
library(scales)

# Load Times New Roman if available (Windows: run font_import() once first).
# Falls back silently on Mac/Linux.
tryCatch(
  loadfonts(device = "win"),
  error = function(e) message("extrafont: skipping font load (non-Windows or fonts not imported).")
)

# -- Path configuration --------------------------------------------------------
# All CSV files should be in the data/ folder relative to the project root.
# If you run scripts from a different working directory, adjust DATA_DIR below.

DATA_DIR <- file.path("data")

# -- File list -----------------------------------------------------------------
# Years 2000-2014 have one file each; 2015 is split by month (01-06).

years_full  <- 2000:2014
months_2015 <- sprintf("%02d", 1:6)

files_full  <- file.path(DATA_DIR,
                 sprintf("DAT_MT_USDCHF_M1_%d.csv", years_full))

files_2015  <- file.path(DATA_DIR,
                 sprintf("DAT_MT_USDCHF_M1_2015%s.csv", months_2015))

all_files   <- c(files_full, files_2015)

# -- Loading function ----------------------------------------------------------
# Each CSV has no header. Columns: V1=date, V2=time, V3=bid, V4=ask.
# Original code reversed row order (order(-index)); preserved here.

load_one <- function(path) {
  df        <- read.csv(path, header = FALSE)
  df$index  <- seq_len(nrow(df))
  df        <- df[order(-df$index), ]
  df
}

# -- Load and combine ----------------------------------------------------------
# Uses lapply + do.call(rbind) — idiomatic R replacement for the manual rbind loop.

message("Loading ", length(all_files), " CSV files from: ", DATA_DIR)

data_list <- lapply(all_files, function(f) {
  if (!file.exists(f)) {
    warning("File not found, skipping: ", f)
    return(NULL)
  }
  load_one(f)
})

# Drop any NULLs (missing files) and combine
data_list <- Filter(Negate(is.null), data_list)
data.1    <- do.call(rbind, data_list)
rownames(data.1) <- NULL

# -- Parse date column ---------------------------------------------------------
data.1$date <- as.Date(strptime(data.1$V1, "%Y.%m.%d"))

# -- Column reference ----------------------------------------------------------
# data.1$V1   : date string  (YYYY.MM.DD)
# data.1$V2   : intraday time (integer, minutes since midnight * some factor)
# data.1$V3   : bid price
# data.1$V4   : ask price
# data.1$date : parsed Date

message("Data loaded. Rows: ", nrow(data.1),
        " | Date range: ", min(data.1$date), " to ", max(data.1$date))
