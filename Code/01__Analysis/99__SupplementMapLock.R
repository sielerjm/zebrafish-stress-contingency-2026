# 99__SupplementMapLock.R
#
# Created by: Michael Sieler
# Date last updated: 2026-04-25
#
# Description:
#   Filesystem-level guardrail for the curated supplementary figure-map YAML files under
#   `Manuscript/Supplementary/`. Useful as an extra protection against accidental overwrites.
#
# Expected input:
#   - CLI: `lock` | `unlock` | `status`
#   - Optional: `--mode-lock=<octal>` / `--mode-unlock=<octal>`
#   - Example:
#       - `Rscript Code/01__Analysis/99__SupplementMapLock.R lock`
#       - `Rscript Code/01__Analysis/99__SupplementMapLock.R unlock`
#       - `Rscript Code/01__Analysis/99__SupplementMapLock.R status`
#
# Expected output:
#   - Prints which YAMLs exist and their file modes; optionally changes modes via chmod.
#

options(stringsAsFactors = FALSE)

script_path_from_cli <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", ca[startsWith(ca, "--file=")])
  if (length(f) == 1L && nzchar(f)) {
    return(normalizePath(f, winslash = "/", mustWork = TRUE))
  }
  NA_character_
}

.sp <- script_path_from_cli()
if (!is.na(.sp)) {
  proj.path <- normalizePath(file.path(dirname(.sp), "..", ".."), winslash = "/", mustWork = TRUE)
} else {
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Install here: install.packages('here')")
  }
  here::i_am("Code/01__Analysis/99__SupplementMapLock.R")
  proj.path <- as.character(here::here())
}

old_wd <- setwd(proj.path)
on.exit(setwd(old_wd), add = TRUE)

args <- commandArgs(trailingOnly = TRUE)
cmd <- args[!grepl("^--", args)][[1L]]
if (is.null(cmd) || !nzchar(cmd) || !cmd %in% c("lock", "unlock", "status")) {
  stop(
    "Usage:\n",
    "  Rscript Code/01__Analysis/99__SupplementMapLock.R lock\n",
    "  Rscript Code/01__Analysis/99__SupplementMapLock.R unlock\n",
    "  Rscript Code/01__Analysis/99__SupplementMapLock.R status\n",
    "\n",
    "Optional:\n",
    "  --mode-lock=0444\n",
    "  --mode-unlock=0644\n",
    call. = FALSE
  )
}

get_flag_value <- function(flag_name, default = NULL) {
  hit <- args[startsWith(args, paste0(flag_name, "="))]
  if (length(hit) < 1L) {
    return(default)
  }
  sub(paste0("^", flag_name, "="), "", hit[[1L]])
}

mode_lock <- get_flag_value("--mode-lock", default = "0444")
mode_unlock <- get_flag_value("--mode-unlock", default = "0644")

path.supp <- file.path(proj.path, "Manuscript", "Supplementary")

curated <- c(
  "supplement_map__01__Diversity.yml",
  "supplement_map__02__Composition.yml",
  "supplement_map__03__DiffAbund.yml",
  "supplement_map__04__DiffGeneExp.yml",
  "supplement_map__05__Mort-Inf.yml",
  "supplement_map__06__Taxon-DEG-Mort.yml",
  "supplement_map__07__FunctionalAnno.yml",
  "supplement_map__08__NeutralModel.yml"
)

paths <- file.path(path.supp, curated)
exists <- file.exists(paths)

mode_string <- function(p) {
  fi <- file.info(p)
  if (nrow(fi) != 1L || is.na(fi$mode[[1L]])) {
    return(NA_character_)
  }
  as.character(fi$mode[[1L]])
}

print_status <- function() {
  message("Supplementary map directory: ", path.supp)
  for (i in seq_along(paths)) {
    p <- paths[[i]]
    if (!file.exists(p)) {
      message("MISSING\t", curated[[i]])
      next
    }
    message(mode_string(p), "\t", curated[[i]])
  }
  invisible(TRUE)
}

if (identical(cmd, "status")) {
  print_status()
} else if (identical(cmd, "lock")) {
  changed <- 0L
  for (p in paths[exists]) {
    ok <- Sys.chmod(p, mode = mode_lock, use_umask = FALSE)
    if (isTRUE(ok)) {
      changed <- changed + 1L
    }
  }
  message("Locked ", changed, " map file(s) with mode ", mode_lock, ".")
  print_status()
} else if (identical(cmd, "unlock")) {
  changed <- 0L
  for (p in paths[exists]) {
    ok <- Sys.chmod(p, mode = mode_unlock, use_umask = FALSE)
    if (isTRUE(ok)) {
      changed <- changed + 1L
    }
  }
  message("Unlocked ", changed, " map file(s) with mode ", mode_unlock, ".")
  print_status()
}

