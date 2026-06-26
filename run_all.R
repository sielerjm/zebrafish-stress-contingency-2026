# run_all.R
# Created by: Michael Sieler
# Date last updated: 2026-06-26
#
# Description: Run analysis drivers 00–08 in order from the repository root.
#   Assumes preprocessing has produced ps.list (see 04__DataPreProcess.R).
#
# Expected input:  Repo root as working directory; dependencies installed.
# Expected output:  Results/*/Figures, Tables, Stats for modules 00–08.

drivers <- c(
  "Code/01__Analysis/00__Overview.R",
  "Code/01__Analysis/01__Diversity.R",
  "Code/01__Analysis/02__Composition.R",
  "Code/01__Analysis/03__DiffAbund.R",
  "Code/01__Analysis/04__DiffGeneExp.R",
  "Code/01__Analysis/05__Mort-Inf.R",
  "Code/01__Analysis/06__Taxon-DEG-Mort.R",
  "Code/01__Analysis/07__FunctionalAnno.R",
  "Code/01__Analysis/08__NeutralModel.R"
)

if (!file.exists("Code/00__Setup/00__InitializeEnvironment.R")) {
  stop(
    "Set working directory to the repository root (contains Code/00__Setup/).\n",
    "Open Sieler2026.Rproj or setwd() to the clone root."
  )
}

message("Sieler2026 — running ", length(drivers), " analysis modules...")
start_time <- Sys.time()

for (script in drivers) {
  message("\n=== ", script, " ===")
  t0 <- Sys.time()
  ok <- tryCatch(
    {
      source(script, local = new.env(parent = globalenv()))
      TRUE
    },
    error = function(e) {
      message("ERROR in ", script, ": ", conditionMessage(e))
      FALSE
    }
  )
  elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 2)
  if (!ok) {
    stop("Pipeline stopped at ", script, " (", elapsed, " min elapsed). Fix the error and re-run from this module.")
  }
  message("Completed in ", elapsed, " min")
}

total <- round(difftime(Sys.time(), start_time, units = "mins"), 2)
message("\nAll modules finished in ", total, " min.")
message("Knit Code/02__Results/*.Rmd for HTML reports; run 98__MainFiguresRefresh.R for manuscript panels.")
