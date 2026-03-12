# S01: correct and harmonize the active COMPADRE snapshot for this project.

# libraries ----
source("code/setup.R")
setup_packages(c("tidyverse", "Rcompadre"))
source("code/functions.R")


# load compadre data ----
compadre <- load_compadre(corrected = FALSE)


# fix ehrlen ----
ehrlen_fix <- which(
  compadre$SpeciesAuthor == "Lathyrus_vernus" &
    compadre$MatrixPopulation == "G"
)

purrr::walk(ehrlen_fix, ~ {
  if (nrow(compadre$mat[[.x]]@matA) == 7) {
    compadre$mat[[.x]]@matA <- compadre$mat[[.x]]@matA[-7, -7]
    compadre$mat[[.x]]@matU <- compadre$mat[[.x]]@matU[-7, -7]
    compadre$mat[[.x]]@matF <- compadre$mat[[.x]]@matF[-7, -7]
    compadre$mat[[.x]]@matC <- compadre$mat[[.x]]@matC[-7, -7]
    compadre$mat[[.x]]@matrixClass <- compadre$mat[[.x]]@matrixClass[-7, ]
  }
})


# fix Aschero ----
aschero_fix <- which(
  compadre$Authors == "Aschero; Morris; Vázquez; Alvarez; Villagra" &
    compadre$MatrixTreatment == "Unmanipulated"
)

if (length(aschero_fix) >= 1) {
  load("data/derived/studies/aschero_U.RData")
  for (i in aschero_fix) {
    if (!isTRUE(all.equal(compadre$mat[[i]]@matU, U, tolerance = 0))) {
      compadre$mat[[i]]@matU <- U
      compadre$mat[[i]]@matA <- compadre$mat[[i]]@matU +
        compadre$mat[[i]]@matF + compadre$mat[[i]]@matC
    }
  }
  rm(U)
}


# fix Lazaro clonal transitions recorded in C rather than U ----
lazaro_fix <- which(
  compadre$SpeciesAuthor == "Dioon_merolae" &
    compadre$MatrixPopulation == "EC"
)

for (i in lazaro_fix) {
  if (any(compadre$mat[[i]]@matC > 0)) {
    compadre$mat[[i]]@matU <- compadre$mat[[i]]@matU + compadre$mat[[i]]@matC
    compadre$mat[[i]]@matC[compadre$mat[[i]]@matC > 0] <- 0
  }
}


# fix Portela clonal transitions recorded in C rather than U ----
portela_fix1 <- which(
  compadre$SpeciesAuthor == "Astrocaryum_aculeatissimum"
)

for (i in portela_fix1) {
  if (any(compadre$mat[[i]]@matC > 0)) {
    compadre$mat[[i]]@matU <- compadre$mat[[i]]@matU + compadre$mat[[i]]@matC
    compadre$mat[[i]]@matC[compadre$mat[[i]]@matC > 0] <- 0
  }
}


# fix Plank Trillium population labels in recent COMPADRE versions ----
trillium_fix <- which(
  compadre$SpeciesAuthor == "Trillium_persistens" &
    compadre$MatrixTreatment == "Unmanipulated"
)

if (length(trillium_fix) > 0) {
  tr_min_mislabeled <- which(
    compadre$SpeciesAuthor == "Trillium_persistens" &
      compadre$MatrixTreatment == "Unmanipulated" &
      compadre$MatrixPopulation == "Moody Creek" &
      seq_along(compadre$mat) %in% trillium_fix
  )

  if (length(tr_min_mislabeled) > 0) {
    tr_min_mislabeled <- tr_min_mislabeled[
      vapply(compadre$mat[tr_min_mislabeled], function(x) abs(x@matF[1, 4] - 0.125) < 1e-8, logical(1))
    ]
  }

  if (length(tr_min_mislabeled) > 0) {
    compadre$MatrixPopulation[tr_min_mislabeled] <- "Moody Creek; Minimum fecundity"
  }
}


# write corrected db to file ----
save(compadre, file = get_compadre_path(corrected = TRUE))
