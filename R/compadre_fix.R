

### libraries
library(tidyverse)
library(Rcompadre)
source("R/functions.R")



### load compadre data
compadre <- cdb_fetch("data/COMPADRE_v.X.X.X.RData")


### Load data from Ellis et al. (2012)
ellis_data <- read.table('data/ellis/Transition_Matrices.txt', sep = '\t',
                         header = TRUE, stringsAsFactors = FALSE) %>%
  as_tibble() %>% 
  mutate(matA = lapply(Mx, stringToMat)) %>% 
  mutate(matU = lapply(Tmx, stringToMat)) %>% 
  mutate(matF = mapply(function(a, b) a - b, matA, matU, SIMPLIFY = F)) %>% 
  mutate(N = lapply(Nx, NxToVec))




## fix typo in A matrix for Eriogonum longifolium (3.420 should be 0.342)
satterthwaite_fix <- which(
  compadre$SpeciesAuthor == 'Eriogonum_longifolium_var._gnaphalifolium_2' &
    compadre$MatrixPopulation == 'Unburned' &
    compadre$MatrixStartYear == 1991
)

compadre$mat[[satterthwaite_fix]]@matA[5,5] <- 0.342


## fix lazaro MatrixComposite
lazaro_fix <- which(
  compadre$SpeciesAuthor == 'Dioon_merolae' &
    (compadre$MatrixEndYear - compadre$MatrixStartYear  == 1)
)

compadre$MatrixComposite[lazaro_fix] <- "Individual"


## fix ehrlen
ehrlen_fix <- which(
  compadre$SpeciesAuthor == "Lathyrus_vernus" &
    compadre$MatrixPopulation == "G"
)

for(i in ehrlen_fix) {
  compadre$mat[[i]]@matA <- compadre$mat[[i]]@matA[-7,-7]
  compadre$mat[[i]]@matU <- compadre$mat[[i]]@matU[-7,-7]
  compadre$mat[[i]]@matF <- compadre$mat[[i]]@matF[-7,-7]
  compadre$mat[[i]]@matC <- compadre$mat[[i]]@matC[-7,-7]
  compadre$mat[[i]]@matrixClass <- compadre$mat[[i]]@matrixClass[-7,]
}


## fix lemke
lemke_fix <- which(compadre$SpeciesAuthor == "Trollius_europaeus")
compadre$MatrixTreatment[lemke_fix] <- "Unmanipulated"


## fix dostalek
dostalek_fix <- which(
  compadre$SpeciesAuthor == "Dracocephalum_austriacum_2" &
    compadre$MatrixTreatment == "Mean"
)

compadre$MatrixComposite[dostalek_fix] <- "Mean"
compadre$MatrixTreatment[dostalek_fix] <- "Unmanipulated"


## fix Aschero
aschero_fix <- which(
  compadre$Authors == "Aschero; Morris; Vázquez; Alvarez; Villagra" &
    compadre$MatrixTreatment == "Unmanipulated"
)

load("data/studies/aschero_U.RData")
compadre$mat[[aschero_fix]]@matU <- U
compadre$mat[[aschero_fix]]@matA <- compadre$mat[[aschero_fix]]@matU +
  compadre$mat[[aschero_fix]]@matF + compadre$mat[[aschero_fix]]@matC
rm(U)

## fix Astragalus_scaphoides_2
assc_fix <- which(
  compadre$SpeciesAuthor == "Astragalus_scaphoides_2" &
    compadre$MatrixPopulation == "McDevitt Creek" &
    compadre$MatrixComposite == "Mean"
)

assc_rep <- ellis_data %>% 
  filter(SPP == "ASSC", POP == "ASSC_mcdevi") %>% 
  .$matF %>% 
  mat_mean()

compadre$mat[[assc_fix]]@matF <- assc_rep
compadre$mat[[assc_fix]]@matA <- compadre$mat[[assc_fix]]@matU +
  compadre$mat[[assc_fix]]@matF + compadre$mat[[assc_fix]]@matC
rm(assc_rep)

# fix Shryock; Esque; Hughes
shryock_fix1 <- which(
  compadre$SpeciesAuthor == "Pediocactus_bradyi" &
    compadre$MatrixPopulation == "Badger Creek" &
    compadre$MatrixStartYear == 1999 &
    compadre$MatrixEndYear == 2000
)

shryock_fix2 <- which(
  compadre$SpeciesAuthor == "Pediocactus_bradyi" &
    compadre$MatrixPopulation == "Soap Creek" &
    compadre$MatrixStartYear == 1997 &
    compadre$MatrixEndYear == 1998
)

i1 <- is.na(compadre$mat[[shryock_fix1]]@matU)
i2 <- is.na(compadre$mat[[shryock_fix2]]@matU)

compadre$mat[[shryock_fix1]]@matU[i1] <- 0
compadre$mat[[shryock_fix2]]@matU[i2] <- 0


## lazaro
lazaro_fix <- which(
  compadre$SpeciesAuthor == "Dioon_merolae" &
    compadre$MatrixPopulation == "EC"
)

for (i in lazaro_fix) {
  compadre$mat[[i]]@matU <- compadre$mat[[i]]@matU + compadre$mat[[i]]@matC
  compadre$mat[[i]]@matC[compadre$mat[[i]]@matC > 0] <- 0
}


## portela
portela_fix1 <- which(
  compadre$SpeciesAuthor == "Astrocaryum_aculeatissimum"
)

for (i in portela_fix1) {
  compadre$mat[[i]]@matU <- compadre$mat[[i]]@matU + compadre$mat[[i]]@matC
  compadre$mat[[i]]@matC[compadre$mat[[i]]@matC > 0] <- 0
}



### write corrected db to file
save(compadre, file = "data/COMPADRE_v.X.X.X_Corrected.RData")

