# Worshop: Pattern Analysis for Evaluating Soil Maps
#  Pedometrics 2024, Las Cruces NM (USA)
# Install required packages

# function
install.pkg <- function(pkg) {
  if (!pkg %in% rownames(installed.packages()))
  {
    install.packages(pkg, dependencies = TRUE)
  } else
  { 
    print(paste(pkg, "already installed"))
  }
}


# list of packages
workshop.pkgs <- c("terra", "raster", "sf", "sp", "gstat",
                   "motif", "philentropy", "sabre",
                   "landscapemetrics", "landscapetools",
                   "ggplot2", "gridExtra", "tidyterra",
                   "dplyr", "RColorBrewer",
                   "soilDB",
                   "glcm", "GLCMTextures",
                   "supercells", "diffeR")

# install them
for (pkg in workshop.pkgs) {
  install.pkg(pkg)
}

# extra -- not in CRAN
if (!"supercells" %in% rownames(installed.packages())) {
  install.packages("supercells", repos = "https://nowosad.r-universe.dev")
} else { print("supercells already installed") }


