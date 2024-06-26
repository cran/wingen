% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resist_gd.R
\name{resist_gd}
\alias{resist_gd}
\title{Create a moving window map of genetic diversity based on resistance}
\usage{
resist_gd(
  gen,
  coords,
  lyr,
  maxdist,
  distmat = NULL,
  stat = "pi",
  fact = 0,
  rarify = FALSE,
  rarify_n = 2,
  rarify_nit = 5,
  min_n = 2,
  fun = mean,
  L = "nvariants",
  rarify_alleles = TRUE,
  sig = 0.05,
  transitionFunction = mean,
  directions = 8,
  geoCorrection = TRUE
)
}
\arguments{
\item{gen}{genetic data either as an object of type vcf or a path to a vcf file (\emph{note:} order matters! The coordinate and genetic data should be in the same order; there are currently no checks for this)}

\item{coords}{coordinates of samples as sf points, a two-column matrix, or a data.frame representing x and y coordinates (see Details for important information about projections)}

\item{lyr}{conductivity layer (higher values should mean greater conductivity) to move window across. Can be either a SpatRaster or RasterLayer.}

\item{maxdist}{maximum cost distance used to define neighborhood; any samples further than this cost distance will not be included (this can be thought of as the neighborhood radius, but in terms of cost distance).
Can either be (1) a single numeric value or (2) a SpatRaster where each pixel is the maximum distance to be used for that cell on the landscape (must be the same spatial scale as \code{lyr}).}

\item{distmat}{distance matrix output from \link[wingen]{get_resdist} (optional; can be used to save time on distance calculations)}

\item{stat}{genetic diversity statistic(s) to calculate (see Details, defaults to \code{"pi"}). Can be a single statistic or a vector of statistics}

\item{fact}{aggregation factor to apply to \code{lyr} (defaults to 0; \emph{note:} increasing this value reduces computational time)}

\item{rarify}{if rarify = TRUE, rarefaction is performed (defaults to FALSE)}

\item{rarify_n}{if rarify = TRUE, number of points to use for rarefaction (defaults to min_n)}

\item{rarify_nit}{if rarify = TRUE, number of iterations to use for rarefaction (defaults to 5). Can also be set to \code{"all"} to use all possible combinations of samples of size \code{rarify_n} within the window.}

\item{min_n}{minimum number of samples to use in calculations (any focal cell with a window containing less than this number of samples will be assigned a value of NA; defaults to 2)}

\item{fun}{function to use to summarize rarefaction results (defaults to mean, must take \code{na.rm = TRUE} as an argument)}

\item{L}{for calculating \code{"pi"}, L argument in \link[hierfstat]{pi.dosage} function. Return the average nucleotide diversity per nucleotide given the length L of the sequence. The wingen default is L = "nvariants", which sets L to the number of variants in the VCF. If L = NULL, returns the sum over SNPs of nucleotide diversity (\emph{note:} L = NULL is the \link[hierfstat]{pi.dosage} default which wingen does not use)}

\item{rarify_alleles}{for calculating \code{"biallelic_richness"}, whether to perform rarefaction of allele counts as in \link[hierfstat]{allelic.richness} (defaults to TRUE)}

\item{sig}{for calculating \code{"hwe"}, significance threshold (i.e., alpha level) to use for hardy-weinberg equilibrium tests (defaults to 0.05)}

\item{transitionFunction}{function to calculate transition values from grid values (defaults to mean)}

\item{directions}{directions in which cells are connected (4, 8, 16, or other), see \link[raster]{adjacent} (defaults to 8)}

\item{geoCorrection}{whether to apply correction to account for local distances (defaults to TRUE). Geographic correction is necessary for all objects of the class Transition that are either: (1) based on a grid in a geographic (lonlat) projection and covering a large area; (2) made with directions > 4 (see \link[gdistance]{geoCorrection} for more details).}
}
\value{
SpatRaster that includes a raster layer of genetic diversity and a raster layer of the number of samples within the window for each cell
}
\description{
Generate a continuous raster map of genetic diversity using resistance distances calculated with a conductivity surface
}
\details{
Coordinates and rasters should be in a Euclidean coordinate system (i.e., UTM coordinates) such that raster cell width and height are equal distances.
As such, longitude-latitude systems should be transformed before using dist_gd. Transformation can be performed using \link[sf]{st_set_crs} for coordinates or \link[terra]{project} for rasters (see vignette for more details).

Coordinates and rasters should be in a projected (planar) coordinate system such that raster cells are of equal sizes.
Therefore, spherical systems (including latitute-longitude coordinate systems) should be projected prior to use.
Transformation can be performed using \link[sf]{st_set_crs} for coordinates or \link[terra]{project} for rasters (see vignette for more details).

Current genetic diversity metrics that can be specified with \code{stat} include:
\itemize{
\item \code{"pi"} for nucleotide diversity (default) calculated using \code{hierfstat} \link[hierfstat]{pi.dosage}. Use the \code{L} argument to set the sequence length (defaults to dividing by the number of variants).
\item \code{"Ho"} for average observed heterozygosity across all sites
\item \code{"allelic_richness"} for average number of alleles across all sites
\item \code{"biallelic_richness"} for average allelic richness across all sites for a biallelic dataset (this option is faster than \code{"allelic_richness"})
\item \code{"hwe"} for the proportion of sites that are not in Hardy–Weinberg equilibrium, calculated using \code{pegas} \link[pegas]{hw.test} at the 0.05 level (other alpha levels  can be specified by adding the sig argument; e.g., \code{sig = 0.10}).
\item \code{"basic_stats"} for a series of statistics produced by \code{hierfstat} \link[hierfstat]{basic.stats} including
mean observed heterozygosity (same as Ho), mean gene diversities within population (Hs),
Gene diversities overall (Ht), and Fis following Nei (1987). Population-based statistics (e.g., FST) normally reported by \link[hierfstat]{basic.stats}
are not included as they are not meaningful within the individual-based moving windows.
}
}
\examples{
load_mini_ex()
rpi <- resist_gd(mini_vcf, mini_coords, mini_lyr, maxdist = 50)
}
