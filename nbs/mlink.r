library("paramlink2")
library("pedprobr")
mlink <- function(x, aff,freq=NULL, model, rho=seq(0,0.45,0.05), liability = NULL, markers = NULL, maxOnly = NA,
               loopBreakers = NULL, peelOrder = NULL, verbose = FALSE) {

  if(model$chrom == "X")
    stop2("X-linked disease models are not implemented at the moment.")

  if (is.singleton(x))
    stop2("This function is not applicable to singleton objects.")

  if(hasInbredFounders(x))
    stop2("Likelihood of linked markers is not implemented in pedigrees with founder inbreeding.\n",
          "(Note that this is usually not well-defined)")

  if(hasSelfing(x))
    stop2("Likelihood of pedigrees with selfing is not implemented.\n",
          "Contact the maintainer if this is important to you.")

  if(!hasMarkers(x))
    stop2("The pedigree does not have any attached markers")

  if(!is.null(markers))
    x = selectMarkers(x, markers)
  else
    markers = seq_len(nMarkers(x))

  if(is.na(maxOnly))
    maxOnly = nMarkers(x) == 1

  if(hasUnbrokenLoops(x)) {
    x = breakLoops(x, loopBreakers = loopBreakers, verbose = verbose)
    markers = x$MARKERS
  }

  # Convert aff into vector of 0-1-2 status
  aff = paramlink2:::fixAff(x, aff, verbose = verbose)

  # Liability classes
  liability = if (!is.null(liability)) liability else rep(1, pedsize(x))

  peel0 = function(dat, nuc) paramlink2:::.peel_MM_AUT(dat, nuc, rho = 0.5)
  peelOrder = peelingOrder(x)
  peelProcess = pedprobr:::peelingProcess

  marker_lod = function(m,mfreq){
    if(!is.null(mfreq))
      model$dfreq = mfreq
    startDat = paramlink2:::startdata_MD_AUT(x, m, aff, model = model, liability = liability)
    # Denominator: Unlinked
    denom = peelProcess(x, m, startdata = startDat, peeler = peel0, peelOrder = peelOrder)
    mlod = vapply(rho, function(r) {
    peel_MD = function(dat, nuc) paramlink2:::.peel_MM_AUT(dat, nuc, rho = r)
    # Numerator
    numer = peelProcess(x, m, startdata = startDat, peeler = peel_MD, peelOrder = peelOrder)
    log10(numer) - log10(denom)
    },
    FUN.VALUE = 0)
    mlod
    }
  lods = mapply(marker_lod,x$MARKERS,freq)
  lods
}
