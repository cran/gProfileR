\name{gorth}
\alias{gorth}
\title{Find orthologs.}
\usage{
  gorth(query, source_organism = "hsapiens",
    target_organism = "mmusculus", region_query = F,
    numeric_ns = "", mthreshold = Inf, filter_na = T,
    df = T)
}
\arguments{
  \item{query}{list of gene IDs to be translated.}

  \item{source_organism}{name of the source organism.}

  \item{target_organism}{name of the target organism.}

  \item{region_query}{interpret query as chromosomal
  ranges.}

  \item{numeric_ns}{namespace to use for fully numeric
  IDs.}

  \item{mthreshold}{maximum number of ortholog names per
  gene to show.}

  \item{filter_na}{logical indicating whether to filter out
  results without a corresponding target name.}

  \item{df}{logical indicating whether the output will be a
  data.frame or list.}
}
\value{
  The output can be either a list or a data.frame. The list
  has an entry for every input gene. The data frame is a
  table closely corresponding to the web interface output.
}
\description{
  Interface to the g:Orth tool. Organism names are
  constructed by concatenating the first letter of the name
  and the family name. Example: human - 'hsapiens', mouse -
  'mmusculus'.
}
\details{
  To alleviate the problem of having many orthologs per
  gene (most of them uninformative) one can set a threshold
  for the number of results. The program tries to find the
  most informative by selecting the most popular ones.
}
\examples{
gorth(c("Klf4","Pax5","Sox2","Nanog"), source_organism="mmusculus", target_organism="hsapiens")
}
\author{
  Raivo Kolde <rkolde@gmail.com>, Juri Reimand
  <juri.reimand@ut.ee>, Tambet Arak <tambet.arak@gmail.com>
}
\references{
  J. Reimand, M. Kull, H. Peterson, J.  Hansen, J. Vilo:
  g:Profiler -- a web-based toolset for functional
  profiling of gene lists from large-scale experiments
  (2007) NAR 35 W193-W200
}

