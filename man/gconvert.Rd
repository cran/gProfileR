\name{gconvert}
\alias{gconvert}
\title{Convert gene IDs.}
\usage{gconvert(query, organism="hsapiens", target="ENSG", region_query=F,
    numeric_ns="", df=T)
}
\description{Convert gene IDs.}
\details{Interface to the g:Convert tool. Organism names are constructed by concatenating the first letter of 
the name and the family name. Example: human - 'hsapiens', mouse - 'mmusculus'.}
\value{The output can be either list or data.frame. List has an entry for every input gene. Data frame is just a two column table with inputs and corrsponding outputs. The input names may be duplicated.}
\references{J. Reimand, M. Kull, H. Peterson, J. Hansen, J. Vilo: g:Profiler - a web-based toolset for functional profiling of gene lists from large-scale experiments (2007) NAR 35 W193-W200}
\author{Juri Reimand <jyri.reimand@ut.ee>, Raivo Kolde <rkolde@gmail.com>, Tambet Arak <tambet.arak@gmail.com>}
\arguments{\item{query}{list of gene IDs.}
\item{organism}{organism name.}
\item{target}{target namespace.}
\item{region_query}{interpret query as chromosomal ranges.}
\item{numeric_ns}{namespace to use for fully numeric IDs.}
\item{df}{logical indicating whether the output will be a data.frame or list.}
}
\examples{gconvert(c("Klf4", "Pax5", "Sox2", "Nanog"), organism = "mmusculus")

# Get all mouse Cell cycle genes 
gconvert(c("GO:0007049"), organism = "mmusculus")}

