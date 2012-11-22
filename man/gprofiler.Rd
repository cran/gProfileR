\name{gprofiler}
\alias{gprofiler}
\title{Annotate gene list functionally.}
\usage{gprofiler(query, organism="hsapiens", ordered_query=F, significant=T,
    exclude_iea=F, region_query=F, max_p_value=1, max_set_size=0,
    correction_method="analytical", hier_filtering="none",
    domain_size="annotated", custom_bg="", numeric_ns="")
}
\description{Annotate gene list functionally.}
\details{Interface to the g:Profiler tool for finding enrichments in gene lists. Organism names are constructed by concatenating the first letter of 
the name and the family name. Example: human - 'hsapiens', mouse - 'mmusculus'.}
\value{Data frame with the enrichment analysis results. If the input consisted of several lists the corresponding list is indicated with a variable 'query number'.}
\references{J. Reimand, M. Kull, H. Peterson, J. Hansen, J. Vilo: g:Profiler - a web-based toolset for functional profiling of gene lists from large-scale experiments (2007) NAR 35 W193-W200}
\author{Juri Reimand <jyri.reimand@ut.ee>, Raivo Kolde <rkolde@gmail.com>, Tambet Arak <tambet.arak@gmail.com>}
\arguments{\item{organism}{organism name.}
\item{query}{list of gene IDs or a list of such lists.}
\item{ordered_query}{in case output gene lists are ranked this option may be used
to get GSEA style p-values.}
\item{significant}{whether all or only statistically significant results should be 
returned.}
\item{exclude_iea}{exclude electronic annotations (IEA).}
\item{region_query}{interpret query as chromosomal ranges.}
\item{max_p_value}{custom p-value threshold, results with a larger p-value are excluded.}
\item{max_set_size}{maximum size of functional category, larger categories are excluded.}
\item{correction_method}{the algorithm used for determining the significance threshold, one of 
"gSCS", "fdr", "bonferroni".}
\item{hier_filtering}{hierarchical filtering strength, one of "none", "moderate", "strong".}
\item{domain_size}{statistical domain size, one of "annotated", "known".}
\item{custom_bg}{vector of gene names to use as a statistical background.}
\item{numeric_ns}{namespace to use for fully numeric IDs.}
}
\examples{gprofiler(c("Klf4", "Pax5", "Sox2", "Nanog"), organism = "mmusculus")}

