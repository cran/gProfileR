
gp_globals = new.env()

gp_globals$version = 
	sessionInfo("gProfileR")$otherPkgs$gProfileR$Version
gp_globals$rcurl_opts =
	RCurl::curlOptions(useragent = paste("gProfileR/", gp_globals$version, sep=""))
gp_globals$png_magic =
	as.raw(c(0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a))
gp_globals$base_url = 
	"http://biit.cs.ut.ee/gprofiler/"

#' Find orthologs.
#' 
#' Interface to the g:Orth tool. Organism names are constructed by concatenating the first letter of 
#' the name and the family name. Example: human - 'hsapiens', mouse - 'mmusculus'.
#' 
#' To alleviate the problem of having many orthologs per gene (most of 
#' them uninformative) one can set a threshold for the number of results. The program 
#' tries to find the most informative by selecting the most popular ones. 
#'
#' @param query list of gene IDs to be translated.
#' @param source_organism name of the source organism.
#' @param target_organism name of the target organism. 
#' @param region_query interpret query as chromosomal ranges.
#' @param numeric_ns namespace to use for fully numeric IDs.
#' @param mthreshold maximum number of ortholog names per gene to show.
#' @param df logical indicating whether the output will be a data.frame or list.
#' @return The output can be either a list or a data.frame. The list has an entry for every input gene.
#' The data frame is just a two column table with inputs and corrsponding outputs. The input names may be duplicated.
#' @references  J. Reimand, M. Kull, H. Peterson, J. Hansen, J. Vilo: g:Profiler -- a web-based toolset for functional profiling of gene lists from large-scale experiments (2007) NAR 35 W193-W200 
#' @author  Raivo Kolde <rkolde@@gmail.com>, Juri Reimand <juri.reimand@@ut.ee>, Tambet Arak <tambet.arak@@gmail.com>
#' @examples
#' gorth(c("Klf4", "Pax5", "Sox2", "Nanog"), source_organism = "mmusculus", target_organism = "hsapiens")
#' @export

gorth <- function(
	query,
	source_organism = "hsapiens",
	target_organism = "mmusculus",
	region_query = F,
	numeric_ns = "",
	mthreshold = 3,
	df = T
){
	my_url = paste(gp_globals$base_url, "gorth.cgi", sep="")
	
	if (length(query) == 0) return(NULL)
	
	raw_query <- RCurl::postForm(my_url, .opts = gp_globals$rcurl_opts,
		output = 'mini', 
		query = paste(query, collapse = " "), 
		organism = source_organism,
		target = target_organism,
		region_query = ifelse(region_query, "1", "0"),
		prefix = numeric_ns
	)
	
	conn <- textConnection(raw_query)
	tab <- read.table(conn, sep = "\t")
	close(conn)
	
	colnames(tab)[2] <- "Source"
	
	res = plyr::dlply(tab, "Source", function(x){
		
		x = x[x$V6 != "N/A",]
		if (nrow(x) == 0){
			return(NULL)
		}
		
		if (length(x$V6) < mthreshold + 1){
			return(as.character(x$V6))
		} 
		else{
			return(names(rev(sort(table(as.character(x$V6)))))[1:mthreshold])
		}
	})
	
	if(df){
		res <- plyr::ldply(res, function(x) data.frame(Target = x, stringsAsFactors = F))
		res$Source = as.character(res$Source)
	}
	
	return(res)
}
 
#' Annotate gene list functionally.
#' 
#' Interface to the g:Profiler tool for finding enrichments in gene lists. Organism names are 
#' constructed by concatenating the first letter of the name and the family name. Example: human - 
#' 'hsapiens', mouse - 'mmusculus'. If requesting PNG output, the request is directed to the g:GOSt 
#' tool in case 'query' is a vector and the g:Cocoa (compact view of multiple queries) tool in case 
#' 'query' is a list. PNG output can fail (return FALSE) in case the input query is too large. In 
#' such case, it is advisable to fall back to a non-image request.
#'
#' @param organism organism name.
#' @param query vector of gene IDs or a list of such vectors.
#' @param ordered_query in case output gene lists are ranked this option may be used
#' to get GSEA style p-values.
#' @param significant whether all or only statistically significant results should be 
#' returned.
#' @param exclude_iea exclude electronic annotations (IEA).
#' @param region_query interpret query as chromosomal ranges. 
#' @param max_p_value custom p-value threshold, results with a larger p-value are excluded.
#' @param max_set_size maximum size of functional category, larger categories are excluded.
#' @param correction_method the algorithm used for determining the significance threshold, one of 
#' "gSCS", "fdr", "bonferroni".
#' @param hier_filtering hierarchical filtering strength, one of "none", "moderate", "strong".
#' @param domain_size statistical domain size, one of "annotated", "known".
#' @param custom_bg vector of gene names to use as a statistical background.
#' @param numeric_ns namespace to use for fully numeric IDs. 
#' @param png_fn request the result as PNG image and write it to png_fn.
#' @return A data frame with the enrichment analysis results. If the input consisted 
#' of several lists the corresponding list is indicated with a variable 'query number'.
#' When requesting a PNG image, either TRUE or FALSE, depending on whether a non-empty result 
#' was received and a file written or not, respectively. 
#' @references  J. Reimand, M. Kull, H. Peterson, J. Hansen, J. Vilo: g:Profiler - a web-based toolset for functional profiling of gene lists from large-scale experiments (2007) NAR 35 W193-W200
#' @author  Juri Reimand <jyri.reimand@@ut.ee>, Raivo Kolde <rkolde@@gmail.com>, Tambet Arak <tambet.arak@@gmail.com>
#' @examples
#' gprofiler(c("Klf4", "Pax5", "Sox2", "Nanog"), organism = "mmusculus")
#' @export

gprofiler <- function(
	query, 
	organism='hsapiens', 	
	ordered_query=F,
	significant=T, 	
	exclude_iea=F,
	region_query = F,
	max_p_value=1.0,
	max_set_size=0,
	correction_method="analytical",
	hier_filtering="none",
	domain_size="annotated",
	custom_bg="",	
	numeric_ns = "",
	png_fn = NULL
) {	
	query_url = ""
	my_url = paste(gp_globals$base_url, "gcocoa.cgi", sep="")
	wantpng = ifelse(is.character(png_fn), T, F)
	output_type = ifelse(wantpng, "mini_png", "mini")
	
	# Query

	if (is.list(query)) {
		qnames = names(query)
		
		for (i in 1:length(query)) {
			cname = paste("> query", i)
			if (!is.null(qnames) && nchar(qnames[i]) > 0)
				cname = paste(">", qnames[i])
			query_url <- paste(sep="\n", query_url, cname, paste(query[[i]], collapse=" "))
		}
	} else if (is.vector(query)) {
		query_url <- paste(query, collapse=" ")		
		if (wantpng)
			my_url = paste(gp_globals$base_url, "index.cgi", sep="")		
	} else {
		warning("Missing query")
		return()
	}
	
	# Significance threshold
	
	if (correction_method == "gSCS")
		correction_method = "analytical"
	if (!correction_method %in% c("analytical", "fdr", "bonferroni"))
		stop("Multiple testing correction method not recognized")
	
	# Hierarchical filtering
	
	if (!hier_filtering %in% c("none", "moderate", "strong"))
		stop("hier_filtering must be one of \"none\", \"moderate\" or \"strong\"")
	if (hier_filtering == "strong") hier_filtering = "compact_ccomp"
	else if (hier_filtering == "moderate") hier_filtering = "compact_rgroups"
	else hier_filtering = ""
	
	# Domain size

	if (!domain_size %in% c("annotated", "known"))
		stop("domain_size must be one of \"annotated\" or \"known\"")
	
	# Custom background
	
	if (is.vector(custom_bg))
		custom_bg = paste(custom_bg, collapse=" ")
	else
		stop("custom_bg must be a vector")

	# Max. set size
	
	if (max_set_size < 0) max_set_size = 0
	
	# HTTP request
	
	raw_query <- RCurl::postForm(my_url, .opts = gp_globals$rcurl_opts,
		organism=organism, 
		query=query_url, 
		output=output_type, 		
		analytical="1", 
		sort_by_structure="1",
		ordered_query = ifelse(ordered_query, "1", "0"),
		significant = ifelse(significant, "1", "0"),
		no_iea = ifelse(exclude_iea, "1", "0"),
		as_ranges = ifelse(region_query, "1", "0"),		
		user_thr = as.character(max_p_value),
		max_set_size = as.character(max_set_size),
		threshold_algo = correction_method,
		hierfiltering = hier_filtering,		
		domain_size_type = domain_size,
		custbg_file = "",
		custbg = custom_bg,
		prefix = numeric_ns
	)
	
	# Requested PNG, write to disk and return

	if (wantpng) {	
		pngv = as.vector(raw_query)
		
		if (length(raw_query) == 0)	
			return(F)
		if (!isTRUE(all.equal(head(pngv, 8), gp_globals$png_magic))) { # check PNG magic number
			warning(paste(
				"Received non-PNG data, but PNG output requested. Please report this error",
				"with the appropriate details to the package maintainer"
			));
			return(F);
		}
		
		conn = file(png_fn, "wb")	
		writeBin(pngv, conn, useBytes = T)
		close(conn)
		return(T);
	}
	
	# Requested text
	
	split_query <- unlist(strsplit(raw_query, split="\n"))
	
	commented_lines <- grep("^#", split_query)
	if (length(commented_lines)>0) {
		split_query <- split_query[-commented_lines]
	}
	
	empty_lines <- which(split_query == "")
	if (length(empty_lines)>0) {
		split_query <- split_query[-empty_lines]
	}
	
	if(length(split_query) > 0){
		conn <- textConnection(paste(split_query, collapse = "\n"))
		split_query <- read.table(conn, sep = "\t", quote = "", stringsAsFactors = F)
		close(conn)
	}
	else{
		split_query <- as.data.frame(matrix(NA, 0, 14))
	}
	
	rownames(split_query) <- NULL
	colnames(split_query) <- c(
		"query.number", "significant", "p.value", 
		"term.size", "query.size", "overlap.size", 
		"precision", "recall", "term.id", 
		"domain", "subgraph.number", "term.name",
		"relative.depth", "intersection"
	)
	split_query$term.name <- gsub("^\\s+", "", split_query$term.name)	
	split_query$significant <- ifelse(split_query$significant == "!", T, F)
	
	if(is.list(query) & !is.null(names(query))){
		split_query$query.number <- names(query)[split_query$query.number]
	}
	
	return(split_query)
}

#' Convert gene IDs.
#' 
#' Interface to the g:Convert tool. Organism names are constructed by concatenating the first letter of 
#' the name and the family name. Example: human - 'hsapiens', mouse - 'mmusculus'.
#'
#' @param query list of gene IDs.
#' @param organism organism name.
#' @param target target namespace.
#' @param region_query interpret query as chromosomal ranges. 
#' @param numeric_ns namespace to use for fully numeric IDs.
#' @param df logical indicating whether the output will be a data.frame or list.
#' @return The output can be either list or data.frame. List has an entry for every input gene. Data frame is just a two column table with inputs and corrsponding outputs. The input names may be duplicated.
#' @references  J. Reimand, M. Kull, H. Peterson, J. Hansen, J. Vilo: g:Profiler - a web-based toolset for functional profiling of gene lists from large-scale experiments (2007) NAR 35 W193-W200
#' @author  Juri Reimand <jyri.reimand@@ut.ee>, Raivo Kolde <rkolde@@gmail.com>, Tambet Arak <tambet.arak@@gmail.com>
#' @examples
#' gconvert(c("POU5F1", "SOX2", "NANOG"), organism = "hsapiens", target="AFFY_HG_U133_PLUS_2")
#' @export

gconvert = function(
	query,
	organism = "hsapiens",
	target = "ENSG",
	region_query = F,
	numeric_ns = "",	
	df = T
) {
	
	url = paste(gp_globals$base_url, "gconvert.cgi", sep="")		
	field = 4

	if (target == "DESC") {
		field = 6
		target = "ENSG"
	}

	raw_query <- RCurl::postForm(url, .opts = gp_globals$rcurl_opts,
		organism=organism, 
		output='mini', 
		query=paste(query, collapse="+"), 
		target=target,
		region_query = ifelse(region_query, "1", "0"),
		prefix = numeric_ns
	)
	
	split_query = unlist(strsplit(raw_query, split="\n"))
	
	commented_lines = grep("^#", split_query)
	if (length(commented_lines)>0) {
		split_query = split_query[-commented_lines]
	}
	empty_lines = which(split_query == "")
	if (length(empty_lines)>0) {
		split_query = split_query[-empty_lines]
	}
	
	split_query = t(as.data.frame(strsplit(split="\t", split_query)))
	rownames(split_query) = c()	
	resulting_mapping = list()

	if (nrow(split_query)) {
		for (id in unique(split_query[,2])) {
			my_map = split_query[split_query[,2]==id,,drop=F][,field]
			my_map = my_map[my_map != "N/A"]
			resulting_mapping[[id]] = my_map
		}
	}
	
	if(df){
		resulting_mapping <- plyr::ldply(resulting_mapping, function(x) data.frame(Target = x, stringsAsFactors = F))
	}
	
	return(resulting_mapping)
}

#' Get current user agent string.
#' 
#' Get the HTTP User-Agent string.
#'
#' @export

get_user_agent = function() {
	gp_globals$rcurl_opts$useragent
}

#' Set custom user agent string.
#' 
#' Set the HTTP User-Agent string. Useful for overriding the default user agent for
#' packages that depend on gProfileR functionality.
#'
#' @param ua the user agent string.
#' @param append logical indicating whether to append the passed string to the default user agent string.
#' @export

set_user_agent = function(ua, append = F) {	
	rco = gp_globals$rcurl_opts
	rco$useragent = ifelse(
		append, paste(rco$useragent, ua, sep=""), ua)
	assign("rcurl_opts", rco, envir=gp_globals)
}
