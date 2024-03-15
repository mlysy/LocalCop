#' Implementation of [rticles::joss_article()] which is compatible with [bookdown::pdf_book()].
#'
#' @param journal,keep_md,latex_engine Same as for [rticles::joss_article()].
#' @param pandoc_args Character vector of arguments to pass to pandoc.
#' @param ... Same as for [rticles::joss_article()], except without `pandoc_args`.
#'
#' @return Same as for [rticles::joss_article()].
joss_article2 <- function(journal = "JOSS", keep_md = TRUE,
                          latex_engine = "xelatex",
                          pandoc_args = NULL, ...) {
  rmarkdown::pandoc_available("2.2", TRUE)
  logo_path <- find_resource("joss", paste0(journal, "-logo.png"))
  journalname <- ifelse(journal == "JOSS",
                        "Journal of Open Source Software",
                        "Journal of Open Source Education")
  pandoc_args <- c(
    pandoc_args,
    "-V", paste0("logo_path=", logo_path),
    "-V", paste0("journal_name=", journalname),
    "-V", "graphics=true"
  )
  pdf_document_format("joss", latex_engine = latex_engine,
                      citation_package = "default",
                      keep_md = keep_md,
                      pandoc_args = pandoc_args, ...)
}

environment(joss_article2) <- asNamespace("rticles")
assignInNamespace("joss_article", joss_article2,
                  ns = "rticles")

rticles::joss_article
