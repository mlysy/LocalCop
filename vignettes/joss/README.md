# LocalCop JOSS Article

This folder contains the source code to reproduce the article accompanying the **LocalCop** package published in the Journal of Open-Source Software (JOSS).  The following are specific instructions for recreating the article submission itself.

1.  Run `rmarkdown::render("_LocalCop-paper.Rmd")` to create the Markdown file `_LocalCop-paper.md`.  This should be identical to `paper.md` in this folder, which is the exact Markdown used for the JOSS submission.

2.  The R Markdown command above generates a PDF output, but this isn't necessarily identical to the JOSS submission.  To obtain the former, one can run the [Open Journals PDF Generator](https://github.com/marketplace/actions/open-journals-pdf-generator) GitHub action.  

	The action should be tweaked to use the correct file `_LocalCop-paper.md`, since it will default to `paper.md` which won't find the figures in the correct location.
	
	Speaking of which, one would also need to un`.gitignore` the figure folder `_LocalCop-paper_files` in order for the GitHub action to properly render the PDF.
