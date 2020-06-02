# Note: Tabs must be used instead of spaces
#       Using spaces results in a "missing separator" error
filename = optimiser-overview
all: clean pdf
pdf: clean
	pdflatex $(filename)    # Creates a bad pdf and lots of temp files
	bibtex $(filename)      # Creates filename.bbl, needed by latex
	pdflatex $(filename)    # Produces a slightly better pdf and temp files
	pdflatex $(filename)    # Produces a correct pdf from the temp files
	make clean
clean:
	rm -f *.aux *.log *.out *.spl *.blg *.bbl *.toc *.maf *.mtc*

