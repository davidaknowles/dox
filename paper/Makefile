NAME := dox_paper
BIB := dox.bib

all: $(NAME).pdf

$(NAME).pdf: $(NAME).tex $(BIB)
	pdflatex $(NAME)
	bibtex $(NAME)
	pdflatex $(NAME)
	pdflatex $(NAME)

.PHONY: clean

clean:
	rm -f *pdf *~ *.aux *.blg *.bbl *.lof *.log *.lot *.out *.toc
