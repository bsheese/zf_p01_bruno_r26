# Manuscript build.
#   make            → manuscript.pdf
#   make clean      → remove the built PDF
#
# Assembles the section files in order, resolves @citekey citations against
# ../lit/references.bib (Chicago author-date via pandoc's built-in citeproc),
# and renders to PDF with xelatex. Figures are read from ../output/figures.
# Requires: pandoc (>=3), a TeX install with xelatex, and the DejaVu Serif font
# (broad Unicode coverage for α, ×, ↔, → used in the text).

PANDOC  := pandoc
BIB     := ../lit/references.bib
OUT     := manuscript.pdf
SRC     := metadata.yaml intro.md methods.md results.md discussion.md references.md

# -implicit_figures: render images in place with manual "Figure N." captions,
#   so figure/table numbers in the text stay fixed (no LaTeX float reordering).
FORMAT  := markdown-implicit_figures

FLAGS := --from=$(FORMAT) \
         --citeproc --bibliography=$(BIB) \
         --pdf-engine=xelatex \
         --number-sections \
         -V mainfont="DejaVu Serif" \
         -V geometry:margin=1in \
         -V fontsize=11pt \
         -V colorlinks=true -V linkcolor=Maroon -V citecolor=Maroon -V urlcolor=NavyBlue \
         -V linkReferences=true

$(OUT): $(SRC) $(BIB)
	$(PANDOC) $(SRC) $(FLAGS) -o $(OUT)
	@echo "Built $(OUT)"

.PHONY: clean
clean:
	rm -f $(OUT)
