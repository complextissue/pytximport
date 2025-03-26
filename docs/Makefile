.DEFAULT_GOAL := build-docs

SPHINXOPTS    ?=
SOURCEDIR     = source
BUILDDIR      = build

help:
	sphinx-build -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# cp ../README.md ./source/introduction.md
build-docs:
	@./install.sh
	make html

.PHONY: help Makefile

reset-build:
	if [ -d ./build/doctrees ]; then \
        rm -r ./build/doctrees; \
    fi && \
    if [ -d ./build/html ]; then \
        rm -r ./build/html; \
    fi && \
    mkdir -p ./build/doctrees && \
    mkdir -p ./build/html

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	make reset-build
	sphinx-build -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
