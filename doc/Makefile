all: manual.html
.PHONY: init-wiki update-wiki

# Documentation
# if HAVE_PANDOC
pages=wiki/Installation.md wiki/Command-Line-Options.md wiki/Output.md

init-wiki:
	git submodule init; git submodule update

update-wiki: init-wiki
	cd wiki; git pull origin master

manual.html: ${pages} knitr.css
	pandoc -t html5 -s -S --toc --toc-depth=2 \
	  -M title='The scrm Quick Reference' \
	  -M author='Paul R. Staab, Joe (Sha) Zhu, Dirk Metzler and Gerton Lunter' \
	  -M date="Version ${VERSION}" \
	  -c knitr.css \
	  --standalone \
	  -o manual.html ${pages}
#endif
