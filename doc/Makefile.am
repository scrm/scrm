all: doc
.PHONY: update-wiki init-wiki doc

pages=wiki/Command-Line-Options.md

doc: init-wiki manual.html

init-wiki:
	cd ..; git submodule init; git submodule update

update-wiki: init-wiki
	cd wiki; git pull origin master

manual.html: ${pages} Makefile knitr.css
	pandoc -s -S --toc --toc-depth=2 \
	  -M title='scrm manual' \
	  -M author='Paul R. Staab, Joe (Sha) Zhu, Dirk Metzler and Gerton Lunter' \
	  -c knitr.css \
	  --standalone \
	  -o manual.html ${pages}
