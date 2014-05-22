all: manual.html

pages=wiki/Command-Line-Options.md

manual.html: ${pages}
	pandoc -s -S --toc -o manual.html ${pages}
