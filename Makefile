.PHONY: all html open clean

all: html

html:
	./render.R

open:
	open src/rapport/rapport.html

clean:
	rm -f src/rapport/rapport.html src/rapport/*.knit.md

