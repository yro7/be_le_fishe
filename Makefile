.PHONY: all html open clean

all: html

html:
	./render.R

open:
	open rapport/rapport.html

clean:
	rm -f rapport/rapport.html rapport/*.knit.md

