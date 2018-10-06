SPiCCAto:
	cd src && $(MAKE)
	cp src/run-SPiCCAto .
clean:
	rm run-SPiCCAto && cd src && $(MAKE) clean
