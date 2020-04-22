# common make routines and rules

# automatic prerequisite generation for c files
# NOTE doesn't include system headers
%.d: %.c
	@set -e; rm -f $@ && \
	$(CC) -MM $(CPPFLAGS) $< > $@.$$$$ && \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@ && \
	rm -f $@.$$$$

