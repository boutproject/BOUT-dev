all: derivs.hxx


derivs.hxx: derivs.hxx.in.py derivs.hxx.in.jinja
# run first dependency and write output to target
# on failure delete output and report error
	@echo "  Generating $@"
	@./$< > $@.tmp || (fail=$?; echo "touch $@ to ignore failed generation" ; exit $fail)
	@mv $@.tmp $@
	@clang-format -i $@ || echo "Formatting failed"
