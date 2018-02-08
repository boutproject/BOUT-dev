
BOUT_TOP = ../..

SOURCEC		= field.cxx field2d.cxx field3d.cxx fieldperp.cxx field_data.cxx \
		  fieldgroup.cxx field_factory.cxx fieldgenerators.cxx \
		  initialprofiles.cxx vecops.cxx vector2d.cxx vector3d.cxx \
		  where.cxx globalfield.cxx generated_fieldops.cxx
SOURCEH		= $(SOURCEC:%.cxx=%.hxx) field_data.hxx
TARGET		= lib


include $(BOUT_TOP)/make.config

# The C++ source file depends on both the Python driver and the Jinja
# template files. If either have changed, run the driver (first
# dependency), dumping the output in the C++ target file. We then try
# to apply clang-format to this file in-place, but don't worry if the
# formatting fails
generated_fieldops.cxx: gen_fieldops.py gen_fieldops.jinja
	@echo "  Generating $@"
	@./$< --filename $@.tmp || (fail=$?; echo "touch $@ to ignore failed generation" ; exit $fail)
	@mv $@.tmp $@
	@clang-format -i $@ || echo "Formatting failed"
