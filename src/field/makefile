
BOUT_TOP = ../..

SOURCEC		= field.cxx field2d.cxx field3d.cxx fieldperp.cxx field_data.cxx \
		  fieldgroup.cxx field_factory.cxx fieldgenerators.cxx \
		  initialprofiles.cxx vecops.cxx vector2d.cxx vector3d.cxx \
		  where.cxx globalfield.cxx generated_fieldops.cxx
SOURCEH		= $(SOURCEC:%.cxx=%.hxx) field_data.hxx
TARGET		= lib


include $(BOUT_TOP)/make.config

generated_fieldops.cxx:gen_fieldops.py
	python3 $< > $@
	@clang-format -i $@ || echo "Formatting failed"
