from pathlib import Path
import os
import re


def get_full_filepath(filepath):

    main_directory = Path(os.path.abspath(__file__)).parent.parent
    return Path(main_directory) / filepath


def update_version_number_in_file(relative_filepath, pattern, new_version_number):

    def get_replacement(match):
        return match[0].replace(match[1], new_version_number.as_string())

    full_filepath = get_full_filepath(relative_filepath)
    with open(full_filepath, "r", encoding='UTF-8') as file:
        file_contents = file.read()
        updated_text = re.sub(pattern, get_replacement, file_contents, flags=re.MULTILINE)
    with open(full_filepath, "w", encoding='UTF-8') as file:
        file.write(updated_text)


def bump_version_numbers(new_version_number):

    short_version_number = ShortVersionNumber(new_version_number.major_version, new_version_number.minor_version)
    bout_next_version_number = VersionNumber(new_version_number.major_version,
                                             new_version_number.minor_version + 1,
                                             new_version_number.patch_version)

    update_version_number_in_file("configure.ac",
                                  r"^AC_INIT\(\[BOUT\+\+\],\[(\d+\.\d+\.\d+)\]", new_version_number)
    update_version_number_in_file("CITATION.cff",
                                  r"^version: (\d+\.\d+\.\d+)", new_version_number)
    update_version_number_in_file("manual/sphinx/conf.py",
                                  r"^version = \"(\d+\.\d+)\"", short_version_number)
    update_version_number_in_file("manual/sphinx/conf.py",
                                  r"^release = \"(\d+\.\d+\.\d+)\"", new_version_number)
    update_version_number_in_file("manual/doxygen/Doxyfile_readthedocs",
                                  r"^PROJECT_NUMBER         = (\d+\.\d+\.\d+)", new_version_number)
    update_version_number_in_file("manual/doxygen/Doxyfile",
                                  r"^PROJECT_NUMBER         = (\d+\.\d+\.\d+)", new_version_number)
    update_version_number_in_file("CMakeLists.txt",
                                  r"^set\(_bout_previous_version \"v(\d+\.\d+\.\d+)\"\)", new_version_number)
    update_version_number_in_file("CMakeLists.txt",
                                  r"^set\(_bout_next_version \"(\d+\.\d+\.\d+)\"\)", bout_next_version_number)
    update_version_number_in_file("tools/pylib/_boutpp_build/backend.py",
                                  r"_bout_previous_version = \"v(\d+\.\d+\.\d+)\"", new_version_number)
    update_version_number_in_file("tools/pylib/_boutpp_build/backend.py",
                                  r"_bout_next_version = \"v(\d+\.\d+\.\d+)\"", bout_next_version_number)


class VersionNumber:

    major_version: int
    minor_version: int
    patch_version: int

    def __init__(self, major_version, minor_version, patch_version):
        self.major_version = major_version
        self.minor_version = minor_version
        self.patch_version = patch_version

    def as_string(self):
        return "%d.%d.%d" % (self.major_version, self.minor_version, self.patch_version)


class ShortVersionNumber(VersionNumber):

    def __init__(self, major_version, minor_version):
        self.major_version = major_version
        self.minor_version = minor_version
        self.patch_version = None

    def as_string(self):
        return "%d.%d" % (self.major_version, self.minor_version)


if __name__ == '__main__':

    bump_version_numbers(new_version_number=VersionNumber(63, 15, 12))
