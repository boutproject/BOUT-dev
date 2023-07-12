from pathlib import Path
import re


def get_full_filepath(filepath):

    return Path(r"C:\git\BOUT-dev") / filepath


def update_version_number_in_file(relative_filepath, pattern, new_version_number):

    def get_replacement(match):
        return match[0].replace(match[1], new_version_number)

    full_filepath = get_full_filepath(relative_filepath)
    with open(full_filepath, "r", encoding='UTF-8') as file:
        file_contents = file.read()
        updated_text = re.sub(pattern, get_replacement, file_contents, flags=re.MULTILINE)
    with open(full_filepath, "w", encoding='UTF-8') as file:
        file.write(updated_text)


def bump_version_numbers(new_version_number):

    short_version_number = new_version_number[0:3]
    new_minor_version_number = str(int(new_version_number[2]) + 1)
    bout_next_version_number = re.sub(r"(?<=\d\.)\d(?=\.\d)", new_minor_version_number, new_version_number)

    update_version_number_in_file("configure.ac",
                                  r"^AC_INIT\(\[BOUT\+\+\],\[(\d\.\d\.\d)\]", new_version_number)
    update_version_number_in_file("CITATION.cff",
                                  r"^version: (\d\.\d\.\d)", new_version_number)
    update_version_number_in_file("manual/sphinx/conf.py",
                                  r"^version = \"(\d\.\d)\"", short_version_number)
    update_version_number_in_file("manual/sphinx/conf.py",
                                  r"^release = \"(\d\.\d\.\d)\"", new_version_number)
    update_version_number_in_file("manual/doxygen/Doxyfile_readthedocs",
                                  r"^PROJECT_NUMBER         = (\d\.\d\.\d)", new_version_number)
    update_version_number_in_file("manual/doxygen/Doxyfile",
                                  r"^PROJECT_NUMBER         = (\d\.\d\.\d)", new_version_number)
    update_version_number_in_file("CMakeLists.txt",
                                  r"^set\(_bout_previous_version \"v(\d\.\d\.\d)\"\)", new_version_number)
    update_version_number_in_file("CMakeLists.txt",
                                  r"^set\(_bout_next_version \"(\d\.\d\.\d)\"\)", bout_next_version_number)
    update_version_number_in_file("tools/pylib/_boutpp_build/backend.py",
                                  r"_bout_previous_version = \"v(\d\.\d\.\d)\"", new_version_number)
    update_version_number_in_file("tools/pylib/_boutpp_build/backend.py",
                                  r"_bout_next_version = \"v(\d\.\d\.\d)\"", bout_next_version_number)


if __name__ == '__main__':

    bump_version_numbers(new_version_number="6.1.0")
