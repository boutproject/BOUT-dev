from pathlib import Path
from unidecode import unidecode
import yaml

authors_from_git = ["A Allen", "Aaron Fisher", "Adam Dempsey", "Andrew Allen", "arkabokshi", "Ben Dudson", "bendudson",
                    "Benjamin Dudson", "Brendan", "Brendan Shanahan", "Brett Friedman", "brey", "BS", "bshanahan",
                    "Chenhao Ma", "Chris MacMackin", "D Dickinson", "David", "David Bold", "David Dickinson",
                    "David Schwörer", "dependabot[bot]", "Dmitry Feliksovich Meyerson", "Dmitry Meyerson", "dschwoerer",
                    "Eric Medwedeff", "Erik Grinaker", "Fabio Riva", "George Breyiannis", "GitHub Merge Button",
                    "github-actions[bot]", "hahahasan", "Haruki Seto", "Haruki SETO", "holger", "Holger Jones",
                    "Hong Zhang", "Ilon Joseph", "Ilon Joseph - x31405", "Jarrod Leddy", "JB Leddy", "j-b-o",
                    "Jed Brown", "Jens Madsen", "John Omotani", "johnomotani", "jonesholger", "Joseph Parker",
                    "Joshua Sauppe", "Kab Seok Kang", "kangkabseok", "Kevin Savage", "Licheng Wang", "loeiten",
                    "Luke Easy", "M Leconte", "M V Umansky", "Marta Estarellas", "Matt Thomas", "Maxim Umansky",
                    "Maxim Umansky - x26041", "Michael Løiten", "Michael Loiten Magnussen", "Minwoo Kim",
                    "Nicholas Walkden", "Nick Walkden", "nick-walkden", "Olivier Izacard", "Pengwei Xi", "Peter Hill",
                    "Peter Naylor", "Qin, Yining", "Sajidah Ahmed", "Sanat Tiwari", "Sean Farley", "seanfarley",
                    "Seto Haruki", "Simon Myers", "Tianyang Xia", "Toby James", "tomc271", "Tongnyeol Rhee",
                    "Xiang Liu", "Xinliang Xu", "Xueqiao Xu", "Yining Qin", "ZedThree", "Zhanhui Wang"]


def parse_cff_file(filename):
    with open(filename, "r", encoding='UTF-8') as stream:
        try:
            return yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)


def get_authors_from_cff_file():
    filename = Path(r"C:\git\BOUT-dev") / "CITATION.cff"
    file_contents = parse_cff_file(filename)
    try:
        return file_contents["authors"]
    except KeyError as key_error:
        print("Failed to find section:", key_error, "in", filename)


class ExistingAuthorNames:

    def __init__(self, existing_authors):
        self._existing_author_names = \
            [(unidecode(a.get("given-names")), unidecode(a.get("family-names"))) for a in existing_authors]

    def last_name_matches_surname_and_first_name_or_first_letter_matches_given_name(self, last_name, first_name):
        matches = [n for n in self._existing_author_names if
                   n[1].casefold() == last_name.casefold()]  # Last name matches surname

        for match in matches:
            if match[0].casefold() == first_name.casefold():  # The given name also matches author first name
                return True
            if match[0][0].casefold() == first_name[0].casefold():  # The first initial matches author first name
                return True

    def first_name_matches_surname_and_last_name_matches_given_name(self, first_name, last_name):
        matches = [n for n in self._existing_author_names if
                   n[1].casefold() == first_name.casefold()]  # First name matches surname

        for match in matches:
            if match[0].casefold() == last_name.casefold():  # The given name also matches author last name
                return True

    def surname_matches_whole_author_name(self, author):

        surname_matches = [n for n in self._existing_author_names if n[1].casefold() == author.casefold()]
        if len(surname_matches) > 0:
            return True

    def given_name_matches_matches_whole_author_name(self, author):
        given_name_matches = [n for n in self._existing_author_names if n[0].casefold() == author.casefold()]
        if len(given_name_matches) > 0:
            return True

    def combined_name_matches_whole_author_name(self, author):

        combined_name_matches = [n for n in self._existing_author_names if
                                 (n[0] + n[1]).casefold() == author.casefold()]
        if len(combined_name_matches) > 0:
            return True

    def combined_name_reversed_matches(self, author):

        combined_name_reversed_matches = [n for n in self._existing_author_names if
                                          (n[1] + n[0]).casefold() == author.casefold()]
        if len(combined_name_reversed_matches) > 0:
            return True

    def author_name_is_first_initial_and_surname_concatenated(self, author):
        first_character = author[0]
        remaining_characters = author[1:]
        matches = [n for n in self._existing_author_names if
                   n[1].casefold() == remaining_characters.casefold()]  # Second part of name matches surname
        for match in matches:
            if match[0][0].casefold() == first_character.casefold():  # The first initial matches author first name
                return True


def author_found_in_existing_authors(author, existing_authors):
    existing_author_names = ExistingAuthorNames(existing_authors)

    names = author.split()
    first_name = unidecode(names[0].replace(",", ""))
    last_name = unidecode(names[-1])

    if existing_author_names.last_name_matches_surname_and_first_name_or_first_letter_matches_given_name(
            last_name, first_name):
        return True

    if existing_author_names.first_name_matches_surname_and_last_name_matches_given_name(first_name, last_name):
        return True

    if existing_author_names.surname_matches_whole_author_name(author):
        return True

    if existing_author_names.given_name_matches_matches_whole_author_name(author):
        return True

    if existing_author_names.combined_name_matches_whole_author_name(author):
        return True

    if existing_author_names.combined_name_reversed_matches(author):
        return True

    if existing_author_names.author_name_is_first_initial_and_surname_concatenated(author):
        return True

    return False


def update_citations():

    existing_authors = get_authors_from_cff_file()

    nonhuman_authors = [a for a in authors_from_git if "github" in a.casefold() or "dependabot" in a.casefold()]

    known_authors = [a for a in authors_from_git if a in KNOWN_AUTHORS]

    human_authors = [a for a in authors_from_git if a not in nonhuman_authors]

    authors_to_search_for = [a for a in human_authors if a not in known_authors]

    unrecognised_authors = [a for a in authors_to_search_for if
                            not author_found_in_existing_authors(a, existing_authors)]

    print("The following authors were not recognised. Add to citations?")
    for author in unrecognised_authors:
        print(author)


KNOWN_AUTHORS = {"bendudson",
                 "brey",
                 "David Schwörer",
                 "dschwoerer",
                 "hahahasan",
                 "Ilon Joseph - x31405",
                 "kangkabseok",
                 "loeiten",
                 "Michael Loiten Magnussen",
                 "Maxim Umansky - x26041",
                 "nick-walkden",
                 "ZedThree",
                 # "tomc271"
                 }

if __name__ == '__main__':
    update_citations()
