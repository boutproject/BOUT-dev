#!/usr/bin/env python3
import argparse
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from ruamel.yaml import YAML
from unidecode import unidecode


@dataclass
class KnownAuthor:
    family_names: str
    given_names: str


TOP_DIR = Path(__file__).parent.parent
CITATION_FILE = TOP_DIR / "CITATION.cff"
NONHUMAN_AUTHORS = ("github", "dependabot")
KNOWN_AUTHORS = {
    "bendudson": KnownAuthor("Dudson", "Benjamin"),
    "brey": KnownAuthor("Breyiannis", "George"),
    "David Schwörer": KnownAuthor("Bold", "David"),
    "dschwoerer": KnownAuthor("Bold", "David"),
    "hahahasan": KnownAuthor("Muhammed", "Hasan"),
    "Ilon Joseph - x31405": KnownAuthor("Joseph", "Ilon"),
    "kangkabseok": KnownAuthor("Kang", "Kab Seok"),
    "loeiten": KnownAuthor("Løiten", "Michael"),
    "Michael Loiten Magnussen": KnownAuthor("Løiten", "Michael"),
    "Maxim Umansky - x26041": KnownAuthor("Umansky", "Maxim"),
    "nick-walkden": KnownAuthor("Walkden", "Nicholas"),
    "ZedThree": KnownAuthor("Hill", "Peter"),
    "tomc271": KnownAuthor("Chapman", "Tom"),
    "j-b-o": KnownAuthor("Bold", "Jessica"),
    "BS": KnownAuthor("Brendan", "Shanahan"),
}


def get_authors_from_git():
    output = subprocess.run(
        ["git", "log", "--format='%aN %aE'"],
        capture_output=True,
        check=True,
        text=True,
    )

    authors_string = output.stdout
    authors_list = authors_string.split("\n")
    authors_without_quotes = [a.strip("'") for a in authors_list]

    distinct_authors = set(authors_without_quotes)
    distinct_authors_list_without_empty_strings = [
        a.rsplit(maxsplit=1) for a in distinct_authors if a
    ]
    authors_with_emails = defaultdict(list)
    for author, email in distinct_authors_list_without_empty_strings:
        authors_with_emails[author].append(email)

    return authors_with_emails


def parse_cff_file():
    yaml = YAML(typ="rt")
    with open(CITATION_FILE, "r", encoding="UTF-8") as stream:
        return yaml.load(stream)


def get_authors_from_cff_file():
    file_contents = parse_cff_file()
    try:
        return file_contents["authors"]
    except KeyError as key_error:
        raise ValueError(f"Failed to find section:{key_error} in {CITATION_FILE}")


class ExistingAuthorNames:
    def __init__(self, existing_authors):
        self._existing_author_names = [
            (unidecode(a.get("given-names")), unidecode(a.get("family-names")))
            for a in existing_authors
        ]

    def last_name_matches_surname_and_first_name_or_first_letter_matches_given_name(
        self, last_name, first_name
    ):
        # Last name matches surname
        matches = [
            n
            for n in self._existing_author_names
            if n[1].casefold() == last_name.casefold()
        ]

        for match in matches:
            # The given name also matches author first name
            if match[0].casefold() == first_name.casefold():
                return True
            # The first initial matches author first name
            if match[0][0].casefold() == first_name[0].casefold():
                return True

    def first_name_matches_surname_and_last_name_matches_given_name(
        self, first_name, last_name
    ):
        # First name matches surname
        matches = [
            n
            for n in self._existing_author_names
            if n[1].casefold() == first_name.casefold()
        ]

        # The given name also matches author last name
        for match in matches:
            if match[0].casefold() == last_name.casefold():
                return True

    def surname_matches_whole_author_name(self, author):
        surname_matches = [
            n
            for n in self._existing_author_names
            if n[1].casefold() == author.casefold()
        ]
        if len(surname_matches) > 0:
            return True

    def given_name_matches_matches_whole_author_name(self, author):
        given_name_matches = [
            n
            for n in self._existing_author_names
            if n[0].casefold() == author.casefold()
        ]
        if len(given_name_matches) > 0:
            return True

    def combined_name_matches_whole_author_name(self, author):
        combined_name_matches = [
            n
            for n in self._existing_author_names
            if (n[0] + n[1]).casefold() == author.casefold()
        ]
        if len(combined_name_matches) > 0:
            return True

    def combined_name_reversed_matches(self, author):
        combined_name_reversed_matches = [
            n
            for n in self._existing_author_names
            if (n[1] + n[0]).casefold() == author.casefold()
        ]
        if len(combined_name_reversed_matches) > 0:
            return True

    def author_name_is_first_initial_and_surname_concatenated(self, author):
        first_character = author[0]
        remaining_characters = author[1:]
        # Second part of name matches surname
        matches = [
            n
            for n in self._existing_author_names
            if n[1].casefold() == remaining_characters.casefold()
        ]
        for match in matches:
            # The first initial matches author first name
            if match[0][0].casefold() == first_character.casefold():
                return True


def author_found_in_existing_authors(author, existing_authors):
    existing_author_names = ExistingAuthorNames(existing_authors)

    names = author.split()
    first_name = unidecode(names[0].replace(",", ""))
    last_name = unidecode(names[-1])

    if existing_author_names.last_name_matches_surname_and_first_name_or_first_letter_matches_given_name(
        last_name, first_name
    ):
        return True

    if existing_author_names.first_name_matches_surname_and_last_name_matches_given_name(
        first_name, last_name
    ):
        return True

    if existing_author_names.surname_matches_whole_author_name(author):
        return True

    if existing_author_names.given_name_matches_matches_whole_author_name(author):
        return True

    if existing_author_names.combined_name_matches_whole_author_name(author):
        return True

    if existing_author_names.combined_name_reversed_matches(author):
        return True

    if existing_author_names.author_name_is_first_initial_and_surname_concatenated(
        author
    ):
        return True

    return False


def is_nonhuman(author: str) -> bool:
    lower_case_author = author.casefold()
    return any(nonhuman in lower_case_author for nonhuman in NONHUMAN_AUTHORS)


def update_citations():
    authors_from_git = get_authors_from_git()
    existing_authors = get_authors_from_cff_file()

    known_authors = [a for a in authors_from_git if a in KNOWN_AUTHORS]
    human_authors = {a: e for a, e in authors_from_git.items() if not is_nonhuman(a)}
    authors_to_search_for = {
        a: e for a, e in human_authors.items() if a not in known_authors
    }

    unrecognised_authors = {
        a: e
        for a, e in authors_to_search_for.items()
        if not author_found_in_existing_authors(a, existing_authors)
    }

    if not unrecognised_authors:
        return

    for author, email in unrecognised_authors.items():
        print(f"{author} (email(s): {email})")

    while True:
        reply = (
            input("\nThe above authors were not recognised. Add to citations? [y/N] ")
            .lower()
            .strip()
        )
        if not reply or reply[0] == "n":
            return
        if reply[0] == "y":
            break

    new_authors = []
    for author in unrecognised_authors:
        first_name, last_name = author.rsplit(maxsplit=1)
        new_authors.append({"family-names": last_name, "given-names": first_name})

    print(f"Adding new authors:\n{new_authors}")

    yaml_file = parse_cff_file()
    yaml_file["authors"].extend(new_authors)

    # Try to preserve indentation and quotes as much as possible
    yaml = YAML()
    yaml.indent(mapping=2, sequence=4, offset=2)
    yaml.preserve_quotes = True
    with open(CITATION_FILE, "w", encoding="UTF-8") as f:
        yaml.dump(yaml_file, f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Update CITATIONS.cff based on git authors"
    )
    parser.parse_args()

    update_citations()
