"""Ask a yes/no question and return the answer.

"""

from builtins import input
import sys


def query_yes_no(question, default="yes"):
    """Ask a yes/no question via input() and return their answer.

    Answers are case-insensitive.

    Probably originally from https://code.activestate.com/recipes/577058/
    via https://stackoverflow.com/a/3041990/2043465

    Parameters
    ----------
    question : str
        Question to be presented to the user
    default : {"yes", "no", None}
        The presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    Returns
    -------
    bool
        True if the answer was "yes" or "y", False if "no" or "n"
    """

    valid = {"yes":True,   "y":True,  "ye":True,
             "no":False,     "n":False,   "No":False,  "N":False }

    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "\
                             "(or 'y' or 'n').\n")
