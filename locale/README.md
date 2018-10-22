
To start a new translation:

    msginit --input=libbout.pot --locale=fr --output=fr/libbout.po

Edit the .po file, then convert to a .mo file:
Note: Change the content line to: "Content-Type: text/plain; charset=UTF-8\n"

    msgfmt --output-file=LC_MESSAGES/libbout.mo libbout.po
