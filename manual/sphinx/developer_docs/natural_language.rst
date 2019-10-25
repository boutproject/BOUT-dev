Natural Language Support
========================

BOUT++ uses `GNU gettext <https://www.gnu.org/software/gettext/>`_. to
provide translations of output strings. Configuration is described in
:ref:`sec-config-nls` and running in :ref:`sec-run-nls`. Currently
only ``fr``, ``zh_TW``, and ``zh_CN`` have been added, but it is quite
easy to add more. See ``locale/README.md`` or below. 


Marking strings for translation
-------------------------------

In the code strings are wrapped with ``_()`` e.g. ``"hello world"``
becomes ``_("hello world")``. Find a string you want to replace (which
can include formatting like ``%d``), surround it with ``_()``. Then in
the locale directory:: 

  make libbout.pot

will update the template file ``libbout.pot`` under
``BOUT_TOP/locale``. The template file should not be edited, but is
used to generate language-specific files (``libbout.po``).
To update these language files see the next section.

Adding translations
-------------------

Adding support for a new language, or improving the translations in
the existing files can be done by:

1. Going to the ``locale`` BOUT++ subdirectory and running::

         make locale-ll

   where ``ll`` is the language code e.g. ``make locale-zh_TW`` or
   ``make locale-de``. This will create a file ``libbout.po`` under a
   ``locale/ll`` subdirectory.
2. Edit the ``locale/ll/libbout.po`` file. Edit the .po file in de
   subdirectory (not the .pot file!), adding the translations. Each
   ``msgid`` entry should have a translated ``msgstr`` entry. If you
   don't want to translate them all, just delete the ones you don't
   translate. Any missing will just revert to the version in the
   code. If you're adding UTF-8 characters, change the content line in
   the .po file to have charset=UTF-8.
3. In the ``locale`` directory run ``make``. This should output
   something like::
     
         Building language:  fr
         Building language:  zh_CN
         Building language:  zh_TW

The new language should now be available (no need to recompile BOUT++).
