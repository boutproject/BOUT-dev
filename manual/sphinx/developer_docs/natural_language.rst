Natural Language Support
========================

BOUT++ uses GNU ``gettext`` to provide translations of output strings.
Configuration is described in :ref:`sec-config-nls` and running in
:ref:`sec-run-nls`.

Adding support for a new language, or improving the translations in
the existing files can be done by:

1. Going to the ``locale`` BOUT++ subdirectory and running::

         make locale-ll

   where ``ll`` is the language code e.g. ``make locale-de``. This
   will create a file ``libbout.po`` under a ``locale/ll``
   subdirectory.
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
