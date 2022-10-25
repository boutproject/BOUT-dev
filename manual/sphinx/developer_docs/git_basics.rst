.. Use bash as the default language for syntax highlighting in this file
.. highlight:: console

.. _sec-gitbasics:

Git crash course
----------------

The BOUT++ distribution is hosted on Github:

https://github.com/boutproject/BOUT-dev

For a full guide to using Git, see the `git website`_ or online
tutorials. This manual just explains some basic ways to use Git, and the
recommended work flow when working with BOUT++.

.. _git website: https://git-scm.com/

If you’re just starting with BOUT++, current developers will want to
check your changes before submitting them to the repository. In this
case you should fork the git repository, make any changes and then
submit a pull request on Github. Fortunately Git makes this process
quite easy: First get a copy of BOUT++::

    $ git clone https://github.com/boutproject/BOUT-dev.git

The BOUT++ repository will now be in a directory called “BOUT-dev”
(sorry - github doesn’t like ’+’ in project names). To get the latest
changes, use::

    $ git pull

To see the status of the repository, commits etc. in a GUI use::

    $ gitk

This is also useful for showing what changes you’ve made which need to
be committed, or which haven’t yet been sent to the main repository.

You can make edits as normal, and commit them using::

    $ git commit -a

which is pretty much the equivalent of ``svn commit`` in that it commits
all changes, though importantly it doesn’t send them to a central
server. To see which changes will be committed, use::

    $ git status

To choose which files you want to commit, use::

    $ git add file1, file2, ...
    $ git commit

(Git can actually only commit selected parts of files if you want). To
make using Git easier, you can create a config file ``$HOME/.gitconfig``
containing:

.. code-block:: cfg

    [user]
        name = A. Developer
        email = a.developer@example.com

    [alias]
        st = status
        ci = commit
        br = branch
        co = checkout
        df = diff
        lg = log -p
        who = shortlog -s --

(though obviously you should change the name and email).

Once you’re done making changes, you should first pull the latest
changes from the server::

    $ git pull

**Read carefully** what git prints out. If there are conflicts then git
will try to resolve them, but in some cases you will have to resolve
them yourself. To see a list of conflicting changes run ``git status``
(or ``git st`` if you’re using the above ``.gitconfig`` file). Once
you’ve finished resolving conflicts, run ``git commit -a`` to commit the
merge.

Accessing github from behind a firewall
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you’re working on a machine which can’t access github directly (such
as grendel, smaug etc. at LLNL), you can still seamlessly access github
by using another machine as a proxy over SSH. To do this, edit your SSH
config file `` /.ssh/config`` and add the following lines:

.. code-block:: aconf

    Host            gh
    HostName        github.com
    User            git
    ProxyCommand    ssh -q -x user@euclid.nersc.gov nc %h %p

where ``euclid.nersc.gov`` can be replaced by any machine you can access
which has netcat (``nc``) installed, and which can access github.com. If
you have set up a github account with SSH keys, you should now be able
to get a copy of BOUT++ by running::

    $ git clone gh:boutproject/BOUT-dev.git

Creating a private repository
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Whilst we would prefer it if improvements to BOUT++ were shared,
sometimes you might want to keep changes private for a while before
publishing them. Creating a private repository with Git is very simple,
because every clone of a repository is itself a repository. Git doesn’t
have the concept of a central repository, which can seem strange coming
from the world of SVN and CVS. What it means is that you can create your
own private repository anywhere you have access to. Sharing it with only
some people means as giving them read or write access to the repository
directory.

The following assumes you have a NERSC account and want to create a
private repository on Franklin. To apply this to a different machine
just replace ``franklin.nersc.gov`` with the machine you want to put the
repository on.

#. SSH to ``franklin.nersc.gov``, or wherever you want your repository::

           $ ssh username@franklin.nersc.gov
         

#. Create a “bare” Git repository by cloning a repository with the
   ``–bare`` option::

           $ cd ~
           $ git clone --bare git@github.com:boutproject/BOUT-dev.git  bout_private
         

   where you can replace ``git@github.com:boutproject/BOUT-dev.git`` with any
   other repository you can access. ``bout_private`` will be the name of
   the directory which will be created. This will make a repository
   without a working version. This means you can’t modify the code in it
   directly, but can pull and push changes to it. If you want to work on
   the code on Franklin, make a clone of your private repository::

           $ git clone bout_private bout
         

   which creates a repository ``bout`` from your private repository.
   Running ``git pull`` and ``git push`` from within this new repository
   will exchange patches with your ``bout_private`` repository.

#. You can now clone, pull and push changes to your private repository
   over SSH e.g.::

           $ git clone username@franklin.nersc.gov:bout_private
         

#. To keep your private repository up to date you may want to pull
   changes from github into your private repository. To do this, you
   need to use a third repository. Log into Franklin again::

           $ cd ~
           $ git clone bout_private bout_tmp
         

   This creates a repository ``bout_tmp`` from your private repository.
   Now cd to the new directory and pull the latest changes from github::

           $ cd bout_tmp
           $ git pull https://github.com/boutproject/BOUT-dev.git
         

   Note: You should be able to access this repository from Franklin, but
   if not then see the previous subsection for how to access github from
   behind a firewall.

#. This pull might result in some conflicts which need to be resolved.
   If so, git will tell you, and running::

           $ git status
         

   will give a list of files which need to be resolved. Edit each of the
   files listed, and when you’re happy commit the changes::

           $ git commit -a
         

#. Your ``bout_tmp`` directory now contains a merge of your private
   repository and the repository on github. To update your private
   repository, just push the changes back::

           $ git push
         

   You can now delete the ``bout_tmp`` repository if you want.

