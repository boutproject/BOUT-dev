Introduction
============

.. toctree::
   :maxdepth: 1
   :caption: Contents:

This is a manual describing the core BOUT++ [1]_, and is intended for
anyone who wants to work on improving BOUT++. It does its best to
describe the details of how BOUT++ works, and assumes that the user is
very comfortable with C++. For a general introduction, and
instructions for using BOUT++ see the users’ guide. The user’s guide
assumes only minimal knowledge of C++, and provides only those details
needed to use BOUT++.

Since BOUT++ is a scientific code, it is constantly changing and
(hopefully) being improved. This provides a moving target for
documentation, and means that describing all the details of how BOUT++
works in one document is probably impossible. This is particularly
true since often our first priority is to write papers and code - not
documentation - and so whatever is documented is likely to be slightly
out of date. A source of up-to-date documentation of the BOUT++ code
is the comments and Doxygen tags: running ``doxygen`` on the source
should produce a set of HTML documentation. See www.doxygen.org for
more details.

.. [1] http://www.sciencedirect.com/science/article/B6TJ5-4VTCM95-3/2/ed200cd23916d02f86fda4ce6887d798

Using the BOUT++ repository
---------------------------

The BOUT++ distribution is hosted on Github:

https://github.com/boutproject/BOUT-dev

For a full guide to using Git, see the `git website`_ or online
tutorials. This manual just explains some basic ways to use Git, and the
recommended work flow when working with BOUT++.

.. _git website: http://git-scm.com/

If you’re just starting with BOUT++, current developers will want to
check your changes before submitting them to the repository. In this
case you should fork the git repository, make any changes and then
submit a pull request on Github. Fortunately Git makes this process
quite easy: First get a copy of BOUT++

.. code-block:: bash

    $ git clone https://github.com/boutproject/BOUT-dev.git

The BOUT++ repository will now be in a directory called “BOUT-dev”
(sorry - github doesn’t like ’+’ in project names). To get the latest
changes, use

.. code-block:: bash

    $ git pull

To see the status of the repository, commits etc. in a GUI use

.. code-block:: bash

    $ gitk

This is also useful for showing what changes you’ve made which need to
be committed, or which haven’t yet been sent to the main repository.

You can make edits as normal, and commit them using

.. code-block:: bash

    $ git commit -a

which is pretty much the equivalent of ``svn commit`` in that it commits
all changes, though importantly it doesn’t send them to a central
server. To see which changes will be committed, use

.. code-block:: bash

    $ git status

To choose which files you want to commit, use

.. code-block:: bash

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
changes from the server:

.. code-block:: bash

    $ git pull

**Read carefully** what git prints out. If there are conflicts then git
will try to resolve them, but in some cases you will have to resolve
them yourself. To see a list of conflicting changes run ``git status``
(or ``git st`` if you’re using the above ``.gitconfig`` file). Once
you’ve finished resolving conflicts, run ``git commit -a`` to commit the
merge.


Developing BOUT++
~~~~~~~~~~~~~~~~~

If you are doing a lot of development of BOUT++, it will probably make
sense for you to push changes directly to the online repository. In this
case you’ll need to sign up for an account on ``github.com``, then
upload an ssh key and ask to be added. The process is then almost
identical except that you clone using SSH:

.. code-block:: bash

    $ git clone https://github.com/boutproject/BOUT-dev.git

and rather than creating a patch, you push changes to the repository:

.. code-block:: bash

    $ git push

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
to get a copy of BOUT++ by running

.. code-block:: bash

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

#. SSH to ``franklin.nersc.gov``, or wherever you want your repository

   .. code-block:: bash

           $ ssh username@franklin.nersc.gov
         

#. Create a “bare” Git repository by cloning a repository with the
   ``–bare`` option:

   .. code-block:: bash

           $ cd ~
           $ git clone --bare git@github.com:boutproject/BOUT-dev.git  bout_private
         

   where you can replace ``git@github.com:boutproject/BOUT-dev.git`` with any
   other repository you can access. ``bout_private`` will be the name of
   the directory which will be created. This will make a repository
   without a working version. This means you can’t modify the code in it
   directly, but can pull and push changes to it. If you want to work on
   the code on Franklin, make a clone of your private repository:

   .. code-block:: bash

           $ git clone bout_private bout
         

   which creates a repository ``bout`` from your private repository.
   Running ``git pull`` and ``git push`` from within this new repository
   will exchange patches with your ``bout_private`` repository.

#. You can now clone, pull and push changes to your private repository
   over SSH e.g.

   .. code-block:: bash

           $ git clone username@franklin.nersc.gov:bout_private
         

#. To keep your private repository up to date you may want to pull
   changes from github into your private repository. To do this, you
   need to use a third repository. Log into Franklin again:

   .. code-block:: bash

           $ cd ~
           $ git clone bout_private bout_tmp
         

   This creates a repository ``bout_tmp`` from your private repository.
   Now cd to the new directory and pull the latest changes from github:

   .. code-block:: bash

           $ cd bout_tmp
           $ git pull git://github.com/boutproject/BOUT-dev.git
         

   Note: You should be able to access this repository from Franklin, but
   if not then see the previous subsection for how to access github from
   behind a firewall.

#. This pull might result in some conflicts which need to be resolved.
   If so, git will tell you, and running

   .. code-block:: bash

           $ git status
         

   will give a list of files which need to be resolved. Edit each of the
   files listed, and when you’re happy commit the changes

   .. code-block:: bash

           $ git commit -a
         

#. Your ``bout_tmp`` directory now contains a merge of your private
   repository and the repository on github. To update your private
   repository, just push the changes back:

   .. code-block:: bash

           $ git push
         

   You can now delete the ``bout_tmp`` repository if you want.

House rules
-----------

BOUT++ consists of about 60,000 lines of C/C++  [2]_ along with 18,500
lines of IDL and 12,000 of Python. Of this, about 40,000 lines is the
core BOUT++ code, and the remainder a mix of pre- and post-processors,
and physics modules. As production codes go, this is not particularly
huge, but it is definitely large enough that keeping the code ‘clean’
and understandable is necessary. This is vital if many people are going
to work on the code, and also greatly helps code debugging and
verification. There are therefore a few house rules to keep in mind when
modifying the BOUT++ code.

When modifying the core BOUT++ code, please keep in mind that this
portion of the code is intended to be general (i.e. independent of any
particular physical system of equations), and to be used by a wide range
of users. Making code clear is also more important in this section than
the physics model since the number of developers is potentially much
greater.

Here are some rules for editing the core BOUT++ code:

- **NO FORTRAN**. EVER. Though it may be tempting for scientific
   programmers to use a little Fortran now and then, please please
   don’t put any into BOUT++.  Use of Fortran, particularly when mixed
   with C/C++, is the cause of many problems in porting and modifying
   codes.

-  If a feature is needed to study a particular system, only include it
   in the core code if it is more generally applicable, or cannot be put
   into the physics module.

Coding conventions
------------------

See CONTRIBUTING.md for guidelines on naming and other coding
conventions used within BOUT++.


.. [2]
   generated using Al Danial’s cloc

