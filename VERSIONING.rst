===============================
Grid Cell repository versioning
===============================

This is a git managed repository. Assuming the central repository is on
BitBucket_, you can obtain and work with the repository in the following way.

1. If you don't know much about git, the best way is to spend an hour or two
   (possibly more :-)) reading the first three chapters of the `Git Book`_.
   This will be useful in the case here, but it is also a useful skill when
   collaborating in any kind of project that contains a human readable source
   code.

2. In case you just want to get the code, you have to install git_ on your
   system (supports Linux, MaxOS and Windows)

3. To download the repository, run the command::

        $ git clone https://<yourname>@bitbucket.org/lsolanka/gridcells.git

   Replace ``<yourname>`` with your Bitbucket_ username. You should also have
   at least a read-only access to the repository. If this runs succesfully, the
   command will create a directory named gridcells with the source code in it.

4. The repository contains a number of branches. The *master* branch is the
   most stable one. After you have cloned the repository, git will
   automatically check out the master branch for you.

5. After this, it is good to crate your own branch. For instance if you are
   planning to add resonance into the model, create a branch named
   ``resonance``::

        $ git checkout -b resonance

   This command creates a new branch and automatically switches your context to
   that branch. You can now work inside this branch and make any changes, and
   they won't be affecting the master branch.

6. I strongly recomment going through the `Git Book`_ and keep versioning the
   changes you have made to the model. It might be useful when for instance we
   want to publish the resonance or any other models that have been created.
   Also versioning allows you to keep track of your parameter values and go
   back when something goes wrong and the model no longer works.

7. The next step is to read README to see what the directory structure of the
   repository is, and INSTALL to see how to run the simulations.


.. _Bitbucket: https://bitbucket.org/
.. _Git Book: http://git-scm.com/book
.. _git: http://git-scm.com/
