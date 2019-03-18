We'd love you to contribute, and here is the mini how-to:
talk, fork, branch, hack, rob & pull request.

Longer version:

1. Before diving into technicalities: if you intend to make major changes,
   beyond fixing bugs and small functionality improvements, please open a Github
   issue first, so we can discuss before coding. Please explain what you intend
   to accomplish and why. That often saves a lot of time and trouble in the long
   run.

   Use the issue to plan your changes. Try to solve only one problem at a time,
   instead of fixing several issues and adding different features in a single
   shot. Small changes are easier to handle, also for the reviewer in the last
   step below.

   Mention in the corresponding issue when you are working on it. "Claim" the
   issue to avoid duplicate efforts.

2. Install Roberto, which is the driver for our CI setup. It can also replicate
   the continuous integration on your local machine, which makes it easier to
   prepare a passable pull request. See TODO FIX URL.

3. Make a fork of the project, using the Github "fork" feature.

4. Clone the original repository on your local machine and enter the directory

   .. code-block:: bash

    git clone git@github.com:theochem/cellcutoff.git
    cd cellcutoff

5. Add your fork as a second remote to your local repository, for which we will
   use the short name `mine` below, but any short name is fine:

   .. code-block:: bash

    git remote add mine git@github.com:<your-github-account>/cellcutoff.git

6. Make a new branch, with a name that hints at the purpose of your
   modification:

   .. code-block:: bash

    git checkout -b new-feature

7. Make changes to the source. Please, make it easy for others to understand
   your code. Also, add tests that verify your code works as intended.
   Rules of thumb:

   - Write transparent code, e.g. self-explaining variable names.
   - Add comments to passages that are not easy to understand at first glance.
   - Write docstrings explaining the API.
   - Add unit tests when feasible.

8. Commit your changes with a meaningful commit message. The first line is a
   short summary, written in the imperative mood. Optionally, this can be
   followed by an empty line and a longer description.

   If you feel the summary line is too short to describe what you did, it
   may be better to split your changes into multiple commits.

9. Run Roberto and fix all problems it reports. Either one of the following
   should work

   .. code-block:: bash

    rob                 # Normal case
    python3 -m roberto  # Only if your PATH is not set correctly

   Style issues, failing tests and packaging issues should all be detected at
   this stage.

10. Push your branch to your forked repository on Github:

    .. code-block:: bash

        git push mine -u new-feature

    A link should be printed on screen, which will take the next step for you.

11. Make a pull request from your branch `new-feature` in your forked repository
    to the `master` branch in the original repository.

12. Wait for the tests on Travis-CI to complete. These should pass. Also
    coverage analysis will be shown, but this is merely indicative. Normally,
    someone should review your pull request in a few days. Ideally, the review
    results in minor corrections at worst. We'll do our best to avoid larger
    problems in step 1.
