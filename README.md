CEM
===

A version of *Coastline Evolution Model* with a *BMI*.

You can find the official version of CEM here:
https://github.com/csdms-contrib/cem.


Install a pre-built version of CEM
----------------------------------

To install a pre-built version of this package with conda run:

    $ conda install -c csdms cem

This will grab the latest version of *CEM* from the CSDMS channel on
[Anaconda.org](https://anaconda.org/csdms/cem) and install it into your
current environment.

Build and install from source
-----------------------------

To build and install from source:

    mkdir _build && cd _build
    cmake .. -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_BUILD_TYPE=Release
    make all -j4
    make install

Create a new release
--------------------

New releases are built and uploaded to
[Anaconda.org](https://anaconda.org/csdms/cem) whenever a new tag that starts
with the letter "v" is
[created and pushed to](https://git-scm.com/book/en/v2/Git-Basics-Tagging)
[GitHub](https://github.com/csdms-contrib/cem). As an example, the following
will cause a new release to be built,

    $ git tag v0.1.1
    $ git push --tags

A new release is created (*v0.1.1*) and the tag pushed to GitHub.
[Travis-CI](https://travis-ci.org/csdms-contrib/cem) notices the tagged commit,
and after building and testing the package, creates a fresh new package that
is uploaded to [Anaconda.org](https://anaconda.org/csdms/cem).

A couple notes about creating a new version:

1.  The version given in the tag name should match that in
    `conda-recipe/meta.yaml`.
1.  If you mess up (forget to update all the versions, for example), you can
    always [delete the tag and recreate it](https://git-scm.com/docs/git-tag).
    To do this, you'll need to delete both the remote tag and the local tag:

        $ git push --delete origin <tagname> # Delete the tag on the remote repository
        $ git tag --delete <tagname> # Delete the tag from the local repository
1.  If your new tag was successfully pushed to GitHub, you will be able to see
    it with the rest of the
    [releases](https://github.com/csdms-contrib/cem/releases) and
    [tags](https://github.com/csdms-contrib/cem/tags).
1.  To see if your new release was created successfully, you can do one or all
    of the following:

    *  Check the logs for the build of your tagged commit on
       [Travis-CI](https://travis-ci.org/csdms-contrib/cem).
    *  Check [Anaconda.org](https://anaconda.org/csdms/cem) to see if your
       releases appear there.
    *  Check if `conda` can see your new release with `conda search cem -c
       csdms`. See the [conda docs](http://conda.pydata.org/docs/using/index.html)
       for a description of `conda` and how to use it, or you can always use
       `conda -h` from the command line.

Helpful links
-------------

1.  [Using conda](http://conda.pydata.org/docs/using/index.html): What `conda`
    is and how to use it.
1.  [git tags](https://git-scm.com/book/en/v2/Git-Basics-Tagging): What git
    tags are and how to create them.
1.  [The git tag command](https://git-scm.com/docs/git-tag): A description
    of all of the options for the `git tag` command (including `git tag
    --delete`).
1.  [CEM on Travis](https://travis-ci.org/csdms-contrib/cem): The latest
    Travis builds of CEM.
1.  [CEM on Anaconda](https://anaconda.org/csdms/cem): The conda packages for
    CEM releases.
