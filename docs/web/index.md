% The MFrontGenericInterfaceSupport project
% Thomas Helfer
% 20/11/2018

This project aims at proving tools (functions, classes, bindings,
etc...) to handle behaviours written using `MFront` generic interface.
For information about `MFront`, see <http://tfel.sourceforge.net>.

Those tools are meant to be used by (`FEM`, `FFT`, etc.) solver
developers. This tools are *not* linked to the `TFEL` libraries.
Permissive licences have been chosen to allow integration in open-source
and proprietary codes.

# The `MFrontGenericInterface` `C++` library

The project is build around the `MFrontGenericInterface` library. This
library provides two main functions:

- the `mgis::behaviour::load` functions loads `MFront` behaviours from
  external shared libraries and retrieve all relevant meta data
  function. Those relevant information are stored in the
  `mgis::behaviour::Behaviour` class.
- the `mgis::behaviour::integrate` integrates the behaviour over one
  time step. The data associated with an integration point are handled
  by the `mgis::behaviour::BehaviourData` class which contains the state
  of the integration point at the beginning and at the end of the time
  step.

The library also supports handling a group of integration points though
the `mgis::behaviour::MaterialStateManager` class.

An introduction to the `C++` library may be found [here](bindings-cxx.html)

# Bindings

## Existing bindings

The following bindings are available:

- `c` binding
- `python` binding: this binding also provide an interface to `FEniCS` which is 
  discussed [here along with a collection of commented demos](mgis_fenics.html)
- `fortran` binding
- `julia` binding

## Future bindings (contributors are welcomed)

The following bindings are under consideration:

- `octave` binding

# Versions, branches

## Branches

- the `master` branch follows the evolution of the `master` branch of
  the `TFEL` project
- the `rliv-1.1` follows the evolution of the 3.3.x series of the `TFEL`
  project.
- the `rliv-1.0` follows the evolution of the 3.2.x series of the `TFEL`
  project. Note that this branch is **not** compatible with
  `TFEL-3.2.0`.

## Release notes

- [Version 1.2](release-notes-1.2.html)
- [Version 1.1](release-notes-1.1.html)

# Installation

The installation process is discussed [here](install.html).

# Contributing to `MGIS`

Contributions to `MGIS` are greatly appreciated.

Please take a moment to review this document in order to make the
contribution process easy and effective for everyone involved.

Following these guidelines helps to communicate that you respect the time of
the developers managing and developing this open source project. In return,
they should reciprocate that respect in addressing your issue or assessing
patches and features.

## Using the issue tracker

The [issue
tracker](https://github.com/thelfer/MFrontGenericInterfaceSupport/issues)
is the preferred channel for [bug reports](#bugs), [features
requests](#features) and [submitting pull requests](#pull-requests), but
please respect the following restrictions:

* Please **do not** use the issue tracker for personal support requests
  (contact directly the authors
  [tfel-contact@cea.fr](mailto:tfel-contact@cea.fr)

## Bug reports

A bug is a _demonstrable problem_ that is caused by the code in the repository.
Good bug reports are extremely helpful - thank you!

Guidelines for bug reports:

1. **Use the GitHub issue search**: check if the issue has already been
   reported.

2. **Check if the issue has been fixed**: try to reproduce it using the
   latest `master` or development branch in the repository.

3. **Isolate the problem**: ideally create a [reduced test
   case].

A good bug report shouldn't leave others needing to chase you up for more
information. Please try to be as detailed as possible in your report. What is
your environment? What steps will reproduce the issue? What compiler(s) and OS
experience the problem? What would you expect to be the outcome? All these
details will help people to fix any potential bugs.

Example:

> Short and descriptive example bug report title
>
> A summary of the issue, versions of `TFEL/MFront` used and the
> OS/compiler environment in which it occurs. If suitable, include the
> steps required to reproduce the bug.
>
> 1. This is the first step
> 2. This is the second step
> 3. Further steps, etc.
>
> Any other information you want to share that is relevant to the issue being
> reported. This might include the lines of code that you have identified as
> causing the bug, and potential solutions (and your opinions on their
> merits).

## Feature requests

Feature requests are welcome. But take a moment to find out whether your idea
fits with the scope and aims of the project. It's up to *you* to make a strong
case to convince the project's developers of the merits of this feature. Please
provide as much detail and context as possible.

## Pull requests

Good pull requests - patches, improvements, new features - are a fantastic
help. They should remain focused in scope and avoid containing unrelated
commits.

**Please ask first** before embarking on any significant pull request (e.g.
implementing features, refactoring code, porting to a different language),
otherwise you risk spending a lot of time working on something that the
project's developers might not want to merge into the project.

Please adhere to the coding conventions used throughout a project (indentation,
accurate comments, etc.) and any other requirements (such as test coverage).

Adhering to the following this process is the best way to get your work
included in the project:

1. [Fork](http://help.github.com/fork-a-repo/) the project, clone your fork,
   and configure the remotes:

   ```bash
   # Clone your fork of the repo into the current directory
   git clone https://github.com/<your-username>/MFrontGenericInterfaceSupport.git
   # Navigate to the newly cloned directory
   cd MFrontGenericInterfaceSupport
   # Assign the original repo to a remote called "upstream"
   git remote add upstream https://github.com/thelfer/MFrontGenericInterfaceSupport.git
   ```

2. If you cloned a while ago, get the latest changes from upstream:

   ```bash
   git checkout master
   git pull upstream master
   ```

3. Create a new topic branch (off the main project development branch) to
   contain your feature, change, or fix:

   ```bash
   git checkout -b <topic-branch-name>
   ```

4. Commit your changes in logical chunks. Please adhere to these [git commit
   message guidelines](http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html)
   or your code is unlikely be merged into the main project. Use Git's
   [interactive rebase](https://help.github.com/articles/interactive-rebase)
   feature to tidy up your commits before making them public.

5. Locally merge (or rebase) the upstream development branch into your topic branch:

   ```bash
   git pull [--rebase] upstream master
   ```

6. Push your topic branch up to your fork:

   ```bash
   git push origin <topic-branch-name>
   ```

7. [Open a Pull Request](https://help.github.com/articles/using-pull-requests/)
    with a clear title and description.

**IMPORTANT**: By submitting a patch, you agree to allow the project owners to
license your work under the the terms of the *LGPL License*.

<!--
This paragraph is merely a copy of the `CONTRIBUTING.md` file of the
html5boilerplate project

Copyright (c) HTML5 Boilerplate

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
-->

# Acknowledgement

This project uses code extracted from the following projects:

- <https://github.com/bitwizeshift/string_view-standalone> by Matthew
  Rodusek
- <https://github.com/mpark/variant>: by Michael Park
- <https://github.com/progschj/ThreadPool> by Jakob Progsch and Václav
  Zeman
- <https://github.com/martinmoene/span-lite> by Martin Moene


