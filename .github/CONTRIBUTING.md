# Contributing Guidelines for `spaghetti`

Thank you for your interest in contributing! We work primarily on Github. Please review the contributing procedures [here](http://pysal.org/getting_started#for-developers) and [here](https://github.com/pysal/pysal/wiki/GitHub-Standard-Operating-Procedures) so that we can accept your contributions! Alternatively, contact someone in the [development chat channel](https://gitter.im//pysal/Spaghetti).


## Contact

1. All questions, comments, & discussions should happen in a public forum, where possible. Please start a [discussion](https://github.com/pysal/spaghetti/discussions) for questions or open an [issue](https://github.com/pysal/spaghetti/issues) if there appears to be a bug. Private messages and emails will not be answered in a substantive manner.


## Style and format

1. At the time of this writing, Python 3.10, 3.11, and 3.12 are the officially supported versions.
2. This project implements the linting and formatting conventions of [`ruff`](https://docs.astral.sh/ruff/) on all incoming Pull Requests. To ensure a PR is properly linted and formatted prior to creating a Pull Request, [install `pre-commit`](https://pre-commit.com/#installation) in your development environment and then [set up the configuration of pre-commit hooks](https://pre-commit.com/#3-install-the-git-hook-scripts). 
3. Import packages, classes, and functions with their full name where possible.
  * For example:
    
    :white_check_mark:
    ```python
    import spaghetti
    from shapely.geometry import Point
    ```
    :x:
    ```python
    import spaghetti as spgh
    from shapely.geometry import Point as pt
    ```


