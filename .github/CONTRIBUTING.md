# Contributing Guidelines for `spaghetti`

Thank you for your interest in contributing! We work primarily on Github. Please review the contributing procedures [here](http://pysal.org/getting_started#for-developers) and [here](https://github.com/pysal/pysal/wiki/GitHub-Standard-Operating-Procedures) so that we can accept your contributions! Alternatively, contact someone in the [development chat channel](https://gitter.im//pysal/Spaghetti).


## Contact

1. All questions, comments, & discussions should happen in a public forum, where possible. Please start a [discussion](https://github.com/pysal/spaghetti/discussions) for questions or open an [issue](https://github.com/pysal/spaghetti/issues) if there appears to be a bug. Private messages and emails will not be answered in a substantive manner.

## Style and format

1. Python 3.7, 3.8, 3.9, and 3.10 are the officially supported versions.
2. This project follows the formatting conventions of [`black`](https://black.readthedocs.io/en/stable/) and utilizes [`pre-commit`](https://pre-commit.com) to format commits prior to pull requests being made. 
    * LJ Miranda provides an [excellent, concise guide](https://ljvmiranda921.github.io/notebook/2018/06/21/precommits-using-black-and-flake8/) on setting up and implementing a `pre-commit` hook for `black`.
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


