---
name: Bug report
about: Create a report to help us improve
title: ''
labels: bug
assignees: jGaboardi

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Additional context**
Add any other context about the problem here.

**System & Versions**
- Platform information:
```python
>>> import os; print(os.name, os.sys.platform);print(os.uname())
```
- Python version: 
```python
>>> import sys; print(sys.version)
```
- spaghetti version:
```python
>>> import spaghetti; print(spaghetti.__version__)
```
- libpysal version:
```python
>>> import libpysal; print(libpysal.__version__)
```
- SciPy version:
```python
>>> import scipy; print(scipy.__version__)
```
- NumPy version:
```python
>>> import numpy; print(numpy.__version__)
```
- GeoPandas verion (if applicable):
```python
>>> import geopandas; print(geopandas.__version__)
```

Also, please upload any relevant data as [a file
attachment](https://help.github.com/articles/file-attachments-on-issues-and-pull-requests/). Please **do not** upload pickled objects, since it's nearly impossible to troubleshoot them without replicating your exact namespace. Instead, provide the minimal subset of the data required to replicate the problem. If it makes you more comfortable submitting the issue, feel free to:

1. remove personally identifying information from data or code
2. provide only the required subset of the full data or code
