===========================
Contributing
===========================

Want to contribute to developement? That's fantastic! Have a discussion with us first via `issues <https://github.com/xraysoftmat/kkcalc-xraysoftmat/issues>`_ to avoid wasted effort, and so we can make our visions align.

To install for development, create your own `fork/branch <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo>`_ of ``kkcalc2`` on ``github`` so you can save and contribute your own changes, and clone that to your local file system using ``Github Desktop`` or ``Git`` via Windows Terminal / CMD.

Follow the instructions in the :ref:`install` tab to set up an editable installation environment from the source code.

.. code-block:: shell

    # Install the library from the source code
    pip install -e . --group dev


Continuous Integration
======================
`kkcalc2` uses a combination of continuous integration (CI) tools to manage the code quality and documentation.

- `pre-commit <https://pre-commit.com/>`_: To fix consistency issues in code before pushing to the repository.
- `ruff <https://github.com/astral-sh/ruff>`_: To provide linting and code quality checks. Standardizes code formatting so the entire repository is consistent and readable. Previously used and style continuing close to `black <https://github.com/psf/black>`_.
- `numpydoc <https://numpydoc.readthedocs.io/>`_: Readable documentation to standard with other scientific python packages.
- `pytest <https://github.com/pytest/pytest>`_: A framework for writing and running tests on example calculations and results, so _unit_ functionality doesn't break between versions.

Some of these tools can be used as hooks - that is quality checks that are run before you commit or push your code. Install these hooks as follows:

.. code-block:: shell

    pre-commit install

Unit-testing
============
To run unit testing, invoke ``pytest``:

.. code-block:: shell

    pytest

Ideally, ensure all tests pass (and write new ones for new code) before creating a pull request.

Documentation
=============

To build the documentation, install the documentation dependencies from the package directory:

.. code-block:: shell

    pip install . --group docs
    # or
    pip install . --group dev # for everything

To build the documentation, run

.. code-block:: shell

    sphinx-build -M html ./docs ./docs/_build

which should create a new set of static `HTML` files in the `kkcalc2/docs/_build` directory.

.. To update the documentation after edits, run:
..     sphinx-apidoc -o ./docs ./kkcalc2
.. in the source directory, then rebuild as above.
