# GitHub Actions CI workflows
Definitions for GitHub Actions (continuous integration) workflows

## Publishing
To publish a new version of the `pymap3d` package to PyPI, create and publish a
release in GitHub (preferrably from a Git tag) for the version; the workflow
will automatically build and publish an sdist and wheel from the tag.

Requires the repo secret `PYPI_API_TOKEN` to be set to a PyPI API token: see
[PyPI's help](https://pypi.org/help/#apitoken) for instructions on how to
generate one.
