# beiko-lab/ARETE: Contributing Guidelines

Hey! Thank you for taking an interest in contributing to ARETE.

We use GitHub for managing issues, contribution requests and everything else. So feel free to communicate with us using new issues and discussions, whatever best fits your idea for your contribution

## Contribution Workflow

The standard workflow for contributing to ARETE is as follows:

1. Check first if there isn't already an issue for your feature request, bug, etc. If there isn't one, **you should create a new [issue](https://github.com/beiko-lab/arete/issues/new/choose) or [discussion](https://github.com/beiko-lab/arete/discussions/new/choose) for your planned contribution before starting working on it**.
2. [Fork](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) the [beiko-lab/ARETE](https://github.com/beiko-lab/arete/) repository to your GitHub account
3. Make the necessary changes or additions within your forked repository. **You should probably create a new branch for your contribution, instead of committing directly to the master branch of your repository**.
4. In case any parameters were added or changed, use `nf-core schema build` to add them to the pipeline JSON schema and `nf-core schema docs --output docs/params.md --force` to update the respective documentation (requires [nf-core tools](https://github.com/nf-core/tools) >= 1.10).
5. Optionally run the pipeline's unit tests locally using [`nf-test`](https://github.com/askimed/nf-test#installation): `nf-test test tests/subworkflows/local/*`.
6. Submit a Pull Request to our `master` branch and wait for your changes to be reviewed and merged by our maintainers.

Our Github Actions workflows should perform a few pipeline tests automatically after receiving your pull request.
Any errors on them should be looked at since they could point to underlying issues in your changes.
