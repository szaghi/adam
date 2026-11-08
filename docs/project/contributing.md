---
title: Contributing
---

# Contributing

This is a FOSS project — anyone interested in using, developing, or contributing
is welcome. The project follows a KISS (Keep It Simple and Stupid) philosophy.

## Reporting Issues

- Open a ticket on the repository **GitHub Issues** page
- Clearly describe the problem, including steps to reproduce for bugs
- Note the earliest version you know has the issue

## Pull Requests

ADAM is developed trunk-based on `develop`: there are no long-lived feature branches.

1. Fork the repository on GitHub and branch from `develop`:
   ```bash
   git checkout -b my-change develop
   ```
2. Run the regression suites of the applications your change can reach, as CI does:
   ```bash
   ./src/tests/prism/regression/run.sh cpu                         # PRISM, CPU
   ./src/tests/flume/regression/run.sh cpu                         # FLUME, CPU
   ./src/tests/prism/regression/run.sh fnl --varset local_nvf      # PRISM, GPU (OpenACC); CI has no GPU
   REGRESSION_CASES="rmf rmf-amr" ./src/tests/prism/regression/run.sh cpu --no-build   # only some cases
   ```
   A change to a numerical result must also update the affected goldens, with the reason in the commit message.
3. Check for unnecessary whitespace: `git diff --check`
   A change to the documentation ships with the code it describes. When it names build modes, check them:
   `python3 scripts/check_docs_modes.py` (CI runs it); a change to the fobos files regenerates the
   [build-modes page](/reference/build-modes) with `python3 scripts/make_build_modes_page.py`.
4. Open a pull request against `develop`, with a [Conventional Commits](#commit-style) message

## Fortran Coding Style

- **Clarity over brevity**: `real :: gas_ideal_air` is better than `real :: gia`
- Single-character variable names only for loop counters
- Name all constants
- `implicit none` in every module and program
- Declare `intent` for all procedure arguments, ordered: pass arg → `inout` → `in` → `out` → optional
- Indent with three spaces (not tabs); keep lines within 132 columns
- No trailing whitespace; blank lines must contain no spaces
- Use `>, <, ==` instead of `.gt., .lt., .eq.`
- Avoid Windows-style CRLF line endings

### Recommended git whitespace settings

```ini
[color]
  ui = true
[color "diff"]
  whitespace = red reverse
[core]
  whitespace = fix,-indent-with-non-tab,trailing-space,cr-at-eol
```

## Commit style

Use [Conventional Commits](https://www.conventionalcommits.org/) so that `CHANGELOG.md` is generated automatically from the git log:

| Prefix | Purpose | Changelog section |
|--------|---------|-------------------|
| `feat:` | New feature or capability | New features |
| `fix:` | Bug fix | Bug fixes |
| `perf:` | Performance improvement | Performance |
| `refactor:` | Code restructuring | Refactoring |
| `docs:` | Documentation only | Documentation |
| `test:` | Tests | Testing |
| `build:` | Build system | Build system |
| `ci:` | CI/CD pipeline | CI/CD |
| `chore:` | Maintenance | Miscellaneous |

Give the part of the code as the scope, `type(scope): description`. Append `!` for breaking changes (`feat!:`,
`fix!:`). Reference issues with `#123`; they are auto-linked.

```
fix(prism): a non-positive output cadence disables that output instead of crashing on mod(it, 0) (#97)
perf(lib/fnl): collapse the RK, eikonal and L2-norm nests like their OpenMP twins (#95 P4)
test(flume): the vortex-rk3 regression case, the only one on the low-storage Runge-Kutta path (#95)
```

---

## Creating a release

Releases are fully automated via `scripts/bump.sh` and GitHub Actions. The only steps needed are:

```bash
# Install git-cliff once (bump.sh needs the git-cliff binary on PATH)
cargo install git-cliff

# Then, to release:
scripts/bump.sh patch   # v1.2.3 → v1.2.4
scripts/bump.sh minor   # v1.2.3 → v1.3.0
scripts/bump.sh major   # v1.2.3 → v2.0.0
scripts/bump.sh v2.1.0  # explicit version
```

`bump.sh` will ask for confirmation, then:

1. Regenerate `CHANGELOG.md` from the git log via [git-cliff](https://git-cliff.org/)
2. Commit with `chore(release): vX.Y.Z`
3. Create an annotated git tag
4. Push commit + tag

Pushing the tag triggers the release workflow (`.github/workflows/release.yml`). It runs no test and builds no
documentation. It:
- packages a versioned tarball `adam-vX.Y.Z.tar.gz` (without `.git`, `exe`, `obj`, `mod` and the docs build);
- extracts the tag's section of `CHANGELOG.md` as the release notes;
- publishes a GitHub release with the tarball and `scripts/install.sh` attached.

The other workflows:

| Workflow | Trigger | What it does |
|---|---|---|
| `ci.yml` | every branch push and pull request (not tags) | the coverage rule (`fobis rule --ex makecoverage-analysis`); the PRISM and FLUME CPU regression suites (`--varset ci_gnu`); the FLUME OpenMP-vs-serial bitwise check; the RIPPLE CPU oracles; the check of the build modes named in the docs |
| `docs.yml` | a push to the default branch, or manual | the coverage rule and the coverage badge, the build-modes page, the API reference, the VitePress build, and the deploy to GitHub Pages |
| `install.yml` | after a successful release | installs the release through `install.sh` with FoBiS, with fpm, and with CMake when the tag has a `CMakeLists.txt` |

Coverage is not uploaded to any external service: `docs.yml` writes it to the badge `docs/public/coverage-badge.json`.
