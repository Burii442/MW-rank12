# MW-rank12

Exact computations accompanying Dominik Burek, *An elliptic Calabi–Yau
threefold of Mordell–Weil rank twelve*.

The program checks the explicit formulas and finite computations used in the
paper. It runs independently of the manuscript. The 24 rows of the node table
are stored in `data/node_table.json`; each row gives a value of `v`, a choice of
`tau`, and sections through the four points with `t = tau, i*tau, -tau, -i*tau`.
Every coefficient is rational, and computations use exact arithmetic, including
arithmetic over `Q(i)`.

## Reproducing the computations

With an existing SageMath installation, open a terminal in this directory and run:

```sh
sage -python verify.py --output ../MW-rank12-results/results.json
```

The checks use SymPy (tested with version 1.14.0). A successful run prints `PASS`; the full
computational report is saved at the path above. No manuscript file is needed.

Alternatively, use Python 3.9 or later with the pinned SymPy dependency:

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python verify.py --output ../MW-rank12-results/results.json
```

A successful run ends with:

```text
PASS: 395 identities, two height matrices, 96 node incidences and nonzero graph Jacobians; e(X)=48.
```

The JSON report records every identity, both height matrices and their principal
minors, all 96 node incidences and graph Jacobians, the numerical invariants, and
SHA-256 hashes of the program and its input data. It also records the Python and
SymPy versions. Reports are generated afresh outside the repository; none are
distributed as precomputed results.

To compare the stored table with the marked table in a manuscript source, use:

```sh
python verify.py --paper /path/to/MW_rank12.tex --output ../MW-rank12-results/results-with-paper.json
```

This additionally checks all table entries against the manuscript and records
the source hash. It does not check that all other formulas printed in the
manuscript agree with the formulas encoded in the program.

The exit status is zero only when every check passes. A failed mathematical or
table-validation check produces a `FAIL` report and a nonzero exit status. Reports
are replaced atomically; a new run first marks its report `RUNNING`, so an earlier
`PASS` cannot be mistaken for a completed new run. Reports must be outside the
repository and cannot overwrite a supplied manuscript. Checks remain active
with Python's `-O` option.

## Relation to the paper

| Part of the paper | Computations |
| --- | --- |
| §2, the elliptic curves and their invariants | The identity `A^3 - C^2 = 108 B^4`. |
| §4, the crepant resolution of the quartic | The lowest weighted forms at the six indicated singularities. |
| §5, the elliptic pencil | The binary-quartic invariants, the pointed quartic-to-Weierstrass identity, and the change of scale. |
| §6, the Weierstrass model | The transition identities for its coefficients. |
| §7, ordinary double points | The factorization and squarefreeness of `J`, coprimality checks, node identities, all 96 node incidences, and nonzero section-graph Jacobians. |
| §8, Hodge numbers | Discriminant identities and intersection, genus and Euler arithmetic using the geometric inputs of the paper. |
| Appendix A | Six generic section identities, their twelve specializations and boundary coefficients; pairwise intersections and height matrices at `v=5,7`, positive leading principal minors, and determinant `4096/9`. |
| Appendix B, Table 1 | All 24 rows, representing 96 nodes; optional comparison with the manuscript source. |

The reported 395 identities are the polynomial equalities checked by the
program. The squarefreeness, positivity, nonvanishing and table-size checks are
additional checks and are not included in that count.

The computations supplement the proofs. The program does not independently
prove the global resolution, flatness, the Calabi–Yau property, the generic
Picard number, or the upper bound for the Mordell–Weil rank. Intersection
corrections and Euler contributions use the geometric arguments in the paper.
In particular, `h11=25` is supplied from the manuscript; `h21=1` is then derived
from `h11` and `e(X)=48`. These are distinguished in the report.

## Repository contents

- `verify.py`: the computation and optional manuscript-table comparison.
- `data/node_table.json`: exact node-table data, with a description of its format.
- `requirements.txt`: the pinned dependency.
- `tests/test_cli.py`: input-validation and failure-report checks.
- `.github/workflows/verify.yml`: automatic verification with Python 3.12.
- `CITATION.cff`: citation metadata for version 1.0.1.

Run the input-validation tests with `python -m unittest discover -s tests`.
The mathematical verification is run separately by `verify.py`.

GitHub Actions runs both on each update. Repository users with write access can
start a fresh run under **Actions → Exact verification → Run workflow**, selecting
`main`. A green **Success** means the complete computational verifier and its
input-validation tests passed. The report remains temporary; the workflow does
not upload it as an artifact.

Only verification code, required mathematical inputs, instructions and citation
metadata are maintained here. Manuscripts, conversations, working notes and saved
execution outputs are excluded.

## Citation

Dominik Burek, *MW-rank12: exact computations for an elliptic Calabi–Yau
threefold of Mordell–Weil rank twelve*, version 1.0.1 (2026).

Repository: https://github.com/Burii442/MW-rank12
