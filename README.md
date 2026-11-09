# Moving Apolar Constraints on Hodge Cycles of the Fermat Sextic Fourfold

[![verify](https://github.com/Li-Hongmin/fermat-sextic-moving-apolar/actions/workflows/verify.yml/badge.svg)](https://github.com/Li-Hongmin/fermat-sextic-moving-apolar/actions/workflows/verify.yml) [![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.24054740-blue)](https://doi.org/10.5281/zenodo.24054740)

**Hongmin Li** (Institute of Science Tokyo; The University of Tokyo) · ORCID [0000-0003-0228-0600](https://orcid.org/0000-0003-0228-0600)

Main result of the manuscript: on the Fermat sextic fourfold, exact second-order calculations along a pencil of primitive classes from two Aoki–Shioda surfaces separate the difference class (obstructed) from the sum class (extends through a formal Hodge jet of order six), with the transported apolar module computed explicitly.

- Manuscript on rStation: <https://rstation.ai/problems/hodge_conjecture/papers/41/v1> (version 1)
- Archived on Zenodo: concept DOI [10.5281/zenodo.24054740](https://doi.org/10.5281/zenodo.24054740) (all versions); this repository mirrors version 1, DOI [10.5281/zenodo.24054741](https://doi.org/10.5281/zenodo.24054741)

## Contents

- `paper/paper.md`, `paper/fermat-sextic-moving-apolar.pdf` — the manuscript (version 1).
- `verification/` — the supplement of that version, unpacked unchanged; its own README describes each file.
- `.github/workflows/verify.yml` — reruns every check below on each push.

## What the programs check

| Check | Command (run inside `verification/`) |
|---|---|
| First-order pencil | `sage calculations/bundles/aoki-sextic-second-order/aoki-sextic-pencil-first-order.sage` |
| Second-order difference | `sage calculations/bundles/aoki-sextic-second-order/aoki-sextic-difference-second-order.sage` |
| Sum jet to order six | `sage calculations/bundles/aoki-sextic-second-order/aoki-sextic-sum-order-six.sage` |
| Moving apolar S2 | `sage calculations/bundles/aoki-sextic-moving-apolar-s2/aoki-sextic-moving-apolar-s2.sage` |
| Hilbert tangent at y0 | `Y0_EXACT=1 sage calculations/bundles/aoki-sextic-y0-hilbert-tangent/aoki-y0-hilbert-tangent.sage` |
| Singular: central union special fibre | `sage -singular -q calculations/bundles/aoki-sextic-ci-special-fiber/aoki-sextic-ci-special-fiber.sing < /dev/null` |

The programs check the finite and exact computations used in the manuscript; the mathematical arguments are those written in the paper.

## Reproducing

```sh
cd verification
sage calculations/bundles/aoki-sextic-second-order/aoki-sextic-pencil-first-order.sage
sage calculations/bundles/aoki-sextic-second-order/aoki-sextic-difference-second-order.sage
sage calculations/bundles/aoki-sextic-second-order/aoki-sextic-sum-order-six.sage
sage calculations/bundles/aoki-sextic-moving-apolar-s2/aoki-sextic-moving-apolar-s2.sage
Y0_EXACT=1 sage calculations/bundles/aoki-sextic-y0-hilbert-tangent/aoki-y0-hilbert-tangent.sage
sage -singular -q calculations/bundles/aoki-sextic-ci-special-fiber/aoki-sextic-ci-special-fiber.sing < /dev/null
```

Requirements: SageMath. Run the commands from a copy if you want to keep the recorded outputs untouched.

## Cite

```bibtex
@misc{li2026fermatsexticmovingapolar,
  author    = {Li, Hongmin},
  title     = {Moving Apolar Constraints on Hodge Cycles of the Fermat Sextic Fourfold},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.24054740},
  url       = {https://doi.org/10.5281/zenodo.24054740}
}
```

## License

Manuscript and data: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Code: MIT (see `LICENSE`).
