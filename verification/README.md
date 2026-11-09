# Moving-apolar constraints on the Fermat sextic

Requirements: SageMath 10.9 and Singular 4.4.1, or compatible versions. All included programs are standalone; SageMath calculations use exact coefficient fields and the central-union Singular check uses the explicitly stated finite field. No files outside this archive, private checkout, unavailable /tmp objects, or network access are needed.

From the extracted archive root, run:

1. `sage calculations/bundles/aoki-sextic-second-order/aoki-sextic-pencil-first-order.sage`
2. `sage calculations/bundles/aoki-sextic-second-order/aoki-sextic-difference-second-order.sage`
3. `sage calculations/bundles/aoki-sextic-second-order/aoki-sextic-sum-order-six.sage`
4. `sage calculations/bundles/aoki-sextic-moving-apolar-s2/aoki-sextic-moving-apolar-s2.sage`
5. `Y0_EXACT=1 sage calculations/bundles/aoki-sextic-y0-hilbert-tangent/aoki-y0-hilbert-tangent.sage`
6. `Singular -q calculations/bundles/aoki-sextic-ci-special-fiber/aoki-sextic-ci-special-fiber.sing`

The moving-apolar-s2 computation was rerun on October 22, 2026; its newly captured output is included. Other output files are historical recorded outputs. The second-order apolar calculation uses streamed coefficient rows; its earlier recorded peak memory was about 450 MiB. The first-order pencil and central Hilbert-tangent computations can take more resources. Run one program at a time.

The finite sum jet reaches order six. It is not an all-orders or convergence certificate. The ramified 261/262 claim of an older working draft is outside this submitted manuscript: its dependency objects have not been recovered and it is not an input to any result here.

Copyright Hongmin Li. Distributed for review and research inspection under the rStation publication terms; no additional software license is asserted.

The Singular program in `aoki-sextic-ci-special-fiber/` constructs only the central union S union S' union V(ell,ell',F), contained in the Fermat sextic. Its degree-18 Hilbert numerator does not certify non-membership of any separate special surface. No such non-membership result is claimed.

The retained historical tangent output uses the label `Hom tangent ambient ... dim 424` for a redundant polynomial-representative kernel. Quotienting by the zero representatives of dimensions 1,8,8,33,33,225 gives the embedded Hilbert tangent dimension 116. The revised source explicitly computes these ideal-piece dimensions from its Groebner basis and reports the quotient. Its 537/538 incidence inconsistency is unchanged. This reporting correction has not been described as a fresh rerun of the full characteristic-zero system.

A quotient-only check on October 23, 2026 computed the six ideal-piece dimensions from the actual Groebner basis, both over GF(109) and over Q(zeta_12,2^(1/3)). Both give [1,8,8,33,33,225] and quotient dimension 116 using the recorded rank 423. Its output is `calculations/bundles/aoki-sextic-y0-hilbert-tangent/hilbert-quotient-check.output.txt`. This small check did not rerun the historical large syzygy or incidence linear systems.
