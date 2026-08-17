# SI3BP rigorous entropy certificate

This package is designed to validate positive topological entropy for the spatial isosceles three-body problem at

\[
\beta=3/100,\qquad e=987/1000,\qquad H=-1,
\]
so that
\[
\alpha=12/97,\qquad \varpi^2=25831/1800,\qquad \beta^2+e^2=0.975069<1.
\]

## What is proved when the program exits with code 0

The code validates:

1. **C1:** a unique symmetric period-four point of the periapsis map and hyperbolicity of `F=P^4`;
2. **C2:** existence and uniqueness of the local energy lift `r=r(p_z,z)` before interval Newton refinement;
3. **C3:** a hyperbolic block, invariant unstable cone, stable-coordinate trapping and unstable exit faces; the graph-transform proposition then gives the true local unstable manifold as a graph;
4. **C4:** opposite signs of the apoapsis `z` coordinate on two tubes containing true unstable-manifold points, regularity of the apoapsis event, and strict negativity of `dp_z/du` along the propagated unstable arc. Reversibility gives a homoclinic point and the Burns--Weiss crossing criterion gives positive topological entropy;
5. **C4+ (optional):** if `dz/du>1e5` is also certified, the homoclinic intersection is transverse and a Smale horseshoe follows directly.

A successful proof run ends with

```text
C1: PASS
C4 (Burns--Weiss crossing): PASS
ALL ENTROPY CERTIFICATES PASSED.
```

`C4+` is stronger but is **not required** for positive entropy.

## Why the energy lift is rigorous

For every local parameter box used by the certificate, the code first checks uniformly

```text
Phi(0.36;p_z,z) > 0,
Phi(0.39;p_z,z) < 0,
H_r < 0 on [0.36,0.39].
```

Thus IVT plus monotonicity gives exactly one energy root for every parameter pair. Only after this existence/uniqueness proof does interval Newton narrow the root enclosure.

## Local build

With CAPD already built:

```bash
make CAPD_CONFIG=/absolute/path/to/CAPD/build/bin/capd-config
./si3bp_capd | tee capd-proof.log
```

The source is written against the public CAPD `IPoincareMap`, `C1Rect2Set` and `C1HORect2Set` interfaces.

## GitHub Actions

The included workflow `.github/workflows/capd-proof.yml` builds a pinned CAPD revision

```text
731079217a9254ea2948d742df2b170895effe7f
```

and then builds and runs the certificate. The log is uploaded as the artifact `si3bp-capd-proof-log`, even if a numerical validation fails. This makes tuning failures distinguishable from mathematical failures.

The long `A o P^24` propagation is subdivided into 32 `u`-cells to reduce wrapping. If C1/C3 validate but C4 does not, the next numerical tuning parameter to increase is this subdivision count; no mathematical statement should be weakened merely to force a pass.
