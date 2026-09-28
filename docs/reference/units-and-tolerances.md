# Units and tolerances

## Units

::: brille._brille.AngleUnit

::: brille._brille.LengthUnit

::: brille._brille.NodeType

## Floating point tolerance

brille decides the zone and mesh geometry with exact predicates (Shewchuk's
[robust predicates](http://www.cs.cmu.edu/~quake/robust.html) and exact
expansion arithmetic), but compares user-provided lattices, symmetry
operations and points with a floating point tolerance, to allow for rounding
in the input and in calculations. Two numbers $a$ and $b$ are equal if
$|a - b| \le t + \epsilon\,|a + b|$, with $t$ an absolute tolerance and
$\epsilon$ machine epsilon.

No one absolute tolerance suits all input: atom positions written as 32-bit
floats, for example, are far less precise than the lattice vectors beside
them. So real-space and reciprocal-space comparisons have separate
tolerances.

The absolute tolerances are $10^{-12}$ Å in real space and $10^{-12}$ Å⁻¹ in
reciprocal space. They can be changed for the whole module,

```python
import brille
brille.real_space_tolerance(1e-10)
brille.reciprocal_space_tolerance(1e-14)
```

or for one call, with an [`ApproxConfig`][brille._brille.ApproxConfig]:

```python
config = brille.ApproxConfig()
config.real_space_tolerance = 1e-10
config.reciprocal_space_tolerance = 1e-14
```

!!! note
    Typical input should not need either. Changing them can cause surprising
    errors elsewhere.

::: brille._brille.ApproxConfig

::: brille._brille.real_space_tolerance

::: brille._brille.reciprocal_space_tolerance
