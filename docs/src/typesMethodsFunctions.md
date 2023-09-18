# Types/Methods/Functions

## Index
```
[NmodController.IonCurrent](#ioncurrent)
```
## Public Interface
```
### NmodController.IonCurrent - Type
"""
    ψ = AMFMcomp(C)
    ψ = AMFMcomp(a, ω, φ)
    ψ = AMFMcomp(a, ω)

Create a 'AMFMcomp' parameterized by a single 'AMFMtriplet'.

# Examples
```@example
 using ISA
 𝐶₀ = AMFMtriplet(t->exp(-t^2),t->2.0,0.0)
 ψ₀ = AMFMcomp(𝐶₀)
```
Another convenient way to create a 'AMFMcomp' is by providing
an *instantenouse amplitude function* `a`,
an *instantaneous frequency function* `ω`, and a *phase reference* `φ`.

Called with two inputs `a, ω`, this is equivalent to `AMFMcomp(a, ω, 0.0)`.

```@example
using ISA
a₀(t) = exp(-t^2)
ω₀(t) = 2.0
φ₀ = 0.0
ψ₀ = AMFMcomp(a₀,ω₀,φ₀)
```
"""
```