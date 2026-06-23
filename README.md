# dwdm-vdr

## Exact DWDM Fiber Propagation Simulator

Exact arithmetic simulation of Dense Wavelength Division Multiplexing fiber optic propagation using the split-step Fourier method. All channel interactions — self-phase modulation and cross-phase modulation — computed with zero accumulated arithmetic error.

Built on [vdr-math](https://pypi.org/project/vdr-math/), an exact arithmetic library where every value is an ordered triple [V, D, R] and the remainder slot catches what conventional systems discard.

## What This Demonstrates

The split-step Fourier method propagates optical fields through fiber in alternating linear (dispersion + attenuation) and nonlinear (Kerr effect) steps. Each step feeds its output into the next. This is the sequential chain structure where floating-point rounding errors accumulate.

dwdm-vdr runs this propagation in exact rational arithmetic alongside an identical float64 implementation. Same algorithm, same constants, same operation order. The only difference is the arithmetic: exact VDR rationals vs IEEE 754 float64.

## Results

### Zero Distance — Exact Identity

Input power equals output power with zero error. Not approximately zero. Exactly zero. The VDR propagation engine introduces no arithmetic artifacts when no physical propagation occurs.

```
Power preserved exactly: True — PASS
```

### Single Channel Attenuation — 1 km SMF-28

One channel at 193.1 THz propagated through 1 km of standard single-mode fiber in 3 steps. Third-order Taylor attenuation model.

```
Pin: 1.000000e-04  Pout: 9.120103e-05  Ratio: 0.912010
```

The output is an exact rational number — the ratio of two integers. Every digit is correct. There is no rounding, no truncation, no accumulated error from the 3 sequential propagation steps.

### Symmetry Preservation — The Key Result

Four channels at 50 GHz spacing, symmetric around the ITU anchor frequency, propagated through 1 km with cross-phase modulation enabled. Channels equidistant from center must produce identical output power by symmetry of the physics.

```
VDR Ch0 vs Ch3: EXACT
VDR Ch1 vs Ch2: EXACT
Float Ch0 vs Ch3: |diff| = 0.00e+00
Float Ch1 vs Ch2: |diff| = 0.00e+00
```

VDR preserves symmetry structurally — not within tolerance, not to machine epsilon, but identically. The concentration at channel i equals the concentration at channel N-1-i because the arithmetic is exact.

### VDR-vs-Float Drift — Machine Epsilon

The per-channel power drift between exact arithmetic and float64 after 2 propagation steps:

```
Ch0 drift: 5.94e-16
Ch1 drift: 8.92e-16
Ch2 drift: 8.92e-16
Ch3 drift: 5.94e-16
```

The drift is symmetric (Ch0 matches Ch3, Ch1 matches Ch2) and sits at machine epsilon (~10⁻¹⁶). This is the expected result for a short propagation chain with a linear-dominant system. The drift itself is real float64 arithmetic error — not an algorithm mismatch, not a modeling difference. Both paths run the identical algorithm. The drift is purely from IEEE 754 rounding.

## Why This Matters for DWDM Engineers

### Simulator Validation

Every production DWDM simulator (VPItransmissionMaker, OptSim, GNPy) uses float64 arithmetic. Each operation introduces a rounding error of approximately 10⁻¹⁶. For a single operation this is negligible. For a fiber propagation simulation — where hundreds of sequential nonlinear operations feed their outputs into the next step — the errors accumulate.

dwdm-vdr provides exact ground truth for any propagation scenario. Run your scenario in dwdm-vdr and in your production simulator. The difference is your simulator's arithmetic noise floor. If you are making engineering decisions where the signal is close to that noise floor, you need to know where it is.

### Channel Spacing Decisions

The industry is pushing toward 12.5 GHz and 25 GHz channel spacings. At these spacings, the nonlinear interactions between channels are stronger and the engineering margins are thinner. The question a DWDM engineer needs to answer is: "Will these two channels interfere destructively at this spacing?"

If the arithmetic noise in your simulator is comparable to the interference signal you are trying to measure, the answer is unreliable. dwdm-vdr lets you establish the exact answer for validation geometries and calibrate your production tools against it.

### Symmetry as a Diagnostic

A symmetric channel configuration must produce symmetric output power. This is a physical invariant — it does not depend on fiber parameters, power levels, or propagation distance. If your simulator breaks this symmetry, the asymmetric component is arithmetic noise, not physics.

dwdm-vdr preserves this symmetry exactly. You can use symmetry violation as a diagnostic: set up a symmetric configuration in your production simulator, measure the asymmetry, and that tells you the magnitude of arithmetic contamination in your specific scenario.

### Where the Drift Grows

The 2-step, 1 km result shown here produces drift at machine epsilon because the system is linear-dominant and the chain is short. The drift grows with:

**More propagation steps.** A transoceanic link with 12,500 steps accumulates more rounding than a datacenter link with 5 steps.

**Stronger nonlinearity.** Higher launch powers, tighter channel spacings, and more channels increase the Kerr phase per step. The small-angle rotation in the nonlinear step does not exactly preserve field magnitude — it grows power by a factor of (1 + φ²) per step. In exact arithmetic this is tracked precisely. In float64 it accumulates asymmetrically.

**Multi-span amplification.** EDFA gain applied after each span multiplies the field (and its arithmetic error) by the gain factor. Over many spans, gain-error multiplication compounds the drift beyond what single-span propagation produces.

## How It Works

### Split-Step Method

The symmetric split-step method applies half-linear, full-nonlinear, half-linear per step:

1. **Half linear step:** Attenuation via third-order Taylor approximation of exp(-αz/2), then chromatic dispersion via precomputed twiddle factors (cos and sin from Taylor series at depth 6).

2. **Full nonlinear step:** Kerr phase rotation via first-order small-angle approximation (re' = re - im·φ, im' = im + re·φ) where φ includes self-phase modulation and cross-phase modulation contributions.

3. **Half linear step:** Same as step 1.

### Exact Arithmetic

All values are plain VDR rationals. No basis projection, no rounding, no truncation. The denominator grows naturally through exact rational operations. Every intermediate value at every channel at every step is an exact ratio of two integers.

Transcendental constants (π for frequency-to-angular-frequency conversion) use the Q335 precomputed value from vdr-math — a 100-digit rational on the 2³³⁵ denominator grid.

The Taylor series for sin and cos (twiddle factors) produce exact rational partial sums. At depth 6, the truncation error for the small dispersion phases in this simulation is below 10⁻³⁰ — far smaller than any physical effect.

### Float Mirror

An identical algorithm runs in float64 in parallel. Same Taylor depths, same small-angle approximation, same operation order. The per-channel drift between VDR and float measures exactly how much error float64 accumulates over the propagation chain.

### Physical Parameters

| Parameter | Value | Representation |
|---|---|---|
| Attenuation | 0.2 dB/km | VDR(46052, 10⁹) per meter |
| Dispersion β₂ | -20.41 ps²/km | VDR(-20410, 10³⁰) s²/m |
| Nonlinear γ | 1.3 /W/km | VDR(13, 10000) /W/m |
| Launch power | -10 dBm | VDR(1, 10000) W |
| Field amplitude | exact √(power) | VDR(1, 100) W^½ |
| Channel spacing | 50 GHz | VDR(50000000000) Hz |
| Grid points | 4 channels | symmetric around ITU anchor |

## Installation

```
pip install vdr-math
python dwdm-vdr.py
```

Requires Python 3.8+ and vdr-math. No other dependencies.

## Extending This Work

### More Channels and Steps

Increase `n_ch` and `steps` in the symmetry test to push float64 harder. The VDR computation time grows with step count and channel count (the rational denominators grow through the chain). The float mirror runs in microseconds regardless. The gap between VDR symmetric (exact) and float symmetric (broken at sufficient chain length) is the result.

### Four-Wave Mixing

Add FWM to the nonlinear step for degenerate channel triplets. The phase-matching efficiency uses a Lorentzian approximation (exact rational). FWM contributions are additive field perturbations. This is straightforward to add and exercises the tightest-spacing regime where arithmetic error matters most.

### Concentration-Dependent Effects

Replace the constant γ with a power-dependent nonlinear coefficient. This creates the nonlinear feedback loop where each step's output modifies the next step's coefficients — the same structure as concentration-dependent diffusion in TCAD. This is where float drift should grow beyond machine epsilon.

### Production Validation Workflow

1. Define your channel plan (frequencies, spacings, powers)
2. Define your link (fiber type, span lengths, amplifier gains)
3. Run in dwdm-vdr to get exact output powers and phases
4. Run the identical scenario in your production simulator
5. Compare per-channel: the difference is your arithmetic noise floor
6. If the noise floor is below your engineering margin, your simulator is trustworthy for that scenario
7. If the noise floor approaches your margin, you need higher precision or shorter step sizes

## Dependencies

- [vdr-math](https://pypi.org/project/vdr-math/) — exact arithmetic library

## License

MIT
