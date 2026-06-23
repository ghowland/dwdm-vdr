#!/usr/bin/env python3
"""
dwdm-vdr — Exact DWDM fiber propagation simulator using VDR arithmetic.

Demonstrates exact arithmetic preserving structural properties (symmetry)
that float64 breaks, using the split-step Fourier method.

All values are plain VDR rationals. Transcendental constants from
vdr.math.transcendental (precomputed Q335). Arithmetic is exact.

Requires: pip install vdr-math
"""

from __future__ import annotations
import time

try:
    from vdr import VDR
    from vdr.math.transcendental import PI, sin_series, cos_series
except ImportError:
    raise ImportError(
        "dwdm-vdr requires vdr-math. Install with: pip install vdr-math"
    )


def _f(x):
    """VDR to Python float."""
    return float(x.to_fraction())


print("Initializing...", flush=True)

_TWO_PI = VDR(2) * PI
_ALPHA = VDR(46052, 10**9)
_BETA2 = VDR(-20410, 10**30)
_GAMMA = VDR(13, 10000)

_ZERO = VDR(0)
_ONE = VDR(1)
_TWO = VDR(2)
_HALF = VDR(1, 2)
_SIXTH = VDR(1, 6)

_FIELD = VDR(1, 100)  # sqrt(1/10000) exact

TWIDDLE_DEPTH = 6

print("Ready.", flush=True)


def itu_grid(n_channels, spacing_hz):
    anchor = VDR(193_100_000_000_000)
    channels = []
    half_n = n_channels - 1
    for i in range(n_channels):
        offset = 2 * i - half_n
        freq = anchor + VDR(offset) * spacing_hz / _TWO
        channels.append({"freq": freq, "re": _FIELD, "im": _ZERO})
    return channels

def precompute_twiddles(freqs, dz):
    n = len(freqs)
    # True center of symmetry: midpoint between the two middle channels
    if n % 2 == 0:
        center = (freqs[n // 2 - 1] + freqs[n // 2]) / _TWO
    else:
        center = freqs[n // 2]
    twiddles = []
    for i in range(n):
        delta_f = freqs[i] - center
        omega = _TWO_PI * delta_f
        omega_sq = omega * omega
        phase = _BETA2 * _HALF * omega_sq * dz
        c = cos_series(phase, TWIDDLE_DEPTH)
        s = sin_series(phase, TWIDDLE_DEPTH)
        twiddles.append((c, s))
    return twiddles


def apply_linear_step(fields_re, fields_im, atten, twiddles):
    for i in range(len(fields_re)):
        re = fields_re[i] * atten
        im = fields_im[i] * atten
        c, s = twiddles[i]
        fields_re[i] = re * c - im * s
        fields_im[i] = re * s + im * c


def apply_nonlinear_step(fields_re, fields_im, dz, include_xpm):
    n = len(fields_re)
    powers = [fields_re[i] * fields_re[i] + fields_im[i] * fields_im[i]
              for i in range(n)]
    total_pw = _ZERO
    if include_xpm:
        for p in powers:
            total_pw = total_pw + p
    for i in range(n):
        phi = _GAMMA * powers[i] * dz
        if include_xpm:
            phi = phi + _TWO * _GAMMA * (total_pw - powers[i]) * dz
        re = fields_re[i]
        im = fields_im[i]
        fields_re[i] = re - im * phi
        fields_im[i] = im + re * phi


def propagate(fields_re, fields_im, freqs, span_m, steps, include_xpm):
    dz = span_m / VDR(steps)
    half_dz = dz * _HALF
    twiddles = precompute_twiddles(freqs, half_dz)
    x = _ALPHA * half_dz
    atten = _ONE - x + x * x * _HALF - x * x * x * _SIXTH
    for step in range(steps):
        print("    step %d/%d" % (step + 1, steps), flush=True)
        apply_linear_step(fields_re, fields_im, atten, twiddles)
        apply_nonlinear_step(fields_re, fields_im, dz, include_xpm)
        apply_linear_step(fields_re, fields_im, atten, twiddles)


def channel_powers(fields_re, fields_im):
    return [fields_re[i] * fields_re[i] + fields_im[i] * fields_im[i]
            for i in range(len(fields_re))]


def propagate_float(f_re, f_im, freqs_f, span_m_f, steps, include_xpm,
                    alpha_f, beta2_f, gamma_f, two_pi_f):
    n = len(f_re)
    dz = span_m_f / steps
    half_dz = dz / 2.0

    if n % 2 == 0:
        center_f = (freqs_f[n // 2 - 1] + freqs_f[n // 2]) / 2.0
    else:
        center_f = freqs_f[n // 2]

    tw_c = []
    tw_s = []
    for i in range(n):
        df = freqs_f[i] - center_f
        omega = two_pi_f * df
        omega_sq = omega * omega
        phase = beta2_f * 0.5 * omega_sq * half_dz
        c_val, term = 1.0, 1.0
        for k in range(1, TWIDDLE_DEPTH):
            term = -term * phase * phase / ((2 * k - 1) * (2 * k))
            c_val += term
        s_val, term = phase, phase
        for k in range(1, TWIDDLE_DEPTH):
            term = -term * phase * phase / ((2 * k) * (2 * k + 1))
            s_val += term
        tw_c.append(c_val)
        tw_s.append(s_val)
    ax = alpha_f * half_dz
    atten = 1.0 - ax + ax * ax * 0.5 - ax * ax * ax / 6.0
    for step in range(steps):
        for i in range(n):
            f_re[i] *= atten; f_im[i] *= atten
            r, m = f_re[i], f_im[i]
            f_re[i] = r * tw_c[i] - m * tw_s[i]
            f_im[i] = r * tw_s[i] + m * tw_c[i]
        powers = [f_re[i] ** 2 + f_im[i] ** 2 for i in range(n)]
        total_p = sum(powers) if include_xpm else 0.0
        for i in range(n):
            phi = gamma_f * powers[i] * dz
            if include_xpm:
                phi += 2.0 * gamma_f * (total_p - powers[i]) * dz
            r, m = f_re[i], f_im[i]
            f_re[i] = r - m * phi; f_im[i] = m + r * phi
        for i in range(n):
            f_re[i] *= atten; f_im[i] *= atten
            r, m = f_re[i], f_im[i]
            f_re[i] = r * tw_c[i] - m * tw_s[i]
            f_im[i] = r * tw_s[i] + m * tw_c[i]
    return [f_re[i] ** 2 + f_im[i] ** 2 for i in range(n)]


def test_zero_distance():
    print("\n--- Zero Distance ---", flush=True)
    channels = itu_grid(4, VDR(50_000_000_000))
    re = [ch["re"] for ch in channels]
    im = [ch["im"] for ch in channels]
    p1 = channel_powers(re, im)
    p2 = channel_powers(re, im)
    ok = all(p1[i] == p2[i] for i in range(4))
    print("  Exact: %s — %s" % (ok, "PASS" if ok else "FAIL"), flush=True)


def test_single_channel():
    print("\n--- Single Channel 1km ---", flush=True)
    channels = itu_grid(1, VDR(50_000_000_000))
    re = [ch["re"] for ch in channels]
    im = [ch["im"] for ch in channels]
    freqs = [ch["freq"] for ch in channels]
    p_in = _f(channel_powers(re, im)[0])
    propagate(re, im, freqs, VDR(1000), 3, False)
    p_out = _f(channel_powers(re, im)[0])
    print("  Pin: %.6e  Pout: %.6e  Ratio: %.6f — PASS" % (
        p_in, p_out, p_out / p_in), flush=True)


def test_symmetry():
    print("\n--- Symmetry: VDR vs Float ---", flush=True)
    n_ch = 4
    steps = 2
    spacing = VDR(50_000_000_000)
    span_m = VDR(1000)

    channels = itu_grid(n_ch, spacing)
    re = [ch["re"] for ch in channels]
    im = [ch["im"] for ch in channels]
    freqs = [ch["freq"] for ch in channels]
    t0 = time.time()
    propagate(re, im, freqs, span_m, steps, True)
    t_vdr = time.time() - t0
    vdr_pw = channel_powers(re, im)

    channels2 = itu_grid(n_ch, spacing)
    fre = [_f(ch["re"]) for ch in channels2]
    fim = [_f(ch["im"]) for ch in channels2]
    freqs_f = [_f(ch["freq"]) for ch in channels2]
    t0 = time.time()
    flt_pw = propagate_float(fre, fim, freqs_f, _f(span_m), steps, True,
                             _f(_ALPHA), _f(_BETA2), _f(_GAMMA), _f(_TWO_PI))
    t_float = time.time() - t0

    print("  VDR: %.2fs  Float: %.4fs" % (t_vdr, t_float), flush=True)

    vdr_sym = True
    for i in range(n_ch // 2):
        j = n_ch - 1 - i
        match = (vdr_pw[i] == vdr_pw[j])
        if not match:
            vdr_sym = False
        print("  VDR Ch%d vs Ch%d: %s" % (i, j, "EXACT" if match else "BROKEN"), flush=True)

    for i in range(n_ch // 2):
        j = n_ch - 1 - i
        diff = abs(flt_pw[i] - flt_pw[j])
        print("  Float Ch%d vs Ch%d: |diff| = %.2e" % (i, j, diff), flush=True)

    max_drift = 0.0
    for i in range(n_ch):
        p_v = _f(vdr_pw[i])
        p_fl = flt_pw[i]
        drift = abs(p_v - p_fl) / p_v if p_v > 0 else 0.0
        max_drift = max(max_drift, drift)
        print("  Ch%d drift: %.2e" % (i, drift), flush=True)

    print("  VDR symmetric: %s  Max drift: %.2e" % (vdr_sym, max_drift), flush=True)
    print("  %s" % ("PASS" if vdr_sym else "FAIL"), flush=True)


def main():
    print("=" * 50)
    print("dwdm-vdr: Exact DWDM Fiber Propagation")
    print("=" * 50, flush=True)
    test_zero_distance()
    test_single_channel()
    test_symmetry()
    print("\n" + "=" * 50)
    print("COMPLETE")
    print("=" * 50)


if __name__ == "__main__":
    main()
