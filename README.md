# hxrgutils

Python utilities for raw data from Teledyne HxRG infrared detectors (H2RG / H4RG),
as used by the NIRPS, SPIRou and SPIP spectrographs: reading individual reads,
reference-pixel correction, up-the-ramp (slope) fitting and noise diagnostics.

```bash
pip install -e .
```

Requires Python >= 3.10 and `numpy`, `scipy`, `astropy`, `matplotlib`, `seaborn`, `tqdm`
and `natsort`; `pip` installs them from `setup.py`.

## Contents

1. [Concepts](#concepts)
2. [Reading an individual read (`hxread`)](#1-reading-an-individual-read-hxread)
3. [CDS with reference-pixel correction (`hxread` only)](#2-cds-with-reference-pixel-correction-hxread-only)
4. [Building a ramp file (`hxramp`)](#3-building-a-ramp-file-hxramp-from-hxread-reads)
5. [Reference-pixel options in `hxramp`](#4-reference-pixel-options-in-hxramp)
6. [Reading a ramp file (`hxramp`): NIRPS vs SPIRou/SPIP](#5-reading-a-ramp-file-hxramp-nirps-vs-spirouspip)

## Concepts

- Arrays are square, `1024 * x` pixels per side: `x=2` for an H2RG (2048), `x=4` for an
  H4RG (4096). The geometry is **always inferred from the data you load**; the `h2rg`
  keyword of `hxramp` is only a hint that the data overrides.
- The detector is read through 32 amplifiers (vertical ribbons of width `w/32`).
- The outer 4 pixels on each edge are reference pixels: the left/right ones track fast
  1/f noise, the top/bottom ones track per-readout drift.
- A full frame takes `(1024*x)**2 * 1e-5 / 32` s (`hxramp.inttime`); this is the time
  step of the ramp fit.
- Default saturation threshold: 45000 ADU.

---

## 1. Reading an individual read (`hxread`)

`hxread` is built empty and then **called** to load data: `__call__` is the loader.

```python
from hxrgutils.type.read import hxread

r = hxread()
r('NIRPS_R01_R01.fits')          # from a FITS file (data in the primary HDU)
r(array)                         # from a 2D numpy array (header optional)
r(array, header=some_header)     # ... with a header attached
```

| Source | Header (`r.H`) | `r.fname` |
| --- | --- | --- |
| FITS path | header of the primary HDU (or `header=` if given) | basename of the file |
| 2D array | `header=` if given, otherwise an **empty** `fits.Header` | `"tmp.fits"` |

Notes:

- A path without `/` is resolved relative to the current directory.
- The data is converted to `float`. A primary HDU without data, or an array that is not
  2D, raises `ValueError` (use `make_ramp()` / `hxramp` for cubes).
- `r.H` is always a real `fits.Header`, never `None`, whichever source was used.
- After loading, `r.im` (alias `r.data`) holds the image and `r.shape`, `r.h`, `r.w`, `r.x`
  are re-derived from it. A size that is not a multiple of 1024 only prints a warning.

Useful members once a read is loaded:

```python
r.im                    # the 2D image (float)
r.H                     # header
r.stats()               # sigma-clipped (mean, median, std); NaN/inf are dropped
r.stats('median')       # a single value: 'mean', 'median' or 'std'
r.get_amp(5)            # amplifier 5 (1..32) as an hxrg_amp
r.unfold()              # flip every other amplifier into readout-time order (in place)
r.copy()                # deep copy (image and header)
r.writeto('out.fits')   # write image + header
r.show()                # zscale image;  r.dist() plots the pixel histogram
r.fft()                 # 1/f / pattern-noise amplitude spectrum (Hz, amplitude)
```

Arithmetic works between reads and with scalars. `a - b`, `a + b`, `a * b`, `a / b` return
a **new** `hxread` (header copied from `a`); `+=`, `-=`, `*=`, `/=` work in place.

---

## 2. CDS with reference-pixel correction (`hxread` only)

A correlated double sampling (CDS) image is the difference of two reads. No dedicated
class is needed: subtract two `hxread` objects, then apply the reference-pixel
correction to the result.

```python
from hxrgutils.type.read import hxread

r1, r2 = hxread(), hxread()
r1('read_01.fits')
r2('read_02.fits')

cds = r2 - r1          # new hxread; r1 and r2 are untouched
cds.refpx()            # full correction, in place on cds.im
cds.writeto('cds.fits')
```

`refpx()` accepts one of two keyword flags:

| Call | What is corrected | Underlying method |
| --- | --- | --- |
| `cds.refpx()` | side reference pixels **+** top and bottom (default) | `refpxcorr()` |
| `cds.refpx(toponly=True)` | side reference pixels + **top** only (no bottom) | `refpxcorrtop()` |
| `cds.refpx(sideonly=True)` | side (left/right) reference pixels only | `refpxcorrside()` |

`refpx()` also appends a `refpixcorr performed [...]` comment to `cds.H`.

For the finer options call the underlying method directly. They accept `oddeven` and
`medfiltersize`, which `refpx()` does not forward:

```python
cds.refpxcorr(oddeven=True)                    # odd/even columns corrected separately
cds.refpxcorrtop(oddeven=True)
cds.refpxcorrside(medfiltersize=21)            # width of the running median on the side pixels (default 15)
```

- `oddeven=True`: odd and even columns inside each amplifier are measured and corrected
  independently (useful for arrays where they behave differently).
- `medfiltersize`: width of the running median applied to the side reference pixels
  before they are subtracted (default 15).

What each step does: the side pixels (columns 0-3 and N-4..N-1) are median-filtered and
averaged into one value per row, which is subtracted from the whole row (removes fast
1/f noise common to all amplifiers). Then, per amplifier, the median of the top and bottom
reference rows defines a linear drift between the bottom and the top, which is subtracted.
For `toponly`, only the top reference rows are used, as a constant offset per amplifier.

> **Caveats**
>
> - The correction works **in place**: it modifies `cds.im` and returns it. A copy is
>   made by `r2 - r1`, so `r1`/`r2` stay intact, but if you load a float64 array with
>   `r(array)` the read may share memory with your array. Pass `np.copy(array)` if you
>   need the original.
> - The correction is **not NaN-safe** on the side pixels: one NaN in a side reference
>   column wipes out whole rows after the median filter. Clean NaNs first.
> - The same methods can be called on a single read (`r.refpx()`), not only on a CDS.

---

## 3. Building a ramp file (`hxramp`) from `hxread` reads

`hxramp` fits a line to each pixel's signal vs. time. Reads are accumulated one at a
time with `<<` (running sums plus a saturation mask), then `fit()` solves the fit.

### Streaming, step by step

```python
from hxrgutils.type.read import hxread
from hxrgutils.type.ramp import hxramp

R = hxramp(h2rg=True, saturation=45000)   # h2rg=True -> 2048x2048, else 4096x4096 (data wins)
r = hxread()
for f in files:                           # files must be in time order
    r(f)
    R << r                                # an hxread; a raw ndarray or a path also works
R.fit()                                   # fills R.a (slope, ADU/s) and R.b (intercept)
R.writeto('ramp.fits')
```

What `fit()` produces:

- `R.a`: slope in ADU/s. `R.b`: intercept (bias added back). `R.n`: number of
  non-saturated reads used per pixel.
- Pixels with fewer than 2 valid reads keep slope 0, not NaN. A NaN pixel never passes the
  saturation mask, so it is dropped **from that read onward**: a NaN in read 3 leaves
  `n == 2` and a slope fitted on the first two reads; only a pixel that is NaN from the
  start ends with `n == 0`. Check `R.n` if your data contains NaNs.
- Without a bias file, the **first read is used as the bias** (self-bias) and is
  subtracted from every read.

### In one call: `make_ramp`

```python
from hxrgutils.type.ramp import make_ramp

R = make_ramp(files,                       # list of paths/arrays, a cube file name, or a 3D array
              errors=True,                 # also compute the slope error (see below)
              nonlin='nl.fits',            # optional non-linearity file
              bias='bias.fits',            # optional super-bias
              h2rg=True)                   # remaining kwargs go to hxramp()
R.writeto('ramp.fits')
```

`make_ramp` creates the `hxramp`, streams every read, calls `fit()`, and optionally
`create_error()`. Pass the list **already in time order**: streaming uses the order
given (only the error step re-sorts file lists with `natsort`).

### Loading reads without a loop

```python
R = hxramp(h2rg=True)
R(cube_array)        # 3D array, or list of 2D arrays, or list of paths: streamed through <<
R('cube.fits')       # FITS file holding a cube of reads (3D primary HDU)
R.fit()              # you still have to call fit()
```

### Optional calibration files

| Option | Method / kwarg | File expected |
| --- | --- | --- |
| Super-bias | `R.upload_bias('bias.fits')` or `make_ramp(bias=...)` | 2D image, same shape as the reads |
| Non-linearity | `R.upload_nonlin('nl.fits')` or `make_ramp(nonlin=...)` | ext 0 = C3 coefficient, ext 1 = C2 coefficient (2D each); applied as `x + C3*x**3 + C2*x**2` on the bias-subtracted read |

Upload them **before** streaming reads. Called directly, `upload_bias()` /
`upload_nonlin()` print a message for a missing file and fall back to self-bias / no
non-linearity correction. **`make_ramp(bias=..., nonlin=...)` is silent**: a path that does
not exist is skipped without any message, so a typo gives a self-biased, uncorrected ramp.
Check the path, or call the `upload_*` methods yourself.

#### Expected shapes

Both files must match the detector format of the reads, as `(height, width)` (row-major):
**2048 x 2048** for an H2RG, **4096 x 4096** for an H4RG.

| File | HDUs | Content | Shape |
| --- | --- | --- | --- |
| Super-bias | primary HDU (`ext 0`) | 2D image of the bias level, in ADU | `(N, N)`, same as one read |
| Non-linearity | `ext 0` | C3 coefficient (2D) | `(N, N)` |
| | `ext 1` | C2 coefficient (2D) | `(N, N)`, same as `ext 0` |

- The bias must be exactly 2D; a cube or 1D array raises `ValueError`.
- Uploading the bias **before any read** also sets the ramp geometry (so `h2rg` is
  not needed); the reads that follow must then match it.
- If the ramp geometry is already known (a read was streamed, or a bias uploaded
  first), a bias or NL file of a different shape raises `ValueError`.
- If the NL file is uploaded first and the reads turn out to have a different shape,
  the first read prints a warning and the non-linearity correction is **disabled**
  (no error), so check the output messages.
- A `UID` keyword in the primary header of either file is copied into the ramp header
  (`BIASUID` / `NLUID`).

### Slope error and output layout

`R.create_error(files_or_cube)` computes the formal error on the slope from the same
reads. It needs at least 3 reads. After it has run, `writeto()` writes the 5-HDU layout
(see [section 5](#5-reading-a-ramp-file-hxramp-nirps-vs-spirouspip)); otherwise it
writes the 4-HDU layout.

---

## 4. Reference-pixel options in `hxramp`

All options are keyword arguments of `hxramp(...)` (and therefore of `make_ramp(...)`).
`fit()` prints which combination was applied.

| Keyword | Default | Effect |
| --- | --- | --- |
| `norefpxcorr` | `False` | Skip the reference-pixel correction entirely. |
| `refpxeachread` | `False` | Correct each read (after bias and non-linearity) as it is streamed, instead of correcting the slope once at the end. |
| `toponly` | `False` | Use the side + **top** reference pixels (`refpxcorrtop`). |
| `notopnobottom` | `False` | Use the **side** reference pixels only (`refpxcorrside`). |
| `oddeven` | `False` | Odd and even columns are corrected separately. |
| `saturation` | `45000` | Reads above this value (ADU) are excluded from the fit for that pixel and all later reads. |
| `h2rg` | `False` | Hint for a 2048x2048 array (the loaded data overrides it). |

Which correction is used (priority: `toponly`, then `notopnobottom`, then default):

| `toponly` | `notopnobottom` | Pixels used |
| --- | --- | --- |
| `False` | `False` | side + top + bottom (full `refpxcorr`) |
| `True` | any | side + top |
| `False` | `True` | side only |

When it is applied:

- **Default** (`refpxeachread=False`): once, on the **slope** `R.a` at the end of `fit()`.
  The intercept `R.b` is not corrected.
- **`refpxeachread=True`**: on every `(read - bias)` as it is streamed, so the
  fit sees corrected reads.
- **`norefpxcorr=True`**: never. (With `refpxeachread=True` as well, `norefpxcorr`
  wins at fit time, but the per-read correction in `<<` is still applied; do not combine them.)

Examples:

```python
R = hxramp(h2rg=True)                                  # full correction on the slope
R = hxramp(h2rg=True, oddeven=True)                    # ... odd/even columns separately
R = hxramp(h2rg=True, refpxeachread=True)              # correct every read
R = hxramp(h2rg=True, toponly=True)                    # side + top only
R = hxramp(h2rg=True, notopnobottom=True)              # side only
R = hxramp(h2rg=True, norefpxcorr=True)                # raw slope, no correction
```

The correction is applied through `hxramp.applyrefpxcorr(array)`, which you can also call
yourself on any 2D array with the options currently set on the ramp.

> **Known limitations.** `applyrefpxcorr()` does not look at `norefpxcorr`, and
> `create_error()` always calls it, so with `norefpxcorr=True` the slope is uncorrected
> but the slope error is still computed on reference-pixel-corrected reads. Also,
> `fit()` prints "Ref. px. corrected on the slope" and "All ref. px. used" even when
> `norefpxcorr=True`; the printed messages are not reliable in that case.

---

## 5. Reading a ramp file (`hxramp`): NIRPS vs SPIRou/SPIP

An already-written ramp is loaded by calling the object with its path. `hxramp` tells
the instrument layout apart by the **number of HDUs**:

| Instrument | HDUs | Layout |
| --- | --- | --- |
| **NIRPS** | 4 | `0` primary (header), `1` slope, `2` intercept, `3` n |
| **SPIRou / SPIP** | 5 | `0` primary (header), `1` slope, `2` intercept, `3` errors (slope error), `4` n |

```python
from hxrgutils.type.ramp import hxramp

R = hxramp()
R('ramp_nirps.fits')        # 4 HDUs -> NIRPS
R('ramp_spirou.fits')       # 5 HDUs -> SPIRou/SPIP (prints "SPIRou/SPIP")

R.a                         # slope (ADU/s)
R.b                         # intercept
R.n                         # number of valid reads per pixel (int16)
R.errslope                  # SPIRou/SPIP only; NaN-filled for NIRPS files
R.H                         # primary header
R.total_integration         # see below
```

What differs:

- **NIRPS (4 HDUs)**: no error extension. `R.err_calulated` stays `False`, and `R.errslope`
  is a NaN array. If you call `R.create_error(reads)` and then `R.writeto(...)`, the
  output becomes the 5-HDU layout.
- **SPIRou/SPIP (5 HDUs)**: the slope error is read from HDU 3 and `R.err_calulated` is
  `True`. `writeto()` re-creates the same 5-HDU layout.
- Any other HDU count raises `ValueError("Shape of ramp file not recongnized ...")`.
- A file with a single HDU, or a 3D primary data array, is **not** a ramp: it is treated
  as a cube of reads and streamed (you then call `fit()`), see
  [section 3](#loading-reads-without-a-loop).
- `total_integration`: for NIRPS files it is taken from the header keyword
  `HIERARCH ESO DET SEQ1 DIT` if present; otherwise it is `inttime * max(n)`.
- The geometry (H2RG vs H4RG) is taken from the file, regardless of the `h2rg` keyword.
- A ramp loaded from a fitted file is read-only for streaming: `R << read` raises
  `RuntimeError`. Create a new `hxramp()` to stream reads.

Other things you can do with a loaded ramp:

```python
a, b = R.get_coeff()      # slope and intercept as hxread objects (header copied)
R.unfold()                # flip odd amplifiers of a, b, n into readout-time order
R.preprocess()            # quick preprocessing of the slope
R.show()                  # plot slope and intercept
R.writeto('copy.fits')    # 4 or 5 HDUs depending on whether errors are present
```

`R.create_error(reads)` also works on a loaded ramp, to add the slope error to a NIRPS
file from its reads.

---

## Notes

- NaN/inf pixels are dropped by `stats()`. `fft()` skips them per sample and interpolates
  samples where all 32 amplifiers are NaN.
- Plotting and optional dependencies (`matplotlib`, `seaborn`, `tqdm`, `scipy.fft`) are
  imported lazily inside the methods that use them.
- `hxrgutils/fits2ramp/fits2ramp.py` is the legacy Python 2.7 standalone script and the
  `*.bkp` files are dead snapshots; they are kept as reference only.
