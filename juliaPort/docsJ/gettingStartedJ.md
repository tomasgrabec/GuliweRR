# Getting started with GuliweRR in Julia

GuliweRR exists in MATLAB and in Julia. Both give the same results (to round-off); the Julia
version is several times faster, so use it for calculations. This page gets you from nothing
to a first plot. After that, continue with `quickStartJ.md` and the tutorials.

## 1. Install Julia

Use **juliaup**, Julia's version manager (Julia 1.13 or newer).

- macOS / Linux, in a terminal:
  ```
  curl -fsSL https://install.julialang.org | sh
  ```
- Windows, in a terminal:
  ```
  winget install --name Julia --id 9NJNWW8PVKMN -e -s msstore
  ```

Open a new terminal and check with `julia --version`.

## 2. Make a work environment (once)

A Julia *environment* is a folder that lists the packages your scripts use — the counterpart
of MATLAB's path. Make one for your GuliweRR work, add GuliweRR and a plotting package:

```
mkdir ~/guliwerrWork
cd ~/guliwerrWork
julia --project=.
```

In Julia:

```julia
using Pkg
Pkg.develop(path = "/path/to/GuliweRR/juliaPort/guliwerrJ")   # your copy of the repository
Pkg.add("GLMakie")       # plots in windows; use "CairoMakie" for files / no display
```

This downloads and compiles the dependencies — a few minutes, once.

## 3. Start Julia and run a first calculation

Always start Julia in your environment **with threads** — Julia's threads are fixed at start-up:

```
cd ~/guliwerrWork
julia --project=. -t auto
```

```julia
using GuliwerrJ, GLMakie

ec  = ElastConsts([170.0, 122.0, 74.0], 2)        # Cu: C11, C12, C44 (GPa), cubic
mat = Material(elConsts = ec, rho = 8.96)          # density g/cm^3; default cut (001)
cp  = CalcParams(relDepth = 20, legOrder = 50)
str = Structure(materIn = mat, calcParamsIn = cp, name = "cu_001", nMaxModes = 500)

sweep = calc(angSweepHS(structure = str, angles = collect(0.0:1.0:90.0),
                        calcMode = :full, retentionLevel = :full))
plotGreenMap(sweep)                                # |G13_eff| over angle and velocity
display(current_figure())                          # show it (needed in scripts)
```

The **first** run of a session takes ~20 s: Julia compiles the code it needs. Running it
again takes under a second.

To run a whole script: `julia --project=. -t auto myScript.jl`, or `include("myScript.jl")`
at the Julia prompt. With GLMakie, use the `include` route: a window opened by a script run
from the terminal closes as soon as the script ends. (With CairoMakie, end the script with
`save("fig.png", current_figure())` instead.)

## 4. An editor

**VS Code** with the **Julia** extension is the supported Julia editor and the closest to the
MATLAB desktop (variables pane, plot pane, debugger):

1. Install VS Code and the extension "Julia" (by julialang).
2. Settings → search `julia.NumThreads` → set to `auto`.
3. Open your `guliwerrWork` folder; click "Julia env" in the bottom bar and choose it.
4. In a `.jl` file, `Shift+Enter` runs the current line or block, like a MATLAB section.

Other options: the plain terminal prompt with any text editor; **Jupyter** notebooks
(`Pkg.add("IJulia")`); **Pluto** notebooks, which update dependent cells automatically
(`Pkg.add("Pluto"); using Pluto; Pluto.run()`).

## 5. MATLAB → Julia in brief

| MATLAB | Julia |
|---|---|
| `run('guliwerr/init.m')` | `using GuliwerrJ` (in your environment) |
| `s = Calc.angSweepHS(str, angles, 'full', 'full', []); s = s.calc();` | `s = calc(angSweepHS(structure = str, angles = angles, calcMode = :full, retentionLevel = :full))` |
| `'full'`, `'ritzOnly'` (strings) | `:full`, `:ritzOnly` (symbols) |
| `0:1:90` | `collect(0.0:1.0:90.0)` |
| `a .* b`, `sin(x)` on arrays | `a .* b`, `sin.(x)` — the dot applies any function element-wise |
| `x(1)`, `x(end)` | `x[1]`, `x[end]` — square brackets, also 1-based |
| `% comment` | `# comment` |
| `[]` for "not given" | `nothing` |

Constructors take **named arguments** (`Material(elConsts = ec, rho = 8.96)`), not positions —
except `ElastConsts`, which is positional (`ElastConsts([170.0, 122.0, 74.0], 2)`).

## Next

- `quickStartJ.md` — the full forward example, classification methods, viewers, inversion.
- `juliaPort/tutorialsJ/` — the tutorials, file for file the same as the MATLAB ones.
- `functionReferenceJ.md` — every function; `pitfallsJ.md` — known traps;
  `matlabVsJuliaJ.md` — every difference from the MATLAB version.
- Help at the prompt: `?` followed by a name, e.g. `?ElastConsts`.
