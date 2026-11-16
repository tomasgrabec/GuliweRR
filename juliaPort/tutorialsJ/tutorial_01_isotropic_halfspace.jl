## Tutorial: isotropic half-space -- minimal SAW velocity calculation (Julia port)
#
# Same physics as examples/tutorials/tutorial_01_isotropic_halfspace.m: shortest path from
# "a material" to a SAW (surface acoustic wave) velocity. Verifies the textbook approximation
# v_SAW ~ 0.92*v_T (transverse bulk velocity) for an isotropic solid, which has no angle
# dependence, so a single propagation angle suffices.
#
# This is tutorial 1 of the GuliweRR-Julia ("GuliwerrJ") series -- it teaches the absolute
# minimum Julia needed to run one calculation. See juliaPort/docsJ/quickStartJ.md for the
# general call sequence, juliaPort/docsJ/pitfallsJ.md for the gotchas specific to the Julia
# port, and docs/pitfalls.md for the legO/relDepth rule and the other gotchas that carry over
# unchanged from MATLAB.
#
# --- How to run this script -------------------------------------------------------------
# Two ways, matching MATLAB's "press Run" habit. Both need the package's dependencies
# downloaded once (see juliaPort/docsJ/gettingStartedJ.md, section 2, for the recommended
# own environment with GuliwerrJ and a plotting package; this tutorial alone can also run in
# the package's own environment after
# `julia --project=juliaPort/guliwerrJ -e "using Pkg; Pkg.instantiate()"`).
#
#   1. From a terminal (like running a MATLAB script with `matlab -batch "run('...')"`):
#        julia --project=juliaPort/guliwerrJ -t auto juliaPort/tutorialsJ/tutorial_01_isotropic_halfspace.jl
#      `--project=...` tells Julia which package ENVIRONMENT to use -- the closest analogue
#      to MATLAB's path: it says "use exactly the package versions recorded for this
#      environment", not whatever happens to be globally installed. The later tutorials plot
#      with GLMakie (or CairoMakie), which the package's own environment does not contain
#      (GuliwerrJ treats plotting as optional); run those in your own environment from
#      gettingStartedJ.md, with `--project=<that environment>` instead.
#      `-t auto` starts Julia's thread-based parallelism -- the analogue of MATLAB's `parfor`
#      -- which is OFF by default (a script without `-t auto`/`-t N` runs single-threaded:
#      correct, just serial), unlike MATLAB, which always has `parfor` available once a
#      parallel pool is started. This tutorial's own sweep is one angle, so it makes no
#      difference here; it does for the LARGER sweeps in later tutorials.
#
#   2. From an already-running Julia session (the closer analogue to MATLAB's "Run" button,
#      keeping variables alive afterwards so you can inspect them):
#        julia> include("/path/to/GuliweRR/juliaPort/tutorialsJ/tutorial_01_isotropic_halfspace.jl")
#      This tutorial has no plotting, so nothing else is needed for this particular script.
#
# ------------------------------------------------------------------------------------------
# Julia concepts introduced in this tutorial:
#   - `using` to load a package (like MATLAB's path, but explicit and per-package)
#   - structs (Material/ElastConsts/CalcParams/Structure are all Julia "structs" -- fixed-shape
#     bundles of named fields, MATLAB's classdef value classes translated one-for-one)
#   - keyword-only constructors, e.g. `Material(; elConsts = ec, rho = 2.70)`
#   - immutability: these input structs (`ElastConsts`, `Material`, `CalcParams`, `Structure`)
#     can never be edited in place; every "setter" returns a NEW object instead. The `Calc*`
#     results are different: they are mutable, like MATLAB handle objects
#   - multiple dispatch replacing classdef methods: `calc(sweep)`, not `sweep.calc()`
#   - Symbols (`:ritzOnly`, `:saw_only`) as Julia's typed-option idiom, replacing MATLAB's
#     quoted option strings
#   - 1-based indexing is UNCHANGED from MATLAB -- `sweep.sawVelMap[1][1]` reads exactly like
#     MATLAB's `sweep.sawVelMap{1}(1)` would, just without the curly/round-bracket distinction
#   - Julia scripts never auto-print a statement's value (unlike a bare MATLAB line without a
#     trailing `;`), so there is nothing to suppress here -- explained inline below
#   - `Printf.@printf`, the direct analogue of MATLAB's `fprintf`
#
# Run guliwerr/init.m once per MATLAB session before the MATLAB tutorials -- GuliwerrJ needs
# no equivalent step: `using GuliwerrJ` below is the whole of "path setup" (see the module
# docstring in src/GuliwerrJ.jl for why: Julia has no `addpath`/`genpath` at all, so there is
# no MATLAB-style path hazard to avoid in the first place -- see docs/pitfalls.md OTH2's MATLAB
# warning, which simply does not apply here).

using GuliwerrJ
# `using GuliwerrJ` loads the GuliwerrJ package and brings every name it exports (Material,
# ElastConsts, CalcParams, Structure, angSweepHS, calc, ...) into scope, unqualified. This is
# the Julia analogue of MATLAB's `run('guliwerr/init.m')` path setup -- except it is a genuine
# language statement (not a script side effect), and it only ever affects the names visible in
# THIS file, never a global path setting that leaks into other sessions.
using Printf
# The standard-library module providing `@printf`, used at the bottom of this script. Standard
# library modules ship with every Julia install but still need their own `using` line -- they
# are not automatically in scope, unlike MATLAB's built-in functions.

## Material -- isotropic Al. ElastConsts' isotropic class takes [c11, c44]
# only (c12 = c11 - 2*c44 is implied) -- see SymmClassLibrary() symmNo=1 (same convention as
# the MATLAB original's SymmClassLibrary.m).
# The values are the Hill average of the single-crystal Al constants from the materials
# library (`ListMaterialLibrary(; print = true)`), the usual isotropic estimate. Both calls
# below return two values (`_` discards the library row; the Hill average is the pair C11, C44).
matAlCrystal, _ = MaterialFromLibrary("Al")
alC11, alC44 = HillAverages(matAlCrystal.elConsts.voigt)   # GPa
alConsts = [alC11, alC44]
alSymm   = 1       # isotropic
alRho    = matAlCrystal.rho   # g/cm^3

ec  = ElastConsts(alConsts, alSymm)
# ElastConsts takes its two required arguments POSITIONALLY (consts, then symmNo) -- there is
# no ambiguity to resolve with keywords here, so the port keeps this one positional, exactly
# mirroring MATLAB's own `ElastConsts(alConsts, alSymm)` call syntax.
mat = Material(; elConsts = ec, rho = alRho)
# Material's constructor, by contrast, is KEYWORD-ONLY: every argument after the semicolon `;`
# must be named. `Material(ec, alRho)` (MATLAB-style positional) does NOT work here -- it is
# not merely unconventional, it will not compile. This is deliberate, not a style preference:
# a value class in this port with several optional constructor arguments is translated to a
# Julia keyword constructor throughout so
# that a call reads unambiguously regardless of which arguments are supplied, exactly the way
# MATLAB's own `Material.m` constructor accepts `elConsts`/`rho`/`dirMat`/`name` in any
# combination via `nargin`-guarded `if` checks -- Julia expresses that same "any subset, in any
# order" flexibility with keywords instead. `dirMat`/`name` are both omitted here and fall back
# to Material's own MATLAB-matching defaults (dirMat = [1 0 0; 0 0 1; 0 -1 0], name = "mater").

## Discretisation -- CalcParams(; relDepth, legOrder), "easy" tier (isotropic/principal plane)
relDepth = 20   # truncation depth of the computational domain, in wavelengths
legOrder = 50   # maximal order of Legendre-polynomial basis functions spanning
                # that depth -- should be >= 2*relDepth or the depth profile is
                # under-resolved (docs/pitfalls.md MTH1)
cp = CalcParams(; relDepth = relDepth, legOrder = legOrder)
# Note the field is called `legOrder` here (Julia), not MATLAB's `legO` positional-argument
# name -- the port keeps MATLAB's PROPERTY names (`legOrder` is `CalcParams.m`'s own property),
# not every local variable name every MATLAB tutorial happened to use for it.

## Structure -- single half-space layer
strName   = "al_isotropic_halfspace"   # Julia strings use double quotes; single quotes are a
                                          # single CHARACTER in Julia, e.g. 'a', not a string
nMaxModes = 500   # max number of Ritz modes the solver will compute/retain
str = Structure(; materIn = mat, calcParamsIn = cp, name = strName, nMaxModes = nMaxModes)
# Structure's constructor takes only the FIRST layer this way; a film-on-substrate stack (later
# tutorials) is built by then calling `AddLayer(str, filmMaterial, filmCalcParams)` -- a plain
# function taking the structure as its first argument and returning a NEW Structure, replacing
# MATLAB's `str = str.AddLayer(matFilm, cpFilm)` method-call syntax. `AddLayer` itself stays
# POSITIONAL (not keyword) because, unlike the constructor, its two arguments are always both
# given and never ambiguous about order.

## Single-angle, Ritz-only sweep (isotropic -> angle choice doesn't matter)
angle          = 0.0            # deg (isotropic: any angle gives the same result)
calcMode       = :ritzOnly      # :ritzOnly = SAW velocity only; :full also
                                 # runs Christoffel + Green's function
retentionLevel = :saw_only      # how much Ritz eigenvector data to keep --
                                 # :saw_only keeps just the classified SAW
                                 # mode(s), :full keeps every mode, :none
                                 # discards eigenvectors after classification
# `angles` is always a VECTOR in GuliwerrJ -- a single angle is a one-element vector, `[angle]`
# (a bare number is rejected; MATLAB's scalar-is-a-1x1-matrix rule has no Julia counterpart).
# `:ritzOnly`/`:saw_only` are Julia SYMBOLS -- a leading colon marks a lightweight, interned
# name, Julia's idiomatic replacement for MATLAB's quoted option strings ('ritzOnly',
# 'saw_only'). A Symbol is not a String: it prints without quotes, compares by identity (very
# cheap), and is the standard way this port spells a fixed menu of named options -- MATLAB's
# validated string. `nVelGreen` (the MATLAB tutorial's 5th sweep argument, left `[]` there
# because it is unused in `:ritzOnly` mode) is simply OMITTED below rather than passed as a
# placeholder: a Julia keyword argument you don't need, you just don't write -- there is no `[]`
# "use the default" convention to remember, because every keyword already has its own default.
sweep = angSweepHS(; structure = str, angles = [angle], calcMode = calcMode,
                      retentionLevel = retentionLevel)
sweep = calc(sweep)
# `calc(sweep)`, not `sweep.calc()`. Julia has no methods attached to objects -- every function
# is "free-standing", and which implementation runs is chosen from the TYPES of every argument
# together (this is "multiple dispatch", introduced properly in later tutorials once there is
# more than one type in play). Here there is only one `calc` method that matches a
# `CalcAngSweepHS` argument, so the effect is the same as MATLAB's `sweep = sweep.calc()` --
# `calc` returns the filled-in `CalcAngSweepHS`, and assigning it back to the same name keeps
# the script reading like MATLAB's `sweep = sweep.calc()`.

vSAW = sweep.sawVelMap[1][1]   # km/s
# `sawVelMap` is a Vector of Vectors -- one inner vector per swept angle (here, just one angle),
# each holding every SAW velocity found at that angle (there can be more than one on a general
# crystal cut). `[1][1]` reads exactly like MATLAB's `{1}(1)` did: first angle, first mode --
# 1-based indexing needs no adjustment at all coming from MATLAB.
vT   = sqrt(alC44 / alRho)     # transverse bulk velocity, km/s

@printf("v_SAW = %.4f km/s, v_T = %.4f km/s, ratio = %.3f\n", vSAW, vT, vSAW / vT)
# `@printf` is a MACRO (the leading `@` marks it as one -- a macro rewrites the code you wrote
# into different code before it runs, rather than being called at runtime like an ordinary
# function; you do not need to understand how that works to use it, just that the `@` is part
# of its name, not a typo). Its format-string mini-language (`%.4f`, `\n`, ...) is the same one
# C and MATLAB's `fprintf` both use, so every `%`-specifier in the original MATLAB tutorials
# translates over unchanged.
#
# Aside on trailing semicolons: MATLAB prints a bare statement's result to the console unless
# you end the line with `;` (`vSAW` above would otherwise echo `vSAW = 5.4321`). Julia scripts
# do NOT do this automatically -- only the interactive REPL prompt echoes an expression's
# value, and only when nothing has been done to suppress that display. A `.jl` file run via
# `julia script.jl`, or brought in via `include("script.jl")` from a running session, never
# auto-prints intermediate results, trailing `;` or not. You may still see a trailing `;` in
# Julia code (some style guides use it to mark "this line is an assignment, ignore its value"
# when working directly at the REPL), but it is not doing the job it does in MATLAB, and this
# tutorial series does not rely on it anywhere.
