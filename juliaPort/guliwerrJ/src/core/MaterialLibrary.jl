# MaterialLibrary.jl -- port of Material.FromLibrary / Material.ListLibrary /
# Material.ReadLibraryFile (guliwerr/core/@Material/).
#
# Both languages read ONE data file, guliwerr/data/materialsLibrary.csv, so the numbers exist
# once and cannot drift between MATLAB and Julia. It is found from this file's own location
# (juliaPort/guliwerrJ/src/core -> repository root -> guliwerr/data); every tier ships
# guliwerr/ and juliaPort/ together, so the file is always there.
#
# Naming: MATLAB's static methods hang off the Material class (Material.FromLibrary); Julia's
# Material is a struct and cannot carry them, so they become free functions with the class
# name folded in -- MaterialFromLibrary, ListMaterialLibrary -- rather than the bare generic
# FromLibrary/ListLibrary.
#
# The CSV is parsed by hand (no CSV.jl dependency), exactly as ReadLibraryFile.m does:
# '#' lines are comments, the first other line is the header, fields are comma-separated and
# contain no commas themselves.

const _LIBRARY_PATH = joinpath(dirname(dirname(dirname(dirname(@__DIR__)))),
                               "guliwerr", "data", "materialsLibrary.csv")

"""
    _ReadLibraryFile() -> Vector{NamedTuple}

Read `guliwerr/data/materialsLibrary.csv` into one `NamedTuple` per data set, with fields
`name, aliases, formula, structure, symmNo, consts, rho, T, variant, default, status, reference, source,
note` (`consts::Vector{Float64}` in GPa, `symmNo::Int`, `rho::Float64` [g/cm^3],
`default::Bool`, the rest `String`). Port of `Material.ReadLibraryFile` (private there too).
"""
function _ReadLibraryFile()
    isfile(_LIBRARY_PATH) || error("Materials library not found: $(_LIBRARY_PATH)")
    lines = filter(l -> !isempty(strip(l)) && !startswith(strip(l), "#"),
                   readlines(_LIBRARY_PATH))
    header = Symbol.(strip.(split(lines[1], ',')))
    rows = NamedTuple[]
    for (r, line) in enumerate(lines[2:end])
        fields = strip.(split(line, ','))
        length(fields) == length(header) || error(
            "Library row $r has $(length(fields)) fields, the header has $(length(header)): $line")
        raw = NamedTuple{Tuple(header)}(Tuple(String.(fields)))
        push!(rows, (name = raw.name, aliases = raw.aliases, formula = raw.formula,
                     structure = raw.structure,
                     symmNo = parse(Int, raw.symmNo),
                     consts = parse.(Float64, split(raw.consts)),
                     rho = parse(Float64, raw.rho), T = raw.T, variant = raw.variant,
                     default = raw.default == "1", status = raw.status,
                     reference = raw.reference, source = raw.source, note = raw.note))
    end
    return rows
end

"""
    ListMaterialLibrary(; print=false) -> Vector{NamedTuple}

All data sets of the materials library (`guliwerr/data/materialsLibrary.csv`), one
`NamedTuple` per row (fields as in [`MaterialFromLibrary`](@ref)'s `entry`). With
`print=true` a short overview is also printed. Port of `Material.ListLibrary` (which prints
when called without an output argument).
"""
function ListMaterialLibrary(; print::Bool = false)
    lib = _ReadLibraryFile()
    if print
        Base.println(rpad("name", 14), rpad("structure", 14), rpad("symm", 6),
                     rpad("consts [GPa]", 45), rpad("T", 6),
                     rpad("rho", 9), rpad("variant", 9), rpad("status", 13), "default")
        for e in lib
            Base.println(rpad(e.name, 14), rpad(e.structure, 14), rpad(string(e.symmNo), 6),
                         rpad(join(_fmtg.(e.consts), " "), 45), rpad(e.T, 6),
                         rpad(_fmtg(e.rho), 9), rpad(e.variant, 9), rpad(e.status, 13),
                         e.default ? "*" : "")
        end
    end
    return lib
end

_fmtg(x::Real) = isinteger(x) ? string(Int(x)) : string(x)

"""
    MaterialFromLibrary(name; dirMat=nothing, T=nothing, variant=nothing) -> (Material, entry)

Build a [`Material`](@ref) from the materials library (`guliwerr/data/materialsLibrary.csv`).

- `name`    -- library name or alias, case-insensitive (`"Cu"`, `"copper"`)
- `dirMat`  -- orientation matrix (see [`Material`](@ref)); `nothing`: Material's default
- `T`       -- temperature, a number in K or `"RT"`; must match a stored row
- `variant` -- variant label of the row ([`ListMaterialLibrary`](@ref) shows them)

Without `T` or `variant` the material's default row is used. Returns the `Material` (named
after the library entry, plus variant and T when not the default row) and the row's
`NamedTuple` (`reference`, `source`, `note`, ...). For rows with status `provisional`
(secondary source only, or a doubt given in the row's note) a one-line note saying why is
logged (`@info`; MATLAB prints it). Port of `Material.FromLibrary`.
"""
function MaterialFromLibrary(name::AbstractString; dirMat = nothing, T = nothing,
                             variant = nothing)
    lib = _ReadLibraryFile()
    key = lowercase(name)
    rows = filter(e -> lowercase(e.name) == key ||
                       any(lowercase.(split(e.aliases, '|')) .== key), lib)
    isempty(rows) && error("No material \"$name\" in the library. Available: ",
                           join(unique(e.name for e in lib), ", "), ".")

    sel = if T === nothing && variant === nothing
        filter(e -> e.default, rows)
    else
        filter(e -> (T === nothing || lowercase(e.T) == lowercase(_fmtT(T))) &&
                    (variant === nothing || lowercase(e.variant) == lowercase(variant)), rows)
    end
    available = join(("T=$(e.T) variant=$(e.variant)" for e in rows), "; ")
    isempty(sel) && error("No \"$(rows[1].name)\" row matches the request. Available rows: ",
                          available, ".")
    length(sel) > 1 && error("Several \"$(rows[1].name)\" rows match; give T and/or variant. ",
                             "Available rows: ", available, ".")
    entry = sel[1]

    matName = entry.default ? entry.name : "$(entry.name) ($(entry.variant), T=$(entry.T))"
    # why it is provisional: no primary reference yet (secondary source only), or a doubt
    # stated in the note
    why = entry.reference == "pending" ? "secondary source only: $(entry.source)" : entry.note
    entry.status == "provisional" && @info "$matName -- provisional library values ($why)"
    mat = Material(elConsts = ElastConsts(entry.consts, entry.symmNo), rho = entry.rho,
                   dirMat = dirMat, name = matName)
    return mat, entry
end

_fmtT(T::Real) = _fmtg(T)
_fmtT(T::AbstractString) = String(T)

export MaterialFromLibrary, ListMaterialLibrary
