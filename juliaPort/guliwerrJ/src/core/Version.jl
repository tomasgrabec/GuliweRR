# Version.jl -- port of guliwerr/guliwerrVersion.m.
#
# MATLAB and Julia ship together, so the port carries the same version as MATLAB: this
# string, `guliwerr/guliwerrVersion.m` and Project.toml's `version` are kept identical. A
# result computed by GuliwerrJ is meant to be the same result MATLAB would give, so it
# carries the same version stamp (MATLAB records it in every saved result as
# meta.guliwerrVersion).

"""
    guliwerrVersion() -> String

GuliweRR version string, `"MAJOR.MINOR.PATCH"` with an optional `-dev` suffix for code
between releases (e.g. `"6.1.0-dev"` is work towards 6.1.0). MAJOR marks a large
restructuring, MINOR new or changed functionality (may change or remove existing
functions/options), PATCH fixes only.

Port of `guliwerr/guliwerrVersion.m`: the MATLAB version this port reproduces.
"""
guliwerrVersion() = "6.0.0"

export guliwerrVersion
