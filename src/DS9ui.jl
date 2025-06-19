module DS9ui

using XPA
using OffsetArrays
using StaticArrays
using FITSIO
using Metaheuristics
using Printf
import REPL
using REPL.TerminalMenus

include("connection.jl")
include("interface.jl")
include("mask.jl")
include("pickobj.jl")

export ds9select, ds9, ds9close, ds9image, ds9image!, ds9delete, ds9pickobj

macro public(ex)
    if VERSION >= v"1.11.0-DEV.469"
        args = ex isa Symbol ? (ex,) : Base.isexpr(ex, :tuple) ? ex.args : error("something informative")
        esc(Expr(:public, args...))
    else
        nothing
    end
end

@public ds9cursor, ds9set, ds9get, ds9wcs, ds9fitprofiles, ds9mask, ds9getregions

end
