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

public ds9cursor, ds9set, ds9get, ds9wcs, ds9fitprofiles, ds9mask, ds9getregions

end
