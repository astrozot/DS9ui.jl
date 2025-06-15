# Connection

The connection to the SAOImage DS9 is performed using the XPA interface. By
default, all routines connect to the first DS9 window found. Optionally,
however, the user can select a different window using [`ds9select`](@ref).

All further operations will be performed by default over the connected window
(i.e., the current `access_point`).

## User interface

```@docs
DS9ui.ds9
DS9ui.ds9select
DS9ui.ds9close
```

## Internal routines and variables

```@docs
DS9ui.access_point
DS9ui.connect
DS9ui._current_ap
```
