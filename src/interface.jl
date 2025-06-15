"""
    ds9set([access_point], "command"; kw...)
    ds9set([access_point], (; key1=val1, key2=val2...))

Transmits one or more commands to DS9 using the XPA interface.

In the first form, the command is entered as a single string and is
transmitted unchanged to DS9. Alternatively, one can pass a NamedTuple with a
series of keys and values. Each keywords result in a single command sent to
DS9. The value of the keyword is directly used without any interpolation. The
only exception is that a Boolean value is converted into "yes"/"no" (following
DS9 notation).
"""
function ds9set(ap, command::AbstractString; kw...)
    r = XPA.set(ap, command; kw...)
    if XPA.has_errors(r)
        m = XPA.get_message(r)
        msg = strip(m[10:end])
        f = findfirst("(DS9:", msg)
        if !isnothing(f)
            msg = msg[begin:first(f)-1]
        end
        @warn "XPA $msg"
    elseif r.replies == 0
        @warn "No replies for command `$cmd`"
    end
end
ds9set(command; kw...) = ds9set(_current_ap(), command; kw...)

function ds9set(ap, commands::AbstractDict; kw...)
    for (k, v) ∈ commands
        w = isa(v, Bool) ? (v ? "yes" : "no") : string(v)
        command = "$k $w"
        ds9set(ap, command; kw...)
    end
end
ds9set(ap, commands::NamedTuple; kw...) = ds9set(ap, pairs(commands); kw...)


"""
    r = ds9get([type [, dims],] [access_point,] command; kw...)

Returns the output of an XPA get command sent to DS9.

For complex outputs, one can use the `type` and optionally the `dims`
parameters, which are directly passed to `XPA.get`.
"""
function ds9get(ap, command; kw...)
    if string(command) != ""
        r = XPA.get(ap, string(command); kw...)
        if XPA.has_errors(r)
            m = XPA.get_message(r)
            msg = strip(m[10:end])
            f = findfirst("(DS9:", msg)
            if !isnothing(f)
                msg = msg[begin:first(f)-1]
            end
            @warn "XPA $msg"
        elseif r.replies == 0
            @warn "No replies for command `$command`"
        end
        data = XPA.get_data(String, r)
        lines = filter(!isempty, split(data, "\n"))
        if isnothing(findfirst(x -> isnothing(tryparse(Int, x)), lines))
            result = tryparse.(Int, lines)
        elseif isnothing(findfirst(x -> isnothing(tryparse(Float64, x)), lines))
            result = tryparse.(Float64, lines)
        else
            result = lines
        end
        if length(result) == 0
            return nothing
        elseif length(result) == 1
            return first(result)
        else
            return result
        end
    end
    nothing
end
ds9get(command; kw...) = ds9get(_current_ap(), command; kw...)

@inline ds9get(::Type{X}, ap, command; kw...) where {X} = XPA.get(X, ap, command; kw...)
@inline ds9get(::Type{X}, command; kw...) where {X} = XPA.get(X, _current_ap(), command; kw...)
@inline ds9get(::Type{X}, dims::Union{Integer,NTuple{N,Integer}}, ap, command; kw...) where {X,N} =
    XPA.get(X, dims, ap, command; kw...)
@inline ds9get(::Type{X}, dims::Union{Integer,NTuple{N,Integer}}, command; kw...) where {X,N} =
    XPA.get(X, dims, _current_ap(), command; kw...)


"""
    ds9cursor([access_point]; coords=:image, event=:button)

Return the coordinates of the cursor in the DS9 window.

- `ident`: the identifier of the DS9 window.
- `coords`: the type of coordinates to returnL can be `:image`, `:physical`,
  `:fk5`, `:galactic`
- `event`: the type of event to capture the cursor position. Can be `:button`,
  `:key`, `:any`

The function returns a tuple `(key, coords)`, where `key` is a string
representing the event and `coords` is an array of coordinates as float
numbers.
"""
function ds9cursor(ap=_current_ap(); coords=:image, event=:button)
    if event ∉ (:button, :key, :any)
        error("Unknown event type $event")
    end
    XPA.set(ap, "raise")
    reply = XPA.get(ap, "imexam $event coordinate $coords")
    data = split(XPA.get_data(String, reply))
    if event != :button
        key = data[1]
        p = parse.(Float64, @view data[2:end])
    else
        key = "<1>"
        p = parse.(Float64, data)
    end
    key, p
end

"""
    ds9wcs([access_point]; useheader=true)

Return the WCS transformation of the current image in the DS9 window.

# Keyword Arguments
- `ident`: the identifier of the DS9 window.
- `useheader`: if `true`, the WCS is extracted from the FITS header. If
  `false`, the WCS is extracted from the DS9 window.
"""
function ds9wcs(ap=_current_ap(); useheader=true)
    if useheader
        reply = XPA.get(ap, "fits header")
        header = XPA.get_data(String, reply)
    else
        path = tempname()
        reply = XPA.set(ap, "wcs save $path")
        header = open(path) do f
            read(f, String)
        end
        rm(path, force=true)
    end
    header
    # FIXME: Gravity.wcs_from_header(Gravity.fits_split_header(replace(header, "\n" => "")))
end

"""
    ds9image([access_point,] image; usefile=true, wcs=nothing)

Display the `image` in the current frame of DS9.
"""
function ds9image(ap, image; usefile=true, wcs=nothing)
    if usefile
        path = tempname()
        try
            # FIXME
            Gravity.writefits(image, path, wcs=wcs, silent=true)
            reply = XPA.set(ap, "file $path")
        finally
            rm(path, force=true)
        end
    else
        factor = isa(image, AbstractMatrix{<:AbstractFloat}) ? -8 : 8
        reply = XPA.set(ap, 
            "array new[xdim=$(size(image, 2)),ydim=$(size(image, 1)),bitpix=$(factor * sizeof(eltype(image)))]";
            data=image)
        if !isnothing(wcs)
            # FIXME
            XPA.set(ap, "wcs replace"; data=Gravity.fits_join_header(wcs_to_header(wcs)))
        end
    end
end
@inline ds9image(image; kw...) = ds9image(_current_ap(), image; kw...)


"""
    ds9getregions([access_point,] name=""; coords=:image, selected=false)

Return the regions defined in the DS9 window.

The optional argument `name` is the name of the group of regions to extract.
If `name` is an empty string, all regions are extracted.

# Keyword Arguments
- `ident`: the identifier of the DS9 window.
- `coords`: the type of coordinates to returnL can be `:image`, `:physical`,
  `:fk5`, `:galactic`

The return value is an array of tuples `(shape, coordinates, properties)`,
where `shape` is a symbol indicating shape of the region, `coordinates` is an
array of coordinates, and `properties` is a dictionary with the properties of
the region.
"""
function ds9getregions(ap, name::String; coords=:image, selected=false)
    function parsevalue(s)
        r = tryparse(Int, s)
        if r === nothing
            r = tryparse(Float64, s)
            if r === nothing
                r = s
            end
        end
        r
    end
    # Query DS9 for the regions and the selections
    group = (name != "") ? "-group $name" : ""
    sel = selected ? "selected" : ""
    reply = XPA.get(ap, "regions $sel -format ds9 $group -system $coords")
    regs = split(XPA.get_data(String, reply), "\n")
    # Extract the global settings
    global_lines = filter(line -> startswith(line, "global "), regs)
    global_props = Dict{String, Any}()
    for line in global_lines
        for m ∈ eachmatch(r"\b(\w+)\b=\b([\w.]+)\b", line)
            global_props[m.captures[1]] = parsevalue(m.captures[2])
        end
    end
    regions = Tuple{Symbol, Vector{Float64}, Dict{String, Any}}[]
    for line ∈ regs
        m = match(r"^\s*(-?circle|-?ellipse|-?box|-?polygon|-?annulus|-?panda|-?epanda)\(([^)]+)\)\s*(#.*)?$", line)
        if m !== nothing
            shape = Symbol(m.captures[1])
            coords = parse.(Float64, split(m.captures[2], ","))
            comment = isnothing(m.captures[3]) ? "" : m.captures[3][2:end]
            local_props = copy(global_props)
            for m′ ∈ eachmatch(r"\b(\w+)\b=\b([\w.]+)\b", comment)
                local_props[m′.captures[1]] = parsevalue(m′.captures[2])
            end
            push!(regions, (shape, coords, local_props))
        end
    end
    return regions
end
@inline ds9getregions(name::String=""; kw...) = ds9getregions(_current_ap(), name; kw...)


