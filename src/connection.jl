"""A string saving the XPA address of the current access point, or "" if no
connection is established.
"""
access_point::String = ""

"""
    _current_ap()

The default access point.

This function just returns [`access_point`](@ref) if this variable is set, or
establishes a new connection using [`connect`](@ref).
"""
function _current_ap()
    global access_point
    if access_point == ""
        try
            access_point = connect()
        catch
            @warn "Failed to connect to SAOImage DS9"
        end
    end
    return access_point
end

"""
    connect(ident="DS9:*"; method="local")

Establish a new connection to DS9 using the given `ident`.

If `method` is not `undefined`, the environment variable `XPA_METHOD` is set
to the corresponding string. Typically the most reliable connection is
`local`, the default value.
"""
function connect(ident::Union{Regex,AbstractString}="DS9:*"; method="local", kw...)
    global access_point
    if !haskey(ENV, "XPA_METHOD") && !isnothing(method)
        ENV["XPA_METHOD"] = string(method)
    end
    apt = XPA.find(ident; kw...)
    apt === nothing && error("no matching SAOImage/DS9 server found")
    rep = XPA.get(apt, "version"; nmax=1)
    if length(rep) != 1 || !XPA.verify(rep)
        error("XPA server at address \"$apt\" is not a valid SAOImage DS9 server")
    end
    addr = XPA.address(apt)
    access_point = addr
end


"""
    ds9select(ident="DS9:*"; method="local"; silent=false, interactive=true)

Select a DS9 window for further interactions.

If multiple DS9 windows match the `ident`, the user can select the correct
window using a simple interface; if `interactive=false` the first matching
window is selected.

!!! warning
    If `method` is not `undefined`, the environment variable `XPA_METHOD` is set
    to the corresponding string. Typically the most reliable connection is
    `local`, the default value.
"""
function ds9select(ident::Union{Regex,AbstractString}="DS9:*"; method="local",
    silent=false, interactive=true)
    global access_point
    if !haskey(ENV, "XPA_METHOD") && !isnothing(method)
        ENV["XPA_METHOD"] = string(method)
    end
    i = findfirst(isequal(':'), ident)
    if i === nothing
        # allow any class
        class = "*"
        name = ident
    else
        class = ident[1:i-1]
        name = ident[i+1:end]
    end
    apts = XPA.list()
    good_apts = filter(a -> (name == "*" || name == a.name) && (class == "*" || class == a.class), apts)
    if length(good_apts) == 0
        if !silent
            @error "No valid XPA access point found (server not reachable?)"
        end
        return nothing
    else
        if length(good_apts) == 1 || !interactive
            choice = 1
        else
            menu = RadioMenu(["$(p.class):$(p.name) (user=$(p.user))" for p ∈ apts])
            choice = request("Please select the correct access point:", menu)
            if choice == -1
                return nothing
            end
        end
        if !silent
            class = good_apts[choice].class 
            name = good_apts[choice].name
            user = good_apts[choice].user
            @info "Connected to the XPA access point $class:$name (user=$user)"
        end
        access_point = XPA.address(good_apts[choice])
        return nothing
    end
end


"""
    ds9([name]; method="local", path="ds9")

Launch the DS9 application.

By default the application will be named using the current PID.

The `method` optional keywords is used to set the XPA communication method:
"local" is the recommended way for local executions.

This function authomaticall sets the default access point, so that all further
requests are forwarded to the newly open DS9 window.
"""
function ds9(name::String=string(getpid()); method="local", path="ds9", silent=false)
    global access_point
    if !haskey(ENV, "XPA_METHOD") && !isnothing(method)
        ENV["XPA_METHOD"] = string(method)
    end
    command = detach(`$path -xpa yes -xpa connect -title $name`)
    run(command; wait=false)
    if !silent
        printstyled("[ Info: "; color=Base.default_color_info, bold=true)
        print("Opening DS9")
    end
    for i ∈ 1:20
        silent || print(".")
        sleep(0.4)
        apt = XPA.find("DS9:$name")
        if !isnothing(apt)
            access_point = XPA.address(apt)
            silent || println(" done")
            return
        end
    end
    silent || println(" failed")
    @warn "Timeout establishing an XPA connection."
end


"""
    ds9close([access_point])

Close the current DS9 session.
"""
function ds9close(ap=_current_ap())
    global access_point
    ds9set(ap, "exit")
    if access_point == ap
        access_point = ""
    end
    nothing
end



