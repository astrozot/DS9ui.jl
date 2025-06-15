"""
    ds9fitprofiles([access_point,] name; modes=[:gaussian],
                   usefits=false, silent=false, max_calls=1_000_000)

Fit a region in DS9 with elliptical profiles.

This function extract a region in DS9 and then fits the data with a
superposition of profiles identified by `modes`. A background value for the
sky is always added. Known profiles are `:gaussian` (the default), `:moffat`,
`:lorentzian`, `:sersic`, `:exponential`, `:devaucouleurs`, and `:king`. Note
that `modes` is an array: all given profiles will be combined in the specified
order for the final fit.

Other options control the way data are retrieved from DS9 (`usefits`), the
output to print (`silent`), and the maximum number of calls to the chi-square
funtion used for the fit (`max_calls`).

The function returns the best-fit parameters (as a vector) and a matrix with
the residuals (data - model).
"""
function ds9fitprofiles(ap, name; modes=[:gaussian], usefits=false,
    silent=false, max_calls = 1_000_000)
    function modelpred!(image_out, pars)
        sky = pars[1]
        image_out .= sky
        i = 1
        xrange = extrema(axes(image_out, 1))
        yrange = extrema(axes(image_out, 2))
        for mode ∈ modes
            norm, xₚ, yₚ, a, p, α = pars[1 + i:6 + i]
            b = a * p
            sinα, cosα = sincos(α)
            xₒ = OffsetArray(xrange[1]:xrange[2], xrange[1] - 1) .- xₚ
            yₒ = OffsetArray(yrange[1]:yrange[2], yrange[1] - 1)' .- yₚ
            xₒs = xₒ .* cosα .+ yₒ .* sinα
            yₒs = -xₒ .* sinα .+ yₒ .* cosα
            r² = (xₒs ./ a) .^ 2 .+ (yₒs ./ b) .^ 2
            if mode == :gaussian
                image_out .+= norm .* exp.(-0.5 * r²)
            elseif mode == :lorentzian
                image_out .+= norm ./ (1 .+ r²)
            elseif mode == :moffat
                image_out .+= norm .* (1 .+ r²).^(-pars[7 + i])
                i += 1
            elseif mode == :sersic
                image_out .+= norm .* exp.(-r².^(0.5/pars[7 + i]))
                i += 1
            elseif mode == :exponential
                image_out .+= norm .* exp.(-r².^0.5)
            elseif mode == :devaucouleurs
                image_out .+= norm .* exp.(-r².^0.125)
            elseif mode == :king
                image_out .+= norm .* (1 .+ r²).^(-1)
            else
                error("Unknown mode $mode")
            end
            i += 6
        end
        image_out
    end
    χ²(pars) = sum(abs2, (image .- modelpred!(image_out, pars)) .* weights)
    # Query DS9 for the regions and the selections
    image, mask = ds9mask(ap, name; usefits=usefits) 
    weights = Float64.(mask)
    image_out = copy(image)
    xrange = extrema(axes(image, 1))
    yrange = extrema(axes(image, 2))
    # Perform the fitting
    lbs = Float64[minimum(image)] # sky
    ubs = Float64[maximum(image)]
    for mode ∈ modes
        # parameters: norm, x, y, a, p=b/a, theta
        lb = [0.0, xrange[1], yrange[1], 0.1, 0.01, 0]
        rrange = sqrt((xrange[2] - xrange[1])^2 + (yrange[2] - yrange[1])^2)
        ub = [maximum(image), xrange[2], yrange[2], rrange, 1.0, π]
        if mode == :moffat
            push!(lb, 0.5)
            push!(ub, 5.0)
        elseif mode == :sersic
            push!(lb, 0.5)
            push!(ub, 10.0)
        end
        append!(lbs, lb)
        append!(ubs, ub)
    end
    bounds = Metaheuristics.BoxConstrainedSpace(lb=lbs, ub=ubs)
    options = Metaheuristics.Options(f_calls_limit=max_calls, f_tol=1e-5)
    algorithm = Metaheuristics.ECA(options=options)
    result = Metaheuristics.optimize(x -> sum(abs2, (image .- modelpred!(image_out, x)) .* weights), 
        bounds, algorithm)
    best = result.best_sol.x
    if !silent
        i = 0
        m_str = length(modes) == 1 ? "" : "   "
        @printf "Sky level:%s  %.6g adu\n" m_str best[1]
        for (m, mode) ∈ enumerate(modes)
            m_str = length(modes) == 1 ? "" : "[$m]"
            t_str = length(modes) == 1 ? "" : "  "
            @printf "Mode%s:%s       %s\n" m_str t_str mode
            @printf "Center%s:     (%.2f, %.2f) pix\n" m_str best[i + 3] best[i + 4]
            if mode == :gaussian
                @printf "%sFWHM%s:       (%.2f, %.2f) pix\n" t_str m_str (best[i + 5]*sqrt(8*log(2))) (best[i + 5]*best[i + 6]*sqrt(8*log(2)))
                @printf "%sSigma%s:      (%.2f, %.2f) pix\n" t_str m_str (best[i + 5]) (best[i + 5]*best[i + 6])
            else
                @printf "%sScales%s:     (%.2f, %.2f) pix\n" t_str m_str (best[i+5]) (best[i+5] * best[i+6])
            end
            @printf "%sAxis ratio%s: %.4f\n" t_str m_str best[i + 6]
            @printf "%sAngle%s:      %.1f°\n" t_str m_str (best[i+7] * 180 / π)
            @printf "%sPeak%s:       %.5g adu\n" t_str m_str best[i + 2]
            @printf "%sFlux%s:       %.5g adu pix²\n" t_str m_str (sum(image) - best[1])
            if mode ∈ (:moffat, :sersic)
                @printf "%sIndex%s:     %.5g\n" t_str m_str best[i+8]
                i += 1
            end
            i += 6
        end
    end
    return best, image .- modelpred!(image_out, best)
end
@inline ds9fitprofiles(name::AbstractString=""; kw...) = ds9fitprofiles(_current_ap(), name; kw...)


"""
    ds9pickobj([access_point]; ident="DS9:*", mode=:gaussian, r=32, 
               shape=:circle, silent=false, usefits=false)

Fits interactively a region in DS9 with an elliptical profile.

This function allows the selection of a circular or quadratic region in DS9
and fits it with a given profile using [`ds9fitprofiles`](@ref).
"""
function ds9pickobj(ap=_current_ap(); mode=:gaussian, r=32, shape=:circle, 
    silent=false, usefits=false)
    function showregion()
        if shape == :circle
            XPA.set(ap, "region",
                data="image; circle($(coords[1]), $(coords[2]), $(round(r) + 0.5)) # tag={ds9pickobj} width=2")
        else
            XPA.set(ap, "region",
                data="image; box($(coords[1]), $(coords[2]), $(round(2*r) + 0.5), $(round(2*r) + 0.5)) # tag={ds9pickobj} width=2")
        end
    end
    function hideregion()
        XPA.set(ap, "region group ds9pickobj delete")
    end

    XPA.set(ap, "raise")
    modenames = [:sky, :stat, :gaussian, :lorentzian, :moffat,
        :sersic, :exponential, :devaucouleurs, :king]
    if !silent
        @info "Mouse: pick object - ?/H: help - Esc/Q: quit"
    end

    local coords
    while true
        key, coords = ds9cursor(ap; event=:any)
        if key ∈ ("minus", "plus", "equal", "underscore", "space", "slash")
            if key == "minus"
                r = max(r / √2, 1)
            elseif key == "underscore"
                r = max(r - 1, 1)
            elseif key == "equal"
                r = r * √2
            elseif key == "plus"
                r = r + 1
            elseif key == "slash"
                if shape == :box
                    shape = :circle
                else
                    shape = :box
                end
            end
            showregion()
            sleep(key ∈ ("underscore", "plus") ? 0.1 : 0.3)
            hideregion()
        elseif key ∈ ("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")
            nmode = parse(Int, key)
            if nmode > length(modenames)
                println("Unknown fitting profile; type ? for help")
            else
                mode = modenames[nmode]
                println("Fitting function: $mode")
            end
        elseif key ∈ ("question", "h", "H")
            println("Plus/Minus/Slash: change size and shape - Space: show fit box")
            println(join(("$(i): $(modenames[i])" for i in 1:length(modenames)), "; "))
        elseif key ∈ ("<1>", "Return")
            break
        elseif key ∈ ("Escape", "q", "Q")
            if !silent
                println("Quitting...")
            end
            return -1
        end
    end

    showregion()
    modes = mode == :sky ? [] : [mode]
    res = ds9fitprofiles(ap, "ds9pickobj"; modes=modes, usefits=usefits)
    hideregion()
    res
end

