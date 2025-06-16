function ray_intersect_line(x, p1, p2)
    δ1 = p2[1] - p1[1]
    δ2 = p2[2] - p1[2]
    if (x[1] - p1[1]) * δ2 == (x[2] - p1[2]) * δ1 && min(p1[1], p2[1]) <= x[1] <= max(p1[1], p2[1])
        -1
    else
        point_inside = ((p1[2] > x[2]) ⊻ (p2[2] > x[2])) && x[1] < δ1 * (x[2] - p1[2]) / δ2 + p1[1]  
        point_inside ? 1 : 0
    end
end

function point_in_poly(point, vertices)
    minx = minimum(first, vertices)
    miny = minimum(last, vertices)
    maxx = maximum(first, vertices)
    maxy = maximum(last, vertices)
    if point[1] < minx || point[1] > maxx || point[2] < miny || point[2] > maxy
        return 0
    else
        c = 0
        for i in 1:length(vertices)
            p1 = vertices[i]
            p2 = vertices[mod1(i + 1, length(vertices))]
            intersection = ray_intersect_line(point, p1, p2)
            if intersection == -1
                return -1
            end
            c += intersection
        end
        return c % 2
    end
end


"""
    ds9mask([access_point,] name::String=""; coords=:image, silent=true, 
            full=false, usefits=false)

Return a mask of the regions defined in the DS9 window.

The optional argument `name` is the name of the group of regions to extract.
If `name` is an empty string, all regions are extracted.

# Keyword Arguments
- `ident`: the identifier of the DS9 window.
- `coords`: the type of coordinates to return: can be `:image`, `:physical`,
  `:fk5`, `:galactic`
- `silent`: if `true`, no warning is issued when unknown regions are found.
- `full`: if `true`, the mask is the full image. If `false`, the mask is
  cropped to the smallest rectangle containing all regions.
- `usefits`: if `true`, the mask is extracted from the FITS file. If `false`,
  the mask is extracted from the DS9 window.

The return value is a tuple `(image, mask)`, where `image` is the image
extracted from the DS9 window (that is, the pixel values), and `mask` is a
boolean `OffsetArray` with the mask (that is, `true` inside the regions and
`false` outside).
"""
function ds9mask(ap, name; coords=:image, selected=false, 
    silent=true, full=false, usefits=false)
    reply = XPA.get(ap, "fits header")
    regions = ds9getregions(ap, name; coords, selected)
    if full
        header = split(XPA.get_data(String, reply), "\n")
        ix = findfirst(startswith("NAXIS1  = "), header)
        iy = findfirst(startswith("NAXIS2  = "), header)
        if !isnothing(ix) && !isnothing(iy)
            nx = parse(Int, header[ix][10:30])
            ny = parse(Int, header[iy][10:30])
            xrange = 1:nx
            yrange = 1:ny
        end
    else
        x₀, x₁ = Inf, -Inf
        y₀, y₁ = Inf, -Inf
        for region ∈ regions
            shape, coordinates, props = region
            if shape == :circle
                x, y, r = coordinates
                x₀ = min(x₀, x - r)
                x₁ = max(x₁, x + r)
                y₀ = min(y₀, y - r)
                y₁ = max(y₁, y + r)
            elseif shape == :ellipse
                x, y = coordinates[1:2]
                a, b, α = coordinates[end-2:end]
                r = max(a, b)
                x₀ = min(x₀, x - r)
                x₁ = max(x₁, x + r)
                y₀ = min(y₀, y - r)
                y₁ = max(y₁, y + r)
            elseif shape == :box
                x, y = coordinates[1:2]
                a, b, α = coordinates[end-2:end]
                r = ceil(sqrt(a^2 + b^2) / 2)
                x₀ = min(x₀, x - r)
                x₁ = max(x₁, x + r)
                y₀ = min(y₀, y - r)
                y₁ = max(y₁, y + r)
            elseif shape == :polygon
                npts = length(coordinates) ÷ 2
                for i in 1:npts
                    xᵢ = coordinates[2i - 1]
                    yᵢ = coordinates[2i]
                    x₀ = min(x₀, xᵢ)
                    x₁ = max(x₁, xᵢ)
                    y₀ = min(y₀, yᵢ)
                    y₁ = max(y₁, yᵢ)
                end
            elseif shape == :annulus
                x, y = coordinates[1:2]
                r = coordinates[end]
                x₀ = min(x₀, x - r)
                x₁ = max(x₁, x + r)
                y₀ = min(y₀, y - r)
                y₁ = max(y₁, y + r)
            end
        end
        if x₀ == Inf
            @error "No valid regions found"
            return nothing
        end
        xrange = floor(Int, x₀):ceil(Int, x₁)
        yrange = floor(Int, y₀):ceil(Int, y₁)
    end
    mask = zeros(Bool, xrange, yrange)
    unmask = ones(Bool, xrange, yrange)
    xrange′ = OffsetArray(xrange, first(xrange) - 1)
    yrange′ = OffsetArray(yrange, first(yrange) - 1)
    unknown_regions = Symbol[]
    for region ∈ regions
        shape, coordinates, props = region
        if shape == :circle || shape == Symbol("-circle")
            x, y, r = coordinates
            if shape == :circle
                mask .= mask .| ((xrange′ .- x) .^ 2 .+ (yrange′' .- y) .^ 2 .< r^2)
            else
                unmask .= unmask .& ((xrange′ .- x) .^ 2 .+ (yrange′' .- y) .^ 2 .>= r^2)
            end
        elseif shape == :ellipse || shape == Symbol("-ellipse")
            sinα, cosα = sincosd(coordinates[end])
            x, y = coordinates[1:2]
            a, b = coordinates[end-2:end-1]
            xₒ = xrange′ .- x
            yₒ = yrange′' .- y
            xₒs = xₒ * cosα .+ yₒ * sinα
            yₒs = -xₒ * sinα .+ yₒ * cosα
            e⁻¹ = a / b
            r² = xₒs .^ 2 .+ (yₒs .* e⁻¹) .^ 2
            inside = r² .< a^2  # Inside the ellipse
            inout = true
            for i ∈ length(coordinates)-4:-2:3
                if inout
                    inside .= inside .& (r² .>= coordinates[i]^2)
                else
                    inside .= inside .| (r² .< coordinates[i]^2)
                end
                inout = !inout
            end
            if shape == :ellipse
                mask .= mask .| inside
            else
                unmask .= unmask .& (.!inside)
            end
        elseif shape == :box || shape == Symbol("-box")
            sinα, cosα = sincosd(coordinates[end])
            x, y = coordinates[1:2]
            a, b = coordinates[end-2:end-1]
            xₒ = xrange′ .- x
            yₒ = yrange′' .- y
            xₒs = xₒ * cosα .+ yₒ * sinα
            yₒs = -xₒ * sinα .+ yₒ * cosα
            inside = (abs.(2 * xₒs) .< a) .& (abs.(2 * yₒs) .< b)
            inout = true
            for i ∈ length(coordinates)-4:-2:3
                a, b = coordinates[i:i+1]
                if inout
                    inside .= inside .& ((abs.(2 * xₒs) .>= a) .| (abs.(2 * yₒs) .>= b))
                else
                    inside .= inside .| ((abs.(2 * xₒs) .< a) .& (abs.(2 * yₒs) .< b))
                end
                inout = !inout
            end
            if shape == :box
                mask .= mask .| inside
            else
                unmask .= unmask .& (.!inside)
            end
        elseif shape == :polygon || shape == Symbol("-polygon")
            vertices = reinterpret(SVector{2,Float64}, coordinates)
            if shape == :polygon
                mask .= mask .| (point_in_poly.(SVector{2}.(xrange′, yrange′'), Ref(vertices)) .== 1)
            else
                unmask .= unmask .& (point_in_poly.(SVector{2}.(xrange′, yrange′'), Ref(vertices)) .== 0)
            end
        elseif shape == :annulus || shape == Symbol("-annulus")
            x, y = coordinates[1:2]
            r² = (xrange′ .- x) .^ 2 .+ (yrange′' .- y) .^ 2
            inside = r² .< coordinates[end]^2
            inout = true
            for i ∈ length(coordinates)-1:-1:3
                if inout
                    inside .= inside .& (r² .>= coordinates[i]^2)
                else
                    inside .= inside .| (r² .< coordinates[i]^2)
                end
                inout = !inout
            end
            if shape == :annulus
                mask .= mask .| inside
            else
                unmask .= unmask .& (.!inside)
            end
        else
            if shape ∉ unknown_regions
                push!(unknown_regions, shape)
            end
        end
    end
    if !silent && !isempty(unknown_regions)
        @warn "Unknown regions found: $unknown_regions"
    end
    # Perform the unmasking
    mask = mask .& unmask
    # Shrink the mask as much as possible
    if !full
        xmask = sum(mask, dims=2)
        ymask = sum(mask, dims=1)
        try
            x₀ = findfirst(>(0), xmask)[1]
            x₁ = findlast(>(0), xmask)[1]
            y₀ = findfirst(>(0), ymask)[2]
            y₁ = findlast(>(0), ymask)[2]
        catch e
            @error "Empty area"
            return nothing
        end
        mask = OffsetArrays.OffsetArray(mask[x₀:x₁, y₀:y₁], (x₀, y₀))
    end
    # Get the image
    if usefits
        local image, path
        try
            path = tempname()
            reply = XPA.set(ap, "save fits $path")
            f = FITS(path, "r")
            fullimage = read(f[1])
            close(f)
            image = fullimage[axes(mask)...]
        catch
            rm(path, force=true)
            @error "Could not create temporary file"
            return nothing
        end
        rm(path, force=true)
    else
        xrange = extrema(axes(mask, 1))
        yrange = extrema(axes(mask, 2))
        reply = XPA.get(ap,
            "data image $(xrange[1]) $(yrange[1]) $(xrange[2] - xrange[1] + 1) $(yrange[2] - yrange[1] + 1) no")
        data = split(XPA.get_data(String, reply), "\n")
        image = zeros(axes(mask)...)
        for v ∈ split.(data, r",|=")
            if length(v) < 3
                continue
            end
            v′ = parse.(Float64, v)
            image[round(Int, v′[1]), round(Int, v′[2])] = v′[3]
        end
    end
    image, mask
end

@inline ds9mask(name=""; kw...) = ds9mask(_current_ap(), name; kw...)