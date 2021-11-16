using VectorFEM

export upper, lower, cell, limits

"""
    limits

Possible intergral limits:
1. `upper`: upper limit;
2. `lower`: lower limit.
"""
@enum limits begin
    upper
    lower
end  # limits

"""
    cellmethod(f::Function, rlimits::Dict{limits, Real}, slimits::Dict{limits, Real}, 
        n_r_sections::Int, n_s_sections::Int)

Multiple integral calculation with cell method. f has to be depended on r and s only. 
So other functions have to be converted to such format before calling this method.
This method is very simple but not effective. This method was implemented for tests, 
so it's recommended to use Gauss integration instead.

# Arguments
- `f`: function that needs to be integrated, need to be depended on 2 variables: 
``f = f(r, s)``;
- `rlimits::Dict{limits, Real}`: limits for the first variable;
- `slimits::Dict{limits, Real}`: limits for the second variable;
- `n_r_sections::Int`: number of sections on ``r`` axis;
- `n_s_sections::Int`: number of sections on ``s`` axis.
"""
function cellmethod(f::Function, rlimits::Dict{limits, Real}, slimits::Dict{limits, Real}, 
        n_r_sections::Int, n_s_sections::Int
)
    hr = (rlimits[upper] - rlimits[lower]) / n_r_sections
    hs = (slimits[upper] - slimits[lower]) / n_s_sections
    rpoints = collect(Float64, rlimits[lower]:hr:rlimits[upper])
    spoints = collect(Float64, slimits[lower]:hs:slimits[upper])
    sum = 0
    for i in rpoints
        buffsum = 0
        for j in spoints
            buffsum += (f(i, j) * hr * hs)
        end
        sum += buffsum
    end
    return sum
end  # cellmethod

"""
    getgausspoints_2d(int_order::Int)

Return Gauss integration points for given integration order.

# Arguments
- `int_order`: integration order.
"""
function getgausspoints_2d(int_order::Int)
    points_coords = Array{Tuple{Float64, Float64}}(undef, int_order^2)
    if int_order == 2
        points_coords[1] = (1 / sqrt(3), 1 / sqrt(3))
        points_coords[2] = (-1 / sqrt(3), 1 / sqrt(3))
        points_coords[3] = (-1 / sqrt(3), -1 / sqrt(3))
        points_coords[4] = (1 / sqrt(3), -1 / sqrt(3))
        return points_coords
    elseif int_order == 3
        points_coords[1] = (sqrt(0.6), sqrt(0.6))
        points_coords[2] = (-sqrt(0.6), sqrt(0.6))
        points_coords[3] = (-sqrt(0.6), -sqrt(0.6))
        points_coords[4] = (sqrt(0.6), -sqrt(0.6))
        points_coords[5] = (0, sqrt(0.6))
        points_coords[6] = (-sqrt(0.6), 0)
        points_coords[7] = (0, -sqrt(0.6))
        points_coords[8] = (sqrt(0.6), 0)
        points_coords[9] = (0, 0)
        return points_coords
    else
        @error("Wrong integration order in getgausspoints_2d()")
        return nothing
    end
end  # getgausspoints_2d

"""
    gaussmethod(f::Function, int_order::Int)

Multiple integral calculation with Gauss method. f has to be depended on r and s only. 
So other functions have to be converted to such format before calling this method.
Interval of integration should be equal to [-1; 1].
Supported integration orders: 2.

# Arguments
- `f::Function`: function that needs to be integrated, need to be depended on 2 variables: 
``f = f(r, s)``;
- `int_order::Int`: integration order.
"""
function gaussmethod(f::Function, int_order::Int)
    rarray = Array  # Array of integration points by r coordinate
    sarray = Array  # Array of integration points by s coordinate
    weights = Array  # Array of integration weights
    if int_order == 2
        r = [-1 / sqrt(3), 1 / sqrt(3)]
        s = [-1 / sqrt(3), 1 / sqrt(3)]
        weights = [1, 1]
    else
        @error("Wrong integration order in gaussmethod()")
        return nothing
    end
    resultsum = 0
    for i in 1:int_order
        for j in 1:int_order
            resultsum += (weights[i] * weights[j] * f(r[i], s[j]))
        end
    end
    return resultsum
end  # gaussmethod

"""
    gaussmethod_matrix_1d(f::Function, int_order::Int)

Integrate matrix of functions depending on 1 variable with Gauss method.

# Arguments
- `f::Function`: function returning matrix of functions depending on 1 variable;
- `int_order::Int`: integration order.
"""
function gaussmethod_matrix_1d(f::Function, int_order::Int)
    x = Array{Float64}(undef, int_order)
    weights = Array{Float64}(undef, int_order)
    if int_order == 2
        # x = [-1 / sqrt(3), 1 / sqrt(3)]
        x = [-0.577350269189626, 0.577350269189626]
        weights = [1.0, 1.0]
    elseif int_order == 3
        x = [-0.774596669241483, 0, 0.774596669241483]
        weights = [0.555555555555556, 0.888888888888889, 0.555555555555556]
    elseif int_order == 4
        x = [-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053]
        weights = [0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454]
    else
        @error("Wrong integration order in gaussmethod_matrix_1d()")
        return nothing
    end
    n_rows = size(f(1))[1]
    if length(size(f(1))) == 1
        n_cols = 1
    else
        n_cols = size(f(1))[2]
    end
    result_matrix = Matrix{Float64}(undef, n_rows, n_cols)
    for k in 1:n_rows
        for l in 1:n_cols
            # Integrating [k, l] element of source matrix
            result_sum = 0
            for i in 1:int_order
                result_sum += (weights[i] * f(x[i])[k, l])
            end
            result_matrix[k, l] = result_sum
        end
    end
    return result_matrix
end  # gaussmethod_matrix_1d

"""
    gaussmethod_matrix_2d(f::Function, int_order::Int)

Integrate matrix of functions depending on 2 variables with Gauss method.

# Arguments
- `f::Function`: functions returning matrix of functions depending on 2 variables;
- `int_order::Int`: integration order.
"""
function gaussmethod_matrix_2d(f::Function, int_order::Int)
    r = Array{Float64}(undef, int_order)  # Array of integration points by r coordinate
    s = Array{Float64}(undef, int_order)  # Array of integration points by s coordinate
    weights = Array{Float64}(undef, int_order)  # Array of integration weights
    if int_order == 2
        r = [-1 / sqrt(3), 1 / sqrt(3)]
        s = [-1 / sqrt(3), 1 / sqrt(3)]
        weights = [1, 1]
    elseif int_order == 3
        r = [-0.774596669241483, 0, 0.774596669241483]
        s = [-0.774596669241483, 0, 0.774596669241483]
        weights = [0.555555555555556, 0.888888888888889, 0.555555555555556]
    elseif int_order == 4
        r = [-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053]
        s = [-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053]
        weights = [0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454]
    else
        @error("Wrong integration order in gaussmethod_matrix_2d()")
        return nothing
    end
    n_rows = size(f(1, 1))[1]
    if length(size(f(1, 1))) == 1
        n_cols = 1
    else
        n_cols = size(f(1, 1))[2]
    end
    result_matrix = Matrix{Float64}(undef, n_rows, n_cols)
    for k in 1:n_rows
        for l in 1:n_cols
            # Integrating [k, l] element of source matrix
            result_sum = 0
            for i in 1:int_order
                for j in 1:int_order
                    result_sum += (weights[i] * weights[j] * f(r[i], s[j])[k, l])
                end
            end
            result_matrix[k, l] = result_sum
        end
    end
    return result_matrix
end  # gaussmethod_matrix


"""
    gaussmethod_matrix_3d(f::Function, int_order::Int)

Integrate matrix of functions depending on 3 variables with Gauss method.

# Arguments
- `f::Function`: functions returning matrix of functions depending on 2 variables;
- `int_order::Int`: integration order.
"""
function gaussmethod_matrix_3d(f::Function, int_order::Int)
    r = Array{Float64}(undef, int_order)  # Array of integration points by r coordinate
    s = Array{Float64}(undef, int_order)  # Array of integration points by s coordinate
    t = Array{Float64}(undef, int_order)  # Array of integration points by t coordinate
    weights = Array{Float64}(undef, int_order)  # Array of integration weights

    # Defining integration points and weights
    if int_order == 2
        r = [-1 / sqrt(3), 1 / sqrt(3)]
        s = [-1 / sqrt(3), 1 / sqrt(3)]
        t = [-1 / sqrt(3), 1 / sqrt(3)]
        weights = [1, 1]
    elseif int_order == 3
        r = [-0.774596669241483, 0, 0.774596669241483]
        s = [-0.774596669241483, 0, 0.774596669241483]
        t = [-0.774596669241483, 0, 0.774596669241483]
        weights = [0.555555555555556, 0.888888888888889, 0.555555555555556]
    elseif int_order == 4
        r = [-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053]
        s = [-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053]
        t = [-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053]
        weights = [0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454]
    else
        @error("Wrong integration order in gaussmethod_matrix_3d()")
        return nothing
    end

    # Creating result matrix
    n_rows = size(f(1, 1, 1))[1]
    if length(size(f(1, 1, 1))) == 1
        n_cols = 1
    else
        n_cols = size(f(1, 1, 1))[2]
    end
    result_matrix = zeros(n_rows, n_cols)

    # Precalculating weights
    n_weights = size(weights)[1]
    weights_matr = Array{Float64}(undef, n_weights, n_weights, n_weights)
    for i in 1:int_order
        for j in 1:int_order
            for k in 1:int_order
                weights_matr[i, j, k] =  weights[i] * weights[j] * weights[k]
            end
        end
    end
    
    # Integrating
    for i in 1:int_order
        for j in 1:int_order
            for v in 1:int_order
                fmatrix = f(r[i], s[j], t[v])
                fmatrix .*= weights_matr[i, j, v]
                for k in 1:n_cols
                    for l in k:n_rows
                        result_matrix[l, k] += fmatrix[l, k]
                        if l != k
                            result_matrix[k, l] += fmatrix[l, k]
                        end
                    end
                end
            end
        end
    end

    return result_matrix
end

"""
    gaussmethod_matrix_2d_nonloc(f::Function, impactfunc::Function, impactdist::Number, 
        loc_to_glob::Function, x_coords::Array{Float64}, y_coords::Array{Float64} 
        int_order::Int)

Integrate matrix of functions depending on 2 variables with Gauss method twice taking
into account impact function.

`f` should depend on 4 variables: local r, s and impact r, s. Impact function is
assumed to be in Gauss form and should depend on 3 variables: normalization factor, impact
distance and current distance.

# Arguments
- `f::Function`: functions returning matrix of functions depending on 6 variables;
- `impactfunc::Function`: impact function;
- `impactdist::Number`: impact distance;
- `loc_to_glob::Function`: function which allows to convert local coordinates to global 
ones;
- `x_coords::Array{Float64}`: x-coords of each node in current element;
- `y_coords::Array{Float64}`: y-coords of each node in current element; 
- `int_order::Int`: integration order.
"""
function gaussmethod_matrix_2d_nonloc(f::Function, impactfunc::Function, impactdist::Number, 
    loc_to_glob::Function, x_source::Array{Float64}, y_source::Array{Float64},
    x_impact::Array{Float64}, y_impact::Array{Float64}, int_order::Int
)
    # Defining integration points and weights
    if int_order == 2
        r = [-1 / sqrt(3), 1 / sqrt(3)]
        s = [-1 / sqrt(3), 1 / sqrt(3)]
        weights = [1, 1]
    elseif int_order == 3
        r = [-0.774596669241483, 0, 0.774596669241483]
        s = [-0.774596669241483, 0, 0.774596669241483]
        weights = [0.555555555555556, 0.888888888888889, 0.555555555555556]
    elseif int_order == 4
        r = [-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053]
        s = [-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053]
        weights = [0.347854845137454, 0.652145154862546, 0.652145154862546, 
                    0.347854845137454]
    else
        @error("Wrong integration order in gaussmethod_matrix_2d_nonloc()")
        return nothing
    end

    # Creating result matrix
    n_rows = size(f(1, 1, 1, 1))[1]
    if length(size(f(1, 1, 1, 1))) == 1
        n_cols = 1
    else
        n_cols = size(f(1, 1, 1, 1))[2]
    end
    result_matrix = zeros(n_rows, n_cols)

    # Precalculating weights
    n_weights = size(weights)[1]
    weights_matr = Array{Float64}(undef, n_weights, n_weights)
    for i in 1:int_order
        for j in 1:int_order
            weights_matr[i, j] =  weights[i] * weights[j]
        end
    end

    # Normalization factor for impact function
    normfact = 1 / (2 * pi * impactdist^2)

    # Integrating
    firstpt_loc = (0, 0)  # Center of element (in r, s coordinates)
    firstpt_glob = loc_to_glob(firstpt_loc[1], firstpt_loc[2], x_source, y_source)
    for i_loc in 1:int_order
        for j_loc in 1:int_order
            for i_imp in 1:int_order
                for j_imp in 1:int_order
                    fmatrix = f(r[i_loc], s[j_loc], r[i_imp], s[j_imp])
                    fmatrix .*= weights_matr[i_loc, j_loc]
                    fmatrix .*= weights_matr[i_imp, j_imp]

                    # Calculating impact function
                    secondpt_loc = (r[i_imp], s[j_imp])
                    secondpt_glob = loc_to_glob(secondpt_loc[1], secondpt_loc[2], x_impact, 
                        y_impact)
                    vec = vecfrompoints_2d(firstpt_glob, secondpt_glob)
                    dist = veclength_3d(vec)
                    fmatrix .*= impactfunc(normfact, impactdist, dist)

                    for k in 1:n_cols
                        for l in k:n_rows
                            result_matrix[l, k] += fmatrix[l, k]
                            if l != k
                                result_matrix[k, l] += fmatrix[l, k]
                            end
                        end
                    end
                end
            end
        end
    end

    return result_matrix
end

"""
    gaussmethod_matrix_3d_nonloc(f::Function, impactfunc::Function, impactdist::Number, 
        int_order::Int)

Integrate matrix of functions depending on 3 variables with Gauss method twice taking
into account impact function.

`f` should depend on 6 variables: local r, s, t and impact r, s, t. Impact function is
assumed to be in Gauss form and should depend on 3 variables: normalization factor, impact
distance and current distance.

# Arguments
- `f::Function`: functions returning matrix of functions depending on 6 variables;
- `impactfunc::Function`: impact function;
- `impactdist::Number`: impact distance;
- `int_order::Int`: integration order.
"""
function gaussmethod_matrix_3d_nonloc(f::Function, impactfunc::Function, impactdist::Number, 
    int_order::Int
)
    # Defining integration points and weights
    if int_order == 2
        r = [-1 / sqrt(3), 1 / sqrt(3)]
        s = [-1 / sqrt(3), 1 / sqrt(3)]
        t = [-1 / sqrt(3), 1 / sqrt(3)]
        weights = [1, 1]
    elseif int_order == 3
        r = [-0.774596669241483, 0, 0.774596669241483]
        s = [-0.774596669241483, 0, 0.774596669241483]
        t = [-0.774596669241483, 0, 0.774596669241483]
        weights = [0.555555555555556, 0.888888888888889, 0.555555555555556]
    elseif int_order == 4
        r = [-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053]
        s = [-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053]
        t = [-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053]
        weights = [0.347854845137454, 0.652145154862546, 0.652145154862546, 
                    0.347854845137454]
    else
        @error("Wrong integration order in gaussmethod_matrix_3d_nonloc()")
        return nothing
    end

    # Creating result matrix
    n_rows = size(f(1, 1, 1, 1, 1, 1))[1]
    if length(size(f(1, 1, 1, 1, 1, 1))) == 1
        n_cols = 1
    else
        n_cols = size(f(1, 1, 1, 1, 1, 1))[2]
    end
    result_matrix = zeros(n_rows, n_cols)

    # Precalculating weights
    n_weights = size(weights)[1]
    weights_matr = Array{Float64}(undef, n_weights, n_weights, n_weights)
    for i in 1:int_order
        for j in 1:int_order
            for k in 1:int_order
                weights_matr[i, j, k] =  weights[i] * weights[j] * weights[k]
            end
        end
    end

    # Normalization factor for impact function
    normfact = 1 / (2 * pi * impactdist^2)

    # Integrating
    firstPt = (0, 0, 0)  # Center of element (in r, s, t coordinates)
    for i_loc in 1:int_order
        for j_loc in 1:int_order
            for v_loc in 1:int_order
                for i_imp in 1:int_order
                    for j_imp in 1:int_order
                        for v_imp in 1:int_order
                            fmatrix = f(r[i_loc], s[j_loc], t[v_loc], r[i_imp], s[j_imp], 
                                t[v_imp])
                            fmatrix .*= weights_matr[i_loc, j_loc, v_loc]
                            fmatrix .*= weights_matr[i_imp, j_imp, v_imp]

                            # Calculating impact function
                            secondPt = (r[i_imp], s[j_imp], t[v_imp])
                            vec = vecfrompoints_3d(firstPt, secondPt)
                            dist = veclength_3d(vec)
                            fmatrix .*= impactfunc(normfact, impactdist, dist)

                            for k in 1:n_cols
                                for l in k:n_rows
                                    result_matrix[l, k] += fmatrix[l, k]
                                    if l != k
                                        result_matrix[k, l] += fmatrix[l, k]
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    return result_matrix
end
