# Additional functions to the SymPy package to compute partial derivatives of 2D functions
function partialDiff1(f::Function, n::Int64=1)
    @vars x y
    diff1 = sympy.diff(f(x, y), x, n)
    return lambdify(diff1, (x, y))
end

function partialDiff2(f::Function, n::Int64=1)
    @vars x y
    diff2 = sympy.diff(f(x, y), y, n)
    return lambdify(diff2, (x, y))
end

# Utility function that transform any constant to a function or does nothing if it is already a function
function transformToFunction(constant)
    # If constant is not a function, build a constant function
    if !(typeof(constant) <: Function)
        constant_function(V) = Float64(constant)
        return constant_function
    # Otherwise just return constant
    else
        return constant
    end
end

# Bisection function
function bisection(f::Function, a::Union{Float64, Int64}, b::Union{Float64, Int64}; tol=1e-6::Float64, max_iter=100::Int64)
    # If no passing through 0, throw error
    if f(a) * f(b) > 0
        error("Function must have different signs at interval endpoints!")
    end
    
    # Start the bisection loop
    iter = 0
    while (b - a) / 2 > tol && iter < max_iter
        c = (a + b) / 2
        if f(c) == 0
            return c
        elseif f(c) * f(a) < 0
            b = c
        else
            a = c
        end
        iter += 1
    end
    
    return (a + b) / 2
end
