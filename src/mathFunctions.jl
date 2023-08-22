# Additional functions to the SymPy package to compute partial derivatives of 2D functions
function partialDiff1(f::Function, n::Int64=1)
    @vars x y
    sympy.diff(f(x, y), x, n)
end

function partialDiff2(f::Function, n::Int64=1)
    @vars x y
    sympy.diff(f(x, y), y, n)
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