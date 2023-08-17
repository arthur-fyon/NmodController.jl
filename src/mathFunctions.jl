# Additional functions to the SymPy package to compute partial derivatives of 2D functions
function partialDiff1(f::Function, n::Int64=1)
    @vars x y
    sympy.diff(f(x, y), x, n)
end

function partialDiff2(f::Function, n::Int64=1)
    @vars x y
    sympy.diff(f(x, y), y, n)
end