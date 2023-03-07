using LinearAlgebra
using Printf
using SparseArrays

function powerMethod(A, esp::Float64=0.0001, b0=nothing)
    #=
        Uses the power method for finding the largest eigenvalue
        A needs to be a square array and can be sparse.
        To get the eigenvalue from the eigenvector it uses Rayleigh quotient
        Occasionaly b0 would not converge so I added a counter i so it will stop eventually
        If the max eigenvalue is 0, this will return NaN
    =#
    err = 1
    i=1
    n = size(A)[1]
    if b0===nothing || length(B)!=n
        b0=ones(n,1)
    end
    while err > esp && i<=1000
        Ab=A*b0
        b1 = Ab / sqrt(dot(Ab,Ab))
        err = sqrt(dot(b1-b0,b1-b0))
        b0=b1
        i=i+1
    end
    #Rayleigh's quotient returns a 1x1 matrix
    eigenvalue = abs.(transpose(b0)*A*b0/(dot(b0,b0)))
    return eigenvalue[1], b0
end
#=
A= sprand(5000,5000,0.5)
#display(A)
println("Julia's standard eigenvalue method")
display(maximum(abs.(eigvals(Array(A)))))
eigenvalue, eigenvector =powerMethod(A, 0.000001)
println("powermethod I wrote")
display(eigenvalue)
#display(eigenvector)
=#