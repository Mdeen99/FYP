using LinearAlgebra
using Printf
using SparseArrays
using Statistics

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


#Playing around as I noticed that the size and sparsity of matix under the sprand function will have simmilar max eigenvalues
#eigenmax is approx n*spartisty/2
#=
means=zeros(9,6)
S=[0.1, 0.2 , 0.3, 0.4 ,0.5, 0.6, 0.7, 0.8, 0.9]
N=[10,50,100,500,1000,5000]
for s=1:size(S)[1]
    for n=1:size(N)[1]
        eigenmax=zeros(50)
        for i=1:size(eigenmax)[1]
            B=sprand(N[n],N[n],S[s])
            eigenmax[i],w=powerMethod(B)
        end
        means[s,n]=mean(eigenmax)
        println("n:",N[n]," s: ",S[s])
        display(means[s,n])
    end
end
=#
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