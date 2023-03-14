using LinearAlgebra
using SparseArrays
include("Power method.jl")

function lanczos(A, m::Int64=size(A)[1])
    #=
        Lanczos algorithm to decompose A into VTV^*
        Note that this function assumes that A is a real matrix
        Not numerically stable in the slightest
    =#
    n=size(A)[1]
    v = [1; zeros(n-1,1)]
    V = zeros(n,m)
    V[:,1] = v
    T=spzeros(m,m)
    w = A*v
    T[1,1] = dot(w,v)
    w = w-T[1,1]*v
    for j=2:m
        beta = norm(w)
        T[j,j-1]=beta
        T[j-1,j]=beta
        if beta != 0
            v=(1/beta)*w
            for k=1:j
                v=v - dot(V[:,j],v)*V[:,j]
            end
            v=v/norm(v)
        else
            #getting a v that is orthogonal to all other v
            #Beta is never 0
            v=transpose(V[:,1:j-1]) \ zeros(j-1,1)
            display(V[:,1:j-1]*v===zeros(j-1,1))
            v= v/norm(v)
        end
        V[:,j] = v
        w = A*v
        T[j,j] = dot(w,v)
        w = w-T[j,j]*v-beta*V[:,j-1]
    end
    return T, V
end


A= sprand(1000,1000,0.5)
T,V= lanczos(A)
eigenvalue, eigenvector =powerMethod(A, 0.000001)
println("powermethod")
display(eigenvalue)
println("Lanczos with power")
eigenvalue, eigenvector =powerMethod(T, 0.000001)
display(eigenvalue)

