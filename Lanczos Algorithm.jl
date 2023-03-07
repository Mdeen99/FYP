using LinearAlgebra
using SparseArrays
function neworthogonalvector(V)
function lanczos(A, m::Int64=size(A)[1])
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
        else
            #getting a v that is orthogonal to all other v
            v=transpose(V[:,1:j-1]) \ zeros(j-1,1)
            println(V[:,1:j-1]*v==zeros(j-1,1))
            v= v/norm(v)
        end
        V[:,j] = v
        w = A*v
        T[j,j] = dot(w,v)
        w = w-T[j,j]*v-beta*V[:,j-1]
    end
    return 1
end
A= sprand(50,50,0.5)
b= lanczos(A)
#display(A)
#println("Julia's standard eigenvalue method")
#display(maximum(abs.(eigvals(Array(A)))))
#eigenvalue, eigenvector =powerMethod(A, 0.000001)
#println("powermethod I wrote")
#display(eigenvalue)
#display(eigenvector)

#println("test")