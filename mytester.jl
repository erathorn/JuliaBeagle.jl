using Pkg
Pkg.activate(".")

using JuliaBeagle
using LinearAlgebra
using LoopVectorization
X = reshape(1:9, 3,3)

y1 = 1:3
y = diagm(y1)
C = zeros(3,3)
@benchmark mygemmturbo_diag!($C, $X, $y1)
@benchmark $X * $y1


function mygemmturbo_diag!(C, A, B)
    @tturbo for m ∈ axes(A, 1), n ∈ axes(B, 1)
        Cmn = zero(eltype(C))
        for k ∈ axes(A, 2)
            t = k == n ? B[n] : 0
            Cmn += A[m, k] * t
        end
        C[m, n] = Cmn
    end
end