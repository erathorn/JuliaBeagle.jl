function maxcbs(S)::Int
    floor(512/S)
end

function threadblocks(R, C, CBS)::Int
    R*ceil(C/CBS)
end


function mygemmturbo!(C, A, B)
    @tturbo for m ∈ axes(A, 1), n ∈ axes(B, 2)
        Cmn = zero(eltype(C))
        for k ∈ axes(A, 2)
            Cmn += A[m, k] * B[k, n]
        end
        C[m, n] = Cmn
    end
end

function mygemmturbo_diag!(C, A, B)
    @tturbo for m ∈ axes(A, 1), n ∈ axes(B, 1)
        Cmn = zero(eltype(C))
        for k ∈ axes(A, 2)
            t = m == n ? exp(B[n]) : zero(eltype(C))
            Cmn += A[m, k] * t
        end
        C[m, n] = Cmn
    end
end