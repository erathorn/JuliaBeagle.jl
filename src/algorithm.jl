cbs = 16
sbs = 4

# Number of states (depends on model)
S = 4
# number of columns
C = 100
# Number of Rate categories
R = 4

function parallel_transition_prob(E::Matrix, D::Vector, EI::Matrix, rates::Vector, t::Vector)
    C = similar(E)
    out = Array{Float64, 3}(undef, (length(r), length(t)*length(D), length(t)*length(D)))
    for r in rates
        mygemmturbo!(C, E, diagm(exp.(D*t*r)))
        mygemmturbo!(out[r, :, :], C, EI)
    end
    out
end


function parallel_partial_likelihoods(R::Int, C::Int, CBS::Int, S::Int)
    ThreadBlocks = threadblocks(R, C, CBS)
    # ToDo: actual transition matrix
    large_t_mat = zeros(4,4,4,4)
    for threadblock in 1:ThreadBlocks
        rate = floor((threadblock-1)/ceil(C/CBS))
        cblock = (threadblock-1)%ceil(C/CBS)
        c_t_mat = large_t_mat[rate, :, :]
        lchild = child_recursion()
        rchild = child_recursion()
        parent = similar(lchild)
        for s in 1:S
            for c in 1:CBS
                sum_product!(parent, lchild, rchild, c_t_mat, rate, c, s)
            end
        end
    end
end

"""
    Calculate the partial likelihoods (in place) for a node given the values of its two children 
"""
function sum_product!(parent, lchild, rchild, transprobs, r::I, c::I, s::I)::Nothing where I<:Int
    tmp = zero(eltypt(parent))
    tmp2 = zero(eltypt(parent))
    tmp += lchild[r, c, s] * transprobs[r, 1, s]
    tmp += lchild[r, c, s] * transprobs[r, 2, s]
    tmp += lchild[r, c, s] * transprobs[r, 3, s]
    tmp += lchild[r, c, s] * transprobs[r, 4, s]
    
    tmp2 += rchild[r, c, s] * transprobs[r, 1, s]
    tmp2 += rchild[r, c, s] * transprobs[r, 2, s]
    tmp2 += rchild[r, c, s] * transprobs[r, 3, s]
    tmp2 += rchild[r, c, s] * transprobs[r, 4, s]
    parent[r, c, s] = tmp * tmp2
    nothing
end

