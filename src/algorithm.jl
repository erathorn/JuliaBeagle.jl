cbs = 16
sbs = 4

# Number of states (depends on model)
S = 4
# number of columns
C = 100
# Number of Rate categories
R = 4

function parallel_transition_prob(E::Matrix, D::Vector, EI::Matrix, rates::Vector, mu::Float64, t::Vector)
    C = similar(E)
    out = Array{Float64, 3}(undef, (length(t)*length(D), length(t)*length(D),length(r)))
    for r in rates
        mygemmturbo!(C, E, diagm(exp.(D*t*r*mu)))
        mygemmturbo!(out[:, :, r], C, EI)
    end
    out
end


function felsenstein(R, S, Data, mu, E, D, EI, rates, tree)
    t = get_branchlength_vector(tree)
    pmat = parallel_transition_prob(E, D, EI, rates, mu, t)

end

function recurser_ll(node, R, S, Data, pmat)
    
    ll = 0.0
    if node.nchild > 0
        lchild = node.children[1]
        rchild = node.children[2]
        tmp_data = Array{Float64, 3}(undef, (lenght(R), CBS, S))
        for rate in 1:R, c in 1:CBS, s in 1:s
            sum_product!(tmp_data, lc_Data, rc_Data, pmat, rate, c, s)
        end
        

    end
end

function parallel_partial_likelihoods(R::Int, C::Int, CBS::Int, S::Int)
    ThreadBlocks = threadblocks(R, C, CBS)
    # ToDo: actual transition matrix
    large_t_mat = zeros(4,4,4,4)
    for threadblock in 1:ThreadBlocks
        rate = floor((threadblock-1)/ceil(C/CBS))
        cblock = (threadblock-1)%ceil(C/CBS)
        c_t_mat = large_t_mat[rate, :, :]
        #lchild = child_recursion()
        #rchild = child_recursion()
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
    tmp += lchild[r, c, s] * transprobs[1, s, r]
    tmp += lchild[r, c, s] * transprobs[2, s, r]
    tmp += lchild[r, c, s] * transprobs[3, s, r]
    tmp += lchild[r, c, s] * transprobs[4, s, r]
    
    tmp2 += rchild[r, c, s] * transprobs[1, s, r]
    tmp2 += rchild[r, c, s] * transprobs[2, s, r]
    tmp2 += rchild[r, c, s] * transprobs[3, s, r]
    tmp2 += rchild[r, c, s] * transprobs[4, s, r]
    parent[r, c, s] = tmp * tmp2
    nothing
end

