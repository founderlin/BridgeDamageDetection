module lhs # module of Latin Hybercube Sampling

using Distributions
export LHS_all

    function LHS_one(rn)
        v0 = zeros(rn, 1);
        dtn = 1/rn # ru -> interval number
        for i = 1:rn
            dtn = 1/rn
            pp = [ item for item = 0:dtn:1-dtn ]
            pa = pp[i]
            pb = pa + dtn
            v = rand(Uniform(pa, pb))
            v0[i] = v0[i] + v
        end
        return v0
    end

    function LHS_all(rn, n)
        v1 = zeros(rn, n)
        for i = 1:n
            v0 = LHS_one(rn)
            v1[:,i] = v1[:,i] + v0
        end
        return v1
    end
end
