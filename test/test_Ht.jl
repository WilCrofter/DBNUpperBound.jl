begin # to localize variable scope

    # Determine absolute path to Python sample output file
    local fpath=joinpath(dirname(@__FILE__),"sample_output_Ht_real.txt") 
    
    # Python sample output from KM; [1] is data, [2] is header
    local tmp = readdlm(fpath,',',header=true)[1]

    # Set random seed for reproducibility.
    # seed = rand(UInt32) # (when test was written)
    # 0xdcfdf996
    srand(0xdcfdf996) 

    # 100 random rows of Python output. There are 7006 total.
    # Sampling to avoid long test runtimes
    local idx = rand(1:size(tmp,1),100)

    # Generate Julia results for comparison; a is a vector
    # of 2-tuples [1] being Ht, [2] being estimated error
    # of quadrature.
    local a=[Ht(tmp[i,1],tmp[i,2],upper_limit=10.0) for i in idx]

    # Function to compare Julia and Python results
    local function f(a,x)
        # Are a[1] and x[3] within estimated errors,
        # a[2] and x[4], of one another?
        return abs(a[1]-x[3]) ≤ a[2]+x[4]
    end

    @test all([f(a[i], tmp[idx[i],:]) for i in 1:length(idx)])

end