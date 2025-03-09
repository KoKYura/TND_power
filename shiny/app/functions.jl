using Random, Distributions, DataFrames, GLM

logit(x) = log(x/(1-x))
invlogit(x) = 1/(1+exp(-x))

function get_prop_g(prob_case::Float64, prop_neg_g::Vector{Float64}, ve_g::Vector{Float64})
    total_num = Matrix{Float64}(undef, length(prop_neg_g), 2)
    prob_pos = invlogit.(logit(prob_case) .+ log.(1 .- ve_g))
    total_num[:, 1] = prop_neg_g .* prob_pos ./ (1 .- prob_pos)
    total_num[:, 2] = prop_neg_g
    prop_g = (total_num[:, 1] .+ total_num[:, 2]) ./ sum(total_num)
    return prop_g
end

# Based on codes of Endo et al. 
# https://github.com/akira-endo/TND-biascorrection

function get_misclassified(y_::Vector{Int}, sens::Float64, spec::Float64)
    true_TD = sum(y_)
    true_ND = length(y_) - true_TD
    
    sample_pos = [
        rand(Binomial(true_TD, sens)),
        rand(Binomial(true_ND, 1 - spec))
    ]
    
    obs_TD = vcat(
        findall(x -> x == 1, y_)[1:sample_pos[1]],
        findall(x -> x == 0, y_)[1:sample_pos[2]]
    )
    
    y_obs_ = zeros(Int, length(y_))
    y_obs_[obs_TD] .= 1
    
    return y_obs_
end

function simulate_data(n::Int, prob_case::Float64, prop_g::Vector{Float64}, ve_g::Vector{Float64}, sens::Float64, spec::Float64)

    n_g = [round(Int, n * p) for p in prop_g]
    
    df_sim = DataFrame(group=String[], y_true=Int[], y_obs=Int[])
    
    for g in 1:7
        group = "group$g"
        prob = invlogit(logit(prob_case) + log(1 - ve_g[g]))
        y_true = rand(Binomial(1, prob), n_g[g])
        y_obs = get_misclassified(y_true, sens, spec)
        
        df_sim_g = DataFrame(group=fill(group, n_g[g]), y_true=y_true, y_obs=y_obs)
        append!(df_sim, df_sim_g)
    end
    
    return df_sim
end

# Based on codes of Endo et al. 
# https://github.com/akira-endo/TND-biascorrection

n_iter = 100
function get_ve_corrected(df_sim::DataFrame, sens::Float64, spec::Float64)
    res_bs = DataFrame(group=String[], est_bs=Float64[])
    
    for iter in 1:n_iter

        indices = sample(1:nrow(df_sim), nrow(df_sim), replace=true)
        df_sim_bs = df_sim[indices, :]
        
        glm_bs = glm(@formula(y_obs ~ group), df_sim_bs, Binomial(), LogitLink())
        
        D_pos = filter(row -> row[:y_obs] == 1, df_sim_bs)
        D_neg = filter(row -> row[:y_obs] == 0, df_sim_bs)
        phat_pos = predict(glm_bs, D_pos)
        phat_neg = predict(glm_bs, D_neg)
        pflip_pos = clamp.((sens * (1 .- phat_pos) ./ phat_pos .- (1 .- sens)) .* ((1 .- spec) ./ (sens .+ spec .- 1)), 0, 1)
        pflip_neg = clamp.((spec * phat_neg ./ (1 .- phat_neg) .- (1 .- spec)) .* ((1 .- sens) ./ (sens .+ spec .- 1)), 0, 1)
        
        n_pos = nrow(D_pos)
        n_neg = nrow(D_neg)
        idflip_pos = rand(n_pos) .< pflip_pos
        idflip_neg = rand(n_neg) .< pflip_neg
        
        copy_d_pos = deepcopy(D_pos)
        copy_d_neg = deepcopy(D_neg)
        copy_d_pos.y_cor = copy_d_pos.y_obs
        copy_d_neg.y_cor = copy_d_neg.y_obs
        
        copy_d_pos[idflip_pos, :y_cor] .= 0
        copy_d_neg[idflip_neg, :y_cor] .= 1
        
        data = vcat(copy_d_pos, copy_d_neg)
        
        glm_cor = glm(@formula(y_cor ~ group), data, Binomial(), LogitLink())
        coef_table = coeftable(glm_cor)  
        coef_df = DataFrame(coef_table)
        
        for g in 2:7
            group_name = "group$g"
            est_bs = 1 - exp(coef_df[g,2]) 
            push!(res_bs, (group_name, est_bs))
        end
    end
    
        df_est_med = combine(groupby(res_bs, :group), 
                            :est_bs => median => :est_bs_med)
        df_est_l = combine(groupby(res_bs, :group), 
                            :est_bs => x -> quantile(x, 0.025))
        rename!(df_est_l, :est_bs_function => :est_bs_l)
        df_est_u = combine(groupby(res_bs, :group), 
                            :est_bs => x -> quantile(x, 0.975))
        rename!(df_est_u, :est_bs_function => :est_bs_u)
        
        df_est_ = leftjoin(df_est_med, df_est_l, on = :group)
        df_est_bs = leftjoin(df_est_, df_est_u, on = :group)
        
        df_est_bs.sdif_bs = (df_est_bs.est_bs_l .> 0) .& (df_est_bs.est_bs_u .> 0)
    
    return df_est_bs
end

function get_est(n::Int, prob_case::Float64, prop_neg_g::Vector{Float64}, ve_g::Vector{Float64}, sens::Float64, spec::Float64, null)
    
    prop_g = get_prop_g(prob_case, prop_neg_g, ve_g)
    df_sim = simulate_data(n, prob_case, prop_g, ve_g, sens, spec)
    
    glm_true = glm(@formula(y_true ~ group), df_sim, Binomial(), LogitLink())
    glm_obs = glm(@formula(y_obs ~ group), df_sim, Binomial(), LogitLink())

    coef_true = DataFrame(coeftable(glm_true))
    coef_obs = DataFrame(coeftable(glm_obs))
    
    df_est_uncorrected = DataFrame()
    for g in 2:7
        group_name = "group$(g)"
        est_true = 1 - exp(coef_true[g,2])
        l_true = 1 - exp(coef_true[g,7])
        u_true = 1 - exp(coef_true[g,6])

        est_obs = 1 - exp(coef_obs[g,2])
        l_obs = 1 - exp(coef_obs[g,7])
        u_obs = 1 - exp(coef_obs[g,6])
        
        df_est_uncorrected_ = DataFrame(n = n,
                                        group = group_name,
                                        est_true = est_true,
                                        l_true = l_true,
                                        u_true = u_true,
                                        est_obs = est_obs,
                                        l_obs = l_obs,
                                        u_obs = u_obs)
        append!(df_est_uncorrected, df_est_uncorrected_)
    end
   
    df_est_corrected = get_ve_corrected(df_sim, sens, spec)
    
    df_est = leftjoin(df_est_uncorrected, df_est_corrected, on = :group)
    df_est.null .= null
    
    return df_est
end