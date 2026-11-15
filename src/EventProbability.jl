"""
    EventProbability{T, M<:AbstractMatrix{T}, A<:AbstractArray{T,3}, V<:AbstractVector{T}}

Contiguous column-major storage for per-event unbinned PID probabilities:
- `yields::M`: Expected photon yields of size `(ntracks, npids)`
- `sig_prob::A`: Normalized signal PDF values of size `(nhits, ntracks, npids)` (`hit` is stride-1)
- `bkg_prob::V`: Background PDF values of size `(nhits,)`
"""
struct EventProbability{
    T<:AbstractFloat,
    M<:AbstractMatrix{T},
    A<:AbstractArray{T,3},
    V<:AbstractVector{T},
}
    ntracks::Int
    nhits::Int
    npids::Int
    yields::M
    sig_prob::A
    bkg_prob::V
end

function EventProbability(
    ntracks::Int,
    nhits::Int,
    npids::Int,
    ::Type{T} = Float64,
) where {T<:AbstractFloat}
    yields = zeros(T, ntracks, npids)
    sig_prob = zeros(T, nhits, ntracks, npids)
    bkg_prob = zeros(T, nhits)
    return EventProbability(ntracks, nhits, npids, yields, sig_prob, bkg_prob)
end

"""
    log_sum(prob::EventProbability{T}, ids::AbstractVector{<:Integer})::T

Zero-allocation evaluation of the unbinned extended log-likelihood for a given vector of
1-based hypothesis indices `ids` (one per track).
"""
@inline function log_sum(
    prob::EventProbability{T},
    ids::AbstractVector{<:Integer},
)::T where {T<:AbstractFloat}
    ntracks = prob.ntracks
    nhits = prob.nhits
    if nhits == 0
        return zero(T)
    end

    inv_nhits = one(T) / T(nhits)
    nsum = zero(T)
    @inbounds for t = 1:ntracks
        nsum += prob.yields[t, ids[t]]
    end
    bg_frac = max(T(nhits) - nsum, one(T)) * inv_nhits

    total_ll = zero(T)
    @inbounds for i = 1:nhits
        p = prob.bkg_prob[i] * bg_frac
        for t = 1:ntracks
            h = ids[t]
            fr = prob.yields[t, h] * inv_nhits
            p = muladd(fr, prob.sig_prob[i, t, h], p)
        end
        total_ll += log(max(p, floatmin(T)))
    end

    return total_ll
end

@inline function _evaluate_coordinate_sweep!(
    log_probs::AbstractMatrix{T},
    event_prob::EventProbability{T},
    best_id::AbstractVector{Int},
) where {T<:AbstractFloat}
    ntracks = event_prob.ntracks
    npids = event_prob.npids
    @inbounds for t = 1:ntracks
        last_id = best_id[t]
        for h = 1:npids
            best_id[t] = h
            log_probs[t, h] = log_sum(event_prob, best_id)
        end
        best_id[t] = last_id
    end
    return nothing
end

"""
    solve_unbinned_likelihoods!(
        log_probs::AbstractMatrix{T},
        event_prob::EventProbability{T},
        particles::AbstractVector{Particle},
        hypotheses::AbstractVector{<:Integer};
        n_iter::Int = 3,
        cheated::Bool = false,
        best_id_buf::Union{Vector{Int},Nothing} = nothing,
    )

Runs coordinate-ascent (or cheated MC truth) multi-track log-likelihood evaluation on
`event_prob`, writing the per-track, per-hypothesis log-likelihoods into `log_probs[t, h]`.
"""
function solve_unbinned_likelihoods!(
    log_probs::AbstractMatrix{T},
    event_prob::EventProbability{T},
    particles::AbstractVector{Particle},
    hypotheses::AbstractVector{<:Integer};
    n_iter::Int = 3,
    cheated::Bool = false,
    best_id_buf::Union{Vector{Int},Nothing} = nothing,
) where {T<:AbstractFloat}
    ntracks = event_prob.ntracks
    npids = event_prob.npids
    fill!(log_probs, zero(T))
    if ntracks == 0 || npids == 0
        return log_probs
    end

    if best_id_buf === nothing
        best_id = ones(Int, ntracks)
    else
        length(best_id_buf) < ntracks && resize!(best_id_buf, ntracks)
        best_id = best_id_buf
        @inbounds for t = 1:ntracks
            best_id[t] = 1
        end
    end

    if cheated
        @inbounds for t = 1:ntracks
            true_pid = abs(particles[t].pid)
            idx = findfirst(==(true_pid), hypotheses)
            best_id[t] = idx === nothing ? 1 : idx
        end
        _evaluate_coordinate_sweep!(log_probs, event_prob, best_id)
    else
        @inbounds for iter = 1:n_iter
            if iter > 1
                for t = 1:ntracks
                    best_h = 1
                    best_val = log_probs[t, 1]
                    for h = 2:npids
                        if log_probs[t, h] > best_val
                            best_val = log_probs[t, h]
                            best_h = h
                        end
                    end
                    best_id[t] = best_h
                end
            end
            _evaluate_coordinate_sweep!(log_probs, event_prob, best_id)
        end
    end

    return log_probs
end
