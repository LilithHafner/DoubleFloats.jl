@inline function add_fpdd_dd(x::T, y::Tuple{T,T}) where {T<:IEEEFloat}
    yhi, ylo = y
    yhi, ylo = two_sum(x, yhi, ylo)
    return yhi, ylo
end

@inline function sub_fpdd_dd(x::T, y::Tuple{T,T}) where {T<:IEEEFloat}
    yhi, ylo = y
    yhi, ylo = two_diff(x, yhi, ylo)
    return yhi, ylo
end

#=
# Algorithm 12 from Tight and rigourous error bounds.  relative error <= 5u²
modified to handle +/- Inf properly
=#
@inline function mul_fpdd_dd(x::T, y::Tuple{T,T}) where T<:IEEEFloat
    mul_ddfp_dd(y, x)
end

#=
function mul_fpdd_dd(x::T, y::Tuple{T,T}) where T<:IEEEFloat
    hi, lo = y
    hihi, hilo = two_prod(x, hi)
    lohi, lolo = two_prod(x, lo)
    hi, lo = two_sum(hihi, hilo, lohi)
    isinf(hihi) ? (hihi, NaN) : (hi, lo)
end
=#	   
=#
@inline function mul_fpdd_dd(x::T, y::Tuple{T,T}) where T<:IEEEFloat
    yhi, ylo = x
    hihi, hilo = two_prod(x, yhi)
    t = x * ylo
    t = t + hilo
    hi, lo = two_hilo_sum(hihi, t)
    isinf(hi) ? (hi, NaN) : isinf(hihi) ? (hihi, NaN) : (hi, lo)
end

@inline function dvi_fpdd_dd(x::T, y::Tuple{T,T}) where {T<:IEEEFloat}
    yhi, ylo = y    
    hixy = x / yhi
    uh, ul = two_prod(hixy, y)
    lo = ((((xhi - uh) - ul) + xlo))/y
    hi,lo = two_hilo_sum(hixy, lo)
    isinf(hi) || isinf(hixy) ? (hi, NaN) : (hi, lo)
end

@inline function dvi_fpdd_dd_fast(x::T, y::Tuple{T,T}) where {T<:IEEEFloat}
    yinv = inv_dd_dd_fast(y)
    zhi, zlo = mul_fpdd_dd(x, yinv)
    return zhi, zlo
end
