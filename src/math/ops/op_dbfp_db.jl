@inline function add_dbfp_db(x::DoubleFloat{T}, y::T) where {T<:IEEEFloat}
    isfinite(HI(x) + HI(y)) && return DoubleFloat{T}(add_dbfp_dd(HILO(x), y))
    add_dbfp_db_nonfinite(x, y)
end

@inline function sub_dbfp_db(x::DoubleFloat{T}, y::T) where {T<:IEEEFloat}
    isfinite(HI(x) - HI(y)) && return DoubleFloat{T}(sub_dbfp_dd(HILO(x), y))
    sub_dbfp_db_nonfinite(x, y)
end

@inline function mul_dbfp_db(x::DoubleFloat{T}, y::T) where {T<:IEEEFloat}
    isfinite(HI(x) + HI(y)) && return DoubleFloat{T}(mul_dbfp_dd(HILO(x), y))
    mul_dbfp_db_nonfinite(x, y)
end

@inline function dvi_dbfp_db(x::DoubleFloat{T}, y::T) where {T<:IEEEFloat}
    isfinite(HI(x) / HI(y)) && return DoubleFloat{T}(dvi_dbfp_dd(HILO(x), y))
    dvi_dbfp_db_nonfinite(x, y)
end

@inline function add_dbfp_db_nonfinite(x::DoubleFloat{T}, y::T) where {T<:IEEEFloat}
    z = HI(x) + y
    return DoubleFloat{T}(z, T(NaN))
end

@inline function sub_dbfp_db_nonfinite(x::DoubleFloat{T}, y::T) where {T<:IEEEFloat}
    z = HI(x) - y
    return DoubleFloat{T}(z, T(NaN))
end

@inline function mul_dbfp_db_nonfinite(x::DoubleFloat{T}, y::T) where {T<:IEEEFloat}
    z = HI(x) * y
    return DoubleFloat{T}(z, T(NaN))
end

@inline function dvi_dbfp_db_nonfinite(x::DoubleFloat{T}, y::T) where {T<:IEEEFloat}
    z = HI(x) / y
    return DoubleFloat{T}(z, T(NaN))
end
