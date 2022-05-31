# more comprehensive tests of corner and extremal cases

@test iszero( (2.0 * Double64(pi)) - (Double64(pi) * 2.0) )
@test iszero( (0.5 * Double64(pi)) - (Double64(pi) * 0.5) )

@test iszero( (Double64(2.0) * Double64(pi)) - 
              (Double64(pi) * Double64(2.0)) )
@test iszero( (Double64(0.5) * Double64(pi)) - 
              (Double64(pi) * Double64(0.5)) )

FT = Float64
DT = DoubleFloat{FT}

dnan, fnan = DT(NaN), FT(NaN)
dinf, finf = DT(Inf), FT(Inf)

d0, f0 = zero(DT), zero(FT)
d1, f1 = one(DT), one(FT)

dmin, fmin = floatmin(DT), HI(floatmin(DT))
dmax, fmax = floatmax(DT), HI(floatmax(DT))
dmin₊, fmin₊ = dmin * (1025/1024), HI(dmin * (1025/1024))
dmax₋, fmax₋ = dmax * (1023/1024), HI(dmax * (1023/1024))

sqrt_dmin, sqrt_fmin = sqrt(dmin), HI(sqrt(dmin))
sqrt_dmax, sqrt_fmax = sqrt(dmax), HI(sqrt(dmax))
cbrt_dmin, cbrt_fmin = cbrt(dmin), HI(cbrt(dmin))
cbrt_dmax, cbrt_fmax = cbrt(dmax), HI(cbrt(dmax))
cbrt2_dmin, cbrt2_fmin = cbrt_dmin^2, HI(cbrt_dmin^2)
cbrt2_dmax, cbrt2_fmax = cbrt_dmax^2, HI(cbrt_dmax^2)

sqrt_dmin₊, sqrt_fmin₊ = sqrt(dmin₊), HI(sqrt(dmin₊))
sqrt_dmax₋, sqrt_fmax₋ = sqrt(dmax₋), HI(sqrt(dmax₋))
cbrt_dmin₊, cbrt_fmin₊ = cbrt(dmin₊), HI(cbrt(dmin₊))
cbrt_dmax₋, cbrt_fmax₋ = sqrt(dmax₋), HI(sqrt(dmax₋))
cbrt2_dmin₊, cbrt2_fmin₊ = cbrt_dmin₊^2, HI(cbrt_dmin₊^2) 
cbrt2_dmax₋, cbrt2_fmax₋ = cbrt_dmax₋^2, HI(cbrt_dmax₋^2)

dvals = (d1, 
        dmin, dmax, dmin₊, dmax₋, 
        sqrt_dmin, sqrt_dmax, sqrt_dmin₊, sqrt_dmax₋, 
        cbrt_dmin, cbrt_dmax, cbrt_dmin₊, cbrt_dmax₋,
        cbrt2_dmin, cbrt2_dmax, cbrt2_dmin₊, cbrt2_dmax₋,)

fvals = (f1, 
        dmin, fmax, fmin₊, fmax₋, 
        sqrt_fmin, sqrt_fmax, sqrt_fmin₊, sqrt_fmax₋, 
        cbrt_fmin, cbrt_fmax, cbrt_fmin₊, cbrt_fmax₋,
        cbrt2_fmin, cbrt2_fmax, cbrt2_fmin₊, cbrt2_fmax₋,)

dcomb = combinations(dvals, 2)
fcomb = combinations(fvals, 2)

dfvals = zip(dcomb, fcomb)

@testset "more addition" begin
  for (d,f) in dfvals
      @test isapprox(+(d...), +(f...))
  end
end

@testset "more subtraction" begin
  for (d,f) in dfvals
      @test isapprox(-(d...), -(f...))
  end
end

@testset "more multiplication" begin
  for (d,f) in dfvals
      @test isapprox(*(d...), *(f...))
  end
end

#=
@testset "more division" begin
  for (d,f) in dfvals
      if !isapprox(/(d...), /(f...))
         println("<<<<<<<<<<<<<<<")
         println((d, f))
         println(">>>>>>>>>>>>>>>")
      end
      @test isapprox(/(d...), /(f...))
  end
end
=#
