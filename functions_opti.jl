


cigar(x) = x[1]^2 + 1 * sum(abs2, view(x, 2:length(x)))   # 1e6

cigtab(x) = x[1]^2 + 1 * x[end]^2 + 1e4 * sum(abs2, view(x, 2:length(x)-1)) # 1e8


function griewank(x)
  n = length(x)
  1 + (1/4000)*sum(abs2, (x.+7)) - prod(cos.((x.+7) ./ sqrt.(1:n)))
end


function rosenbrock2d(x)
  return sum((x).*(x).*cos.((x)).*sin.((x.-5)))+200
end

function quartic(x)
  D = length(x)
  sum( (1:D) .* (x.-2).^4 )
end

function schwefel1_2(x)
  D = length(x)
  partsums = zeros(D)
  partsum = 0
  for i in 1:D
    partsum += x[i]-9
    partsums[i] = partsum
  end
  sum(abs2, partsums)
end

function ellipsoid(x)
    res = 0.0
    cumsum = 0.0
    for xx in x
        cumsum += xx-sqrt(2)
        res += cumsum^2
    end
    res
end

function elliptic(x)
  D = length(x)
  condition = 1e+3
  if D==1
      coefficients=1
  else
      coefficients = condition .^ range(0,stop=1,length=D)
  end
  sum(coefficients .* (x.+3/2).^2)
end

function rastrigin(x)
  D = length(x)
  10 * D + sum(abs2, (x.+0.7)) - 10 * sum(xx -> cos.(2π * xx), (x.+0.7))
end

function sphere(x)
  sum(abs2, (x.-1.3))
end


function sin_01_x_2(x)
  sum(sin.((x.+0.7))+0.01.*(x.+0.7).^2)
end

function x_i(x)
  ress=0
  for i=1:length(x)
    ress+=(x[i]-i-2.1)^2
  end
  return ress
end

function x_5_sq(x)
  sum((x.-5).^2).-5
end

function x_plus_y(x)
    ress=0
    for i=1:length(x)
      ress+=x[i]
    end
    return ress
end

#


function rand_poly(x)
    return sum(rand(length(x)).*x)
end
    
# rand_poly(rand(10))
# x=-10:0.001:10
# y=sin.(10rand()x).+(rand()/1000)x.^2
# plot(x,y)

# function limited_domain(x::Vector)::Float64
#   sum( sqrt.(1 ./ (x.+1) .+ (x.+1) ) )  ## min is at 1,1,1...
# end#function
# limited_domain(rand(9))


# x=-10:0.01:10

# y=rand()rastrigin.(x)-.010rand()elliptic.(x)
# scatter(x,y)


function alpine_1(x)
  return sum(abs.(x.*sin.(x)+0.1x))
end