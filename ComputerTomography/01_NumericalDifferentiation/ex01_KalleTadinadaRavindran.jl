#Group: Kalle, Ravindran, Tadinada

using Random, Distributions, Plots
divz = 0.01
x = -1.0:divz:1.0

g(x) = exp(x)

gᵥ = g.(x)

dₙ = Normal(0.0, 0.01)

#Noise generated with seed, to ensure same results
noise = rand(MersenneTwister(314), dₙ, length(gᵥ))

#Noisy data
gᵈ = gᵥ + noise

plot(x, g.(x), label="g(x)", legend=:bottomright, xlabel="x", linealpha=0.8)
plot!(x, gᵈ, label="gᵈ(x)", linealpha=0.8)

dgᵈ(h,z) = (1/(2*h))*(gᵈ[Int(z+((h/divz)+1))] - gᵈ[Int(z-((h/divz)+1))])

dgᵈᵥ = zeros(length(gᵈ))

h = divz*15
hᵢ = findall(y->y==h, 0:divz:1.0)[1]

for i in 1:hᵢ
    global dgᵈᵥ[i] = gᵈ[i]
end

for i in (length(dgᵈᵥ)-hᵢ):length(dgᵈᵥ)
    global dgᵈᵥ[i] = gᵈ[i]
end

for i in (hᵢ+1):(length(gᵥ)-hᵢ)
    global dgᵈᵥ[i] = dgᵈ.(h,i)
end

plot(x, gᵥ, legend=:topleft, label="g'")
plot!(x[(hᵢ+1):(length(gᵥ)-hᵢ)],dgᵈᵥ[(hᵢ+1):(length(gᵥ)-hᵢ)], label="dgᵈ")
