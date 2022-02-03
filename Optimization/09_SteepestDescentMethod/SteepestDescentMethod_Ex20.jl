using LinearAlgebra, Plots

x0 = 2.0
fk = 15
f(x) = (0.125)*(x^2)
∇f(x)= 0.25*x
plot(xlim=(1,fk), ylim=(0,2.05))
b = 1

while x0>0.0
    print(x0)
    k = 1
    Ws = [x0]

    γ = 0.5

    τ = 0.5
    while k<fk
        w = Ws[k]
        sk = ((-0.5)^(float(k)))*∇f(w)
        σₖ=1.0
        a = 0.0
        while a<10000.0
            σₖ = τ*σₖ
            j = f(w+(σₖ*sk)) - f(w)
            α = γ*σₖ*(∇f(w)*sk)
            if j <= α
                a = 10000.0 + 10
            end

            a=a+1
        end
        w_new = w + (σₖ*sk)
        Ws=[Ws;[w_new]]
        if ∇f(w_new)==0
             k = fk+10
        end
         k = k+1
    end

    plot!(f.(Ws), marker=true, markersize = 0.5, linealpha=0.75, color=palette(:tab10)[b])
    plot!(repeat([0.5*x0], fk), linestyle=:dot, color=palette(:tab10)[b])
    global x0 = x0 - 0.25
    global b = b+1
end

plot!(legend=:false)
