using LinearAlgebra, Plots, DataFrames, CSV

#Inital guess
x0=[1, -0.5]

#Rosenbrock function
f(x) = ((1-x[1])^2) + (100*((x[2] - (x[1]^2))^2))

#Gradient of f(x)
∇f(x)= 2*[-x[1]*(1+(200*(x[2]-(x[1]^2)))),
            100*(x[2] - (x[1]^2))]

k = 1

#Ws is an array that will contain the history of points towards to the required solution.
Ws = [x0]

Sk= [[0,1], [1,0], [1,1],
    [0,-1], [-1,0], [-1,-1],
    [1,-1],[-1,1]]

γ = 0.5
τ = 0.5

fk = 100
while k<fk
    if k == 1
        w = Ws[k]
    else w = Ws[k]'
    end
    for i in 1:length(Sk)
        s_temp = dot(∇f(w), Sk[i])
        if s_temp<0
            global sk = Sk[i]
        end
    end

    σₖ=2.0  #initating Armijo step size [it is 2 because in the loop, we take σₖ = τ*σₖ and the first guess will become 1.0]
    a = 0.0
    while a<10000.0
    #to stop the possible infinite search for an Armijo step size, we restrict the loop to 10000 iterations
        σₖ = τ*σₖ
        # defining the RHS and LHS of the Armijo condition
        j = f(w+(σₖ*sk)) - f(w)
        α = γ*σₖ*dot(∇f(w), sk)
        if j <= α   #Amrijo condition
            #if the Amrijo condition is fullfilled, we redefine 'a' to break the search loop
            a = 10000.0 + 10
        end
        a=a+1
    end

    w_new = w + (σₖ*sk)

    # appending new point to Ws
    global Ws=[Ws;[w_new]']
    if ∇f(w_new)==0
        global k = fk+10
    end
    global k = k+1
end

#Changing the orientation of the elements of Ws to easily view it
for i in 2:length(Ws)
    Ws[i] = Ws[i]'
end

#Creating a .csv to clearly view the output in Spreadsheet softwares
Wf = Float64.(zeros(length(Ws), 3))
for i in 1:length(Ws)
    Wf[i, 1] = Ws[i][1]
    Wf[i, 2] = Ws[i][2]
    Wf[i, 3] = f(Ws[i])
end
CSV.write("WfEx19.csv", DataFrame(Wf, :auto))
