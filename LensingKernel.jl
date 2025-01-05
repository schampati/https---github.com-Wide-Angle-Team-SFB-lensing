using Pkg

Pkg.add("SpecialFunctions")
Pkg.add("Distributions")
Pkg.add("Plots")
Pkg.add("SphericalFourierBesselDecompositions")
Pkg.add("Dierckx")
using SpecialFunctions
using Plots
using Distributions
using SphericalFourierBesselDecompositions
using Dierckx
const SFB = SphericalFourierBesselDecompositions

#%%

# Defines the appropriate constants and plots the Hubble rate in terms of z, the redshift.

# Defines Ωs for the Friedmann Equations
Ω_K = 0
Ω_R = 5.5*10^(-5)
Ω_M = 0.3111
Ω_Λ = 1-Ω_R-Ω_M
H_0 = 67.7


# Defines the Constants for the Leading Coefficient of the Lensing Kernel
G = 6.67*10^(-11)
c = 2.998*10^(5)
ρ_c = (3*H_0^2)/(8*pi*G)
ρ_bar=Ω_M*ρ_c

# Parameters for 
rmin = 200.0
rmax = 1000.0
kmax=0.05

#%%


H(z) = H_0 * sqrt((Ω_R*(1/(1+z))^(-4)+Ω_M*(1/(1+z))^(-3)+Ω_K*(1/(1+z))^(-2)+Ω_Λ))
H_2(a) = H_0 * sqrt((Ω_R*(a)^(-4)+Ω_M*(a)^(-3)+Ω_K*(a)^(-2)+Ω_Λ))
plot(H,0,3,xlabel="z",ylabel="H (km/s/Mpc)",title="Hubble Rate versus redshift")

#%%

#Defines a function to do an integration using the trapezoid method. First on a function f

function trapezoid_integral(f::Function,x_i,x_f,N)

    xs = range(x_i,x_f,N)
    Δx = xs[2]-xs[1]
    A = Δx*0.5*(2*sum(f.(xs))-f(xs[1])-f(xs[end]))
    return A
end



# N needs to be equal to the length of fxs. Fxs is the array of function values.
function trapezoid_integral(f::Function,g::Function,x_i,x_f,y_i,y_f,N)
    xs = range(x_i,x_f,N)
    ys = range(y_i,y_f,N)
    Δx = xs[2]-xs[1]
    Δy = ys[2]-ys[1]
    for m1=1:N
        for m2=1:N

        end
    end
    return A
end



#%%


i1(z) = c/H(z)
r1(z) = trapezoid_integral(i1,0,z,1000)
plot(r1,0.00001,3,xlabel="z",ylabel="r (Mpc)",title="Distance versus Redshift")
#%%
i2(a) = (-c/a^2)*(1/H((1-a)/a))
r2(a) = trapezoid_integral(i2,1,a,1000)
plot(r2,0.000001,1,xlabel="a",ylabel="r (Mpc)",title="Distance versus Expansion Rate")
#%%
print(r1(0.25))
print(r2(0.8))
#%%
#This section calculates the kernel K_nl (k) for the lensing magnification term.

amodes = SFB.AnlmModes(kmax,rmin,rmax)
rr = range(200.0,1000.0,1000)
glnr = @.amodes.basisfunctions.gnl[4,3](rr)
plot(rr,glnr,xlabel="r",ylabel="g_nl",title="G_nl, n = 4, l = 2")
#%%

# how to get the j_l(kr) functions? in the specialfunctions module they don't have the k parameter.

amodes = SFB.AnlmModes(kmax,rmin,rmax) # take outside function


function OuterFunction(n,l,rmin,rmax,kmax)

    gnl_r(r) = amodes.basisfunctions.gnl[n,l+1](r)


    f(r)=l*(l+1)*r*gnl_r(r)

    return f

end

OuterFunction(1,1,rmin,rmax,kmax)(1)


#%%


function a_r_L(r2::Function,zmin,zmax,N)

    amax = 1/(1+zmin) # rmin/rmax
    amin = 1/(1+zmax) # rmin/rmax

    as = range(amax,amin,N)
    rs = r2.(as)

    chi_a = Dierckx.Spline1D(rs,as)

    return chi_a 

end

print(a_r_L(r2,0,5,1000)(6000))
#%%

function LensingCoefficient()
    return ((4*pi*G*ρ_bar)/c)
end

LensingCoefficient()

#%%

function InnerFunction(k,l)

    y(r_L,r) = ((r-r_L)/(r_L))*sqrt(2/pi)*k*SpecialFunctions.sphericalbesselj(l,k*r_L)*(k^(-2))*(a_r_L(r2,0,5,1000)(r_L))^2
    return y

end

InnerFunction(1,1)(1,1)

#%%

function Integrand(k,l,n,rmin,rmax,kmax)

    gnl_r(r) = amodes.basisfunctions.gnl[n,l+1](r)

    f(r,r_L) =  ((4*pi*G*ρ_bar)/c)*l*(l+1)*r*gnl_r(r)*((r-r_L)/(r_L))*sqrt(2/pi)*k*SpecialFunctions.sphericalbesselj(l,k*r_L)*(k^(-2))*(a_r_L(r2,0,5,1000)(r_L))^2
    return f 
end

f(x,y) = Integrand(1,1,1,100.5,200.6,0.05)(x,y)
print(Integrand(1,1,1,rmin,rmax,kmax)(1.0,0.00000000001))
f(1,2)*5

#%%

function LensingKernel(k,l,n ; rmin=rmin,rmax=rmax,kmax=kmax, N_r = 10)

# @time

    gnl_r(r) = amodes.basisfunctions.gnl[n,l+1](r)

    f(r,r_L) =  ((4*pi*G*ρ_bar)/c)*l*(l+1)*r*gnl_r(r)*((r-r_L)/(r_L))*sqrt(2/pi)*k*SpecialFunctions.sphericalbesselj(l,k*r_L)*(k^(-2))*(a_r_L(r2,0,5,1000)(r_L))^2


    print(f(1.111,6.6666))


    rs = range(rmin,rmax,length=N_r) 
    Δr = rs[2]-rs[1]
    As = zeros(N_r) #pre-allocate array



    for r in rs

        A = 0 
        r_Ls = range(0.00001,r,step = Δr/10)
        Δr_L = r_Ls[2]-r_Ls[1]

        for r_L in r_Ls

            A = A + 2*f(r,r_L)
        end




        A = A - f(r,0.00001)-f(r,r)
        A = A*Δr_L
        append!(As,A)
    end

    A = Δr*0.5*(2*sum(As) - A[1]-A[end])

    print(A,"that was A")

end

#%%

LensingKernel(5,2,2)

#%%

N_r = 10

f(r,r_L) = Integrand(5,2,2,rmin,rmax,kmax)(r,r_L)

rs = range(0,rmax,length=N_r) #Do we start with rmin as the integration bounds, or just do from zero to rmax

for r in rs
    A = 0 
    r_Ls = range(0,r,10)

    for r_L in r_Ls
        print(f(r,r_L))
    end

end


#%%
y = range(1,10,step=0.5)

y[1]
y[2]
y[3]
y[4]

