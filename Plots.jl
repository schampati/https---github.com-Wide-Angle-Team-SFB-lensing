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
print(Distributions.Uniform([0,5],6))
#%%

# Defines the appropriate constants and plots the Hubble rate in terms of z, the redshift.


Ω_K = 0
Ω_R = 5.5*10^(-5)
Ω_M = 0.3111
Ω_Λ = 1-Ω_R-Ω_M
H_0 = 67.7

G = 6.67*10^(-11)
c = 2.998*10^(5)
ρ_c = (3*H_0^2)/(8*pi*G)
ρ_bar=Ω_M*ρ_c
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
glnr = @.amodes.basisfunctions.gnl[2,3](rr)
plot(glnr)
#%%

# how to get the j_l(kr) functions? in the specialfunctions module they don't have the k parameter.

function OuterFunction(n,l,rmin,rmax,kmax)

    amodes = SFB.AnlmModes(kmax,rmin,rmax)

    gnl_r(r) = amodes.basisfunctions.gnl[n,l+1](r)


    f(r)=l*(l+1)*r*gnl_r(r)

    return f

end

OuterFunction(1,1,rmin,rmax,kmax)(1)

#%%


function a_r_L(r2::Function,zmin,zmax,N)

    amax = 1/(1+zmin)
    amin = 1/(1+zmax)

    as = range(amax,amin,N)
    rs = r2.(as)

    chi_a = Dierckx.Spline1D(rs,as)

    return chi_a 

end

print(spline(r2,0,5,1000)(6000))
#%%

function LensingCoefficient()
    return ((4*pi*G*ρ_bar)/c)
end

#%%

function InnerFunction(k,l)

    y(r_L,r) = ((r-r_L)/(r_L))*sqrt(2/pi)*k*SpecialFunctions.sphericalbesselj(l,k*r_L)*(k^(-2))*(a_r_L(r2,0,5,1000)(r_L))^2
    return y

end

InnerFunction(1,1)(1,1)

#%%

function Integrand(k,l,n,rmin,rmax,kmax)

    f(r,r_L) = LensingCoefficient()*OuterFunction(n,l,rmin,rmax,kmax)(r)*InnerFunction(k,l)(r_L,r)
    return f 
end

Integrand(1,1,1,rmin,rmax,kmax)(2,1)

#%%

function LensingKernel(k,l,n,rmin=rmin,rmax=rmax,kmax=kmax,N_r = 10)

    f(r,r_L) = Integrand(k,l,n,rmin,rmax,kmax)(r,r_L)
    rs = range(rmin,rmax,length=N_r)
    Δr = rs[2]-rs[1]
    As=[]
    for r in rs
        A = 0 
        r_Ls = range(0,r,10)
        Δr_L = r_Ls[2]-r_Ls[1]
        # There's a problem here. It gets less precise with each step.
        for r_L in r_Ls
            A = A + (2*f(r,r_L))
        end
        A = A - f(r,0)-f(r,r)
        A = A*Δr_L
        append!(As,A)
        print(A)
    end

    A = Δr*0.5*(2*sum(As) - A[1]-A[end])

    return A

end

#%%

print(LensingKernel(5,2,2))




