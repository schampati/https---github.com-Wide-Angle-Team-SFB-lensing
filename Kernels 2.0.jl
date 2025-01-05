using Pkg

Pkg.add("SpecialFunctions")
Pkg.add("Distributions")
Pkg.add("Plots")
Pkg.add("SphericalFourierBesselDecompositions")
Pkg.add("Dierckx")
Pkg.add("DelimitedFiles")
Pkg.add("LaTeXStrings")
using LaTeXStrings

using SpecialFunctions
using Plots
using Distributions
using SphericalFourierBesselDecompositions
using Dierckx
using DelimitedFiles
const SFB = SphericalFourierBesselDecompositions


#%%

# Defines the appropriate constants and plots the Hubble rate in terms of z, the redshift.

# Defines Ωs for the Friedmann Equations
Ω_K = 0
Ω_R = 5.5*10^(-5)
Ω_M = 0.3111
Ω_Λ = 1-Ω_R-Ω_M
H_0 = 100.0


# Defines the Constants for the Leading Coefficient of the Lensing Kernel
G = 6.67*10^(-11)
c = 2.998*10^(5)
ρ_c = (3*H_0^2)/(8*pi*G)
ρ_bar=Ω_M*ρ_c

# Paramters for the bessel functions
rmin = 200.0
rmax = 1000.0
kmax= 0.05

#This section defines the modes that we will use for the bessel functions and harmonics.

amodes = SFB.AnlmModes(kmax,rmin,rmax)


# Parameters for general kernel
b = 1.5
S = 0

#Number of intervals for Total Kernel

N_r = 100




#Defines a function to do an integration using the trapezoid method. First on a function f

function trapezoid_integral(f::Function,x_i,x_f,N)

    xs = range(x_i,x_f,N)
    Δx = xs[2]-xs[1]
    A = Δx*0.5*(2*sum(f.(xs))-f(xs[1])-f(xs[end]))
    return A
end


# Defines the Hubble rate and absolute distance in terms of redshift and the expansion metric.



H(z) = H_0 * sqrt((Ω_R*(1/(1+z))^(-4)+Ω_M*(1/(1+z))^(-3)+Ω_K*(1/(1+z))^(-2)+Ω_Λ))
H_2(a) = H_0 * sqrt((Ω_R*(a)^(-4)+Ω_M*(a)^(-3)+Ω_K*(a)^(-2)+Ω_Λ))

i1(z) = c/H(z)
r1(z) = trapezoid_integral(i1,0,z,1000)

i2(a) = (-c/a^2)*(1/H((1-a)/a))
r2(a) = trapezoid_integral(i2,1,a,1000)


# Creates a spline for the expansion rate in terms of the absolute distance.

amax = 1
amin = 0.1

as = range(amax,amin,1000)
rs = r2.(as)

chi_a = Dierckx.Spline1D(rs,as)




#Creates a matter power spectrum from data, to compare against generated power spectra.

P_m,header = readdlm("pk.tsv",header=true)
i=1
while 5==5
    if P_m[i,1]>kmax*1.1
        break
    end
    i=i+1
end
max_ind = i 


P_spline = Dierckx.Spline1D(P_m[:,1],P_m[:,2])

#%%
Δr = (rmax-rmin)/N_r
rs = range(rmin + Δr/2,rmax-Δr/2,length=N_r)




n_prime = 3
n= 2
l = 1

Δk = (kmax-0.0001)/max_ind
ks = range(0.0001+Δk/2,kmax-Δk/2,length=max_ind)




for m in range(1,max_ind) 
    k = ks[m]

    PSvalue = 0.0

    Lensing_kernel = 0.0
    Clustering_kernel = 0.0

    Lensing_kernel2 = 0.0
    Clustering_kernel2 = 0.0



    for r in rs


        A = 0.0
        A_prime = 0.0
    
        B = 0.0
        B_prime = 0.0

        r_Lmin = 0.00001
        Δr_L = (r-r_Lmin)/N_r
        r_Ls = range(r_Lmin+Δr_L/2,r-Δr_L/2,length=N_r)

        for r_L in r_Ls
            
            A = A + ((4*pi*G*ρ_bar)/(c^2))*l*(l+1)*r*amodes.basisfunctions.gnl[n,l+1](r)*((r-r_L)/(r_L))*sqrt(2/pi)*k*SpecialFunctions.sphericalbesselj(l,k*r_L)*(k^(-2))*(chi_a(r_L))^2
            
            A_prime = A_prime +  ((4*pi*G*ρ_bar)/(c^2))*l*(l+1)*r*amodes.basisfunctions.gnl[n_prime,l+1](r)*((r-r_L)/(r_L))*sqrt(2/pi)*k*SpecialFunctions.sphericalbesselj(l,k*r_L)*(k^(-2))*(chi_a(r_L))^2


        end

        B = B + r^2 * sqrt(2/pi) * k * amodes.basisfunctions.gnl[n,l+1](r) * SpecialFunctions.sphericalbesselj(l,k*r)

        B_prime = B_prime + r^2 * sqrt(2/pi) * k * amodes.basisfunctions.gnl[n_prime,l+1](r) * SpecialFunctions.sphericalbesselj(l,k*r)

        A = A*Δr_L

        A_prime = A*Δr_L

        Lensing_kernel = Lensing_kernel + A 

        Lensing_kernel2 = Lensing_kernel2 + A_prime

        Clustering_kernel = Clustering_kernel + B

        Clustering_kernel2 = Clustering_kernel2 + B_prime

    end

    Lensing_kernel = Lensing_kernel*Δr 
    Clustering_kernel = Clustering_kernel*Δr

    Lensing_kernel2 = Lensing_kernel2*Δr 
    Clustering_kernel2 = Clustering_kernel2*Δr


    TK= b*Clustering_kernel - 2*(1-2.5*S)*Lensing_kernel

    TK_prime = b*Clustering_kernel2 - 2*(1-2.5*S)*Lensing_kernel2


    PSvalue = PSvalue + TK_prime*TK*P_spline(k)*Δk

    print(PSvalue,"\n")
end





#%%
N_r = 1000
#%%


Δr = (rmax-rmin)/N_r
rs = range(rmin + Δr/2,rmax-Δr/2,length=N_r)



k = 4
n= 2
l = 1



Clustering_kernel = 0.0


for r in rs
        
    B = 0.0

    B = B + r^2 * sqrt(2/pi) * k * amodes.basisfunctions.gnl[n,l+1](r) * SpecialFunctions.sphericalbesselj(l,k*r)

    Clustering_kernel = Clustering_kernel + B

end

Clustering_kernel = Clustering_kernel*Δr

print(Clustering_kernel)

#%%


#LENSING KERNEL


Δr = (rmax-rmin)/N_r
rs = range(rmin + Δr/2,rmax-Δr/2,length=N_r)

k = 4
n= 2
l = 1



Lensing_kernel = 0.0

for r in rs


    A = 0.0


    r_Lmin = 0.00001
    Δr_L = (r-r_Lmin)/N_r
    r_Ls = range(r_Lmin+Δr_L/2,r-Δr_L/2,length=N_r)

    for r_L in r_Ls
            
        A = A + ((4*pi*G*ρ_bar)/(c^2))*l*(l+1)*r*amodes.basisfunctions.gnl[n,l+1](r)*((r-r_L)/(r_L))*sqrt(2/pi)*k*SpecialFunctions.sphericalbesselj(l,k*r_L)*(k^(-2))*(chi_a(r_L))^2
            

    end

    
    A = A*Δr_L

    Lensing_kernel = Lensing_kernel + A 

end

Lensing_kernel = Lensing_kernel*Δr 
   



print(Lensing_kernel)


#%%


f(r) = r^2 * sqrt(2/pi) * k * amodes.basisfunctions.gnl[n,l+1](r) * SpecialFunctions.sphericalbesselj(l,k*r)


prs = range(rmin,rmax,100)


plot(prs,f.(prs))





#%%