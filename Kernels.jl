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

N_r = 400




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


#%%
rr = range(rmin,rmax,1000)

#%%
plot(H,0,3,xlabel="z",ylabel="H (km/s/Mpc)",title="Hubble Rate versus redshift")

#%%
plot(r1,0.00001,3,xlabel="z",ylabel="r (Mpc)",title="Distance versus Redshift")
#%%
plot(r2,0.000001,1,xlabel="a",ylabel="r (Mpc)",title="Distance versus Expansion Rate")
#%%


#%%
glnr = @.amodes.basisfunctions.gnl[4,3](rr)
plot(rr,glnr,xlabel="r",ylabel="g_nl",title="G_nl, n = 4, l = 2")
#%%





















#%%

# Creates a Function for the Lensing Kernel (for plotting purposes)

function LensingKernel(k,n,l ; rmin=rmin,rmax=rmax,kmax=kmax, N_r = 400) # use kmax to find N_r


    gnl_r(r) = amodes.basisfunctions.gnl[n,l+1](r)

    f(r,r_L) =  ((4*pi*G*ρ_bar)/(c^2))*l*(l+1)*r*gnl_r(r)*((r-r_L)/(r_L))*sqrt(2/pi)*k*SpecialFunctions.sphericalbesselj(l,k*r_L)*(k^(-2))*(chi_a(r_L))^2

    rs = range(rmin,rmax,length=N_r) 
    Δr = rs[2]-rs[1]
    As = zeros(N_r) #pre-allocate array

    for i in range(1,length(rs))

        r = rs[i]
        A = 0 
        r_Ls = range(0.00001,r,step = Δr/10)
        Δr_L = r_Ls[2]-r_Ls[1]

        for r_L in r_Ls
            A = A + 2*f(r,r_L)
        end


        A = A - f(r,0.00001)-f(r,r)
        A = A*Δr_L*0.5
        As[i]=A #take out

    end


    Lensing = Δr*0.5*(2*sum(As) - As[1]-As[end])

    return Lensing

end


#%%

LensingKernel(4,2,1)


#%%

Clustering(4,2,1)

#%%
function Clustering(k,n,l,rmin=rmin,rmax=rmax)

    f(r) = r^2 * sqrt(2/pi) * k * amodes.basisfunctions.gnl[n,l+1](r) * SpecialFunctions.sphericalbesselj(l,k*r) #take r part out


    xs = range(rmin,rmax,N_r)
    Δx = xs[2]-xs[1]
    A = Δx*0.5*(2*sum(f.(xs))-f(xs[1])-f(xs[end]))
    Clustering = A

    return Clustering

end

#%%

#%%

function Total_Kernel(k,n,l,rmin=rmin,rmax=rmax,S=S) # Original version of Total_Kernel, leaving gnl_r as a function and evaluating it for each loop iteration.

    gnl_r(r) = amodes.basisfunctions.gnl[n,l+1](r)

    f(r,r_L) =  ((4*pi*G*ρ_bar)/(c^2))*l*(l+1)*r*gnl_r(r)*((r-r_L)/(r_L))*sqrt(2/pi)*k*SpecialFunctions.sphericalbesselj(l,k*r_L)*(k^(-2))*(chi_a(r_L))^2

    rs = range(rmin,rmax,length=N_r) 
    Δr = rs[2]-rs[1]
    As = zeros(N_r) #pre-allocate array

    for i in range(1,length(rs))

        r = rs[i]
        A = 0 
        r_Ls = range(0.00001,r,step = Δr/10)
        Δr_L = r_Ls[2]-r_Ls[1]

        for r_L in r_Ls

            A = A + 2*f(r,r_L)
        end

        A = A - f(r,0.00001)-f(r,r)
        A = A*Δr_L*0.5
        As[i]=A #take out

    end

    Lensing = Δr*0.5*(2*sum(As) - As[1]-As[end])
    f(r) = r^2 * sqrt(2/pi) * k * gnl_r(r) * SpecialFunctions.sphericalbesselj(l,k*r) #take r part out
    Clustering =  trapezoid_integral(f,rmin,rmax,N_r)
    print("Clustering:",Clustering)
    Kernel = b*Clustering - 2*(1-2.5*S)*Lensing
    return Kernel

end



#%%
Total_Kernel(4,2,1)
#%%
@time(Total_Kernel(1,1,1))
@time(Total_KernelBeta(1,1,1))
@time(Total_KernelNumber3(1,1,1))
#%%


#%%

function LensingPowerSpectrum(n,n_prime,l,max_ind=max_ind) #First version of LensingPowerSpectrum, uses Total_Kernel Function

    A = 0

    for m in range(1,max_ind-1) 

        A = A + (P_m[m+1,1]-P_m[m,1]) * (  (Total_Kernel(P_m[m+1,1],n,l) *Total_Kernel(P_m[m+1,1],n_prime,l))  *P_m[m+1,2]  + (Total_Kernel(P_m[m,1],n,l) * Total_Kernel(P_m[m,1],n_prime,l))  *P_m[m,2]  )/2
        
    end

    return A

end
#%%




#%%
@time(LensingPowerSpectrum(2,3,1))
#%%
@time(LensingPowerSpectrumnofunctions(1,1,1))


#%%
@time(LensingPowerSpectrumBETA(1,1,1))
@time(LensingPowerSpectrumArray(1,1,1))
#%%



function NoLensingPowerSpectrum(n,n_prime,l,max_ind=max_ind,S=S)

    f(k) = b^2*Clustering(k,n,l)*Clustering(k,n_prime,l)

    A=0

    for m in range(1,max_ind)

        A = A + (P_m[m+1,1]-P_m[m,1]) * (f(P_m[m+1,1])*P_m[m+1,2]+f(P_m[m,1])*P_m[m,2])/2
    
    end
    
    return A
    
end


#%%

#The next few sections plot the power spectra (no lensing, lensing, and matter power spectrum)

xs = amodes.knl[ .! isnan.(amodes.knl) ]
N_r = 400
l_test=1
n_max_test=amodes.nmax_l[l_test+1]
n_range_test = 1:n_max_test

#%%
k_nl_test = amodes.knl[n_range_test,l_test+1]
length(k_nl_test)
#%%

@time(C_lnn_test = NoLensingPowerSpectrum.(n_range_test,n_range_test,l_test))
@time(C_lnn_test2 = LensingPowerSpectrum.(n_range_test,n_range_test,l_test))

#%%
plot(k_nl_test,C_lnn_test2,xlabel=L"$k (Mpc^{-1})$",ylabel=L"C_lnn (h^-1 Mpc)^3",label = L"$\tilde{C}_{lnn}(k)$",title=L"Power Spectra, $l=1, n = n'$")
plot!(P_m[1:max_ind,1],b*b*P_m[1:max_ind,2],label=L"$b^2P(k)$")
plot!(k_nl_test,C_lnn_test,ylims = (30000,59000),xlabel=L"$k (Mpc^{-1})$",ylabel=L"C_{lnn} (h^{-1} Mpc)^3",label=L"C_{lnn}(k)")

#%%
plot(k_nl_test,C_lnn_test-C_lnn_test2,xlabel=L"$k (Mpc^{-1})$",ylabel=L"C_{lnn}- \tilde{C}_{lnn} (h^{-1} Mpc)^3",label =L"C_{lnn}- \tilde{C}_{lnn}",title = L"Power Spectrum Difference, $l=1, n = n'$" )
#%%

LensingPowerSpectrum(12,12,l_test)
#%%








# This sections plots the total kernel, and the two components of the kernel (clustering, and lensing) separately.


#%%
n_test = 2
l_test = 1

krange = range(0,kmax,50)
y1 = Total_Kernel.(krange,n_test,l_test)
y2 = LensingKernel.(krange,n_test,l_test)


#%%
y3 = Clustering.(krange,n_test,l_test)
#%%
plot(krange,y1, title = L"$\tilde{W}_{2,1}(k)$",xlabel=L"$k (Mpc^{-1})$",ylabel = L"$\tilde{W}_{2,1}(k) \ [h^{-1} Mpc]^{\frac{1}{2}}$",label= L"\tilde{W}_{2,1}(k)" )
vline!([amodes.knl[3,2+1]],label=L"$k_{nl}$ mode")
#%%
plot(krange,y2,title = L"$K_{2,1}(k)$" ,xlabel=L"$k (Mpc^{-1})$",ylabel = L" $K_{2,1}(k) \ [h^{-1} Mpc]^{\frac{1}{2}}$", label = L"K_{2,1}(k)")
vline!([amodes.knl[3,2+1]],label=L"$k_{nl}$ mode")
#%%
plot(krange,y3,title=L"$W_{2,1}(k)$",xlabel=L"$k (Mpc^{-1})$",ylabel = L"$W_{2,1}(k) \ [h^{-1} Mpc]^{\frac{1}{2}}$",label = L"W_{2,1}(k)")
vline!([amodes.knl[3,2+1]],label=L"$k_{nl}$ mode")
#%%



























#%%


function Total_KernelNumber3(k,n,l,rmin=rmin,rmax=rmax,S=S) # A re-work of Total_Kernel, writing pre-defining g_nl for all the r values instead of creating a function and running it each time.

    rs = range(rmin,rmax,length=N_r) 
    Δr = rs[2]-rs[1]
    As = zeros(N_r) #pre-allocate array

    for i in range(1,length(rs))

        r = rs[i]
        A = 0.0 
        r_Ls = range(0.00001,r,step = Δr/10)
        Δr_L = r_Ls[2]-r_Ls[1]

        for i in range(1,length(r_Ls)-1)
            r_L = 0.5*(r_Ls[i]+r_Ls[i+1])
            A = A + 2*((4*pi*G*ρ_bar)/(c^2))*l*(l+1)*r*amodes.basisfunctions.gnl[n,l+1](r)*((r-r_L)/(r_L))*sqrt(2/pi)*k*SpecialFunctions.sphericalbesselj(l,k*r_L)*(k^(-2))*(chi_a(r_L))^2
        end

        A = A*Δr_L*0.5
        As[i]=A #take out

    end

    Lensing = Δr*0.5*(2*sum(As) - As[1]-As[end])
    f(r) = r^2 * sqrt(2/pi) * k * amodes.basisfunctions.gnl[n,l+1](r) * SpecialFunctions.sphericalbesselj(l,k*r) #take r part out
    Clustering = Δr*0.5*(2*sum(rs.^2 .* sqrt(2/pi) .* k .* amodes.basisfunctions.gnl[n,l+1].(rs) .* SpecialFunctions.sphericalbesselj.(l,k.*rs)) - (rmin^2 * sqrt(2/pi) * k * amodes.basisfunctions.gnl[n,l+1](rmin) * SpecialFunctions.sphericalbesselj(l,k*rmin))-(rmax^2 * sqrt(2/pi) * k * amodes.basisfunctions.gnl[n,l+1](rmax) * SpecialFunctions.sphericalbesselj(l,k*rmax)))
    Kernel = b*Clustering - 2*(1-2.5*S)*Lensing
    return Kernel

end


#%%


function Total_KernelBeta(k,n,l,rmin=rmin,rmax=rmax,S=S) # A re-work of Total_Kernel, writing pre-defining g_nl for all the r values instead of creating a function and running it each time.

    rs = range(rmin,rmax,length=N_r) 

    gnl_rarray = amodes.basisfunctions.gnl[n,l+1].(rs)

    f2(r_L) = (rs.*gnl_rarray.*((rs/r_L)[:].-1))
    f1(r_L) =  (((4*pi*G*ρ_bar)/(c^2))*l*(l+1)*sqrt(2/pi)*k*SpecialFunctions.sphericalbesselj(l,k*r_L)*(k^(-2))*(chi_a(r_L))^2)

    Δr = rs[2]-rs[1]
    As = zeros(N_r) #pre-allocate array

    for i in range(1,length(rs))

        r = rs[i]
        A = 0 
        r_Ls = range(0.00001,r,step = Δr/10)
        Δr_L = r_Ls[2]-r_Ls[1]

        for r_L in r_Ls
            # Insert f1 and f2
            A = A + 2*f1(r_L)*f2(r_L)[i]
        end


        A = A - f1(0.00001)*f2(0.00001)[i]-f1(r)*f2(r)[i]
        A = A*Δr_L*0.5
        As[i]=A #take out

    end

    Lensing = Δr*0.5*(2*sum(As) - As[1]-As[end])

    f(r) = r^2 * sqrt(2/pi) * k * amodes.basisfunctions.gnl[n,l+1](r) * SpecialFunctions.sphericalbesselj(l,k*r) #take r part out

    Clustering =  trapezoid_integral(f,rmin,rmax,N_r)

    Kernel = b*Clustering - 2*(1-2.5*S)*Lensing
    
    return Kernel

end

#%%

function LensingPowerSpectrumArray(n,n_prime,l,max_ind=max_ind) #Third version of LensingPowerSpectrum, defining the gnl_r functions for each value of r.
    rs = range(rmin,rmax,length=N_r)
    
    gnl1_r = amodes.basisfunctions.gnl[n,l+1].(rs)
    gnl2_r = amodes.basisfunctions.gnl[n_prime,l+1].(rs)
    PowerSpectrumValue = 0 
    
    for m in range(1,max_ind-1) 

        k1 = P_m[m+1,1]

        k2 = P_m[m,1]
        
        

        f1(r_L) =  ((4*pi*G*ρ_bar)/(c^2))*l*(l+1)*(rs.*gnl1_r.*((rs.-r_L)/(r_L)))*sqrt(2/pi)*k1*SpecialFunctions.sphericalbesselj(l,k1*r_L)*(k1^(-2))*(chi_a(r_L))^2
        f2(r_L) = ((4*pi*G*ρ_bar)/(c^2))*l*(l+1)*(rs.*gnl2_r.*((rs.-r_L)/(r_L)))*sqrt(2/pi)*k1*SpecialFunctions.sphericalbesselj(l,k1*r_L)*(k1^(-2))*(chi_a(r_L))^2
        f3(r_L) = ((4*pi*G*ρ_bar)/(c^2))*l*(l+1)*(rs.*gnl1_r.*((rs.-r_L)/(r_L)))*sqrt(2/pi)*k2*SpecialFunctions.sphericalbesselj(l,k2*r_L)*(k2^(-2))*(chi_a(r_L))^2
        f4(r_L) = ((4*pi*G*ρ_bar)/(c^2))*l*(l+1)*(rs.*gnl2_r.*((rs.-r_L)/(r_L)))*sqrt(2/pi)*k2*SpecialFunctions.sphericalbesselj(l,k2*r_L)*(k2^(-2))*(chi_a(r_L))^2
 

        Δr = rs[2]-rs[1]
        As = zeros(N_r) #pre-allocate array
        Bs = zeros(N_r)
        Cs = zeros(N_r) #pre-allocate array
        Ds = zeros(N_r)


        for i in range(1,length(rs))

            r = rs[i]
            A = 0 
            B = 0
            C = 0
            D = 0
            r_Ls = range(0.00001,r,step = Δr/10)
            Δr_L = r_Ls[2]-r_Ls[1]

            for r_L in r_Ls

                A = A + 2*f1(r_L)[i]
                B = B + 2*f2(r_L)[i]
                C = C + 2*f3(r_L)[i]
                D = D + 2*f4(r_L)[i]

            end


            A = A - f1(0.00001)[i]-f1(r)[i]
            A = A*Δr_L*0.5
            As[i]=A #take out

            B = B - f2(0.00001)[i]-f2(r)[i]
            B = B*Δr_L*0.5
            Bs[i]=B #take out

            C = C - f3(0.00001)[i]-f3(r)[i]
            C = C*Δr_L*0.5
            Cs[i]=C #take out

            D = D - f4(0.00001)[i]-f4(r)[i]
            D = D*Δr_L*0.5
            Ds[i]=D #take out

        end

        Lensing1 = Δr*0.5*(2*sum(As) - As[1]-As[end])
        cint1 = rs.^2 .* sqrt(2/pi) .* k1 .* gnl1_r .* SpecialFunctions.sphericalbesselj.(l,k1*rs) #take r part out
        Clustering1 = Δr*0.5*(2*sum(cint1)-cint1[1]-cint1[end])
        Kernel1 = b*Clustering1 - 2*(1-2.5*S)*Lensing1

        Lensing2 = Δr*0.5*(2*sum(Bs) - Bs[1]-Bs[end])
        cint2 = rs.^2 .* sqrt(2/pi) .* k1 .* gnl2_r .* SpecialFunctions.sphericalbesselj.(l,k1*rs) #take r part out
        Clustering2 = Δr*0.5*(2*sum(cint2)-cint2[1]-cint2[end])
        Kernel2 = b*Clustering2 - 2*(1-2.5*S)*Lensing2

        Lensing3 = Δr*0.5*(2*sum(Cs) - Cs[1]-Cs[end])
        cint3 = rs.^2 .* sqrt(2/pi) .* k2 .* gnl1_r .* SpecialFunctions.sphericalbesselj.(l,k2*rs) #take r part out
        Clustering3 = Δr*0.5*(2*sum(cint3)-cint3[1]-cint3[end])
        Kernel3 = b*Clustering3 - 2*(1-2.5*S)*Lensing3

        Lensing4 = Δr*0.5*(2*sum(Ds) - Ds[1]-Ds[end])
        cint4 = rs.^2 .* sqrt(2/pi) .* k2 .* gnl2_r .* SpecialFunctions.sphericalbesselj.(l,k2*rs) #take r part out
        Clustering4 = Δr*0.5*(2*sum(cint4)-cint4[1]-cint4[end])
        Kernel4 = b*Clustering4 - 2*(1-2.5*S)*Lensing4
        
        
        PowerSpectrumValue = PowerSpectrumValue + (P_m[m+1,1]-P_m[m,1]) * (  (Kernel1 * Kernel2)  *P_m[m+1,2]  + (Kernel3 * Kernel4)  *P_m[m,2]  )/2

    end
    return PowerSpectrumValue


end





#%%

function LensingPowerSpectrumBETA(n,n_prime,l,max_ind=max_ind) #Second version of Lensing, doesn't use the function Total_Kernel, but still defines gnl1_r as a function and passes through values in each iteration.

    
    gnl1_r(r) = amodes.basisfunctions.gnl[n,l+1](r)
    gnl2_r(r) = amodes.basisfunctions.gnl[n_prime,l+1](r)
    PowerSpectrumValue = 0 
    
    for m in range(1,max_ind-1) 

        k1 = P_m[m+1,1]

        k2 = P_m[m,1]
        


        f1(r,r_L) =  ((4*pi*G*ρ_bar)/(c^2))*l*(l+1)*r*gnl1_r(r)*((r-r_L)/(r_L))*sqrt(2/pi)*k1*SpecialFunctions.sphericalbesselj(l,k1*r_L)*(k1^(-2))*(chi_a(r_L))^2
        f2(r,r_L) = ((4*pi*G*ρ_bar)/(c^2))*l*(l+1)*r*gnl2_r(r)*((r-r_L)/(r_L))*sqrt(2/pi)*k1*SpecialFunctions.sphericalbesselj(l,k1*r_L)*(k1^(-2))*(chi_a(r_L))^2
        f3(r,r_L) = ((4*pi*G*ρ_bar)/(c^2))*l*(l+1)*r*gnl1_r(r)*((r-r_L)/(r_L))*sqrt(2/pi)*k2*SpecialFunctions.sphericalbesselj(l,k2*r_L)*(k2^(-2))*(chi_a(r_L))^2
        f4(r,r_L) = ((4*pi*G*ρ_bar)/(c^2))*l*(l+1)*r*gnl2_r(r)*((r-r_L)/(r_L))*sqrt(2/pi)*k2*SpecialFunctions.sphericalbesselj(l,k2*r_L)*(k2^(-2))*(chi_a(r_L))^2
 

        rs = range(rmin,rmax,length=N_r)
        Δr = rs[2]-rs[1]
        As = zeros(N_r) #pre-allocate array
        Bs = zeros(N_r)
        Cs = zeros(N_r) #pre-allocate array
        Ds = zeros(N_r)

        for i in range(1,length(rs))

            r = rs[i]
            A = 0 
            B = 0
            C = 0
            D = 0
            r_Ls = range(0.00001,r,step = Δr/10)
            Δr_L = r_Ls[2]-r_Ls[1]

            for r_L in r_Ls


                A = A + 2*f1(r,r_L)
                B = B + 2*f2(r,r_L)
                C = C + 2*f3(r,r_L)
                D = D + 2*f4(r,r_L)

            end

           

            A = A - f1(r,0.00001)-f1(r,r)
            A = A*Δr_L*0.5
            As[i]=A #take out

            B = B - f2(r,0.00001)-f2(r,r)
            B = B*Δr_L*0.5
            Bs[i]=B #take out

            C = C - f3(r,0.00001)-f3(r,r)
            C = C*Δr_L*0.5
            Cs[i]=C #take out

            D = D - f4(r,0.00001)-f4(r,r)
            D = D*Δr_L*0.5
            Ds[i]=D #take out

        end

        Lensing1 = Δr*0.5*(2*sum(As) - As[1]-As[end])
        cint1(r) = r^2 * sqrt(2/pi) * k1 * gnl1_r(r) * SpecialFunctions.sphericalbesselj(l,k1*r) #take r part out
        Clustering1 =  trapezoid_integral(cint1,rmin,rmax,N_r)
        Kernel1 = b*Clustering1 - 2*(1-2.5*S)*Lensing1

        Lensing2 = Δr*0.5*(2*sum(Bs) - Bs[1]-Bs[end])
        cint2(r) = r^2 * sqrt(2/pi) * k1 * gnl2_r(r) * SpecialFunctions.sphericalbesselj(l,k1*r) #take r part out
        Clustering2 =  trapezoid_integral(cint2,rmin,rmax,N_r)
        Kernel2 = b*Clustering2 - 2*(1-2.5*S)*Lensing2

        Lensing3 = Δr*0.5*(2*sum(Cs) - Cs[1]-Cs[end])
        cint3(r) = r^2 * sqrt(2/pi) * k2 * gnl1_r(r) * SpecialFunctions.sphericalbesselj(l,k2*r) #take r part out
        Clustering3 =  trapezoid_integral(cint3,rmin,rmax,N_r)
        Kernel3 = b*Clustering3 - 2*(1-2.5*S)*Lensing3

        Lensing4 = Δr*0.5*(2*sum(Ds) - Ds[1]-Ds[end])
        cint4(r) = r^2 * sqrt(2/pi) * k2 * gnl2_r(r) * SpecialFunctions.sphericalbesselj(l,k2*r) #take r part out
        Clustering4 =  trapezoid_integral(cint4,rmin,rmax,N_r)
        Kernel4 = b*Clustering4 - 2*(1-2.5*S)*Lensing4
        
        
        PowerSpectrumValue = PowerSpectrumValue + (P_m[m+1,1]-P_m[m,1]) * (  (Kernel1 * Kernel2)  *P_m[m+1,2]  + (Kernel3 * Kernel4)  *P_m[m,2]  )/2

    end


    return PowerSpectrumValue


end

#%%

function LensingPowerSpectrumnofunctions(n,n_prime,l,max_ind=max_ind) #First version of LensingPowerSpectrum, uses Total_Kernel Function

    A = 0.0

    for m in range(1,max_ind-1) 

        A = A + (P_m[m+1,1]-P_m[m,1]) * (  (Total_KernelNumber3(P_m[m+1,1],n,l) *Total_KernelNumber3(P_m[m+1,1],n_prime,l))  *P_m[m+1,2]  + (Total_KernelNumber3(P_m[m,1],n,l) * Total_KernelNumber3(P_m[m,1],n_prime,l))  *P_m[m,2]  )/2
        
    end
    return A

end





#%%





































