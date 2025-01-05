function Total_Kernel(k,n,l,rmin=rmin,rmax=rmax) 

    gnl_r(r) = amodes.basisfunctions.gnl[n,l+1](r)
    f(r,r_L) =  l*(l+1)*gnl_r(r)*k*SpecialFunctions.sphericalbesselj(l,k*r_L)*(k^(-2))
    As = zeros(N_r) 



    for i in range(1,length(rs))

        r = rs[i]
        A = 0 
        r_Ls = range(0.00001,r,step = Δr/10) #use center of each bin
        Δr_L = r_Ls[2]-r_Ls[1]

        for d in range(1,length(r_Ls))
            r_L=r_Ls[d]

            A = A + 2*f(r,r_L)*Rarray[i,d]
        end

        A = A - f(r,0.00001)*Rarray[i,1]-f(r,r)*Rarray[i,length(r_Ls)]
        A = A*Δr_L*0.5
        As[i]=A #take out
    end

    Lensing = Δr*0.5*(2*sum(As) - As[1]-As[end])

    f(r) = r^2 * sqrt(2/pi) * k * gnl_r(r) * SpecialFunctions.sphericalbesselj(l,k*r) #take r part out

    Clustering =  trapezoid_integral(f,rmin,rmax,N_r)

    Kernel = b*Clustering - 2*(1-2.5*S)*Lensing
    
    return Kernel

end