

#%%

using Pkg

Pkg.add("SpecialFunctions")
Pkg.add("Plots")
using SpecialFunctions
using Plots

real_axis = range(-10,10,100)
radial_axis = range(0,10,100)
im_axis = range(-10,10,100)*im

#%%

plot(real_axis,besselj1.(real_axis))

plot(real_axis,besselj0.(real_axis))

plot(real_axis,sphericalbesselj.(0,radial_axis))

plot(real_axis,sphericalbesselj.(1,radial_axis))


#%%

G = 10
rho_bar = 100
a = 10
#%%
#Calculate the 


#%%

'''




'''

