using Pkg
using ForwardDiff
using Polynomials
using LinearAlgebra
using Printf
#Pkg.add("Sundials")
using Sundials  # for solving differential equations
#Pkg.add("PlotThemes")
using PlotThemes

# NAMING EXAMPLE v1prim MEANS v1' (derivative)

function fun(x, y, yprim)
    #extract mass from args
    mass0, mass1, mass2 = y[13], y[14], y[15]
   
    # y contains blocks data for every body
    # take position(vec) and speed(v) from y array
    vec0, v0 = y[1:2], y[3:4]
    vec1, v1 = y[5:6], y[7:8]
    vec2, v2 = y[9:10], y[11:12]
    
    # positions's derivative is speed
    vec0prim, vec1prim, vec2prim = v0, v1, v2
    
    #distances between three bodies
    d0  = (vec2 - vec1) / ( norm(vec2 - vec1)^3.0 )
    d1  = (vec0 - vec2) / ( norm(vec0 - vec2)^3.0 )
    d2  = (vec1 - vec0) / ( norm(vec1 - vec0)^3.0 )    
    
    # speed derivative == acceleration
    v0prim = mass1*d2 - mass2*d1
    v1prim = mass2*d0 - mass0*d2
    v2prim = mass0*d1 - mass1*d0
    
    yprim[:] = [vec0prim; v0prim; vec1prim; v1prim; vec2prim; v2prim;0;0;0]  # y size > because of three mass values so here also
end

using Plots

theme(:dark)
random = rand(50)
random1 = rand(50)
random2 = rand(50)
random3 = rand(50)


function calculate(m0, m1, m2, time)    

    # init positions and speeds 
    init0 = [ 1.0; -1.0; 0.0; 0.0]
    init1 = [ 1.0; 3.0; 0.0; 0.0]
    init2 = [ -2.0; -1.0; 0.0; 0.0]
    
    time_steps   = collect(range(0.0,stop = time, length = round(Int,time*500)))
    
    # to solve problem with type of m0, m1, m2 widgets outside calculate function 
    # I am passing masses through the initial values
    y0  = [init0; init1; init2; m0; m1; m2]  

    solution = Sundials.cvode(fun, y0, time_steps, reltol=1e-10)
   
    # take the positions and speeds from result of cvode function
    # of course not taking three zeros at the array end
    pos0, v0, pos1, v1, pos2, v2 = solution[:,1:2], solution[:,3:4], solution[:,5:6], solution[:,7:8], solution[:,9:10], solution[:,11:12]
    
    # calc mass centrum
    mass_centrum_x = [(pos0[i,1]*m0 + pos1[i,1]*m1 + pos2[i,1]*m2) / (m0 + m1 + m2) for i=1:length(time_steps)]
    mass_centrum_y = [(pos0[i,2]*m0 + pos1[i,2]*m1 + pos2[i,2]*m2) / (m0 + m1 + m2) for i=1:length(time_steps)]
    
    t1 = 1
    t2 = round(Int,(length(time_steps)-1) * time/time) + 1
    
    X,Y = 1,2
    
    #stars_x = vcat(random*xmin, random1*xmax, random*xmin, random1*xmax)
    #stars_y = vcat(random2*ymin, random3*ymax, random3*ymax, random2*ymin)
    
    
    if(t2-300 >= t1)
        plot(pos0[t2-300:t2-2,X], pos0[t2-300:t2-2,Y], size=(1000,500), legend= :outertopright, color="green", label="Earth trail",linewidth=3,ribbon =0.05) 
        #scatter!(stars_x, stars_y, markersize=0.5)
        scatter!(pos0[t2-2:t2,X], pos0[t2-2:t2,Y], label="Earth", color="green", markersize=m0)
        plot!(pos1[t2-300:t2-2,X], pos1[t2-300:t2-2,Y], color="grey", label="Jupiter trail",linewidth=3,ribbon =0.05)
        scatter!(pos1[t2-2:t2,X], pos1[t2-2:t2,Y], label ="Jupiter", color="grey", markersize=m1, markershape=:heptagon)
        plot!(pos2[t2-300:t2-2,X], pos2[t2-300:t2-2,Y], color="orange", label="Star trail",linewidth=3,ribbon =0.05)
        scatter!(pos2[t2-2:t2,X], pos2[t2-2:t2,Y], label ="Star", color="orange", markersize=m2, markershape=:star4)
        scatter!(mass_centrum_x[t2-1:t2], mass_centrum_y[t2-1:t2], label="centre of mass",linewidth=3)
    else
        plot(pos0[t1:t2-2,X], pos0[t1:t2-2,Y], size=(1000,500), legend= :outertopright, color="green", label="Earth trail",linewidth=3,ribbon =0.05)
        #scatter!(stars_x, stars_y, markersize=0.5)
        scatter!(pos0[t2-2:t2,X], pos0[t2-2:t2,Y], label="Earth", color="green", markersize=m0)
        plot!(pos1[t1:t2-2,X], pos1[t1:t2-2,Y], color="grey", label="Jupiter trail",linewidth=3,ribbon =0.05)
        scatter!(pos1[t2-2:t2,X], pos1[t2-2:t2,Y], label ="Jupiter", color="grey", markersize=m1, markershape=:heptagon)
        plot!(pos2[t1:t2-2,X], pos2[t1:t2-2,Y], color="orange", label="Star trail",linewidth=3,ribbon =0.05)
        scatter!(pos2[t2-2:t2,X], pos2[t2-2:t2,Y], label ="Star", color="orange", markersize=m2, markershape=:star4)
        scatter!(mass_centrum_x[t2-1:t2], mass_centrum_y[t2-1:t2], label="centre of mass",linewidth=3)
    end
end



using WebIO
#WebIO.install_jupyter_nbextension()
using Interact

@manipulate for m0 = widget(5:20, label="EARTH MASS"),
                m1 = widget(5:20, label="JUPITER MASS"),
                m2 = widget(5:20, label="STAR MASS"),
                s_len = widget(0.1:0.01:10, label="SIMULATION TIME")
    calculate(m0,m1,m2, s_len)
end


using WebIO
#WebIO.install_jupyter_nbextension()
using Interact

timer = Observable(0.0)
@async while true
    sleep(0.1)
    timer[]=timer[]+0.1
end


@manipulate for m0 = widget(5:20, label="EARTH MASS", color="green"),
                m1 = widget(5:20, label="JUPITER MASS"),
                m2 = widget(5:20, label="STAR MASS"),
                s_len = timer
    calculate(m0,m1,m2, s_len)
end
