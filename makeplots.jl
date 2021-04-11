using Plots
include("isolatedmodel.jl")
function getn(c)
    if length(size(c))==3
        return size(c)[1],size(c)[3]
    else
        return size(c)
    end
end

function singleconcplot(r,c)
    n_spatial,ntime = getn(c)
    c10min= c[:,ntime÷6]
    c30min=c[:,ntime÷2]
    c60min=c[:,ntime]
    concentrationplot= plot(r,[c10min,c30min,c60min],
    xlabel="Dimensionless Distance",
    ylabel="Dimensionless Concentration",
    title="Concentration vs Radial Postion",
    label=["10min" "30min" "1hr"],
    legend=:topleft,
    reuse=false
    )
    return concentrationplot
end
function doubleconcplot(r,c1,c2,custlabel)
    n_spatial1,ntime1 = getn(c1)
    n_spatial2,ntime2 = getn(c2)
    c10min= c1[:,ntime÷6]
    c30min=c1[:,ntime÷2]
    c60min=c1[:,ntime]
    c10min2= c2[:,ntime2÷6]
    c30min2=c2[:,ntime2÷2]
    c60min2=c2[:,ntime2]

    concentrationplot= plot(r,[c10min,c30min,c60min,c10min2,c30min2,c60min2],
    xlabel="Dimensionless Distance",
    ylabel="Dimensionless Concentration",
    title="Concentration vs Radial Postion",
    label=custlabel,
    legend=:topleft,
    reuse=false
    )
    return concentrationplot
end
function errorconcplot(r,c1,c2)
    n_spatial1,ntime1 = getn(c1)
    n_spatial2,ntime2 = getn(c2)
    c10min= c1[:,ntime÷6]
    c30min=c1[:,ntime÷2]
    c60min=c1[:,ntime]
    c10min2= c2[:,ntime2÷6]
    c30min2=c2[:,ntime2÷2]
    c60min2=c2[:,ntime2]
    c10minerror=broadcast(abs,(c10min .- c10min2)./(c10min2))
    c30minerror=broadcast(abs,(c30min .- c30min2)./(c30min2))
    c60minerror=broadcast(abs,(c60min .- c60min2)./(c60min2))

    errorplot= plot(r,[c10minerror,c30minerror,c60minerror],
    xlabel="Dimensionless Distance",
    ylabel="ϵ",
    title="ϵ vs Radial Postion",
    label=["10min-error" "30min-error" "1hr-error" ],
    legend=:topleft,
    reuse=false
    )
    return errorplot
end
function singleaccumplot(c)
    n_spatial,n_time = getn(c)
    t,accum =Accumulation_Model(c,n_spatial,n_time)
    println("yeet",n_spatial," ",n_time)
    accumplot= plot(t,accum,
    xlabel="Time(min)",
    ylabel="Dimensionless Concentration Accumulation",
    title="Accumulation vs Time",
    legend=:topleft,
    reuse=false
    )
    return accumplot
end
function doubleaccumplot(c1,c2,custlabel)
    n_spatial1,n_time1 = getn(c1)
    n_spatial2,n_time2 = getn(c2)
    t,accum1=Accumulation_Model(c1,n_spatial1,n_time1)
    t,accum2=Accumulation_Model(c2,n_spatial2,n_time2)

    accumplot= plot(t,[accum1,accum2],
    xlabel="Time(min)",
    ylabel="Dimensionless Concentration Accumulation",
    title="Accumulation vs Time",
    label=custlabel,
    legend=:topleft,
    reuse=false
    )
    return accumplot
end
function Accumerrorplot(c1,c2)
    n_spatial1,n_time1 = getn(c1)
    n_spatial2,n_time2 = getn(c2)
    t,accum1=Accumulation_Model(c1,n_spatial1,n_time1)
    t,accum2=Accumulation_Model(c2,n_spatial2,n_time2)
    caccum_error=broadcast(abs,(accum1 .- accum2)./(accum1))
    accumeplot= plot(t,caccum_error,
    xlabel="Time(min)",
    ylabel="ϵ",
    title="Accumulation Error vs Time",
    label=["ϵ"],
    legend=:topleft,
    reuse=false
    )
    return accumeplot
end
function velocityplot(r,v)
    vplot = plot(r,v,
    xlabel="Dimensionless Radial Position",
    ylabel="Dimensionless Velocity",
    title="Velocity vs Radial Position",
    label="Velocity",
    legend=:topleft,
    reuse = false)
    return vplot
end

function pressureplot(r,pressure)
    pplot= plot(r,pressure,
    xlabel="Dimensionless Radial Position",
    ylabel="Dimensionless Pressure",
    title="Pressure vs Radial Position",
    label="Pressure",
    legend=:bottomleft
    )
    return pplot
end
