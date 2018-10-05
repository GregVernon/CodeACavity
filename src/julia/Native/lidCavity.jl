module lidCavity
import Plots

export main
export plotFlow

function main(N,tmax,tol,max_tstep,doPlot)
    NX = N;
    xmin = 0.;
    xmax = 1.;
    x = range(xmin,stop=xmax,length=NX);
    Δx = x[2] - x[1];

    NY = N;
    ymin = 0.;
    ymax = 1.;
    y = range(ymin,stop=ymax,length=NY);
    Δy = y[2] - y[1];

    Visc = 0.01;
    Δt_max = tmax / 10000.;

    Ψ = zeros(Float64, NY, NX);
    Ω = zeros(Float64, NY, NX);
    U = zeros(Float64, NY, NX);
    V = zeros(Float64, NY, NX);
    VELOCITY = zeros(Float64, NY, NX);
    cREYNOLDS = zeros(Float64, NY, NX);

    for ii = 1:NX
        Ω[NY,ii]   = -2.0 / Δy;
        Ω[NY-1,ii] = -2.0 / Δy;
    end


    plotData = [];
    t = 0.;
    pIter = 100;
    pCount = 0;
    tstep = 0;
    while t < tmax && tstep < max_tstep # start the time integration
        pCount += 1;
        tstep += 1;
        # Compute timestep
        Δt = compute_Δt(VELOCITY, cREYNOLDS, Δt_max, Δx, Δy, NX, NY);
        # Compute streamfunction
        Ψ = compute_Ψ(Ψ, Ω, tol, Δx, Δy, NX, NY);
        # Apply vorticity boundary conditions
        Ω = applyBC_Ω(Ψ, Ω, Δx, Δy, NX, NY);
        # Compute Velocity
        U, V, VELOCITY, cREYNOLDS = compute_VELOCITY(Ψ, Δx, Δy, NX, NY);
        # Compute vorticity
        Ω = compute_Ω(Ψ, Ω, U, V, Δt, Δx, Δy, NX, NY, Visc);
        # Increment time value by timestep
        t += Δt;
        # Plot
        if doPlot == true
            if pCount == pIter || tstep == 1
                pCount = 0;
                println("Time: ", t)
                plt = plotFlow(Ψ, Ω, U, V, VELOCITY, cREYNOLDS, x, y, plotData);
            end
        end
    end

    return Ψ,Ω,U,V
end

function plotFlow(Ψ, Ω, U, V, VELOCITY, cREYNOLDS, x, y, plt)
    if isempty(plt) == true
        # Initialize plots
        Plots.pyplot(); # Set PyPlot as backend
        plt = Plots.heatmap(Plots.heatmap(Ψ), Plots.heatmap(Ω), Plots.heatmap(VELOCITY))
        Plots.gui(plt)
    else
        Plots.heatmap!(plt,Ψ,Ω,VELOCITY)
    end
    return plt
end

function compute_Δt(VELOCITY, cREYNOLDS, Δt_max, Δx, Δy, NX, NY);
    vMax = 0.;
    crMax = 0.;
    for ii = 1:NX
        for jj = 1:NY
            if VELOCITY[jj,ii] > vMax
                vMax = VELOCITY[jj,ii]
            end
            if cREYNOLDS[jj,ii] > crMax
                crMax = cREYNOLDS[jj,ii]
            end
        end
    end

    Δt1 = (Δx*Δy) / (4*vMax);
    Δt2 = (2*vMax) / crMax;
    if isfinite(Δt1)
        Δt = Δt1;
        if isfinite(Δt2) && Δt2 < Δt
            Δt = Δt2;
        end
        if Δt_max < Δt
            Δt = Δt_max;
        end
    elseif isfinite(Δt2)
        Δt = Δt2;
        if Δt_max < Δt
            Δt = Δt_max;
        end
    else
        Δt = Δt_max;
    end
    return Δt
end

function compute_Ψ(Ψ, Ω, tol, Δx, Δy, NX, NY);
    err = Inf;
    while err > tol
        err = 0;
        Ψ0 = Ψ;
        for ii=2:NX-1
            for jj=2:NY-1
                Ψ[jj,ii]=0.25*(Ψ[jj+1,ii]+Ψ[jj-1,ii]+Ψ[jj,ii+1]+Ψ[jj,ii-1]+Δx*Δy*Ω[jj,ii]);
                err += natjul_abs(Ψ0[jj,ii] - Ψ[jj,ii]) / (NX*NY);
            end
        end
    end
    return Ψ
end

function compute_Ω(Ψ, Ω, U, V, Δt, Δx, Δy, NX, NY, Visc);
    Ω0 = Ω;
    for ii = 2:NX-1
        for jj = 2:NY-1
            Ω[jj,ii] = Ω0[jj,ii] + Δt * (-(U[jj,ii] * ((Ω0[jj,ii+1] - Ω0[jj,ii-1]) / (2*Δx))) + -(V[jj,ii] * ((Ω0[jj+1,ii] - Ω0[jj-1,ii]) / (2*Δy))) + Visc*(((Ω0[jj,ii-1] - 2*Ω0[jj,ii] + Ω0[jj,ii+1])/(Δx*Δx)) + ((Ω0[jj-1,ii] - 2*Ω0[jj,ii] + Ω0[jj+1,ii])/(Δy*Δy))));
        end
    end
    return Ω
end

function compute_VELOCITY(Ψ, Δx, Δy, NX, NY);
    U = zeros(Float64,NY,NX);
    V = zeros(Float64,NY,NX);
    VELOCITY = zeros(Float64,NY,NX);
    cREYNOLDS = zeros(Float64,NY,NX);

    for ii = 2:NX-1
        for jj = 2:NY-1
            U[jj,ii] =  (Ψ[jj+1,ii] - Ψ[jj-1,ii]) / (2*Δy);
            V[jj,ii] = -(Ψ[jj,ii+1] - Ψ[jj,ii-1]) / (2*Δx);
            VELOCITY[jj,ii] = natjul_sqrt(U[jj,ii]*U[jj,ii] + V[jj,ii]*V[jj,ii]);
            cREYNOLDS[jj,ii] = natjul_abs(U[jj,ii]) + abs(V[jj,ii]);
        end
    end
    return U, V, VELOCITY, cREYNOLDS
end

function applyBC_Ω(Ψ, Ω, Δx, Δy, NX, NY);
    for ii = 2:NX-1;
        Ω[1,ii] = -2.0*Ψ[2,ii]/(Δy*Δy); # vorticity on bottom wall
        Ω[NY,ii]= -2.0*Ψ[NY-1,ii]/(Δy*Δy) - 2.0/Δy; # vorticity on top wall
    end
    for jj = 2:NY-1;
        Ω[jj,1] = -2.0*Ψ[jj,2]/(Δx*Δx);    # vorticity on left wall
        Ω[jj,NX]= -2.0*Ψ[jj,NX-1]/(Δx*Δx); # vorticity on right wall
    end
    return Ω
end

function natjul_abs(v)
    v = v * ((v>0.) - (v<0.));
    return v
end

function natjul_sqrt(y0)
    tol = 1e-14;
    if y0 == 0
        x = 0.;
        return x
    elseif y0 < 1
        x = 1;
    else
        x = y0;
    end

    err = x * x;
    while err >= tol
        x = x - (x * x - y0) / (2 * x);
        err = natjul_abs((x * x) - y0) / y0;
    end
    return x
end

end
