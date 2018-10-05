module lidCavity
import IterativeSolvers
import SparseArrays
import Plots

export main
export plotFlow
export compute_Δt
export compute_Ψ
export compute_Ω
export compute_VELOCITY
export applyBC_Ω
export assembleCoeffMatrix
export assembleRHS

function main(N,tmax,method,tol,max_tstep,doPlot)
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

    Ω[[NY NY-1],:] = fill(-2.0 / Δy,2,NY);

    FDM = assembleCoeffMatrix(Δx,Δy,NX,NY);

    plotData = [];
    t = 0.;
    pIter = 100;
    pCount = 0;
    tstep = 0;
    while t < tmax && tstep < max_tstep # start the time integration
        pCount += 1;
        tstep += 1;
        # Compute timestep
        Δt = compute_Δt(VELOCITY, cREYNOLDS, Δt_max, Δx, Δy);
        # Compute streamfunction
        Ψ = compute_Ψ(Ψ, Ω, FDM, method, tol, Δx, Δy, NX, NY);
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

function compute_Δt(VELOCITY, cREYNOLDS, Δt_max, Δx, Δy);
    vMax = maximum(VELOCITY);
    crMax = maximum(cREYNOLDS);

    Δt1 = (Δx*Δy) / (4*vMax);
    Δt2 = (2*vMax) / crMax;
    Δt = minimum(filter(x->isfinite(x),[Δt1 Δt2 Δt_max]));

    return Δt
end

function compute_Ψ(Ψ, Ω, FDM, method, tol, Δx, Δy, NX, NY);
    b = assembleRHS(Ψ,Ω,NX,NY,Δx,Δy);
    x = reshape(Ψ,NX*NY,1);
    if lowercase(method) == "direct"
        x = FDM \ b;
    elseif lowercase(method) == "cg"
        IterativeSolvers.cg!(x,FDM,b;tol=tol);
    end
    Ψ = collect(reshape(transpose(x),NY,NX));
    return Ψ
end

function compute_Ω(Ψ, Ω, U, V, Δt, Δx, Δy, NX, NY, Visc);
    Ω0 = Ω;
    Threads.@threads for ii = 2:NX-1
        for jj = 2:NY-1
            Ω[jj,ii] = Ω0[jj,ii] + Δt * (-(U[jj,ii] * ((Ω0[jj,ii+1] - Ω0[jj,ii-1]) / (2.0*Δx))) + -(V[jj,ii] * ((Ω0[jj+1,ii] - Ω0[jj-1,ii]) / (2.0*Δy))) + Visc*(((Ω0[jj,ii-1] - 2.0*Ω0[jj,ii] + Ω0[jj,ii+1])/(Δx^2.)) + ((Ω0[jj-1,ii] - 2.0*Ω0[jj,ii] + Ω0[jj+1,ii])/(Δy^2.))));
        end
    end
    return Ω
end

function compute_VELOCITY(Ψ, Δx, Δy, NX, NY);
    U = zeros(Float64,NY,NX);
    V = zeros(Float64,NY,NX);
    VELOCITY = zeros(Float64,NY,NX);
    cREYNOLDS = zeros(Float64,NY,NX);

    Threads.@threads for ii = 2:NX-1
        for jj = 2:NY-1
            U[jj,ii] =  (Ψ[jj+1,ii] - Ψ[jj-1,ii]) / (2.0*Δy);
            V[jj,ii] = -(Ψ[jj,ii+1] - Ψ[jj,ii-1]) / (2.0*Δx);
            VELOCITY[jj,ii] = sqrt(U[jj,ii]^2. + V[jj,ii]^2.);
            cREYNOLDS[jj,ii] = abs(U[jj,ii]) + abs(V[jj,ii]);
        end
    end
    return U, V, VELOCITY, cREYNOLDS
end

function applyBC_Ω(Ψ, Ω, Δx, Δy, NX, NY);
    for ii = 2:NX-1;
        Ω[1,ii] = -2.0*Ψ[2,ii]/(Δy^2.); # vorticity on bottom wall
        Ω[NY,ii]= -2.0*Ψ[NY-1,ii]/(Δy^2) - 2.0/Δy; # vorticity on top wall
    end
    for jj = 2:NY-1;
        Ω[jj,1] = -2.0*Ψ[jj,2]/(Δx^2.);    # vorticity on left wall
        Ω[jj,NX]= -2.0*Ψ[jj,NX-1]/(Δx^2.); # vorticity on right wall
    end
    return Ω
end

function assembleCoeffMatrix(dx,dy,NX,NY)
    # Generate A for a domain from 1:N
    sRow = zeros(Int64, 5 * NY * NX)
    sCol = zeros(Int64, 5 * NY * NX)
    sVal = fill(NaN, 5 * NY * NX)
    sIdx = 0;

    I = 1;
    C = -2*(dy^2. + dx^2.); #-4;
    L = 1 * dy^2.; #1;
    R = 1 * dy^2.; #1;
    D = 1 * dx^2.; #1;
    U = 1 * dx^2.; #1;

    nNodes = (NX)*(NY);
    for node = 1:nNodes
        jj = Int64(rem(node-1, NY)) + 1;
        ii = Int64((node - jj)/(NY)) + 1;
        if ii == 1 && jj == 1
            # Bottom Left Node
            sIdx = sIdx[end] + 1;
            sRow[sIdx] = node;
            sCol[sIdx] = jj + (ii - 1)*(NY);
            sVal[sIdx] = I;
        elseif ii == NX && jj == 1
            # Bottom Right Node
            sIdx = sIdx[end] + 1;
            sRow[sIdx] = node;
            sCol[sIdx] = jj + (ii - 1)*(NY);
            sVal[sIdx] = I;
        elseif ii == 1 && jj == NY
            # Top Left Node
            sIdx = sIdx[end] + 1;
            sRow[sIdx] = node;
            sCol[sIdx] = jj + (ii - 1)*(NY);
            sVal[sIdx] = I;
        elseif ii == NX && jj == NY
            # Top Right Node
            sIdx = sIdx[end] + 1;
            sRow[sIdx] = node;
            sCol[sIdx] = jj + (ii - 1)*(NY);
            sVal[sIdx] = I;
        elseif ii == 1
            # Left Node
            sIdx = sIdx[end] + 1;
            sRow[sIdx] = node;
            sCol[sIdx] = jj + (ii - 1)*(NY);
            sVal[sIdx] = I;
        elseif ii == NX
            # Right Node
            sIdx = sIdx[end] + 1;
            sRow[sIdx] = node;
            sCol[sIdx] = jj + (ii - 1)*(NY);
            sVal[sIdx] = I;
        elseif jj == 1
            # Bottom Node
            sIdx = sIdx[end] + 1;
            sRow[sIdx] = node;
            sCol[sIdx] = jj + (ii - 1)*(NY);
            sVal[sIdx] = I;
        elseif jj == NY
            # Top Node
            sIdx = sIdx[end] + 1;
            sRow[sIdx] = node;
            sCol[sIdx] = jj + (ii - 1)*(NY);
            sVal[sIdx] = I;
        else
            # Interior Node
            sIdx = sIdx[end] + 1 : sIdx[end] + 5;
            sRow[sIdx] = [node, node, node, node, node];
            sCol[sIdx] = [jj, jj, jj, jj+1, jj-1] + ([ii, ii+1, ii-1, ii, ii] .- 1).*(NY);
            sVal[sIdx] = [C, R, L, U, D];
        end
    end
    filter!(x->x≠0,sRow);
    filter!(x->x≠0,sCol);
    filter!(x->!isnan(x),sVal);
    A = SparseArrays.sparse(sRow,sCol,sVal,(NX*NY),(NX*NY));
    return A
end

function assembleRHS(Ψ,Ω,NX,NY,Δx,Δy)
    ## Generate b for a domain from 2:N-1
    nNodes = NX * NY;
    b = zeros(Float64,nNodes,1);
    Threads.@threads for node = 1:nNodes
        jj = Int64(rem(node-1, NY)) + 1;
        ii = Int64((node - jj)/NY) + 1;
        bidx = jj + (ii - 1)*(NY);
        ## Check to see if node is an interior Node
        if ii≠1 && ii≠NY && jj≠1 && jj≠NY
            b[bidx] = -(Δx^2.0 * Δy^2.0 * Ω[bidx]);
        end
    end
    return b
end

end
