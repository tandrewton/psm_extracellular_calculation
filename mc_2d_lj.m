%% monte carlo simulation of 2D LJ discs with non-periodic rectangular boundaries
close all; clear;
rng(1)
% N - number of particles
% T - reduced temperature
% Nsteps - number of steps
% maxdr - maximum particle displacement
% initialConfig - initial configuration of particles (2 by N matrix)
% initialU - initial energy of the configuration
% rCutoff - the cutoff distance for the energy
N = 50;
T = 0.1;
rho = 0.55;
Nsteps = 1000;
radius = 0.5;
maxdr = 2*radius;
rCutoff = 2.5;
L = sqrt(N/rho);
moveCount = 0;
maxMoves = 25;
movedParticle = 0;

% create initial configuration
initialConfig = (L-2*radius)*rand(2,N) + radius;
initialDistances = sqrt(bsxfun(@(x1,x2) (x1-x2).^2 ,...
       initialConfig(1,:),initialConfig(1,:)')...
       +bsxfun(@(x3,x4) (x4-x3).^2 ,...
       initialConfig(2,:),initialConfig(2,:)'));

dist = initialDistances;
U = 4*sum(sum(dist.^(-12)-dist.^(-6),"omitnan"),"omitnan");
particlesPosition = initialConfig;

dUList = [];

for step=1:Nsteps
    if (moveCount > maxMoves)
        break;
    end
    % begin trial move for Metropolis algorithm

    % choose particle to move 
    movedParticle = movedParticle + 1;
    if movedParticle == N + 1
        movedParticle = 1;
    end
    % choose displacement:
    displacex = maxdr*rand - (maxdr/2);
    displacey = maxdr*rand - (maxdr/2);
    displace = sqrt(displacex^2 + displacey^2);

    % calculate trial move properties
    newParticlesPosition = particlesPosition;
    newParticlesPosition(:,movedParticle) = particlesPosition(:,movedParticle) + [displacex; displacey];
    newDist = reCalcDist(dist, movedParticle, newParticlesPosition, N,[]);

    % calculate change in energy
    dU = Uchange(movedParticle, dist, newDist, N, rCutoff);
    dU = dU + wallEnergy(movedParticle, newParticlesPosition, L, radius);

    % calculate probability of accepting move
    if dU < 75 * T % trial move is larger than lower bound on move probability
        if dU < 0 % automatically accept energy reducing move
            U = U + dU;
            dist = newDist;
            particlesPosition = newParticlesPosition;
            moveCount = moveCount + 1;
            plotConfiguration(particlesPosition);
            yline(0); yline(L); xline(0); xline(L);
            dUList(end+1) = dU;
        else
            % probability of keeping new state corresponds to a boltzmann factor
            % otherwise reject move, keep original configuration
            if rand < exp(-dU/T)
                U = U + dU;
                dist = newDist;
                particlesPosition = newParticlesPosition;
                moveCount = moveCount + 1;
                plotConfiguration(particlesPosition);
                yline(0); yline(L); xline(0); xline(L);
                dUList(end+1) = dU;
            end
        end
    end
end

finalU = U;
finalConfiguration = particlesPosition;
finalDistances = dist;
plotConfiguration(finalConfiguration);
yline(0); yline(L); xline(0); xline(L);

function newdist = reCalcDist(dist,movedParticles,...
                newParticlesPosition,N,nlist)
    for i = 1:length(movedParticles)
        movedP = movedParticles(i);
        % recalculates pair distances after moving a particle    
        xi = newParticlesPosition(1,movedP);
        yi = newParticlesPosition(2,movedP);
        newdist = dist;
        % recalculate the relevent row elements in dist matrix
        if movedP > 1
            if ~isempty(nlist)
                neiInd =...
                    nlist.neighborsindy(nlist.neighborsindx == movedP);
                newdist(movedP,neiInd) =...
                distNoPBC(xi,yi,newParticlesPosition(:,neiInd));
            else
                newdist(movedP,1:(movedP-1)) =...
                    distNoPBC(xi,yi,newParticlesPosition(:,1:(movedP-1)));
            end
        end
        % recalculate the relevent column elements in dist matrix
        if movedP < N
            if ~isempty(nlist)
                neiInd =...
                    nlist.neighborsindx(nlist.neighborsindy == movedP);
                newdist(neiInd,movedP) =...
                distNoPBC(xi,yi,...
                    newParticlesPosition(:,neiInd));
            else
                newdist((movedP + 1):N,movedP) =...
                distNoPBC(xi,yi,newParticlesPosition(:,(movedP+1):N));
            end
        end
    end   
end

function dist = distNoPBC(x,y,allPositions)
    distx = abs(x - allPositions(1,:));
    disty = abs(y - allPositions(2,:));
    dist = sqrt(distx.^2 + disty.^2);
end

function dU = Uchange(movedParticle,dist,newDist,N,rCutoff)
    % calculates the change in energy after a particle has moved
    % calculate the old energy for the relevant particle pairs
    
    if movedParticle > 1
        oldUrow = ...
            pairU(dist(movedParticle,1:(movedParticle - 1)),rCutoff);
    else 
        oldUrow = 0;
    end
    
    if movedParticle < N
        oldUcol = ...
            pairU(dist((movedParticle + 1):N,movedParticle),rCutoff);
    else 
        oldUcol = 0;
    end
    oldU = oldUrow + oldUcol;
    
    % calculate the new energy for the relevant particle pairs
    
    if movedParticle > 1
        newUrow = pairU(newDist...
                (movedParticle,1:(movedParticle - 1)),rCutoff);
    else 
        newUrow = 0;
    end
    
    if movedParticle < N
        newUcol = pairU(newDist...
            ((movedParticle + 1):N,movedParticle),rCutoff);
    else 
        newUcol = 0;
    end
    
    newU = newUrow + newUcol;
    
    % calculate change in energy
    dU = newU - oldU;
end

function U = pairU(dist,rCutoff)
    % calculates the reduced energy according to the pair
    % potantial, only pair closer than rCutoff are regarded.
    % input: dist is a row vector of all pair distances
    % output: U is the total energy 
    dist_rCutoff = dist(dist < rCutoff);
    % u is the energies of each pair
    invR = 1./dist_rCutoff;
    u = 4*((invR.^12)-(invR.^6)); 
    U = sum(u);
end

function wallU = wallEnergy(particle, pos, L, radius)
    if (pos(1,particle) < radius || pos(1,particle) > L-radius)
        wallU = inf;
        %fprintf("(x) = %f\n", pos(1,particle))
        %fprintf("%d %d\n", pos(1,particle)<radius, pos(1,particle)>L-radius)
        %fprintf("%f %f\n", radius, L-radius)
    elseif (pos(2,particle) < radius || pos(2,particle) > L-radius)
        wallU = inf;
        %fprintf("(y) = %f\n", pos(2,particle))
        %fprintf("%d %d\n", pos(2,particle)<radius, pos(2,particle)>L-radius)
    else 
        wallU = 0;
    end
end

function plotConfiguration(config)
    figure(); clf; hold on;
    radius = 1/2;
    for ii=1:length(config(1,:))
        circle2(config(1, ii), config(2,ii), radius);
    end
    axis equal;
    axis off;
end

function h = circle2(x,y,r)
    d = r*2;
    px = x-r;
    py = y-r;
    h = rectangle('Position',[px py d d],'Curvature',[1,1], 'FaceColor', 'k');
    daspect([1,1,1])
end