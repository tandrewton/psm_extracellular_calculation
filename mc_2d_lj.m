%% monte carlo simulation of 2D LJ discs with non-periodic rectangular boundaries
close all; clear;
set(0,'DefaultFigureWindowStyle','docked')
rng(1)
% N - number of particles
% T - reduced temperature
% Nsteps - number of steps
% maxdr - maximum particle displacement
% initialConfig - initial configuration of particles (2 by N matrix)
% initialU - initial energy of the configuration
% rCutoff - the cutoff distance for the energy
N = 100;
%radius = normrnd(1.0, 0.1, N, 1)/2; % polydisperse normal distribution 1,0.1
%assert(isempty(radius(radius <= 0)));
T = 0.5;
phi = 0.55;
Nsteps = 100000000;
radius = 0.5;
maxdr = radius;
rCutoff = 2.5;
L = sqrt(N*pi*radius^2/phi);
moveCount = 0;
numFrames = 1;
numCyclesPrint = 10;
isSaveFigure = true; % determines whether we save the figures we print
%numCyclesPrint = 0; % print every n MC cycles = N*n accepted moves
%minMoves = 15000;
minMoves = 15000;
maxMoves = minMoves + numCyclesPrint*N*numFrames;

movedParticle = 0;

% create initial configuration
initialConfig = (L-2*radius)*rand(2,N) + radius;
initialDistances = sqrt(bsxfun(@(x1,x2) (x1-x2).^2 ,...
       initialConfig(1,:),initialConfig(1,:)')...
       +bsxfun(@(x3,x4) (x4-x3).^2 ,...
       initialConfig(2,:),initialConfig(2,:)'));

dist = initialDistances;
%U = 4*sum(sum(dist.^(-12)-dist.^(-6),"omitnan"),"omitnan");
U = pairU(dist, rCutoff);
% use uniform size for first configuration, equilibrate using general
% equation for polydisperse radii
particlesPosition = initialConfig;
plotConfiguration(particlesPosition);
yline(0); yline(L); xline(0); xline(L);

dUList = [];
UList = [];

for step=1:Nsteps
    if (step == Nsteps)
        disp("error: step hit max counter, aborting")
        assert(step < Nsteps)
    end
    if moveCount >= maxMoves
        moveCount
        break;
    elseif moveCount == minMoves && ~isempty(dUList)
        % after enough moves, consider the initial configuration melted and
        % restart energy counting for better numerical stability
        U = pairU(dist, rCutoff);
        dUList = [];
        UList = [];
        plotConfiguration(particlesPosition);
        yline(0); yline(L); xline(0); xline(L);
        disp("Melted, restarting energy counter\n")
    end
    % begin trial move for Metropolis algorithm

    isMoveAccepted = false;

    % choose particle to move 
    movedParticle = mod(movedParticle, N)+1;

    % choose displacement:
    displacex = maxdr*rand - (maxdr/2);
    displacey = maxdr*rand - (maxdr/2);
    displace = sqrt(displacex^2 + displacey^2);

    % calculate trial move properties
    newParticlesPosition = particlesPosition;
    newParticlesPosition(:,movedParticle) = particlesPosition(:,movedParticle) + [displacex; displacey];
    newDist = reCalcDist(dist, movedParticle, newParticlesPosition, N);

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
            %plotConfiguration(particlesPosition);
            %yline(0); yline(L); xline(0); xline(L);
            dUList(end+1) = dU;
            UList(end+1) = U;
            isMoveAccepted = true;
        else
            % probability of keeping new state, otherwise reject move
            if rand < exp(-dU/T)
                U = U + dU;
                dist = newDist;
                particlesPosition = newParticlesPosition;
                moveCount = moveCount + 1;
                %plotConfiguration(particlesPosition);
                %yline(0); yline(L); xline(0); xline(L);
                dUList(end+1) = dU;
                UList(end+1) = U;
                isMoveAccepted = true;
            end
        end
    end

    if (isMoveAccepted && mod(moveCount, numCyclesPrint*N) == 0)
        %every nth MC cycle
        if (moveCount >= minMoves)
            % if past moveCount which is high enough to consider our system equilibrated,
            % plot and save the configuration
            filename = "mc_simulation_frames/" + ...
                "MC_cycle" + (moveCount-minMoves)/N/numCyclesPrint;
            exportFrameGraphic(filename, L, particlesPosition, isSaveFigure);
        end
        moveCount-minMoves
    end
end

finalU = U;
finalConfiguration = particlesPosition;
finalDistances = dist;
plotConfiguration(finalConfiguration);
yline(0); yline(L); xline(0); xline(L);

function newdist = reCalcDist(dist,movedParticles,...
                newParticlesPosition,N)
    for i = 1:length(movedParticles)
        movedP = movedParticles(i);
        % recalculates pair distances after moving a particle    
        xi = newParticlesPosition(1,movedP);
        yi = newParticlesPosition(2,movedP);
        newdist = dist;
        % recalculate the relevent row elements in dist matrix
        newdist(movedP,:) =...
        distNoPBC(xi,yi,newParticlesPosition);
        % recalculate the relevent column elements in dist matrix
        newdist(:,movedP) =...
        distNoPBC(xi,yi,newParticlesPosition);
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

% function U = pairU(dist,rCutoff)
%     % calculates the reduced energy according to the pair potential
%     % input: dist is a row vector of all pair distances
%     % output: U is the total energy 
%     dist_rCutoff = dist(dist < rCutoff);
%     invR = 1./dist_rCutoff;
%     u = 4*((invR.^12)-(invR.^6));  % pair energies, LJ
%     %u = 1e10 * (dist_rCutoff <= 1); % pair energies, hard sphere
%     U = sum(u,"omitnan");
% end

function U = pairU(dist, rCutoff)
    % calculates the reduced energy according to the pair potential
    % input: dist is a row vector of all pair distances
    % output: U is the total energy 
    epsilon_adh = 1;
    %dist_rCutoff = dist(dist < rCutoff);
    dist_hard_contact = dist(dist < 1.0); % r < 1
    dist_rCutoff = dist_hard_contact(dist_hard_contact < rCutoff); % 1 < r < 2.5

    uCut = (1-(1/2.5).^6).^2; % reference: virrueta, o'hern, regan Proteins 2016

    uHard = hardRepulsive(dist_hard_contact) - epsilon_adh * uCut;

    uAttractive = epsilon_adh*(hardRepulsive(dist_rCutoff) - uCut);

    U = sum(uHard+uAttractive);
end

function U = hardRepulsive(dist)
    U = (1 - (1./dist).^6).^2; % pair energy for repulsive hard disc
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

function h = circle2(x,y,r)
    d = r*2;
    px = x-r;
    py = y-r;
    h = rectangle('Position',[px py d d],'Curvature',[1,1], 'FaceColor', 'k');
    daspect([1,1,1])
end

function plotConfiguration(config, fignum)
    if exist('fignum','var')
        figure(fignum); clf; hold on; % specified figure
    else
        figure(); clf; hold on; % new figure
    end

    radius = 1/2;
    for ii=1:length(config(1,:))
        circle2(config(1, ii), config(2,ii), radius);
        %scatter(config(1,ii), config(2,ii), 'k')
    end
    axis equal;
    axis off;
end

function exportFrameGraphic(filename, L, particlesPosition, isPrint)
    plotConfiguration(particlesPosition, 2)
    if (isPrint)
        xlim([0,L]*1.1 - 0.05*L)
        ylim([0,L]*1.1 - 0.05*L)
        yline(0); yline(L); xline(0); xline(L);
        exportgraphics(gcf, filename+".tif", 'Resolution', 100);
    
        clf; hold on;
        axis equal; axis off;
        xlim([0,L]*1.1 - 0.05*L)
        ylim([0,L]*1.1 - 0.05*L)
        yline(0); yline(L); xline(0); xline(L);
        exportgraphics(gcf, filename+"_bd.tif", 'Resolution', 100);
    end
end
