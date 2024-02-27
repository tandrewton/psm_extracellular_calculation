 %% for DP simulations
close all; clear; 
if ismac
    set(0,'DefaultFigureWindowStyle','docked')
    folder = "/Users/AndrewTon/Documents/YalePhD/projects/dpm/output/cells/psm/";
elseif ispc
    set(0,'DefaultFigureWindowStyle','docked')
    folder = "C:/Users/atata/projects/dpm/output/cells/psm/";
else
    % for cluster runs
    folder = "/gpfs/gibbs/pi/ohern/at965/dpm/psm/output/";
end

rng(1)
calA0 = "1.0";
k_off = "1.0";
ka_arr = ["5.0"];
kb_arr = ["0.01"];
v0_arr = ["0.1"];
att_arr = ["0.0" "0.001" "0.005" "0.01" "0.05" "0.1"];
att2_arr = ["0.0" "0.001" "0.005" "0.01" "0.05" "0.1"];
t_stress_arr = ["10000.0"];
gamma_arr = ["0" "0.001" "0.1"];
numSeeds = 12;
phi = "0.8";
Duration = "500";
N = "40";
startSeed = "1";

cellDiam = 50; % in pixels
windowLengthPixels = round(2*cellDiam);

for ka=ka_arr
    disp("loop over ka!")
    for kb=kb_arr
        % for the original image
        figure(1); 
        % for the histograms
        figure(2); hold on; tiledlayout(1, length(windowLengthPixels));
        
        for v0_ind=1:length(v0_arr)
            v0 = v0_arr(v0_ind);
            for t_stress=t_stress_arr
                for gamma=gamma_arr
                    for att=att_arr
                        for att2=att2_arr
                            windowPackingFractions = [];
                            windowN = [];
                            randWindowPackingFractions = [];
                            randWindowN = [];
                            for seed=1:numSeeds
                                firstFrame = 25; % make sure to adjust these as necessary
                                lastFrame = 75;
                                packingFractionPerSeed = [];
                                for ff=firstFrame:lastFrame
                                    %filename = "psmtest/psmpsmtest8_seed1fr"+ff;
                                    filename = "psm_calA0"+calA0+"_phi"+phi+"_tm"+t_stress+"_v0"+v0+"_t_abp1.0_gamma"+gamma+"_kl1.0_ka"+ka+"_kb"+kb+"/"...
                                        + "psm_N"+N+"_dur"+Duration+"_att"+att+"_att2"+att2+"_sd"+string(seed)+"fr"+ff;
                                    raw_file = filename+".tif";
                                    bd_file = filename+"_bd.tif";
                                    try
                                        I_color = imread(folder+raw_file);
                                    catch ME
                                        disp("failed to read file: "+folder+raw_file);
                                        continue;
                                    end
                                    I = imbinarize(im2gray(I_color), 0.8);
                                    %I = imdilate(I, strel("disk",1)); % let edges bleed by 1 pixel
                                    psm_mask = ~imbinarize(im2gray(imread(folder+bd_file)));
                                    if (ff == lastFrame)
                                        figure(1); hold on;
                                        imshow(I_color);
                                        figure(2); hold on;
                                        imshow(I);
                                    end
                                    % find indices of largest area blob in psm_mask
                                    psm_mask = bwareafilt(psm_mask,1);
                                
                            %         [pct, n] = windowImageCountsFractional(I, windowLengthPixels, psm_mask);
                            %         windowPackingFractions = [windowPackingFractions (1-pct)];
                            %         windowN = [windowN n];
                                    % pct = fraction of positive pixels in window
                                    % n = measured window size as a fraction of max size
                                    
                                    [pct, n] = randWindowImageCountsFractional(I, windowLengthPixels, psm_mask, 50);
                                    randWindowPackingFractions = [randWindowPackingFractions (1-pct)];
                                    packingFractionPerSeed = [packingFractionPerSeed (1-pct)];
                                    randWindowN = [randWindowN n];
                                end
                                % write (seed, randWindowN', v0, att, k_ecm) to file
                
                                % to speed up the writing process, I could move this out of
                                % the seed loop, write randWindowN' instead of
                                % packingFractionPerSeed', and then use rowSeed instead of repmat(seed, rowDim)
                                % where rowSeed = []; rowSeed.append(repmat(seed,
                                % packingFractionPerSeed')) for each seed
                                rowDim = size(packingFractionPerSeed');
                                A = table(packingFractionPerSeed', repmat(att, rowDim), ...  
                                    repmat(v0, rowDim), repmat(att2, rowDim), ...
                                    repmat(t_stress, rowDim), repmat(gamma, rowDim), ...
                                    repmat(seed, rowDim));
                    
                                A.Properties.VariableNames = ["phi", "att", "v0", "att2", "tm", "gamma", "seed"];
                                writetable(A, 'windowedPhiDataFrame_calA'+calA0+'_phi'+phi+'.txt', 'WriteMode', 'append');
                            end
                
                            binedges = 0:0.02:1;
                        
                        %     % transform to weighted bin counts
                        %     [b,~] = discretize(windowPackingFractions, binedges);
                        %     c = accumarray(b', windowN');
                        %     windowBinCounts=[c;zeros(numel(binedges)-1-numel(c),1)];
                        
                            [b,~] = discretize(randWindowPackingFractions, binedges); % get bincounts a for each bin, hist inds b for each element
                            c = accumarray(b', randWindowN'); % into inds b, accumulate weights N
                            randWindowBinCounts = [c;zeros(numel(binedges)-1-numel(c),1)]; % zero pad to right size
                        
                            figure(3+v0_ind);
                            nexttile; hold on;
                        
                            %histogram('binedges', binedges, 'bincounts', windowBinCounts','normalization', 'pdf');
                            histogram('binedges', binedges, 'bincounts', randWindowBinCounts','normalization', 'pdf');
                            text(0.4, 12, "a="+att+", att2="+att2)
                                
                            %histogram(windowPackingFractions, binedges, 'normalization', 'pdf')
                            %histogram(randWindowPackingFractions, binedges, 'normalization', 'pdf')
                            xlim([0 1])
                            xticks([0, 0.4, 0.8, 1])
                            %xlabel("\phi")
                            %ylabel("P(\phi)")
                            set(gca,"Fontsize", 10)
                            axis square
                            %title(string(windowLengthPixels/cellDiam))
                        end
                    end
                end
            end
            figure(3+v0_ind); linkaxes
            ylim([0 18.0])
        end
    end
end
disp("done!")
