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
ka = "2.5";
kl_arr = ["0.2"];   
%kl_arr = ["1.0"];
kb_arr = ["0.01"];
v0_arr = ["0.1"];
<<<<<<< Updated upstream
%att_arr = ["0.001" "0.005" "0.01" "0.03" "0.05" "0.1"];
%att2_arr = ["0.0005" "0.001" "0.005" "0.01" "0.05" "0.1"];
att_arr = ["0.015" "0.02" "0.025" "0.03", "0.035"];
att2_arr = ["0" "0.01" "0.015" "0.02" "0.025"];
=======
att_arr = ["0.001" "0.005" "0.01" "0.03" "0.05" "0.1"];
att2_arr = ["0.0005" "0.001" "0.005" "0.01" "0.05" "0.1"];
%att2_arr = ["0.0"];
>>>>>>> Stashed changes
t_stress_arr = ["10000.0"];
gamma_arr = ["0"];
kon_arr=["1.0"];
koff_arr=["0.1"];
%kecm_arr=["0.001", "0.005", "0.01", "0.05"];
numSeeds = 20;
phi = "0.8"; 
Duration = "300";
N = "40";
startSeed = "1";

cellDiam = 75; % in pixels
windowLengthPixels = round(2*cellDiam);

for kl=kl_arr
    disp("loop over kl!")
    for kb=kb_arr
        % for the original image
        figure(1); 
        % for the histograms
        figure(2); hold on; tiledlayout(1, length(windowLengthPixels));
        
        for v0_ind=1:length(v0_arr)
            v0 = v0_arr(v0_ind);
            for t_stress=t_stress_arr
                for gamma=gamma_arr
                    for k_on=kon_arr
                        for k_off=koff_arr
                            %for k_ecm=kecm_arr
                            for att=att_arr
                                for att2=att2_arr
                                    k_ecm=att2;
                                    windowPackingFractions = [];
                                    windowN = [];
                                    randWindowPackingFractions = [];
                                    randWindowN = [];
                                    for seed=1:numSeeds
                                        firstFrame = 25; % make sure to adjust these as necessary
                                        lastFrame = 60;
                                        packingFractionPerSeed = [];
                                        for ff=firstFrame:lastFrame
                                            %filename = "psmtest/psmpsmtest8_seed1fr"+ff;
                                            filename = "psm_calA0"+calA0+"_phi"+phi+"_tm"+t_stress+"_v0"+v0+"_t_abp1.0_gamma"+gamma+...
                                                "_k_on_"+k_on+"_k_off_"+k_off+"_k_ecm_"+k_ecm+"_kl"+kl+"_ka"+ka+"_kb"+kb+"/"...
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
                                            repmat(kl, rowDim), repmat(gamma, rowDim), ...
                                            repmat(k_on, rowDim), repmat(k_off, rowDim), ...
                                            repmat(k_ecm, rowDim), repmat(seed, rowDim));
                            
                                        A.Properties.VariableNames = ["phi", "att", "v0", "att2", "kl", "gamma", "k_on", "k_off", "k_ecm", "seed"];
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
                            %end
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
