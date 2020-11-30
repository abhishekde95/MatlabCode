function [plot_counter,nbins,nbins1] = compute_projection_val(plot_counter,mode,subunits,nbins2,nbins3,flag,ch)

% Computes projection values on the basis vectors derived from the subunits
out = STPROJmod('return');
projs = out{1};
Lspike = out{2};
clear STPROJmod;
clear out;

if (flag == 1)
    labels = ['STA'; 'PC1'; 'PC2'; 'PC3'; 'PCi'];
    X_label = labels(ch(1),:);
    Y_label = labels(ch(2),:);
else
    X_label = 'S1 proj';
    Y_label = 'S2 proj';
end
if (mod(mode,3)==2)
    min_val = floor(min(projs(:))*100)/100;
    max_val = ceil(max(projs(:))*100)/100;
    num_bins = 15;
    bin_interval = (max_val-min_val)/num_bins; 
    nbins1 = min_val:bin_interval:max_val;
    % choose the bins such that they lie between 5 and 95 percentile
    nbins = linspace(prctile(projs(:),5), prctile(projs(:),95),numel(nbins1)+2);
    if (size(projs,2) == 2)
        figure(plot_counter); scatter(projs(:,1),projs(:,2)); hold on;
        scatter(projs(Lspike>0,1),projs(Lspike>0,2),'r'); 
        xlabel(X_label); ylabel(Y_label);       
        grid on;
        axis([min(projs(:,1))-0.1 max(projs(:,1))+0.1 min(projs(:,2))-0.1 max(projs(:,2))+0.1]);
        hold off;
        plot_counter = plot_counter + 1;
    end
else % in the current siuation if mode == 1
    nbins1 = nbins3;
    nbins = nbins2;
end
xmin = min(nbins); xmax = max(nbins);

if (size(projs,2) == 3)
    % using multivariate histogram if u have more than 2 subunits
    % Not perfect yet. Need to fix it.
    [count_raw,~,~,~] = histcn(projs,nbins,nbins,nbins);
    [count_st,~,mid_st,~] = histcn(projs(Lspike>0,:), nbins,nbins,nbins);
    count_raw(count_raw == 0) = 1;
    non_lin = count_st./count_raw;
    
    figure(plot_counter);
    subplot(1,subunits,1); Z = squeeze(sum(non_lin,3)); Z = Z/max(Z(:)); surfl(Z);
    xlabel('S2 proj'), ylabel('S1 proj'), zlabel('Non-linearity(S1-S2)');
    subplot(1,subunits,2); Z = squeeze(sum(non_lin,2)); Z = Z/max(Z(:)); surfl(Z);
    xlabel('S3 proj'), ylabel('S1 proj'), zlabel('Non-linearity(S1-S3)');
    subplot(1,subunits,3); Z = squeeze(sum(non_lin,1)); Z = Z/max(Z(:)); surfl(Z);
    xlabel('S3 proj'), ylabel('S2 proj'), zlabel('Non-linearity(S2-S3)');
    plot_counter = plot_counter + 1;
    
elseif (size(projs,2) == 2)
    % The program will run this code block if there are 2 subunits.
    % Converting the histogram into probability distribution and Estimating the non-linearity
    new_nbins = [nbins(1)-mean(diff(nbins)), nbins, mean(diff(nbins)) + nbins(end)];
    [n_spike,out_spike] = hist3([projs(Lspike>0,1), projs(Lspike>0,2)],{new_nbins,new_nbins});
    [n_raw,out_raw] = hist3([projs(:,1), projs(:,2)],{new_nbins,new_nbins});
    n_spike = n_spike(2:end-1,2:end-1);
    n_raw = n_raw(2:end-1,2:end-1);
    
    figure(plot_counter);
    % left hand side image for the spike trigerred ensemble image
    subplot(221); hist3([projs(Lspike>0,1), projs(Lspike>0,2)],{nbins1,nbins1}); xlabel('S2 Proj'); ylabel('S1 Proj');
    title('Spike Triggered Ensemble'); set(gcf,'renderer','opengl'); set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
    subplot(222); hist3([projs(:,1), projs(:,2)],{nbins1,nbins1}); xlabel(X_label); ylabel(Y_label);
    title('Raw Ensemble'); set(gcf,'renderer','opengl'); set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
    subplot(224); imagesc([min(new_nbins) max(new_nbins)],[min(new_nbins) max(new_nbins)], n_raw); xlabel(X_label); ylabel(Y_label); title('Raw Ensemble');
    subplot(223); imagesc([min(new_nbins) max(new_nbins)],[min(new_nbins) max(new_nbins)], n_spike); xlabel(X_label); ylabel(Y_label); title('Spike Triggered Ensemble');
    plot_counter = plot_counter + 1;
    n_raw(n_raw==0) = 1; % You can either enter 'NaN'or '1' to avoid division by zero
    non_lin = n_spike./n_raw;
    non_lin_pa = padarray(non_lin,[2 2],'replicate'); % pad the array with the border elements
    filt = fspecial('gaussian',5,1.0); % building a gaussian filter
    non_lin_blurred = conv2(non_lin_pa,filt,'same'); % convolving a gaussian filter with the firing rate map
    non_lin_blurred = non_lin_blurred(3:end-2,3:end-2);
        
    figure(plot_counter);
    subplot(2,2,1), surfl(nbins',nbins,non_lin); set(gca,'xlim',[xmin xmax],'ylim',[xmin xmax],'zlim',[min(non_lin(:)) max(non_lin(:))]);
    xlabel(Y_label), ylabel(X_label); title('Firing  Rate');
    subplot(2,2,2), plot(nbins,sum(non_lin,1)/max(sum(non_lin,1)),'LineWidth',2);   % plotting the normalised marginal non-linearity
    xlabel(Y_label), ylabel('Firing Rate'); title('Projection onto S2 vector');
    axis([nbins(1) nbins(end) 0 1]);
    subplot(2,2,3),plot(nbins,sum(non_lin,2)/max(sum(non_lin,2)),'LineWidth',2);
    ylabel('Firing Rate'), xlabel(X_label); title('Projection onto S1 vector');
    axis([nbins(1) nbins(end) 0 1]); %view([-90 90]);
    subplot(2,2,4), imagesc([xmin xmax],[xmin xmax],non_lin); 
    xlabel(Y_label), ylabel(X_label); title('2-D Firing  Rate Intensity Plot');
    plot_counter = plot_counter + 1;
    
    figure(plot_counter);
    subplot(231), imagesc([xmin xmax],[xmin xmax],non_lin);
    xlabel(Y_label), ylabel(X_label); title('FR w/o the filter');
    subplot(232), imagesc([xmin xmax],[xmin xmax],non_lin_blurred); 
    xlabel(Y_label), ylabel(X_label); title('FR after the filter');
    subplot(234), surfl(nbins',nbins,non_lin); set(gca,'xlim',[xmin xmax],'ylim',[xmin xmax],'zlim',[min(non_lin(:)) max(non_lin(:))]);
    xlabel(Y_label), ylabel(X_label); 
    subplot(235), surfl(nbins',nbins,non_lin_blurred); set(gca,'xlim',[xmin xmax],'ylim',[xmin xmax],'zlim',[min(non_lin_blurred(:)) max(non_lin_blurred(:))]);
    xlabel(Y_label), ylabel(X_label);
    subplot(233),contour(nbins',nbins,non_lin_blurred);set(gca,'YDir','Reverse');
    xlabel(Y_label), ylabel(X_label); title('Contour plot');
    plot_counter = plot_counter + 1;
    
elseif (size(projs,2) == 1)
    
    % The program will run this code block if there is only 1 subunit.
    figure(plot_counter);
    % left hand side image for the spike trigerred ensemble image
    subplot(1,2,1); hist(projs(Lspike>0),nbins1);
    hold on;
    xlabel('Proj values'); ylabel('frequency'); title('Spike Triggered Ensemble'); hold off;
    subplot(1,2,2);hist(projs ,nbins1)
    hold on;
    xlabel('Proj values'); ylabel('frequency'); title('Raw Ensemble'); hold off;
    plot_counter = plot_counter + 1;
    
    % Converting the histogram into probability distribution and Estimating the non-linearity
    [n_spike,out_spike] = hist(projs(Lspike>0),nbins);
    [n_raw,out_raw] = hist(projs,nbins);
    n_raw(n_raw==0) = 1; % to avoid division by zero
    non_lin = n_spike./n_raw;
    
    figure(plot_counter);
    plot(nbins,sum(non_lin,1)/max(sum(non_lin,1)),'LineWidth',2);   % plotting the normalised marginal non-linearity
    xlabel('Proj values'), ylabel('Firing Rate');
    axis([nbins(1) nbins(end) 0 1]);
    plot_counter = plot_counter + 1;
end

%*****************************************************************************
%             Trying out the fminsearch to check how well the data fits a
%             linear model. Also calculating the goodness of fit value to
%             check how well the model fits to the experimental data.
%*****************************************************************************
% Change this value if u wanna fit the model
if (mod(mode,3) == 2)
    prompt = ['Which kind of fitting would you prefer? Choose one of the options if there are more than 1 subunit:(1,2,3,4,5,6,7)\n '...
        '1) Fit the firing rate surface to a weighted sum of the S1 and S2 projection values\n '...
        '2) Fit the firing rate surface to the dot product of the marginal S1 and S2 firing rates\n '...
        '3) Fit the firing rate surface to a weighted sum of marginal S1 and S2 firing rates\n '...
        '4) Fit the firing rate surface to the dot product of the "estimated" marginal S1 and S2 firing rates\n '...
        '5) Predict the firing rate surface from the projection values onto the vector pointing towards the STA\n '...
        '6) Predict the firing rate surface from the "estimated" marginals of the projection values onto the vector pointing towards the STA\n '...
        '7) Estimate the projection vector using MLE and then perform the same analysis as 6)\n'];
    str = input(prompt,'s');
    
    if (size(projs,2)==1)
        % fitting procedure when there is only subunit present
        if (strcmp('1',str))
            x1 = nbins;
            myfun = @(w)sum((non_lin-w*x1).^2);
            [w_val,~] = fminsearch(myfun,[-1.2]);
            model_fit = w_val*x1;
            r_value = corrcoef(non_lin(:),model_fit(:));
            fprintf('Goodness of fit value for subunit-stimulus response is %d \n',r_value(1,2));
        else
        end
        
    elseif(size(projs,2)==2)
        % fitting implemented when 2 subunits are present
        % Implementing 5 different kinds of fitting procedure
        if (strcmp('1',str))
            % Fitting a exponential-linear weighted sum model to the firing rate surface
            [x1,x2] = meshgrid(nbins',nbins);
            myfun = @(w)sum(sum((non_lin-exp(w(1)*x1+w(2)*x2)).^2));
            [w_val,~] = fminsearch(myfun,[-1.2, 1]);
            model_fit = exp(w_val(1)*x1 + w_val(2)*x2);
            %             model_fit = model_fit*max(non_lin(:))/max(model_fit(:));
            r_value = corrcoef(non_lin(:),model_fit(:));
            fprintf('Goodness of fit value for subunit-stimulus response is %d \n',r_value(1,2));
            
            figure(plot_counter);
            subplot(221),surfl(x1,x2,non_lin); set(gca,'xlim',[xmin xmax],'ylim',[xmin xmax],'zlim',[min(non_lin(:)) max(non_lin(:))]);
            xlabel(Y_label), ylabel(X_label); title('Original firing rate surface');
            subplot(222),surfl(x1,x2,model_fit); set(gca,'xlim',[xmin xmax],'ylim',[xmin xmax],'zlim',[min(model_fit(:)) max(model_fit(:))])
            xlabel(Y_label), ylabel(X_label); title('Model-Predicted Firing Rate');
            subplot(224),surfl(x1,x2,non_lin-model_fit); set(gca,'xlim',[xmin xmax],'ylim',[xmin xmax],'zlim',[min(non_lin(:)-model_fit(:)) max(non_lin(:)-model_fit(:))]);
            xlabel(Y_label), ylabel(X_label); title('Residual Firing Rate Error');
            plot_counter = plot_counter + 1;
            
        elseif(strcmp('2',str))
            % Fitting the dot product of the marginal firing rates to the joint firing rate surface
            [x1,x2] = meshgrid(nbins',nbins);
            S1 = sum(non_lin,2)/max(sum(non_lin,2)); % projection onto the S1 vector
            S2 = sum(non_lin,1)/max(sum(non_lin,1)); % projection onto the S2 vector
            model_fit = S1*S2;
            %             model_fit = model_fit*max(non_lin(:))/max(model_fit(:));
            r_value = corrcoef(non_lin(:),model_fit(:));
            fprintf('Goodness of fit value for subunit-stimulus response is %d \n',r_value(1,2));
            
            figure(plot_counter);
            subplot(221),surfl(x1,x2,non_lin); set(gca,'xlim',[xmin xmax],'ylim',[xmin xmax],'zlim',[min(non_lin(:)) max(non_lin(:))]);
            xlabel(Y_label), ylabel(X_label); title('Original firing rate surface');
            subplot(222),surfl(x1,x2,model_fit); set(gca,'xlim',[xmin xmax],'ylim',[xmin xmax],'zlim',[min(model_fit(:)) max(model_fit(:))]);
            xlabel(Y_label), ylabel(X_label); title('Model-Predicted Firing Rate');
            subplot(223),plot(nbins',S1,'LineWidth',2); title('Vector 1'); set(gca,'xlim',[xmin xmax]);
            subplot(224),plot(nbins,S2,'LineWidth',2); title('Vector 2'); set(gca,'xlim',[xmin xmax]);
            plot_counter = plot_counter + 1;
            
            figure(plot_counter);
            imagesc([xmin xmax],[xmin xmax],non_lin);
            xlabel(Y_label), ylabel(X_label); title('2-D Firing  Rate Intensity Plot');
            plot_counter = plot_counter + 1;
            
            figure(plot_counter);
            surfl(x1,x2,non_lin-model_fit); set(gca,'xlim',[xmin xmax],'ylim',[xmin xmax],'zlim',[min(non_lin(:)-model_fit(:)) max(non_lin(:)-model_fit(:))]);
            xlabel(Y_label), ylabel(X_label); title('Residual Firing Rate Error');
            plot_counter = plot_counter + 1;
            
        elseif(strcmp('3',str))
            % Fitting a weighted sum of the marginal firing rate to the firing rate surface
            [x1,x2] = meshgrid(nbins',nbins);
            S1 = sum(non_lin,2)/max(sum(non_lin,2)); % projection onto the S1 vector
            S2 = sum(non_lin,1)/max(sum(non_lin,1)); % projection onto the S2 vector
            [S1_mesh,S2_mesh] = meshgrid(S1,S2);
            myfun = @(w)sum(sum((non_lin-w(1)*S1_mesh-w(2)*S2_mesh).^2));
            [w_val,~] = fminsearch(myfun,[-1.2, 1]);
            model_fit = w_val(1)*S1_mesh + w_val(2)*S2_mesh;
            %             model_fit = model_fit*max(non_lin(:))/max(model_fit(:));
            r_value = corrcoef(non_lin(:),model_fit(:));
            fprintf('Goodness of fit value for subunit-stimulus response is %d \n',r_value(1,2));
            
            figure(plot_counter);
            subplot(221),surfl(x1,x2,non_lin); set(gca,'xlim',[xmin xmax],'ylim',[xmin xmax],'zlim',[min(non_lin(:)) max(non_lin(:))]);
            xlabel(Y_label), ylabel(X_label); title('Original firing rate surface');
            subplot(222),surfl(x1,x2,model_fit); set(gca,'xlim',[xmin xmax],'ylim',[xmin xmax],'zlim',[min(model_fit(:)) max(model_fit(:))]);
            xlabel(Y_label), ylabel(X_label); title('Model-Predicted Firing Rate');
            subplot(223),plot(nbins',S1,'LineWidth',2); title('Vector 1'); set(gca,'xlim',[xmin xmax]);
            subplot(224),plot(nbins,S2,'LineWidth',2); title('Vector 2'); set(gca,'xlim',[xmin xmax]);
            plot_counter = plot_counter + 1;
            
            figure(plot_counter);
            surfl(x1,x2,non_lin-model_fit); set(gca,'xlim',[xmin xmax],'ylim',[xmin xmax],'zlim',[min(non_lin(:)-model_fit(:)) max(non_lin(:)-model_fit(:))]);
            xlabel(Y_label), ylabel(X_label); title('Residual Firing Rate Error');
            plot_counter = plot_counter + 1;
            
        elseif(strcmp('4',str))
            % Fitting a weighted sum of the 'estimated' marginal firing rate (obtained using
            % a polynomial fitting to the firing rate surface)
            [x1,x2] = meshgrid(nbins',nbins);
            S1 = sum(non_lin,2)/max(sum(non_lin,2)); % projection onto the S1 vector
            S2 = sum(non_lin,1)/max(sum(non_lin,1)); % projection onto the S2 vector
            myfun1 = @(w)sum((S1-w(1)*exp(nbins'*w(2))- w(3)).^2);
            [w1,~] = fminsearch(myfun1,[1,1,1]);
            myfun2 = @(u)sum((S2-u(1)*exp(nbins*u(2))- u(3)).^2);
            [w2,~] = fminsearch(myfun2,[1,1,1]);
            vec_est1 = w1(1)*exp(nbins'*w1(2)) + w1(3);
            vec_est2 = w2(1)*exp(nbins*w2(2)) + w2(3);
            model_fit = vec_est1 * vec_est2;
            %             model_fit = model_fit*max(non_lin(:))/max(model_fit(:));
            r_value = corrcoef(non_lin(:),model_fit(:));
            fprintf('Goodness of fit value for subunit-stimulus response is %d \n',r_value(1,2));
            
            figure(plot_counter);
            subplot(221),surfl(x1,x2,non_lin); set(gca,'xlim',[xmin xmax],'ylim',[xmin xmax],'zlim',[min(non_lin(:)) max(non_lin(:))]);
            xlabel(Y_label), ylabel(X_label); title('Original firing rate surface');
            subplot(222),surfl(x1,x2,model_fit); set(gca,'xlim',[xmin xmax],'ylim',[xmin xmax],'zlim',[min(model_fit(:)) max(model_fit(:))]);
            xlabel(Y_label), ylabel(X_label); title('Model-Predicted Firing Rate');
            subplot(223),plot(nbins',S1,'LineWidth',2); hold on; plot(nbins',vec_est1,'r','LineWidth',2); set(gca,'xlim',[xmin xmax]);
            legend('Original MP','Estimated MP'); title('Vector 1');
            subplot(224),plot(nbins,S2,'LineWidth',2); hold on; plot(nbins',vec_est2,'r','LineWidth',2); set(gca,'xlim',[xmin xmax]);
            legend('Original MP','Estimated MP'); title('Vector 2');
            plot_counter = plot_counter + 1;
            
            figure(plot_counter);
            surfl(x1,x2,non_lin-model_fit); set(gca,'xlim',[xmin xmax],'ylim',[xmin xmax],'zlim',[min(non_lin(:)-model_fit(:)) max(non_lin(:)-model_fit(:))]);
            xlabel(Y_label), ylabel(X_label); title('Residual Firing Rate Error');
            plot_counter = plot_counter + 1;
            
        elseif(strcmp('5',str))
            % Fitting procedure 5 - need to think of a name for it
            [spike_count,~] = hist3([projs(Lspike>0,1), projs(Lspike>0,2)],{nbins1,nbins1});
            [raw_count,~] = hist3([projs(:,1), projs(:,2)],{nbins1,nbins1});
            p1 = mean(projs(Lspike>0,1)) - mean(projs(:,1)); % based on the actual projection values
            p2 = mean(projs(Lspike>0,2)) - mean(projs(:,2));
            vec = [p1 p2]; vec = vec/sqrt(vec*vec'); % defining and normalizing the vector here
            spike_proj = [projs(Lspike>0,1) projs(Lspike>0,2)]*vec'; % dot product of the spike triggered data onto the vector
            raw_proj = projs*vec'; % dot product of the raw data onto the vector
            ub = ceil(max([max(spike_proj) max(raw_proj)]));
            lb = floor(min([min(spike_proj) min(raw_proj)]));
            hist_bins = lb:mean(diff(nbins)):ub;
            [spike_hist,~] = hist(spike_proj,hist_bins);
            [raw_hist,~] = hist(raw_proj,hist_bins); raw_hist(raw_hist==0) = 1;
            non_lin_1D = spike_hist./raw_hist;
            [X,Y] = meshgrid(nbins',nbins);
            model_fit = zeros(size(X));
            
            for i=1:size(X,1)
                for j=1:size(X,2)
                    model_fit(i,j) = non_lin_1D(find(hist_bins>(vec*[Y(i,j);X(i,j)]),1));
                end
            end
            %             model_fit = model_fit*max(non_lin(:))/max(model_fit(:));
            r_value = corrcoef(non_lin(:),model_fit(:));
            fprintf('Goodness of fit value for subunit-stimulus response is %d \n',r_value(1,2));
            
            figure(plot_counter);
            subplot(221),surfl(X,Y,non_lin); set(gca,'xlim',[xmin xmax],'ylim',[xmin xmax],'zlim',[min(non_lin(:)) max(non_lin(:))]);
            xlabel(Y_label), ylabel(X_label); title('Original firing rate surface');
            subplot(222),surfl(X,Y,model_fit); set(gca,'xlim',[xmin xmax],'ylim',[xmin xmax],'zlim',[min(model_fit(:)) max(model_fit(:))]);
            xlabel(Y_label), ylabel(X_label); title('Model-Predicted Firing Rate');
            subplot(223),plot(hist_bins,non_lin_1D,'LineWidth',2); set(gca,'xlim',[min(hist_bins) max(hist_bins)]);
            legend('1D Non-linearity'); title('Vector 1');
            subplot(224),surfl(X,Y,non_lin-model_fit);
            set(gca,'xlim',[xmin xmax],'ylim',[xmin xmax],'zlim',[min(non_lin(:)-model_fit(:)) max(non_lin(:)-model_fit(:))]);
            xlabel(Y_label), ylabel(X_label); title('Residual Firing Rate Error');
            plot_counter = plot_counter + 1;
            
        elseif(strcmp('6',str))
            % Fitting procedure same as 6, with a fit to the 1-D nonlinearity
            [spike_count,~] = hist3([projs(Lspike>0,1), projs(Lspike>0,2)],{nbins1,nbins1});
            [raw_count,~] = hist3([projs(:,1), projs(:,2)],{nbins1,nbins1});
            p1 = mean(projs(Lspike>0,1)) - mean(projs(:,1)); % based on the actual projection values
            p2 = mean(projs(Lspike>0,2)) - mean(projs(:,2));
            vec = [p1 p2]; vec = vec/sqrt(vec*vec'); % defining and normalizing the vector here
            spike_proj = [projs(Lspike>0,1) projs(Lspike>0,2)]*vec'; % dot product of the spike triggered data onto the vector
            raw_proj = projs*vec'; % dot product of the raw data onto the vector
            ub = ceil(max([max(spike_proj) max(raw_proj)]));
            lb = floor(min([min(spike_proj) min(raw_proj)]));
            hist_bins = lb:mean(diff(nbins)):ub;
            [spike_hist,~] = hist(spike_proj,hist_bins);
            [raw_hist,~] = hist(raw_proj,hist_bins); raw_hist(raw_hist==0) = 1;
            non_lin_1D = spike_hist./raw_hist;
            
            [X,Y] = meshgrid(nbins',nbins);
            model_fit = zeros(size(X));
            max_proj_val = vec(2)*X + vec(1)*Y;
            ind = find(hist_bins > max(max_proj_val(:))+0.1,1);
            non_lin_1D_cut = non_lin_1D(hist_bins>=-1 & hist_bins<= hist_bins(ind));
            hist_bins_cut = hist_bins(hist_bins>=-1 & hist_bins<=hist_bins(ind));
            
            prompt = 'Which fit would u prefer? (mle/lse)';
            str1 = input(prompt,'s');
            
            if(strcmp('mle',str1))
                % maximum likelihood estimate
                spike_hist_cut = spike_hist(hist_bins>=-1 & hist_bins<=hist_bins(ind));
                raw_hist_cut = raw_hist(hist_bins>=-1 & hist_bins<=hist_bins(ind));
                options = optimset('MaxFunEvals',1e10,'MaxIter',1e10);
                
                % exponential fit
                [w,fval] = fminsearch(@(params) mlefit(params,hist_bins_cut,raw_hist_cut,spike_hist_cut), [0.1,-2,0.05],options);
                vec_est = w(1)*exp(w(2)*(hist_bins_cut+ w(3)));
                val = [fval];
                
                %quadratic fit, 
                [w,fval] = fminsearch(@(params) mleqfit(params,hist_bins_cut,raw_hist_cut,spike_hist_cut), [1,1,1],options);
                vec_est = [vec_est; w(1)*hist_bins_cut.^2 + w(2)*hist_bins_cut+ w(3)];
                val = [val; fval];
                disp(val); % displaying the log likelihood values for the exponential and quadratic parametric fits 
                
                % Checking which one of the MLE is better by comparing the R values, 
                %find out the best fit and use that vec_est
                [~,I] = max(val);
                vec_est = vec_est(I,:);
                          
            elseif(strcmp('lse',str1))
                % least square estimate
                myfun1 = @(w)sum((non_lin_1D_cut-w(1)*exp(hist_bins_cut*w(2))- w(3)).^2);
                [w,~] = fminsearch(myfun1,[1,1,1]);
                vec_est = w(1)*exp(hist_bins_cut*w(2)) + w(3);
            end
            
            [r_value,model_fit] = get_R_value(non_lin,vec_est,vec,X,Y,hist_bins_cut);
            fprintf('Goodness of fit value for subunit-stimulus response is %d \n',r_value);
            
            figure(plot_counter);
            subplot(241); imagesc([min(nbins1) max(nbins1)],[min(nbins1) max(nbins1)], raw_count); xlabel(X_label); ylabel(Y_label); title('Raw Ensemble');
            subplot(245); imagesc([min(nbins1) max(nbins1)],[min(nbins1) max(nbins1)], spike_count); xlabel(X_label); ylabel(Y_label); title('Spike Triggered Ensemble');
            subplot(242),surfl(X,Y,non_lin); set(gca,'xlim',[xmin xmax],'ylim',[xmin xmax],'zlim',[min(non_lin(:)) max(non_lin(:))]);
            xlabel(Y_label), ylabel(X_label); title('Original firing rate surface');
            subplot(246), imagesc([xmin xmax],[xmin xmax],non_lin); xlabel(Y_label), ylabel(X_label); title('Original firing rate surface');
            subplot(243),surfl(X,Y,model_fit); set(gca,'xlim',[xmin xmax],'ylim',[xmin xmax],'zlim',[min(model_fit(:)) max(model_fit(:))]);
            xlabel(Y_label), ylabel(X_label); title('Model-Predicted Firing Rate');
            subplot(247), imagesc([xmin xmax],[xmin xmax],model_fit); xlabel(Y_label), ylabel(X_label); title('Model-Predicted Firing Rate');
            subplot(244),plot(hist_bins_cut,non_lin_1D_cut,'LineWidth',2); hold on; plot(hist_bins_cut,vec_est,'r','LineWidth',2); hold off;
            set(gca,'xlim',[min(hist_bins_cut) max(hist_bins_cut)]); legend('Data','Fit'); title('1-D Non-linearity'); xlabel('Projection values onto p_sta');
            subplot(248),imagesc([xmin xmax],[xmin xmax],non_lin-model_fit);
            xlabel(Y_label), ylabel(X_label); title('Residual Firing Rate Error');
            plot_counter = plot_counter + 1;
            disp(vec);
            
        elseif(strcmp('7',str))
            
            % Estimate the projection vector i.e. the vector that points to the direction of the cloud of the STA 
            % Fitting procedure same as option 6, with a fit to 1-D nonlinearity
            [spike_count,~] = hist3([projs(Lspike>0,1), projs(Lspike>0,2)],{nbins1,nbins1});
            [raw_count,~] = hist3([projs(:,1), projs(:,2)],{nbins1,nbins1});
            
            % This is the place where u introduce the function 'estimate_vector' calculates the MLE of the vector
            % which maximally distinguishes the raw ensemble from the spike triggered ensemble
            options = optimset('MaxFunEvals',1e7,'MaxIter',1e7);
            mle_mode = {'exp' 'quad'}; w = []; val = [];
            % optimization technique not quite working - something FISHY
            % here - ABHISHEK mend it
            iter = 20;
            for i = 1:numel(mle_mode)
                for j = 1:iter
                    [wt,fval] = fminsearch(@(params) estimate_vector(params,projs,Lspike,nbins,mle_mode{1,i}), [rand(1,3),randi(10),rand(1)],options);
                    w = [w; wt]; val = [val; fval];
                end      
            end
            [~,I] = max(val(1:50,:));
            w = w(I,:);
            
            vec = [w(1) w(2)]; vec = vec/sqrt(vec*vec'); % defining and normalizing the vector here
            spike_proj = [projs(Lspike>0,1) projs(Lspike>0,2)]*vec'; % dot product of the spike triggered data onto the vector
            raw_proj = projs*vec'; % dot product of the raw data onto the vector
            ub = ceil(max([max(spike_proj) max(raw_proj)]));
            lb = floor(min([min(spike_proj) min(raw_proj)]));
            hist_bins = lb:mean(diff(nbins)):ub;
            [spike_hist,~] = hist(spike_proj,hist_bins);
            [raw_hist,~] = hist(raw_proj,hist_bins); raw_hist(raw_hist==0) = 1;
            non_lin_1D = spike_hist./raw_hist;
            
            [X,Y] = meshgrid(nbins',nbins);
            model_fit = zeros(size(X));
            max_proj_val = vec(2)*X + vec(1)*Y;
            ind = find(hist_bins > max(max_proj_val(:))+0.1,1);
            non_lin_1D_cut = non_lin_1D(hist_bins>=-1 & hist_bins<= hist_bins(ind));
            hist_bins_cut = hist_bins(hist_bins>=-1 & hist_bins<=hist_bins(ind));
            
            if (I<=iter)
                vec_est = w(3)*exp(w(4)*(hist_bins_cut + w(5)));
            else
                vec_est = w(3)*hist_bins_cut.^2 + w(4)*hist_bins_cut + w(5);
            end
            
            [r_value,model_fit] = get_R_value(non_lin,vec_est,vec,X,Y,hist_bins_cut);
            fprintf('Goodness of fit value for subunit-stimulus response is %d \n',r_value);
            
            figure(plot_counter);
            subplot(241); imagesc([min(nbins1) max(nbins1)],[min(nbins1) max(nbins1)], raw_count); xlabel(X_label); ylabel(Y_label); title('Raw Ensemble');
            subplot(245); imagesc([min(nbins1) max(nbins1)],[min(nbins1) max(nbins1)], spike_count); xlabel(X_label); ylabel(Y_label); title('Spike Triggered Ensemble');
            subplot(242),surfl(X,Y,non_lin); set(gca,'xlim',[xmin xmax],'ylim',[xmin xmax],'zlim',[min(non_lin(:)) max(non_lin(:))]);
            xlabel(Y_label), ylabel(X_label); title('Original firing rate surface');
            subplot(246), imagesc([xmin xmax],[xmin xmax],non_lin); xlabel(Y_label), ylabel(X_label); title('Original firing rate surface');
            subplot(243),surfl(X,Y,model_fit); set(gca,'xlim',[xmin xmax],'ylim',[xmin xmax],'zlim',[min(model_fit(:)) max(model_fit(:))]);
            xlabel(Y_label), ylabel(X_label); title('Model-Predicted Firing Rate');
            subplot(247), imagesc([xmin xmax],[xmin xmax],model_fit); xlabel(Y_label), ylabel(X_label); title('Model-Predicted Firing Rate');
            subplot(244),plot(hist_bins_cut,non_lin_1D_cut,'LineWidth',2); hold on; plot(hist_bins_cut,vec_est,'r','LineWidth',2); hold off;
            set(gca,'xlim',[min(hist_bins_cut) max(hist_bins_cut)]); legend('Data','Fit'); title('1-D Non-linearity'); xlabel('Projection values onto p_sta');
            subplot(248),imagesc([xmin xmax],[xmin xmax],non_lin-model_fit);
            xlabel(Y_label), ylabel(X_label); title('Residual Firing Rate Error');
            plot_counter = plot_counter + 1;
            disp(val);
            disp(vec);
        end
        
    elseif(size(projs,2)==3)
        % fitting implemented when 3 subunits are present
        if (strcmp('1',str))
            bins = mid_st{1,2};
            x1 = repmat(bins,[numel(bins) 1 numel(bins)]);
            x2 = repmat(bins',[1 numel(bins) numel(bins)]);
            tmp = zeros(1,1,numel(bins));
            tmp(1,1,:) = bins';
            x3 = repmat(tmp,[numel(bins) numel(bins) 1]);
            myfun = @(w)sum(sum(sum((non_lin-w(1)*x1-w(2)*x2-w(3)*x3).^2)));
            [w_val,~] = fminsearch(myfun,[-1.2, 1, 1]);
            model_fit = w_val(1)*x1 + w_val(2)*x2 + w_val(3)*x3;
            r_value = corrcoef(non_lin(:),model_fit(:));
            fprintf('Goodness of fit value for subunit stimulus-response is %d \n',r_value(1,2));
            
        elseif(strcmp('2',str))
            % CODE BROKEN, HACK HERE - ABHISHEK
        elseif(strcmp('3',str))
        elseif(strcmp('4',str))
        elseif(strcmp('5',str))
        elseif(strcmp('6',str))
        elseif(strcmp('7',str))
        end
    end
end

end

