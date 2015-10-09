%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulations for Figures: 3-6 
% "Using a novel source-localized phase regressor technique for evaluation
% of the vascular contribution to semantic category area localization in BOLD fMRI"
% By: An T Vu and Jack Gallant
% Frontiers in Neuroscience 2015
%

% generate stimlus time course
trsPerBlock = 16;
stim = repmat([0 1 0 1 0 1 0 1 0 1 0 1 0 1],[trsPerBlock 1]);
stim = stim(:);
stimOn = logical(stim);
stimOff = ~stimOn;
n = length(stim);

% noise standard deviation for phase and magnitude time-courses
p_sig = .1;
m_sig = .1;


for figNumber = 3:6
    
    switch figNumber
        case 3
            m_scales = .5; % magnitude BOLD signal strength
            p_scales = .5; % phase BOLD signal strength
        case 4
            m_scales = .5;
            p_scales = 0;
        case 5
            m_scales = 0;
            p_scales = .5;
        case 6
            m_scales = 0:.01:1;
            p_scales = 0:.01:1;
    end
    
    t_est1 = zeros(length(p_scales),length(m_scales));
    t_est2 = t_est1;
    t_est3 = t_est1;
    
    m_cnt = 1;
    for m_scale = m_scales
        p_cnt = 1;
        for p_scale = p_scales
            
            m = (stim*m_scale+randn(n,1)*m_sig);
            p = (stim*p_scale+randn(n,1)*p_sig);
            
            m_var = var(m);
            p_var = var(p);
            
            m_bar = mean(m);
            p_bar = mean(p);
            p2 = p-p_bar;
            
            % Cost function: OLS
            b1_1 = (n*sum(p.*m)-sum(m)*sum(p))/(n*sum(p.^2)-(sum(p)^2));
            b0_1 = m_bar - b1_1*p_bar;
            
            % Cost function: ChiSq as implemented in this manuscript (Vu 2015)
            r = sign(corr2(m,p));
            if r == 0; r= -1; end % handel zero case
            
            b1_3a = r*1/2*(p_var/m_var)*sum(p.^2-p_bar^2);
            b1_3b = (sum((p-p_bar).^2) - r*(p_var/m_var)*p_bar*sum(m-m_bar));
            b1_3c = -(sum(p.*(m-m_bar)) + r*1/2*(p_var/m_var)*sum((m-m_bar).^2));
            b1_3 = (-b1_3b+sqrt(b1_3b.^2-4*(b1_3a*b1_3c)))/(2*b1_3a);
            b0_3 = m_bar - b1_3*p_bar;
            
            
            % Estimation of macrovascular timecourses
            m_est1 = b0_1 + b1_1*p;
            m_est3 = b0_3 + b1_3*p;
            
            
            % Resultant phase regressor (microvascular) timecourses
            z1 = m-m_est1;
            z3 = m-m_est3;
            
            % t-values for magnitude and phase timecourses
            t_m = mean(m(stimOn)-m(stimOff))/((std(m(stimOn))+std(m(stimOff)))/2);
            t_p = mean(p(stimOn)-p(stimOff))/((std(p(stimOn))+std(p(stimOff)))/2);
            
            % t-values for phase regressor timecourses
            t_est1(p_cnt,m_cnt) = mean(z1(stimOn)-z1(stimOff))/((std(z1(stimOn))+std(z1(stimOff)))/2);
            t_est3(p_cnt,m_cnt) = mean(z3(stimOn)-z3(stimOff))/((std(z3(stimOn))+std(z3(stimOff)))/2);
            
            
            p_cnt = p_cnt+1;
        end
        m_cnt = m_cnt+1;
    end
    
    
    if length(m_scales) == 1
        
        t = 1:length(m);
        gl = -.3;
        yl = [-.4 1];
        
        h=figure;
        subplot(2,2,1)
        plot((p)*10,'g'); hold on
        
        ylim(yl*10); xlim([0 length(m)]);
        legend('phase data')
        xlabel('time'); ylabel('phase (deg)');
        title(sprintf('phase data, CNR=%.2f',t_p))
        for k = 1:7; h1 = line([1+(2*k-1)*16 2*k*16],[gl gl]*10); set(h1,'color',[.5 .5 .5]); end
        
        
        subplot(2,2,2)
        plot(m,'r'); hold on
        plot(m_est1,'c');
        plot(m_est3,'b');
        
        ylim(yl); xlim([0 length(m)])
        legend('magnitude data','phase-ChiSq fit','phase-OLS fit')
        xlabel('time'); ylabel('magnitude')
        title(sprintf('magnitude data, CNR=%.2f',t_m))
        for k = 1:7; h1 = line([1+(2*k-1)*16 2*k*16],[gl gl]); set(h1,'color',[.5 .5 .5]); end
        
        
        subplot(2,2,3);
        plot(m-m_est3,'b'); hold on
        tmp = ones(size(stim));
        tmp(stimOn) = mean(z3(stimOn));
        tmp(stimOff) = mean(z3(stimOff));
        plot(tmp,'k');
        
        ylim(yl); xlim([0 length(m)])
        legend('PR (PR-ChiSq)','model fit');
        xlabel('time'); ylabel('magnitude');
        title(sprintf('PR-ChiSq, CNR=%.2f',t_est3(1)))
        for k = 1:7; h1 = line([1+(2*k-1)*16 2*k*16],[gl gl]); set(h1,'color',[.5 .5 .5]); end
        
        
        subplot(2,2,4);
        plot(m-m_est1,'c'); hold on
        tmp = ones(size(stim));
        tmp(stimOn) = mean(z1(stimOn));
        tmp(stimOff) = mean(z1(stimOff));
        plot(tmp,'k')
        
        ylim(yl); xlim([0 length(m)])
        legend('sPR (PR-OLS)','model fit')
        xlabel('time'); ylabel('magnitude')
        title(sprintf('PR-OLS, CNR=%.2f',t_est1(1)))
        for k = 1:7; h1 = line([1+(2*k-1)*16 2*k*16],[gl gl]); set(h1,'color',[.5 .5 .5]); end
    end
    
    
    % for figure 6 only
    ind = [1 21 41 61 81 101];
    clim = [-5 5];
    if length(m_scales) > 1
        h=figure; subplot(1,2,2)
        imagesc(mean(t_est1,3)); caxis(clim); colormap(rwb(256)); colorbar
        set(gca,'YTick',ind,'YTickLabel',p_scales(ind)./p_sig,'XTick',ind,'XTickLabel',m_scales(ind)./m_sig)
        ylabel('phase fSNR'); xlabel('magnitude fSNR'); title('PR-OLS'); axis square
        subplot(1,2,1);
        imagesc(mean(t_est3,3)); caxis(clim); colormap(rwb(256)); colorbar
        set(gca,'YTick',ind,'YTickLabel',p_scales(ind)./p_sig,'XTick',ind,'XTickLabel',m_scales(ind)./m_sig)
        ylabel('phase fSNR'); xlabel('magnitude fSNR'); title('PR-ChiSq'); axis square
    end
end
