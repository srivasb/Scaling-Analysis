function result_plot_scaling_analysis()
    %> plot diode efficiency vs. VZ-split (or magnetic field)
    
    % Open data and save important variables
    job_num_vec = (1:6)' ;
    num_cores_vec = zeros(length(job_num_vec),1);
    Time0_vec = zeros(length(job_num_vec),1);
    tEnd_vec = zeros(length(job_num_vec),1);
    
    % load data
    for ii=1:length(job_num_vec)
        
        job_id = num2str(job_num_vec(ii));
        result_dir_name = ['Results_scaling_analysis_',job_id];
        main_path = fullfile('.',result_dir_name,['data_scaling_analysis_job',num2str(job_num),'.mat']);
        
        % save data
        if exist(main_path,'file')
            num_cores_temp = load(main_path,'num_cores');
            Time0_temp = load(main_path,'Time0');
            tEnd_temp = load(main_path,'tEnd');
            num_cores_vec(ii) = num_cores_temp.num_cores;
            Time0_vec(ii) = Time0_temp.Time0;
            tEnd_vec(ii) = tEnd_temp.tEnd; 
        else
            warning('Data File not found: %s', main_path)
        end
        
    end
    
    % Important plotting quantities
    %> Speedup
    SpeedUp_vec = Time0_vec(1)./tEnd_vec;
    %> Parallel Efficiency
    ParallelEfficiency_vec = (SpeedUp_vec ./ num_cores_vec );
    
    % plot
    hp = plot_scaling_analysis(num_cores_vec,SpeedUp_vec,ParallelEfficiency_vec) ;
    
    try
        close(hp.Parent) ;
    catch ME
        warning('Plot cannot be closed!')
        fprintf('Failed close plot. Reason: %s\n', ME.message);
    end
    
end

function plot_scaling_analysis(nCPU,S,E)


% Inputs:
% nCPU : Number of CPUs
% S    : Parallel Speedup
% E    : Parallel Efficiency (0-1)

close all

hfig = figure;

hfig.Color   = 'w';
hfig.Visible = 'on';

hfig.Units = 'centimeters';
hfig.Position = [2 2 11.5 10.5];

%% Upper axes : Speedup


hax1 = axes('Parent',hfig);

set(hax1,...
    'Units','centimeters',...
    'Position',[1.4 6.0 9.5 3.5],...
    'Box','on',...
    'FontName','Times',...
    'FontSize',11,...
    'LineWidth',1,...
    'TickDir','out');

hold(hax1,'on')

% Measured speedup
hp(1) = plot(hax1,nCPU,S,...
    '-o',...
    'LineWidth',1.5,...
    'MarkerSize',6,...
    'DisplayName','Measured');

% Ideal speedup
hp(2) = plot(hax1,nCPU,nCPU,...
    '--k',...
    'LineWidth',1.2,...
    'DisplayName','Ideal');

grid(hax1,'on')

ylabel(hax1,'Parallel Speedup, S','FontSize',11)

xlim(hax1,[min(nCPU) max(nCPU)])
ylim(hax1,[0 1.1*max([S(:);nCPU(:)])])

legend(hax1,...
    'Location','northwest',...
    'Box','off')


%% Lower axes : Efficiency

hax2 = axes('Parent',hfig);

set(hax2,...
    'Units','centimeters',...
    'Position',[1.4 1.2 9.5 3.5],...
    'Box','on',...
    'FontName','Times',...
    'FontSize',11,...
    'LineWidth',1,...
    'TickDir','out');

hold(hax2,'on')

hp(3) = plot(hax2,nCPU,100*E,...
    '-s',...
    'LineWidth',1.5,...
    'MarkerSize',6);

grid(hax2,'on')

xlabel(hax2,'Number of CPUs','FontSize',11)

ylabel(hax2,'Parallel Efficiency (%)','FontSize',11)

xlim(hax2,[min(nCPU) max(nCPU)])
ylim(hax2,[0 105])


%% Export

set(hfig,'PaperUnits','centimeters');
set(hfig,'PaperSize',[11.5 10.5]);
set(hfig,'PaperPosition',[0 0 11.5 10.5]);

print(hfig,'ScalingAnalysis','-dpdf','-painters');

% Optional PNG (600 dpi)
print(hfig,'ScalingAnalysis','-dpng','-r600');


end


%{
function hp = plot_scaling_analysis(num_cores,SpeedUp,ParallelEfficiency)


    hfig = figure;
    hfig.Visible = 0;
    hax = gca;
    hax.Parent = hfig;
    hax.Box = 'on';
    hax.Title.FontWeight = 'normal';
    hax.FontName = 'Times';
    % hax.FontName = 'Helvetica';
    %{a
    figureRatio = 0.45;
    defaultFigureWidth = 11.5; % textwidth for margin = 2.5cm
    hfig.Units = 'centimeters';
    hfig.Position = [0 0 defaultFigureWidth figureRatio*defaultFigureWidth];
    %Axes Position
    AxesRatio = 0.42;
    defaultAxesWidth = 10;
    defaultAxisHeigth = defaultAxesWidth*AxesRatio;
    hax.Units = 'centimeters';
    hax.Position = [1.1 0.8 defaultAxesWidth defaultAxisHeigth]; 

    % Colors
    blue_Color_gradient_vec = ['#A7C7E7';'#6699CC';'#336699';'#1F4E79';'#0B2E59'];

    

    % sub plots:
    hp = plot(num_,ParallelEfficiency,'-',...
                            'LineWidth',1,...
                            'Color',blue_Color_gradient_vec(3,:),...
                            'Marker','o',...
                            'MarkerSize',8,...
                            'MarkerEdgeColor',blue_Color_gradient_vec(5,:),...
                            'MarkerFaceColor',blue_Color_gradient_vec(5,:));
                        
    % names
    title('Rashba=0, $VZ=2\,10^{-3}$','Interpreter','latex')
    xlabel('$\Delta E_{Zeeman}\,[eV]\,(10^{-4})$','Interpreter','latex')
    ylabel('$\eta$','Interpreter','latex')
   
    % Extra congifuration
    hax.Box = 'on';
    hax.Units = 'points';
    hax.TickLabelInterpreter = 'latex';
    hax.Title.FontSize = 10;
    hax.XLabel.FontSize = 11;	
    hax.YLabel.FontSize = 11;
    hax.XAxis.FontSize = 10;
    hax.YAxis.FontSize = 10;
    hax.XRuler.Axle.LineWidth = 0.6;
    hax.YRuler.Axle.LineWidth = 0.6;
    hax.XRuler.TickLabelGapOffset = -1; % positive -> down / negative -> up
    hax.YRuler.TickLabelGapOffset = 1; % positive -> left / negative -> right
    hax.TickDir = 'in';
    hax.YMinorTick = 'on';
    
    % save fig:
    pos = get(hfig,'Position');
    hfig.PaperPositionMode = 'Auto';
    hfig.PaperUnits = 'centimeters';
    hfig.PaperSize = [pos(3), pos(4)];
    
    output_filename = 'results_diode_VZ.pdf';
    try
        exportgraphics(hfig, output_filename, 'ContentType', 'vector');
        fprintf('Vector PDF plot successfully created: %s\n', output_filename);
    catch ME
        warning('Failed to export vector PDF. Saving PNG image copy instead.');
        fprintf('Failed to pdf plot. Reason: %s\n', ME.message);
        
        exportgraphics(hfig, 'results_diode_VZ.png', 'Resolution', 300);
    end

    % close fig
    close(hfig);

end

%}


