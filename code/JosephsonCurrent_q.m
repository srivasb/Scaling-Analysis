function [I0,I0_q0,crit_current,Diode_efficiency,tEnd] = JosephsonCurrent_q( filenum, height, num_cores, kpoint_chunk_size, varargin)

    if ~exist('filenum', 'var')
        filenum = 1;
    end

    if ~exist('height', 'var')
        height = 1626;
    end
    
    if ~exist('num_cores','var')
        num_cores = 4;
    end

    if ~exist('kpoint_chunk_size','var')
        % use for scaling analysis only
        kpoint_chunk_size = 500; % number of k-point to compute the cpr
    end

    %> Parser
    p = inputParser;
    % Global Inputs
    p.addParameter('resultsdir', 'results',@isstr); %
    p.addParameter('T', 0.04, @isscalar); % Temperature [K]
    p.addParameter('EF', 0.02, @isscalar); %The doping level in eV
    p.addParameter('pair_potential',0.0012); %eV
    % Spatial Potetial
    p.addParameter('surfacePotential', 0, @isscalar); %The doping potential inside the SNS junction in eV
    p.addParameter('lambda', 186.9612, @isscalar); %
	p.addParameter('y1', 0.28*height, @isscalar); %
    p.addParameter('vargamma_sc', []);
    p.addParameter('eta', 0); %
    % SOC
    p.addParameter('Zeeman_split', 2e-4, @isscalar); 
    p.addParameter('Rashba', 2.97e-3, @isscalar); %in eV
    p.addParameter('ValleyZeeman', 2.97e-3, @isscalar); % eV
    
    p.parse(varargin{:});
    
    % values
    T = p.Results.T;
    EF = p.Results.EF;
    vargamma_sc = p.Results.vargamma_sc;
    PairPotential = p.Results.pair_potential;
    surfacePotential = p.Results.surfacePotential;
    lambda = p.Results.lambda;
	y1 = p.Results.y1;
    eta = p.Results.eta;
    % SOC
    zeeman_split = p.Results.Zeeman_split;
    Rashba = p.Results.Rashba;
    ValleyZeeman = p.Results.ValleyZeeman;
    
    resultsdir = p.Results.resultsdir;

    % vec to store U=f(coord) potential
	% surfacePotential_vec = [];

	disp(['calculating Josephson Current with file extension ', num2str(filenum)]);

    filename = mfilename('fullpath');
    [directory, fncname] = fileparts( filename );
    
    %> filename containing the output XML
    outfilename = [fncname, '_',num2str( filenum )];
    %> Output XML file
    outputXML = [outfilename,'.xml'];
    %> the output folder
    outputdir   = [];
    
    %> creating output directories
    setOutputDir();
    
    %> filename containing the input XML    
    inputXML = 'Basic_Input_Graphene_SOC_Lattice.xml';
    %> Parsing the input file and creating data structures
    [Opt, param] = parseInput( fullfile( directory, inputXML ) );
    
    % number of cores
    Opt.workers = num_cores;

    % Imaginary Part of site-energy
    %{a
    if abs(eta) > 1e-6
		param.Leads{1}.epsilon = param.Leads{1}.epsilon + 1i*eta;
		param.Leads{2}.epsilon = param.Leads{2}.epsilon + 1i*eta;
    end
    %}

    %{
    if abs(eta) > 1e-6
		param.Leads{1}.epsilon = param.Leads{1}.epsilon - EF + 1i*eta;
		param.Leads{2}.epsilon = param.Leads{2}.epsilon - EF + 1i*eta;
    end
    %}

   % Hopping between scat. region and SC Leads
    if ~isempty( vargamma_sc )
		param.Leads{1}.vargamma_sc = vargamma_sc;
		param.Leads{2}.vargamma_sc = vargamma_sc;
    end
    
    % staggered A-B potential
    param.Leads{1}.deltaAB = 0*param.Leads{1}.vargamma;
    param.Leads{2}.deltaAB = 0*param.Leads{1}.vargamma;
    param.scatter.deltaAB = 0*param.Leads{1}.vargamma;
    
    % setting the Rashba-type SOC
    param.Leads{1}.Rashba = 0.*param.Leads{1}.vargamma;
    param.Leads{2}.Rashba = 0.*param.Leads{1}.vargamma;
    param.scatter.Rashba = Rashba;%0.001*param.Leads{1}.vargamma;
    
    % setting the Intrinsic SOC
    param.Leads{1}.Intrinsic = 0.00*param.Leads{1}.vargamma;
    param.Leads{2}.Intrinsic = 0.00*param.Leads{1}.vargamma;
    param.scatter.Intrinsic = 0.00*param.Leads{1}.vargamma;
    
    % setting the Valley Zeeman SOC
    param.Leads{1}.ValleyZeeman = 0.00*param.Leads{1}.vargamma;   
    param.Leads{2}.ValleyZeeman = 0.00*param.Leads{1}.vargamma;   
    param.scatter.ValleyZeeman = ValleyZeeman; 
   
    % Zeeman y potential strength
    %zeeman_split = 2e-4;
       
    %  setting the Fermi level and the height of the npn potential
    %EF = 1e-2;
    
    % superconducting pairing potential
    %Delta = 1e-3;
    Delta = abs(PairPotential);
    param.scatter.pair_potential = 0;
    param.Leads{1}.pair_potential = Delta;
    param.Leads{2}.pair_potential = Delta;
    
    
    % Planck contant
    h = 6.626e-34;
    hbar = h/(2*pi);
    % The charge of the electron
    qe = 1.602e-19;
    % atomic distance
    rCC = 1.42*1e-10; %In Angstrom
    % flux quanta
    flux0 = h/(2*qe);
    
    % the array of the phase differences
    DeltaPhi_vec = [];
    Edb = 1111;%2111;%bemenet(2);%885; 
    
    setvectors();
        
    total_width = 23500; %total width
    width_uc = 2; %width of the unit cell
	                    
%%   
    % for graphene infinite mass boundary
    % zigzag edged ribbon, q in units of 1/(3r_{cc})
    W = (total_width/2*2 + (total_width/2-1) ); %in units of r_CC
    m = 0:1:total_width/width_uc;
    interpq = transpose( (m+0.5)*pi/(W/3) ); % for graphene infinite mass boundary
    
    % Extract a representative chunk subset of k-points
    interpq_chunk = setvectors_kpoints_chunk(interpq,kpoint_chunk_size);
    
    Opt.Silent = false;
    Continue = 1;
    adaptive_moments = adaptiveQ( Opt, interpq_chunk, 'outFileName', [outputdir,'/',outfilename, '_adaptiveq.mat'], 'resume', Continue );
    Opt.Silent = true;
   
    
    % time loop parallel computation
    tStart = tic;
    
    parallelmanager = Parallel( Opt );
    parallelmanager.openPool();    
    
    % test to calculate the CPR fo a given q
    % currentVSDeltaPhi_q( 0 )
    
    currentvec_1range_vec = [];    
    
    adaptive_moments.runAdaptiveIterations( @currentVSDeltaPhi_q );
    parallelmanager.closePool();
    
    % time loop parallel end
    tEnd = toc(tStart);
    
    %> current result
    current_matrix = adaptive_moments.Read( 'interpy' );
    currentvec_q0 = current_matrix(:,1);
    
    currentvec_1range_vec = summation_q( current_matrix );
	phi_interp = min(DeltaPhi_vec):(max(DeltaPhi_vec)-min(DeltaPhi_vec))/500:max(DeltaPhi_vec);
	currentvec = interp1( DeltaPhi_vec, currentvec_1range_vec, phi_interp, 'spline');
    crit_current = max(abs(currentvec));

    % Capture I anomalous: I(DPhi=0)=I0
    I0 = currentvec_1range_vec(1); %*qe^2/hbar*1e9; % [nA]
    I0_q0 = currentvec_q0(1);

    % Diode factor
    IcUp = max(currentvec);
    IcDown = abs(min(currentvec));
    Diode_efficiency = (IcUp-IcDown)/(IcUp+IcDown);
    
    % save
    save([outputdir,'/',outfilename,'.mat'], 'currentvec_1range_vec', 'currentvec_q0','DeltaPhi_vec',...
        'crit_current','I0','I0_q0','IcUp','IcDown','Diode_efficiency',...
        'T','height','total_width','EF','zeeman_split','Rashba','ValleyZeeman');

    % plot
    % EgeszAbra();
	% save([outputdir,'/',outfilename,'.mat']);    
    
    
%% set vectors    
    function setvectors()
        phi_db = 20;%bemenet(3);%40;
        phi_max = 2*pi;%pi/5;
        deltaPhi = phi_max/(phi_db+1);     
        DeltaPhi_vec = 0:deltaPhi:phi_max;      
        %phi2vec = pi/3;              
    end

%% set k-pints representative, chunk size
    function ret = setvectors_kpoints_chunk(q_points,kpoint_chunk_size)
        total_num_kpoints = length(q_points);
        % generate index for getting random numbers
        indexes_chunk_kpoints = int16( sort( (total_num_kpoints-1)*rand([kpoint_chunk_size,1])+ 1 ) );
        % get the chunk of kpoints
        chunk_q_points = q_points(indexes_chunk_kpoints);
        % ensure q=0 is in chunk
        if q_points(1) == chunk_q_points(1)
            ret = chunk_q_points;
            return
        else
            aux_index_vec = zeros(size(indexes_chunk_kpoints));
            aux_index_vec(1) = 1;
            for ii=2:kpoint_chunk_size
                aux_index_vec(ii) = indexes_chunk_kpoints(ii-1);
            end
            
            ret = q_points(aux_index_vec);
        end
    end

    
%% currentVSDeltaPhi
    function current_1range = currentVSDeltaPhi_q( qvec )
        
        current_1range = zeros(length(qvec), length(DeltaPhi_vec)  );
    
        for idx = 1:length(qvec)       
%         for idx = 1:1 % test
            q = qvec(idx);
        			
   			hRibbon = Ribbon_tweak( 'Opt', Opt, 'param', param, 'width', width_uc, 'height', height, 'EF', EF, 'silent', true, 'filenameOut', [outputdir, filesep, outputXML], 'q', q );
            
            hZeemanPot = @ZeemanPot;
            SNSJosephson_handles = SNSJosephson( Opt, 'T', T, 'junction', hRibbon, 'gfininvfromHamiltonian', true, 'scatterPotential', hZeemanPot);

            %bandWidth = SNSJosephson_handles.getBandWidth();
            % Evec = createIntegrationContour(); % Original
            % 23 Aug
            %{a
            bandWidth = SNSJosephson_handles.getBandWidth();
            Emin = -max([bandWidth.Scatter.Emax, bandWidth.Lead.Emax ]);
            Emid = -bandWidth.Lead.Emin*10;
            Emax = 0;
            Evec = createIntegrationContour();
            % save('Evec.mat','Evec');
            % return
            %}
            
            current_tmp = SNSJosephson_handles.CurrentCalc_discrete( DeltaPhi_vec, 'Evec', Evec );
            current_1range(idx,:) = current_tmp;
            
            %> To save partial results
            save([outputdir,'/',outfilename,'_partial','.mat'],'current_1range',...
                'T','height','total_width','EF','zeeman_split','Rashba','ValleyZeeman');
            
        end
    
        %----------------------------------------
        % 23 Aug
        function Evec = createIntegrationContour( )
          
            %{a
            phi0 = 0.25*pi; %default
            
            radius1 = abs( Emax-Emid) /2;
            R1 = radius1/cos(phi0);
            
            radius2 = abs( Emid-Emin) /2;
            R2 = radius2/cos(phi0);
            
            Edb1 = fix(Edb*0.25);
            Edb2 = Edb-Edb1;
            
            phivec1 = phi0 +  0.5*(1.0 - cos(0:pi/(Edb1+1):pi) ) * (pi-2*phi0); %(0.5 + 0.5*sin( -0.5*pi:pi/(Edb1+1):0.5*pi )) * DeltaPhi + (pi-DeltaPhi)/2;
            phivec2 = phi0 +  0.5*(1.0 - cos(0:pi/(Edb2+1):pi) ) * (pi-2*phi0); %(0.5 + 0.5*sin( -0.5*ppi/(Edb2+1):0.5*pi )) * DeltaPhi + (pi-DeltaPhi)/2;
            
            Evec1 = R1*(cos(phivec1)-cos(phi0)) + 1i*R1*(sin(phivec1)-sin(phi0)) + Emax; %R1*(cos(phivec1) + 1i*sin(phivec1)) - radius1 + Emax - 1i*R1*sin((pi-DeltaPhi)/2);
            Evec2 = R2*(cos(phivec2)-cos(phi0)) + 1i*R2*(sin(phivec2)-sin(phi0)) + Emid; %R2*(cos(phivec2) + 1i*sin(phivec2)) - radius2 + Emid - 1i*R2*sin((pi-DeltaPhi)/2);
            
            Evec = [Evec1, Evec2];
            %}

        end

    end


%% ZeemanPot
% Zeeman potential applied in the scattering region
    function ret = ZeemanPot( CreateH, Energy )
        
        % getting coordinates
        coordinates = CreateH.Read('coordinates');
        
        % getting the scatter Hamiltonian
        Hscatter = CreateH.Read('Hscatter');
        
                
        % indices of spin up sites
        spinup = coordinates.spinup;
        
        indexes_spinup = 1:length(coordinates.x);
        indexes_spinup = indexes_spinup(spinup);

        indexes_spindown = 1:length(coordinates.x);
        indexes_spindown = indexes_spindown(~spinup);

        
        % allocate matrix for zeeman potential
        pot = sparse([], [], [], size(Hscatter,1), size(Hscatter, 2));
        
        % Bx spin down -> spin up 
        pot = pot + sparse( indexes_spinup, indexes_spindown, zeeman_split, length(coordinates.x), length(coordinates.x) );

        % Bx spin up -> spin down 
        pot = pot + sparse( indexes_spindown, indexes_spinup, zeeman_split, length(coordinates.x), length(coordinates.x) );
                
        % add the potential to the Hamiltonian
        Hscatter = Hscatter + pot;
        
        % save the modified Hamiltonain 
        CreateH.Write('Hscatter', Hscatter);
        
        %> use when No Spatial Potential
        ret = zeros(size(coordinates.x));
        
        % Spatial Potential
        %{
        if isempty(surfacePotential)
            ret = zeros(size(coordinates.x));
            return;
        end		
        
        %lead_doping = [param.Leads{1}.epsilon, param.Leads{2}.epsilon] - EF;
        lattice_const = norm(coordinates.a);
        
        %print(lattice_const);
        
        lambda_num = lambda*lattice_const; %transition length of the pn potential in units of rCC
        y1_num     = y1*lattice_const;
        y2_num     = height*lattice_const - y1_num;
        
        x = coordinates.x;
        y = coordinates.y;    

        if ~isempty(surfacePotential_vec) && ( norm(size( surfacePotential_vec) - size(y) ) == 0)
			ret = surfacePotential_vec;
			return;
        end
        
        % Original
        %ret = surfacePotential*( tanh((y-y1_num)/lambda_num) - tanh((y-y2_num)/lambda_num) )/2;
        % Nov 10
        %ret = (surfacePotential+EF)*( tanh((y-y1_num)/lambda_num) - tanh((y-y2_num)/lambda_num) )/2;
        
        % Jan 27
        lambda_L = lambda_num;
        lambda_R = lambda_L;
        ret = surfacePotential*( tanh((y-y1_num)/lambda_L) - tanh((y-y2_num)/lambda_R) )/2;
        
        surfacePotential_vec = ret;
        %}
            
    end



%--------------------------------------
    function current = summation_q( current_mtx)
        current = NaN( size(current_mtx,2),1);
        for iidx = 1:size(current_mtx,2)
            current(iidx) = sum( current_mtx(:,iidx));
        end
          
    end


%% plotfunction
    function EgeszAbra()
        
        figure1 = figure();
        
        fontsize = 16;
               
        Position = [0.18 0.18 0.53 0.53];
        axes1 = axes('Parent',figure1, 'Position', Position,...
                'Visible', 'on',...
                'FontSize', fontsize,...                 'xlim', x_lim,...                    'ylim', y_lim,...                'XTick', XTick,...                'YTick', YTick,...
                'Box', 'on',...
                'FontName','Times New Roman');
        hold on; 
        
        % Create xlabel
        %xlabel_position = [0 0 0];
        xlabel('\Delta\phi/\pi','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes1);
        %xlabel_handle = get(axes_cond,'XLabel');  
        %set(xlabel_handle, 'Position', get(xlabel_handle, 'Position') + xlabel_position);

        % Create ylabel
        %ylabel_position = [0 0 0];
        ylabel('I [nA]','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes1);
        %ylabel_handle = get(axes_cond,'YLabel'); 
        %set(ylabel_handle, 'Position', get(ylabel_handle, 'Position') + ylabel_position);
        
        legend_labels = [];
         
        if ~isempty( currentvec_1range_vec )
            indexes = ~isnan( currentvec_1range_vec );
            if ~isempty(DeltaPhi_vec(indexes))
                plot(DeltaPhi_vec(indexes)/pi, currentvec_1range_vec(indexes)*qe^2/hbar*1e9, 'c-', 'Parent', axes1)
            end
        end
        
           
        try
            print('-dpdf',[outputdir,'/',outfilename,'.pdf'],'-painters','-fillpage')
            disp('pdf figure was created successfully.')
        catch
            print('-depsc2', [outputdir,'/',outfilename,'.eps'])
            disp('eps figure was created successfully.')
            system(['convert ', outputdir,'/',outfilename,'.eps', ' ', outputdir,'/',outfilename,'.jpg'])
        end
	close(figure1)
        
        
    end



%% sets the output directory
    function setOutputDir()
        % resultsdir = 'results_Rashba';
        mkdir(resultsdir);
        outputdir = resultsdir;        
    end

        
        
    
end
