function Job_JosephsonCurrent_q(job_num,num_cores,kpoints_chunk_size)

    if ~exist('job_num','var')
        job_num='1';
    end
    job_num = str2double(job_num);

    if ~exist('num_cores','var')
        num_cores = '4';
    end
    num_cores = str2double(num_cores);
    
    if ~exist('kpoints_chunk_size','var')
        kpoints_chunk_size = '500';
    end
    kpoints_chunk_size = str2double(kpoints_chunk_size);
    
    % results dir 
    result_dir_name = ['Results_scaling_analysis_',num2str(job_num)];

    %> Global Parameters: Nano Letter Side contact parameters
    height = 518; % 1016; % ~249.64 nm % 1628; % ~400.41nm % sites
    T = 0.0400; % K
    EF = 0.350; % eV
    pair_potential = 0.0012; %1e-3; %eV
    vargamma_sc = 2.97*0.67;
    
    % Magnetic Field (Zeeman Split Value)
    num_Zemman_split_runs = 4; 
    Zeeman_split_vec = zeros(num_Zemman_split_runs,1);
    FillZeemanVec()
    
    %> Job inputs
    eta = 0 ; % Img site-energy % pair_potential*0.17 ;
    idx_zeeman_split = 1;
    Zeeman_split_val = Zeeman_split_vec(idx_zeeman_split);
    Rashba_val = 0 ; %2e-3 ; % 2.97e-3; % eV
    ValleyZeeman_val = 2e-3; % eV
    
    rcc=0.142; %nm
    Lnm = @(x) (x-1)*sqrt(3)*rcc;
    height_nm = Lnm(height);
    
    fprintf('\n')
    fprintf('num_cores = %d\n', num_cores);
    fprintf('Chunk size = %d\n',kpoints_chunk_size);
    
    fprintf('\n')
    fprintf('Job: %d | height = %.2f nm | T = %.4f K | EF = %.3f eV | Delta = %.4f eV | vargamma = %.3f | eta = %.4f |Evz = %.4f eV | Rashba = %.4f eV | Ezs = %.3fe-5 eV\n', ...
    job_num, height_nm, T, EF, pair_potential, vargamma_sc, eta, ValleyZeeman_val, Rashba_val, Zeeman_split_val*1e5);

    fprintf('\n');
    fprintf('NOT USING SPATIAL POTENTIAL!');
    fprintf('\n');
   
    
%% Results

    % File numbering
    aux_file_num = job_num;
    file_num = aux_file_num;
    
    % Runs:
    % measure run-time
    tStart = tic;
    
    %> Iterate over Zeeman_split_vec elements
    [I0_val,I0q0_val,crit_current_val,Diode_efficiency_val,Time0] = JosephsonCurrent_q(file_num,height,num_cores,kpoints_chunk_size,...
        'resultsdir',result_dir_name,...
        'T',T,'EF',EF,...
        'pair_potential',pair_potential,...
        'vargamma_sc',vargamma_sc,...
        'eta',eta,...
        'zeeman_split',Zeeman_split_val,...
        'Rashba',Rashba_val,...
        'ValleyZeeman',ValleyZeeman_val);

    % time end
    tEnd = toc(tStart);
    fprintf('Total Runtime = %.2f h \n',toc(tStart)/(60*60));
    
    % save
    save(['./',result_dir_name,'/data_scaling_analysis_job',num2str(job_num),'.mat'],...
            'I0_val','I0q0_val','crit_current_val','Diode_efficiency_val',...
            'Rashba_val','Zeeman_split_val','ValleyZeeman_val',...
            'height','T','EF','eta','num_cores','tEnd','Time0');
        
    fprintf('------------------------\n')
    fprintf('CODE ENDED CORRECTLY!\n')
    fprintf('------------------------\n')
        

    % plot
    % plot_DiodeEffi_result();
    
    
%% Extra Functions

    function FillZeemanVec()
        
        % Bohr Magneton
        % mu_B = qe*hbar/(2*me);
        mu_B = 0.579e-4; %eV/T
        
        % Zeeman Split Energy
        hZeemanSplit_func = @(B) mu_B*B;
        
        % Magnetic Field vec in T
        B_vec = 0:1/(num_Zemman_split_runs-1):1;
        
        % Fill vec
        Zeeman_split_vec = hZeemanSplit_func(B_vec)';
        
    end

    function ret_doping = doping_density(EF,Usurf)
        % Doping density in [cm]^-1 units respect to Fermy energy EF [eV] 
        % and Surface Doping Usurf [eV]

        if ~exist('EF', 'var')
            EF = 0.35; %[eV]
        end

        % constants
        qe = 1.602e-19; %[C]
        hbar = 6.582e-16;  %[eVs]
        hbar_J = 1.0546e-34; %[Js]
        vF = 1e6; %[m/s]

        % Use Js units of hbar
        % density(idx) = sign(EF-surfacePotential).*(EF-surfacePotential).^2*qe^2/(hbar_J^2*vF^2)/pi/1e4/1e11; 
        ret_doping = sign(EF-Usurf).*(EF-Usurf).^2*qe^2/(hbar_J^2*vF^2)/pi/1e4/1e11; %[C]^2/[cm]^-2

    end

    
end