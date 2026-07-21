%%    Eotvos Quantum Transport Utilities - Ribbon
%    Copyright (C) 2009-2015 Peter Rakyta, Ph.D.
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see http://www.gnu.org/licenses/.
%
%> @addtogroup utilities Utilities
%> @{
%> @file Ribbon.m
%> @brief A class for calculations on a ribbon of finite width for equilibrium calculations mostly in the zero temperature limit.
%> @image html Ribbon_structure.jpg
%> @image latex Ribbon_structure.jpg
%> @} 
%> @brief A class for calculations on a ribbon of finite width for equilibrium calculations mostly in the zero temperature limit.
%> <tr class="heading"><td colspan="2"><h2 class="groupheader"><a name="avail"></a>Structure described by the class</h2></td></tr>
%> @image html Ribbon_structure.jpg
%> @image latex Ribbon_structure.jpg
%> The drawing represents a two-terminal structure made of two leads, of a scattering region and of two interface regions between the leads and scattering center.
%> Each rectangle describes a unit cell including singular and non-singular sites. 
%> The scattering center is also described by a set of identical unit cells, but arbitrary potential can be used.
%> Arrows indicate the hopping direction stored by the attributes in the corresponding classes (see attributes #CreateLeadHamiltonians.H1 and #InterfaceRegion.Hcoupling for details).
%> The orientation of the lead is +1 if the lead is terminated by the interface region in the positive direction, and -1 if the lead is terminated by the interface region in the negative direction.
%> (see attribute #CreateLeadHamiltonians.Lead_Orientation for details)
classdef Ribbon_tweak < NTerminal 

   
    properties ( Access = public )
        %> An instance of class #Lead (or its subclass) describing the unit cell of the scattering region
        Scatter_UC
        %> A function handle pot = f( #Coordinates ) or pot=f( #CreateLeadHamiltonians, Energy) of the transverse potential applied in the lead. (Instead of #CreateLeadHamiltonians can be used its derived class)
        transversepotential
        %> width of the scattering region (number of the nonsingular atomic sites in the cross section)        
        width
        %> height (length) of the scattering region (number of unit cells)        
        height
        %> the shift of the coordinates of the sites (two component vector)
        shift     
    end
    
    
methods ( Access = public )
%% Contructor of the class
%> @brief Constructor of the class.
%> @param varargin Cell array of optional parameters. For details see #InputParsing.
%> @return An instance of the class
    function obj = Ribbon_tweak( varargin )             
        obj = obj@NTerminal();     
        
        
        obj.param = [];        
        obj.Scatter_UC         = [];
        obj.transversepotential = [];      
        obj.width               = [];    
        obj.height              = [];
        obj.shift               = 0;      
        
        
    if strcmpi( class(obj), 'Ribbon_tweak') 
        % processing the optional parameters
        obj.InputParsing(varargin{:})                 

        
        obj.display(['EQuUs:Utils:', class(obj), ':Ribbon: Creating a Ribbon object'])
        
        % create the shape of the scattering region
        obj.createShape();
        
        % set the Fermi level
        obj.setFermiEnergy();
        
        % exporting calculation parameters into an XML format
        createOutput( obj.filenameOut, obj.Opt, obj.param );
        
        % create class intances and initializing class attributes
        obj.CreateHandles();
        
        % create Hamiltonians and coordinates of the unit cells 
        obj.CreateRibbon();
    end        
        

        
    end
    
%% Transport
%> @brief Calculates the transport through the two terminal setup on two dimensional lattices. Use for development pupose only.
%> @param Energy The energy value.
%> @param varargin Cell array of optional parameters (https://www.mathworks.com/help/matlab/ref/varargin.html):
%> @param 'constant_channels' Logical value. Set true (default) to keep constant the number of the open channels in the leads for each energy value, or false otherwise.
%> @param 'gfininvfromHamiltonian' Logical value. Set true calculate the surface Greens function of the scattering region from the Hamiltonaian of the scattering region, or false (default) to calculate it by the fast way (see Phys. Rev. B 90, 125428 (2014) for details).
%> @param 'decimateDyson' Logical value. Set true (default) to decimate the sites of the scattering region in the Dyson equation.
%> @param 'PotInScatter' Obsolete parameter. Use 'scatterPotential' instead.
%> @param 'scatterPotential' A function handle pot=f( #Coordinates ) or pot=f( #CreateHamiltonians, Energy) for the potential to be applied in the Hamiltonian (used when FiniteGreensFunctionFromHamiltonian=true).
%> @param 'selfEnergy' Logical value. Set true to use the self energies of the leads in the Dyson equation, or false (default) to use the surface Green function instead.
%> @param 'Smatrix' Set true (default) to calculate the conductance by using the scattering matrix via #Transport_Interface.Conduktance, or false to use function #Transport_Interface.Conductance2. (The latter one also works with complex energies.)
%> @return [1] Conductivity The calculated conductivity.
%> @return [2] aspect_ratio The aspect ratio W/L of the junction.
%> @return [3] Conductance The conductance tensor
%> @return [4] ny Array of the open channel in the leads.
%> @return [5] DeltaC Error of the unitarity.
%> @return [6] S The scattering matrix.
    function [Conductivity,aspect_ratio,Conductance,ny,DeltaC,S] = Transport(obj, Energy, varargin)
        
        p = inputParser;
        p.addParameter('constant_channels', false);
        p.addParameter('gfininvfromHamiltonian', false);
        p.addParameter('decimateDyson', true);
        p.addParameter('PotInScatter', []) %OBSOLETE use scatterPotential instead
        p.addParameter('scatterPotential', []) %NEW overrides optional argument 'PotInScatter'
        p.addParameter('selfEnergy', false)
		p.addParameter('Smatrix', true) %logcal value to use the S-matrix method to calculate the conductance or not.
        p.parse(varargin{:});
        constant_channels     = p.Results.constant_channels;
        gfininvfromHamiltonian = p.Results.gfininvfromHamiltonian;
        
        scatterPotential               = p.Results.PotInScatter;        
        if ~isempty( p.Results.scatterPotential )
            scatterPotential               = p.Results.scatterPotential;
        end
        
        decimateDyson          = p.Results.decimateDyson;
        selfEnergy             = p.Results.selfEnergy; 
		Smatrix                = p.Results.Smatrix;
        
        obj.CalcSpectralFunction(Energy, 'constant_channels', constant_channels, 'gfininvfromHamiltonian', gfininvfromHamiltonian, ...
            'decimateDyson', decimateDyson, 'scatterPotential', scatterPotential, 'SelfEnergy', selfEnergy);

        if strcmpi(obj.Opt.Lattice_Type, 'Graphene' ) && strcmp(obj.param.scatter.End_Type, 'A')
            aspect_ratio = (obj.width*3/2)/((obj.height+3)*sqrt(3));
        elseif strcmpi(obj.Opt.Lattice_Type, 'Graphene' ) && strcmp(obj.param.scatter.End_Type, 'Z')
            aspect_ratio = ((obj.width-0.5)*sqrt(3))/((obj.height+3)*3);
        elseif strcmpi(obj.Opt.Lattice_Type, 'Square' )
            aspect_ratio = obj.width/obj.height;
        else
            aspect_ratio = 1;
        end
        
        if Smatrix
        	[S,ny] = obj.FL_handles.SmatrixCalc();
        
        	norma = norm(S*S'-eye(sum(ny)));
            
        	if norma >= 1e-3
        	    obj.display( ['error of the unitarity of S-matrix: ',num2str(norma)] )
        	end
        
        	if ny(1) ~= ny(2)
        	   obj.display( ['openchannels do not match: ',num2str(ny(1)), ' ', num2str(ny(2))] )
        	end
        

        
        	conductance = obj.FL_handles.Conduktance();
        	conductance = abs([conductance(1,:), conductance(2,:)]);
        	DeltaC = std(conductance);
        
        	C = mean(conductance);
        	Conductivity = C/aspect_ratio;
        	Conductance = C;
        
        	obj.display( ['aspect ratio = ', num2str(aspect_ratio), ' conductance = ', num2str(C), ' open_channels= ', num2str(ny(1)), ' conductivity = ', num2str(Conductivity)])  

		else      
			Conductance = obj.FL_handles.Conductance2();
			DeltaC = NaN;
			Conductivity = Conductance/aspect_ratio;
			ny = NaN;
			S = [];
			obj.display( ['aspect ratio = ', num2str(aspect_ratio), ' conductance = ', num2str(Conductance), ' conductivity = ', num2str(Conductivity)])  
        end
        
    end

%% getCoordinates
%> @brief Gets the coordinates of the central region
%> @return [1] Coordinates of the central region.
%> @return [2] Coordinates of the interface region.
function [coordinates, coordinates_interface] = getCoordinates( obj )
    try
        if isempty( obj.CreateH ) || ~obj.CreateH.Read('HamiltoniansCreated')
            obj.CreateScatter()
        end
        
        coordinates = obj.CreateH.Read('coordinates');
        non_singular_sites = obj.CreateH.Read('kulso_szabfokok');
        indexes = false( size(coordinates.x) );
        indexes(non_singular_sites) = true;
        
        coordinates = coordinates.KeepSites( indexes );
        
    catch errCause
        err = MException('EQuUs:Utils:Ribbon:getCoordinatesFromRibbon', 'Error occured when retriving the coordinates of the ribbon.');
        err = addCause(err, errCause);
        save('Error_Ribbon_getCoordinatesFromRibbon.mat');
        throw( err );
    end
    
    
    try
        if isempty( obj.Interface_Regions )
            coordinates_interface = [];
            return
        end
        
        coordinates_interface = cell( length(obj.Interface_Regions), 1);
        for idx = 1:length( obj.Interface_Regions )
            coordinates_interface{idx} = obj.Interface_Regions{idx}.Read('coordinates');
        end
        
    catch errCause
        err = MException('EQuUs:Utils:Ribbon:getCoordinatesFromRibbon', 'Error occured when retriving the coordinates of the interface region from a Ribbon interface');
        err = addCause(err, errCause);
        save('Error_Ribbon_getCoordinatesFromRibbon.mat');
        throw( err );
    end    
        
end

%% ShiftCoordinates
%> @brief Shifts the coordinates of the sites in the ribbon by an integer multiple of the lattice vector. The coordinates of the Leads are automatically adjusted later
%> @param shift An integer value.
    function ShiftCoordinates( obj, shift )  
        obj.Scatter_UC.ShiftCoordinates( shift );
        if ~isempty( obj.Interface_Regions )
            for idx = 1:length(obj.Interface_Regions)
                obj.Interface_Regions{idx}.ShiftCoordinates( shift );
            end
        end
        
        obj.shift = obj.shift + shift;
    end

%% CreateScatter
%> @brief Initializes class #CreateHamiltonians for storing and manipulate the Hamiltonian of the the scattering region. The created object is stored in attribute #CreateH.
    function Scatter_UC = CreateScatter( obj ) 
        Scatter_UC = obj.Scatter_UC.CreateClone('empty', true);
        
        % create Hamiltonian of one unit cell of the scattering region
        Scatter_UC.CreateHamiltonians( 'toSave', 0);
        Scatter_UC.ShiftCoordinates( obj.shift ); 
        
        %applying transverse potential
        obj.ApplyTransversePotential( Scatter_UC )
       
        
       	% apply magnetic field in the unit cell of the scattering region
        % can be applied if the vector potential is identical in each unit cells
        if ~isempty( obj.PeierlsTransform_Scatter ) && obj.Opt.magnetic_field_trans_invariant %for finite q the vector potential must be parallel to q, and perpendicular to the unit cell vector
           	obj.display(['EQuUs:Utils:',class(obj),':CreateScatter: Applying magnetic field in the unit cell of the scattering region']);
           	obj.PeierlsTransform_Scatter.PeierlsTransformLeads( Scatter_UC );
        end
        
        % Create the Hamiltonian of the scattering region
        createH = CreateHamiltonians(obj.Opt, obj.param, 'q', obj.q);
        createH.CreateScatterH( 'Scatter_UC', Scatter_UC );    
        
       	% apply magnetic field in the whole Hamiltonian of the scattering region
        % can be applied for non-translational invariant vector potentials
        if ~isempty( obj.PeierlsTransform_Scatter ) && ~obj.Opt.magnetic_field_trans_invariant 
           	obj.display(['EQuUs:Utils:',class(obj),':CreateScatter: Applying magnetic field in the whole Hamiltonian of the scattering region']);
           	obj.PeierlsTransform_Scatter.PeierlsTransform( createH );
        end          
        
        
        obj.CreateH = createH;    
        
        
        Scatter_UC = obj.Scatter_UC.CreateClone('empty', true);
        
        params = Scatter_UC.Read( 'params' );
        params.ValleyZeeman = 0;
        Scatter_UC.Write( 'params', params );
        
        Scatter_UC.CreateHamiltonians( 'toSave', 0);
        
        % obtaining the sites of the scattering region that are coupled to the leads
        H0 = Scatter_UC.Read('H0');
        H1 = Scatter_UC.Read('H1');
        
        [rows, cols] = find( H1 );
        rows = unique(rows);   % cols identical to non_singular_sites
        cols = unique(cols);   % cols identical to non_singular_sites

        % identify sites that would be directly connected to the leads
        non_singular_sites = [reshape(cols, 1, length(cols)), (obj.height-1)*size(H0,1)+reshape(rows, 1, length(rows))];   
        createH.Write('kulso_szabfokok', non_singular_sites );

        
    end

%% setEnergy
%> @brief Sets the energy for the calculations
%> @param Energy The value of the energy in the same units as the Hamiltonian.
    function setEnergy( obj, Energy )
        
        setEnergy@NTerminal( obj, Energy );
        
        if ~isempty( obj.Scatter_UC ) && strcmpi(class(obj.Scatter_UC), 'Lead')  
            obj.Scatter_UC.Reset();
        end
        
        % recreate the Hamiltonian of the scattering region
        if ~isempty(obj.Scatter_UC)
            obj.CreateRibbon();
        end
        
    end    


%% CustomDysonFunc
%> @brief Custom Dyson function for a two terminal arrangement on a two dimensional lattice.
%> @param varargin Cell array of optional parameters (https://www.mathworks.com/help/matlab/ref/varargin.html):
%> @param 'gfininv' The inverse of the Greens function of the scattering region. For default the inverse of the attribute #G is used.
%> @param 'constant_channels' Logical value. Set true (default) to keep constant the number of the open channels in the leads for each energy value, or false otherwise.
%> @param 'onlyGinverz' Logical value. Set true to calculate only the inverse of the total Green operator, or false (default) to calculate #G as well.
%> @param 'recalculateSurface' A vector of the identification numbers of the lead surfaces to be recalculated.
%> @param 'decimate' Logical value. Set true (default) to eliminate all inner sites in the Greens function and keep only the selected sites. Set false to omit the decimation procedure.
%> @param 'kulso_szabfokok' Array of sites to be kept after the decimation procedure. (Use parameter 'keep_sites' instead)
%> @param 'selfEnergy' Logical value. Set true to use the self energies of the leads in the Dyson equation, or false (default) to use the surface Green function instead.
%> @param 'keep_sites' Name of sites to be kept in the resulted Green function (Possible values are: 'scatter', 'interface', 'lead').
%> @param 'UseHamiltonian' Set true if the interface region is matched to the whole Hamiltonian of the scattering center, or false (default) if the surface Green operator of the scattering center is used in the calculations.
%> @return [1] The calculated Greens function. 
%> @return [2] The inverse of the Green operator. 
%> @return [3] An instance of structure #junction_sites describing the sites in the calculated Green operator.
    function [Gret, Ginverz, junction_sites] = CustomDysonFunc( obj, varargin ) %NEW output
        
    p = inputParser;
    p.addParameter('gfininv', []);
    p.addParameter('constant_channels', true);
    p.addParameter('onlyGinverz', false );
    p.addParameter('recalculateSurface', [1 2] );
    p.addParameter('decimate', true );
    p.addParameter('kulso_szabfokok', []); %The list of sites to be left after the decimation procedure
    p.addParameter('SelfEnergy', false); %set true to calculate the Dyson equation with the self energy
    p.addParameter('keep_sites', 'lead'); %Name of sites to be kept (scatter, interface, lead)
    p.addParameter('UseHamiltonian', false); %true if the interface region is matched to the whole Hamiltonian of the scattering center, false if the surface Green operator of the scattering center is used in the calculations.
    p.parse(varargin{:});
    gfininv     = p.Results.gfininv;
    constant_channels = p.Results.constant_channels;
    onlyGinverz        = p.Results.onlyGinverz;
    recalculateSurface = p.Results.recalculateSurface;  
    decimate           = p.Results.decimate; 
    kulso_szabfokok = p.Results.kulso_szabfokok;
    useSelfEnergy         = p.Results.SelfEnergy;
    keep_sites        = p.Results.keep_sites;
    UseHamiltonian    = p.Results.UseHamiltonian;
    

    if ~isempty(recalculateSurface)
    
        % creating interfaces for the Leads
        if constant_channels
            shiftLeads = ones(length(obj.param.Leads),1)*obj.E;
        else
            shiftLeads = ones(length(obj.param.Leads),1)*0;
        end
    
        % creating Lead instaces and calculating the retarded surface Green operator/self-energy
        coordinates_shift = [-2, obj.height+1] + obj.shift;             
        obj.FL_handles.LeadCalc('coordinates_shift', coordinates_shift, 'shiftLeads', shiftLeads, 'transversepotential', obj.transversepotential, ...
            'SelfEnergy', useSelfEnergy, 'SurfaceGreensFunction', ~useSelfEnergy, 'gauge_field', obj.gauge_field, 'leads', recalculateSurface, 'q', obj.q, ...
            'leadmodel', obj.leadmodel);
        
        for idx = 1:length(recalculateSurface)
            obj.CreateInterface( recalculateSurface(idx), 'UseHamiltonian', UseHamiltonian );
        end
        
    end
    
    [Gret, Ginverz, junction_sites] = CustomDysonFunc@NTerminal(obj, 'gfininv', gfininv, ...
                                    'onlyGinverz', onlyGinverz, ...
                                    'recalculateSurface', [], ...
                                    'decimate', decimate, ...
                                    'kulso_szabfokok', kulso_szabfokok, ...
                                    'SelfEnergy', useSelfEnergy, ...
                                    'keep_sites', keep_sites);                                 

        
    end

%% CalcFiniteGreensFunction
%> @brief Calculates the Green operator of the scattering region by the fast way (see PRB 90, 125428 (2014)).
%> @param varargin Cell array of optional parameters identical to #NTerminal.CalcFiniteGreensFunction.
    function CalcFiniteGreensFunction( obj, varargin )
    
    	p = inputParser;
        p.addParameter('gauge_trans', false); % logical: true if want to perform gauge transformation on the Green's function and Hamiltonians   
        p.addParameter('onlyGinv', false);
        p.parse(varargin{:});
        gauge_trans     = p.Results.gauge_trans;
        onlyGinv        = p.Results.onlyGinv;   
        
        
        if obj.Opt.magnetic_field && ~obj.Opt.magnetic_field_trans_invariant
           warning( ['EQuUs:Utils:', class(obj), ':CalcFiniteGreensFunction'], 'The vector potential should be translational invariant in this calculation!' )
        end
        
        % recreate The Hamiltonian of the unit cell 
      	obj.CreateRibbon(); 
        
        if ~obj.Scatter_UC.Read('OverlapApplied')
            obj.Scatter_UC.ApplyOverlapMatrices(obj.E); 
        end
        
        % getting the Hamiltonian for the edge slabs
        [K0, K1, K1adj] = obj.Scatter_UC.qDependentHamiltonians();
       
        
        obj.display( ['EQuUs:Utils:', class(obj), ':CalcFiniteGreensFunction: Solving the eigenproblem in the scattering region'] )        
            
        % Trukkos sajatertek
        obj.display(['EQuUs:Utils:', class(obj), ':CalcFiniteGreensFunction: Eigenvalues of the scattering region'])
        obj.Scatter_UC.TrukkosSajatertekek(obj.E);
        % group velocity
        obj.display(['EQuUs:Utils:', class(obj), ':CalcFiniteGreensFunction: Group velocity for the scattering region'])
        obj.Scatter_UC.Group_Velocity();
        
        
        % transforming the Hamiltonians by SVD if necessary
        if obj.Scatter_UC.Read('is_SVD_transformed')
            V = obj.Scatter_UC.Get_V();
            K0 = V'*K0*V;
            K1 = V'*K1*V;
            K1adj = V'*K1adj*V;
        end
       

        %> the first two and last two slabs are added manualy at the end
        if ~obj.Scatter_UC.Read('is_SVD_transformed')
            z1 = 1;
            z2 = obj.height-2;
        else
            z1 = 2;
            z2 = obj.height-3;            
        end
        
        obj.G = [];
        obj.Ginv = [];
        
        obj.Scatter_UC.FiniteGreenFunction(z1,z2, 'onlygfininv', true);        
        obj.Ginv = obj.Scatter_UC.Read( 'gfininv' ); 
        
        % adding to the scattering region the first set of transition layers   
        Neff = obj.Scatter_UC.Get_Neff();
        non_singular_sites = obj.Scatter_UC.Read( 'kulso_szabfokok' ); 
        if isempty( non_singular_sites )
            non_singular_sites = 1:Neff;            
        end
        
        non_singular_sites_edges = [non_singular_sites, size(K0,1)+1:2*size(K0,1)];        
        ginv_edges = -[K0, K1; K1adj, K0];
        ginv_edges = obj.DecimationFunction( non_singular_sites_edges, ginv_edges );
        
        %Ginv = [first_slab, H1, 0;
        %             H1',  invG, H1;
        %              0  , H1', last_slab];        
        
        obj.Ginv = [ginv_edges(1:Neff,1:Neff), [ginv_edges(1:Neff, Neff+non_singular_sites), zeros(Neff, Neff+size(K0,2))]; ...
            [ginv_edges(Neff+non_singular_sites, 1:Neff); zeros(Neff, Neff)], obj.Ginv, [zeros(Neff, size(ginv_edges,2)-Neff); ginv_edges(1:Neff, Neff+1:end)]; ...
                            [zeros(size(K0,2),2*Neff), ginv_edges(Neff+1:end, 1:Neff)], ginv_edges(Neff+1:end, Neff+1:end)];  
                        
        [rows, ~] = find( K1 );
        rows = unique(rows);   % cols identical to non_singular_sites

        %non_singular_sites_Ginv = [1:length(non_singular_sites), size(obj.Ginv,1)-size(K0,1)+1:size(obj.Ginv,1)];
        non_singular_sites_Ginv = [1:length(non_singular_sites), size(obj.Ginv,1)-size(K0,1)+reshape(rows, 1, length(rows))];        
        obj.Ginv = obj.DecimationFunction( non_singular_sites_Ginv, obj.Ginv );
        
             
        % terminate the scattering region by the first and last slabs and transform back to the normal space from the SVD representation    
        if obj.Scatter_UC.Read('is_SVD_transformed')                           
            %> adding the first and last slab to 
            %Ginv = [first_slab, H1, 0;
            %             H1',  invG, H1;
            %              0  , H1', last_slab];        
        
            obj.Ginv = [K0, [K1(:, 1:Neff), zeros(size(K0,1), 2*size(K0,2))]; ...
                [K1adj(1:Neff,:); zeros(size(K0))], obj.Ginv, [zeros(Neff, size(K0,2)); K1]; ...
                            zeros(size(K0)), [zeros(size(K0,1), Neff), K1adj], K0];                        
        

            non_singular_sites_Ginv = [1:size(K0,1), size(obj.Ginv,1)-size(K0,1)+1:size(obj.Ginv,1)];
            obj.Ginv = obj.DecimationFunction( non_singular_sites_Ginv, obj.Ginv ); 
            
            % transform back to the normal space
            V_tot = [V, zeros(size(K0)); zeros(size(K0)), V];
            obj.Ginv = V_tot*obj.Ginv*V_tot';
        end
        
                     
%disp( obj.Ginv )           
        % gauge transformation of the vector potential in the effective Hamiltonians
        if gauge_trans
            try
            	% gauge transformation on Green's function
                if ~isempty(obj.Ginv) && ~isempty(obj.PeierlsTransform_Scatter) && isempty(obj.q) && ~isempty(obj.gauge_field)
                    coordinates_scatter = obj.getCoordinates();
                    %> gauge transformation on the inverse Green's function
                    obj.Ginv = obj.PeierlsTransform_Scatter.gaugeTransformation( obj.Ginv, coordinates_scatter, obj.gauge_field );
                end
            catch  errCause  
                err = MException(['EQuUs:Utils:', class(obj), ':CalcFiniteGreensFunction'], 'Unable to perform gauge transformation');
                err = addCause(err, errCause);
                save('Error_Ribbon_CalcFiniteGreensFunction.mat')
                throw(err);
            end
        end 
        
        if ~onlyGinv
            rcond_Ginv = rcond(obj.Ginv);
            if isnan(rcond_Ginv ) || abs( rcond_Ginv ) < 1e-15
                obj.display( ['EQuUs:Utils:', class(obj), ':CalcFiniteGreensFunction: Regularizing Ginv by SVD'], 1);
                obj.G = obj.inv_SVD( obj.Ginv );
            else
                obj.G = inv(obj.Ginv);
            end       
            obj.Ginv = [];
        else
            obj.G = [];
        end
        
    end

%% CalcFiniteGreensFunctionFromHamiltonian
%> @brief Calculates the Green operator of the scattering region from the whole Hamiltonian.
%> @param varargin Cell array of optional parameters (https://www.mathworks.com/help/matlab/ref/varargin.html):
%> @param 'gauge_trans' Logical value. Set true to perform gauge transformation on the Green operator and on the Hamiltonians.
%> @param 'onlyGinv' Logical value. Set true to calculate only the inverse of the surface Greens function #Ginv, or false (default) to calculate #G as well. In the latter case the attribute #Ginv is set to empty at the end.
%> @param 'PotInScatter' Obsolete parameter. Use 'scatterPotential' instead.
%> @param 'scatterPotential' A function handle pot=f( #Coordinates ) or pot=f( #CreateHamiltonians, Energy) for the potential to be applied in the Hamiltonian (used when FiniteGreensFunctionFromHamiltonian=true).
    function CalcFiniteGreensFunctionFromHamiltonian( obj, varargin )
        
    	p = inputParser;
        p.addParameter('gauge_trans', false); % logical: true if want to perform gauge transformation on the Green's function and Hamiltonians 
        p.addParameter('onlyGinv', false);
        p.addParameter('PotInScatter', []) %OBSOLETE use scatterPotential instead
        p.addParameter('scatterPotential', []) %NEW overrides optional argument 'PotInScatter', might be a cell array of function handles
        p.parse(varargin{:});
        gauge_trans              = p.Results.gauge_trans;
        onlyGinv                 = p.Results.onlyGinv;
        
        scatterPotential               = p.Results.PotInScatter;        
        if ~isempty( p.Results.scatterPotential )
            scatterPotential               = p.Results.scatterPotential;
        end
        
        % creating the Hamiltonian of the scattering region
        obj.CreateScatter();
        CreateH = obj.CreateH.CreateClone();

        % obtaining the Hamiltonian of the scattering region
        Hscatter = CreateH.Read('Hscatter');
        Hscatter_transverse = CreateH.Read('Hscatter_transverse');        
        
        % apply custom potential in the scattering center
        if ~isempty(scatterPotential)        
            if iscell( scatterPotential )
                for idx = 1:length( scatterPotential )
                    obj.ApplyPotentialInScatter( CreateH, scatterPotential{idx} );
                end
            else
                obj.ApplyPotentialInScatter( CreateH, scatterPotential); 
            end
            Hscatter = CreateH.Read('Hscatter');
        end          
        
        non_singular_sites_Ginv = CreateH.Read('kulso_szabfokok');
        
        % apply the periodic boundary condition in the transverse direction
        q = CreateH.Read('q');
        if ~isempty( q ) && ~CreateH.Read('HamiltoniansDecimated')
            Hscatter = Hscatter + Hscatter_transverse*diag(exp(1i*q)) + Hscatter_transverse'*diag(exp(-1i*q)); 
        end
         
        
        % reordering the sites of the scattering region to group the external point into the top corner of the matrix
        obj.display(['EQuUs:Utils:', class(obj), ':CalcFiniteGreensFunctionFromHamiltonian: Calculating the surface Green function of the scattering region.'])
        Hscatter = (sparse(1:size(Hscatter,1), 1:size(Hscatter,1), obj.E, size(Hscatter,1),size(Hscatter,1))   - Hscatter);  
        indexes = false( size(Hscatter,1), 1);
        indexes(non_singular_sites_Ginv) = true;
        Hscatter = [ Hscatter( ~indexes, ~indexes) , Hscatter(~indexes, indexes);
                     Hscatter(indexes, ~indexes), Hscatter(indexes, indexes)];

        obj.G = obj.partialInv( Hscatter, length(non_singular_sites_Ginv) );

        % recreate the Hamiltonian of the unit cell
        obj.CreateRibbon();
        
%disp( inv(obj.G) )         
        % gauge transformation of the vector potential in the effective Hamiltonians
        if gauge_trans
            try                
            	% gauge transformation on Green's function
                if ~isempty(obj.G) && ~isempty(obj.PeierlsTransform_Scatter) && isempty(obj.q) && ~isempty(obj.gauge_field)
                    surface_coordinates = obj.getCoordinates();
                    %> gauge transformation on the inverse Green's function
                    obj.G = obj.PeierlsTransform_Scatter.gaugeTransformation( obj.G, surface_coordinates, obj.gauge_field );
                end
            catch  errCause  
                err = MException(['EQuUs:Utils:', class(obj), ':CalcFiniteGreensFunctionFromHamiltonian', 'Unable to perform gauge transformation']);
                err = addCause(err, errCause);
                save('Error_Ribbon_CalcFiniteGreensFunctionFromHamiltonian.mat')
                throw(err);
            end
        end
              
        if onlyGinv
            rcond_G = rcond(obj.G);
            if isnan(rcond_G ) || abs( rcond_G ) < 1e-15
                obj.display( ['EQuUs:Utils:', class(obj), ':CalcFiniteGreensFunctionFromHamiltonian: Regularizing Ginv by SVD'],1);
                obj.Ginv = obj.inv_SVD( obj.G );
            else
                obj.Ginv = inv(obj.G);
            end       
            obj.G = [];
        else
            obj.Ginv = [];
        end
    end

%% CreateRibbon
%> @brief Creates the Hamiltonians of the unit cell of the ribbon shaped scattering region.
    function CreateRibbon( obj ) 
                
        % create the Hamiltonians for a unit cell of the scattering region
        if ~obj.Scatter_UC.Read( 'HamiltoniansCreated' )
        	obj.display( ['EQuUs:Utils:', class(obj), ':CreateRibbon: Creating ribbon Hamiltonian.'] )
        	obj.Scatter_UC.CreateHamiltonians( 'toSave', 0);
        	obj.Scatter_UC.ShiftCoordinates( obj.shift ); 
            
            %applying transverse potential
            obj.ApplyTransversePotential( obj.Scatter_UC )
        end     
        
       
        
       	% apply magnetic field in the unit cell of the scattering region
        if ~isempty( obj.PeierlsTransform_Scatter ) && ~obj.Scatter_UC.Read('MagneticFieldApplied') && obj.Opt.magnetic_field_trans_invariant
           	obj.display(['EQuUs:Utils:', class(obj), ':CreateRibbon: Applying magnetic field in ribbon Hamiltonians'])
           	obj.PeierlsTransform_Scatter.PeierlsTransformLeads( obj.Scatter_UC ); 
        end    
               
        
	
    end
    

%% CreateInterface
%> @brief Creates the Hamiltonians for the interface regions between the leads and scattering center.
%> @param idx Identification number of the interface region. 
%> @param varargin Cell array of optional parameters (https://www.mathworks.com/help/matlab/ref/varargin.html):
%> @param 'UseHamiltonian' Logical value. Set true if the interface region should be created to match to the whole Hamiltonian of the scattering center, false (default) if only the surface Green operator of the scattering center is used in the calculations.
	function CreateInterface( obj, idx, varargin )
        
        p = inputParser;
        p.addParameter('UseHamiltonian', false); %true if the interface region is matched to the whole Hamiltonian of the scattering center, false if the surface Green operator of the scattering center is used in the calculations.
        p.parse(varargin{:});
        UseHamiltonian     = p.Results.UseHamiltonian;
        
        %> Hamiltoninans of the interface region		
        Interface_Region = obj.Interface_Regions{idx};
        
        % The regularization of the interface is performed according to the Leads
        Leads = obj.FL_handles.Read( 'Leads' );
        Lead = Leads{idx};
        
        Interface_Region.Write( 'coordinates2', obj.getCoordinates() );
        Interface_Region.Write( 'K0', Lead.Read('K0'));
        Interface_Region.Write( 'K1', Lead.Read('K1'));
        Interface_Region.Write( 'K1adj', Lead.Read('K1adj'));
        Interface_Region.Write( 'K1_transverse', Lead.Read('K1_transverse'));
        Interface_Region.Write( 'K1_skew_left', Lead.Read('K1_skew_left'));
        Interface_Region.Write( 'K1_skew_right', Lead.Read('K1_skew_right'));
        Interface_Region.Write( 'coordinates', Lead.Read('coordinates'));
        Interface_Region.Write( 'kulso_szabfokok', Lead.Read('kulso_szabfokok'));
        Interface_Region.Write( 'OverlapApplied', true);
        
        coordinates_shift = [1, -1 ]; %relative to the leads
        Interface_Region.ShiftCoordinates( coordinates_shift(idx) );  
        
        %> coupling between the interface and the scattering region
        Surface_sc = obj.createSurface_sc( idx ); 
        Surface_sc.ApplyOverlapMatrices(obj.E);
        

        Lead_Orientation = Surface_sc.Read('Lead_Orientation');
        [~, K1, K1adj] = Surface_sc.qDependentHamiltonians();
        
        %K1 = Surface_sc.Read('K1');
        %K1adj = Surface_sc.Read('K1adj');
        [rows, cols] = find( K1 );
        rows = unique(rows);       
        cols = unique(cols);  % identical as non_singular_sites
        if Lead_Orientation == 1
            
            if UseHamiltonian
                Hscatter = obj.CreateH.Read('Hscatter');
                if isempty( Hscatter)
                    error(['EQuUs:utils:', class(obj), ':CreateInterface: Hamiltonian of the scattering region needs to be constructed first.'])
                end
                Neff = size(Hscatter,1)-size(K1,2);
                if ~Lead.Read('is_SVD_transformed')
                    % simple decimation                
                    Kcoupling = [K1, sparse([], [], [], size(K1,1), Neff )];
                    Kcouplingadj = [K1adj; sparse([], [], [], Neff, size(K1,1))];                 
                else
                    % regularization with SVD
                    Kcoupling = [K1, sparse([], [], [], size(K1,1), Neff )];
                    Kcouplingadj = [K1adj; sparse([], [], [], Neff, size(K1,1))]; 
                end                
            else
                Neff = length(rows); %number of coupling sites at Lead_Oriantation=-1
                if ~Lead.Read('is_SVD_transformed')
                    % simple decimation                
                    Kcoupling = [K1(:,cols), zeros( size(K1,1), Neff )];
                    Kcouplingadj = [K1adj(cols,:); zeros(Neff, size(K1,1))];                 
                else
                    % regularization with SVD
                    Kcoupling = [K1(:,cols), zeros(size(K1,1), Neff )];
                    Kcouplingadj = [K1adj(cols,:); zeros(Neff, size(K1,1))]; 
                end                                
            end
            

            
        elseif Lead_Orientation == -1
            
            if UseHamiltonian
                Hscatter = obj.CreateH.Read('Hscatter');
                if isempty( Hscatter)
                    error('EQuUs:utils:Ribbon:CreateInterface: Hamiltonian of the scattering region needs to be constructed first.')
                end
                Neff = size(Hscatter,1)-size(K1adj,2);
                if ~Lead.Read('is_SVD_transformed')
                    % simple decimation                
                    Kcoupling = [sparse([], [], [], length(cols), Neff ), K1adj(cols,:)];
                    Kcouplingadj = [sparse([], [], [], Neff, length(cols)); K1(:,cols)];                    
                else 
                    % regularization with SVD
                    Kcoupling = [sparse([], [], [], size(K1adj,1), Neff ), K1adj];
                    Kcouplingadj = [sparse([], [], [], Neff, size(K1,1) ); K1];                 
                end                   
            else
                Neff = length(cols); %number of coupling sites at Lead_Oriantation=1
                if ~Lead.Read('is_SVD_transformed')
                    % simple decimation                
                    Kcoupling = [zeros(length(cols), Neff), K1adj(cols,rows)];
                    Kcouplingadj = [zeros(Neff, length(cols)); K1(rows,cols)];                    
                else 
                    % regularization with SVD
                    Kcoupling = [zeros(size(K1,1), Neff), K1adj(:,rows)];
                    Kcouplingadj = [zeros(Neff, size(K1,2)); K1(rows,:)];                 
                end                                
            end
            

        else
            error('EQuUs:Utils:Ribbon:CreateInterface', 'Unknown lead orientation');            
        end
        
        
               
        Interface_Region.Write('Kcoupling', Kcoupling);
        Interface_Region.Write('Kcouplingadj', Kcouplingadj);  
        
        % method to adjust the interface region and coupling to the scattering region by an external function.
        if ~isempty( obj.interfacemodel )
            obj.interfacemodel( Interface_Region ); 
        end
        
        Interface_Region.Calc_Effective_Hamiltonians( 0, 'Lead', Lead );
        
    end

%% CreateClone
%> @brief Creates a clone of the present object.
%> @return Returns with the cloned object.
    function ret = CreateClone( obj )
        
        ret = Ribbon( 'width', obj.width, ...
            'height', obj.height, ...
            'filenameIn', obj.filenameIn, ...
            'filenameOut', obj.filenameOut, ...
            'E', obj.E, ...
            'EF', 0, ...
            'phi', obj.phi, ...
            'silent', obj.silent, ...
            'transversepotential', obj.transversepotential, ...
			'Opt', obj.Opt, ...
			'param', obj.param, ...
            'q', obj.q, ...
            'leadmodel', obj.leadmodel, ...
            'interfacemodel', obj.interfacemodel);
        
        ret.EF = obj.EF;
        ret.CreateH = obj.CreateH.CreateClone();
        ret.FL_handles = obj.FL_handles.CreateClone();
        ret.Scatter_UC = obj.Scatter_UC.CreateClone();
        ret.Interface_Regions = cell(size(obj.Interface_Regions));
        for idx = 1:length(obj.Interface_Regions)
            ret.Interface_Regions{idx} = obj.Interface_Regions{idx}.CreateClone();
        end
        if ~isempty( obj.PeierlsTransform_Scatter )
            ret.PeierlsTransform_Scatter = obj.PeierlsTransform_Scatter.CreateClone();
        end
        
        if ~isempty( obj.PeierlsTransform_Leads )
            ret.PeierlsTransform_Leads   = obj.PeierlsTransform_Leads.CreateClone();    
        end
        ret.gauge_field              = obj.gauge_field;  % function handle for the scalar field to transform the vector potential from Landauy to Landaux
        
    end  


end % methods public

methods ( Access = protected )

%% setFermiEnergy
%> @brief Sets the Fermi energy on the atomic sites for the calculations (use the same units as the elements of the Hamiltonian).
    function setFermiEnergy( obj )
        if ~isempty(obj.EF)
            obj.param.scatter.epsilon = obj.param.scatter.epsilon - obj.EF;
            for idx = 1:length(obj.param.Leads)
                obj.param.Leads{idx}.epsilon = obj.param.Leads{idx}.epsilon - obj.EF;
            end
        end      
    end

%% createSurface_sc
%> @brief Creates the copuling Hamiltonians between the scattering and interface region
%> @param idx The identification number of the interface region. (Integer value.)
%> @return An instance of class #Lead describing the copuling between the scattering and interface region
    function Surface_sc = createSurface_sc( obj, idx )
        Surface_sc = obj.Scatter_UC.CreateClone('empty', true); 
        
        params = Surface_sc.Read( 'params' );
        params.ValleyZeeman = 0;
        Surface_sc.Write( 'params', params );

        if ~isempty( obj.param.Leads{idx}.vargamma_sc )
            params = Surface_sc.Read( 'params' );
            params.vargamma = obj.param.Leads{idx}.vargamma_sc;
            Surface_sc.Write( 'params', params );
        end
        
        Surface_sc.CreateHamiltonians( 'toSave', 0);            
        if idx == 1
            coordinates_shift = 0 + obj.shift ;
        elseif idx == 2
            coordinates_shift = obj.height-1 + obj.shift;
        end
        Surface_sc.ShiftCoordinates( coordinates_shift )
        Surface_sc.Write('Hanyadik_Lead', idx);
        Surface_sc.Write('Lead_Orientation', obj.Interface_Regions{idx}.Read('Lead_Orientation'));
            
        %> applying transverse potential
        obj.ApplyTransversePotential( Surface_sc )        
        
        %> applying magnetic field
        if ~isempty( obj.PeierlsTransform_Leads )
            %> In superconducting lead one must not include nonzero magnetic
            %> field.
            %> Hamiltonians in transverse computations must remain
            %> traslational invariant. 
            if ~Surface_sc.isSuperconducting()
                obj.display('EQuUs:Utils:Ribbon:createSurface_sc: Applying magnetic field in the Hamiltonians')
               	obj.PeierlsTransform_Leads.PeierlsTransformLeads( Surface_sc );  
            else
                obj.display('EQuUs:Utils:Ribbon:createSurface_sc: Applying gauge transformation in the Hamiltonians')
               	obj.PeierlsTransform_Leads.gaugeTransformationOnLead( Surface_sc, obj.gauge_field );
            end
        end       
        
        
    end

%% ApplyTransversePotential
%> @brief Apply the tranvesre potential in the Hamiltonians
%> @param Scatter_UC An instance of class #Lead containing the Hamiltonians.
    function ApplyTransversePotential( obj, Scatter_UC )
       if ~isempty(obj.transversepotential) && isempty(obj.q) %In transverse computations no transverse potential can be applied
			coordinates = Scatter_UC.Read('coordinates');
       		if nargin( obj.transversepotential ) == 1
            	potential2apply = obj.transversepotential( coordinates );
	        elseif nargin( obj.transversepotential ) == 2
    	        potential2apply = obj.transversepotential( Scatter_UC, obj.E );
    	    else
    	        error('EQuUs:Utils:Ribbon:ApplyTransversePotential', 'To many input arguments in function handle scatterpotential');
    	    end

           	if isprop(coordinates, 'BdG_u')
				if size( potential2apply, 1) == 1 || size( potential2apply, 2) == 1
	                potential2apply(~coordinates.BdG_u) = -potential2apply(~coordinates.BdG_u);
				else
					potential2apply(~coordinates.BdG_u, ~coordinates.BdG_u) = -conj(potential2apply(~coordinates.BdG_u, ~coordinates.BdG_u));
				end
           	end
           	Scatter_UC.AddPotential( potential2apply );
        end 
    end




%% CreateHandles
%> @brief Initializes the attributes of the class.
    function CreateHandles( obj )  
        
        CreateHandles@NTerminal( obj )
        
        obj.Scatter_UC = obj.FL_handles.SurfaceGreenFunctionCalculator([], 'createCore', 1, 'q', obj.q);
 
    end

%% calculate_lead_attach_points
%> @brief Determines the site indexes at which the leads are connected to the scattering center.
    function calculate_lead_attach_points( obj )
        for idx = 1:length(obj.param.Leads)
            obj.param.Leads{idx}.M = obj.param.scatter.shape.width;
        end
    end

%% createShape
%> @brief Creates the geometry data of the ribbon shaped scattering region.
    function createShape( obj )
        
        if ~isempty( obj.width ) && ~isempty( obj.height )
            obj.param.scatter.shape.width = obj.width;
            obj.param.scatter.shape.height = obj.height;
        end
        
        if ~isempty(obj.param.scatter.shape.width)
            obj.width = obj.param.scatter.shape.width;
        else
            err = MException(['EQuUs:utils:', class(obj), ':createShape'], 'Shape is not given correctly, width is missing');
            throw(err)
        end
        
        if ~isempty(obj.param.scatter.shape.height)
            obj.height = obj.param.scatter.shape.height;
        else
            err = MException(['EQuUs:utils:', class(obj), ':createShape'], 'Shape is not given correctly, height is missing');
            throw(err)
        end
        
        obj.calculate_lead_attach_points();
       
    end
    
end % protected methods
    
    
methods (Access=protected)        


%% InputParsing
%> @brief Parses the optional parameters for the class constructor.
%> @param varargin Cell array of optional parameters (https://www.mathworks.com/help/matlab/ref/varargin.html):
%> @param 'width' Integer. The number of the nonsingular atomic sites in the cross section of the ribbon.
%> @param 'height' Integer. The height of the ribbon in units of the lattice vector.
%> @param 'filenameIn' The input filename containing the computational parameters. (Use parameters 'Op' and 'param' instead)
%> @param 'filenameOut' The output filename to export the computational parameters.
%> @param 'WorkingDir' The absolute path to the working directoy.
%> @param 'E' The energy value used in the calculations (in the same units as the Hamiltonian).
%> @param 'EF' The Fermi energy in the same units as the Hamiltonian. Attribute #E is measured from this value. (Use for equilibrium calculations in the zero temperature limit. Overrides the one comming from the external source)
%> @param 'silent' Set true to suppress output messages.
%> @param 'transversepotential' A function handle pot = f( #Coordinates ) or pot=f( #CreateLeadHamiltonians, Energy) of the transverse potential applied in the lead. (Instead of #CreateLeadHamiltonians can be used its derived class)
%> @param 'leadmodel' A function handle #Lead=f( idx, E, varargin ) of the alternative lead model with equivalent inputs and return values as #Transport_Interface.SurfaceGreenFunctionCalculator and with E standing for the energy.
%> @param 'interfacemodel' A function handle f( #InterfaceRegion ) to manually adjus the interface regions. (Usefull when 'leadmodel' is also given. For example see @InterfaceModel)
%> @param 'Opt' An instance of the structure #Opt.
%> @param 'param' An instance of the structure #param.
%> @param 'q' The transverse momentum quantum number.
    function InputParsing( obj, varargin )
    
        p = inputParser;
        p.addParameter('width', obj.width);
        p.addParameter('height', obj.height);
        p.addParameter('filenameIn', obj.filenameIn, @ischar);
        p.addParameter('filenameOut', obj.filenameOut, @ischar);
        p.addParameter('WorkingDir', obj.WorkingDir, @ischar);
        p.addParameter('E', obj.E, @isscalar);
        p.addParameter('EF', obj.EF);
        p.addParameter('silent', obj.silent);   
        p.addParameter('transversepotential', obj.transversepotential);
        p.addParameter('leadmodel', obj.leadmodel); %individual physical model for the contacts
        p.addParameter('interfacemodel', obj.interfacemodel); %individual physical model for the interface regions
        p.addParameter('Opt', obj.Opt);
        p.addParameter('param', obj.param);
        p.addParameter('q', obj.q);
        
        p.parse(varargin{:});
        
        InputParsing@NTerminal( obj, 'filenameIn', p.Results.filenameIn, ...
                            'filenameOut', p.Results.filenameOut, ...
                            'WorkingDir', p.Results.WorkingDir, ...
                            'E', p.Results.E, ...
                            'EF', p.Results.EF, ...
                            'silent', p.Results.silent, ...
                            'leadmodel', p.Results.leadmodel, ...
                            'interfacemodel', p.Results.interfacemodel, ...
                            'Opt', p.Results.Opt, ...
                            'param', p.Results.param, ...
                            'q', p.Results.q);

       
        obj.width        = p.Results.width;
        obj.height       = p.Results.height;
        obj.transversepotential    = p.Results.transversepotential;      
        
    end
    
end % methdos protected



end
