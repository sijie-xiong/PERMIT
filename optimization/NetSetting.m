classdef NetSetting
    
    properties
        synthetic
        device
        where
        input_typ
        data
        
        n_dev
        n_eps
        n_rho
        i_eps
        i_rho
        eps_all
        rho_all
        min_rho
        
        in_sizes
        in_rates
        arrival_rate
        num_slots_perpkt
        in_byte_rate
        coeff_var
        
        shaper
        scaled_alphabet
        increased_alphabet
        how_to_increase
        
        out_sizes
        same_out_sizes
        N
        hN
        
        pmf
        scale
        
        eps_pool = [0.01, 0.1, 1, 2, 5, 10];
    end
    
    methods
        %% Constructor
        function self = NetSetting(varargin)
            % Num of args must be even.
            if rem(nargin, 2) ~= 0
                error('NetSetting needs propertyName/propertyValue pairs!');
            end
            % Set default values from a struct.
            args = struct('where','cluster','synthetic',false,...
                          'device','SleepCamSwitchMerged','input_typ','iid',...
                          'scale',1,'pmf','Zipf',...
                          'shaper','DPS_VSCI',...
                          'n_dev',2,'n_eps',5,'n_rho',6,...
                          'same_out_sizes',true,...
                          'scaled_alphabet',false,...
                          'increased_alphabet',false,...
                          'how_to_increase','same_bin');
            args_names = fieldnames(args);
            for attr = reshape(args_names,1,[])
                name = attr{1};
                self.(name) = args.(name); % .(name) = getattr(name) 
            end
            % Change attribute values passed in by user.
            for pair = reshape(varargin, 2, []) % name-value pairs
                name = pair{1};
                if any(strcmp(name, args_names))
                    % set the property values.
                    self.(name) = pair{2};
                else
                    error('NetSetting does not have %s property!', name);
                end
            end
            
            if self.n_eps > length(self.eps_pool)
            	error('Not enough eps values in the pool!');
            end
            self.eps_all = self.eps_pool(1:self.n_eps);
            
            if ~self.synthetic
                self = self.set_data();
            end
            self.in_sizes = self.get_in_sizes();
            if self.scaled_alphabet
                self.in_sizes = 3*self.in_sizes;
            end
            self.in_rates = self.get_in_rates();
            self.in_byte_rate = self.get_mu();
            self.coeff_var = self.get_cv();
            self.arrival_rate = 1 - self.in_rates(1);
            self.num_slots_perpkt = 1 / self.arrival_rate;
            
            if self.same_out_sizes
                self.out_sizes = self.in_sizes;
                min_rho = self.in_byte_rate / max(self.out_sizes);
                if strcmp(self.shaper, 'DPS_VSVI') || strcmp(self.shaper, 'VSVI_CSVI')
                    min_rho = min_rho / self.arrival_rate;
                end
                self.rho_all = linspace(0.98, min_rho+0.01, self.n_rho);
                self.min_rho = min_rho;
            else
                min_rho = self.in_byte_rate / max(self.in_sizes) / self.arrival_rate;
                rho_all_tmp = linspace(0.98, min_rho+0.01, self.n_rho);
                
                out_sizes = ceil( self.in_byte_rate ./ rho_all_tmp / self.arrival_rate );
                self.rho_all = self.in_byte_rate ./ out_sizes / self.arrival_rate;
                self.out_sizes = [0 out_sizes];
            end
            self.N = length(self.in_sizes);
            self.hN = length(self.out_sizes);
            
            if self.increased_alphabet
                n_dev = 4;
                self = self.get_new_ns(n_dev);
            end
        end
        
        %% Get (eps, rho) values
        function [eps, rho] = get_eps_rho(self, i_eps, i_rho)
            if (i_eps > self.n_eps || i_eps < 1 || i_rho > self.n_rho || i_rho < 1)
                error("i_eps & i_rho must be in [1,%d] and [1,%d] respectively!",...
                      self.n_eps, self.n_rho);
            end
            eps = self.eps_all(i_eps);
            rho = self.rho_all(i_rho);
        end
        
        %% Optimization
        function res = optim_shaper(self, eps, rho)
            switch self.shaper
                case {'DPS_VSCI', 'DPS_VSVI'}
                    nx = self.N * self.hN;
                    
                    % Objective functions
                    fun_Qsz = @(x)self.WHF(x, 1e-4);
                    fun_Bout = reshape( diag(self.in_rates) * kron(ones(self.N,1), self.out_sizes), nx, 1);
                    
                    % Linear equality constraint - byte rate (in_byte_rate/rho = out_byte_rate)
                    Aeq_byterate = fun_Bout';
                    beq_byterate = self.in_byte_rate / rho;

                    % Linear equality constraint - rows sum to ones
                    Aeq_rowsum = kron(ones(1,self.hN), eye(self.N));
                    beq_rowsum = ones(self.N, 1);
                    
                    % Linear inequality constraints - privacy
                    % priv constr for each column
                    if strcmp(self.shaper, 'DPS_VSCI')
                        [eps_s, eps_t] = deal(eps, eps); 
                        tt1 = [diag(ones(self.N,1)*exp(eps_t/2));...
                               kron(ones(self.N-1,1), diag([exp(eps_t/2); ones(self.N-1,1)*exp(eps_s)]))];
                        tt2 = kron(eye(self.N), ones(self.N,1)) - tt1; 
                    else % 'DPS_VSVI'
                        tt1 = kron(ones(self.N-1,1), diag([0; ones(self.N-1,1)*exp(eps)]));
                        I = eye(self.N);
                        tt2 = kron(I(2:end,:), [0; ones(self.N-1,1)]) - tt1;
                    end
                    % repeat priv constr for all columns
                    A_priv = kron( eye(self.hN), tt2 ); 
                    A_priv(isnan(A_priv)) = 0;
                    [n_privcon, ~] = size(A_priv);
                    b_priv = zeros(n_privcon, 1);
                    %     AA = kron( eye(ns.hN), kron(eye(ns.N), ones(ns.N,1)) - ...
                    %                           exp(eps) * kron(ones(ns.N,1), eye(ns.N)) );
                case {'VSCI_CSCI', 'VSVI_CSVI'}
                    nx = self.hN;
                    
                    if strcmp(self.shaper, 'VSCI_CSCI')
                        fun_Qsz = @(x)self.WHF(ones(self.N,1)*x', 1e-4);
                        Aeq_byterate = self.out_sizes;
                    else
                        I = eye(self.hN);
                        top_row = I(1,:);
                        fun_Qsz = @(x)self.WHF([top_row; ones(self.N-1,1)*x'], 1e-4);
                        Aeq_byterate = self.arrival_rate * self.out_sizes;
                    end
                    beq_byterate = self.in_byte_rate / rho;
                    
                    % Linear equality constraint - rows sum to ones
                    Aeq_rowsum = ones(1,self.hN);
                    beq_rowsum = 1;
                    
                    A_priv = []; b_priv = [];
                otherwise
                    error("Choose 'shaper' from 'DPS_VSCI', 'DPS_VSVI', 'VSCI_CSCI' or 'VSVI_CSVI'!");
            end
            % Concatenate equality constraints
            Aeq = [Aeq_byterate; Aeq_rowsum];
            beq = [beq_byterate; beq_rowsum];
            
            % Lower and upper bounds on the entries of q or \mu or \nu
            lb = zeros(nx,1); ub = ones(nx,1);
            
            % Constrained optimization
            exitflag = -2;
            while exitflag ~= 1
                % Random initial guess
                switch self.shaper
                    case {'DPS_VSCI', 'DPS_VSVI'}
                        tt1 = triu(ones(self.N,self.hN)).*rand(self.N,self.hN);
                        tt1 = tt1./(sum(tt1,2)*ones(1,self.hN));
                        tt1(isnan(tt1)) = 0;
                        if self.N > self.hN
                            tt1(self.hN+1:self.N,self.hN) = 1;
                        end
                        x0 = tt1(:);
                    case {'VSCI_CSCI', 'VSVI_CSVI'}
                        tt1 = rand(self.hN-1,1);
                        tt1 = tt1/sum(tt1);
                        x0 = [0;tt1];
                end
             
                try
                    options = optimoptions('fmincon','Display','iter','Algorithm','sqp',...
                                           'MaxFunctionEvaluations',200*self.N*self.hN,...
                                           'ConstraintTolerance',1e-10,...
                                           'UseParallel',true);
                    [x,fval,exitflag] = fmincon(fun_Qsz,x0,A_priv,b_priv,Aeq,beq,lb,ub,[],options);
                    %     [xx,ffval,eexitflag] = fmincon(fun,x0,AA,b,Aeq,beq,lb,ub,[],options);
                catch
                    warning("There's an optimization iterate preventing WHF from converging!");
                    continue;
                end
            end
            
            switch self.shaper
                case {'DPS_VSCI', 'DPS_VSVI'}
                    Q_opt = reshape(x,[self.N,self.hN]);
                    [eps_s_actual, eps_t_actual] = self.eps_actual(x);
                    rho_est = self.rho_est(x);

                    disp('Q_opt = ');
                    disp(Q_opt);
                    fprintf('eps = %f, eps_s_actual = %f, eps_t_actual = %f\n', ...
                            eps, eps_s_actual, eps_t_actual);
                    fprintf('rho = %f, rho_est = %f\n', rho, rho_est);
                    fprintf('qsz_est = %f\n', fval);

                    res.eps_s_actual = eps_s_actual;
                    res.eps_t_actual = eps_t_actual;
                    res.rho_est = rho_est;
                    res.qsz_est = fval;
                    res.Q_opt = Q_opt;
                case {'VSCI_CSCI', 'VSVI_CSVI'}
                    mu_opt = x;
                    if strcmp(self.shaper, 'VSCI_CSCI')
                        rho_est = self.rho_est(ones(self.N,1)*x');
                    else
                        rho_est = self.rho_est([top_row; ones(self.N-1,1)*x']);
                    end

                    disp('mu_opt = ');
                    disp(mu_opt);
                    fprintf('rho = %f, rho_est = %f\n', rho, rho_est);
                    fprintf('qsz_est = %f\n', fval);

                    res.rho_est = rho_est;
                    res.qsz_est = fval;
                    res.mu_opt = mu_opt;
            end
        end
        
        %% WHF
        function queue_sz_est = WHF(self, x, delta)
            Q_opt = x;
            if any(size(Q_opt) ~= [self.N, self.hN])
                Q_opt = reshape( Q_opt, [self.N, self.hN] )
            end
            
            rho_est = self.rho_est(Q_opt);
            
            if (any(x(:)<-1e-6) || any(x(:)>1+1e-6) || ...
                any(sum(Q_opt, 2) > 1+1e-6) || ...
                any(sum(Q_opt, 2) < 1-1e-6) || ...
                rho_est >= 1 || rho_est < 0)
                queue_sz_est = Inf;
            else
                ladder_heights = (kron(self.in_sizes', ones(1,self.hN)) - ...
                                  kron(self.out_sizes, ones(self.N,1)));
                tt = kron(self.in_rates', ones(1,self.hN));
                ladder_heights_prob = tt .* Q_opt;
                
                g = abs(min(ladder_heights(:)));
                h = max(ladder_heights(:));
                
                u_dict = containers.Map('KeyType','double','ValueType','any');
                for j = (-g):h
                    u_dict(j) = 0;
                end
                
                for i = 1:self.N
                    for j = 1:self.hN
                        ladder = ladder_heights(i,j);
                        u_dict(ladder) = u_dict(ladder)+ladder_heights_prob(i,j);
                    end
                end
                
                a_dict = containers.Map('KeyType','double','ValueType','any');
                for i = 1:h
                    a_dict(i) = 0;
                end
                cnt = 0;
                
                b_dict = containers.Map('KeyType','double','ValueType','any');
                while true
                    cnt = cnt + 1;
                    if cnt > 5000
                        disp('Utilization: ')
                        disp(rho_est)
                        disp('Q_opt: ')
                        disp(Q_opt)
                        disp('Row sum: ')
                        disp(sum(Q_opt,2))
%                         pause;
                        error('This optimization iterate prevents WHF from converging!');
                    end
                
                    for j = g:-1:0
                        b_dict(j) = u_dict(-j);
                        for i = 1:h
                            if isKey(b_dict, j+i)
                                b_dict(j) = b_dict(j) + a_dict(i) * b_dict(j+i);
                            else
                                break;
                            end
                        end
                    end
                    b_dict(0) = min(1, b_dict(0));

                    S = sum(cell2mat(values(b_dict))) - b_dict(0);

                    max_diff = -inf;
                    for i = h:-1:1
                        ai = u_dict(i);
                        for j = 1:g
                            if isKey(a_dict, i+j)
                                ai = ai + a_dict(i+j) * b_dict(j);
                            else
                                break;
                            end
                        end
                        ai = ai/S;

                        diff = abs(ai-a_dict(i));                
                        a_dict(i) = ai;

                        if max_diff < diff
                            max_diff = diff;
                        end
                    end

                    if max_diff < delta
                        break;
                    end
                end
                
                w0 = 1-sum(cell2mat(values(a_dict)));
                queue_sz_est = 0;
                for i = 1:h
                    queue_sz_est = queue_sz_est + i*a_dict(i);
                end
                queue_sz_est = queue_sz_est/w0;
                
                fprintf('%d iterations for WHF to converge.\n', cnt);
%             disp(rho_est)
%             disp(queue_sz_est)
            end
        end
        
        %% Utilization level
        function rho = rho_est(self, x)
            Q_opt = x;
            if any(size(Q_opt) ~= [self.N, self.hN])
                Q_opt = reshape( Q_opt, [self.N, self.hN] );
            end
            out_byte_rate = self.in_rates * Q_opt * self.out_sizes';
            rho = self.in_byte_rate / out_byte_rate;
        end
        
        %% Privacy level
        function varargout = eps_actual(self, x)
            if nargout > 2
                error("Can't have more than 2 outputs!");
            end
            varargout{1} = eps_s_actual(self,x);
            varargout{2} = eps_t_actual(self,x); % eps_t = 2*eps_o
        end
        function eps_s = eps_s_actual(self, x)
            Q_opt = round(x,15);
            if any(size(Q_opt) ~= [self.N, self.hN])
                Q_opt = reshape( Q_opt, [self.N, self.hN] );
            end
            tt1 = kron( eye(self.N-1), ones(self.N-1,1) ) * Q_opt(2:end,:);
            tt2 = kron( ones(self.N-1,1), eye(self.N-1) ) * Q_opt(2:end,:);
            tmp = tt1 ./ tt2;
            eps_s = log( nanmax(tmp(:)) );
        end
        function eps_t = eps_t_actual(self, x)
            Q_opt = round(x,15);
            if any(size(Q_opt) ~= [self.N, self.hN])
                Q_opt = reshape( Q_opt, [self.N, self.hN] );
            end
            I = eye(self.N);
            tt1 = [kron(ones(self.N,1),I(1,:)); I(2:end,:)] * Q_opt;
            tt2 = [I; kron(ones(self.N-1,1),I(1,:))] * Q_opt;
            tmp = tt1 ./ tt2;
            eps_o = log( nanmax(tmp(:)) );
            eps_t = 2*eps_o;
        end
        
        %% Input byte rate and coefficient of variation
        function mu = get_mu(self, varargin)
            if nargin > 1
                [sizes, rates] = deal(varargin{1}, varargin{2});
            else
                [sizes, rates] = deal(self.in_sizes, self.in_rates);
            end
            if sum(rates) ~= 1
                rates = rates / sum(rates);
            end
            mu = sum( sizes.*rates );
        end
        function cv = get_cv(self, varargin)
            if nargin > 1
                [sizes, rates] = deal(varargin{1}, varargin{2});
            else
                [sizes, rates] = deal(self.in_sizes, self.in_rates);
            end
            mu = self.get_mu(sizes, rates);
            if sum(rates) ~= 1
                rates = rates / sum(rates);
            end
            cv = sqrt(sum((sizes - mu).^2.*rates)) / mu;
        end
        
        %% Find the scale for larger alphabet with the same in_byte_rate
        function ns = get_new_ns(self, n_dev)
            n = 5000;
            new_in_sizes = self.get_in_sizes(n_dev);
            scale_all = linspace(0.01,10,n);
            Bin_all = zeros(n,1);
            CV_all = zeros(n,1);
            for i_scale = 1:length(scale_all)
                new_scale = scale_all(i_scale);
                new_in_rates = self.get_in_rates(new_in_sizes, new_scale);
                Bin_all(i_scale) = self.get_mu(new_in_sizes, new_in_rates);
                CV_all(i_scale) = self.get_cv(new_in_sizes, new_in_rates);
            end
            if strcmp(self.how_to_increase, 'same_bin')
                [~, idx] = min(abs(self.in_byte_rate-Bin_all));
            elseif strcmp(self.how_to_increase, 'same_cv')
                [~, idx] = min(abs(self.coeff_var-CV_all));
            else
                error("'how_to_increase' need to be either 'same_bin' or 'same_cv'!");
            end
%             figure;
%             plot(scale_all, in_byte_rate_all);
%             disp(err);
            closest_scale = scale_all(idx);
            ns = NetSetting('n_dev',n_dev, 'n_eps',self.n_eps, 'n_rho',self.n_rho,... 
                            'pmf',self.pmf, 'scale',closest_scale,...
                            'same_out_sizes',self.same_out_sizes,...
                            'scaled_alphabet',self.scaled_alphabet,...
                            'increased_alphabet',false,...  % temporarily set as false to avoid infinite loop
                            'how_to_increase',self.how_to_increase);
            ns.increased_alphabet = true; % set it back for validation purpose from other files.
            ns.rho_all = self.rho_all;
        end
        
        %% Construct filename for optimization result
        % Get result path
        function rp = get_rp(self, subfolder)
            switch self.where
                case 'local'
                    if self.synthetic
                        rp = strcat('/Users/Sijie/GoogleDrive/PERMIT/3-Simulation Code',...
                                    '/packet size obfuscation/python_code/run_on_cluster/',...
                                    'DiscreteTime/matlab_code_new/res/linspace_rho_all/', subfolder);
                    else
                        rp = strcat('/Users/Sijie/GoogleDrive/PERMIT/3-Simulation Code',...
                                    '/packet size obfuscation/python_code/run_on_cluster/',...
                                    'DiscreteTime/matlab_code_new/res/', subfolder);
                    end
                case 'cluster'
                    rp = '';
            end             
        end
        % Create filename
        function fn = get_fn(self, res_prefix, how)
            if self.synthetic
                fmt1 = strcat(res_prefix,'_%d_dev_%.2f_scale_%2.2f_ibr_%.2f_ar_%.2f_cv');
                fn1 = sprintf(fmt1,self.n_dev,self.scale,self.in_byte_rate,...
                             self.arrival_rate,self.coeff_var);
            else
                fmt1 = strcat(self.device,'_', res_prefix,'_%d_szs_%2.2f_ibr_%.2f_ar_%.2f_cv');
                fn1 = sprintf(fmt1,self.N,self.in_byte_rate,self.arrival_rate,...
                             self.coeff_var);
            end
            fmt2 = '_%d_%d.mat';
            if strcmp(how, 'agg')
                fn2 = sprintf(fmt2, self.n_eps, self.n_rho);
            elseif strcmp(how, 'ind')
                fn2 = sprintf(fmt2, self.i_eps, self.i_rho);
            end
            fn = strcat(fn1, fn2);
        end
        % Create result name = result path + filename
        function rn = get_rn(self, subfolder, res_prefix)
            rp = self.get_rp(subfolder);
            fn = self.get_fn(res_prefix);
            rn = strcat(rp, '/', fn);
        end
        
        %% Get input pkt sizes and rates
        function in_sizes = get_in_sizes(self, varargin)
            if self.synthetic
                log2_minsz = 5;
                if nargin > 1
                    nn_dev = varargin{1};
                else
                    nn_dev = self.n_dev;
                end
                in_sizes = [0,2.^linspace(log2_minsz,log2_minsz+nn_dev-1,nn_dev)];
            else
                in_sizes = self.data.sizes;
            end
        end
        function in_rates = get_in_rates(self, varargin)
            if self.synthetic 
                if nargin > 1
                    [sizes, scale] = deal(varargin{1}, varargin{2});
                else
                    [sizes, scale] = deal(self.in_sizes, self.scale);
                end
                switch self.pmf
                    case 'Zipf'
                        in_rates = Zipf(sizes, scale);
                    case 'Poisson'
                        in_rates = Poisson(sizes, scale);
                    otherwise
                        error('Choose PMF from "Zipf" or "Poisson"!');
                end
            else
                in_rates = self.data.pmf;
            end
        end
        function self = set_data(self)
            path = strcat('/Users/Sijie/GoogleDrive/PERMIT/3-Simulation Code',...
                          '/packet size obfuscation/python_code/other/SmartHomeIoT/',...
                          'matfiles/Filtered/');
            data_fn = strcat(self.device,'.mat');
            switch self.where
                case 'local'
                    self.data = load(strcat(path,data_fn));
                case 'cluster'
                    self.data = load(data_fn);
            end
        end
    end
end

%% Zipf
function pmf = Zipf(sizes, s)
    n_sizes = length(sizes);
    denom = 0;
    for j = 1:n_sizes
        denom = denom + (1/j)^s;
    end
    pmf = 1 ./ (1:n_sizes).^s / denom;
end

%% Poisson
function pmf = Poisson(sizes, r)
    n_sizes = length(sizes);
    k = 1:n_sizes;
    pmf = r.^k * exp(-r) ./ factorial(k);
    pmf = pmf / sum(pmf);
end
