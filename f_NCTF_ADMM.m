function out = f_NCTF_ADMM(T)
%% NCTF_ADMM
% NCTF based on ADMM
% input: T- data information (struct)
% required
%        T.subject(cell) - data of tensors or matrices
%        T.R (vector)    - component number or rank of each data
%        T.C (vector)    - coupled component number of data
% default
%        T.tol           - stopping tolerance
%        T.eps           - nonnegative papameter
%        T.maxT          - maximum of execution time-10e10;
%        T.maxIter       - maximum of Iteration Number-1000;
%        T.islargescale  - largescale or not (1)
%        T.isnorm        - factor matrix normalize or not(0)
% output: out - results (struct)
%         out.U---------factor matrices of tensor and matrix
%         out.G---------core tensors or matrices/features
%         out.obj-------object function value
%         out.relerr----relative error
%         out.tensorfit-tensorfit
%         out.iter------number of iteration
%         out.Executime-Execution time of current algorithm
% author: Wang xiulin
% last modified by March 10,2021
% note: this code is used for joint analysis of tensor/matrix data
if nargin < 1
    demo_NCTF_ADMM;
    return;
end
t1 = clock;
out = NCTF_ADMM(T);
t2 = clock;
out.RunningTime = etime(t2,t1);
end

function demo_NCTF_ADMM
%% generate data tensor
clear;
clc;
close all;
startup;
%% generation of tensor
for iter = 1:1
    % Params
    S = 2;                 % the number of tensors
    R = zeros(S,1);
%     for s = 1:S
%         R(s) = 1;          % the number of components
%     end
    R(1) = 3;
    R(2) = 4;
    size_tens{1} = [64 130 100 19]; % the size of tensors
    size_tens{2} = [64 130 100 20]; % the size of tensors
    C = [2 2 2 0];            % the coupled number on each mode
    %- coupled tensor generation
    N = numel(size_tens{1});
    UC = cell(N,1);
    for n = 1:N
        UC{n} = max(0, rand(size_tens{1}(n), C(n))); % couple information
    end
    U0 = cell(S,1);
    for p = 1:S
        for n = 1:N
            U0{p}{n} = max(0, rand(size_tens{p}(n), R(p))); % factor matrix-randn
            U0{p}{n}(:,1:C(n)) = UC{n}(:,1:C(n));
        end
        Z.object{p} = full(ktensor(U0{p}));
    end
    Z.R      = R;
    Z.C      = C;
    Z.tol  = 1e-6;
    Z.eps  = 1e-16;
    Z.maxT   =  10e10;
    Z.maxIter = 1000;
    Z.U0 = U0;
    Z.isnorm = 0;
    %%
    tic
    out = NCTF_ADMM(Z);
    out.iter
    out.obj(end)
    out.relerr(end)
    out.tensorfit(end)
    %     out.RunningTime
end
% semilogy(pi)
end
function out = NCTF_ADMM(Z)
% NCTF based on ADMM
% author: Wang xiulin
%% parse parameters
iter = 0;
T         = Z.object; % tensor
for p = 1:numel(T)   
    size_tens{p} = size(T{p}); % size of tensors
end
N         = numel(size_tens{1}); % order of tensors
R         = Z.R;
C         = Z.C;
if isfield(Z,'tol');     tol = Z.tol;         else tol = 1e-8;         end % stopping tolerance
if isfield(Z,'eps');     eps = Z.eps;         else eps = 1e-16;        end % nonnegative papameter
if isfield(Z,'maxIter'); maxIter = Z.maxIter; else maxIter = 1000;     end % max # of iterations
if isfield(Z,'maxT');    maxT = Z.maxT;       else maxT = 10e10;       end % max time in seconds
if isfield(Z,'islargescale') islargescale = Z.islargescale; else islargescale = 1; end; % largescale or not
if isfield(Z,'isnorm'); isnorm = Z.isnorm;    else isnorm = 0;         end % factor matrix normalize or not
%% initialization
obj0 = 0;
Mnrm = cell(numel(T),1);
MaT  = cell(numel(T),N);
U    = cell(numel(T),1);
U_dual = cell(numel(T),1);
U_aux  = cell(numel(T),1);
for p = 1:numel(T)
    Mnrm{p} = norm(T{p});
    obj0 = obj0 + 0.5*Mnrm{p}.^2;
    for n = 1:N
        tensm    = tenmat(T{p},n);
        MaT{p,n} = tensm.data;
        U{p}{n} = max(0, rand(size_tens{p}(n), R(p))); % factor matrix-randn
        U{p}{n}  = U{p}{n}/norm(U{p}{n},'fro')*Mnrm{p}^(1/(N)); % normalize and average
        U_dual{p}{n} = zeros(size(U{p}{n}));
        U_aux{p}{n}  = zeros(size(U{p}{n}));
    end
end
relerr0 = 1;
%% iteration
krrr  = cell(numel(T),1);
krrn  = cell(numel(T),1);
mttkr = cell(numel(T),1);
rho   = zeros(numel(T),1);
start_time = tic;
fprintf('Iteration:     ');
while iter<maxIter
    iter = iter + 1;
    fprintf('\b\b\b\b\b%5d', iter);
    %% Updating factor matrices
    for n = 1:N
        for p = 1:numel(T)
            krrr{p} = khatrirao(U{p}([end:-1:n+1 n-1:-1:1]));
            mttkr{p} = MaT{p,n}*krrr{p};
%           mttkr{p} = mttkrp(tensor(T{p}),U{p},n)*G{p};
            if islargescale == 1
                krrn{p} = Largescalekrrn(U{p},n,N);
            else
                krrn{p} = krrr{p}'*krrr{p};
            end
            rho(p) = trace(krrn{p})/R(p);
            rho(p) = max(rho(p), 1e-12);
        end
        if C(n)
            %% Primal variables-coupled parts of common components
            upC1  = mttkr{1}(:,1:C(n));
            upC2  = U{1}{n}(:,C(n)+1:end)*krrn{1}(C(n)+1:end,1:C(n));
            upC3  = U_dual{1}{n}(:,1:C(n));
            upC4  = rho(1)*U_aux{1}{n}(:,1:C(n));
            downC = krrn{1}(1:C(n),1:C(n)) + rho(1)*eye(C(n));
            for p = 2:numel(T)
                upC1 = upC1 + mttkr{p}(:,1:C(n));
                upC2  = upC2 + U{p}{n}(:,C(n)+1:end)*krrn{p}(C(n)+1:end,1:C(n));
                upC3  = upC3 + U_dual{p}{n}(:,1:C(n));
                upC4  = upC4 + rho(p)*U_aux{p}{n}(:,1:C(n));
                downC = downC + krrn{p}(1:C(n),1:C(n)) + rho(p)*eye(C(n));
            end
            upC = upC1-upC2-upC3 + upC4;
            U{1}{n}(:,1:C(n)) = upC*(pinv(downC));
            if n ~= N && isnorm 
                U{1}{n}(:,1:C(n)) = bsxfun(@rdivide,U{1}{n}(:,1:C(n)),sqrt(sum(U{1}{n}(:,1:C(n)).^2)));
            end
        end
        %% Primal variables-individual parts of common components or individual factors
        for p = 1:numel(T)
            U{p}{n}(:,1:C(n)) = U{1}{n}(:,1:C(n));
            downI = krrn{p}(C(n)+1:end,C(n)+1:end) + rho(p)*eye(R(p)-C(n));
            upI   = mttkr{p}(:,C(n)+1:end)-U{p}{n}(:,1:C(n))*krrn{p}(1:C(n),C(n)+1:end)-U_dual{p}{n}(:,C(n)+1:end) + rho(p)*U_aux{p}{n}(:,C(n)+1:end);
            U{p}{n}(:,C(n)+1:end) = upI*(pinv(downI));
            if n ~= N && isnorm
                U{p}{n}(:,C(n)+1:end) = bsxfun(@rdivide,U{p}{n}(:,C(n)+1:end),sqrt(sum(U{p}{n}(:,C(n)+1:end).^2)));
            end
        end
        %% Auxiliary variables-coupled parts of common components
        if C(n)
            U_aux{1}{n}(:,1:C(n)) = U{p}{n}(:,1:C(n)) + upC3./sum(rho);
        end
        %% Auxiliary variables-individual parts of common components or individual factors
        for p = 1:numel(T)
            U_aux{p}{n}(:,1:C(n))     = U_aux{1}{n}(:,1:C(n));
            U_aux{p}{n}(:,C(n)+1:end) = U{p}{n}(:,C(n)+1:end) + U_dual{p}{n}(:,C(n)+1:end)./rho(p);
            U_aux{p}{n} = max(eps, U_aux{p}{n}); %%%%
            if n ~= N && isnorm 
                U_aux{p}{n} = bsxfun(@rdivide,U_aux{p}{n},sqrt(sum(U_aux{p}{n}.^2)));
            end
        end
        %% Dual variables
        for p = 1:numel(T)
            U_dual{p}{n} = U_dual{p}{n} + rho(p)*(U{p}{n}-U_aux{p}{n});
            U{p}{n} = U_aux{p}{n};
        end
    end
    %% compute the target function value-relative error-tensorfit-
    obj = 0;
    relerr = 0;
    tensorfit = 0;
    for p = 1:numel(T)
        % Computation speed of target function value 1 is faster than 2 for
        % larger-scale problem
        if islargescale == 1
            Usp = U{p}{N}'*U{p}{N};
            obj1 = 0.5*(sum(sum(Usp.*krrn{p}))-2*sum(sum(U{p}{N}.*mttkr{p}))+Mnrm{p}^2);% target function value 1
        else
            Mn = MaT{p,N} - U{p}{N}*krrr{p}';
            obj1 = 0.5*(norm(Mn,'fro').^2); % target function value 2
        end
        relerr1 = (2*obj1)^.5/Mnrm{p}; % relative error
        tensorfit1 = 1 - relerr1;
        obj = obj + obj1;
        relerr = relerr + relerr1;
        tensorfit = tensorfit + tensorfit1;
    end
    out.Executime(iter) = toc(start_time);
    out.obj(iter)       = obj;
    out.relerr(iter)    = relerr;
    out.tensorfit(iter) = tensorfit./numel(T);
    %-stopping criterion-
    relerr2 = abs(relerr-relerr0);  % retative error chage
    if relerr2 < tol || out.Executime(iter) > maxT; break; end
    relerr0 = relerr;
end
out.iter = iter;
out.U = U;
end
function krn = Largescalekrrn(U,n,N)
% kr'*kr
krn = ones(size(U{1},2));
for nn = N:-1:1
    if nn ~= n
        krn = krn.*(U{nn}'*U{nn});
    end
end
end




