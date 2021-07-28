function [h,familycouplings]=DCAparameters(inputfile,stype)
% Direct Coupling Analysis (DCA)
%
% function dca(inputfile,stype)
% ~
% INPUTS: 
%   inputfile  - file containing the FASTA alignment
%   stype      - species type: 1 for proteins
%			       2 for RNA and DNA
%
% OUTPUTS:
%   h 		    - q x N matrix with the average local fields obtained
%                     from the ones calculated pairwise.
%   Familycouplings - (q*N) x (q*N) matrix with the eij(alpha,beta)
%                     couplings. Each submatrix i,j (size q x q) contains
%                     the couplings for the i,j pair.
%   align - Translated alignment Aminoacid --> Numbers

%
% SOME RELEVANT VARIABLES:
%   N        number of residues in each sequence (no insert)
%   M        number of sequences in the alignment
%   Meff     effective number of sequences after reweighting
%   q        equal to 21 (20 aminoacids + 1 gap) for stype 1
%            equal to 4 for stype 2 (5 if given an alignment with gaps)
%   align    M x N matrix containing the alignmnent
%   Pij_true N x N x q x q matrix containing the reweigthed frequency
%            counts.
%   Pij      N x N x q x q matrix containing the reweighted frequency 
%            counts with pseudo counts.
%   C        N(q-1) x N(q-1) matrix Pij_true,Pi_truecontaining the covariance matrix.
%
%   NOTE: The gauge chosen for this implementation is such that local fields and
%	  couplings of q are fixed to be zero. Zeros are explicitly included on the
%	  corresponding matrices.
%
% Copyright for this implementation: 
%             2011/12 - Andrea Pagnani and Martin Weigt
%                       andrea.pagnani@gmail.com 
%                       martin.weigt@upmc.fr
% 
% Permission is granted for anyone to copy, use, or modify this
% software and accompanying documents for any uncommercial
% purposes, provided this copyright notice is retained, and note is
% made of any changes that have been made. This software and
% documents are distributed without any warranty, express or
% implied. All use is entirely at the user's own risk.
%
% Any publication resulting from applications of DCA should cite:
%
%     F Morcos, A Pagnani, B Lunt, A Bertolino, DS Marks, C Sander, 
%     R Zecchina, JN Onuchic, T Hwa, M Weigt (2011), Direct-coupling
%     analysis of residue co-evolution captures native contacts across 
%     many protein families, Proc. Natl. Acad. Sci. 108:E1293-1301.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    tic
    pseudocount_weight = 0.5; % relative weight of pseudo count   
    theta = 0.2;              % threshold for sequence id in reweighting

    
    [N,M,q,align] = return_alignment(inputfile,stype);


	if (M>80000)
		 tic
		[Pij_true,Pi_true, Meff]=newCompute_True_Frequencies(align,M,N,q, theta);
		t_stop =toc;
		fprintf ( 1, 'Elapsed CPU time frequencies= %f\n', t_stop );
	else
		tic
		[Pij_true,Pi_true, Meff]=oldCompute_True_Frequencies(align,M,N,q, theta);
		t_stop =toc;
		fprintf ( 1, 'Elapsed CPU time frequencies= %f\n', t_stop );
	end

    fprintf('### N = %d M = %d Meff = %.2f q = %d\n', N,M,Meff,q); 
    [Pij,Pi] = with_pc(Pij_true,Pi_true,pseudocount_weight,N,q);
    C = Compute_C(Pij,Pi,N,q);
    invC = inv(C);
    familycouplings=nicematrix(-invC,q);
    Pairwisehfield=Compute_Results(Pi, invC, N, q);
    Pairwisehfield=symmetriclocal(Pairwisehfield,N,q);
    [h,~]=averagehfield(Pairwisehfield);
    t_stop= toc;
    fprintf ( 1, 'Elapsed CPU time= %f\n', t_stop);
end

function [N,M,q,Z] = return_alignment(inputfile,stype)
% reads alignment from inputfile, removes inserts and converts into numbers

    align_full = fastaread(inputfile);
    M = length(align_full);
    ind = align_full(1).Sequence ~= '.' & ...
        align_full(1).Sequence == upper( align_full(1).Sequence );
    N = sum(ind);
    Z = zeros(M,N);
            for i=1:M
                counter = 0;
                for j=1:length(ind)
                    if( ind(j) )
                        counter = counter + 1;
                        Z(i,counter)=letter2number( align_full(i).Sequence(j),stype );
                    end
                end
            end      
    q = max(max(Z));
end

function Pairwisehfield=Compute_Results(Pi,invC, N,q)
% computes and prints the mutual and direct informations

    Pairwisehfield=zeros(N*q,2*N);
    for i=1:N
        for j=(i+1):N
            % direct information from mean-field
            W_mf = ReturnW(invC,i,j,q); 
            Pairwisehfield(((i-1)*q+1):i*q,(2*j-1):2*j)= bp_link(i,j,W_mf,Pi,q);
         end
    end
end

function [Pij_true,Pi_true,Meff] = newCompute_True_Frequencies(align,M,N,q,theta)
% computes reweighted frequency counts
    
    W = ones(1,M);
    if( theta > 0.0 )
        parfor seq=1:M
            for comparing=1:M
                W(seq)=W(seq)+(pdist([align(seq);align(comparing)],'hamm')<theta);
            end
        end
        W = (1./W);
    end
    Meff=sum(W);

    Pij_true = zeros(N,N,q,q);
    Pi_true = zeros(N,q);

    for j=1:M
        for i=1:N
            Pi_true(i,align(j,i)) = Pi_true(i,align(j,i)) + W(j);
        end
    end
    Pi_true = Pi_true/Meff;

    for l=1:M
        for i=1:N-1
            for j=i+1:N
                Pij_true(i,j,align(l,i),align(l,j)) = Pij_true(i,j,align(l,i),align(l,j)) + W(l);
                Pij_true(j,i,align(l,j),align(l,i)) = Pij_true(i,j,align(l,i),align(l,j));
            end
        end
    end
    Pij_true = Pij_true/Meff;

    scra = eye(q,q);
    for i=1:N
        for alpha=1:q
            for beta=1:q
                Pij_true(i,i,alpha,beta) = Pi_true(i,alpha) * scra(alpha,beta);
            end
        end
    end
end

function [Pij_true,Pi_true,Meff] = oldCompute_True_Frequencies(align,M,N,q,theta)
% computes reweighted frequency counts

    W = ones(1,M);
    if( theta > 0.0 )
        W = (1./(1+sum(squareform(pdist(align,'hamm')<theta))));
    end
    Meff=sum(W);

    Pij_true = zeros(N,N,q,q);
    Pi_true = zeros(N,q);

    for j=1:M
        for i=1:N
            Pi_true(i,align(j,i)) = Pi_true(i,align(j,i)) + W(j);
        end
    end
    Pi_true = Pi_true/Meff;

    for l=1:M
        for i=1:N-1
            for j=i+1:N
                Pij_true(i,j,align(l,i),align(l,j)) = Pij_true(i,j,align(l,i),align(l,j)) + W(l);
                Pij_true(j,i,align(l,j),align(l,i)) = Pij_true(i,j,align(l,i),align(l,j));
            end
        end
    end
    Pij_true = Pij_true/Meff;

    scra = eye(q,q);
    for i=1:N
        for alpha=1:q
            for beta=1:q
                Pij_true(i,i,alpha,beta) = Pi_true(i,alpha) * scra(alpha,beta);
            end
        end
    end
end

function x=letter2number(a,stype)
	switch (stype)
	case 1
    switch(a)
        % full AA alphabet
        case '-'
             x=1;
        case 'A'    
            x=2;    
        case 'C'    
            x=3;
        case 'D'
            x=4;
        case 'E'  
            x=5;
        case 'F'
            x=6;
        case 'G'  
            x=7;
        case 'H'
            x=8;
        case 'I'  
            x=9;
        case 'K'
            x=10;
        case 'L'  
            x=11;
        case 'M'
            x=12;
        case 'N'  
            x=13;
        case 'P'
            x=14;
        case 'Q'
            x=15;
        case 'R'
            x=16;
        case 'S'  
            x=17;
        case 'T'
            x=18;
        case 'V'
            x=19;
        case 'W'
            x=20;
        case 'Y'
            x=21;
        otherwise
            x=1;
    end
	case 2
    switch(a)
        % full AA alphabet
        case 'A'
             x=1;
        case 'C'    
            x=2;    
        case 'G'    
            x=3;
        case 'T'
            x=4;
        case 'U'
            x=4;
        case '-'
            x=5;
        otherwise
            x=5;
    end
    end
end

function [Pij,Pi] = with_pc(Pij_true, Pi_true, pseudocount_weight,N,q)
% adds pseudocount

    Pij = (1.-pseudocount_weight)*Pij_true + pseudocount_weight/q/q*ones(N,N,q,q);
    Pi = (1.-pseudocount_weight)*Pi_true + pseudocount_weight/q*ones(N,q);

    scra = eye(q);

    for i=1:N
        for alpha = 1:q
            for beta = 1:q
               Pij(i,i,alpha,beta) =...
                   (1.-pseudocount_weight)*Pij_true(i,i,alpha,beta) + ...
                   pseudocount_weight/q*scra(alpha,beta);
            end
        end
    end 
end

function C = Compute_C(Pij,Pi,N,q)
% computes correlation matrix

    C=zeros(N*(q-1),N*(q-1));
    for i=1:N
        for j=1:N
            for alpha=1:q-1
                for beta=1:q-1
                     C(mapkey(i,alpha,q),mapkey(j,beta,q)) = ...
                         Pij(i,j,alpha,beta) - Pi(i,alpha)*Pi(j,beta);
                end
            end
        end
    end
end

function A=mapkey(i,alpha,q)
%index mapping to the submatrices
    A = (q-1)*(i-1)+alpha;
end

function W=ReturnW(C,i,j,q)
% extracts coupling matrix for columns i and j

    W = ones(q,q);
    W(1:q-1,1:q-1) = exp( -C(mapkey(i,1:q-1,q),mapkey(j,1:q-1,q)) );

end

function hihj = bp_link(i,j,W,P1,q)
% Adjustment to the correct gauge.

    [mu1, mu2] = compute_mu(i,j,W,P1,q);
    mu1=mu1/mu1(q);
    mu2=mu2/mu2(q);
    hihj=[log(mu1'), log(mu2')];
    return;
end

function [mu1,mu2] = compute_mu(i,j,W,P1,q)
%computes e^hi for each pair using message passing algorithm
    epsilon=1e-4;
    diff =1.0;
    mu1 = ones(1,q)/q;
    mu2 = ones(1,q)/q;
    pi = P1(i,:);
    pj = P1(j,:);

    while ( diff > epsilon )

        scra1 = mu2 * W';
        scra2 = mu1 * W;

        new1 = pi./scra1;
        new1 = new1/sum(new1);

        new2 = pj./scra2;
        new2 = new2/sum(new2);

        diff = max( max( abs( new1-mu1 ), abs( new2-mu2 ) ) );

        mu1 = new1;
        mu2 = new2;

    end
end

function symm=symmetriclocal(Localfield,N,q)
%Copy the information to the bottom part
    symm=Localfield;
    for i=1:N
        for j=(i+1):N
            symm(((j-1)*q+1):j*q,(2*i-1):2*i)=[Localfield(((i-1)*q+1):i*q,2*j),...
                Localfield(((i-1)*q+1):i*q,(2*j-1))];
         end
    end
end

function [hi,sigma]=averagehfield(Pairwisehfield)

%   INPUT:
%       Pairwisehfield  -  q*M x 2*M matrix with the Pairwise h local
%                          fields given by the function DCAparameters.
%   OUTPUT:
%           hi           -  (q-1) x M matrix with the average local fields
%                           by site. WITHOUT GAPS (residues displaced one
%                           place in the numeration i-->(i-1)

    N=size(Pairwisehfield,2)/2;
    q=size(Pairwisehfield,1)/N;
    hi=zeros(q,N);
    sigma=zeros(q,N);
    for i=1:N
        if (i==1)
            hi(1:(q),i)=mean(Pairwisehfield(1:q,3:2:2*N),2);
            sigma(1:(q),i)=std(Pairwisehfield(1:q,3:2:2*N),0,2);
        else if (i==N)
                hi(1:(q),i)=mean(Pairwisehfield(((N-1)*q+1):N*q,1:2:(2*N-3)),2);
                sigma(1:(q),i)=std(Pairwisehfield( ((N-1)*q+1):N*q ,1:2:(2*N-3)),0,2);
            else
                hi(1:(q),i)=mean([Pairwisehfield(( ((i-1)*q+1):((i-1)*q+21) ),1:2:(2*(i-1))),...
                    Pairwisehfield( ((i-1)*q+1):i*q ,(2*i+1):2:2*N)],2);
                sigma(1:(q),i)=std([Pairwisehfield(( ((i-1)*q+1):i*q ),1:2:(2*(i-1))),...
                    Pairwisehfield( ((i-1)*q+1):i*q ,(2*i+1):2:2*N)],0,2);
            end
        end
    end
end

function coupling=nicematrix(familycouplinmgs,q)
    n=size(familycouplinmgs,1)/(q-1);
    coupling=zeros(q*n);
    for i=1:n
        for j=1:n
            ii=((i-1)*(q-1)+1):i*(q-1);
            jj=((j-1)*(q-1)+1):j*(q-1);
            newii=((i-1)*q+1):(i*q-1);
            newjj=((j-1)*q+1):(j*q-1);
            coupling(newii,newjj)=familycouplinmgs(ii,jj);
        end
    end
end
