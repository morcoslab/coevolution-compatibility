function SF=StructureFitness(inputfile,couplings,localfields,Htype)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   AUTHORS Xian-Li Jiang, Rey P. Dimas, Clement T. Y. Chan, Faruck Morcos
%	General Compatibility Score
%	INPUT:
%		inputfile   - file containing the FASTA alignment.
%		couplings   - Corresponding coupling matrix (familycouplings) 
%                 as returned from the DCAparameters function of the 
%			      DCA framework.'eij' in example parameter mat file.
%		localfields - Corresponding h fields matrix as returned
%			      from the DCAparameters function of the DCA
%			      framework.'lf' in example parameter mat file.
%		Htype 	  - Score type: type 1 sums couplings over
%			      pairs across two species while type 2
%			      evaluates the complete Compatibility for the
%			      complete sequence as in the traditional Potts
%			      model.
%
%	OUTPUT:
%		SF 	    - Vector of the Structural Fitness Scores of each sequences.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    InterPairs=load('1lbg_monomer.txt');
    InterPairs=InterPairs(:,:);
    UniqRes=unique(InterPairs);
    [N,M,q,Sequences] = return_alignment(inputfile);
    
    SF=zeros(M,1);
    %Couplings
    switch Htype
        case 1
            %Compute across 2 different modules
              for seq=1:M
                  for res=1:length(InterPairs)
                      iindex=q*(InterPairs(res,1)-1)+Sequences(seq,InterPairs(res,1));
                      jindex=q*(InterPairs(res,2)-1)+Sequences(seq,InterPairs(res,2));
                      SF(seq)=couplings(iindex,jindex)+SF(seq); 
                  end                                       
              end
        case 2
            %Complete Hamiltonian
            for seq=1:M
                  for res=1:N
                      if (res<N)
                          for pair=(res+1):N
                              iindex=q*(res-1)+Sequences(seq,res);
                              jindex=q*(pair-1)+Sequences(seq,pair);
                              SF(seq)=couplings(iindex,jindex)+SF(seq);
                          end
                      end
                  end
            end
    end
    SF=-SF;
end

function [N,M,q,Z] = return_alignment(inputfile)
% reads alignment from inputfile, removes inserts and converts into numbers

    align_full = fastaread(inputfile);
    M = length(align_full);
    ind = align_full(1).Sequence ~= '.' & align_full(1).Sequence == upper( align_full(1).Sequence );
    N = sum(ind);
    Z = zeros(M,N);

    for i=1:M
        counter = 0;
        for j=1:length(ind)
            if( ind(j) )
                counter = counter + 1;
                Z(i,counter)=letter2number( align_full(i).Sequence(j) );
            end
        end
    end
    q = max(max(Z));
end

function x=letter2number(a)
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
end

