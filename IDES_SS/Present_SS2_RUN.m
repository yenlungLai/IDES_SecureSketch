warning('off')
% 'inner' code

k=20; n=k; %dimension of G1

k2=n+20; n2=k2+10;
% 'outer' code k2>=n
G2=randi([0,1],n2,k2);
% 'inner' code
G1=G2(1:n,(k2-n)+1:k2);

delta=1; % this parameter control the hamming distance between w and w'


MSG=randi([0 1],1,k);MSG=gf(MSG);% given an input w=MSG
N_string=randi([1 length(MSG)],1,n2); 
ss=sketching(MSG,N_string,G1,G2);% sketching

% e_vector_dec2=randerr(1,k,1:floor(k/2));
e_vector_dec2=randerr(1,k,delta);
MSG2=MSG+e_vector_dec2; % adversary can sample any random w_e(i)=MSG2
[Recovered_codeword,Noisy_string]=recover(ss,MSG2,N_string,G1,G2);



%%%%%%%%%%% Sketching %%%%%%%%%%%
function [ss]=sketching(MSG,N_string,G1,G2)
n=size(G1,1);
k2=size(G2,2);
msg=ones(1,n); % create galois field codeword of ones
c_star=gf(msg);  %Our proof relies on w_e(i)=w_e(j), means vsyn=c_star
c_star=c_star+MSG;

concat_vstar=zeros(1,k2-n);
vstar=[concat_vstar c_star]; % v^* =0||vsyn,  (size k2)
c_xi = G2*(vstar'); % encoding v^* with 'outer' code
RV1=repmat(HAMMING_HASH(N_string,MSG),1,1); % phi = RV1

ss=c_xi+RV1'; % output sketch SS
end





%%%%%%%%%%% Recovery %%%%%%%%%%%
function [Recovered_codeword,MSG2_e]=recover(ss,MSG2,N_string,G1,G2)

n=size(G1,1);
k=size(G1,2); k2=size(G2,2);

error_weight=0;
while error_weight<=floor(k/2)
    
    ii=1;
    
    if error_weight==0
        error_vct=zeros(1,k);
        nstop=1;
    else
        error_vct=[zeros(1,k-error_weight),ones(1,error_weight)];
        nstop=nchoosek(k,error_weight);
    end
    
    while ii<=nstop
        error_list=[error_vct' nextperms(error_vct,ii)];% double the error parameter to 2eps over chi'
        e_vector_i=error_list(:,ii);
        aa=strcat("Decoding ",num2str(ii)," out of ",num2str(nstop));
        disp(aa)
        
        MSG2_e=MSG2+e_vector_i';%w_e_i=w'xor e, this allow EC up to t2/n2+eps
        RV2=HAMMING_HASH(N_string,MSG2_e); % hamming hash for w_e_i
        
        c_xi_prime=ss-RV2'; % compute the curropted  c'
        c_xi=double(c_xi_prime.x);
        [Decoded_vstar,vld]= gflineq(G2,c_xi); %first decoding
      
        
        if vld==1 && sum(Decoded_vstar(1:k2-n)==0)==k2-n % if first k-n^* elements are all zeros
            vstar_prime=Decoded_vstar(end-(n-1):end); % set the last n^* value of v^*
            vstar_prime=gf(vstar_prime);
            Recovered_codeword=(vstar_prime'-MSG2_e)';
            d9=strcat("1st Decoding successful recovered c*=", int2str(Recovered_codeword.x)');
            
            
            disp(d9)
            
            G1_kn=G1(:,1:(k2-n)); %set desired k2-n (size) for G1
           
            [Decoded_MSG,vld2]= gflineq(G1,double(Recovered_codeword.x)); %try second decoding
            error_weight=floor(k/2);  ii=nstop;
       
            if vld2==0
                d6=strcat("Second decoding got NULL solution");
                d62=strcat("Square mat solution got NULL");
                disp(d6)
                disp(d62)
                
                
            else
                Decoded_MSG=gf(Decoded_MSG);
                d8=strcat("Second decoding get solution w=", int2str(Decoded_MSG.x)');
          
              
                disp(d8)
            
            end    
        elseif vld==0
            Recovered_codeword=0;
            
        else
            Recovered_codeword=0;
        end
        
        ii=ii+1;
    end
    error_weight=error_weight+1;
    
   
end
 if Recovered_codeword==0
        disp('no solution found, n too small')
 end
end




