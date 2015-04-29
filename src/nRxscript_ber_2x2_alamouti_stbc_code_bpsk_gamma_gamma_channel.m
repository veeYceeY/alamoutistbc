
function[simBer,Eb_N0_dB,theoryBer_nRx1,theoryBerMRC_nRx2,theoryBerAlamouti_nTx2_nRx1]=nRxscript_ber_2x2_alamouti_stbc_code_bpsk_gamma_gamma_channel(alpha,gamma,N,nRx)

Eb_N0_dB = [0:25]; 
for ii = 1:length(Eb_N0_dB)

    
    ip = rand(1,N)>0.5; 
    s = 2*ip-1; 

   
    sCode = 1/sqrt(2)*kron(reshape(s,2,N/2),ones(1,2)) ;

    % channel
    nch = 1/sqrt(2)*[randGamma(alpha,gamma,N) + j*randGamma(alpha,gamma,N)];
  
    for ki=1:nRx
    h(ki,:)=nch;
    end
    n = 1/sqrt(2)*[randn(nRx,N) + j*randn(nRx,N)];
    y = zeros(nRx,N);
    yMod = zeros(nRx*2,N);
    hMod = zeros(nRx*2,N);
    for kk = 1:nRx

        hMod = kron(reshape(h(kk,:),2,N/2),ones(1,2)); 
        hMod = kron(reshape(h(kk,:),2,N/2),ones(1,2));
        temp = hMod;
        hMod(1,[2:2:end]) = conj(temp(2,[2:2:end])); 
        hMod(2,[2:2:end]) = -conj(temp(1,[2:2:end]));

        
        y(kk,:) = sum(hMod.*sCode,1) + 10^(-Eb_N0_dB(ii)/20)*n(kk,:);

        
        yMod([2*kk-1:2*kk],:) = kron(reshape(y(kk,:),2,N/2),ones(1,2));
    
        
        hEq([2*kk-1:2*kk],:) = hMod;
        hEq(2*kk-1,[1:2:end]) = conj(hEq(2*kk-1,[1:2:end]));
        hEq(2*kk,  [2:2:end]) = conj(hEq(2*kk,  [2:2:end]));

    end

    
    hEqPower = sum(hEq.*conj(hEq),1);
    yHat = sum(hEq.*yMod,1)./hEqPower; 
    yHat(2:2:end) = conj(yHat(2:2:end));

    
    ipHat = real(yHat)>0;

    
    nErr(ii) = size(find([ip- ipHat]),2);

end

simBer = nErr/N; 
EbN0Lin = 10.^(Eb_N0_dB/10);
theoryBer_nRx1 = 0.5.*(1-1*(1+1./EbN0Lin).^(-0.5)); 

p = 1/2 - 1/2*(1+1./EbN0Lin).^(-1/2);
theoryBerMRC_nRx2 = p.^2.*(1+2*(1-p)); 

pAlamouti = 1/2 - 1/2*(1+2./EbN0Lin).^(-1/2);
theoryBerAlamouti_nTx2_nRx1 = pAlamouti.^2.*(1+2*(1-pAlamouti)); 


