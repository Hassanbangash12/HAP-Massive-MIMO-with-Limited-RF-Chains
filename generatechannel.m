function H = generatechannel(fr,rk,B,K,P_t,Ka,c,Ms,Ns)


%%Rician Fading chennal model 

 
%lamda=c/fr;
lamda=1;
dc=2*lamda;       %% spacing 
dr=dc;

startNumber= -pi;
endNumber=pi;
startNumber1=0;
endNumber1=pi/2;

fullH=zeros(K,Ms*Ns);


for k=1:K
        theta=startNumber+(endNumber-startNumber)*rand();
        phi=startNumber1+(endNumber1-startNumber1)*rand();

        %disp(l)
        alpha = ((4*pi*rk)/lamda)^-2;

    d_h = dr*cos(phi)*(sin(theta/lamda));    %%%% take transpose for size problem
    d_v = dc*cos(phi)*cos(theta/lamda);      %%%% take transpose for size problem

    a=zeros(Ns,1);
    a(1,1)=1;
    b=zeros(Ns,1);
    b(1,1)=1;
    for i=1:Ns-1
        a(i+1,1)=exp(1i*i*2*pi*d_h);
        b(i+1,1)=exp(1i*i*2*pi*d_v);
    end

    h_bar=kron(a,b);

    R=zeros(Ns,Ms);    %%%%% first when i do code it was 100 according to paper 

    startNumber2= -pi;
    endNumber2=pi;
    meu=startNumber2+(endNumber2-startNumber2)*rand();
    phi_not=30*(pi/180);
    sym_k=5;
    sigma=10*(pi/180);

    for p=1:Ms*Ns
        for q=1:Ms*Ns
            d1= (p-q)*dr*cos(phi)*sin(theta);
            d2= (p-q)*dc*cos(phi)*cos(theta);

            expo=exp(j*2*pi/0.1250)*(d1+d2);
            %R(p,q)=expo;
            %f_phi=exp(-sqrt(2)*abs(phi-0.5236)/0.1745); %%% abs
            %f_theta=exp(5*cos(theta-meu))/2*pi*(besseli(0,5));
            ftheta=@(theta) exp(5*cos(theta-meu))/2*pi*(besseli(0,5));
            fphi=@(phi) exp(-sqrt(2)*abs(phi-0.5236)/0.1745);
            fthetaint=integral(ftheta,-pi,pi);
            fphiint=integral(fphi,0,pi/2);

            R(p,q)=fthetaint*fphiint*expo;
        end
    end

    W = randn(Ms*Ns,1) + i*randn(Ms*Ns,1);
    h_NLOS=sqrt(R)*W;
    h_one=(sqrt(alpha))*((sqrt(Ka/1+Ka)*h_bar + sqrt(Ka/1+Ka)*h_NLOS));

  fullH(k,:)=h_one';
  
end

H=fullH;

