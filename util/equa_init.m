function [sig_rigid] = equa_init(vars, mode_idx)
    % Generates the rigid-lid solution given the tropospheric
    % non-dimensional parameters below. Used to initialize the
    % stratospheric solution. mode_idx is the Nth fastest growing mode,
    % (1 means the fastest growing mode, 2 the second, etc.).
    % Adapted from Kerry Emanuel.

    % The parameters below correspond to those in the paper (Khairoutdinov 
    % and Emanuel, 2018), with Greek letter spelled out.
    alpha=vars.alpha;       % Default=1.5
    chi=vars.chi;           % Default=1.5
    C=vars.C;               % Default=0.8
    gamma=vars.gamma;       % Default=1. This should not be smaller than about 0.2.
    D=vars.D;               % Default=1.5
    d=vars.d;               % Default=0.02
    kappa=vars.kappa;       % Default=1
    G=vars.G;               % Default=0.1
    delta=vars.delta;       % Default=30
    ep=vars.ep;             % Precipitation efficiency
    nm=vars.nm;
    km=vars.km;

    %------------------------------------------------------------------------
    N=[-1 0 1 2 3]; % Mode number. Note: One can add more modes (N>3) but
                    % associated Hermite polynomials will have to be added
                    % below if eigenfunction plots are desired.

    close all
    clear hp
    qn=max(size(N));
    k=0:10;
    qk=max(size(k));
    ai=complex(0,1);

    % Below are factors that are common to all the modes. In the notation of
    % equation (9) of the paper: 
    %    a1 = f1*sig^2 + f2*sig + f3
    %    a2 = f4*sig + f5
    %    a3 = f6*sig + f7
    %
    f1=gamma;  
    f2=d.*k.^2-kappa.*C+chi*gamma;
    f3=D*(1+C)+chi*(d.*k.^2-kappa*C);
    f4=gamma*alpha;
    f5=alpha*(d.*k.^2-kappa*C+1+C);
    f6=gamma;
    f7=d.*k.^2-kappa*C+G*(1+C);

    sigsort=zeros(qn,2*qk-1,8);
    sigsortg=zeros(qn,2*qk-1,4);
    sigsortwtg=zeros(qn,2*qk-1);
    hp=zeros(qn,1);
    [~, jm] = min(abs(nm-N));

    % Loop through meridional modes for which N > -1 (The latter are dealt with
    %     separately.)
    for m=jm
        sigwtg=-(N(m).*f5+ai.*k.*f7)./(N(m).*f4+ai.*k.*f6); % WTG dispersion relation
        %
        sig=zeros(qk,8);
        b1=zeros(qk,8);
        b2=zeros(qk,8);
        %
        sigg=zeros(qk,4);
        b1g=zeros(qk,4);
        b2g=zeros(qk,4);
        %
        for i=1:qk
            % Here we use MATLAB's equation solver owing to the algebraic
            % complexity. Note that this is pretty slow. 
            syms x;
            eqn=[((N(m)+0.5)^2)*((f4*x+f5(i))^2+4*x*(f1*x^2+f2(i)*x+f3(i))*(f6*x+f7(i))) == ...
                (ai*k(i)*(f6*x+f7(i))-0.5*(f4*x+f5(i))-(x/delta)*(x*(f1*x^2+f2(i)*x+f3(i)) ...
                +ai*k(i)*(f4*x+f5(i))+k(i)^2*(f6*x+f7(i))))^2]; %#ok<NBRAK>
            S=solve(eqn,x);
            sig(i,:)=vpa(S);
            %
            % This loop tests all 8 roots for viability. Non-viable roots set
            % to zero.
            %
            for j=1:8
                a1=f1*sig(i,j)^2+f2(i)*sig(i,j)+f3(i);
                a2=f4*sig(i,j)+f5(i);
                a3=f6*sig(i,j)+f7(i);
                det=sqrt(a2^2+4*a1*sig(i,j)*a3);
                b1(i,j)=(-a2+det)/(4*sig(i,j)*a3);
                b2(i,j)=(-a2-det)/(4*sig(i,j)*a3);
                if N(m) > -1
                    fac=(sig(i,j)/delta)*(a1*sig(i,j)+ai*k(i)*a2+a3*k(i)^2);
                else 
                    fac=0;
                end    
                res1=ai*k(i)*a3-0.5*a2-(N(m)+0.5)*det-fac;
                res2=ai*k(i)*a3-0.5*a2+(N(m)+0.5)*det-fac;
                %
                if (real(b1(i,j)) <= 0 || (real(b1(i,j)) > 0 && abs(res1) > 1e-4)) && ...
                    (real(b2(i,j)) <= 0 || (real(b2(i,j))> 0 && abs(res2) > 1e-4))
                    sig(i,j)=-1;
                end  
                if N(m) == 0 && abs(sig(i,j)*a1+ai*k(i)*a2+k(i)^2*a3) < 1e-8 % Handle special case when N=0
                    sig(i,j)=-1;
                end    
            end
            %
            % Now calculate 4 roots for geostrophic modes
            %
            A(1)=4*f1.*f6.*(N(m)+0.5).^2; %  A's are coefficients of quartic dispersion relation
            A(2)=4.*(f1.*f7(i)+f2(i).*f6).*(N(m)+0.5).^2;
            A(3)=f4.^2.*N(m)*(N(m)+1)+ai.*k(i).*f4.*f6+f6.^2.*k(i).^2+4.*(f2(i).* ...
                f7(i)+f3(i).*f6).*(N(m)+0.5).^2;
            A(4)=2.*f4.*f5(i).*N(m)*(N(m)+1)+ai.*k(i).*(f4.*f7(i)+f5(i).*f6)+ ...
                2.*k(i).^2.*f6.*f7(i)+4.*f3(i).*f7(i).*(N(m)+0.5).^2;
            A(5)=f5(i).^2.*N(m)*(N(m)+1)+ai.*k(i).*f5(i).*f7(i)+k(i).^2.*f7(i).^2;
            %
            sigg(i,:)=roots(A);   % Note:  This is pretty fast.
            %
            % This loop tests all 4 roots for viability. Non-viable roots set
            % to zero.
            %
            for j=1:4
                a1=f1.*sigg(i,j).^2+f2(i).*sigg(i,j)+f3(i);
                a2=f4.*sigg(i,j)+f5(i);
                a3=f6.*sigg(i,j)+f7(i);
                det=sqrt(a2.^2+4.*a1.*sigg(i,j).*a3);
                b1g(i,j)=(-a2+det)./(4.*sigg(i,j).*a3);
                b2g(i,j)=(-a2-det)./(4.*sigg(i,j).*a3);
                res1=ai.*k(i)*a3-0.5*a2-(N(m)+0.5)*det;
                res2=ai.*k(i)*a3-0.5*a2+(N(m)+0.5)*det;
                %
                if (real(b1g(i,j)) <= 0 || (real(b1g(i,j)) > 0 && abs(res1) > 1e-4)) && ...
                    (real(b2g(i,j)) <= 0 || (real(b2g(i,j))> 0 && abs(res2) > 1e-4))
                    sigg(i,j)=-1;
                end
            end        
        end   
        %
        % Next remap positive k and two-signed frequency into two-signed
        % frequency and two-signed k. (We only use aboslute value of frequency
        % in plots, though.)
        %
        k2=-max(k):max(k);
        qk2=max(size(k2));
        sig2=zeros(qk2,8);
        sig2g=zeros(qk2,4);
        sigwtg2=zeros(1,qk2);
        for i=1:qk
            for j=1:8
              if imag(sig(i,j)) > 0
                  sig2(qk-i+1,j)=sig(i,j);
              else
                  sig2(i+qk-1,j)=sig(i,j);
              end
            end
            %
            for j=1:4
              if imag(sigg(i,j)) > 0
                  sig2g(qk-i+1,j)=sigg(i,j);
              else
                  sig2g(i+qk-1,j)=sigg(i,j);
              end
            end
            %
            if imag(sigwtg(i)) > 0
              sigwtg2(qk-i+1)=sigwtg(i);
            else
              sigwtg2(i+qk-1)=sigwtg(i);
            end  
        end    
        %
        % Sort roots by real part of growth rate
        %
        for i=1:qk2
            sigsort(m,i,:)=sort(sig2(i,:),'descend','ComparisonMethod','real');
            sigsortg(m,i,:)=sort(sig2g(i,:),'descend','ComparisonMethod','real');
            sigsortwtg(m,i)=sigwtg2(i);
        end 
    end
    % Treat v=0 (N=-1) mode seperately
    p(1)=f1;
    for i=1:qk 
        p(2)=f2(i);
        p(3)=f3(i)+f6.*k(i)^2+ai*k(i)*f4;
        p(4)=f7(i)*k(i)^2+ai*k(i)*f5(i);
        sig=roots(p);
        sig(imag(sig)>0)=0; % Roots must have negative imaginary part fo satisfy BCs
        sigsort(1,i+qk-1,1:3)=sort(sig,'descend','ComparisonMethod','real');
    end

    sigsortg(1,:,1:3)=sigsort(1,:,1:3);
    sigsortwtg(1,qk:qk2)=-(-1*f5+ai.*k.*f7)./(-1*f4+ai.*k.*f6); % WTG dispersion relation
    sigsort(sigsort==-1)=NaN.*(1+ai);
    sigsortg(sigsortg==-1)=NaN.*(1+ai);
    sigsortwtg(sigsortwtg==-1)=NaN.*(1+ai);
    msizefac(1,1)=60./max(max(real(sigsort(:,:,1)))); % Set marker size factor according to fastest growing mode

    qn=max(size(N));
    qk2=max(size(k2));
    plotmode='f';
    close all
    clear Ad

    [~,jm]=min(abs(nm-N)); % Calculate index for n value
    [~,im]=min(abs(km-k2)); % Calculate index for k2 value
    [~,imk]=min(abs(abs(km)-k)); % Calculate index for k value

    sigs = squeeze(sigsort(jm,im,:));
    sigs(real(sigs) <= 0) = NaN + 1i * NaN;
    sigs(abs(sigs) < 1e-3) = NaN + 1i * NaN;
    sigs = sigs(~isnan(sigs));

    % Calculate and plot eigenfunctions
    x=vars.x; % Full wavelenth in x
    y=vars.y;   % Note choice of range of y to plot here
    qy=max(size(y));
    qx=max(size(x));
        
    if numel(sigs) == 0
        % No solutions!
        v_xy = zeros(qx,qy);
        u_xy = zeros(qx,qy);
        s_xy = zeros(qx,qy);
        sm_xy =zeros(qx,qy);

        % Save rigid-lid solution into a state-matrix.
        uBT_y = zeros(1, vars.Ny);
        uBC_y = interp1(vars.y, u_xy(1, :), vars.y, 'linear', 0);
        vBT_y = zeros(size(uBT_y));    
        vBC_y = interp1(vars.y, v_xy(1, :), vars.y, 'linear', 0);
        sT_y = interp1(vars.y, s_xy(1, :), vars.y, 'linear', 0);
        smT_y = interp1(vars.y, sm_xy(1, :), vars.y, 'linear', 0);
        uS_yz = zeros(vars.Nz, vars.Ny);
        vS_yz = zeros(vars.Nz, vars.Ny);
        phiS_yz = zeros(vars.Nz, vars.Ny);
        eta_y = zeros(1, vars.Ny);
        qT_y = zeros(1, vars.Ny);
        state_matrix = state_to_matrix(uBT_y, uBC_y, vBT_y, vBC_y, sT_y, smT_y, uS_yz, vS_yz, phiS_yz, eta_y, qT_y);
        sig_rigid = 0 + 0i;
        save(vars.init_file, 'vars', 'state_matrix', 'sig_rigid');
        fprintf('Warning: No rigid-lid growing solution...\nSaved %s!\n', vars.init_file);
    end
    
    for s_idx = 1:numel(sigs)
        sig=sigsort(jm,im,s_idx);
        if km <= -1
            sig = real(sig) + 1i * (-imag(sig));
        end
        %km=abs(km);

        % Recalculate meridional decay factor using chosen root for sig. Then
        % figure out which root of meridional decay factor to use. 
        a1=f1.*sig.^2+f2(imk).*sig+f3(imk);
        a2=f4.*sig+f5(imk);
        a3=f6.*sig+f7(imk);
        det=sqrt(a2^2+4.*a1*sig*a3);
        b1=(-a2+det)/(4*sig*a3);
        b2=(-a2-det)/(4*sig*a3);

        if nm > -1
            fac=(sig/delta)*(a1*sig+ai*km*a2+a3*km^2);
        else 
            fac=0;
        end    
        res1=ai.*km*a3-0.5*a2-(nm+0.5)*det-fac;
        res2=ai.*km*a3-0.5*a2+(nm+0.5)*det-fac;
        if nm > -1 && min(abs(res1),abs(res2)) > 1e-4  % (This should never be invoked!)
            disp(fprintf('Warning: Residual not small, Residual = %f\n', min(abs(res1),abs(res2))));
            continue
        end
        b=b1;
        if abs(res2) < abs(res1) && b2 > 0
            b=b2;
        end 

        if km > 0
            X=exp(ai.*x); % We plot one full wavelength in x regardless of k value
        elseif km < 0
            X=exp(-ai.*x);
        else
            X=zeros(1,qx)+1;
        end    

        Fm=zeros(qx,qy);

        % Now calculate eigenfunctions
        if nm == -1  
            Y=exp(0.5*y.^2.*ai*km/sig);
            for i=1:qx
                for j=1:qy
                    Fm(i,j)=X(i).*Y(j);
                end
            end
            %
            V=zeros(qx,qy);
            U=real(Fm);
            W=real(-ai.*km.*Fm);
            S=real(sig.*Fm./(ai.*km));
            Sm=real(((chi+sig).*(sig.*Fm./(ai.*km))+ai.*km.*Fm+alpha.*Fm)/(1+C));
            Vor=real(-ai.*km.*Fm./sig).*repmat(y,[qx,1]);
            M=W+((1-ep)/ep)*(Sm-chi*S-alpha*U);
        else
            nm1=nm+1;
            Y=exp(-b.*y.^2);
            %
            % First calculate Hermite polynomials and their first derivatives
            %
            H=zeros(4,qy);
            HD=zeros(4,qy);
            H(1,1:qy)=1;
            H(2,1:qy)=y;
            if strcmp(plotmode,'f')
                H(3,1:qy)=0.5.*y.^2+sig*a3/(0.5*a1/b-ai*abs(km)*a3);
                c2=2.*sig*a3/(6*b*sig*a3+2.*a2-ai*km*a3);
            else
                H(3,1:qy)=0.5.*y.^2+sig*a3/(2.*b*sig*a3-ai*km*a3+a2+(sig/delta)*( ...
                    a1*sig+ai*km*a2+km^2*a3));
                c2=2.*sig*a3/(6*b*sig*a3+2.*a2-ai*km*a3+(sig/delta)*( ...
                    a1*sig+ai*km*a2+km^2*a3));
            end    
            H(4,1:qy)=(y.^3)/3+c2.*y;
            HD(1,1:qy)=0;
            HD(2,1:qy)=1;
            HD(3,1:qy)=y;
            HD(4,1:qy)=y.^2+c2;

            % Now calculate eigenfunctions
            sraw=zeros(qx,qy);
            uraw=zeros(qx,qy);
            vorraw=zeros(qx,qy);
            wraw=zeros(qx,qy);
            for i=1:qx
                for j=1:qy
                    Fm(i,j)=X(i)*Y(j)*H(nm1,j);
                    if nm == 0
                        sraw(i,j)=y(j)*X(i)*Y(j)*sig/delta;
                        uraw(i,j)=(ai*km*sraw(i,j)+y(j)*Fm(i,j))/sig;
                        wraw(i,j)=2*b*y(j)*Fm(i,j)-ai*km*uraw(i,j);
                        vorraw(i,j)=-X(i)*Y(j)*(1-4*b*y(j)^2+(4/3)*b^2*y(j)^4)/sig;
                        vorraw(i,j)=vorraw(i,j)+ai*km*Fm(i,j)/delta;
                    else    
                        sraw(i,j)=X(i)*Y(j)*(a3*(HD(nm1,j)-2*b*y(j)*H(nm1,j))+ ...
                            y(j)*H(nm1,j)*(ai*km*a3-a2)/sig)/(a1+ai*km* ...
                            (a2-ai*km*a3)/sig);
                        uraw(i,j)=(ai*km.*sraw(i,j)+y(j)*Fm(i,j))/sig;
                        wraw(i,j)=-(a2*uraw(i,j)+a1*sraw(i,j))/a3;
                        vorraw(i,j)=-(ai*km*y(j)*uraw(i,j)+Fm(i,j)+y(j)*X(i)*Y(j).* ...
                            (HD(nm1,j)-2.*b.*y(j).*H(nm1,j)))/sig;
                        vorraw(i,j)=vorraw(i,j)+ai*km*Fm(i,j)/delta;
                    end
                end
            end   
            S=real(sraw);
            U=real(uraw);
            V=real(Fm);
            Vor=real(vorraw);
            W=real(wraw);
            Sm=real((chi+sig)*sraw+alpha*uraw-(a2*uraw+a1.*sraw)/a3)/(1+C);
            M=W+((1-ep)/ep).*(Sm-chi*S-alpha*U);
        end    

        if s_idx == mode_idx
            if nm == -1
                v_xy = zeros(qx,qy);
                u_xy = Fm;
                s_xy = sig.*Fm./(ai.*km);
                sm_xy =((chi+sig).*(sig.*Fm./(ai.*km))+ai.*km.*Fm+alpha.*Fm)/(1+C); k_mode=km;
            else
                s_xy=sraw; u_xy=uraw; v_xy=Fm; sm_xy=((chi+sig)*sraw+alpha*uraw+wraw)/(1+C); k_mode=km;
            end
            if max(abs(sig*u_xy - (1i * km * s_xy + y.*v_xy)), [], 'all') > 1e-5
                fprintf('Solution error!');
            end

            % Save rigid-lid solution into a state-matrix.
            uBT_y = zeros(1, vars.Ny);
            uBC_y = interp1(vars.y, u_xy(1, :), vars.y, 'linear', 0);
            vBT_y = zeros(size(uBT_y));    
            vBC_y = interp1(vars.y, v_xy(1, :), vars.y, 'linear', 0);
            sT_y = interp1(vars.y, s_xy(1, :), vars.y, 'linear', 0);
            smT_y = interp1(vars.y, sm_xy(1, :), vars.y, 'linear', 0);
            uS_yz = zeros(vars.Nz, vars.Ny);
            vS_yz = zeros(vars.Nz, vars.Ny);
            phiS_yz = zeros(vars.Nz, vars.Ny);
            eta_y = zeros(1, vars.Ny);
            qT_y = zeros(1, vars.Ny);
            state_matrix = state_to_matrix(uBT_y, uBC_y, vBT_y, vBC_y, sT_y, smT_y, uS_yz, vS_yz, phiS_yz, eta_y, qT_y);
            sig_rigid = sig;
            save(vars.init_file, 'vars', 'state_matrix', 'sig_rigid');
            fprintf('Saved %s!\n', vars.init_file);
        end
    end
end
