% indexToF, S�ren Schmidt Jan 28th 2015.
function [Nmax,Umax_out,rmax_out,omega_calc,lambda_calc]=18_indexToF(omega_meas_in,lambda_meas_in,mode,Umax_in,rmax_in)

% Mode 1: find U_max, mode 3: improve results, sampling a smaller region of the
% Rodrigues space around the results previously found, mode 2: plot the fitting
% curves
readeqmat; % for cubic
B=2*pi/2.8665*[1 0 0;0 1 0;0 0 1]; % Fe
h=[1 1 0];
%; 2 0 0; 2 1 1]; % 2 2 0; 3 1 0; 2 2 2; 3 2 1; 4 0 0; 4 2 0; 3 3 2]; % bcc

lambda_min = 2.1;
lambda_max = 4.2;

[idx idy]=find(lambda_meas_in>=lambda_min & lambda_meas_in<=lambda_max);
lambda_meas = lambda_meas_in(idx);
omega_meas = omega_meas_in(idx);

% Set up hkls
n=0;

omega_calc=0;
lambda_calc=0;

for i=1:size(h,1)
    for m=1:24
        E=reshape(Eqmat(m,:),3,3);
        s=E*h(i,:)';
        same=1;
        for j=1:n
           if norm(s'-squeeze(v(j,:)))==0
               same=0;
           end
        end
        if same==1
            n=n+1;
            v(n,:)=s;
        end
    end
end
v=(B*v')';
% 1/v^2
inv = (1./(v(:,1).^2+v(:,2).^2+v(:,3).^2))';

% Sample orientation space

if mode == 3
    r1=(rmax_in(1)-0.1):0.005:(rmax_in(1)+0.1);
    r2=(rmax_in(2)-0.1):0.005:(rmax_in(2)+0.1);
    r3=(rmax_in(3)-0.1):0.005:(rmax_in(3)+0.1);
    % Improve the results sampling a small region of the Rodrigues space

    [X,Y,Z] = meshgrid(r1,r2,r3);
    rx=reshape(X,1,size(X,1)^3);
    ry=reshape(Y,1,size(Y,1)^3);
    rz=reshape(Z,1,size(Z,1)^3);

else
    a=sqrt(2.)-1.;
    %r=-a:0.015:a;
    r=-a:0.02:a;

    [X,Y,Z] = meshgrid(r,r,r);
    rx=reshape(X,1,size(X,1)^3);
    ry=reshape(Y,1,size(Y,1)^3);
    rz=reshape(Z,1,size(Z,1)^3);

end

min(rx)
max(rx)

min(ry)
max(ry)

min(rz)
max(rz)

% define omegas
omega=0:3:177;
%Om =[cosd(omega); -sind(omega)]';  % Right hand
Om =[cosd(omega); sind(omega)]';    % Left hand

size(Om)

Nmax=0;
Smin=1E10;
Umax_out=0;
rmax_out=0;

if mode==1 || mode==3
    nr=size(rx,2);
else
   nr=1;
end

for i=1:nr
    %i
    %size(rx,2)
    %[rx(i) ry(i) rz(i)]
    U=r2U([rx(i) ry(i) rz(i)]);

    if mode==1 || mode==3
        w=U*v';
    else
        w=Umax_in*v';
    end

    p =-4*pi*[inv.*w(1,:); inv.*w(2,:)];

    if mode==2
        figure;
        lambda=Om*p;
        [idx idy]=find(lambda>=lambda_min & lambda<=lambda_max);
        for j=1:size(idx,1)
            hold on
            plot(3*(idx(j)),lambda(idx(j),idy(j)),'.r');
            %plot(3*(idx(j)-1),lambda(idx(j),idy(j)),'.r');
            lambda_calc(j)=lambda(idx(j),idy(j))';
            omega_calc(j) =3*(idx(j))';
        end
        lambda_calc=lambda_calc';
        omega_calc=omega_calc';
        plot(omega_meas,lambda_meas,'.');
    else
        N=0;
        S=0;
        for j=1:size(omega,2)
            lambda=Om(j,:)*p; %prediction
            id=0;
            idp=find(lambda>=lambda_min & lambda<=lambda_max);
            idm=find(omega_meas==omega(j)); %3*(j-1));
            %lambda(idp)
            %lambda_meas(idm)
            if size(lambda(idp),1)>0 && size(lambda_meas(idm),1)>0
                ii=0; jj=0;
                %[ii jj]=find(pdist2(lambda(idp)', lambda_meas(idm))<0.1);
                [ii jj]=find(pdist2(lambda(idp)', lambda_meas(idm))<0.1);

                %pdist2(lambda(idp)', lambda_meas(idm))
                %min(pdist2(lambda(idp)', lambda_meas(idm)))

                s=sum(min(pdist2(lambda(idp)', lambda_meas(idm))));

                %return
                % evaluate both omega AND lambda!

                N=N+size(ii,1);

                S=S+s;
            end
        end
        if S<Smin %N>Nmax
            Nmax=N;
            Smin=S;
            Umax_out=U;
            rmax_out=[rx(i) ry(i) rz(i)];
        end
    end
end
