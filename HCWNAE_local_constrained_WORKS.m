
function [R,J] = HCWNAE_local_constrained_WORKS(x,BC,n_var,D,N,Neq,ns,M,ip,J,J_local,Elb,Erb,SLb,SRb,R,LJ_local,BCtype,xc,yc,zc)
%% State, Costates and BCs
x = reshape(x,[12*N,M]);
%ns = 1;
x1 = x(1:N,:);
x2 = x(N+1:2*N,:);
y1 = x(2*N+1:3*N,:);
y2 = x(3*N+1:4*N,:);
z1 = x(4*N+1:5*N,:);
z2 = x(5*N+1:6*N,:);

L1 = x(6*N+1:7*N,:);
L2 = x(7*N+1:8*N,:);
L3 = x(8*N+1:9*N,:);
L4 = x(9*N+1:10*N,:);
L5 = x(10*N+1:11*N,:);
L6 = x(11*N+1:12*N,:);

x10 = BC(1);
x20 = BC(2);
y10 = BC(3);
y20 = BC(4);
z10 = BC(5);
z20 = BC(6);
x1f = BC(7);
x2f = BC(8);
y1f = BC(9);
y2f = BC(10);
z1f = BC(11);
z2f = BC(12);
q11 = 0;
q22 = 0;
q33 = 0;
q44 = 0;
q55 = 0;
q66 = 0;


%% System of Equations
for k = 1 : M % for each cover
      
      if k ==2
         %Residuals   
        R(1:N,k)       = D(:,:,k)*x1(:,k) - x2(:,k); 
        R(N+1:2*N,k)   = D(:,:,k)*x2(:,k) + 3*n_var*n_var*x1(:,k);
        R(2*N+1:3*N,k) = D(:,:,k)*y1(:,k) - y2(:,k);
        R(3*N+1:4*N,k) = D(:,:,k)*y2(:,k);    
        R(4*N+1:5*N,k) = D(:,:,k)*z1(:,k) - z2(:,k);
        R(5*N+1:6*N,k) = D(:,:,k)*z2(:,k) - n_var*n_var*z1(:,k);

        R(6*N+1:7*N,k) = D(:,:,k)*L1(:,k) + q11*x1(:,k) - 3*n_var*n_var*L2(:,k);
        R(7*N+1:8*N,k) = D(:,:,k)*L2(:,k) + q22*x2(:,k) + L1(:,k);
        R(8*N+1:9*N,k) = D(:,:,k)*L3(:,k) + q33*y1(:,k);
        R(9*N+1:10*N,k) = D(:,:,k)*L4(:,k) + q44*y2(:,k) + L3(:,k);
        R(10*N+1:11*N,k) = D(:,:,k)*L5(:,k) + q55*z1(:,k) + n_var*n_var*L6(:,k);
        R(11*N+1:12*N,k) = D(:,:,k)*L6(:,k) + q66*z2(:,k) + L5(:,k);
        %first attempt
        % R(1:N,k)       = D(:,:,k)*x1(:,k) - x2(:,k); 
        % R(N+1:2*N,k)   = D(:,:,k)*x2(:,k) - x2(:,k).^2./(xc-x1(:,k));
        % R(2*N+1:3*N,k) = D(:,:,k)*y1(:,k) - y2(:,k);
        % R(3*N+1:4*N,k) = D(:,:,k)*y2(:,k) - y2(:,k).^2./(yc-y1(:,k));    
        % R(4*N+1:5*N,k) = D(:,:,k)*z1(:,k) - z2(:,k);
        % R(5*N+1:6*N,k) = D(:,:,k)*z2(:,k) - z2(:,k).^2./(zc-z1(:,k));
        % 
        % R(6*N+1:7*N,k) = D(:,:,k)*L1(:,k) + q11*x1(:,k) + L2(:,k).*x2(:,k).^2./(xc-x1(:,k)).^2;
        % R(7*N+1:8*N,k) = D(:,:,k)*L2(:,k) + q22*x2(:,k) + L1(:,k) + 2*L2(:,k).*x2(:,k)./(xc-x1(:,k));
        % R(8*N+1:9*N,k) = D(:,:,k)*L3(:,k) + q33*y1(:,k) + L4(:,k).*y2(:,k).^2./(yc-y1(:,k)).^2;
        % R(9*N+1:10*N,k) = D(:,:,k)*L4(:,k) + q44*y2(:,k) + L3(:,k) +2*L4(:,k).*y2(:,k)./(yc-y1(:,k));
        % R(10*N+1:11*N,k) = D(:,:,k)*L5(:,k) + q55*z1(:,k) + L6(:,k).*z2(:,k).^2./(zc-z1(:,k)).^2;
        % R(11*N+1:12*N,k) = D(:,:,k)*L6(:,k) + q66*z2(:,k) + L5(:,k) + 2*L6(:,k).*z2(:,k)./(zc-z1(:,k));

        %Jacobian

        %good
        dR1dx1 = D(:,:,k);              dR1dx2 = -eye(N);               dR1dy1 = zeros(N); dR1dy2 = zeros(N);         dR1dz1 = zeros(N);           dR1dz2 = zeros(N);
            dR1dl1 = zeros(N);              dR1dl2 = zeros(N);              dR1dl3 = zeros(N); dR1dl4 = zeros(N);         dR1dl5 = zeros(N);           dR1dl6 = zeros(N);
        %good
        dR2dx1 = 3*n_var*n_var*eye(N);     dR2dx2 = D(:,:,k);              dR2dy1 = zeros(N); dR2dy2 = zeros(N);  dR2dz1 = zeros(N);           dR2dz2 = zeros(N);
            dR2dl1 = zeros(N);              dR2dl2 = zeros(N);            dR2dl3 = zeros(N); dR2dl4 = zeros(N);         dR2dl5 = zeros(N);           dR2dl6 = zeros(N);
        %good
        dR3dx1 = zeros(N);              dR3dx2 = zeros(N);              dR3dy1 = D(:,:,k); dR3dy2 = -eye(N);          dR3dz1 = zeros(N);           dR3dz2 = zeros(N);
            dR3dl1 = zeros(N);              dR3dl2 = zeros(N);              dR3dl3 = zeros(N); dR3dl4 = zeros(N);         dR3dl5 = zeros(N);           dR3dl6 = zeros(N);
        %good
        dR4dx1 = zeros(N);              dR4dx2 =zeros(N);               dR4dy1 = zeros(N);  dR4dy2 = D(:,:,k);         dR4dz1 = zeros(N);           dR4dz2 = zeros(N);
            dR4dl1 = zeros(N);              dR4dl2 = zeros(N);              dR4dl3 = zeros(N); dR4dl4 = zeros(N);       dR4dl5 = zeros(N);           dR4dl6 = zeros(N);
        %good
        dR5dx1 = zeros(N);              dR5dx2 = zeros(N);              dR5dy1 = zeros(N); dR5dy2 = zeros(N);         dR5dz1 = D(:,:,k);           dR5dz2 = -eye(N);
            dR5dl1 = zeros(N);              dR5dl2 = zeros(N);              dR5dl3 = zeros(N); dR5dl4 = zeros(N);         dR5dl5 = zeros(N);           dR5dl6 = zeros(N);
        %good
        dR6dx1 = zeros(N);              dR6dx2 = zeros(N);              dR6dy1 = zeros(N); dR6dy2 = zeros(N);         dR6dz1 = -n_var*n_var*eye(N); dR6dz2 = D(:,:,k);
            dR6dl1 = zeros(N);              dR6dl2 = zeros(N);              dR6dl3 = zeros(N); dR6dl4 = zeros(N);         dR6dl5 = zeros(N);           dR6dl6 = zeros(N);
        %good
        dR7dx1 = q11*eye(N);              dR7dx2 = zeros(N);              dR7dy1 = zeros(N); dR7dy2 = zeros(N);         dR7dz1 = zeros(N);           dR7dz2 = zeros(N);
            dR7dl1 = D(:,:,k);              dR7dl2 = -3*n_var*n_var*eye(N); dR7dl3 = zeros(N); dR7dl4 = zeros(N);         dR7dl5 = zeros(N);           dR7dl6 = zeros(N);
        %good
        dR8dx1 = zeros(N);              dR8dx2 = q22*eye(N);              dR8dy1 = zeros(N); dR8dy2 = zeros(N);         dR8dz1 = zeros(N);           dR8dz2 = zeros(N);
            dR8dl1 = eye(N);               dR8dl2 = D(:,:,k);              dR8dl3 = zeros(N); dR8dl4 = zeros(N);   dR8dl5 = zeros(N);           dR8dl6 = zeros(N);
        %good
        dR9dx1 = zeros(N);              dR9dx2 = zeros(N);              dR9dy1 = q33*eye(N);    dR9dy2 = zeros(N);         dR9dz1 = zeros(N);           dR9dz2 = zeros(N);
            dR9dl1 = zeros(N);              dR9dl2 = zeros(N);              dR9dl3 = D(:,:,k); dR9dl4 = zeros(N);         dR9dl5 = zeros(N);           dR9dl6 = zeros(N);
        %
        dR10dx1 = zeros(N);             dR10dx2 = zeros(N);             dR10dy1 = zeros(N); dR10dy2 = q44*eye(N);        dR10dz1 = zeros(N);          dR10dz2 = zeros(N);
            dR10dl1 = zeros(N);             dR10dl2 = zeros(N);      dR10dl3 = eye(N); dR10dl4 = D(:,:,k);               dR10dl5 = zeros(N);          dR10dl6 = zeros(N);
        %
        dR11dx1 = zeros(N);             dR11dx2 = zeros(N);             dR11dy1 = zeros(N);dR11dy2 = zeros(N);        dR11dz1 = q55*eye(N);          dR11dz2 = zeros(N);
            dR11dl1 = zeros(N);             dR11dl2 = zeros(N);             dR11dl3 = zeros(N);dR11dl4 = zeros(N);        dR11dl5 = D(:,:,k);          dR11dl6 = n_var*n_var*eye(N);
        %
        dR12dx1 = zeros(N);             dR12dx2 = zeros(N);             dR12dy1 = zeros(N);dR12dy2 = zeros(N);        dR12dz1 = zeros(N);          dR12dz2 = q66*eye(N);
            dR12dl1 = zeros(N);             dR12dl2 = zeros(N);             dR12dl3 = zeros(N);dR12dl4 = zeros(N);        dR12dl5 = eye(N);           dR12dl6 = D(:,:,k);


        % %first attempt
        % dR1dx1 = D(:,:,k);              dR1dx2 = -eye(N);               dR1dy1 = zeros(N); dR1dy2 = zeros(N);         dR1dz1 = zeros(N);           dR1dz2 = zeros(N);
        %     dR1dl1 = zeros(N);              dR1dl2 = zeros(N);              dR1dl3 = zeros(N); dR1dl4 = zeros(N);         dR1dl5 = zeros(N);           dR1dl6 = zeros(N);
        % %good
        % dR2dx1 = - x2(:,k).^2/(xc-x1(:,k)).^2; dR2dx2 = D(:,:,k) - 2*- x2(:,k)/(xc-x1(:,k));              dR2dy1 = zeros(N); dR2dy2 = zeros(N);  dR2dz1 = zeros(N);           dR2dz2 = zeros(N);
        %     dR2dl1 = zeros(N);              dR2dl2 = zeros(N);            dR2dl3 = zeros(N); dR2dl4 = zeros(N);         dR2dl5 = zeros(N);           dR2dl6 = zeros(N);
        % %good
        % dR3dx1 = zeros(N);              dR3dx2 = zeros(N);              dR3dy1 = D(:,:,k); dR3dy2 = -eye(N);          dR3dz1 = zeros(N);           dR3dz2 = zeros(N);
        %     dR3dl1 = zeros(N);              dR3dl2 = zeros(N);              dR3dl3 = zeros(N); dR3dl4 = zeros(N);         dR3dl5 = zeros(N);           dR3dl6 = zeros(N);
        % %good
        % dR4dx1 = zeros(N);              dR4dx2 =zeros(N);               dR4dy1 = - y2(:,k).^2/(yc-y1(:,k)).^2; dR4dy2 = D(:,:,k) - -2*y2(:,k)/(yc-y1(:,k));         dR4dz1 = zeros(N);           dR4dz2 = zeros(N);
        %     dR4dl1 = zeros(N);              dR4dl2 = zeros(N);              dR4dl3 = zeros(N); dR4dl4 = zeros(N);       dR4dl5 = zeros(N);           dR4dl6 = zeros(N);
        % %good
        % dR5dx1 = zeros(N);              dR5dx2 = zeros(N);              dR5dy1 = zeros(N); dR5dy2 = zeros(N);         dR5dz1 = D(:,:,k);           dR5dz2 = -eye(N);
        %     dR5dl1 = zeros(N);              dR5dl2 = zeros(N);              dR5dl3 = zeros(N); dR5dl4 = zeros(N);         dR5dl5 = zeros(N);           dR5dl6 = zeros(N);
        % %good
        % dR6dx1 = zeros(N);              dR6dx2 = zeros(N);              dR6dy1 = zeros(N); dR6dy2 = zeros(N);         dR6dz1 = - z2(:,k).^2/(zc-z1(:,k)).^2; dR6dz2 = D(:,:,k) - 2*z2(:,k)/(zc-z1(:,k));
        %     dR6dl1 = zeros(N);              dR6dl2 = zeros(N);              dR6dl3 = zeros(N); dR6dl4 = zeros(N);         dR6dl5 = zeros(N);           dR6dl6 = zeros(N);
        % %good
        % dR7dx1 = q11*eye(N) + 2*L2(:,k).*x2(:,k).^2/(xc-x1(:,k)).^3;              dR7dx2 = 2*L2(:,k) .* x2(:,k)/(xc-x1(:,k)).^2;              dR7dy1 = zeros(N); dR7dy2 = zeros(N);         dR7dz1 = zeros(N);           dR7dz2 = zeros(N);
        %     dR7dl1 = D(:,:,k);              dR7dl2 = x2(:,k).^2/(xc-x1(:,k)).^2; dR7dl3 = zeros(N); dR7dl4 = zeros(N);         dR7dl5 = zeros(N);           dR7dl6 = zeros(N);
        % %good
        % dR8dx1 = 2*L2(:,k).*x2(:,k)/(xc-x1(:,k)).^2;              dR8dx2 = q22*eye(N) + 2*L2(:,k)/(xc-x1(:,k));              dR8dy1 = zeros(N); dR8dy2 = zeros(N);         dR8dz1 = zeros(N);           dR8dz2 = zeros(N);
        %     dR8dl1 = eye(N);               dR8dl2 = D(:,:,k) + 2*x2(:,k)/(xc-x1(:,k));              dR8dl3 = zeros(N); dR8dl4 = zeros(N);   dR8dl5 = zeros(N);           dR8dl6 = zeros(N);
        % %good
        % dR9dx1 = zeros(N);              dR9dx2 = zeros(N);              dR9dy1 = q33*eye(N) + 2*L4(:,k).*y2(:,k).^2/(yc-y1(:,k)).^3; dR9dy2 = 2*L4(:,k).*y2(:,k)/(yc-y1(:,k)).^2;         dR9dz1 = zeros(N);           dR9dz2 = zeros(N);
        %     dR9dl1 = zeros(N);              dR9dl2 = zeros(N);              dR9dl3 = D(:,:,k); dR9dl4 = y2(:,k).^2/(yc-y1(:,k)).^2;         dR9dl5 = zeros(N);           dR9dl6 = zeros(N);
        % %
        % dR10dx1 = zeros(N);             dR10dx2 = zeros(N);             dR10dy1 = 2*L4(:,k).*y2(:,k)/(yc-y1(:,k)).^2; dR10dy2 = q44*eye(N) + 2*L4(:,k)/(yc-y1(:,k));        dR10dz1 = zeros(N);          dR10dz2 = zeros(N);
        %     dR10dl1 = zeros(N);             dR10dl2 = zeros(N);      dR10dl3 = eye(N); dR10dl4 = D(:,:,k) + 2*y2(:,k)/(yc-y1(:,k));               dR10dl5 = zeros(N);          dR10dl6 = zeros(N);
        % %
        % dR11dx1 = zeros(N);             dR11dx2 = zeros(N);             dR11dy1 = zeros(N);dR11dy2 = zeros(N);        dR11dz1 = q55*eye(N) + 2*L6(:,k).*z2(:,k).^2/(zc-z1(:,k)).^3;          dR11dz2 = 2*L6(:,k).*z2(:,k)/(zc-z1(:,k)).^2;
        %     dR11dl1 = zeros(N);             dR11dl2 = zeros(N);             dR11dl3 = zeros(N);dR11dl4 = zeros(N);        dR11dl5 = D(:,:,k);          dR11dl6 = z2(:,k).^2/(zc-z1(:,k)).^2;
        % %
        % dR12dx1 = zeros(N);             dR12dx2 = zeros(N);             dR12dy1 = zeros(N);dR12dy2 = zeros(N);        dR12dz1 = 2*L6(:,k).*z2(:,k)/(zc-z1(:,k)).^2;          dR12dz2 = q66*eye(N) + 2*L6(:,k)/(zc-z1(:,k));
        %     dR12dl1 = zeros(N);             dR12dl2 = zeros(N);             dR12dl3 = zeros(N);dR12dl4 = zeros(N);        dR12dl5 = eye(N);           dR12dl6 = D(:,:,k) + 2*z2(:,k)./(zc-z1(:,k));
     else

                 %Residuals    
        R(1:N,k)       = D(:,:,k)*x1(:,k) - x2(:,k); 
        R(N+1:2*N,k)   = D(:,:,k)*x2(:,k) - (3*n_var*n_var*x1(:,k)+2*n_var*y2(:,k)-L2(:,k));
        R(2*N+1:3*N,k) = D(:,:,k)*y1(:,k) - y2(:,k);
        R(3*N+1:4*N,k) = D(:,:,k)*y2(:,k) + (2*n_var*x2(:,k)+L4(:,k));    
        R(4*N+1:5*N,k) = D(:,:,k)*z1(:,k) - z2(:,k);
        R(5*N+1:6*N,k) = D(:,:,k)*z2(:,k) + (n_var*n_var*z1(:,k)+L6(:,k));
        
        R(6*N+1:7*N,k) = D(:,:,k)*L1(:,k) + q11*x1(:,k) + 3*n_var*n_var*L2(:,k);
        R(7*N+1:8*N,k) = D(:,:,k)*L2(:,k) + q22*x2(:,k) + L1(:,k)-2*n_var*L4(:,k);
        R(8*N+1:9*N,k) = D(:,:,k)*L3(:,k) + q33*y1(:,k);
        R(9*N+1:10*N,k) = D(:,:,k)*L4(:,k) + q44*y2(:,k) +2*n_var*L2(:,k)+L3(:,k);
        R(10*N+1:11*N,k) = D(:,:,k)*L5(:,k) + q55*z1(:,k) - n_var*n_var*L6(:,k);
        R(11*N+1:12*N,k) = D(:,:,k)*L6(:,k) + q66*z2(:,k) + L5(:,k);

        %Jacobian
        dR1dx1 = D(:,:,k);              dR1dx2 = -eye(N);               dR1dy1 = zeros(N); dR1dy2 = zeros(N);         dR1dz1 = zeros(N);           dR1dz2 = zeros(N);
            dR1dl1 = zeros(N);              dR1dl2 = zeros(N);              dR1dl3 = zeros(N); dR1dl4 = zeros(N);         dR1dl5 = zeros(N);           dR1dl6 = zeros(N);
        dR2dx1 = -3*n_var*n_var*eye(N); dR2dx2 = D(:,:,k);              dR2dy1 = zeros(N); dR2dy2 = -2*n_var*eye(N);  dR2dz1 = zeros(N);           dR2dz2 = zeros(N);
            dR2dl1 = zeros(N);              dR2dl2 = eye(N);            dR2dl3 = zeros(N); dR2dl4 = zeros(N);         dR2dl5 = zeros(N);           dR2dl6 = zeros(N);
        dR3dx1 = zeros(N);              dR3dx2 = zeros(N);              dR3dy1 = D(:,:,k); dR3dy2 = -eye(N);          dR3dz1 = zeros(N);           dR3dz2 = zeros(N);
            dR3dl1 = zeros(N);              dR3dl2 = zeros(N);              dR3dl3 = zeros(N); dR3dl4 = zeros(N);         dR3dl5 = zeros(N);           dR3dl6 = zeros(N);
        dR4dx1 = zeros(N);              dR4dx2 = 2*n_var*eye(N);        dR4dy1 = zeros(N); dR4dy2 = D(:,:,k);         dR4dz1 = zeros(N);           dR4dz2 = zeros(N);
            dR4dl1 = zeros(N);              dR4dl2 = zeros(N);              dR4dl3 = zeros(N); dR4dl4 = eye(N);       dR4dl5 = zeros(N);           dR4dl6 = zeros(N);
        dR5dx1 = zeros(N);              dR5dx2 = zeros(N);              dR5dy1 = zeros(N); dR5dy2 = zeros(N);         dR5dz1 = D(:,:,k);           dR5dz2 = -eye(N);
            dR5dl1 = zeros(N);              dR5dl2 = zeros(N);              dR5dl3 = zeros(N); dR5dl4 = zeros(N);         dR5dl5 = zeros(N);           dR5dl6 = zeros(N);
        dR6dx1 = zeros(N);              dR6dx2 = zeros(N);              dR6dy1 = zeros(N); dR6dy2 = zeros(N);         dR6dz1 = n_var*n_var*eye(N); dR6dz2 = D(:,:,k);
            dR6dl1 = zeros(N);              dR6dl2 = zeros(N);              dR6dl3 = zeros(N); dR6dl4 = zeros(N);         dR6dl5 = zeros(N);           dR6dl6 = eye(N);
        dR7dx1 = q11*eye(N);              dR7dx2 = zeros(N);              dR7dy1 = zeros(N); dR7dy2 = zeros(N);         dR7dz1 = zeros(N);           dR7dz2 = zeros(N);
            dR7dl1 = D(:,:,k);              dR7dl2 = 3*n_var*n_var*eye(N); dR7dl3 = zeros(N); dR7dl4 = zeros(N);         dR7dl5 = zeros(N);           dR7dl6 = zeros(N);
        dR8dx1 = zeros(N);              dR8dx2 = q22*eye(N);              dR8dy1 = zeros(N); dR8dy2 = zeros(N);         dR8dz1 = zeros(N);           dR8dz2 = zeros(N);
            dR8dl1 = eye(N);               dR8dl2 = D(:,:,k);              dR8dl3 = zeros(N); dR8dl4 = -2*n_var*eye(N);   dR8dl5 = zeros(N);           dR8dl6 = zeros(N);
        dR9dx1 = zeros(N);              dR9dx2 = zeros(N);              dR9dy1 = q33*eye(N); dR9dy2 = zeros(N);         dR9dz1 = zeros(N);           dR9dz2 = zeros(N);
            dR9dl1 = zeros(N);              dR9dl2 = zeros(N);              dR9dl3 = D(:,:,k); dR9dl4 = zeros(N);         dR9dl5 = zeros(N);           dR9dl6 = zeros(N);
        dR10dx1 = zeros(N);             dR10dx2 = zeros(N);             dR10dy1 = zeros(N); dR10dy2 = q44*eye(N);        dR10dz1 = zeros(N);          dR10dz2 = zeros(N);
            dR10dl1 = zeros(N);             dR10dl2 = 2*n_var*eye(N);      dR10dl3 = eye(N); dR10dl4 = D(:,:,k);               dR10dl5 = zeros(N);          dR10dl6 = zeros(N);
        dR11dx1 = zeros(N);             dR11dx2 = zeros(N);             dR11dy1 = zeros(N);dR11dy2 = zeros(N);        dR11dz1 = q55*eye(N);          dR11dz2 = zeros(N);
            dR11dl1 = zeros(N);             dR11dl2 = zeros(N);             dR11dl3 = zeros(N);dR11dl4 = zeros(N);        dR11dl5 = D(:,:,k);          dR11dl6 = -n_var*n_var*eye(N);
        dR12dx1 = zeros(N);             dR12dx2 = zeros(N);             dR12dy1 = zeros(N);dR12dy2 = zeros(N);        dR12dz1 = zeros(N);          dR12dz2 = q66*eye(N);
            dR12dl1 = zeros(N);             dR12dl2 = zeros(N);             dR12dl3 = zeros(N);dR12dl4 = zeros(N);        dR12dl5 = eye(N);           dR12dl6 = D(:,:,k);
 
      end


      %end
  
    % Jacobian boundary conditions
    if k == 1
                dR1dx1(1,:) = [1,zeros(1,N-1)]; dR1dx2(1,:) = zeros(1,N);       dR1dy1(1,:) = zeros(1,N); dR1dy2(1,:) = zeros(1,N); dR1dz1(1,:) = zeros(1,N); dR1dz2(1,:) = zeros(1,N); ... 
            dR1dl1(1,:) = zeros(1,N); dR1dl2(1,:) = zeros(1,N);dR1dl3(1,:) = zeros(1,N);dR1dl4(1,:) = zeros(1,N);dR1dl5(1,:) = zeros(1,N);dR1dl6(1,:) = zeros(1,N); 

        dR2dx1(1,:) = zeros(1,N);       dR2dx2(1,:) = [1,zeros(1,N-1)]; dR2dy1(1,:) = zeros(1,N); dR2dy2(1,:) = zeros(1,N); dR2dz1(1,:) = zeros(1,N); dR2dz2(1,:) = zeros(1,N);...
            dR2dl1(1,:) = zeros(1,N); dR2dl2(1,:) = zeros(1,N);dR2dl3(1,:) = zeros(1,N);dR2dl4(1,:) = zeros(1,N);dR2dl5(1,:) = zeros(1,N);dR2dl6(1,:) = zeros(1,N);

        dR3dx1(1,:) = zeros(1,N); dR3dx2(1,:) = zeros(1,N);       dR3dy1(1,:) = [1,zeros(1,N-1)]; dR3dy2(1,:) = zeros(1,N); dR3dz1(1,:) = zeros(1,N); dR3dz2(1,:) = zeros(1,N);...
            dR3dl1(1,:) = zeros(1,N); dR3dl2(1,:) = zeros(1,N);dR3dl3(1,:) = zeros(1,N);dR3dl4(1,:) = zeros(1,N);dR3dl5(1,:) = zeros(1,N);dR3dl6(1,:) = zeros(1,N);

        dR4dx1(1,:) = zeros(1,N);       dR4dx2(1,:) = zeros(1,N); dR4dy1(1,:) = zeros(1,N); dR4dy2(1,:) = [1,zeros(1,N-1)]; dR4dz1(1,:) = zeros(1,N); dR4dz2(1,:) = zeros(1,N);...
            dR4dl1(1,:) = zeros(1,N); dR4dl2(1,:) = zeros(1,N);dR4dl3(1,:) = zeros(1,N);dR4dl4(1,:) = zeros(1,N);dR4dl5(1,:) = zeros(1,N);dR4dl6(1,:) = zeros(1,N);

        dR5dx1(1,:) = zeros(1,N); dR5dx2(1,:) = zeros(1,N);       dR5dy1(1,:) = zeros(1,N); dR5dy2(1,:) = zeros(1,N); dR5dz1(1,:) = [1,zeros(1,N-1)]; dR5dz2(1,:) = zeros(1,N);...
            dR5dl1(1,:) = zeros(1,N); dR5dl2(1,:) = zeros(1,N);dR5dl3(1,:) = zeros(1,N);dR5dl4(1,:) = zeros(1,N);dR5dl5(1,:) = zeros(1,N);dR5dl6(1,:) = zeros(1,N);

        dR6dx1(1,:) = zeros(1,N);       dR6dx2(1,:) = zeros(1,N); dR6dy1(1,:) = zeros(1,N); dR6dy2(1,:) = zeros(1,N); dR6dz1(1,:) = zeros(1,N); dR6dz2(1,:) = [1,zeros(1,N-1)]; ...
            dR6dl1(1,:) = zeros(1,N); dR6dl2(1,:) = zeros(1,N);dR6dl3(1,:) = zeros(1,N);dR6dl4(1,:) = zeros(1,N);dR6dl5(1,:) = zeros(1,N);dR6dl6(1,:) = zeros(1,N);

     end
    if k == M
        %         dR7dx1(N,:) = [zeros(1,N-1),1]; dR7dx2(N,:) = zeros(1,N);       dR7dy1(N,:) = zeros(1,N); dR7dy2(N,:) = zeros(1,N); dR7dz1(N,:) = zeros(1,N); dR7dz2(N,:) = zeros(1,N); ... 
        %     dR7dl1(N,:) = zeros(1,N); dR7dl2(N,:) = zeros(1,N);dR7dl3(N,:) = zeros(1,N);dR7dl4(N,:) = zeros(1,N);dR7dl5(N,:) = zeros(1,N);dR7dl6(N,:) = zeros(1,N); 
        % 
        % dR8dx1(N,:) = zeros(1,N);       dR8dx2(N,:) = [zeros(1,N-1),1]; dR8dy1(N,:) = zeros(1,N); dR8dy2(N,:) = zeros(1,N); dR8dz1(N,:) = zeros(1,N); dR8dz2(N,:) = zeros(1,N);...
        %     dR8dl1(N,:) = zeros(1,N); dR8dl2(N,:) = zeros(1,N);dR8dl3(N,:) = zeros(1,N);dR8dl4(N,:) = zeros(1,N);dR8dl5(N,:) = zeros(1,N);dR8dl6(N,:) = zeros(1,N);
        % 
        % dR9dx1(N,:) = zeros(1,N); dR9dx2(N,:) = zeros(1,N);       dR9dy1(N,:) = [zeros(1,N-1),1]; dR9dy2(N,:) = zeros(1,N); dR9dz1(N,:) = zeros(1,N); dR9dz2(N,:) = zeros(1,N);...
        %     dR9dl1(N,:) = zeros(1,N); dR9dl2(N,:) = zeros(1,N);dR9dl3(N,:) = zeros(1,N);dR9dl4(N,:) = zeros(1,N);dR9dl5(N,:) = zeros(1,N);dR9dl6(N,:) = zeros(1,N);
        % 
        % dR10dx1(N,:) = zeros(1,N);       dR10dx2(N,:) = zeros(1,N); dR10dy1(N,:) = zeros(1,N); dR10dy2(N,:) = [zeros(1,N-1),1]; dR10dz1(N,:) = zeros(1,N); dR10dz2(N,:) = zeros(1,N);...
        %     dR10dl1(N,:) = zeros(1,N); dR10dl2(N,:) = zeros(1,N);dR10dl3(N,:) = zeros(1,N);dR10dl4(N,:) = zeros(1,N);dR10dl5(N,:) = zeros(1,N);dR10dl6(N,:) = zeros(1,N);
        % 
        % dR11dx1(N,:) = zeros(1,N); dR11dx2(N,:) = zeros(1,N);       dR11dy1(N,:) = zeros(1,N); dR11dy2(N,:) = zeros(1,N); dR11dz1(N,:) = [zeros(1,N-1),1]; dR11dz2(N,:) = zeros(1,N);...
        %     dR11dl1(N,:) = zeros(1,N); dR11dl2(N,:) = zeros(1,N);dR11dl3(N,:) = zeros(1,N);dR11dl4(N,:) = zeros(1,N);dR11dl5(N,:) = zeros(1,N);dR11dl6(N,:) = zeros(1,N);
        % 
        % dR12dx1(N,:) = zeros(1,N);       dR12dx2(N,:) = zeros(1,N); dR12dy1(N,:) = zeros(1,N); dR12dy2(N,:) = zeros(1,N); dR12dz1(N,:) = zeros(1,N); dR12dz2(N,:) = [zeros(1,N-1),1]; ...
        %     dR12dl1(N,:) = zeros(1,N); dR12dl2(N,:) = zeros(1,N);dR12dl3(N,:) = zeros(1,N);dR12dl4(N,:) = zeros(1,N);dR12dl5(N,:) = zeros(1,N);dR12dl6(N,:) = zeros(1,N);
        % 

                dR6dx1(1,:) = [1,zeros(1,N-1)]; dR6dx2(1,:) = zeros(1,N);       dR6dy1(1,:) = zeros(1,N); dR6dy2(1,:) = zeros(1,N); dR6dz1(1,:) = zeros(1,N); dR6dz2(1,:) = zeros(1,N); ... 
            dR6dl1(1,:) = zeros(1,N); dR6dl2(1,:) = zeros(1,N);dR6dl3(1,:) = zeros(1,N);dR6dl4(1,:) = zeros(1,N);dR6dl5(1,:) = zeros(1,N);dR6dl6(1,:) = zeros(1,N); 

        dR7dx1(1,:) = zeros(1,N);       dR7dx2(1,:) = [1,zeros(1,N-1)]; dR7dy1(1,:) = zeros(1,N); dR7dy2(1,:) = zeros(1,N); dR7dz1(1,:) = zeros(1,N); dR7dz2(1,:) = zeros(1,N);...
            dR7dl1(1,:) = zeros(1,N); dR7dl2(1,:) = zeros(1,N);dR7dl3(1,:) = zeros(1,N);dR7dl4(1,:) = zeros(1,N);dR7dl5(1,:) = zeros(1,N);dR7dl6(1,:) = zeros(1,N);

        dR8dx1(1,:) = zeros(1,N); dR8dx2(1,:) = zeros(1,N);       dR8dy1(1,:) = [1,zeros(1,N-1)]; dR8dy2(1,:) = zeros(1,N); dR8dz1(1,:) = zeros(1,N); dR8dz2(1,:) = zeros(1,N);...
            dR8dl1(1,:) = zeros(1,N); dR8dl2(1,:) = zeros(1,N);dR8dl3(1,:) = zeros(1,N);dR8dl4(1,:) = zeros(1,N);dR8dl5(1,:) = zeros(1,N);dR8dl6(1,:) = zeros(1,N);

        dR9dx1(1,:) = zeros(1,N);       dR9dx2(1,:) = zeros(1,N); dR9dy1(1,:) = zeros(1,N); dR9dy2(1,:) = [1,zeros(1,N-1)]; dR9dz1(1,:) = zeros(1,N); dR9dz2(1,:) = zeros(1,N);...
            dR9dl1(1,:) = zeros(1,N); dR9dl2(1,:) = zeros(1,N);dR9dl3(1,:) = zeros(1,N);dR9dl4(1,:) = zeros(1,N);dR9dl5(1,:) = zeros(1,N);dR9dl6(1,:) = zeros(1,N);

        dR10dx1(1,:) = zeros(1,N); dR10dx2(1,:) = zeros(1,N);       dR10dy1(1,:) = zeros(1,N); dR10dy2(1,:) = zeros(1,N); dR10dz1(1,:) = [1,zeros(1,N-1)]; dR10dz2(1,:) = zeros(1,N);...
            dR10dl1(1,:) = zeros(1,N); dR10dl2(1,:) = zeros(1,N);dR10dl3(1,:) = zeros(1,N);dR10dl4(1,:) = zeros(1,N);dR10dl5(1,:) = zeros(1,N);dR10dl6(1,:) = zeros(1,N);

        dR11dx1(1,:) = zeros(1,N);       dR11dx2(1,:) = zeros(1,N); dR11dy1(1,:) = zeros(1,N); dR11dy2(1,:) = zeros(1,N); dR11dz1(1,:) = zeros(1,N); dR11dz2(1,:) = [1,zeros(1,N-1)]; ...
            dR11dl1(1,:) = zeros(1,N); dR11dl2(1,:) = zeros(1,N);dR11dl3(1,:) = zeros(1,N);dR11dl4(1,:) = zeros(1,N);dR11dl5(1,:) = zeros(1,N);dR11dl6(1,:) = zeros(1,N);



        % dR1dx1(N,:) = [zeros(1,N-1),1]; dR1dx2(N,:) = zeros(1,N);       dR1dy1(N,:) = zeros(1,N); dR1dy2(N,:) = zeros(1,N); dR1dz1(N,:) = zeros(1,N); dR1dz2(N,:) = zeros(1,N); ... 
        %     dR1dl1(N,:) = zeros(1,N); dR1dl2(N,:) = zeros(1,N);dR1dl3(N,:) = zeros(1,N);dR1dl4(N,:) = zeros(1,N);dR1dl5(N,:) = zeros(1,N);dR1dl6(N,:) = zeros(1,N); 
        % 
        % dR2dx1(N,:) = zeros(1,N);       dR2dx2(N,:) = [zeros(1,N-1),1]; dR2dy1(N,:) = zeros(1,N); dR2dy2(N,:) = zeros(1,N); dR2dz1(N,:) = zeros(1,N); dR2dz2(N,:) = zeros(1,N);...
        %     dR2dl1(N,:) = zeros(1,N); dR2dl2(N,:) = zeros(1,N);dR2dl3(N,:) = zeros(1,N);dR2dl4(N,:) = zeros(1,N);dR2dl5(N,:) = zeros(1,N);dR2dl6(N,:) = zeros(1,N);
        % 
        % dR3dx1(N,:) = zeros(1,N); dR3dx2(N,:) = zeros(1,N);       dR3dy1(N,:) = [zeros(1,N-1),1]; dR3dy2(N,:) = zeros(1,N); dR3dz1(N,:) = zeros(1,N); dR3dz2(N,:) = zeros(1,N);...
        %     dR3dl1(N,:) = zeros(1,N); dR3dl2(N,:) = zeros(1,N);dR3dl3(N,:) = zeros(1,N);dR3dl4(N,:) = zeros(1,N);dR3dl5(N,:) = zeros(1,N);dR3dl6(N,:) = zeros(1,N);
        % 
        % dR4dx1(N,:) = zeros(1,N);       dR4dx2(N,:) = zeros(1,N); dR4dy1(N,:) = zeros(1,N); dR4dy2(N,:) = [zeros(1,N-1),1]; dR4dz1(N,:) = zeros(1,N); dR4dz2(N,:) = zeros(1,N);...
        %     dR4dl1(N,:) = zeros(1,N); dR4dl2(N,:) = zeros(1,N);dR4dl3(N,:) = zeros(1,N);dR4dl4(N,:) = zeros(1,N);dR4dl5(N,:) = zeros(1,N);dR4dl6(N,:) = zeros(1,N);
        % 
        % dR5dx1(N,:) = zeros(1,N); dR5dx2(N,:) = zeros(1,N);       dR5dy1(N,:) = zeros(1,N); dR5dy2(N,:) = zeros(1,N); dR5dz1(N,:) = [zeros(1,N-1),1]; dR5dz2(N,:) = zeros(1,N);...
        %     dR5dl1(N,:) = zeros(1,N); dR5dl2(N,:) = zeros(1,N);dR5dl3(N,:) = zeros(1,N);dR5dl4(N,:) = zeros(1,N);dR5dl5(N,:) = zeros(1,N);dR5dl6(N,:) = zeros(1,N);
        % 
        % dR6dx1(N,:) = zeros(1,N);       dR6dx2(N,:) = zeros(1,N); dR6dy1(N,:) = zeros(1,N); dR6dy2(N,:) = zeros(1,N); dR6dz1(N,:) = zeros(1,N); dR6dz2(N,:) = [zeros(1,N-1),1]; ...
        %     dR6dl1(N,:) = zeros(1,N); dR6dl2(N,:) = zeros(1,N);dR6dl3(N,:) = zeros(1,N);dR6dl4(N,:) = zeros(1,N);dR6dl5(N,:) = zeros(1,N);dR6dl6(N,:) = zeros(1,N);

       % Alternative
        % dR2dx1(1,:) = [1,zeros(1,N-1)]; dR2dx2(1,:) = zeros(1,N);       dR2dy1(1,:) = zeros(1,N); dR2dy2(1,:) = zeros(1,N); dR2dz1(1,:) = zeros(1,N); dR2dz2(1,:) = zeros(1,N); ... 
        %     dR2dl1(1,:) = zeros(1,N); dR2dl2(1,:) = zeros(1,N);dR2dl3(1,:) = zeros(1,N);dR2dl4(1,:) = zeros(1,N);dR2dl5(1,:) = zeros(1,N);dR2dl6(1,:) = zeros(1,N); 
        % 
        % dR3dx1(1,:) = zeros(1,N);       dR3dx2(1,:) = [1,zeros(1,N-1)]; dR3dy1(1,:) = zeros(1,N); dR3dy2(1,:) = zeros(1,N); dR3dz1(1,:) = zeros(1,N); dR3dz2(1,:) = zeros(1,N);...
        %     dR3dl1(1,:) = zeros(1,N); dR3dl2(1,:) = zeros(1,N);dR3dl3(1,:) = zeros(1,N);dR3dl4(1,:) = zeros(1,N);dR3dl5(1,:) = zeros(1,N);dR3dl6(1,:) = zeros(1,N);
        % 
        % dR4dx1(1,:) = zeros(1,N); dR4dx2(1,:) = zeros(1,N);       dR4dy1(1,:) = [1,zeros(1,N-1)]; dR4dy2(1,:) = zeros(1,N); dR4dz1(1,:) = zeros(1,N); dR4dz2(1,:) = zeros(1,N);...
        %     dR4dl1(1,:) = zeros(1,N); dR4dl2(1,:) = zeros(1,N);dR4dl3(1,:) = zeros(1,N);dR4dl4(1,:) = zeros(1,N);dR4dl5(1,:) = zeros(1,N);dR4dl6(1,:) = zeros(1,N);
        % 
        % dR5dx1(1,:) = zeros(1,N);       dR5dx2(1,:) = zeros(1,N); dR5dy1(1,:) = zeros(1,N); dR5dy2(1,:) = [1,zeros(1,N-1)]; dR5dz1(1,:) = zeros(1,N); dR5dz2(1,:) = zeros(1,N);...
        %     dR5dl1(1,:) = zeros(1,N); dR5dl2(1,:) = zeros(1,N);dR5dl3(1,:) = zeros(1,N);dR5dl4(1,:) = zeros(1,N);dR5dl5(1,:) = zeros(1,N);dR5dl6(1,:) = zeros(1,N);
        % 
        % dR6dx1(1,:) = zeros(1,N); dR6dx2(1,:) = zeros(1,N);       dR6dy1(1,:) = zeros(1,N); dR6dy2(1,:) = zeros(1,N); dR6dz1(1,:) = [1,zeros(1,N-1)]; dR6dz2(1,:) = zeros(1,N);...
        %     dR6dl1(1,:) = zeros(1,N); dR6dl2(1,:) = zeros(1,N);dR6dl3(1,:) = zeros(1,N);dR6dl4(1,:) = zeros(1,N);dR6dl5(1,:) = zeros(1,N);dR6dl6(1,:) = zeros(1,N);
        % 
        % dR7dx1(1,:) = zeros(1,N);       dR7dx2(1,:) = zeros(1,N); dR7dy1(1,:) = zeros(1,N); dR7dy2(1,:) = zeros(1,N); dR7dz1(1,:) = zeros(1,N); dR7dz2(1,:) = [1,zeros(1,N-1)]; ...
        %     dR7dl1(1,:) = zeros(1,N); dR7dl2(1,:) = zeros(1,N);dR7dl3(1,:) = zeros(1,N);dR7dl4(1,:) = zeros(1,N);dR7dl5(1,:) = zeros(1,N);dR7dl6(1,:) = zeros(1,N);


     end
    J_local(:,:,k) = [dR1dx1, dR1dx2, dR1dy1, dR1dy2 dR1dz1, dR1dz2, dR1dl1, dR1dl2, dR1dl3, dR1dl4, dR1dl5, dR1dl6; ...
                      dR2dx1, dR2dx2, dR2dy1, dR2dy2 dR2dz1, dR2dz2, dR2dl1, dR2dl2, dR2dl3, dR2dl4, dR2dl5, dR2dl6; ...
                      dR3dx1, dR3dx2, dR3dy1, dR3dy2 dR3dz1, dR3dz2, dR3dl1, dR3dl2, dR3dl3, dR3dl4, dR3dl5, dR3dl6;...
                      dR4dx1, dR4dx2, dR4dy1, dR4dy2 dR4dz1, dR4dz2, dR4dl1, dR4dl2, dR4dl3, dR4dl4, dR4dl5, dR4dl6; ...
                      dR5dx1, dR5dx2, dR5dy1, dR5dy2 dR5dz1, dR5dz2, dR5dl1, dR5dl2, dR5dl3, dR5dl4, dR5dl5, dR5dl6; ...
                      dR6dx1, dR6dx2, dR6dy1, dR6dy2 dR6dz1, dR6dz2, dR6dl1, dR6dl2, dR6dl3, dR6dl4, dR6dl5, dR6dl6; ...
                      dR7dx1, dR7dx2, dR7dy1, dR7dy2 dR7dz1, dR7dz2, dR7dl1, dR7dl2, dR7dl3, dR7dl4, dR7dl5, dR7dl6; ...
                      dR8dx1, dR8dx2, dR8dy1, dR8dy2 dR8dz1, dR8dz2, dR8dl1, dR8dl2, dR8dl3, dR8dl4, dR8dl5, dR8dl6; ...
                      dR9dx1, dR9dx2, dR9dy1, dR9dy2 dR9dz1, dR9dz2, dR9dl1, dR9dl2, dR9dl3, dR9dl4, dR9dl5, dR9dl6; ...
                      dR10dx1, dR10dx2, dR10dy1, dR10dy2 dR10dz1, dR10dz2, dR10dl1, dR10dl2, dR10dl3, dR10dl4, dR10dl5, dR10dl6; ...
                      dR11dx1, dR11dx2, dR11dy1, dR11dy2 dR11dz1, dR11dz2, dR11dl1, dR11dl2, dR11dl3, dR11dl4, dR11dl5, dR11dl6; ...
                      dR12dx1, dR12dx2, dR12dy1, dR12dy2 dR12dz1, dR12dz2, dR12dl1, dR12dl2, dR12dl3, dR12dl4, dR12dl5, dR12dl6];
    
    % Indices of states and segments
    Elb(k) = (k-1) * LJ_local+1;        %left boundary of each segment (Element left boundary)
    Erb(k) = Elb(k) + LJ_local- 1;      %right boundary of each segment
    for ii = 1 : Neq
        SLb(ii,:,k) = Elb(k)+ (ii-1)*N; %left boundary of each state ii for each element k (State left boundary)
        SRb(ii,:,k) = Elb(k) + ii*N-1;
        slb(ii) = (ii-1)*N+1;           %Left boundary of each state ii for the "local matrix"
        srb(ii) = ii*N;                 %Right boundary of each state ii for the "local matrix"
    end
    
    %Global Jacobian matrix
    J(Elb(k):Erb(k),Elb(k):Erb(k)) = J_local(:,:,k); %Global matrix

    % Continuity Conditions of the Jacobian Matrix
    if k > 1
        for ii = 1 : Neq
            J(SLb(ii,:,k):SLb(ii,:,k)+ns,:)                           = zeros(1+ns,length(J));
            J(SLb(ii,:,k):SLb(ii,:,k)+ns,SLb(ii,:,k-1):SRb(ii,:,k-1)) = [zeros(1,N-1),-1];
            J(SLb(ii,:,k):SLb(ii,:,k)+ns,SLb(ii,:,k):SRb(ii,:,k)) = [1,zeros(1,N-1)];
        end
        %     J(SLb(ii,:,k):SLb(ii,:,k)+ns,SLb(ii,:,k-1):SRb(ii,:,k-1)) = [zeros(1,N-1),-1;zeros(1,N-1),-1];
        %     J(SLb(ii,:,k):SLb(ii,:,k)+ns,SLb(ii,:,k):SRb(ii,:,k)) = [1,zeros(1,N-1);1,zeros(1,N-1)];
        % end
    end

        if k > 1 
        R(1,k)     = x1(1,k) - x1(end,k-1);
        R(N+1,k)   = x2(1,k) - x2(end,k-1);
        R(2*N+1,k) = y1(1,k) - y1(end,k-1);
        R(3*N+1,k) = y2(1,k) - y2(end,k-1);
        R(4*N+1,k) = z1(1,k) - z1(end,k-1);
        R(5*N+1,k) = z2(1,k) - z2(end,k-1);
        R(6*N+1,k) = L1(1,k) - L1(end,k-1);
        R(7*N+1,k) = L2(1,k) - L2(end,k-1);
        R(8*N+1,k) = L3(1,k) - L3(end,k-1);
        R(9*N+1,k) = L4(1,k) - L4(end,k-1);
        R(10*N+1,k) = L5(1,k) - L5(end,k-1);
        R(11*N+1,k) = L6(1,k) - L6(end,k-1);
        end
end
%% Boundary conditions
switch BCtype
    case "fixed" %fixed final state
 
        R(1,1)   = x1(1,1) - x10;
        R(N+1,1) = x2(1,1) - x20;
        R(6*N+1,1)   = x1(end,end) - x1f;
        R(7*N+1,1)   = x2(end,end) - x2f;
        
        R(2*N+1,1)   = y1(1,1) - y10;
        R(3*N+1,1) = y2(1,1) - y20;
        R(8*N+1,1)   = y1(end,end) - y1f;
        R(9*N+1,1)   = y2(end,end) - y2f;

        R(4*N+1,1)   = z1(1,1) - z10;
        R(5*N+1,1) = z2(1,1) - z20;
        R(10*N+1,1)   = z1(end,end) - z1f;
        R(11*N+1,1)   = z2(end,end) - z2f;

        % R(N+1,1)   = x1(1,1) - x10;
        % R(2*N+1,1) = x2(1,1) - x20;
        % R(2*N,end)   = x1(end,end) - x1f;
        % R(3*N,end)   = x2(end,end) - x2f;
        % 
        % R(4*N+1,1)   = y1(1,1) - y10;
        % R(5*N+1,1) = y2(1,1) - y20;
        % R(5*N,end)   = y1(end,end) - y1f;
        % R(6*N,end)   = y2(end,end) - y2f;
        % 
        % R(7*N+1,1)   = z1(1,1) - z10;
        % R(8*N+1,1) = z2(1,1) - z20;
        % R(8*N,end)   = z1(end,end) - z1f;
        % R(9*N,end)   = z2(end,end) - z2f;
%         R(1,1)   = x1(1,1) - x10;
%         R(N+1,1) = x2(1,1) - x20;
%         R(N,end)   = x1(end,end) - x1f;
%         R(2*N,end)   = x2(end,end) - x2f;
% 
%         J(SLb(1,:,1),SLb(1,:,1):SRb(1,:,1)) = [1,zeros(1,N-1)];
%         J(SLb(2,:,1),SLb(2,:,1):SRb(2,:,1)) = [1,zeros(1,N-1)];
%         
%         
%       
%       J(SRb(1,:,M),SLb(1,:,M):SRb(1,:,M)) = [zeros(1,N-1),1];
%         
%         
%         J(SRb(2,:,M),SLb(2,:,M):SRb(2,:,M)) = [zeros(1,N-1),1];
%         
% % 
%         

    case "free" %free final state
end

%% Output Residual Vector
% Function output
 R = R(:);  %[all state segment1; all states segment2; ... ; all states last segment];

if ip == 1
    R = 1/2*(R'*R);
end


