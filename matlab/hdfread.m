
warning('off','all')

shouldPause = false;
holdOn1dConc = true;

[filen,filep] = uigetfile(['../Data/*.h5'], 'Select file');
filename = strcat(filep,filen);

WI = [0 0 3]; %[45   30];

hinfo = hdf5info(filename);
Names = hinfo.GroupHierarchy.Datasets;
Name = natsortfiles({Names.Name});
H_names = Name(~cellfun('isempty',strfind(Name,'H_')));
Ex_names = Name(~cellfun('isempty',strfind(Name,'Ex')));
Ez_names = Name(~cellfun('isempty',strfind(Name,'Ez')));
AbsE2_names = Name(~cellfun('isempty',strfind(Name,'|E|2')));
N_names = Name(~cellfun('isempty',strfind(Name,'N')));
T_names = Name(~cellfun('isempty',strfind(Name,'Time')));
EpsX_names = Name(~cellfun('isempty',strfind(Name,'EpsX')));
EpsZ_names = Name(~cellfun('isempty',strfind(Name,'EpsZ')));
PertName = Name(~cellfun('isempty',strfind(Name,'Pert')));
Heat_names = Name(~cellfun('isempty',strfind(Name,'Heat')));
Fill_names = Name(~cellfun('isempty',strfind(Name,'Fill')));
Rat_names = Name(~cellfun('isempty',strfind(Name,'Sz')));

OptName = Name(~cellfun('isempty',strfind(Name,'Opt')));
OptMat = h5read(filename,char(OptName));
Nx = round(OptMat(1));
Nz = round(OptMat(2));
Zmin_Wavelength = OptMat(6);
Zmax_Wavelength = OptMat(7);
Xmin_Wavelength = OptMat(8);
Xmax_Wavelength = OptMat(9);

% addpath ../.
% config
% rmpath ../.

Aratio = [12,1];
z = linspace(Zmin_Wavelength,Zmax_Wavelength,Nz);
x = linspace(Xmin_Wavelength,Xmax_Wavelength,Nx);
[X,Z] = meshgrid(x,z);
zlm = [-8 3];
%%
Pert = h5read(filename,char(PertName));
figure(10)
figure(5)
clf;
surf(Z',X',Pert(1:end-1,:)','linestyle','none');
set(gca,'DataAspectRatio',[Aratio max(max(abs(Pert+1e-9)))])
colormap(jet);
view([0 0 10])
colorbar;
ylim([x(1) x(end)]);
xlim([z(1) z(end)]);

if filen(1) == 'D'
    PertType = 'electron density';
else
    PertType = 'medium permittivity';
end
title(['Perturbation of ',PertType]);
%%
maxH = 0;
maxEx = 0;
maxAbsE2 = 0;
 indMin = 2;
%  indMin = numel(H_names);

indStep = 2;
indStepFindMax = 1;
indMax = numel(H_names);

%indStep = round(numel(H_names)/2)

%%


for i = indMin:indStepFindMax:indMax
    %for i = numel(H_names)
    
    AbsE2 = h5read(filename,char(AbsE2_names(i)));
    AbsE22D = reshape(AbsE2(1:end-Nx),[Nx,Nz]);
    if (max(max(abs(AbsE22D)))> maxAbsE2)
        maxAbsE2 = max(max(abs(AbsE22D)));
    end
    
    H = h5read(filename,char(H_names(i)));
    H2D = reshape(H.real + 1i*H.imag,[Nx,Nz]);
    if (max(max(abs(H2D)))> maxH)
        maxH = max(max(abs(H2D)));
    end
        
    Ex = h5read(filename,char(Ex_names(i)));
    Ex2D = reshape(Ex.real(1:end-Nx) + 1i*Ex.imag(1:end-Nx),[Nx,Nz]);
    if (max(max(abs(Ex2D)))> maxEx)
        maxEx = max(max(abs(Ex2D)));
    end
    
end

disp(['----//----', newline,...
    'max Hy = ', num2str(maxH), newline,...
    'max Ex = ', num2str(maxEx), newline,...
    'max |E|^2 = ', num2str(maxAbsE2)]);
S = 0;
% %%
% fig = uifigure('Name','My Figure');
% pnl = uipanel(fig);
% btn = uibutton(pnl); 
%%
for i = indMin:indStep:indMax
    %for i = numel(H_names)
    T = h5read(filename,char(T_names(i)));
    
    AbsE2 = h5read(filename,char(AbsE2_names(i)));
    AbsE22D = reshape(AbsE2(1:end-Nx),[Nx,Nz]);
    
    Heat = h5read(filename,char(Heat_names(i)));
    Heat2D = reshape(Heat(1:end-Nx),[Nx,Nz]);
    
    H = h5read(filename,char(H_names(i)));
    H2D = reshape(H.real + 1i*H.imag,[Nx,Nz]);
    
    Ex = h5read(filename,char(Ex_names(i)));
    Ex2D = reshape(Ex.real(1:end-Nx) + 1i*Ex.imag(1:end-Nx),[Nx,Nz]);
    
    Ez = h5read(filename,char(Ez_names(i)));
    Ez2D = reshape(Ez.real + 1i*Ez.imag,[Nx,Nz]);
    
    N = h5read(filename,char(N_names(i)));
    N2D = reshape(N(1:end-Nx),[Nx,Nz]);
    
    EpsX = h5read(filename,char(EpsX_names(i)));
    EpsX2D = reshape(EpsX.real(1:end-Nx) + 1i*EpsX.imag(1:end-Nx),[Nx,Nz]);
    
    EpsZ = h5read(filename,char(EpsZ_names(i)));
    EpsZ2D = reshape(EpsZ.real(1:end-Nx) + 1i*EpsZ.imag(1:end-Nx),[Nx,Nz]);
    
    F = h5read(filename,char(Fill_names(i)));
    F2D = reshape(F(1:end-Nx),[Nx,Nz]);

    apb = h5read(filename,char(Rat_names(i)));
    apb2D = reshape(apb(1:end-Nx),[Nx,Nz]);
    
    if Nx > 4
        
        figure(1)
        clf
        surf(Z',X',abs(AbsE22D),'linestyle','none');
        colormap(jet)
        colorbar;
        view(WI);
        title(['|E|^2; Max = ',num2str(max(max(abs(AbsE22D)))),'; |E_{th}|^2 = ',num2str(OptMat(5))]);
        %     set(gca,'DataAspectRatio',[Aratio max(max(abs(AbsE22D)))+0.01])
        ylim([x(1) x(end)]);
        xlim([z(1) z(end)]);
        xlim(zlm)
        hold on
        C = max(max(abs(AbsE22D)));
        contour3(Z',X',real(EpsX2D)+C,[1 1]*C,'color',[1 1 1]);
        %     set(gca,'DataAspectRatio',[1 1 1]);
        
        figure(11)
        clf
        surf(Z',X',abs(Heat2D),'linestyle','none');
        colormap(jet)
        colorbar;
        view(WI);
        title(['Heat; Max = ',num2str(max(max(abs(Heat2D))))]);
        %     set(gca,'DataAspectRatio',[Aratio max(max(abs(AbsE22D)))+0.01])
        ylim([x(1) x(end)]);
        xlim([z(1) z(end)]);
        xlim(zlm)
        hold on
        C = max(max(abs(Heat2D)));
        contour3(Z',X',real(EpsX2D)+C,[1 1]*C,'color',[1 1 1]);
        %     set(gca,'DataAspectRatio',[1 1 1]);
        
        figure(12)
        clf
        surf(Z',X',F2D,'linestyle','none');
        colormap(jet)
        colorbar;
        view(WI);
        title(['f; Max = ',num2str(max(max(abs(F2D))))]);
        %     set(gca,'DataAspectRatio',[Aratio max(max(abs(AbsE22D)))+0.01])
        ylim([x(1) x(end)]);
        xlim([z(1) z(end)]);
        xlim(zlm)
        hold on
        C = max(max(abs(F2D)));
        contour3(Z',X',real(EpsX2D)+C,[1 1]*C,'color',[1 1 1]);

        figure(13)
        clf
        surf(Z',X',apb2D,'linestyle','none');
        colormap(jet)
        colorbar;
        view(WI);
        title(['a/b; Max = ',num2str(max(max(abs(apb2D))))]);
        %     set(gca,'DataAspectRatio',[Aratio max(max(abs(AbsE22D)))+0.01])
        ylim([x(1) x(end)]);
        xlim([z(1) z(end)]);
        xlim(zlm)
        hold on
        C = max(max(abs(apb2D)));
        contour3(Z',X',real(EpsX2D)+C,[1 1]*C,'color',[1 1 1]);
        
        %     set(gca,'DataAspectRatio',[1 1 1]);

        figure(2)
        clf
        surf(Z',X',abs(H2D),'linestyle','none');
        colormap(jet)
        colorbar;
        view(WI);
        title(['|H_y| ; Max = ',num2str(max(max(abs(H2D))))]);
        %     set(gca,'DataAspectRatio',[Aratio max(max(abs(H2D)))+0.01])
        ylim([x(1) x(end)]);
        xlim([z(1) z(end)]);
        xlim(zlm)
        hold on
        C = max(max(abs(H2D)));
        contour3(Z',X',real(EpsX2D)+C,[1 1]*C,'color',[1 1 1]);
        
        figure(3)
        clf
        surf(Z',X',abs(Ex2D),'linestyle','none');
        colormap(jet)
        colorbar;
        view(WI);
        title(['|E_x| ; Max = ',num2str(max(max(abs(Ex2D))))]);
        %     set(gca,'DataAspectRatio',[Aratio max(max(abs(Ex2D)))+0.01])
        ylim([x(1) x(end)]);
        xlim([z(1) z(end)]);
        xlim(zlm)
        hold on
        C = max(max(abs(Ex2D)));
        contour3(Z',X',real(EpsX2D)+C,[1 1]*C,'color',[1 1 1]);
        
        figure(4)
        clf
        surf(Z',X',abs(Ez2D),'linestyle','none');
        colormap(jet)
        colorbar;
        view(WI);
        title(['|E_z| ; Max = ',num2str(max(max(abs(Ez2D))))]);
        %     set(gca,'DataAspectRatio',[Aratio max(max(abs(Ez2D)))+0.01])
        ylim([x(1) x(end)]);
        xlim([z(1) z(end)]);
        xlim(zlm)
        hold on
        C = max(max(abs(Ez2D)));
        contour3(Z',X',real(EpsX2D)+C,[1 1]*C,'color',[1 1 1]);
        
        figure(5)
        clf
        surf(Z',X',abs(N2D),'linestyle','none');
        colormap(jet)
        colorbar;
        view(WI);
        title(['N ; Max =  ',num2str(max(max(abs(N2D)))),'; t = ',num2str(T*2.67),' fs']);
        %     set(gca,'DataAspectRatio',[Aratio max(max(abs(AbsE22D)))+0.01])
        ylim([x(1) x(end)]);
        xlim([z(1) z(end)]);
        xlim(zlm)
        hold on
        C = max(max(abs(N2D)));
        contour3(Z',X',real(EpsX2D)+C,[1 1]*C,'color',[1 1 1]);
        %     set(gca,'DataAspectRatio',[1 1 1]);
        
        
        figure(6)
        clf
        subplot(2,1,1);
        surf(Z',X',real(EpsX2D),'linestyle','none');
        colormap(jet)
        colorbar;
        view(WI);
        title(['Re \epsilon_x / \epsilon_s;', ' t = ',num2str(T*2.67),' fs']);
        %     set(gca,'DataAspectRatio',[Aratio max(max(abs(N)))+0.01])
        ylim([x(1) x(end)]);
        xlim([z(1) z(end)]);
        xlim(zlm)
        hold on
        C = 0;
        contour3(Z',X',real(EpsX2D)+C,[1 1]*C,'color',[1 1 1]);
        subplot(2,1,2);
        surf(Z',X',imag(EpsX2D),'linestyle','none');
        colormap(jet)
        colorbar;
        view(WI);
        title(['Im \epsilon_x / \epsilon_s;', ' t = ',num2str(T*2.67),' fs']);
        %     set(gca,'DataAspectRatio',[Aratio max(max(abs(N)))+0.01])
        ylim([x(1) x(end)]);
        xlim([z(1) z(end)]);
        xlim(zlm)
        hold on
        C = 0;
        contour3(Z',X',real(EpsX2D)+C,[1 1]*C,'color',[1 1 1]);
        
        figure(7)
        clf
        subplot(2,1,1);
        surf(Z',X',real(EpsZ2D),'linestyle','none');
        colormap(jet)
        colorbar;
        view(WI);
        title(['Re \epsilon_z / \epsilon_s;', ' t = ',num2str(T*2.67),' fs']);
        %     set(gca,'DataAspectRatio',[Aratio max(max(abs(N)))+0.01])
        ylim([x(1) x(end)]);
        xlim([z(1) z(end)]);
        xlim(zlm)
        hold on
        C = 0;
        contour3(Z',X',real(EpsZ2D)+C,[1 1]*C,'color',[1 1 1]);
        subplot(2,1,2);
        surf(Z',X',imag(EpsZ2D),'linestyle','none');
        colormap(jet)
        colorbar;
        view(WI);
        title(['Im \epsilon_x / \epsilon_s;', ' t = ',num2str(T*2.67),' fs']);
        %     set(gca,'DataAspectRatio',[Aratio max(max(abs(N)))+0.01])
        ylim([x(1) x(end)]);
        xlim([z(1) z(end)]);
        xlim(zlm)
        hold on
        C = 0;
        contour3(Z',X',real(EpsX2D)+C,[1 1]*C,'color',[1 1 1]);

        
        figure(10)
        clf
        Pert = N2D - sum(N2D)/Nx;
        surf(Z',X',Pert,'linestyle','none');
%       set(gca,'DataAspectRatio',[Aratio max(max(abs(Pert+1e-9)))])
        colormap(jet);
        view([0 0 10])
        colorbar;
        ylim([x(1) x(end)]);
        xlim([z(1) z(end)]);
        xlim(zlm)
        hold on
        C = max(max(Pert));
        contour3(Z',X',real(EpsX2D)+C,[1 1]*C,'color',[1 1 1]);
        
    else
        
        indx=find(abs(x - (x(end)-x(1))/2) == min(abs(x - (x(end)-x(1))/2)));
                
        figure(6)
        plot(z,real(EpsX2D(indx,:)),'r',z,imag(EpsX2D(indx,:)),'b')
        %ylim([0 0.35])
        xlim([z(1) z(end)]);
        title('\epsilon_x / \epsilon_0')
        grid on
        
        figure(7)
        plot(z,real(EpsZ2D(indx,:)),'r',z,imag(EpsZ2D(indx,:)),'b')
        %ylim([0 0.35])
        xlim([z(1) z(end)]);
        title('\epsilon_z / \epsilon_0')
        grid on

        figure(4)
        plot(z,real(EpsZ2D(indx,:)),'r',z,imag(EpsZ2D(indx,:)),'b')
        %ylim([0 0.35])
        xlim([z(1) z(end)]);
        title('\epsilon_z / \epsilon_0')
        grid on
        
        figure(2)
        plot(z,real(H2D(indx,:)),'r',z,imag(H2D(indx,:)),'b',z,abs(H2D(indx,:)),'k')
        ylim([-1 ,1]*(maxH+0.05));
        xlim([z(1) z(end)]);
        title('H_y');
        grid on
        
        figure(3)
        plot(z,real(Ex2D(indx,:)),'r',z,imag(Ex2D(indx,:)),'b',z,abs(Ex2D(indx,:)),'k')
        ylim([-1 ,1]*(maxEx+0.05));
        xlim([z(1) z(end)]);
        title('E_x')
        grid on
                
        figure(5)
        plot(z,N2D(indx,:),'k')
        title('N/\epsilon_0 N_c')
        grid on
        if holdOn1dConc
            hold on;
            plot(z,z-z+1.87,'r--');
        end
        
        figure(1)
        plot(z,sqrt(AbsE22D(indx,:)),'r',z([1,end]),[1,1]*sqrt(OptMat(5)),'k')
        title('|E|')
        grid on
        %%
        figure(11)
        plot(z,Heat2D(indx,:),'r')
        title('|Heat|')
        grid on

        figure(12)
        plot(z,F2D(indx,:),'r')
        title('Volume fraction')
        grid on
        figure(13)
        plot(z,apb2D(indx,:),'r')
        title('a / b')
        grid on

    end
    %%
    if shouldPause
        pause();
    end
    drawnow;
end

