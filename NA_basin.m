%% Neighborhood algorithm implemented to be applied to Geological basin inversion, based on: 

% Sambridge, M. (1999). Geophysical inversion with a neighborhood algorithm–Searching a parameters space: Geophysical Journal International, 138, 479–494, doi:10.1046/j.1365-246X.1999.00876.x.
% Sambridge, M. (2001). Finding acceptable models in nonlinear inverse problems using a neighborhood algorithm: Inverse Problems, 17, 387–403, doi:10.1088/0266-5611/17/3/302.
% Sambridge, M., and K. Mosegaard (2002). Monte Carlo methods in geophysical inverse problems: Reviews of Geophysics, 40, 1–29. doi: 10.1029/2000RG000089.

% The code is given with a set of standard test functions('Multi','Holder',
% 'Goldstein-Price','Booth') but can be applied to custom functions
% provided that the user defines:
% - ra : range of each variable, ex. ra{1}=[min_x1, max_x1];
% - fo : fitness function, written as anonymous function
% - tol_cal: the fo value under which the corresponding variable set is considered as calibrated 
% - header : the header of the output file 'population.txt'
% - formato_out : the printing format of the output file 'population.txt'



%% Licence
%  Copyright (c) 2019, Eloisa Salina Borello & Costanzo Peter
% All rights reserved.
% 
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.


function [best_samples,samples,misfit]=NA_basin()
close all
clc
stream = RandStream('twister','Seed',0);
RandStream.setGlobalStream(stream);


%% Test functions
[header, formato_out,ra,fo,tol_cal]= input_test('Multi');
%[header, formato_out,ra,fo,tol_cal]= input_test('Holder');
%[header, formato_out,ra,fo,tol_cal]= input_test('Goldstein-Price');
%[header, formato_out,ra,fo,tol_cal]= input_test('Booth');

%% Parameters
maxiter = 50;
tol_similarity=1e-5;
np=100;% numero suddivisioni per trovare i limiti della griglia di Voronoy. Più è fitto e più la ricerca è raffinata;
alpha0=1/4;
alpha= alpha0;
Ncal=20;

ndim = length(ra);
nS  = round(ndim^3.5);%% From Sambridge 2002
nP = round(nS/2); % numero di campioni per iterazione
nT1 = ceil(nP*alpha);
nT2 = nP-nT1;


%% sampling uniforme iniziale (sampling uniforme n-dimensionale di ns elementi)
misfit=[];
new_samples1=initial_sampling(nS,ndim,ra);
new_samples2=[];
samples=[];


hbar=waitbar(0,['Runing iter 1/',num2str(maxiter)],'WindowStyle','modal');
fout= fopen('population.txt','w');
fprintf(fout,header);
best_fo = zeros(maxiter,1);
ncal_iter= zeros(maxiter,1);
n_calibrated =0;
ind_calibrated=[];
iter=1;

while iter<=maxiter && n_calibrated < Ncal
    
    [samples,new_samples]=avoid_doubles([new_samples1;new_samples2],samples);
    
    [misfit_new]=fo(new_samples);
    misfit=[misfit;misfit_new];
    
    best_fo(iter) = min(misfit);
    
    for inds=1:size(new_samples,1)
        fprintf(fout,formato_out,[num2str(iter),'_iter'],[new_samples(inds,:),misfit_new(inds,:)].');
    end
    
    
    
    
    %% verifica presenza di modelli calibrati
    [ind_calibrated,n_calibrated]=verify_calibration(misfit);
    ncal_iter(iter) = n_calibrated;
    % sampling uniforme di ns/nT elementi per ogni cella di Voronoy selezionata
    
    if iter < maxiter && n_calibrated < Ncal
        % rankin e selezione di nT elementi
        [ ind_sel,ind_rand] = seleziona(misfit,nT1,nT2,ind_calibrated);
        new_samples1 = Gibbs(samples,ind_sel,nT1,nT1,ndim,ra);
        new_samples2 = Gibbs(samples,ind_rand,nT2,nT2,ndim,ra);
    end
    
    waitbar(iter/maxiter,hbar,[' Completed iter ',num2str(iter),'/',num2str(maxiter)]);
    iter=iter+1;
    
end
fclose (fout);
delete(hbar);

best_samples = samples(ind_calibrated,:);
mis_calibrated = misfit(ind_calibrated);


voronoi_evolution('population.txt',ndim,1);

figure_teo(ra,fo,samples,best_samples,misfit,mis_calibrated,'2D');
figure_teo(ra,fo,samples,best_samples,misfit,mis_calibrated,'3D');

figure_convergency(maxiter, best_fo);



function    [ ind_sel,ind_rand] = seleziona(misfit,nT1,nT2,ind_calibrated)
%% escludo da ulteriore ricampionamento i calibrati
[~,ind_sort]= sort(misfit);
ind_sort_ncal = setdiff(ind_sort,ind_calibrated,'stable');
ind_sel = ind_sort_ncal(1:nT1);
ind_bin = ind_sort_ncal(nT1+1:end);
ir = randperm(length(ind_bin),nT2).';
ind_rand = ind_bin(ir);
end

function figure_teo(ra,fo,samples,best_samples,misfit,misfit_cal,mode)

%% visualizzazione grafica teorico
figure;
xx=linspace(ra{1}(1),ra{1}(2),500).';
yy=linspace(ra{2}(1),ra{2}(2),500).';
[XX,YY] = meshgrid(xx,yy);
ZZ= zeros(size(XX));
for i1=1:length(xx)
    for j2=1:length(yy)
        ZZ(j2,i1)= fo([xx(i1),yy(j2)]);
    end
end


colormap('jet');
switch mode
    case '2D'
        
        contour(XX,YY,ZZ,30);
        hold on
        h1=plot(samples(:,1),samples(:,2),'k.');
        h2=plot(best_samples(:,1),best_samples(:,2),'or');
        caxis([min(misfit),max(misfit)]);
        hc=colorbar();
        title(hc,'fo');
        xlabel('x1')
        ylabel('x2')
        title('fo visualization 2D') ;
        legend([h1,h2],{'samples','best samples'})
        
    case '3D'
        
        
        if (log10( max(misfit))-log10( min(misfit)))>2
            
            
            T = real2rgb(log(ZZ), 'jet');
            % Render the surface
            surf(XX, YY, ZZ, T,'linestyle','none');
            set(gca,'Zscale','log','Clim',[min(ZZ(:)) max(ZZ(:))]);
            % Generate the colorbar
            
            
        else
            surf(XX,YY,ZZ,'linestyle','none');
            caxis([min(misfit),max(misfit)]);
            hc=colorbar();
            title(hc,'fo');
        end
        hold on
        h1=scatter3(samples(:,1),samples(:,2),misfit,'.k');
        h2=scatter3(best_samples(:,1),best_samples(:,2),misfit_cal,'or');
        h1.MarkerEdgeColor='k';
        h2.MarkerEdgeColor='r';
        xlabel('x1');
        ylabel('x2');
        zlabel('fo');
        title('fo visualization 3D') ;
        legend([h1,h2],{'samples','best samples'});
        
        
end
 pause(0.2);
end



function figure_convergency(niter,pop_best)
fig=figure;
plot([1:niter].',pop_best,'-');
xlabel('iter');
ylabel('best fo');
title('fo evolution');
end



function new_samples = Gibbs(samples,ind_sel,nS,nT,ndim,ra)
n_new = round(nS/nT);
new_samples = zeros(n_new*nT,ndim);
for ic=1:length(ind_sel) %% voronoi selected cells
    for is = 1:n_new %% new samples
        rand_walk = samples(ind_sel(ic),:);% il cammino inizia dal punto della cella
        for ia= 1: ndim %% samples dimensions
            rs = voronoi_range(samples,ind_sel(ic),ia,ra{ia},rand_walk);% TODO
            if (rs(2)-rs(1)) == 0
                %% zona sovracampionata, mi sposto in mododel tutto random
                rand_walk = initial_sampling(1,ndim,ra);
                break;
            else
                rand_walk(1, ia) = uniform(rs,1);
            end
        end
        new_samples(is+(ic-1)*n_new, :) = rand_walk;
    end
end
end

function rn=uniform(ra,n)
rn= rand(n,1)*(ra(2)-ra(1))+ra(1);
end


function samples=initial_sampling(nS,ndim,ra)
samples = zeros(nS,ndim);
for is = 1:nS %% new samples
    for ia= 1: ndim %% samples dimensions
        samples(is, ia) = uniform(ra{ia},1);
    end
end
end

function [samples_u,new_samples_u]=avoid_doubles(new_samples,samples)
%% TODO: posso anche usare unique con una soglia, ovvero differenze minime non mi interessano
if ~isempty(samples)
    new_samples_u=[];
    for ii=1:size(new_samples,1);
        for is=1: size(samples,1)
            dis=norm(samples(is,:)-new_samples(ii,:));
            if dis <= tol_similarity
                break;
            end
        end
        % fare un ciclo for per calcolare la NORMA!!!!
        if dis >tol_similarity
            new_samples_u=[new_samples_u;new_samples(ii,:)];
        end
    end
else
    new_samples_u=new_samples;
end
samples_u=[samples;new_samples_u];
%bisogna conservare i vecchi samples se no i misfit vecchi non sono
%più corrisponenti
end


function rn=voronoi_range(samples,ind_sel,ia,ra, randwalk)
xmin=ra(1);
xmax=ra(2);
rn = zeros(length(ind_sel),2);


xp = linspace(xmin,xmax,np);
nn= zeros(1,np);
for is=1:length(ind_sel)
    for ip=1:np
        point = randwalk;% parto sempre dall'ultima posizione del random walk
        point(1,ia) = xp(ip);
        dist=zeros(size(samples,1),1);
        for ig = 1:size(samples,1)
            dist(ig) = norm(samples(ig,:) - point);
        end
        [~,nn(ip)]=min(dist);
    end
    
    ips=find(nn == ind_sel(is));
    if ~isempty(ips)
        rn(is,1)=xp(min(ips));
        rn(is,2)=xp(max(ips));
    else
        rn(is,:)=[samples(ind_sel(is),ia),samples(ind_sel(is),ia)]; % già
        % molto campionato in zona, mantengo quella x -> introduco un caso
        % da 0
        %     elseif isequal(modo,'refine')
        %         rn(is,:)=[samples(ind_sel(is),ia)-abs(normrnd(0,0.1)),samples(ind_sel(is),ia)+abs(normrnd(0,0.1))]; % passo random vicino al valore perchè sono prossima al buono
    end
end
end

function  [header, formato_out,ra,fo,tol_cal]= input_test(test)
switch test
    case 'Multi'
        ra{1}=[-2,2];% first variable range 
        ra{2}=[-2,2];% second variable range
        header=['case\tx1\tx2\tmisfit\n'];
        formato_out=['%s\t%.3f\t%.3f\t%.3f\n'];
        fo = @(x) x(:,1).*sin(4*pi.*x(:,1))-x(:,2).*sin(4*pi.*x(:,2)+pi)+1;
        tol_cal=-2; % fo minimum
    case 'Holder' 
        ra{1}=[-10,10];
        ra{2}=[-10,10];
        header=['case\tx1\tx2\tmisfit\n'];
        formato_out=['%s\t%.3f\t%.3f\t%.3f\n'];
        fo = @(x) -abs(sin(x(:,1)).*cos(x(:,1)).*exp(abs(1-(sqrt(x(:,1).^2+x(:,2).^2)./pi))));
        tol_cal=-18;
    case 'Goldstein-Price'
        ra{1}=[-2,2];
        ra{2}=[-2,2];
        header=['case\tx1\tx2\tmisfit\n'];
        formato_out=['%s\t%.3f\t%.3f\t%.3f\n'];
        fo = @(x) (1+(x(:,1)+x(:,2)+1).^2.*(19-14*x(:,1)+3*x(:,1).^2-14*x(:,2)+6.*x(:,1).*x(:,2)+3*x(:,2).^2)).*(30+(2*x(:,1)-3*x(:,2)).^2.*(18-32*x(:,1)+12*x(:,1).^2+48*x(:,2)-36*x(:,1).*x(:,2)+27*x(:,2).^2));
        tol_cal=4;
    case 'Booth'
        ra{1}=[-10,10];
        ra{2}=[-10,10];
        header=['case\tx1\tx2\tmisfit\n'];
        formato_out=['%s\t%.3f\t%.3f\t%.3f\n'];
        fo = @(x) (x(:,1)+2*x(:,2)-7).^2+(2*x(:,1)+x(:,2)-5).^2;
        tol_cal=0.5;
end
end



function [calibrated,n_calibrated]=verify_calibration(misfit)
calibrated=find(misfit <tol_cal);
n_calibrated=length(calibrated);
end


function voronoi_evolution(file_name,ndim,header)%nS si può omettere, e volendo anche niter


    var_name ={};
    fp=fopen(file_name,'r');
    formato=['%s '];
    for ii=1:ndim
        formato = [formato,'%f '];
        var_name{ii} =['x',num2str(ii)];
    end
    formato=[formato,'%f \n'];


    C = textscan(fp, formato,'HeaderLines',header, 'Delimiter', '\t');
    npop=length(C{1});
    it_ind = strfind(C{1}{npop},'_');

    niter =str2double(C{1}{npop}(1:it_ind-1));

    pop_best=zeros(niter,1);
    ibest=zeros(niter,1);
    imax =zeros(niter,1);

    fo_val = C{ndim+2};

    %RICONOSCIMENTO AUTOMATICO ITERAZIONI
    i_iter=zeros(niter,2);
    i_iter(1,1)=1;
    id_iter=1;
    for ii=2:npop
        v=C{1}{ii};
        if(isequal(v(1),num2str(id_iter+1)))|| isequal(v(1:2),num2str(id_iter+1))
            i_iter(id_iter,2)= ii-1;
            i_iter(id_iter+1,1)= ii;
            id_iter=id_iter+1;
        end
    end
    i_iter(id_iter,2)=npop;



    for ii = 1:niter
        %     imin(ii) = (ii-1)*nS +1;
        imax(ii) =  i_iter(ii,2);%(ii-1)*nS +nS;
        [pop_best(ii),ibest(ii)]=min(fo_val(1:imax(ii)));
    end




    %% evoluzione fo vs iter
    %nfig = binomial(n, 2);
    for ix= 1 : ndim
        ra{ix}=[min(C{ix+1}), max(C{ix+1})];
    end

    for ix= 1 : ndim
        for iy = ix+1 : ndim
            resp=questdlg('See the evolution ?','Voronoi Evolution ','Yes','No','Yes');
            if isequal(resp,'Yes')
            iter_step =1;% max(1,round(niter/10));
              figg=figure();
              hold on;
              for it = 1:iter_step:niter
                   
                    set_colors(fo_val,figg);
                    x=C{ix+1}(1:imax(it));
                    y=C{iy+1}(1:imax(it));
                    fo_iter=fo_val(1:imax(it));
                    plot_vor_iter(x,y,fo_iter,ibest(it),ra{ix},ra{iy},var_name{ix},var_name{iy});

                    title(['voronoy evolution - iter ',num2str(it)]) ;
                    pause(0.01);
                end


            end
                figg=figure();
                set_colors(fo_val,figg);
                x=C{ix+1}(1:imax(niter));
                y=C{iy+1}(1:imax(niter));
                fo_iter=fo_val(1:imax(niter));
                plot_vor_iter(x,y,fo_iter,ibest(niter),ra{ix},ra{iy},var_name{ix},var_name{iy});

                title(['voronoy map at final iter']) ;
                 pause(0.2);
            
        end
    end

end

function set_colors(fo_val,fig)
figure(fig);
colormap(jet(256));
caxis([min(fo_val),max(fo_val)]);
hc=colorbar;
set(hc,'Limits',[min(fo_val),max(fo_val)], 'Ticks', linspace( min(fo_val),max(fo_val),10));
title(hc,'fo');
end


function plot_vor_iter(x,y,fo_iter,ibest,rax,ray,var_namex,var_namey)

[v,c] = voronoin([x(:) y(:)]);
hold on
color_unbounded(x,y,c,v,fo_iter);
h1=plot(x,y,'k.');
h2=plot(x(ibest), y(ibest),'*r');
set(gca,'xlim',rax,'ylim',ray);
xlabel(var_namex);
ylabel(var_namey)
legend([h1,h2],{'samples','best sample'});
end

function color_unbounded(x,y,c,v,col)
%% color the unbounded cells
h = voronoi(x,y);
v1 = shiftdim(reshape([h(2).XData;h(2).YData],2,3,[]),2); % Arranged one edge per row, one vertex per slice in the third dimension
nUnbounded = sum(cellfun(@(ic)ismember(1,ic),c));
v1Unbounded = v1(end-(nUnbounded-1):end,:,:);
[~,iBounded] = min(pdist2(v,v1Unbounded(:,:,1))); % Index of the bounded vertex
vUnbounded = v1Unbounded(:,:,2); % Displayed coordinate of the unbounded end of the cell edge

for icv=1:length(c)
    cPatch = c{icv}; % List of vertex indices
    vPatch = v(cPatch,:); % Vertex coordinates which may contain Inf
    idx = find(cPatch==1); % Check if cell has unbounded edges
    if idx
        cPatch = circshift(cPatch,[-idx,0]); % Move the 1 to the end of the list of vertex indices
        vPatch = [vPatch(1:idx-1,:)
            vUnbounded(iBounded == cPatch(end-1),:)
            vUnbounded(iBounded == cPatch(1),:)
            vPatch(idx+1:end,:)]; % Replace Inf values at idx with coordinates from the unbounded edges that meet the two adjacent finite vertices
    end
    if ~isnan(col(icv))
        patch(vPatch(:,1),vPatch(:,2),col(icv));
    else
        patch(vPatch(:,1),vPatch(:,2),[0.6,0.6,0.6]);
    end
end
end


end
