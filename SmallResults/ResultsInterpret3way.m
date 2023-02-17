%% Francesca Middleton 27/09/2022

%%
f1=figure(1);
f1.Position = [10,10,400, 300];
left_color = [0 0 0];
right_color = [0 0 0];
set(f1,'defaultAxesColorOrder',[left_color; right_color]);
load('colorblind_colormap.mat')
clf
index1plot=1;
indexplot2 =length(fns);
yyaxis left
semilogy(fns(index1plot:indexplot2),smse(index1plot:indexplot2),'.','Color', colorblind(1,:), 'MarkerSize',14)
hold on 
semilogy(fns(index1plot:indexplot2),wsmse(index1plot:indexplot2),'o','Color', colorblind(2,:), 'MarkerSize',6)
ylabel('Error (J/mol)')
xlabel('Number of factors')

yyaxis right
semilogy(fns(index1plot:indexplot2),aard(index1plot:indexplot2),'*','Color', colorblind(4,:), 'MarkerSize',5)
hold on
semilogy(fns(index1plot:indexplot2),waard(index1plot:indexplot2),'^','MarkerFaceColor', [.4660 .6740 .1880],'Color',  [.4660 .6740 .1880], 'MarkerSize',5)
hold off
ylabel('AARD (%)')
legend('SMSE','wSMSE', 'AARD', 'wAARD')
figname = strcat(num2str(Temperature),'K-',fillmethod,'-ErrorMetrics');
if small ==1
    figname = strcat('Small-t=',num2str(Thresh),'%-',figname);
end 
saveas(f1,strcat(figname,'.jpg'))
%% Pure he predictions 
clc 
clear
concen = 1;
postprocess = 0;
load('2waySVD-20comps-threshold90par-LOOCV-Xnone-T=298.15-fillmethod=uni-20-Jan-2023.mat')
Temperature =298.15; 
fillmethod = 'uni';
Thresh = '70';
small = 1;

concentrations = conc_interval;
load(filename, 'comps')
mixtures = zeros(size(comps,1)^2,4);
index = 0;
for i = 1:length(comps)
    for j = 1:length(comps)
        index = index+1;
        mixtures(index,:) = [comps(i,:) comps(j,:)];
    end
end 
Xtemp = tril(Xs(:,:,1),-1)+triu(nan(size(Xs(:,:,1))));
Xtemp(Xtemp==0)=nan;
filled_ind = find(~isnan(Xtemp));
Xtemp = nan(size(Xs,[1,2]));
Xtemp(filled_ind)=1;
[row,col] =  find(~isnan(Xtemp));

load(strcat('heUNIFACforT=',num2str(Temperature),'.mat'), 'mixture', 'conc_interval', 'he') %change the temperature here
mixunind = find(ismember(mixture',mixtures(filled_ind,:), 'rows')); % hopefully this extracts what we want 
mixturepred = mixtures(filled_ind,:);
mixunind2 = find(ismember(mixture', mixtures(filled_ind,[3,4,1,2]),'rows'));
mixuniInd = union(mixunind,mixunind2);
mixtureuni = mixture(:,mixuniInd)';
[La,Locb]=ismember(mixtureuni,mixturepred, 'rows');
%mixture
%mixunind = mixunind(Locb,:);
heunifacpred = he(5:5:96, mixuniInd);
heunifacpred(:, length(mixunind):length(mixuniInd)) = flip(heunifacpred(:, length(mixunind):length(mixuniInd)));
%Xm_boot(Xm_boot==0)=nan;
Xm_boot = Xm_boot(~isnan(Xm_boot));
Xm_boot = reshape(Xm_boot, length(concentrations),length(fns),[]);
%Xm_boot2(Xm_boot2==0)=nan;
Xm_boot2 = Xm_boot2(~isnan(Xm_boot2));
Xm_boot2 = reshape(Xm_boot2, length(concentrations),length(fns),[]);
[compoundnames, codesNames] = findnames([mixturepred(:,[1,2]); mixturepred(:,[3,4])]);
[bacnames,bacgroups] = findBac(mixturepred);

whichside = 1:2;
for fnind =4%1:length(fns)% 
    hepred = [Xm_boot(:,fnind,:) Xm_boot2(:,fnind,:)]; 
    hepred = permute(hepred, [1 3 2]);
    Truth = zeros(size(hepred));
    if postprocess ==1
        [heprednew] = postprocesshe(hepred(:,:,whichside),concentrations);
        if whichside == 2
            hepred(:,:,whichside) = heprednew(:,:);
        elseif whichside ==1 
            hepred(:,:,1) = heprednew(:,:,1);% 1 for polynomials, whichside for just removing outliers
        else 
            hepred = heprednew;
        end 
    end 

    for counter = 1:size(row,1)
        Truth(:,counter,1) = X(row(counter),col(counter),:);
        Truth(:,counter,2) = X(col(counter),row(counter),:);
    end 
    Truth(isnan(Truth)) =0;
    indicespred = find(~isnan(hepred(concen,:,1)));
    %indicespred2 = find((hepred(concen,:,1))~=0);
    %indicespred = intersect(indicespred2,indicespred1);
    errorHE(:,:,whichside) = abs(Truth(:,indicespred,whichside)) - abs(hepred(:,indicespred,whichside));
    
    errorHE(isnan(errorHE)) = 0;
    errorHEsystem{fnind} = mean(abs(errorHE),1);
    errorHEcomp{fnind} = mean(errorHE,2);
    errorUNI(:,:) = Truth(:,1:size(heunifacpred,2),1) - heunifacpred; 
    errorUNIsystem = mean(abs(errorUNI),1);
    %hepred = hepred(:,indicespred,:);
    aardHE(:,:,whichside) = abs(errorHE(:,:,whichside))./abs(Truth(:,indicespred,whichside));
    %aardHE(:,:,2) = abs(errorHE(:,:,2))./abs(reshape(Truth(:,indicespred,2),size(errorHE(:,:,2))));
    aardHE(isinf(aardHE))  =1;
    aardHE(isnan(aardHE))=1;
    for c = 1:size(X,3)
        smsec(fnind,c) = sqrt(sum(sum(errorHE(c,:,whichside).^2))/length(indicespred));
        wsmsec(fnind,c) = sqrt(find_wmse_error(errorHE(c,:,whichside),length(indicespred)*2));
        aardc(fnind,c) = sum(sum(abs(aardHE(c,:,whichside))))/length(aardHE(c,:,1))*100;
        waardc(fnind,c) = sqrt(find_wmse_error(aardHE(c,:,whichside),length(indicespred)*2))*100;
        percless(fnind,c,1) = sum(sum((abs(aardHE(c,:,whichside))<0.05)))/length(indicespred)/2*100;
        percless(fnind,c,2) = sum(sum((abs(aardHE(c,:,whichside))<0.15)))/length(indicespred)/2*100;
        percless(fnind,c,3) = sum(sum((abs(aardHE(c,:,whichside))<0.25)))/length(indicespred)/2*100;
        percless(fnind,c,4) = sum(sum((abs(aardHE(c,:,whichside))<0.50)))/length(indicespred)/2*100;
        percless(fnind,c,5) = sum(sum((abs(aardHE(c,:,whichside))<0.75)))/length(indicespred)/2*100;
    end 
    wsmse(fnind) = sqrt(find_wmse_error(errorHE,numel(errorHE)));
    smse(fnind) = sqrt(sum(sum(sum(errorHE.^2)))/numel(errorHE));
    aard(fnind) = mean(aardHE(:))*100;
    waard(fnind) = sqrt(find_wmse_error(aardHE,numel(aardHE)))*100;
end 
%    metrictbl = table(fns',smse',wsmse', aard',waard', 'VariableNames', ["Rank", "SMSE (J/mol)", "wSMSE (J/mol)", "AARD (%)", "wAARD (%)"]);
%% Classification of correct predictions per system 

%fnind =9; %can change this here 
load('colorblind_colormap.mat')
for fnind = fnind

    errorHEsystem1 = errorHEsystem{fnind};
    errorHEsystem1(:,length(filled_ind)+1:end,:)=[];
    errorHEsystem1 = reshape(errorHEsystem1, length(filled_ind),2);
    %errorUNIsystem = reshape(errorUNIsystem, length(filled_ind),2);
    leftcorrect = (abs(errorHEsystem1(:,1))<abs(errorHEsystem1(:,2)));% logical output 
    errorHEbest = errorHEsystem1(:,1).*(leftcorrect) + errorHEsystem1(:,2).*(1-leftcorrect);

    unifaccorrect = (abs(errorHEbest)>abs(errorUNIsystem')); %logical output 
    numleftcorrect = sum(find(leftcorrect));
    numunicorrect = sum(find(unifaccorrect));
    %figure out x and y later 
    leftC = nan(size(X,[1,2]));
    uniC = nan(size(X,[1,2]));
    for i = 1:length(filled_ind)
        leftC(row(i),col(i)) = leftcorrect(i);
        uniC(row(i),col(i)) = unifaccorrect(i);
    end 
    leftC(find(triu(leftC)))=2;
    indices = 1:3:length(compoundnames);

    f8 = figure(8);
    clf
    figname = strcat(num2str(Temperature),'K-',fillmethod,'-',num2str(fnind),'-UppervsLower');
    if small ==1
        figname = strcat('Small-t=',num2str(Thresh),'%-',figname);
    end 
    f8.Name = figname;
    f8.NumberTitle = 'off';
    f8.Position=[10 80 600 450];
    ax = gca;
    ax.XColor='k';
    ax.YColor='k';
    clf
    [x,y] = meshgrid(1:size(X,1), 1:size(X,2));
    pcolor(x, y,leftC)
    xlabel('Compound 1', 'FontSize',14,'Color','k')
    ylabel('Compound 2', 'FontSize',14,'Color','k')
    colormap([colorblind(6,:); colorblind(8,:); colorblind(4,:)])
    c = colorbar('Ticks', [0,1,2], 'TickLabels', ["Lower"; "Upper"; "N/A"],'Color','k', 'FontSize',12); %upper triangular = left, 0 is false
    xticks(indices)
    yticks(indices)
    xticklabels(compoundnames(indices))
    yticklabels(compoundnames(indices))
    saveas(f8,strcat(figname,'.jpg'))

    f9 =figure(9);
    figname = strcat(num2str(Temperature),'K-',fillmethod,'-',num2str(fnind),'-UNIFACvsMCM');
    if small ==1
        figname = strcat('Small-t=',num2str(Thresh),'%-',figname);
    end 
    f9.Name = figname;
    f9.NumberTitle = 'off';
    f9.Position=[700 80 600 450];
    clf
    [x,y] = meshgrid(1:size(X,1), 1:size(X,2));
    uniC(find(triu(uniC)))=2;
    pcolor(x,y,uniC)
    xlabel('Compound 1', 'FontSize',14,'Color','k')
    ylabel('Compound 2', 'FontSize',14,'Color','k')
    colormap([colorblind(1,:); colorblind(5,:); colorblind(4,:)])
    c = colorbar('Ticks', [0,1,2], 'TickLabels', ["MCM"; "UNIFAC"; "N/A"],'Color','k','FontSize',12); % 0 is false
    xticks(indices)
    yticks(indices)
    xticklabels(compoundnames(indices))
    yticklabels(compoundnames(indices))
    saveas(f9,strcat(figname,'.jpg'))
end
%% Define mixture you wish to predict
f13=figure(13);
f13.Position = [100,80,1000, 600];
load('colorblind_colormap.mat')
R2thresh = 0.91;
func_group = {1,2,3:5, [6,7,11], 8:9, 10, 12:17,18:20, 21, 22};
label_func = {'Alkane', 'Primary alcohol', 'Other alcohol','Cycloalkanes', 'Ketone', 'Alkene', 'Ester', 'Amine','Acid','Aldehyde'};
load('colorblind_colormap.mat')
load(filename, 'conc_original', 'HE_original', 'mixtureT') 
mixture_exp = mixtureT;
concentrations = 0.05:0.05:0.95;
for p = 0:2
    clf
    figname = strcat(num2str(Temperature),'K-',fillmethod,'-',num2str(fnind),'-Mixtures',num2str(p));
    if small ==1
        figname = strcat('Small-t=',num2str(Thresh),'%-',figname);
    end 
    t= tiledlayout(4,4);
    t.XLabel.String = 'Composition of compound 1 (%)';
    t.XLabel.FontSize = 12;
    t.XLabel.Color='k';
    t.YLabel.String = 'Excess enthalpy (J/mol)';
    t.YLabel.FontSize = 12;
    t.YLabel.Color='k';
    nexttile
    clc
    for counter =1:15

        index = counter+p*15;
        mixpredind = filled_ind(index);%change the mixture number here
        mixpred = mixtures(mixpredind,:); % might need to change this depending on the results 

        mixpred = mixturepred(index,:);%match mixpred to mixture_exp 
        [~,indexpred] = ismember(mixpred,mixture_exp,'rows');

        if indexpred==0
            [~,indexpred] = ismember([mixpred(3:4) mixpred(1:2)], mixture_exp,'rows');
            if indexpred == 0
                indexpred = 1;
                disp('Mixture does not exist in experimental data')
                disp(mixpred)
            end 
        end 

        conc_exp = conc_original(:,indexpred);
        heexp = HE_original(:,indexpred);
        conc_exp = conc_exp(~isnan(conc_exp));
        heexp = heexp(~isnan(heexp));
        [conc_exp,indicessort]=sort(conc_exp);
        heexp=heexp(indicessort);
        xvar = [concentrations'; conc_exp];
        %xvar = [ones(size(xvar)) xvar];
        yvar = [Truth(:,index,2); heexp];
        [fit,S] = polyfit(xvar,yvar,3);
        y = polyval(fit,xvar);
        R2=corrcoef(y, yvar);
        R2=R2(1,2);
        disp(R2)
        if R2< R2thresh
            conc_exp = 1-(conc_exp);
            %disp(R2)
        end 
        disp(mixpred)


        nexttile
        plot(conc_exp*100, heexp,'^','Color', colorblind(6,:), 'MarkerSize',4,'MarkerFaceColor',colorblind(6,:))
        hold on
        plot(concentrations*100, (Truth(:,index,2)), 'x','Color',colorblind(7,:), 'MarkerSize',6)
        %predictions made by the model 
        if length(whichside)>1
            heplot1 = (hepred(:,index,1));
            heplot2 = (hepred(:,index,2));
            plot(concentrations*100, flip(heplot1),'.', 'Color',colorblind(8,:), 'MarkerSize',12)      
            hold on 
            plot(concentrations*100, (heplot2),'.', 'Color',colorblind(1,:), 'MarkerSize',12)
        else
            heplot1 = (hepred(:,index, whichside));
            if whichside ==1
                heplot1 = flip(heplot1);
            end 
            plot(concentrations*100, heplot1,'.', 'Color',colorblind(8,:), 'MarkerSize',12)
            hold on 
        end 

        %extract experimental predictions 


        %extract unifac predictions 
        load(strcat('heUNIFACforT=',num2str(Temperature),'.mat'), 'mixture', 'conc_interval', 'he') %change the temperature here
        [~,mixunind] = ismember(mixpred,mixture', 'rows');
        if mixunind ~=0 
            heuniplot = he(:,mixunind);
            plot(conc_interval*100, (heuniplot),'k--')
        else 
            [~, mixunind] = ismember([mixpred(3:4) mixpred(1:2)],mixture','rows');
            if mixunind ~=0
                heuniplot = he(:,mixunind);
                plot(conc_interval*100, (heuniplot),'k--')
            else
                disp(mixpred)
                disp('Not found')
            end 
        end 
        %create the title 
        [~,index1]= ismember(mixpred(:,[1,2]),codesNames,'rows');
        [~,index2]= ismember(mixpred(:,[3,4]),codesNames,'rows');
        compoundname1 = compoundnames{index1};
        compoundname2 = compoundnames{index2};
        title(strcat(compoundname1, " & ", compoundname2))% work on making this words - dictionary of compounds and mixture numeric? 
        %title(strcat(label_func(find(strcmp(func_group, {(mixpred(1,1))}))), ' ', num2str(mixpred(1,2)), '+', label_func(find(strcmp(func_group, {num2str(mixpred(1,3))}))), ' ',num2str(mixpred(1,4))))) 

        hold off
    end 
    delete(nexttile(1))
    if length(whichside)>1 
        legend('Experimental','MCM upper', 'MCM lower','MCM average','UNIFAC (Do)', 'Location', 'southeast', 'FontSize',9)
    else
        legend('Experimental','Predictions',  'UNIFAC (Do)')
    end 
    lgd = legend;
    lgd.Layout.Tile = 1;
    saveas(f13,strcat(figname,'.jpg'))
end 
%%  Errors per type of mixtures 
dim = size(Truth);
type = zeros(1,dim(2));
for counter = 1:dim(2) % for all systems 
    %extract the data needed 
    temp = reshape(Truth(:,counter,whichside),dim(1),[]);
    maxVal = max(temp);
    minVal = min(temp);
    if length(whichside)>1
        maxVal = maxVal(1);
        minVal = minVal(1);
    end 
    if maxVal<0 || minVal<0 
        % passes through the origin or is negative 
        if sign(maxVal)~=sign(minVal) %passes through origin, only one is negative 
            type(counter) = 1;
        else 
            type(counter) = 2;
        end 
    else 
        %positive curve 
        if maxVal > 1000 % adjustable ceilings for this classification
            type(counter) = 3; % very positive 
        elseif maxVal >200
            type(counter) = 4; % moderately positive 
        else 
            type(counter) = 5; % less than 200 at the maximum 
        end 
    end 
end
%errors per type 
for counter = 1:5
    temp = errorHE(:,find(type==counter),:);
    temp2 = aardHE(:,find(type==counter),:);
    temp(isnan(temp))=[];
    temp2(isnan(temp2))=[];
   
    errortype{counter} =errorHE(:,find(type==counter),:);
    smsetype(counter) = sqrt(sum(sum(sum(temp.^2)))/numel(temp));
    wsmsetype(counter) = sqrt(find_wmse_error(temp,prod(size(temp))));
    aardtype(counter) = sum(sum(sum(abs(temp2))))/prod(size(temp2))*100;
    waardtype(counter) =sqrt(find_wmse_error(temp2,numel(temp2)))*100;
    nomix(counter) = length(find(type==counter));
end 
tabletypes = [1:5; nomix;smsetype; wsmsetype; aardtype; waardtype]';
disp(tabletypes)
%% Plots
% Plot BAC interests
f12=figure(12);
f12.Position = [100,80,1000, 300];
figname = strcat(num2str(Temperature),'K-',fillmethod,'-',num2str(fnind),'-MixturesInterest');
if small ==1
    figname = strcat('Small-t=',num2str(Thresh),'%-',figname);
end 
load('colorblind_colormap.mat')
load(filename, 'conc_original', 'HE_original', 'mixtureT', 'HE_data_sparse','comps') 
Xs = HE_data_sparse;
Xtemp = tril(Xs(:,:,1),-1)+triu(nan(size(Xs(:,:,1))));
Xtemp(Xtemp==0)=nan;
%filled_ind = find(~isnan(Xtemp));
mixtures = zeros(size(comps,1)^2,4);
index = 0;
for i = 1:length(comps)
    for j = 1:length(comps)
        index = index+1;
        mixtures(index,:) = [comps(i,:) comps(j,:)];
    end
end 
[bacmixtures,indices] = findBacPlot(mixturepred);

t= tiledlayout(2,4);
t.XLabel.String = 'Composition of compound 1 (%)';
t.XLabel.FontSize = 12;
t.XLabel.Color='k';
t.YLabel.String = 'Excess enthalpy (J/mol)';
t.YLabel.FontSize = 12;
t.YLabel.Color='k';
nexttile
count = 1;
R2thresh = 0.89;
for i =1:9
    index = indices{i};
    for j = (index)'
        mixpred = mixturepred(j,:);
        [~,indexpred] = ismember(mixpred,mixtureT,'rows');%mixtureT is mixture_exp - different to mixtures order
        if indexpred==0
            [~,indexpred] = ismember([mixpred(3:4) mixpred(1:2)], mixtureT,'rows');
            if indexpred == 0
                disp('Mixture does not exist in experimental data')
                disp(mixpred)
            end 
        end 
        if indexpred~=0 
            conc_exp = conc_original(:,indexpred);
            heexp = HE_original(:,indexpred);
            conc_exp = conc_exp(~isnan(conc_exp));
            heexp = heexp(~isnan(heexp));
            [conc_exp,indicessort]=sort(conc_exp);
            heexp=heexp(indicessort);
            xvar = [concentrations'; conc_exp];
            %xvar = [ones(size(xvar)) xvar];
            yvar = [Truth(:,j,2); heexp];
            [fit,S] = polyfit(xvar,yvar,3);
            y = polyval(fit,xvar);
            R2=corrcoef(y, yvar);
            R2=R2(1,2);
            disp(R2)
            if R2< R2thresh
                conc_exp = 1-(conc_exp);
                %disp(R2)
            end 
            disp(mixpred)
            count = count +1;
            ax = nexttile;
            ax.XColor = [0 0 0];
            ax.YColor = [0 0 0];
            plot(conc_exp*100, heexp,'^','Color', colorblind(6,:), 'MarkerSize',4,'MarkerFaceColor',colorblind(6,:))
            hold on
            %plot(concentrations*100, (Truth(:,j,2)), 'x','Color',colorblind(7,:), 'MarkerSize',6)
            %predictions made by the model 
            if length(whichside)>1
                heplot1 = (hepred(:,j,1));
                heplot2 = (hepred(:,j,2));
                heplot3 = mean([flip(heplot1),heplot2],2);
                %disp(heplot3)
                plot(concentrations*100, flip(heplot1),'.', 'Color',colorblind(8,:), 'MarkerSize',12)      
                hold on 
                plot(concentrations*100, (heplot2),'.', 'Color',colorblind(1,:), 'MarkerSize',12)
                plot(concentrations*100, (heplot3),'s', 'Color',colorblind(9,:), 'MarkerSize',5)
            
            else
                heplot1 = (hepred(:,j, whichside));
                if whichside ==1
                    heplot1 = flip(heplot1);
                end 
                plot(concentrations*100, heplot1,'.', 'Color',colorblind(8,:), 'MarkerSize',12)
                hold on 
            end 

            %extract experimental predictions 


            %extract unifac predictions 
            load(strcat('heUNIFACforT=',num2str(Temperature),'.mat'), 'mixture', 'conc_interval', 'he') %change the temperature here
            [~,mixunind] = ismember(mixpred,mixture', 'rows');
            if mixunind ~=0 
                heuniplot = he(:,mixunind);
                plot(conc_interval*100, (heuniplot),'k--')
            else 
                [~, mixunind] = ismember([mixpred(3:4) mixpred(1:2)],mixture','rows');
                if mixunind ~=0
                    heuniplot = he(:,mixunind);
                    plot(conc_interval*100, (heuniplot),'k--')
                else
                    disp(mixpred)
                    disp('Not found')
                end 
            end 
            %create the title 
            [~,index1]= ismember(mixpred(:,[1,2]),codesNames,'rows');
            [~,index2]= ismember(mixpred(:,[3,4]),codesNames,'rows');
            compoundname1 = compoundnames{index1};
            compoundname2 = compoundnames{index2};
            title(strcat(compoundname1, " & ", compoundname2, " (",num2str(i),")"), 'FontSize',10)% work on making this words - dictionary of compounds and mixture numeric? 
            %title(strcat(label_func(find(strcmp(func_group, {(mixpred(1,1))}))), ' ', num2str(mixpred(1,2)), '+', label_func(find(strcmp(func_group, {num2str(mixpred(1,3))}))), ' ',num2str(mixpred(1,4))))) 

            hold off
            %xlabel('Composition of compound 1 (%)', 'FontSize',9)
            %ylabel('Excess enthalpy (J/mol)', 'FontSize',9)
        end 
    end 
end 
delete(nexttile(1))
if length(whichside)>1 
    legend('Experimental','MCM upper', 'MCM lower','MCM average','UNIFAC (Do)', 'Location', 'southeast', 'FontSize',9)
else
    legend('Experimental','Predictions',  'UNIFAC (Do)')
end 
lgd = legend;
lgd.Layout.Tile = 1;
saveas(f12,strcat(figname,'.jpg'))
% Error plot within the 3-way arrays 
f10 = figure(10);
f10.Position = [10 20 800 500];
figname = strcat(num2str(Temperature),'K-',fillmethod,'-',num2str(fnind),'-Err-3way');
if small ==1
    figname = strcat('Small-t=',num2str(Thresh),'%-',figname);
end 
clf
TempIndex = 1;
errorplot = nan(size(X(:,:,:,TempIndex)));
for i =1:length(col)
    errorplot(row(i),col(i),:) = errorHE(:,i,1);
    errorplot(col(i),row(i),:) = errorHE(:,i,2);
end 
[X,Y,Z] = meshgrid(1:size(errorplot,1), 1:size(errorplot,1), 5:5:95);
xslice = 1:1:size(errorplot,1);    % location of y-z planes
yslice = 1:1:size(errorplot,1);     % location of x-z plane
zslice = 5:5:95;         % location of x-y planes
slice(X,Y,Z,errorplot,xslice,yslice,zslice)
%pbaspect([1 1 0.5])
%[caz,cel] = view % find the camera angle 
view(-33, 55) %33 comps -30.1688, 23.6838; 38 comps -26.2693, 23.6777; 56 comps -27.0897, 18.4455
colormap(parula)
c = colorbar;
c.Label.String = 'Error (J/mol)';
c.Label.FontSize = 12;
xlabel('Compound 1')
ylabel('Compound 2')
zlabel('Composition of compound 1 (%)', 'FontSize',12)
indices = 1:3:length(compoundnames);
xticks(indices)
yticks(indices)
xticklabels(compoundnames(indices))
yticklabels(compoundnames(indices))
saveas(f10,strcat(figname,'.jpg'))
% Errors per type of BAC groups in the mixture 
for i = 1:9
    index = find(bacgroups==i);
    if index
        errtemp = errorHE(:,index,whichside);
        ardtemp = aardHE(:,index,whichside);
        SMSE_BAC(i) = sqrt(sum(errtemp(:).^2)/numel(errtemp));
        wSMSE_BAC(i) = sqrt(find_wmse_error(errtemp(:),numel(errtemp)));
        AARD_BAC(i) = mean(ardtemp(:))*100;
        wAARD_BAC(i) = 100*sqrt(find_wmse_error(ardtemp(:),numel(ardtemp)));
        disp('BAC Group')
        disp(i)
        %disp(index)
        comp_type1{i} = (bacnames{index(1),1});
        comp_type2{i} = bacnames{index(1),2};
    end
end 
indexkeep = find((wAARD_BAC)~=0);
tbl_BAC = table(comp_type1(indexkeep)', comp_type2(indexkeep)',SMSE_BAC(indexkeep)',AARD_BAC(indexkeep)',wSMSE_BAC(indexkeep)',wAARD_BAC(indexkeep)', 'VariableNames', ["Compound type 1", "Compound type 2", "SMSE (J/mol)","AARD (%)","wSMSE (J/mol)","wAARD (%)" ]);
label_BAC = {'NA', 'HA','SA','HD'};
[~,ib,~] = intersect(label_BAC,(comp_type1(indexkeep)));
[~,ia,~] = intersect(label_BAC,unique(comp_type2(indexkeep)));
ia = sort(ia);
ib = sort(ib);
% BAC wAARD
f8 = figure(8);
f8.Position = [10 20 450 300];
figname = strcat(num2str(Temperature),'K-',fillmethod,'-',num2str(fnind),'-BAC-wAARD');
if small ==1
    figname = strcat('Small-t=',num2str(Thresh),'%-',figname);
end 
clf
h = heatmap(tbl_BAC, 'Compound type 1', 'Compound type 2','ColorVariable', 'wAARD (%)', 'ColorMethod', 'none');
h.YDisplayData = flip(label_BAC(ia));
h.XDisplayData = label_BAC(ib);
h.Colormap=winter;
h.FontSize = 14;
saveas(f8,strcat(figname,'.jpg'))
% BAC AARD
f9 = figure(9);
f9.Position = [10 20 450 300];
figname = strcat(num2str(Temperature),'K-',fillmethod,'-',num2str(fnind),'-BAC-AARD');
if small ==1
    figname = strcat('Small-t=',num2str(Thresh),'%-',figname);
end 
clf
h = heatmap(tbl_BAC, 'Compound type 1', 'Compound type 2','ColorVariable', 'AARD (%)', 'ColorMethod', 'none');
h.YDisplayData = flip(label_BAC(ia));
h.XDisplayData = label_BAC(ib);
h.Colormap=winter;
h.FontSize = 14;
saveas(f9,strcat(figname,'.jpg'))
% BAC wSMSE 
f10 = figure(10);
f10.Position = [10 20 450 300];
figname = strcat(num2str(Temperature),'K-',fillmethod,'-',num2str(fnind),'-BAC-wSMSE');
if small ==1
    figname = strcat('Small-t=',num2str(Thresh),'%-',figname);
end 
clf
h = heatmap(tbl_BAC, 'Compound type 1', 'Compound type 2','ColorVariable', 'wSMSE (J/mol)', 'ColorMethod', 'none');
h.YDisplayData = flip(label_BAC(ia));
h.XDisplayData = label_BAC(ib);
h.Colormap=winter;
h.FontSize = 14;
saveas(f10,strcat(figname,'.jpg'))
%BAC SMSE
f11 = figure(11);
f11.Position = [10 20 450 300];
figname = strcat(num2str(Temperature),'K-',fillmethod,'-',num2str(fnind),'-BAC-SMSE');
if small ==1
    figname = strcat('Small-t=',num2str(Thresh),'%-',figname);
end 
clf
h = heatmap(tbl_BAC, 'Compound type 1', 'Compound type 2','ColorVariable', 'SMSE (J/mol)', 'ColorMethod', 'none');
h.YDisplayData = flip(label_BAC(ia));
h.XDisplayData = label_BAC(ib);
h.Colormap=winter;
h.FontSize = 14;
saveas(f11,strcat(figname,'.jpg'))
% error distribution 2-way
f6 = figure(6);
f6.Position = [10 20 400 300];
figname = strcat(num2str(Temperature),'K-',fillmethod,'-',num2str(fnind),'-2Hist');
if small ==1
    figname =strcat('Small-t=',num2str(Thresh),'%-',figname);
end 
clf
errorHE(errorHE==0)=nan; 
histogram(errorHE(:), 20, 'FaceColor', colorblind(6,:)) %, 'BinLimits', [-1000 1000]
ylabel('Instances','FontSize',13, 'Color', 'k')
xlabel('MCM error (J/mol)', 'FontSize',13, 'Color', 'k')
saveas(f6,strcat(figname,'.jpg'))
% 3way histogram 
load('colorblind_colormap.mat')
f7 = figure(7);
f7.Position = [10 20 450 300];
figname = strcat(num2str(Temperature),'K-',fillmethod,'-',num2str(fnind),'-3Hist');
if small ==1
    figname = strcat('Small-t=',num2str(Thresh),'%-',figname);
end 
clf
errorHE(errorHE==0)=nan; 
edges = {-1000:100:1000; 0:5:100}; % bin edges
data = [errorHE(:) abs(aardHE(:)).*100]; 
hist3(data, 'Edges', edges,'CdataMode','auto','FaceColor','interp')
zlabel('Instances','FontSize',13, 'Color', 'k')
xlabel('Error (J/mol)', 'FontSize',13, 'Color', 'k')
ylabel('ARD (%)' ,'FontSize', 13, 'Color', 'k')
colormap(colorblind(flip([9,1,8,6,3]),:))
saveas(f7,strcat(figname,'.jpg'))
% parity plot  
f5 = figure(5);
f5.Position=[10 10 500 400];
figname = strcat(num2str(Temperature),'K-',fillmethod,'-',num2str(fnind),'-Parity');
if small ==1
    figname = strcat('Small-t=',num2str(Thresh),'%-',figname);
end 
clf
%choose between plot and loglog 
plotx =Truth(:,:,whichside);
indremove = find(all(plotx(:,:,1)==0,1));
plotx(:,indremove,:)=[];
ploty = hepred(:,:,whichside);
ploty(:,indremove,:)=[];
plot(plotx(:),ploty(:),'.','Color',colorblind(1,:), 'MarkerSize',6)
hold on 
%Import all the UNIFAC predictions 
plotuni = 1;
load(strcat('heUNIFACforT=',num2str(Temperature),'.mat'), 'mixture', 'conc_interval', 'he') %change the temperature here
mixtureT = mixtures(filled_ind,:);
%mixtureT(indremove,:)=[];
for i = 1:size(mixtureT,1)
    mixpred = mixtureT(i,:);
    [~,mixunind] = ismember(mixpred,mixture', 'rows');
    if mixunind ~=0
        heuni(:,i) = he(5:5:96,mixunind); % check these indices 
        indexuni(i)=1;
    else
        [~,mixunind] = ismember([mixpred(3:4) mixpred(1:2)],mixture', 'rows');
        if mixunind ~=0
            heuni(:,i) = he(5:5:96,mixunind); % check these indices 
            indexuni(i)=1;
        else 
            indexuni(i) = 0;
        end 
        
    end 
    
end
if plotuni
    if length(whichside) >1
        plot(plotx(:), [heuni(:); heuni(:)],'.','Color', 'k', 'MarkerSize',6)
    else 
        plot(plotx(:), [heuni(:)],'.','Color', 'k', 'MarkerSize',6)
    end    
end 
%prediction interval
PIx = min(plotx(:)):100:max(plotx(:));
syx = std(plotx(:));
hi = 1/numel(plotx(:)) + (PIx-mean(PIx)).^2/sum((PIx-mean(PIx)).^2);
PIy1 = PIx - 1.96.*syx.*sqrt(1+hi);
PIy2 = PIx + 1.96.*syx.*sqrt(1+hi);
plot(PIx,PIy1, '--', 'Color', colorblind(7,:), 'LineWidth', 1) %make this pink 
plot(PIx,PIy2, '--', 'Color', colorblind(7,:), 'LineWidth', 1) %make this pink 
%y=x line 
plot([min(plotx(:)),max(plotx(:))], [min(plotx(:)),max(plotx(:))],'-', 'Color',colorblind(4,:), 'LineWidth', 1.2)
hold off 
xlabel('Experimental (J/mol)','FontSize',13, 'Color', 'k')
ylabel('Prediction (J/mol)','FontSize',13, 'Color', 'k')
if plotuni
    legend('MCM', 'UNIFAC','PI (95%)', 'FontSize',8,'Location', 'northwest');
else 
    legend('MCM','PI (95%)', 'FontSize',8,'Location', 'northwest');
end 
saveas(f5,strcat(figname,'.jpg'))
% Functional group errors
load(filename, 'comps', 'mixtureT')
load(filenamesave)
mixtureT=mixtures(filled_ind,:);
func_group = {1,2,3:5, [6,7,11], 8:9, 10, 12:17,18:20, 21, 22};
label_func = {'Alkane', 'Primary alcohol', 'Other alcohol','Cycloalkanes', 'Ketone', 'Alkene', 'Ester', 'Amine','Acid','Aldehyde'};
err_funcgroup = cell(length(func_group),length(func_group));
ard_funcgroup = cell(length(func_group),length(func_group));
smse_funcgroup = zeros(length(func_group),length(func_group));
waard_funcgroup = zeros(length(func_group),length(func_group));
wsmse_funcgroup = zeros(length(func_group),length(func_group));
AARD_funcgroup = zeros(length(func_group),length(func_group));
no_mix = zeros(length(func_group),length(func_group));
for func = 1:length(func_group)
    for func2 = 1:length(func_group)
        %extract arrays of groups 
        funcgroup1 = func_group{func};
        funcgroup2 = func_group{func2};
        %find mixtures with this combination
        indexloop = 1;
        index=[];
        for i =1:length(funcgroup1)
            indextemp = find(mixtureT(:,1)==funcgroup1(i));
            if ~isempty(indextemp) %&& index~=0
                index(indexloop:length(indextemp)+indexloop-1,1) = indextemp;
                indexloop = length(index)+indexloop;
            end
        end 
        indexloop = 1;
        index2=[];
        for i =1:length(funcgroup2)
            indextemp = find(mixtureT(:,3)==funcgroup2(i));
            if ~isempty(indextemp) %&& index2~=0
                index2(indexloop:length(indextemp)+indexloop-1,1) = indextemp;
                indexloop = length(index2)+indexloop;
            end
        end 

        indices = intersect(index,index2);
        indices(indices==0) =[];
        %export errors and AARDs
        if ~isempty(indices)
            err_funcgroup{func,func2} = errorHE(:,indices,whichside); 
            errtemp = errorHE(:,indices,whichside);
            smse_funcgroup(func,func2) = sqrt(sum(errtemp(:).^2)/length(errtemp(:)));
            wsmse_funcgroup(func,func2) = sqrt(find_wmse_error(errtemp(:),length(errtemp(:))));
            ardtemp = aardHE(:,indices,whichside);
            waard_funcgroup(func,func2) = sqrt(find_wmse_error(ardtemp(:),length(ardtemp(:))));
            ard_funcgroup{func,func2} = aardHE(:,indices,whichside);
            AARD_funcgroup(func,func2) = mean(ardtemp(:));  
            no_mix(func,func2) = length(ardtemp(:)/9);
        end 
    end 
end 
% Heatmap of the errors per functional group 
SMSEplot = smse_funcgroup;
AARDplot=AARD_funcgroup.*100;
triuSMSE = triu(SMSEplot)';
trilSMSE = tril(SMSEplot);
tempvar = (trilSMSE+triuSMSE);
tempvar(tempvar==0)=nan;
%create a table 
[row,col] = find(~isnan(tempvar));
clear smse_tbl aard_tbl wsmsetbl funcgroup1tbl funcgroup2tbl
for counter = 1:length(row)
    funcgroup1tbl{counter} = label_func{row(counter)};
    funcgroup2tbl{counter} = label_func{col(counter)};
    if ~(any(isnan([SMSEplot(row(counter), col(counter)),  SMSEplot(col(counter),row(counter))])))
        smse_tbl(counter) = mean([SMSEplot(row(counter), col(counter)),SMSEplot(col(counter),row(counter))], 'all');
        aard_tbl(counter) = mean([AARDplot(row(counter),col(counter)),AARDplot(col(counter),row(counter))], 'all');
        wsmsetbl(counter) = mean([wsmse_funcgroup(row(counter),col(counter)) , wsmse_funcgroup(col(counter),row(counter))], 'all');
        waardtbl(counter) = mean([waard_funcgroup(row(counter),col(counter)) , waard_funcgroup(col(counter),row(counter))], 'all')*100;
    elseif  isnan(SMSEplot(row(counter), col(counter)))
        smse_tbl(counter) = SMSEplot(col(counter),row(counter));
        aard_tbl(counter) = AARDplot(col(counter),row(counter));
        wsmsetbl(counter) = wsmse_funcgroup(col(counter),row(counter));
        waardtbl(counter) = waard_funcgroup(col(counter),row(counter))*100;
    else 
        smse_tbl(counter) = SMSEplot(row(counter),col(counter));
        aard_tbl(counter) = AARDplot(row(counter),col(counter));
        wsmsetbl(counter) = wsmse_funcgroup(row(counter),col(counter));
        waardtbl(counter) = waard_funcgroup(row(counter),col(counter))*100; 
    end 
end 
% heatmap 
tbl = table(funcgroup1tbl', funcgroup2tbl',smse_tbl',aard_tbl',wsmsetbl',waardtbl', 'VariableNames', ["Functional group 1", "Functional group 2", "SMSE (J/mol)","AARD (%)","wSMSE (J/mol)","wAARD (%)" ]);
% heatmap 
f1 = figure(1);
f1.Position = [50, 50,700, 360];  %larger figuures 700, 360  small figures  500, 250
figname = strcat(num2str(Temperature),'K-',fillmethod,'-',num2str(fnind),'-Func-SMSE');
if small ==1
    figname = strcat('Small-t=',num2str(Thresh),'%-',figname);
end 
h = heatmap(tbl, 'Functional group 1', 'Functional group 2','ColorVariable', 'SMSE (J/mol)', 'ColorMethod', 'none');
[~,ia,~] = intersect(label_func,unique(funcgroup1tbl));
[~,ib,~] = intersect(label_func,unique(funcgroup2tbl));
ia = sort(ia);
ib = sort(ib);
h.YDisplayData = flip(label_func(ib));
 h.XDisplayData = label_func(ia);
h.Colormap=winter;
h.FontSize = 14;
saveas(f1,strcat(figname,'.jpg'))
% heatmap 
f2 = figure(2);
f2.Position = [50, 50, 700, 360]; %larger figuures 700, 360 
figname = strcat(num2str(Temperature),'K-',fillmethod,'-',num2str(fnind),'-Func-AARD');
if small ==1
    figname =strcat('Small-t=',num2str(Thresh),'%-',figname);
end 
clf
h = heatmap(tbl, 'Functional group 1', 'Functional group 2','ColorVariable', 'AARD (%)', 'ColorMethod', 'none');
h.YDisplayData = flip(label_func(ib));
 h.XDisplayData = label_func(ia);
h.Colormap=winter;
h.FontSize = 14;
saveas(f2,strcat(figname,'.jpg'))
% heatmap 
f3= figure(3);
f3.Position = [50, 50, 700, 360]; %larger figuures 700 360 
figname = strcat(num2str(Temperature),'K-',fillmethod,'-',num2str(fnind),'-Func-wSMSE');
if small ==1
    figname = strcat('Small-t=',num2str(Thresh),'%-',figname);
end 
clf
h = heatmap(tbl, 'Functional group 1', 'Functional group 2','ColorVariable', 'wSMSE (J/mol)', 'ColorMethod', 'none');
h.YDisplayData = flip(label_func(ib));
 h.XDisplayData = label_func(ia);
h.Colormap=winter;
h.FontSize = 14;
saveas(f3,strcat(figname,'.jpg'))
% heatmap 
f4= figure(4);
f4.Position = [50, 50, 700, 360]; %larger figuures 700 360 
figname = strcat(num2str(Temperature),'K-',fillmethod,'-',num2str(fnind),'-Func-wAARD');
if small ==1
    figname = strcat('Small-t=',num2str(Thresh),'%-',figname);
end 
clf
h = heatmap(tbl, 'Functional group 1', 'Functional group 2','ColorVariable', 'wAARD (%)', 'ColorMethod', 'none');
h.YDisplayData = flip(label_func(ib));
 h.XDisplayData = label_func(ia);
h.Colormap=winter;
h.FontSize = 14;
saveas(f4,strcat(figname,'.jpg'))

%% functions 
function [compoundnames,codesout] = findnames(codesin)
    compounds = {'Methane','Ethane', 'Propane', 'Butane', 'Pentane', 'Hexane', 'Heptane', 'Octane', 'Nonane', 'Decane', 'Dodecane', 'Methanol', 'Ethanol', '1-Propanol', '1-Butanol', '1-Pentanol', '1-Hexanol', '1-Heptanol', '1-Octanol', '1-Nonanol', '1-Decanol', '2-Propanol', '2-Butanol', '2-Pentanol', '2-Hexanol','2-Octanol', 'Isobutanol', 'Tertbutylalcohol', 'Cyclopentane', 'Cyclohexane','Cycloheptane', 'Cyclooctane','Benzene', 'Toluene', 'Ethanal', 'Propanal', 'Butanal', 'Pentanal', 'Hexanal', 'Ethene', '1-Pentene', '1-Hexene', '1-Heptene', '1-Octene', '1-Nonene', '1-Decene', 'Acetone', 'Butanone', '2-Propanone', '2-Hexanone', '2-Heptanone', '2-Octanone','Formic acid', 'Acetic acid', 'Propionic acid', 'Butyric acid', 'Pentanoic acid', 'Hexanoic acid', 'Heptanoic acid', 'Octanoic acid', 'Nonanoic acid', 'Decanoic acid', '3-Pentanone', '3-Hexanone', '3-Heptanone', '3-Octanone', 'Methyl formate', 'Methyl acetate', 'Methyl propionate', 'Methyl butyrate', 'Methyl pentanoate', 'Methyl hexanoate', 'Methyl benzoate', 'Ethyl benzoate', 'Ethyl formate', 'Ethyl acetate', 'Ethyl propionate', 'Ethyl butyrate', 'Ethyl pentanoate', 'Ethyl hexanoate', 'Ethyl heptanoate', 'Ethyl octanoate' , 'Propyl formate','Propyl acetate', 'Propyl propionate', 'Propyl butyrate', 'Butyl formate', 'Butyl acetate', 'Butyl butyrate', 'Pentyl acetate', 'Propylamine', 'Butylamine', 'Pentyalamine', 'Hexylamine', 'Heptylamine', 'Octylamine', 'Nonylamine', 'Decylamine', 'Aniline', 'Benzylamine'};
    codes = [ones(11,1), ([1:10,12])'; 2*ones(10,1), (1:10)'; 3*ones(5,1), ([3:6,8])'; 4, 4; 5, 4;11*ones(4,1), (5:8)'; 6, 0; 7, 0; 22*ones(5,1), (2:6)'; 10*ones(7,1), ([2,5:10])'; 8*ones(6,1), (3:8)'; 21*ones(10,1), (1:10)'; 9*ones(4,1), (5:8)';12*ones(6,1), (1:6)'; 17*ones(2,1), (1:2)'; 13*ones(8,1), (1:8)'; 14*ones(4,1), (1:4)'; 15*ones(3,1), ([1,2,4])'; 16, 2; 18*ones(8,1), (3:10)'; 19, 0; 20, 0];
    indices = find(ismember(codes, codesin,'rows'));
    for i = 1:length(indices)
        compoundnames{i} = compounds{indices(i)};
        codesout(i,:) =codes(indices(i),:);
    end 
end 
function [bacmixtures,indices] = findBacPlot(mixtures)
    bacmixtures = cell(9,1);  
    indices = cell(9,1);
    bac{1} = [ones(6,1), [ones(5,1); 2], ones(6,1), [2,3,4,5,7,7]'; 6*ones(3,1), zeros(3,1), [11, 7, 1]', [6,0,7]'; 11,5,11,8];
    bac{2} = [8*ones(6,1), [3*ones(4,1);4;4], [1;6;11;1;6;1], [5;0;6;6;0;7]; 13, 2, 6, 0];
    bac{4} = [8*ones(3,1), [3;3;4], [8;13;13], [4;2;2]];
    bac{5} = [2*ones(9,1), [ones(4,1);2;2;3;3;3], [1;1;6;1;6;1;1;6;1], [2;3;0;6;0;7;2;0;6] ; 3,4,11,6];
    bac{8} = [2*ones(6,1), [1;1;1;2;2;2], [8; 13; 13;8;8;13],[3;2;4;3;4;2]; 2,1,12,2;2,1,15,2];
    bac{8}(1,:)=[];
    bac{9} = [2*ones(7,1), [ones(3,1);2;1;2;3], [2;2;2;21;12;12;21], [2;3;4;2;1;2;2]];
    for i = 1:9
        bacloop = bac{i};
        if length(bacloop)>1
            index = find(ismember(mixtures,bacloop,'rows'));
            index2 = find(ismember(mixtures, [bacloop(:,3:4), bacloop(:,1:2)],'rows'));
            index = [index;index2];
            bacmixtures{i} = mixtures(index,:);
            indices{i} = index;
        end 
         
    end 
end 
function [bacnames,bacgroups] = findBac(mixtures)
    bacNamegroups = {'NA', 'HA', 'SA', 'HD'};
    %codes = [ones(11,1), ([1:10,12])'; 2*ones(10,1), (1:10)'; 3*ones(5,1), ([3:6,8])'; 4, 4; 5, 4;11*ones(4,1), (5:8)'; 6, 0; 7, 0; 22*ones(5,1), (2:6)'; 10*ones(7,1), ([2,5:10])'; 8*ones(6,1), (3:8)'; 21*ones(10,1), (1:10)'; 9*ones(4,1), (5:8)';12*ones(6,1), (1:6)'; 17*ones(2,1), (1:2)'; 13*ones(8,1), (1:8)'; 14*ones(4,1), (1:4)'; 15*ones(3,1), ([1,2,4])'; 16, 1; 18*ones(8,1), (3:10)'; 19, 0; 20, 0];
    for i = 1:size(mixtures,1)
        if any(mixtures(i,1) == [1,6,7,10,11,19,20]) 
            bacnames{i,1} = bacNamegroups{1};
            if any(mixtures(i,3) == [1,6,7,10,11,19,20])
                bacnames{i,2} = bacNamegroups{1};
                bacgroups(i) = 1;
            elseif any(mixtures(i,3) ==[8,9,12:18,22])
                bacnames{i,2} = bacNamegroups{2};
                bacgroups(i) = 2;
            elseif any(mixtures(i,3) == [2,3,4,5,21])
                bacnames{i,2} = bacNamegroups{3};
                bacgroups(i) = 5;
            end 
        elseif any(mixtures(i,1) == [8,9,12:18,22]) 
             bacnames{i,1} = bacNamegroups{2};
            if any(mixtures(i,3) == [8,9,12:18,22])
                bacnames{i,2} = bacNamegroups{2};
                bacgroups(i)=4;
            elseif any(mixtures(i,1) == [2,3,4,5,21])
                bacnames{i,2} = bacNamegroups{3};
                bacgroups(i)=8;
            end 
        elseif any(mixtures(i,1) == [2,3,4,5,21]) && any(mixtures(i,1) == [2,3,4,5,21])
            bacnames{i,1} = bacNamegroups{3};
            bacnames{i,2} = bacNamegroups{3};
            bacgroups(i) = 9;
        else 
            bacnames{i,1} ='0';
        end 
    end 
end 

function [wmse]=find_wmse_error(errors, count)
    errors = reshape(errors,prod(size(errors)),1);
    perc5=prctile(errors,5,'all');
    perc95=prctile(errors,95,'all');
    errors(errors<perc5)=perc5;
    errors(errors>perc95)=perc95;
    wmse = (sum((errors).^2))/count;
end 

function [heprednew] = postprocesshe(hepred,concentrations)
    heprednew = zeros(size(hepred));
    
    uncertainty = zeros(length(size(hepred,2)),length(concentrations));
    threshold = 0.2*(size(hepred,1)*size(hepred,3)-1)/sqrt(size(hepred,1)*size(hepred,3));
    for index = 1:(size(hepred,2))
        temppred = hepred(:,index,:);
        xvar = [concentrations,concentrations];
        temppred = temppred(:);
        %remove abnormal values and zeros 
        temppred(abs(temppred)>2e3)=nan;
        % Use the threshold to remove outliers 
        tempscore = (temppred-mean(temppred(~isnan(temppred))))./std(temppred(~isnan(temppred)));
        temppred(abs(tempscore)>threshold)=nan;
        temppred(isnan(temppred)) = median(temppred(~isnan(temppred)));
        temp = temppred;
%         temppred = temppred(~isnan(temppred));
%         xvar = xvar(~isnan(temppred));
%         xvar = [xvar,0,1];
%         temppred = [temppred; 0;0];
%         %disp(xvar)
%         %disp(temppred)
%         %fit a polynomial 
%         [p,S,mu] = polyfit(xvar, temppred,3);
%         % create polynomial predictions 
%         [temp,uncertainty(index,:)] = polyval(p,concentrations,S,mu);
        heprednew(:,index,:) = reshape(temp, size(hepred,1),[]);
    end 
    
end 
function [X_removed]=remove_nan2(X)
    %check for rows and columns containing only nan values in a matrix and remove them 
    % Input 
    % X = matrix with rows or columns containing nan values
    % Output 
    % X_removed = matrix without rows or columns containing only nan values
    % 
    row_nan = find(all(isnan(X),2));
    col_nan = find(all(isnan(X),1));
    X_removed = X;
    if row_nan
        X_removed(row_nan,:)=[];
    end 
    if col_nan 
        X_removed(:,col_nan)=[];
    end
end 