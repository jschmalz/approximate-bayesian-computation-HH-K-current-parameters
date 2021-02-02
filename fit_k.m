clearvars
clear -all
clc


% load processed data so it is easier to handle
[tk,K,thigh,Khigh] = k_data_return();


% parameters
a1 = 0.01; 
a2 = 10;
a3 = 10;
b1 = 0.125;
b2 = 80;
g_k = 24;
theta_p = [a1,a2,a3,b1,b2,g_k];

% switch to log scale
P_start = log10(theta_p);


% different constant membrane potential
V = -[6,19,26,32,38,51,63,76,88];

% ABC set-up
N = 1000;
epsilon = [0.3,0.22,0.17];

theta_1 = zeros(size(P_start,2),N);
theta_2 = zeros(size(P_start,2),N);
theta_3 = zeros(size(P_start,2),N);

w_1 = ones(size(P_start,2),N);
w_2 = zeros(size(P_start,2),N);
w_3 = zeros(size(P_start,2),N);

error_1 = zeros(1,N);
error_2 = zeros(1,N);
error_3 = zeros(1,N);




% ABC-smc


%% 1
serror = 100;
disp('-#1-')
tic
parfor i = 1:N 
    error = epsilon(1)*10;
    theta_s = zeros(size(P_start));
    
    while abs(real(error)) > epsilon(1)
        % sample parameters
        theta_s = P_start + 0.5*randn(size(P_start));

        % generate profile for all constant V vlaues
        [t,g] = k_conductance_fun(V,10.^theta_s);
        % calculate distance
        error = dist_K(t,g,tk(2:end),K);
        
        
        if 1 == isnan(error)
            error = 10^7;
        end
        
        
    end
    theta_1(:,i) = theta_s;
    [1,i]
    error
    error_1(:,i) = error;
end
% calculate weight
w_1 = w_1/N;
toc

%% 2

tic
disp('-#2-')
parfor i = 1:N
    error = epsilon(2)*10;
    theta_ss = zeros(size(P_start,2),1);
    theta_s = theta_1(:,i)';
    while abs(real(error)) > epsilon(2)*(size(P_start,2) == sum(unifpdf(theta_s,P_start-(0.5*3),P_start+(0.5*3)) ~= 0 )) 
        % sample parameters
        for j = 1:size(P_start,2)
            theta_ss = randsample( theta_1(j,:), 1, true, w_1(j,:) ); % sample
            theta_s(j) = 0.8*theta_ss + 0.4*theta_ss*rand(1); % perturb
        end
        
        % generate profile for all constant V vlaues
        [t,g] = k_conductance_fun(V,10.^theta_s);
        % calculate distance
        error = dist_K(t,g,tk(2:end),K);
        
        if 1 == isnan(error)
            error = 10^7;
        end
            
    end
    % store accepted value
    theta_2(:,i) = theta_s;
    % update weight
    w_2(:,i) = normpdf(theta_s,P_start,0.5*ones(size(P_start)))./sum(  w_1.*unifpdf(ones(size(theta_1)).*theta_s',theta_1 - 0.2*sign(theta_1).*theta_1,theta_1 + 0.2*sign(theta_1).*theta_1  ),2 )'; %ones(1,N).*
    [2,i]
    error_2(:,i) = error;
end

% calculate weight
w_2 = w_2./(sum(w_2));
 
toc

%% 2

tic
disp('-#3-')
parfor i = 1:N
    error = epsilon(3)*10;
    theta_ss = zeros(size(P_start,2),1);
    theta_s = theta_2(:,i)';
    while abs(real(error)) > epsilon(3)*(size(P_start,2) == sum(unifpdf(theta_s,P_start-(0.5*3),P_start+(0.5*3)) ~= 0 )) 
        % sample parameters
        for j = 1:size(P_start,2)
            theta_ss = randsample( theta_2(j,:), 1, true, w_2(j,:) ); % sample
            theta_s(j) = 0.8*theta_ss + 0.4*theta_ss*rand(1); % perturb
        end
        
        % generate profile for all constant V vlaues
        [t,g] = k_conductance_fun(V,10.^theta_s);
        % calculate distance
        error = dist_K(t,g,tk(2:end),K);
        
        if 1 == isnan(error)
            error = 10^7;
        end
            
    end
    % store accepted value
    theta_3(:,i) = theta_s;
    % update weight
    w_3(:,i) = normpdf(theta_s,P_start,0.5*ones(size(P_start)))./sum(  w_2.*unifpdf(ones(size(theta_2)).*theta_s',theta_2 - 0.2*sign(theta_2).*theta_2,theta_2 + 0.2*sign(theta_2).*theta_2  ),2 )'; %ones(1,N).*
    [3,i]
    error_3(:,i) = error;
end

%calculate weight
w_3 = w_3./(sum(w_3));
% 
toc

%%
% save('theta_1.mat','theta_1')
% save('w_1.mat','w_1')
% save('error_1.mat','error_1')
% 
% save('theta_2.mat','theta_2')
% save('w_2.mat','w_2')
% save('error_2.mat','error_2')
% 
% save('theta_3.mat','theta_3')
% save('w_3.mat','w_3')
% save('error_3.mat','error_3')

%% theta plots

figure('DefaultAxesFontSize',14,'DefaultTextFontName','Calibri')
%a1
subplot(3,3,1)
histogram(10.^theta_1(1,:))
xlabel('a_1')
title('t = 0')
ylabel('Frequency')
xlim([0,0.05])

subplot(3,3,4)
histogram(10.^theta_2(1,:))
xlabel('a_1')
title('t = 1')
ylabel('Frequency')
xlim([0,0.05])

subplot(3,3,7)
histogram(10.^theta_3(1,:))
xlabel('a_1')
title('t = 2')
ylabel('Frequency')
xlim([0,0.05])

% a2
subplot(3,3,2)
histogram(10.^theta_1(2,:))
xlabel('a_2')
title('t = 0')
ylabel('Frequency')
xlim([0,100])

subplot(3,3,5)
histogram(10.^theta_2(2,:))
xlabel('a_2')
title('t = 1')
ylabel('Frequency')
xlim([0,100])

subplot(3,3,8)
histogram(10.^theta_3(2,:))
xlabel('a_2')
title('t = 2')
ylabel('Frequency')
xlim([0,100])

%a3
subplot(3,3,3)
histogram(10.^theta_1(3,:))
xlabel('a_3')
title('t = 0')
ylabel('Frequency')
xlim([0,40])

subplot(3,3,6)
histogram(10.^theta_2(3,:))
xlabel('a_3')
title('t = 1')
ylabel('Frequency')
xlim([0,40])

subplot(3,3,9)
histogram(10.^theta_3(3,:))
xlabel('a_3')
title('t = 2')
ylabel('Frequency')
xlim([0,40])


figure('DefaultAxesFontSize',14,'DefaultTextFontName','Calibri')
%b1
subplot(3,3,1)
histogram(10.^theta_1(4,:))
xlabel('b_1')
title('t = 0')
ylabel('Frequency')
xlim([0,0.3])

subplot(3,3,4)
histogram(10.^theta_2(4,:))
xlabel('b_1')
title('t = 1')
ylabel('Frequency')
xlim([0,0.3])

subplot(3,3,7)
histogram(10.^theta_3(4,:))
xlabel('b_1')
title('t = 2')
ylabel('Frequency')
xlim([0,0.3])

% b2
subplot(3,3,2)
histogram(10.^theta_1(5,:))
xlabel('b_2')
title('t = 0')
ylabel('Frequency')
xlim([0,1500])

subplot(3,3,5)
histogram(10.^theta_2(5,:))
xlabel('b_2')
title('t = 1')
ylabel('Frequency')
xlim([0,1500])

subplot(3,3,8)
histogram(10.^theta_3(5,:))
xlabel('b_2')
title('t = 2')
ylabel('Frequency')
xlim([0,1500])

%a3
subplot(3,3,3)
histogram(10.^theta_1(6,:))
xlabel('g_K')
title('t = 0')
ylabel('Frequency')
xlim([0,40])

subplot(3,3,6)
histogram(10.^theta_2(6,:))
xlabel('g_K')
title('t = 1')
ylabel('Frequency')
xlim([0,40])

subplot(3,3,9)
histogram(10.^theta_3(6,:))
xlabel('g_K')
title('t = 2')
ylabel('Frequency')
xlim([0,40])

%% plot
theta_p = 10.^median(theta_3,2)
[t,g] = k_conductance_fun(V,theta_p);
error = dist_K(t,g,tk(2:end),K)



% plots data with model profiles
figure('DefaultAxesFontSize',14,'DefaultTextFontName','Calibri')
plot(t,g(:,1),'k',t,g(:,2),'b',t,g(:,3),'r',...
    t,g(:,4),'g',t,g(:,5),'m',t,g(:,6),'k--',...
    t,g(:,7),'b--',t,g(:,8),'r--',t,g(:,9),'g--',...
    tk(2:end),K(:,1),'ks',tk(2:end),K(:,2),'bs',tk(2:end),K(:,3),'rs',...
    tk(2:end),K(:,4),'gs',tk(2:end),K(:,5),'ms',tk(2:end),K(:,6),'k^',...
    tk(2:end),K(:,7),'b^',tk(2:end),K(:,8),'r^',tk(2:end),K(:,9),'g^',...
    'MarkerSize',6,'LineWidth',2)


legend('V = 6 mV',...
    'V = 19 mV','V = 26 mV',...
    'V = 32 mV','V = 38 mV',...
    'V = 51 mV','V = 63 mV',...
    'V = 76 mV','V = 88 mV')
xlim([0,12])
ylim([0,22])
xlabel('Time (ms)')
ylabel('g_K - (mm mho/cm^2)')


%% HH parameters

a1 = 0.01; 
a2 = 10;
a3 = 10;
b1 = 0.125;
b2 = 80;
g_k = 24;
theta_p = [a1,a2,a3,b1,b2,g_k];

[t,g] = k_conductance_fun(V,theta_p);
error = dist_K(t,g,tk(2:end),K)



% plots data with model profiles
figure('DefaultAxesFontSize',14,'DefaultTextFontName','Calibri')
plot(t,g(:,1),'k',t,g(:,2),'b',t,g(:,3),'r',...
    t,g(:,4),'g',t,g(:,5),'m',t,g(:,6),'k--',...
    t,g(:,7),'b--',t,g(:,8),'r--',t,g(:,9),'g--',...
    tk(2:end),K(:,1),'ks',tk(2:end),K(:,2),'bs',tk(2:end),K(:,3),'rs',...
    tk(2:end),K(:,4),'gs',tk(2:end),K(:,5),'ms',tk(2:end),K(:,6),'k^',...
    tk(2:end),K(:,7),'b^',tk(2:end),K(:,8),'r^',tk(2:end),K(:,9),'g^',...
    'MarkerSize',6,'LineWidth',2)
legend('V = 6 mV',...
    'V = 19 mV','V = 26 mV',...
    'V = 32 mV','V = 38 mV',...
    'V = 51 mV','V = 63 mV',...
    'V = 76 mV','V = 88 mV')
xlim([0,12])
ylim([0,22])
xlabel('Time (ms)')
ylabel('g_K - (mm mho/cm^2)')


