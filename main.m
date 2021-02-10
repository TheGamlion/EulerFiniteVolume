clear variables;
% parameter
global n dx gamma Energy pressure c;
n =50;
dx = 1/n;
gamma = 1.4;
courant =0.6;
t_end = 0.25;
slowdown = 1;
example = 1;


% Energy EQ of state
Energy = @(rho,p,v) p/(gamma-1) + rho*v^2/2;
pressure = @(q) (gamma-1)*(q(3)-q(1)*(q(2)/q(1))^2/2);
c = @(q) sqrt(gamma*pressure(q)/q(1));


% q1 entry : density
% q2 entry : momentum
% q3 entry : Energy

%%  fill initial values

switch example
    %example 1
    case 1
        rho_l = 1;
        rho_r = 0.125;
        velocity_l = 0;
        velocity_r = 0;
        pressure_l = 1;
        pressure_r = 0.1;
    case 2
        % example 2
        rho_l = 1;
        rho_r = 1;
        velocity_l = -1;
        velocity_r = 1;
        pressure_l = 0.4;
        pressure_r = 0.4;
    case 3
        rho_l = 1;
        rho_r = 1;
        velocity_l = 0;
        velocity_r = -1;
        pressure_l = 0.4;
        pressure_r = 0.4;
    otherwise
        disp('Example not implemented');
        quit;
end

q = zeros(n,1);
for i=1:n
    if i <= n/2
        q(i,1) = rho_l;
        q(i,2) = velocity_l*rho_l;
        q(i,3) = Energy(rho_l,pressure_l,velocity_l);
    else
        q(i,1) = rho_r;
        q(i,2) = velocity_r*rho_r;
        q(i,3) = Energy(rho_r,pressure_r,velocity_r);
    end
end
% save all the simulation data
data = cell(1,4);
data(1,1) = {q(:,1)}; % density
data(1,2) = {q(:,2)}; % momentum
data(1,3) = {q(:,3)}; % energy
data(1,4) = {0}; % time

%% main loop

t = 0;
step =1 ;
while t<t_end
    % getting dt
    dt = maxdt(q);
    dt = courant* dt;
    
    for i=1:n
        if(i ==1)
%                         q(i,:) = q(i,:) - dt/dx*(F(q(i,:),q(i+1,:))-F(q(end,:),q(i,:)));
        elseif(i ==n)
%                         q(i,:) = q(i,:) - dt/dx*(F(q(i,:),q(1,:))-F(q(i-1,:),q(i,:)));
        else
            q(i,:) = q(i,:) - dt/dx*(F(q(i,:),q(i+1,:))-F(q(i-1,:),q(i,:)));
        end
        q(1,:) = q(2,:);
        q(end,:) = q(end-1,:);
    end
    
    % loop updates and copying of data
    step = step + 1;
    t = t + dt;
    disp(t);
    data(step,1) = {q(:,1)};
    data(step,2) = {q(:,2)};
    data(step,3) = {q(:,3)};
    data(step,4) = {t};
end

%% check for conservation


for i=1:step
    rho(i,:) = cell2mat(data(i,1)); %#ok<*SAGROW>
    momentum(i,:) = cell2mat(data(i,2));
    energy(i,:) = cell2mat(data(i,3));
    sum_rho(i) = dx*norm(rho(i,:),1);
    sum_mom(i) = dx*norm(momentum(i,:),1);
    sum_energy(i) = dx*norm(energy(i,:),1);
end


%%  plot
x_line = linspace(0,1,n);
t_line = cell2mat(data(:,4));
figure('Renderer', 'painters', 'Position', [10 10 1600 1200])

for i=1:step
    rho(i,:) = cell2mat(data(i,1));
    momentum(i,:) = cell2mat(data(i,2));
    energy(i,:) = cell2mat(data(i,3));
end

subplot(2,3,1)
plot1 = pcolor(x_line,t_line,rho);
plot1.LineWidth = 0.1;
colorbar
title('density');
subplot(2,3,2)
plot1 = pcolor(x_line,t_line,momentum);
plot1.LineWidth = 0.1;
colorbar
title('momentum');
subplot(2,3,3)
plot1 = pcolor(x_line,t_line,energy);
plot1.LineWidth = 0.1;
colorbar
title('energy');


x_line = linspace(0,1,n);
t_line = cell2mat(data(:,4));


for i=1:step
    subplot(2,3,[4,5,6]);
    p = plot(x_line,rho(i,:),x_line,momentum(i,:),x_line,energy(i,:));
    set(gca,'fontsize',18)
    % plot attributes
    p(1).LineWidth = 4;
    p(2).LineWidth = 4;
    p(3).LineWidth = 4;
    title('1D Euler Equasion',t_line(i));
    xlabel('x');
    refline([0 0]);
    legend('density','momentum','energy','Location','southeast');
    ylim([-2 3]);
    if i == step
        pause(1)
    else
        pause(slowdown*(t_line(i+1)-t_line(i)));
    end
end



%just animation
% figure('Renderer', 'painters', 'Position', [10 10 1600 1200])
% for i=1:step
%     rho(i,:) = cell2mat(data(i,1));
%     momentum(i,:) = cell2mat(data(i,2));
%     energy(i,:) = cell2mat(data(i,3));
% end
% 
% x_line = linspace(0,1,n);
% t_line = cell2mat(data(:,4));
% 
% 
% for i=1:step
%     rho(i,:) = cell2mat(data(i,1));
%     momentum(i,:) = cell2mat(data(i,2));
%     energy(i,:) = cell2mat(data(i,3));
%     p = plot(x_line,rho(i,:),x_line,momentum(i,:),x_line,energy(i,:));
%     set(gca,'fontsize',18)
%     % plot attributes
%     p(1).LineWidth = 4;
%     p(2).LineWidth = 4;
%     p(3).LineWidth = 4;
%     title('1D Euler Equasion',t_line(i));
%     xlabel('x');
%     refline([0 0]);
%     legend('density','momentum','energy','Location','southeast');
%     ylim([-2 3]);
%     if i == step
%         pause(1)
%     else
%         pause(slowdown*(t_line(i+1)-t_line(i)));
%     end
% end


%% functions needed

function flux = F(q1,q2)
[EW1,~] = ewev(q1);
[EW2,~] = ewev(q2);
al = min(min(EW1),min(EW2));
ar = max(max(abs(EW1)),max(abs(EW2)));

if 0 <= al
    flux = f(q1);
elseif (al <= 0) && (0<= ar)
    flux = 1/(ar-al)*(ar*f(q1)-al*f(q2)+al*ar*(q2-q1));
elseif 0 >= ar
    flux = f(q2);
end
end

function v = f(q)
global pressure;
v(1) = q(2);
v(2) = q(1)*(q(2)/q(1))^2 + pressure(q);
v(3) = q(2)/q(1)*(q(3)+ pressure(q));
end


function dt = maxdt(q)
global n dx;

max_speed =0;

for i=1:n
    [EW,~] = ewev(q(i,:));
    if max_speed < max(abs(EW))
        max_speed = max(abs(EW));
    end
    dt = dx/max_speed;
end

end

function [EW,R] = ewev(q)
global gamma pressure c; %#ok<*NUSED>
EW(1) = q(2)/q(1) -c(q);
EW(2) = q(2)/q(1);
EW(3) = q(2)/q(1) + c(q);

R =[1 1 1; q(2)/q(1)-c(q) q(2)/q(1) q(2)/q(1)+c(q);(q(3)+pressure(q))/q(1)-q(2)/q(1)*c(q) 1/2*(q(2)/q(1))^2 (q(3)+pressure(q))/q(1)+q(2)/q(1)*c(q)];
end


