%% Question 9.1

index_counter = 1; % Used for indexing the compliance and delta P vectors

%9.1a
for delta_V = 0.1:0.1:0.4
   
time_delay = 20.*delta_V; % Used to calculate the duration of each input injection
t = 0:0.1:180; %time axis values
r = 0.05*(heaviside(t-15)-heaviside(t-(15+time_delay))); %square unit pulse
tau = 60; %Calculated from given graph
k = 5000; % After playing with the values of k at 100,1000,2000,...5000, this seems to be the closest value to the graph
% with the peak spike at 43.67 compared to 43.71 on the homework pdf

% Parameters for state-space model. Representative of delta P

A = -1/tau;
B = k/tau;
C = 1;
D = 0;

% Calculating and plotting the system using the state space model, the
% input and the time vector
sys = ss(A,B,C,D);
P = lsim(sys,r,t); 
plot(t,P+12.5) % Must increase the plot by P_0 in order to mimic graph on pdf
hold on

plot(t,r);
hold on

P_max = max(P);
dP(index_counter) = P_max; % I only shifted my graph during the plotting, so there is no need to subtract P_0
Compliance(index_counter) = delta_V/dP(index_counter);

index_counter = index_counter + 1;
end

title('Pressure vs Time');
legend ('dV = 0.1 mL','dV = 0.2 mL','dV = 0.3 mL','dV = 0.4 mL');
xlabel('Time [s]');
ylabel('P(t) [mmHg]');

disp('The figure shown is the plot for pressure vs time for each saline injection amount')
%Question 9.1b

num = 5000;
denom = [tau 1];

disp('Shown below is the transfer function used in the system')
tfsys = tf(num,denom) % This is my transfer function

%Question 9.1c

disp('Compliance and delta P are in order of 0.1,0.2,0.3,0.4 V descending')
Compliance'
dP



