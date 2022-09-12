clear all;
clc;
%{

    Name = Aman Kumar
    Roll No. = 214103404
    Specialization = Machine Design

Problem: Use a Binary Coded GA to minimize the function f(x1,x2)=x1 + x2
-2x1^2-x2^2+x1x2, in the range of 0.0<=x1,x2<=0.5. Using a random
population size N=6,a single point crossover with probability pc =
1.0,neglect mutation and assume 5 bits for each variable.
%}

%{ 
See if the function is minmization or maximization , if minimization then
convert it to a maximization problem.
%}
m=1; %minimization problem
gen = 1;% count of generations

%Define the bounds of each variable.
N_var = 2;%Number of variables
bound = [0 0.5;0 0.5];%first row for x1 , second for x2.

%Choose initial population
N = 6;
bits =[5;5];
length = sum(bits);
population = zeros(N,length);
p=zeros(1,length);
r = binornd(1,0.6);
for i =1:N
    for j=1:length
        p(j)=binornd(1,0.5);
    end
    population(i,:) = p(1,:);
end
error =10;
while(error>0||gen<50)

%to make real variables out of the selected population
real_pop = zeros(N,N_var);
for i =1:N
    k=1;
    for j=1:N_var
        real_pop(i,j) = reall(population(i,k:(bits(j)+k-1)),bound(j,1),bound(j,2),bits(j));
        k = k+bits(j);
    end
end

%evaluation of fitness value and formulation of mating pool
fitness = zeros(N,1);
for i=1:N
    fitness(i) = objfn(real_pop(i,1),real_pop(i,2),m);
end
total = sum(fitness);
prob = zeros(N,1);
for i=1:N
prob(i)=fitness(i)/total;
end
cumProb = zeros(N,1);
for i=1:N
    if(i==1)
        cumProb(i)=prob(i);
    else cumProb(i)=cumProb(i-1)+prob(i);
    end
end
randomNo = rand(N,1);
matingPool = zeros(N,length);
for i=1:N
    k=1;
    while(randomNo(i)>cumProb(k))
        k = k+1;
    end
    matingPool(i,:)=population(k,:);
end

%Single Point CrossOver
children = zeros(N,length);
count=0;
while(count<N)
    Parent1(1:length) = matingPool(randi([1,N],1,1));
    Parent2(1:length) = matingPool(randi([1,N],1,1));
    CrossSite = randi([1,length-1],1,1);
    children(count+1,1:CrossSite)= Parent1(1:CrossSite);
    children(count+1,CrossSite+1:length)= Parent2(CrossSite+1:length);
    children(count+2,1:CrossSite)= Parent2(1:CrossSite);
    children(count+2,CrossSite+1:length)= Parent1(CrossSite+1:length);
    count = count+2;
end
population = children;
error = mse(children,fitness,N,N_var,bits,bound,m);
gr(gen)=error;
gen = gen+1;
end

final= zeros(N,N_var);
for i =1:N
    k=1;
    for j=1:N_var
        final(i,j) = reall(population(i,k:(bits(j)+k-1)),bound(j,1),bound(j,2),bits(j));
        k = k+bits(j);
    end
end

%% Part 5 Writing output in the output file
fileOD = fopen('output.txt','w');
fprintf(fileOD,'%s\t\t\t\t%s\t\t%s\n','x1','x2','f(x1,x2)');
for i=1:N
    fprintf(fileOD,'%f \t\t %f \t\t %f\n',final(i,1),final(i,2),objfn(final(i,1),final(i,2),0));
end

fileMSE = fopen('MSE.txt','w');
fprintf(fileMSE,'%s %s\n','MSE  ','   Iterations');
for i=1:gen-1
    fprintf(fileMSE,'%f \t %d \n',gr(i),i);
end
plot(1:gen-1,gr);
%Function for MSE
function err = mse(x,xold,N,N_var,bits,bound,m)
temp = zeros(N,N_var);
val = zeros(N,1);
 for i =1:N
    k=1;
    for j=1:N_var
        temp(i,j) = reall(x(i,k:(bits(j)+k-1)),bound(j,1),bound(j,2),bits(j));
        k = k+bits(j);
    end
    val(i) = objfn(temp(i,1),temp(i,2),m);
 end
err = (sum(val.*val)/N)-(sum(xold.*xold)/N);
end


%function to decode real values from the binary values
function realVal = reall(x,xmin,xmax,len)

sum = 0;
for i= 1:len
    sum = sum + 2^(len-i)*x(i);
end

realVal = xmin + ((xmax-xmin)/(2^len-1))*sum;

end


%function for objective function evaluation(fitness value)
function fvalue = objfn(x1,x2,m)
value = x1+x2-2*x1*x1-x2*x2+x1*x2;
if(m==0)
    fvalue = value; % in case of maximization(m==0)
else fvalue = 1/(1+value*value); %in case of minimization(m!=0)
end

end


