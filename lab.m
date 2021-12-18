graphics_toolkit('gnuplot')
load data.mat
t = therm(7:8, :);
line_1 = therm(7, :);
line_2 = therm(8, :);
s = size(line_1);
x = 1:s(2);

figure(1, 'position',[0, 0, 800, 600]);
grid on
hold on
plot(x, line_1)
plot(x, line_2)
plot([1000, 1000], [0, 300], "k")

plot([500,500], [0, 300], "k")
xlim([0, 2200])



left = 500
right = 1000


data_1 = line_1(left:right);
data_2 = line_2(left:right);
indexes = x(left:right);
left = 500
right = 1000


pkg load interval

## ���������� ���� � �������� ���������� ������������ ���������
addpath(genpath('./m'))
n = 5
start_i = 1
step = 80
end_i = step * n
inds = [start_i:step:end_i]

figure
plot(indexes(start_i:end_i), data_1(start_i:end_i), "-")
hold on
plot(indexes(inds), data_1(inds), "o")



## ���������� ������ ���������� ������������ ��������� 
##     y = X * beta = beta1 + beta2 * x 
## � ������������ beta2 >= 0

x = indexes(inds);        # ���������� ������������ �������
y = data_1(inds);        # ����� �������������� ����

epsilon = zeros(1, n) + 0.1;  # ������� ������� ������ ��� y_i

display(x)
display(y)
display(epsilon)

x = reshape(x, size(x)(2), size(x)(1));
y = reshape(y, size(y)(2), size(y)(1));
epsilon = reshape(epsilon, size(epsilon)(2), size(epsilon)(1));

X = [ x.^0 x ];                               # ������� �������� ���������� ��� beta1 � beta2
lb = [-inf 0];                                # ������ ������� beta1 � beta2
irp_temp = ir_problem(X, y, epsilon, lb);    # �������� ����������, ���������� �������� ������ 
                                              #               ���������� ������������ ���������

display(x)
display(y)
display(epsilon)


## ������ ������������ ���������
figure('position',[0, 0, 800, 600]);

ir_scatter(irp_temp);   
grid on
set(gca, 'fontsize', 12)

## �������� ��������� ���

b_lsm = (X \ y)'  
MNK_line = [b_lsm(1) + b_lsm(2) * min(x), b_lsm(1) + b_lsm(2) * max(x)];


figure('position',[0, 0, 800, 600]);

plot(x, y, "o")
hold on
plot([min(x),  max(x)], MNK_line)
grid on
# xlim([0 22])
# ylim([-10 260])
set(gca, 'fontsize', 12)
#title('Steam generator performance');
#xlabel('Fuel consumption');
#ylabel('Steam quantity');
legend("", "Prediction")


## ����������� ������������� ��������������� ���������

#���
eps = epsilon
m = size(x)(1)
C = zeros(1, m + 2);
for i = 1:m
  C(i) = 1;
end
display(C)
A = zeros(2*m, m+2);

for i = 1:m
  A(2 * i - 1, i) = eps(i);
  A(2 * i, i) = eps(i);

  A(2 * i - 1, m + 1) = 1;
  A(2 * i, m + 1) = -1;

  A(2 * i - 1, m + 2) = x(i);
  A(2 * i, m + 2) = -x(i);

end

display(A)

B = zeros(1, 2*m);
for i = 1:m
  B(2 * i - 1) = y(i);
  B(2 * i) = -y(i);
end

display(B)

lb = zeros(1, m+2);
for i = 1:m
  lb(i) = 1;
end

lb(m+2) = -inf;

display(lb)

ctype = "";
for i = 1:2 * m
  ctype(i) = 'L';
end

display(ctype)

vartype = "";
for i = 1:m + 2
  vartype(i) = 'C';
end

display(vartype)
sense = 1
w = glpk(C,A,B,lb,[],ctype,vartype,sense)


scale = max(w(1:n))
for i = 1:n
  eps(i) = epsilon(i) * scale;
end

display(x)
display(y)
display(eps)


X = [ x.^0 x ];                               # ������� �������� ���������� ��� beta1 � beta2
lb = [-inf 0];                                # ������ ������� beta1 � beta2
irp_temp = ir_problem(X, y, eps, lb);    # �������� ����������, ���������� �������� ������ 
                                              #               ���������� ������������ ���������



## ������ ������������ ���������
figure('position',[0, 0, 800, 600]);

ir_scatter(irp_temp);   

set(gca, 'fontsize', 12)



## ������� ��������������� ��������� ������ ���������� ������������ ���������
vertices = ir_beta2poly(irp_temp)

## ������� � �������� ��������� ������� ��������������� ��������� 
[rhoB, b1, b2] = ir_betadiam(irp_temp)


## ������� ������������ ������ ���������� ������ y = beta1 + beta2 * x 
b_int = ir_outer(irp_temp)

## �������� ������ ���������� 
b_maxdiag = (b1 + b2) / 2    # ��� �������� ���������� ��������� ��������������� ���������

b_gravity = mean(vertices)   # ��� ����� ������� ��������������� ��������� 

b_lsm = (X \ y)'             # ������� ���������� ���������

## ����������� ������������� ������� ������������ ������ ��������������� ���������
figure('position',[0, 0, 800, 600]);
ir_plotbeta(irp_temp)
hold on
ir_plotrect(b_int,'r-')
grid on
set(gca, 'fontsize', 12)
xlabel('\beta_1')
ylabel('\beta_2')
title('Information set')

## �������� ������
plot(b_maxdiag(1), b_maxdiag(2), 'ro')
plot(b_gravity(1), b_gravity(2), 'k+')
plot(b_lsm(1), b_lsm(2), 'gx')
legend("", "", "enclosure", "maxdiag",  "gravity", "lsm")


## ����������� ������������� �������� ���������� ������������ ��� ������ y = beta1 + beta2 * x
figure('position',[0, 0, 800, 600]);
xlimits = [400 1000];
ir_plotmodelset(irp_temp, xlimits)     # ������� ���������� ������������

hold on
ir_scatter(irp_temp,'bo')              # ������������ ���������
ir_plotline(b_maxdiag, xlimits, 'r-')   # ����������� � �����������, ���������� ��� ����� ���������� ��������� ��
ir_plotline(b_gravity, xlim, 'b--')     # ����������� � �����������, ���������� ��� ����� ������� ��  
# ir_plotline(b_lsm, xlim, 'b--')         # ����������� � �����������, ���������� ���
# ir_scatter(ir_problem(Xp,ypmid,yprad),'ro')

grid on
set(gca, 'fontsize', 12)

## ����������� ������������� �������� ���������� ������������ ��� ������ y = beta1 + beta2 * x
figure('position',[0, 0, 800, 600]);
xlimits = [400 1000];
ir_plotmodelset(irp_temp, xlimits)     # ������� ���������� ������������

hold on
ir_scatter(irp_temp,'bo')              # ������������ ���������
ir_plotline(b_maxdiag, xlimits, 'r-')   # ����������� � �����������, ���������� ��� ����� ���������� ��������� ��
ir_plotline(b_gravity, xlim, 'b--')     # ����������� � �����������, ���������� ��� ����� ������� ��  
grid on

grid on
set(gca, 'fontsize', 12)

xlim([450, 600])
ylim([70, 95])

## ����������� ������������� �������� ���������� ������������ ��� ������ y = beta1 + beta2 * x
figure('position',[0, 0, 800, 600]);
xlimits = [400 1000];
ir_plotmodelset(irp_temp, xlimits)     # ������� ���������� ������������

hold on
ir_scatter(irp_temp,'bo')              # ������������ ���������
ir_plotline(b_maxdiag, xlimits, 'r-')   # ����������� � �����������, ���������� ��� ����� ���������� ��������� ��
ir_plotline(b_gravity, xlim, 'b--')     # ����������� � �����������, ���������� ��� ����� ������� ��  


grid on
set(gca, 'fontsize', 12)


xlim([480, 520])
ylim([78, 80])

## �������� y, ������������� � ������� ������ y = beta1 + beta2 * x � ������ ������������
yp0 = ir_predict(irp_temp, X)       # ������������ ������� �������� y � ������ x

yp0mid = mean(yp0,2)                 # ������� �������� ���������� ����������
yp0rad = 0.5 * (yp0(:,2) - yp0(:,1)) # ������ ���������� ����������

yp0rad_rel = 100 * yp0rad ./ yp0mid  # ������������� �������� ���������������� ��������� � ���������

## �������� y, ������������� � ������� ������ y = beta1 + beta2 * x
xp = [250; 450; 600; 950; 1800]
Xp = [xp.^0 xp];

yp = ir_predict(irp_temp, Xp)         # ������������ ������� �������� y � ������ xp
ypmid = mean(yp,2)                     # ������� �������� ���������� ����������
yprad = 0.5 * (yp(:,2) - yp(:,1))      # ������ ���������� ����������

yprad_relative = 100 * yprad ./ ypmid  # ������������� �������� ���������������� ��������� � ���������


## ����������� ������������� �������� ���������� ������������ ��� ������ y = beta1 + beta2 * x
figure('position',[0, 0, 800, 600]);
xlimits = [200 2000];
ir_plotmodelset(irp_temp, xlimits)     # ������� ���������� ������������

hold on
ir_scatter(irp_temp,'bo')              # ������������ ���������
ir_plotline(b_maxdiag, xlimits, 'r-')   # ����������� � �����������, ���������� ��� ����� ���������� ��������� ��
ir_plotline(b_gravity, xlim, 'b--')     # ����������� � �����������, ���������� ��� ����� ������� ��  
ir_scatter(ir_problem(Xp,ypmid,yprad),'ro')

grid on
set(gca, 'fontsize', 12)



# ����� ��������� �����
MY_EPS = 0.00001;
## �������� y, ������������� � ������� ������ y = beta1 + beta2 * x � ������ ������������
irp_temp;
cur_x = irp_temp.y;
cur_eps = irp_temp.epsilon;

yp0 = ir_predict(irp_temp, X); 

for i = 1:n
x_top = cur_x(i) + cur_eps(i);
x_bot = cur_x(i) - cur_eps(i);

y_top = yp0(i, 2);
y_bot = yp0(i, 1);

if abs(y_top - x_top) < MY_EPS
    display(i)
    
end

if abs(y_bot - x_bot) < MY_EPS
    display(i)
end


end
i = 1
i = 2
i = 5

border_x = [-5, 5];
border_y = [-0.3, 0.3];

for i = 1:m

cur_point_x = x(i);
cur_point_y = y(i);

figure('position',[0, 0, 800, 600]);
xlimits = [400 1000];
ir_plotmodelset(irp_temp, xlimits)     # ������� ���������� ������������
hold on
ir_scatter(irp_temp,'bo')              # ������������ ���������
grid on
set(gca, 'fontsize', 12)

xlim([cur_point_x + border_x(1), cur_point_x + border_x(2)]);
ylim([cur_point_y + border_y(1), cur_point_y + border_y(2)]);


end









