2.1 The matlab code of task1
clc;clear;close all;
syms C1 C2 C3 S1 S2 S3 L1 L2 L3
T01 = [C1 -S1 0 0
       S1 C1 0 0
       0 0 1 0
       0 0 0 1]
T12 = [C2 -S2 0 L1
       0 0 -1 0
       S2 C2 0 0
       0 0 0 1]
T23 = [C3 -S3 0 L2
       S3 C3 0 0
       0 0 1 0
       0 0 0 1]
T34 = [1 0 0 L3
       0 1 0 0
       0 0 1 0
       0 0 0 1]
T04 = simplify(T01*T12*T23*T34)
A1 = [
      L3*(-S1*C2*C3+S1*S2*S3)-S1*C2*L2-S1*L1
      L3*(C1*C2*C3-C1*S2*S3)+C1*C2*L2+C1*L1
      0
     ]
A2 = [
      L3*(-C1*S2*C3-C1*C2*S3)-C1*S2*L2
      L3*(-S1*S2*C3-S1*C2*S3)-S1*S2*L2
      L3*(C2*C3-S3*S2)+C2*L2
      ]
A3 = [
      L3*(-C1*C2*S3-C1*S2*C3)
      L3*(-S1*C2*S3-S1*S2*C3)
      L3*(-S2*S3+C3*C2)
      ]
M1 = T04(1:3, 1:3).'
AA = [A1 A2 A3]
J4theta = simplify(  M1 * AA)


2.2 The matlab code of task2
clc;clear;close all;
syms g s1 c1 N11 IC11 W11d W11 m1 w1 h1 l1 theta1 theta1d theta1dd
% task 2
IC11 = [m1/12*(h1^2+l1^2) 0 0;
        0 m1/12*(w1^2+h1^2) 0;
        0 0 m1/12*(l1^2+w1^2)]
W11 = [0;0;theta1d]
W11d = [0;0;theta1dd]
N11 = simplify(IC11*W11d+cross(W11,(IC11*W11)))
V11d = [g*s1;g*c1;0]
V1C1d = [(-1/2)*l1*theta1d^2+g*s1;(1/2)*l1*theta1dd+g*c1;0]
F11 = m1*V1C1d

syms F22 l2 c2 s2 R21 theta2 theta2d theta2dd m2 w2 h2 l2
R21 = [c2 0 s2;-s2 0 c2;0 -1 0]
W22 = simplify(R21*W11+[0;0;theta2d])
W22d = simplify(R21*W11d + cross(R21*W11,[0;0;theta2d])+[0;0;theta2dd])
V22d = simplify(R21*(cross(W11d,[l1;0;0])+cross(W11,cross(W11,[l1;0;0]))+V11d))
V2C2d = simplify(cross(W22d,[(1/2)*l2;0;0])+cross(W22,cross(W22,[(1/2)*l2;0;0]))+V22d)
F22 = m2*V2C2d
IC22 = [m2/12*(h2^2+l2^2) 0 0;
        0 m2/12*(w2^2+h2^2) 0;
        0 0 m2/12*(l2^2+w2^2)]
N22 = simplify(IC22*W22d+cross(W22,(IC22*W22)))

syms F33 l3 c3 s3 R32 theta3 theta3d theta3dd m3 w3 h3 l3
R32 = [c3 s3 0;-s3 c3 0;0 0 1]
W33 = simplify(R32*W22+[0;0;theta3d])
W33d = simplify(R32*W22d + cross(R32*W22,[0;0;theta3d])+[0;0;theta3dd])
V33d = simplify(R32*(cross(W22d,[l2;0;0])+cross(W22,cross(W22,[l2;0;0]))+V22d))
V3C3d = simplify(cross(W33d,[(1/2)*l3;0;0])+cross(W33,cross(W33,[(1/2)*l3;0;0]))+V33d)
F33 = m3*V3C3d
IC33 = [m3/12*(h3^2+l3^2) 0 0;
        0 m3/12*(w3^2+h3^2) 0;
        0 0 m3/12*(l3^2+w3^2)]
N33 = simplify(IC33*W33d+cross(W33,(IC33*W33)))

syms 
R23 = R32.'
R12 = R21.'
f33 = F33
n33 = simplify(N33+cross([(1/2)*l3;0;0],F33))
f22 = R23*f33+F22
n22 = simplify(N22+R23*n33+cross([(1/2)*l2;0;0],F22)+cross([l2;0;0],R23*f33))
f11 = simplify(R12*f22+F11)
n11 = simplify(N11+R12*n22+cross([(1/2)*l1;0;0],F11)+cross([l1;0;0],R12*f22))

t1 = simplify(n11(3,:))
t2 = simplify(n22(3,:))
t3 = simplify(n33(3,:))
2.3 The matlab code of task3
% task3_1
% M(theta)
% t3 dd
[uthetat33dd,L33dd] = coeffs(t3, theta3dd)
[uthetat32dd,L32dd] = coeffs(t3, theta2dd)
[uthetat31dd,L31dd] = coeffs(t3, theta1dd)
% t2 dd
[uthetat23dd,L23dd] = coeffs(t2, theta3dd)
[uthetat22dd,L22dd] = coeffs(t2, theta2dd)
[uthetat21dd,L21dd] = coeffs(t2, theta1dd)
% t1 dd
[uthetat13dd,L13dd] = coeffs(t1, theta3dd)
[uthetat12dd,L12dd] = coeffs(t1, theta2dd)
[uthetat11dd,L11dd] = coeffs(t1, theta1dd)
%
Mtheta = [uthetat11dd(1,1) 0 0
          0 uthetat22dd(1,1) uthetat23dd(1,1)
          0 uthetat32dd(1,1) uthetat33dd(1,1)]
% V(theta)
% t3 d d^2
[uthetat33d,L33d] = coeffs(t3, theta3d)
[uthetat32d,L32d] = coeffs(t3, theta2d)
[uthetat31d,L31d] = coeffs(t3, theta1d)
% t2 d d^2
[uthetat23d,L23d] = coeffs(t2, theta3d)
[uthetat22d,L22d] = coeffs(t2, theta2d)
[uthetat21d,L21d] = coeffs(t2, theta1d)
% t1 d d^2
[uthetat13d,L13d] = coeffs(t1, theta3d)
[uthetat12d,L12d] = coeffs(t1, theta2d)
[uthetat11d,L11d] = coeffs(t1, theta1d)
% Vtheta
Vtheta = [uthetat11d(1,1)*theta1d
         uthetat21d(1,1)*theta1d^2 + uthetat22d(1,1)*theta2d^2 + uthetat22d(1,2) * theta2d + uthetat23d(1,1)*theta3d^2
         uthetat31d(1,1)*theta1d^2 + uthetat32d(1,1)*theta2d^2]
% G(theta)
Gtheta = [simplify(t1 - Mtheta(1,1)*theta1dd - Vtheta(1,:))
     simplify(t2 - Mtheta(2,:)*[theta1dd;theta2dd;theta3dd] - Vtheta(2,:))
     simplify(t3 - Mtheta(3,:)*[theta1dd;theta2dd;theta3dd] - Vtheta(3,:))]

% task3_2
J4theta = [0 s3*l2 0
            0 c3*l2 l3
            -l1-l2*c2-l3*c2*c3+l3*s2*s3 0 0]
J4thetaD = [0 c3*l2*theta3d 0
            0 -s3*l2*theta3d 0
            l2*s2*theta2d+l3*s2*c3*theta2d+l3*c2*s3*theta3d+l3*c2*s3*theta2d+l3*s2*c3*theta3d 0 0]
% Mxtheta
Mxtheta = (J4theta^(-1)).' * Mtheta * J4theta^(-1)
% Vxtheta
Vxtheta = (J4theta^(-1)).' * simplify(Vtheta-Mtheta*J4theta^(-1)*J4thetaD*[theta1d;theta2d;theta3d])
% Gxtheta
Gxtheta = (J4theta^(-1)).' * Gtheta


