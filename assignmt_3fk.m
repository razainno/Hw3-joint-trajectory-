

syms q1 q2  d3 d1 a1
a1= FK(q1, d1, 0, -90);
a2= FK(q2, 0, 1, 90);

a3= FK(0, d3, 0, 0);
 
k=[0;0;1];
R01=[cos(q1),0,sin(q1);sin(q1),0,-cos(q1);0,-1,0];
R12=[cos(q2),0,sin(q2);sin(q2),0,cos(q2);0,1,0];
R23=[1,0,0;0,1,0;0,0,1];
R02=R01*R12;

 %*******************angular jacobian***************************%     
 jw1=k;
 jw2=R01*k;
 jw3=[0;0;0];
 jw=[jw1,jw2,jw3];

%*****************************linear jacobian********************?%
  O01=[0;0;1];
  O12=[cos(q2);sin(q2);0];
  O23=[0;0;d3];
  O03=O01+O12+O23;
   
  jv1=cross(k,O03);
  jv2=cross(R01*k,O01);
  jv3=R02*k;
  jv=[jv1,jv2,jv3];
  j=[jv;jw];
  
  %************************************sigularity**********************************************%
  
  D=det(jv);
    
  
  %*******************************velocities**********************************************%
   syms t  
  q1=sin(t);
  q2=cos(2*t);
  d3=sin(3*t);
  O01=[0;0;1];
  O12=[cos(q2);sin(q2);0];
  O23=[0;0;d3];
  O03=O01+O12+O23;
   
  R011=[cos(q1),0,sin(q1);sin(q1),0,-cos(q1);0,-1,0];
  R121=[cos(q2),0,sin(q2);sin(q2),0,cos(q2);0,1,0];
  R231=[1,0,0;0,1,0;0,0,1];
  R021=R011*R121;
  jw1=k;
  jw2=R011*k;
  jw3=[0;0;0];
  jw=[jw1,jw2,jw3];
  jv1=cross(k,O03);
  jv2=cross(R011*k,O01);
  jv3=R021*k;
  jv=[jv1,jv2,jv3];
  j=[jv;jw]
  dq= [diff(q1);diff(q2);diff(d3)];
   v = j*dq;
   vx(t)= v(1,:);
   vy(t)= v(2,:);
   vz(t)= v(3,:);
   
   
   subplot(4,1,1);
fplot(@(t) vx(t))

  subplot(4,1,2);
fplot(@(t) vy(t))

  subplot(4,1,3);
fplot(@(t) vz(t))


subplot(4,1,4);
hold on
fplot(@(t) vx(t));
fplot(@(t) vy(t));
fplot(@(t) vz(t));
hold off

syms q1 q2  d3 d1 a1


foward_kinamatics = a1*a2*a3;
function a=FK(theta,d,a,alpha)
a=[cos(theta)  -sin(theta)*cos(alpha) sin(theta)*sind(alpha)  a*cos(theta);
   sin(theta)  cos(theta)*cosd(alpha)  sind(alpha)*cos(theta)  a*sin(theta);
   0                sind(alpha)              cosd(alpha)            d      ;
   0                   0                      0                  1    ];

end

