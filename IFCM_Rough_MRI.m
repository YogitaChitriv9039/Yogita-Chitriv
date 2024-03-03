clc;clear all;close all;

%f=fopen('t1_icbm_normal_3mm_pn3_rf40.rawb');
f=fopen('t1_icbm_normal_1mm_pn3_rf20.rawb');

a=fread(f);
q=reshape(a,181,217,181);
q1=q(:,:,91);
q1=q1';
q1=q1(end:-1:1,:);
I=q1;

double_th = doubleThreshold(I,25,150);
level = graythresh(I);
otsu_img = singleThreshold(I,level*255);
%Combine the results

out_img1 = logical(otsu_img + double_th);
se = strel('disk',8);
%dilated_img = imdilate(out_img1,se);
eroded_img = imerode(out_img1,se);
dilated_img = imdilate(eroded_img,se);
filled_img = imfill(dilated_img,'holes');
out_img =(double(filled_img) .*I);
o_img=im2double(uint8(out_img));
subplot(1,2,1),imshow(I,[])
subplot(1,2,2),imshow(out_img,[])
for l1=0:255
count1(l1+1)=sum(sum(out_img==l1));
end

expanse=15;

[nr,nc]=size(out_img);
%%Fuzzy Image
%out_img=im2double(uint8(out_img));
Y=out_img/max(max(out_img));
Y_h=(1-Y-((1-Y)./(1+2*Y)));
 Y_n=1-Y-Y_h;
 
[Hg,x]=imhist(o_img,256);
%Hg=[0; Hg(2:256)];
Hn=Hg;
f1=1/9.0*[1 1 1;1 1 1;1 1 1];
% yf=conv2(Y,f1);
% yf=yf(2:nr+1,2:nc+1);
% Y_n=1-(Y + Y_h);   %non-membership
% 
% Y_h=abs(Y-yf); 
% Y_h=(1-Y).*Y_h.^0.4;     
%mx=max(max(Y_h));
for i=1:nr
    for j=1:nc
        if out_img(i,j)==0;
            Y_n(i,j)=0;
        else
            Y_n(i,j)=((1-Y(i,j))./(1+2*Y(i,j)));%Sugeno
            %Y_h(i,j)=Y_h(i,j)/mx;     %Mushrif sir
            %Y_n(i,j)=(1-Y(i,j).^2).^(1/2);%Yagers
        end
    end
end
for i=1:nr
    for j=1:nc
        if out_img(i,j)==0;
            Y_h(i,j)=0;
        else
            Y_h(i,j)=1-Y(i,j)-Y_n(i,j);
        end
    end
end
% 
for i=2:nr-1
    for j=2:nc-1
        kr=i-1:i+1;
        kc=j-1:j+1;
        sm1=0;

        P_m=Y(i,j); 
        P_n=Y_n(i,j); 
        P_h=Y_h(i,j); 
        
        Q_m=Y(kr,kc); 
        Q_n=Y_n(kr,kc); 
        Q_h=Y_h(kr,kc); 
        tmp=(P_m-Q_m).^2 + (P_n-Q_n).^2 + (P_h-Q_h).^2 ;
        sm1=1/2*(sqrt(sum(sum(tmp))));
        expnss(i,j)=sm1;


    end
end
Hn=double(Hn);
sdv=std2(expnss);
D=exp(-0.5*(expnss/sdv).^2);
for i=2:nr-1
    for j=2:nc-1

  
 % D(i,j)=exp(-0.5*(sm1/40)^2); %25
        % fHn(hx(i,j)+1)=fHn(hx(i,j)+1)+D(i,j);
          Hn(out_img(i,j)+1)=Hn(out_img(i,j)+1)+D(i,j);
          end
end

x1=0:255;
xI=0:255;
[HY,HC]=imhist(Y,256);
Hg(1)=0;
Hn(1)=0;
count1(1)=0;
f2=fspecial('gaussian',[20,1],5);
smt_Hn=conv(Hn,f2,'same');

smt_HY=conv(Hg,f2,'same');
figure,plot(x1,Hn,x1,Hg)    
rI=1-((abs(Hg')+eps)./(abs(Hn')+eps));
 figure,plot(xI,rI)
% val=1;

%figure,plot(x1,smt_HY,x1,smt_Hn) 
%rI1=1-((abs(HY))./(abs(Hn)+eps));
%figure,plot(x1,smt_HY)
[pks,locs]=findpeaks(rI,'minpeakdistance',40);
 
 v = locs';
 v=v/255;
 v=[0.01;v];
%v = [0.01;0.28;0.65;0.9];
  C=length(v);
  v_h=(1-v-((1-v)./(1+2*v)));
v_n=1-v_h-v;
% fl=1/3.0*[1 1 1];
% vf=conv(v,fl,'same');
% v_h=abs(v-vf); 
% v_h=(1-v).*v_h.^0.4;     %hesitency
% mh=max(v_h);
% v_h=v_h/mh;
% v_n=1-(v + v_h);   %non-membership
Y_h=(1-Y-((1-Y)./(1+2*Y)));
 Y_n=1-Y-Y_h;
 
 



[M N]=size(out_img);
U(:,:,1)=zeros(size(Y));
U(:,:,2)=Y;
U(:,:,3)=Y_h;
U(:,:,4)=Y_n;

Uold(:,:,1)=zeros(size(Y));
Uold(:,:,2)=zeros(size(Y));
Uold(:,:,3)=zeros(size(Y));
Uold(:,:,4)=zeros(size(Y));
epsilon = 0.0001;
 p=2;
 v_old=zeros(size(v));
 itt=1;
while((sum((v-v_old).^2))>=epsilon) 
%while((abs(v(1)-v_old(1))>=epsilon)&&(abs(v(2)-v_old(2))>=epsilon)&&(abs(v(3)-v_old(3))>=epsilon)&&(abs(v(4)-v_old(4)))>=epsilon)
%while((abs(U-Uold))>=epsilon), 
    itt=itt+1;
    
    disp(['itteration ' num2str(itt)]);
    
      nom=zeros(C,1); den=zeros(C);

      for m=1:M
          for n=1:N
%              Uold=U;
             for i=1:C
                  x=1/2*sqrt((Y(m,n)-v(i)).^2+(Y_h(m,n)-v_h(i)).^2+(Y_n(m,n)-v_n(i)).^2);
                  %x=dist_ifs(Y(m,n),v(i),Y_h(m,n),v_h(i),Y_n(m,n),v_n(i),'ykd_h');
                  dent=0;
                  %x=D(m,n);
                  for j=1:C
                      y=1/2*sqrt((Y(m,n)-v(j)).^2+(Y_h(m,n)-v_h(j)).^2+(Y_n(m,n)-v_n(j)).^2);
                      % y=dist_ifs(Y(m,n),v(j),Y_h(m,n),v_h(j),Y_n(m,n),v_n(j),'ykd_h');
                      if(abs(y)<eps),
                          y=eps; 
                      end
                      dent=dent+(x/y)^(2/(p-1));
                  end
                  if(abs(dent)<eps), dent=eps; end
                  U(m,n,i) = (1/ dent);
             end
             
          end
      end
     
         
       Up = (U).^p; 
       
    for m=1:M
          for n=1:N
              v_old=v;
             for i=1:C
                
                 nom(i,:)=nom(i,:)+Up(m,n,i).*Y(m,n);
                 den(i)=den(i)+Up(m,n,i);   
                 if(abs(den(i))<eps), den(i)=eps; end
                 v(i,:)=nom(i,:)/den(i);
             end
          end
    end
      
    end

 
 
    Up1=U(:,:,1);
      Up2=U(:,:,2);
        Up3=U(:,:,3);
       Up4=U(:,:,4);
%    
%     
    
%     subplot(2,2,3), imshow(B,[]), title('Estimated biasfield');
%     subplot(2,2,4), imshow(Y-B), title('Corrected image');

figure,subplot(2,2,1),imshow(Up1)
subplot(2,2,2),imshow(Up2)
subplot(2,2,3),imshow(Up3)
subplot(2,2,4),imshow(Up4)


U1=U;
 

for r=2:M-1
    for c=2:N-1
        
        for k=1:C
            l=[U(r,c,1),U(r,c,2),U(r,c,3),U(r,c,4)];
            lar=max(l);
            if U(r,c,k)>=lar
                U1(r,c,k)=1;
            else
                 U1(r,c,k)=0;
            end
        end
    end
end


figure,subplot(2,2,1),imshow(U1(:,:,1))
subplot(2,2,2),imshow(U1(:,:,2))
subplot(2,2,3),imshow(U1(:,:,3))
subplot(2,2,4),imshow(U1(:,:,4))
P81=U1(:,:,1);
P82=U1(:,:,2);
P83=U1(:,:,3);
P84=U1(:,:,4);

F=P81+P82*2+P83*3+P84*4;
kk=U(:,:,1)+U(:,:,2)+U(:,:,3)+U(:,:,4);

 figure, 
    subplot(1,2,1), imshow(out_img,[]), title('input image');
    subplot(1,2,2), imshow(F,[]), title('Segmented Image')
bck_pixels=0;
for i=2:M-1
    for j=2:N-1
        if F(i,j)==1
            bck_pixels=bck_pixels+1;
        end
    end
end
bck_pixels

CSF_pixels=0;
for i=2:M-1
    for j=2:N-1
        if F(i,j)==2
            CSF_pixels=CSF_pixels+1;
        end
    end
end
CSF_pixels

GM_pixels=0;
for i=2:M-1
    for j=2:N-1
        if F(i,j)==3
            GM_pixels=GM_pixels+1;
        end
    end
end
GM_pixels

WM_pixels=0;
for i=2:M-1
    for j=2:N-1
        if F(i,j)==4
            WM_pixels=WM_pixels+1;
        end
    end
end
WM_pixels


JS_CSF_IFCM=100*(min(2919,CSF_pixels)/max(2919,CSF_pixels));
JS_GM_IFCM=100*(min(7108,GM_pixels)/max(7108,GM_pixels));
JS_WM_IFCM=100*(min(9480,WM_pixels)/max(9480,WM_pixels));


DC_CSF_IFCM=100*2*(min(2919,CSF_pixels)/(2919+CSF_pixels));
DC_GM_IFCM=100*2*(min(7108,GM_pixels)/(7108+GM_pixels));
DC_WM_IFCM=100*2*(min(9480,WM_pixels)/(9480+WM_pixels));
%
 Result=[JS_CSF_IFCM  DC_CSF_IFCM ;JS_GM_IFCM DC_GM_IFCM ;JS_WM_IFCM DC_WM_IFCM ]
% 


