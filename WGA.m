clear
for iijj=1:1
    clc
    disp('Wild Geese Algorithm: A novel algorithm for large scale optimization based on the natural life and death of wild geese');
    
    npop=120;
    nvar=1000;
    D=nvar;
    In=npop;
    Fi=30;
    Cr=0.5;
    V0=1;
    
    xmin=-600;
    xmax=-xmin;
    dx=xmax-xmin;
    
    % vmax=0.5*dx;
    % vmin=-vmax;
    
    FEs=5000000;
    VarSize = [1 nvar];
    M=[];
    o=[];
    A1=[];
    %%%%%%%%%%%
    it=1;
    while it<=FEs
        % for it=1:maxit
        if it==1
            gbestcost(1)=inf;
            for i=1:npop
                velocity(i,:)=zeros(1,nvar);
                position(i,:)=xmin+(xmax-xmin)*rand(1,nvar);
                cost(i)=Cost(position(i,:));
                pbest(i,:)=position(i,:);
                pbestcost(i)=cost(i);
                vbest(i,:)=velocity(i,:);
                
                delta(i,:)=pi*rand;
                if pbestcost(i)<gbestcost(it)
                    gbest(1,:)=pbest(i,:);
                    gbestcost(it)=pbestcost(i);
                    vG(1,:)= vbest(i,:);
                end
            end
        else
            gbest(1,:)=gbest(1,:);
            gbestcost(it)=gbestcost(it-npop);
            [hh, gg]=sort(pbestcost);
            
            npop=(In-1)-((In-Fi)*((it)/FEs));
            npop=round(npop+1);
            npop=max(npop,Fi);
            npop=min(In,npop);
            B6=In-Fi;
            
            for eee=1:npop
                
                if B6==0
                    i=eee;
                else
                    
                    i=gg(eee);
                end
                [f1, f2]=find(gg==i);
                
                %%%BEDDER
                if f2==npop
                    f2=0;
                end
                jj1=gg(1,f2+1);
                %%%%BETTER
                [f1, f2]=find(gg==i);
                tt=1;
                if f2==1
                    f2=npop+1;
                    tt=-1;
                end
                jj2=gg(1,f2-1);
                if f2==2
                    f2=npop+2;
                    
                end
                jj3=gg(1,f2-2);
                jjj=gg(1,1);
                ff1=gg(1,end);
                
                
                velocity(i,:)= (rand(1,nvar).*velocity(i,:)+rand(1,nvar).*(velocity(jj2,:)-velocity(jj1,:)))+rand(1,nvar).*(pbest(i,:)-position(jj1,:))+rand(1,nvar).*(pbest(jj2,:)-position(i,:))-rand(1,nvar).*(pbest(jj1,:)-position(jj3,:))+rand(1,nvar).*(pbest(jj3,:)-position(jj2,:));%%ORIGINAL
                
                
                BB=(pbestcost(jj2))/(pbestcost(i));
                GG=(pbestcost(jjj))/(pbestcost(i));
                
                position(i,:)=pbest(i,:)+rand(1,nvar).*rand(1,nvar).*(( pbest(jj2,:)+gbest(1,:)-2*pbest(i,:))+(velocity(i,:)));
                
                f1=(gbestcost(it))/(pbestcost(i)+gbestcost(it));
                f0=(pbestcost(jj2))/(pbestcost(jj2)+pbestcost(i));
                
                
                DE1(1,:)=((pbest(jj2,:)-pbest(i,:)));
                for ww=1:nvar
                    if rand<Cr
                        position(i,ww)=pbest(i,ww)+rand*rand*(DE1(1,ww));
                    end
                end
                DE1(1,:)=[];
                position(i,:)=min(max(position(i,:),xmin),xmax);
                
                cost(i)=Cost(position(i,:));
                
                
                if cost(i)<pbestcost(i)
                    pbest(i,:)=position(i,:);
                    pbestcost(i)=cost(i);
                    vbest(i,:)=velocity(i,:);
                    if pbestcost(i)<gbestcost(it)
                        gbest(1,:)=pbest(i,:);
                        gbestcost(it)=pbestcost(i);
                        vG(1,:)= vbest(i,:);
                    end
                end
                
            end
            
        end
        
        
        disp(['Iteration ' num2str(it) ':   Best Cost = ' num2str(gbestcost(it))]);
        
        it=it+npop;
        
    end
    
    plot(log(gbestcost),'g-.', 'LineWidth',2.5); hold
    hold on;
    Cost_Rsult(1, iijj)=gbestcost(end);
    if iijj==1
        save New1
    end
    if iijj==2
        save New2
    end
    if iijj==3
        save New3
    end
    if iijj==4
        save New4
    end
    if iijj==5
        save New5
    end
    if iijj==6
        save New6
    end
    if iijj==7
        save New7
    end
    if iijj==8
        save New8
    end
    if iijj==9
        save New9
    end
    if iijj==10
        save New10
    end
end
Mean=mean(Cost_Rsult)
Std=std(Cost_Rsult)
