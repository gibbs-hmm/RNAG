function y=hier_clus_octave(Dir_data,num_seq,burn_in,iter,indicator_ref, gamma_choice)
%indicator_ref = 0 if there is no reference given, and = 1 if there is a
%reference

%gamma_choice = 0 for default gamma = 1, gamma_choice = 1 to get all output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mainfolder = pwd;
if(Dir_data(end)~='/')
    Dir_data=[Dir_data,'/'];
end
eval(['cd ' Dir_data])
% set(0,'RecursionLimit',1000000)  % comment this out for ocatve

ppv = zeros(1,16);
sen = zeros(1,16);

%overwrite any existing files

if indicator_ref == 1
    dlmwrite('bias_variance','')
else
    dlmwrite('variance','')
end
if gamma_choice == 1 && indicator_ref == 1
    dlmwrite('PPV_SEN','')
end

for ind_seq=1:num_seq

    dlmwrite(['centroid_',num2str(ind_seq)],'')
    dlmwrite(['en_centroid_',num2str(ind_seq)],'')
    
    if indicator_ref == 1
        truth=textread(['tru_str_',num2str(ind_seq)],'%s');
        truth=truth{1};
        %N=length(truth);
    end
    %1)transform from bracket to binary matrix
    tmp_st=read_project(['project_',num2str(ind_seq),'.str']);
    N = length(tmp_st{3});
    stru=cell(1,1+iter+burn_in);
    if indicator_ref == 1
        stru{1}=bra2list(truth);
    else
        stru{1} = [];
    end
    for i=1:iter+burn_in
        stru{i+1}=bra2list(tmp_st{3*i});
    end
    iteration=length(tmp_st)/3;
    data=iteration-burn_in;



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    %2)Hierarchical Clustering
    %burn in the first burn_in
    dis2=find_same(stru(burn_in+2:end),1);
    Z = linkage(squareform(dis2,'tovector'),'average');  % comment for ocatve
    %c=dendrogram(Z);
    %title('Linkage average')
    %saveas(gcf,['tree',num2str(ind_seq),'.jpg']); 
    %close;
    %dstnme=sprintf('stu_%d.mat',ind_seq);  
    %eval(['save ',dstnme,' stru dis2 Z dis']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    %3)CH-index to determine number of clusters 2:20
    %dstnme=sprintf('stu_%d.mat',ind_seq);  
    %eval(['load ',dstnme]);
    best=-1;
    for i=2:20
	tmp_cl= cluster(Z,'maxclust',i); 
%        tmp_cl = myKmeans(dis2, i);  % for Octave
        enen=zeros(1,N^2);
        pp=cell(1,i);
        W=0;
        for ii=1:i
            centroid=zeros(1,N^2);
            cl_ind=find(tmp_cl==ii);
            len=length(cl_ind);
            for iii=cl_ind'
                centroid(stru{burn_in+1+iii})=centroid(stru{burn_in+1+iii})+1;
            end
            enen=enen+centroid;
            centroid=centroid/len;
            [str,pair]=find_cen(centroid,N,1);
            pp{ii}=pair';
            tmp=find_same({pair',stru{burn_in+1+cl_ind}},2);
            W=W+sum(tmp.^2);
        end
        len=data;
        centroid=enen/len;
        [str,pair]=find_cen(centroid,N,1);
        tmp=find_same([pair',pp],2);
        B=sum(tmp.^2);
        tmp_best=B*(data-i)/((i-1)*W);
        if(best<tmp_best)
            cl=tmp_cl;
            cl_no=i;
            best=tmp_best;
        end
    end
    
    
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %4)multiple centroid and 95 credibility limit
    %sort it a little bit
    T=zeros(1,cl_no);
    for i=1:cl_no
        T(i)=sum(cl==i);
    end
    [ord,oo]=sort(T,'descend');
    tmp_cl=zeros(1,data);
    for i=1:cl_no
        tmp_cl(cl==oo(i))=i;
    end
    cl=tmp_cl;
    enen=zeros(1,N^2);
    for i=1:cl_no
        centroid=zeros(1,N^2);
        cl_ind=find(tmp_cl==i);
        len=length(cl_ind);
        for iii=cl_ind
            centroid(stru{burn_in+1+iii})=centroid(stru{burn_in+1+iii})+1;
        end
        enen=enen+centroid;
        centroid=centroid/len;
        dlmwrite(['centroid_',num2str(ind_seq)],['>>centroid ',num2str(i),' :(',num2str(len),')'],'delimiter','','-append')
        if gamma_choice == 1
            for gamma=[2.^[-5:1:10],6]
                [str,pair]=find_cen(centroid,N,gamma);
                tmp_dis=sort(find_same({pair',stru{burn_in+1+cl_ind}},2),'ascend');
                wrt=['g=',num2str(gamma),',CL_95:',num2str(tmp_dis(ceil(.95*len)))];
                dlmwrite(['centroid_',num2str(ind_seq)],[str,wrt],'delimiter','','-append')
            end
        else
            [str,pair] = find_cen(centroid,N,1);
            tmp_dis = sort(find_same({pair',stru{burn_in+1+cl_ind}},2),'ascend');
            wrt=['CL_95:',num2str(tmp_dis(ceil(.95*len)))];
            dlmwrite(['centroid_',num2str(ind_seq)],[str,wrt],'delimiter','','-append')
        end

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5-
    %5)ensemble centroid and 95 credibility limit
    
    len=data;
    centroid=enen/len;
    dlmwrite(['en_centroid_',num2str(ind_seq)],['>>centroid ',num2str(1),' :(',num2str(len),')'],'delimiter','','-append')
    if gamma_choice == 1
        gamma_en_centroids = cell(1,16);
        gammalist = [2.^[-5:1:2],6,2.^[3:1:10]];
        for gammaind = 1:16
            [str,pair]=find_cen(centroid,N,gammalist(1,gammaind));
            gamma_en_centroids{gammaind} = pair;
            tmp_dis=sort(find_same({pair',stru{burn_in+1+(1:len)}},2),'ascend');
            wrt=['g=',num2str(gammalist(1,gammaind)),',CL_95:',num2str(tmp_dis(ceil(.95*len)))];
            dlmwrite(['en_centroid_',num2str(ind_seq)],[str,wrt],'delimiter','','-append')
        end
    else
        [str,pair] = find_cen(centroid,N,1);
        tmp_dis = sort(find_same({pair',stru{burn_in+1+(1:len)}},2),'ascend');
        wrt = ['CL_95:',num2str(tmp_dis(ceil(.95*len)))];
        dlmwrite(['en_centroid_',num2str(ind_seq)],[str,wrt],'delimiter','','-append')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %6)bias and variance
    if indicator_ref == 1
    
        bias_dis=find_same(stru,2);
        bias = mean(bias_dis(burn_in+1:end));
        [str,pair] = find_cen(centroid,N,1);
        var_dis = find_same({pair',stru{burn_in+1:end}},2);
        var_sq = var_dis.^2;
        variance = mean(var_sq);
        %bias-variance
        dlmwrite('bias_variance',['>> sequence ', num2str(ind_seq), ': '],'delimiter','','-append')
        dlmwrite('bias_variance',['Bias: ',num2str(bias)],'delimiter','','-append')
        dlmwrite('bias_variance',['Variance: ',num2str(variance)],'delimiter','','-append')
    else
        [str,pair] = find_cen(centroid,N,1);
        var_dis = find_same({pair',stru{burn_in+1:end}},2);
        var_sq = var_dis.^2;
        variance = mean(var_sq);
        %variance
        dlmwrite('variance',['>> sequence ', num2str(ind_seq), ': '],'delimiter','','-append')
        dlmwrite('variance',['Variance: ',num2str(variance)],'delimiter','','-append')
    end
    
    %7) PPV-SEN
    if indicator_ref == 1 && gamma_choice == 1
        gammalist = [2.^[-5:1:2],6,2.^[3:1:10]];
        for gammaind = 1:16
            [tp,fp,fn] = comparestructures(gamma_en_centroids{gammaind},stru{1});
            ppv(1,gammaind) = ppv(1,gammaind) + tp/(tp+fp);
            sen(1,gammaind) = sen(1,gammaind) + tp/(tp+fn);
        end
    end
        

end
if indicator_ref == 1 && gamma_choice == 1
    ppv = ppv./num_seq;
    sen = sen./num_seq;
    plot(ppv,sen,'*-')
    xlabel('PPV')
    ylabel('SEN')
    title('Positive Predictive Value vs Sensitivity')
    print('-dpdf', 'PPV_SEN', '-r150');
    gammalist = [2.^[-5:1:2],6,2.^[3:1:10]];
    for gammaind = 1:16
        dlmwrite('PPV_SEN',['gamma=',num2str(gammalist(1,gammaind)),': PPV=',num2str(ppv(1,gammaind)),', SEN=',num2str(sen(1,gammaind))],'delimiter','','-append')
    end
end
eval(['cd ',mainfolder])
exit()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ch_ind=ch(num_cl,diss,cluster,T)
%total sum of squares T
%diss=dis2.^2;
%T=sum(sum(diss));
%sum of squares within each cluster W
    WW=0;
    for k=1:num_cl
        ind=find(cluster==k);
        WW=WW+sum(sum(diss(ind,ind)));
    end
    BB=T-WW;
    ch_ind=BB*(size(diss,1)-num_cl)/((num_cl-1)*WW);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=bra2list(bracket)
    topstack = 0;
    openstack = [];
    y=[];
    NN=length(bracket);
    for k=1:NN
        if (bracket(k)=='(')
            topstack = topstack + 1;
            openstack(topstack) = k; 
        elseif (bracket(k)==')')
            y=[y,openstack(topstack)*NN+k];
            topstack = topstack-1;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find_same makes hamming distance matrices (1 = pairwise, 2 = wrt
%reference)
function y=find_same(stru,mode)
    if(mode==1)
    %dist matrix
    y=zeros(length(stru));
    for k=1:length(stru)
        for j=k+1:length(stru)        
            y(k,j)=dist_ham(stru{k},stru{j});
        end
    end
    y=y+y';
    elseif(mode==2)
        y=zeros(1,length(stru)-1);
        for k=2:length(stru)
            y(k-1)=dist_ham(stru{1},stru{k});
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=dist_ham(t1,t2)
    he=unique([t1,t2]);
    y=2*length(he)-length(t1)-length(t2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [str,pair]=find_cen(pij,N,gamma)
    % gamma~[0,1]:direct
    if(gamma<=1)
        str=repmat('.',1,N);
        pair=[];
        for k=find(pij>1/(1+gamma))
            po1=floor(k/N);
            po2=k-N*po1;
            if(po2==0)
                po1=po1-1;
                po2=N;
            end
            str(po1)='(';
            str(po2)=')';
            pair=[pair;k];
        end

    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gamma~[1,...]:DP
        %1) Recursion:
        DP=zeros(N);
        for L=2:N
            for k=1:N+1-L
                j=k+L-1;
                DP(k,j)=max([DP(k+1,j),DP(k,j-1),DP(k+1,j-1)+(1+gamma)*pij(k*N+j)-1,DP(k,[k+1:j-1])+DP(1+[k+1:j-1],j)']);
            end
        end
        %2) Backtrack:
        stack=[1,N];
        pair=[];
        while ~isempty(stack)
            k=stack(1);
            j=stack(2);
            stack=stack(3:end);
            if(k<j)
                if(DP(k+1,j)==DP(k,j))
                    stack=[stack,k+1,j];
                elseif(DP(k,j-1)==DP(k,j))
                    stack=[stack,k,j-1];
                elseif(DP(k+1,j-1)+(1+gamma)*pij(k*N+j)-1==DP(k,j))
                    stack=[stack,k+1,j-1];
                    pair=[pair;[k,j]];
                else
                    for kk=k+1:j-1
                        if(DP(k,kk)+DP(1+kk,j)==DP(k,j))
                            stack=[stack,kk+1,j,k,kk];
                            break;
                        end
                    end
                end
            end
        end
        %3) return:
        str=repmat('.',1,N);
        if(~isempty(pair))
            str(pair(:,1))='(';
            str(pair(:,2))=')';
            pair=pair(:,1)*N+pair(:,2);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [TP,FP,FN] = comparestructures(test,ref)
    TP = 0;
    FP = 0;
    FN = 0;
    for j = 1:length(test)
        tempnum = test(j);
        found = 0;
        k = 1;
        while (found == 0 && k <= length(ref))
            if ref(k) == tempnum
                found = 1;
            end
            k = k+1;
        end
        if found == 1
            TP = TP + 1;
        end
    end
    FP = length(test) - TP;
    FN = length(ref) - TP;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read a project file to match Matlab input
% needed because Octave textread differs from Matlab and differs among
% octave versions
function lines = read_project(filename)
    fid = fopen(filename, 'r');
    lines = {};
    count = 1;
    structure_count = 1;
    while ~feof(fid)
      line = fgetl(fid);
      lines{count} = '>Structure';
      lines{count + 1} = num2str(count);
      line = fgetl(fid);
      lines{count + 2} = line ;
      count = count + 3;
      structure_count = structure_count + 1;
    end
    fclose(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
