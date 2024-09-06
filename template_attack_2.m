function template_attack_2()
succ=0;
num_of_pro=160000;
num_of_poi=8;
num_of_attack=100;
SNR=0.05;
for repeat=1:1000
    state=randi(256,num_of_pro,1)-1;
    %signal part
    weight=randi(40,num_of_poi,8);
    signal=zeros(num_of_pro,num_of_poi);
    for i=1:num_of_pro
        for j=1:num_of_poi
            for k=1:8
                signal(i,j)=signal(i,j)+bitget(state(i,1),k)*weight(j,k);
            end
        end
    end
    var_signal=var(signal);
    %noise part
    var_noise=var_signal/SNR;
    rho=0.4;
    mean_vector=zeros(1,num_of_poi);
    covariance_matrix=zeros(num_of_poi,num_of_poi);
    for i=1:num_of_poi
        for j=1:num_of_poi
            if i==j
                covariance_matrix(i,j)=var_noise(1,i);
            else
                covariance_matrix(i,j)=rho*sqrt(var_noise(1,i)*var_noise(1,j));
            end
        end
    end
    noise=mvnrnd(mean_vector,covariance_matrix,num_of_pro);
    %profile trace
    profile_leakage=zeros(num_of_pro,num_of_poi);
    for i=1:num_of_pro
        for j=1:num_of_poi
            profile_leakage(i,j)=signal(i,j)+noise(i,j);
        end
    end
    %characterization
    mean_vector_c=zeros(256,num_of_poi);
    covariance_matrix_c=zeros(num_of_poi,num_of_poi);
    num_of_group=zeros(256,1);
    leakage_of_group=zeros(256,num_of_poi);
    for i=1:num_of_pro
        leakage_of_group(state(i,1)+1,:)=leakage_of_group(state(i,1)+1,:)+profile_leakage(i,:);
        num_of_group(state(i,1)+1,1)=num_of_group(state(i,1)+1,1)+1;
    end
    for i=1:256
        if num_of_group(i,1)~=0
            mean_vector_c(i,:)=leakage_of_group(i,:)/num_of_group(i,1);
        end
    end
    noise=zeros(num_of_pro,num_of_poi);
    for i=1:num_of_pro
        noise(i,:)=profile_leakage(i,:)-mean_vector_c(state(i,1)+1,:);
    end
    for i=1:num_of_poi
        for j=1:num_of_poi
            temp=cov(noise(:,i),noise(:,j));
            covariance_matrix_c(i,j)=temp(1,2);
        end
    end
    %
    state=randi(256,1)-1;
    %signal part
    signal=zeros(num_of_attack,num_of_poi);
    for i=1:num_of_attack
        for j=1:num_of_poi
            for k=1:8
                signal(i,j)=signal(i,j)+bitget(state,k)*weight(j,k);
            end
        end
    end
    %noise part
    noise=mvnrnd(mean_vector,covariance_matrix,num_of_attack);
    %attack trace
    attack_leakage=zeros(num_of_attack,num_of_poi);
    for i=1:num_of_attack
        for j=1:num_of_poi
            attack_leakage(i,j)=signal(i,j)+noise(i,j);
        end
    end
    %attack
    result=zeros(1,256);
    for guess=0:255
        for j=1:num_of_attack
            result(1,guess+1)=result(1,guess+1)+(attack_leakage(j,:)-mean_vector_c(guess+1,:))/covariance_matrix_c*(attack_leakage(j,:)-mean_vector_c(guess+1,:))';
        end
    end
    [~,b]=sort(result,'ascend');
    if b(1,1)==state+1
        succ=succ+1;
    end
    [repeat,succ]
end
end