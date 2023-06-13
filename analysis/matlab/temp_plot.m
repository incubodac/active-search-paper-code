load(    '/home/juank/Desktop/EEGEYE_EJN2022/Unfold2022/sac_models/Inter+face/Results.mat')

t = unfoldResults.E01.times;
ch = 28;
FFNN = fieldnames(unfoldResults);

fea = 1; I1 = []; for i=1:length(FFNN); I1 = [I1; squeeze(unfoldResults.(FFNN{i}).beta_dc(ch,:,fea))]; end
fea = 3; I2 = []; for i=1:length(FFNN); I2 = [I2; squeeze(unfoldResults.(FFNN{i}).beta_dc(ch,:,fea))]; end

figure; 
    hold on
        plot(t,mean(I1),'r-','linewidth',2)
        plot(t,mean(I2),'b-','linewidth',2)
    hold off
    length('Intercept 1','Intercept 2')
