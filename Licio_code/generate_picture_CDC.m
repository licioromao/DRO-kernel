% Script to generate plot for the CDC paper

load('CDC_2.mat')
figure 
hold on
box on
grid on
load(results{1}.paths_results(2).full_path)
plot(state_partition.get_values.partition.grid_x,value_func_kernel_KME.value_function(:,1),'--b','LineWidth',1.5) % radius equal to zero
load(results{1}.paths_results(3).full_path)
plot(state_partition.get_values.partition.grid_x,value_func_kernel_KME.value_function(:,1),'--*b','LineWidth',1.5) % radius equal to 0.01
load(results{1}.paths_results(4).full_path)
plot(state_partition.get_values.partition.grid_x,value_func_kernel_KME.value_function(:,1),'-.b','LineWidth',1.5) % radius equal to 0.05
load(results{1}.paths_results(5).full_path)
plot(state_partition.get_values.partition.grid_x,value_func_moment.value_function(:,1),'r','LineWidth',2) % radius mean = 0.2, radius variance = 1.5

ll = legend({'Kernel ambiguity $(\\\epsilon = 0)$','Kernel ambiguity ($\\\epsilon = 0.01$)',...
                'Kernel ambiguity ($\\\epsilon$ = 0.05)','Moment ambiguity ($b = 0.1, c = 1.5$, see [28])'},...
                    'Interpreter','latex','FontSize',14,'Orientation','vertical','Position',[0.57,0.20,0.1,0.1]);
set(ll,'color','none');

title_text = {'Comparison of safe probability'};

title(title_text,'Interpreter','latex','FontSize',18)

title_x = {'Initial temperature'};
xlabel(title_x,'Interpreter','latex','FontSize',18)

title_y = {'Safety probability'};
ylabel(title_y,'Interpreter','latex','FontSize',18)


path_to = '/Users/licioromao/OneDrive - Nexus365/Postdoc/Papers/Ashish_collaboration/CDC-paper/ACC23_AshishLicio_KernelSafety/';

kernel_versus_moment = gcf;
temp = strcat(path_to,'kernel-versus-moment.pdf');
exportgraphics(kernel_versus_moment,temp)
