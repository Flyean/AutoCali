%% 预定义变量
clear
c = 3e8;
fc = 5.31e9;                          %  center frequency of the channel
d = 0.026;                            %  antenna interval
f_space = 312.5e3;                    %  frequency interval
fc_space = 4*f_space;                 %  frequency interval of logged CSI
sub_carrier_index = [-58:4:58];
K = 30;
N = 3;
fc_wifi = zeros(K,1);
% num_packets = 1000;

for q = 1:K
    fc_wifi(q,1) = fc+fc_space*sub_carrier_index(q);
end

d_grid = (-400:40:400)*1e-9;
a_grid = -90:0.2:90;
s_grid = zeros(N*K,length(d_grid)*length(a_grid));

aoa_grid = exp(-1i*2*pi*fc*(0:2).'*d*sind(a_grid)/c);
tof_grid = exp(-1i*2*pi*(0:29).'*d_grid*fc_space);
s_grid = kron(tof_grid,aoa_grid);

% 已知位置B408
% x1 = 2; y1 = -1.1;
% x2 = 6; y2 = 3;
% x3 = 2; y3 = 6.3;
% 已知位置B401
x1 = 2; y1 = -0.6;
x2 = 4.6; y2 = 3.4;
x3 = 2; y3 = 6.6;
% 计算目标位置的概率分布
[X,Y] = meshgrid(-1.7:0.1:4.6, -0.6:0.1:6.6);
P = zeros(size(X));

Pi_offset = exp(-1i * pi * ones(1, 30));
%% 读取相位偏移
for i = 1:3
    devicename = num2str(i+3);
    loadpath = strcat('./phaseoffset_ref/sRE', devicename, '/000/*.mat');
    namelist = dir(loadpath);

    % 读取后namelist 的格式为
    % name -- filename
    % date -- modification date
    % bytes -- number of bytes allocated to the file
    % isdir -- 1 if name is a directory and 0 if not

    len = length(namelist);
    for ii = 1:len
        file_name{ii}=namelist(ii).name;
        x{ii} = load(strcat('./phaseoffset_ref/sRE', devicename, '/000/' , file_name{ii}));
    end
    mean_phaseshift12 = x{1}.mean_phaseshift12;
    mean_phaseshift13 = x{2}.mean_phaseshift13;
    phase_shift_12_neg(i,1,:) = mean_phaseshift12;
    phase_shift_12_neg(i,2,:) = mean_phaseshift12 + pi;
    phase_shift_13_neg(i,1,:) = mean_phaseshift13;
    phase_shift_13_neg(i,2,:) = mean_phaseshift13 + pi;
end
phase_shift_12_neg = exp(-1i*phase_shift_12_neg);
phase_shift_13_neg = exp(-1i*phase_shift_13_neg);

%% 读取数据
for i = 1:3
    data_path = strcat('./data/WiFi_Localization_Dataset_v0/B401/B401_sRE', num2str(i+4), '_user4_wo.mat');
    load(data_path);
    num_packets = size(features_csi, 1);
    expr = ['CSI_sRE' num2str(i+4) '=repmat(features_csi,4,1);'];
    expr1 = ['CSI_sRE' num2str(i+4) '(num_packets+1:2*num_packets, 31:60)=features_csi(:, 31:60).* Pi_offset;'];
    expr2 = ['CSI_sRE' num2str(i+4) '(2*num_packets+1:3*num_packets, 61:90)=features_csi(:, 61:90) .* Pi_offset;'];
    expr3 = ['CSI_sRE' num2str(i+4) '(3*num_packets+1:4*num_packets, 31:60)=features_csi(:, 31:60).* Pi_offset;'];
    expr4 = ['CSI_sRE' num2str(i+4) '(3*num_packets+1:4*num_packets, 61:90)=features_csi(:, 61:90) .* Pi_offset;'];
    eval([expr expr1 expr2 expr3 expr4]);
end


%% 处理数据
for i = 1: num_packets
    for j = 1:3
        for k = 1:4
            switch j
                case 1
                    rcscsi = CSI_sRE5((k-1)*num_packets + i, :);
                case 2
                    rcscsi = CSI_sRE6((k-1)*num_packets + i, :);
                case 3
                    rcscsi = CSI_sRE7((k-1)*num_packets + i, :);
            end
            rcscsi = reshape(rcscsi,[30, 3]).';
            rcscsi = reshape(rcscsi,[90, 1]);
            s_hat = s_grid'*rcscsi/(N*K);
            P = sum(s_hat.*conj(s_hat), 2)/num_packets;

            for jj  = 1:1
                tmp       = repmat(P.',N*K, 1);
                Ap        = s_grid.*tmp;
                R         = Ap*s_grid'+eye(N*K)*1e-8;
                clear tmp s_hat;
                R_inv     = inv(R);
                R_invA    = R_inv*s_grid;
                part1     = R_invA'*rcscsi;
                part2     = sum(conj(R_invA).*s_grid, 1).';
                s_hat     = part1./part2;
                P         = sum(s_hat.*conj(s_hat),2)/num_packets;
            end
            PP = reshape(P,[length(a_grid),length(d_grid)]);
            [row,col] = find(PP==max(max(PP)));
            tof = (-400+(col-1)*40)*1e-9;
            aoa = -90+(row-1)*0.2;
            pdf_col = PP(:, col);
            pdf(i, j, k, :) = pdf_col./sum(pdf_col);
        end
    end
    fprintf('\b\b\b\b\b\b\b\b\b\b%d/%d', i, num_packets);
end

%%
for jj = 1:num_packets
    tag = 1;
    for i = 1:4
        for j = 1:4
            for k = 1:4
                P_loc = zeros(size(X));
                for ii = 1:numel(X)
                    x = X(ii);
                    y = Y(ii);
%                     B408
%                     theta1 = atan2d(x-x1,y-y1);
%                     theta2 = atan2d(y-y2,x2-x);
%                     theta3 = atan2d(x3-x,y3-y);
                    theta1 = atan2d(x-x1,y-y1);
                    theta2 = atan2d(y-y2,x2-x);
                    theta3 = atan2d(x-x3,y3-y);
                    [id1, ~] = find_nearest_angle(theta1);
                    [id2, ~] = find_nearest_angle(theta2);
                    [id3, ~] = find_nearest_angle(theta3);
                    p1 = pdf(jj, 1, i, id1);
                    p2 = pdf(jj, 2, j, id2);
                    p3 = pdf(jj, 3, k, id3);
                    P_loc(ii) = p1*p2*p3;
                end
                maximum(tag) = max(max(P_loc));
%                 if tag == 1
%                     % 绘制概率分布图
%                     figure;
%                     r = [-1.7 6];
%                     v = [-1.1 6.3];
%                     imagesc('XData',r,'YData',v,'CData',abs(P_loc));
%                     xlabel('X');
%                     ylabel('Y');
% %                     colormap(flipud(hot));
%                     axis([-1.7 6 -1.1 6.3]);
%                     axis off
%                     set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','bold','LineWidth',1.5)
%                 end
                tag = tag + 1;
            end
        end
    end
    maximum_list(jj, :) = maximum;
    index(jj) = find(maximum==max(maximum));
    fprintf('\b\b\b\b\b\b\b\b\b\b%d/%d', jj, num_packets);
end
