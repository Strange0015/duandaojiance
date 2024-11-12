function data_cleaned = remove_near_duplicates(data, tolerance)
    % ȷ���ݲ������
    if tolerance <= 0
        error('Tolerance must be a positive number.');
    end
 
    % ��ʼ�����������ݺ͵�ǰֵ
    data_cleaned = data(1,:);
    current_value = data(1,1);
 
    for i = 2:length(data)
        if abs(data(i,1) - current_value) > tolerance
            data_cleaned = [data_cleaned; data(i,:)];
            current_value = data(i,1);
        end
    end
end