function all_data = sample_experiment(ndata,experiment,salt,replicate,gene_id)
        if length(gene_id) == 1
            % load second experimental replicate data
            if salt==0.2
                load('../mat_files/Analyzed_FISH_data_2016_06.mat')
                 data = RNAFISHexp(1).gene(gene_id).rep(replicate).RNAhisttot';
              

            elseif salt==0.4
                load('../mat_files/Analyzed_FISH_data_2016_06.mat')
                data = RNAFISHexp(2).gene(gene_id).rep(replicate).RNAhisttot';
                if replicate == 3
                    data = [ones(1,size(data,2)); data]; 
                end
                        

            end

            all_data = zeros(length(experiment.tvec),size(data,2),ndata);
            size(data)
            for k=1:ndata         
                for i=1:length(experiment.tvec)
                    all_data(i,:,k) = sample_d(data(i,:),experiment.Nc(i));
                end
            end
        else
            % all_data = zeros(length(experiment.tvec),size(data,2),ndata);
            for k=1:ndata         
                for i=1:length(experiment.tvec)
                    if salt==0.2
                        load('../mat_files/Analyzed_FISH_data_2016_06.mat')
                        if replicate==1
                            data_i = RNAFISHexp(1).gene(gene_id(i)).rep(1).RNAhisttot';
                        elseif replicate==2
                            data_i = RNAFISHexp(1).gene(gene_id(i)).rep(2).RNAhisttot';
                        end
                    elseif salt==0.4
                        load('../mat_files/Analyzed_FISH_data_2016_06.mat')
                        if replicate==1
                            data_i = RNAFISHexp(2).gene(gene_id(i)).rep(1).RNAhisttot';
                        elseif replicate==2
                            data_i = RNAFISHexp(2).gene(gene_id(i)).rep(2).RNAhisttot';
                        end
                    end
                    all_data(i,:,k) = sample_d(data_i(i,:),experiment.Nc(i));
                end
            end
        end
end