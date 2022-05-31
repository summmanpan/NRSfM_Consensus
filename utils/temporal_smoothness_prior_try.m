
% Order varying temporal regularization
frames = 9;

%% FIRST-ORDER
%only 2 frames are consider to impulse the temporal constraint, 
% 2 entries per columns
L1 = -eye(frames,frames-1); % n_frames x n_frames-1 %eliminate last colum
L1(2:frames+1:end)=1
L1_s = sparse(L1);
% spy(L1_s)

%% SECOND-ORDER
L2 = eye(frames)*2;
L2(2:frames+1:end)=-1;
L2(frames+1:frames+1:end)=-1;
% add boundary conditions to the 1st and last entries
L2(1,1)=1;
L2(frames,frames)=1;

%% FOURTH-ORDER
L4 = eye(frames)*-30;

% L4(2:frames+1:end)= 16;
% L4(3:frames+1:end)= -1;
% L4(frames+1:frames+1:end)= 16;
% L4((frames*2)+1:frames+1:end)= -1;
% add boundary conditions to the 1sts and lasts entries
variable = 0;
for i= 1:2
    if i~=1
        variable = -16;
    end
    L4(i+1:frames+1:end)= 16/(2-i+variable); 
    L4((frames*i)+1:frames+1:end)= 16/(2-i+variable);

    % boundary conditions
    L4(i,i)=variable;
    L4(frames-(i-1),frames-(i-1))=variable;
    L4(i,3-i)=1;
    L4(frames-(i-1),frames-(2-i))=1;
end
L4(1,frames)=0;
% L4(1,1)=0;
% L4(frames,frames)=0;
% L4(2,2)=-16;
% L4(frames-1,frames-1)=-16;
% L4(2,1)=1;L4(1,2)=1;
% L4(frames,frames-1)=1;
% L4(frames-1,frames)=1;

%%
%% FOURTH-ORDER
% L4 = eye(frames)*-30;

% L4(2:frames+1:end)= 16;
% L4(3:frames+1:end)= -1;
% L4(frames+1:frames+1:end)= 16;
% L4((frames*2)+1:frames+1:end)= -1;
% add boundary conditions to the 1sts and lasts entries
% L4(1,1)=0;
% L4(frames,frames)=0;
% L4(2,2)=-16;
% L4(frames-1,frames-1)=-16;
% L4(2,1)=1;L4(1,2)=1;
% L4(frames,frames-1)=1;
% L4(frames-1,frames)=1;


%%
Ll = get_L(5,4)

function L = get_L(frames,order)
% return a F × F matrix encoding temporal smoothness priors.
    switch(order)
        case 1
            % FIRST-ORDER
            % only 2 frames are consider to impulse the temporal constraint, 
            % 2 entries per columns
            L = -eye(frames,frames-1); % n_frames x n_frames-1 %eliminate last colum
            L(2:frames+1:end)=1;

        case 2
            %SECOND-ORDER
            L = eye(frames)*2;
            L(2:frames+1:end)=-1;
            L(frames+1:frames+1:end)=-1;
            % add boundary conditions to the 1st and last entries
            L(1,1)=1;
            L(frames,frames)=1;
        case 4
            % FOURTH-ORDER
            L = eye(frames)*-30;
            variable = 0;
            for i= 1:2
                if i~=1
                    variable = -16;
                end
                L(i+1:frames+1:end)= 16/(2-i+variable); 
                L((frames*i)+1:frames+1:end)= 16/(2-i+variable);
                % boundary conditions
                L(i,i)=variable;
                L(frames-(i-1),frames-(i-1))=variable;
                L(i,3-i)=1;
                L(frames-(i-1),frames-(2-i))=1; 
            end
            L(1,frames)=0; 
        otherwise
            disp('no exist')
            L=0;
    end
end