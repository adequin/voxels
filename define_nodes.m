%define nodes

function n = define_nodes(varargin)
    % inputs: vox_id, or lattice, row,col,slice
    % output: node struct that contains all nodes per face
    %         and includes each faces normal, up, and right vectors. To access
    %         nodes use this table 

    %define UV table
    % XYZ: X points right, Y points in, and Z points up
    % FACE NORMAL,    UP,     RIGHT,    NUMBER
    % -Y              +Z      +X        1    
    % +X              +Z      +Y        2
    % +Y              +Z      -X        3
    % -X              +Z      -Y        4
    % +Z              +Y      +X        5
    % -Z              -Y      -X        6

    % then access using n.("f" + NUMBER).(direction from face center + left, middle, or right)
    % e.g. n.f1.lb -> node on face 1 left bottom relative to the uv table

    narginchk(1,4) % enforce that only 1 or 3 inputs are valid
    if nargin == 1
        vox_id = varargin{1};
        % do your "from id" logic here

        
    elseif nargin == 4
        lat = varargin{1};
        r = varargin{2};
        c = varargin{3};
        s = varargin{4};
        vox_id = lat.id(r,c,s);
    end


        
    %define canonical nodes
    n = struct();
    off = (vox_id-1)*42;

    for i = 1:12
        for j = 1:3
            n.("e"+i+j) = j + (i-1)*3 + off; %ei is the edge and j is the slot (1->3) in the direction of the edge point pair
        end
    end
    for k = 1:6
        n.("fn"+k) = 36 + k + off; %add 6 face nodes at the end
    end

    n.f1.normal = [0;-1;0];
    n.f2.normal = [1;0;0];
    n.f3.normal = [0;1;0];
    n.f4.normal = [-1;0;0];
    n.f5.normal = [0;0;1];
    n.f6.normal = [0;0;-1];

    n.f1.right = [1;0;0];
    n.f2.right = [0;1;0];
    n.f3.right = [-1;0;0];
    n.f4.right = [0;-1;0];
    n.f5.right = [1;0;0];
    n.f6.right = [-1;0;0];

    n.f1.up = [0;0;1];
    n.f2.up = [0;0;1];
    n.f3.up = [0;0;1];
    n.f4.up = [0;0;1];
    n.f5.up = [0;1;0];
    n.f6.up = [0;-1;0];

    face_nodes = ["bl","bm","br", ...   % bottom edge (left→mid→right)
                  "rb","rm","rt", ...   % right edge (bottom→mid→top)
                  "tr","tm","tl", ...   % top edge (right→mid→left)
                  "lt","lm","lb", ...   % left edge (top→mid→bottom)
                  "c"];                 % center

    f1 = [ n.e11, n.e12, n.e13, ...   % bottom (E1: 1–2, left→mid→right)
           n.e61, n.e62, n.e63, ...   % right  (E6: 2–6, bottom→mid→top)
           n.e93, n.e92, n.e91, ...   % top    (E9: 5–6, right→mid→left)
           n.e53, n.e52, n.e51, ...   % left   (E5: 1–5, top→mid→bottom)
           n.fn1 ];

    f2 = [ n.e21, n.e22, n.e23, ...
           n.e71, n.e72, n.e73, ...
           n.e103, n.e102, n.e101, ...
           n.e63, n.e62, n.e61, ...
           n.fn2 ];

    f3 = [ n.e31, n.e32, n.e33, ...
           n.e81, n.e82, n.e83, ...
           n.e113, n.e112, n.e111, ...
           n.e73, n.e72, n.e71, ...
           n.fn3 ];

    f4 = [ n.e41, n.e42, n.e43, ...
           n.e51, n.e52, n.e53, ...
           n.e123, n.e122, n.e121, ...
           n.e83, n.e82, n.e81, ...
           n.fn4 ];

    f5 = [ n.e91,  n.e92,  n.e93, ...
           n.e101, n.e102, n.e103, ...
           n.e111, n.e112, n.e113, ...
           n.e121, n.e122, n.e123, ...
           n.fn5 ]; %all increasing

    f6 = [ n.e33, n.e32, n.e31, ...
           n.e23, n.e22, n.e21, ...
           n.e13, n.e12, n.e11, ...
           n.e43, n.e42, n.e41, ...
           n.fn6 ]; %all flipped

    n.face_names = ["f1", "f2", "f3", "f4", "f5", "f6"];
    n.face_nodes = face_nodes;
    face_vals = [f1; f2; f3; f4; f5; f6];
    for i = 1:numel(n.face_names)
        fname = char(n.face_names(i));
        for j = 1:numel(face_nodes)
            alias = char(face_nodes(j));
            n.(fname).(alias) = face_vals(i,j);
        end
    end

    assert(n.f1.tl == n.f5.bl);
    assert(n.f1.tm == n.f5.bm);
    
end 

