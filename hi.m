% 定义节点结构
classdef Node
    properties
        points
        children
        parent
        interactionList
        u
        d
    end
    methods
        function obj = Node(points)
            obj.points = points;
            obj.children = [];
            obj.parent = [];
            obj.interactionList = [];
            obj.u = 0;
            obj.d = 0;
        end
    end
end

% 生成示例树
function tree = generateTree()
    % 创建根节点
    root = Node(rand(10, 2));
    % 细分根节点
    root.children = [Node(rand(5, 2)), Node(rand(5, 2))];
    root.children(1).parent = root;
    root.children(2).parent = root;
    tree = root;
end

% 前序遍历和细分节点
function tree = subdivideTree(tree, p)
    if length(tree.points) > p
        tree.children = [Node(rand(5, 2)), Node(rand(5, 2))];
        tree.children(1).parent = tree;
        tree.children(2).parent = tree;
        for i = 1:length(tree.children)
            tree.children(i) = subdivideTree(tree.children(i), p);
        end
    end
end

% 构建 U-list 和 V-list
function [U_list, V_list] = constructLists(tree)
    U_list = {};
    V_list = {};
    if isempty(tree.children)
        U_list{end+1} = tree;
        V_list{end+1} = tree;
    else
        for i = 1:length(tree.children)
            [U, V] = constructLists(tree.children(i));
            U_list = [U_list, U];
            V_list = [V_list, V];
        end
    end
end

% 向上传递
function tree = upwardPass(tree)
    if isempty(tree.children)
        tree.u = sum(tree.points(:));
    else
        for i = 1:length(tree.children)
            tree.children(i) = upwardPass(tree.children(i));
            tree.u = tree.u + tree.children(i).u;
        end
    end
end

% 向下传递
function tree = downwardPass(tree)
    if ~isempty(tree.parent)
        tree.d = tree.parent.d + sum(tree.points(:));
    end
    for i = 1:length(tree.children)
        tree.children(i) = downwardPass(tree.children(i));
    end
end

% 收集结果
function E = gatherResults(tree)
    E = [];
    if isempty(tree.children)
        E = sum(tree.points(:)) + tree.d;
    else
        for i = 1:length(tree.children)
            E = [E; gatherResults(tree.children(i))];
        end
    end
end

% 主函数
function main()
    % 生成树
    tree = generateTree();
    % 细分树
    tree = subdivideTree(tree, 5);
    % 构建 U-list 和 V-list
    [U_list, V_list] = constructLists(tree);
    % 向上传递
    tree = upwardPass(tree);
    % 向下传递
    tree = downwardPass(tree);
    % 收集结果
    E = gatherResults(tree);
    disp('结果:');
    disp(E);
end

% 运行主函数
main();

