function [s] = dbn_wrapper(args)

    % parse args
    inPath = args.inPath;
    outPath = args.outPath;

    % this object needs to contain all the useful stuff for the DBN code
    % D, max_indegree, prior_graph, lambas, reg_mode, stdise, silent
    load(inPath);

    fprintf('##### entering network inference routine..\n');
    [e i l] = dynamic_network_inference(D, max_indegree, prior_graph, lambdas, reg_mode, stdise, silent);
    fprintf('##### network inference completed.\n');
    save(outPath, 'e', 'i', 'l');

    s = 1;
end
