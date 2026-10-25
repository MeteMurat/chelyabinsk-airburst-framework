function report = r1c6e_native_preflight(repositoryRoot, auditRoot)
% R1-C6E R2: native MATLAB code analysis and numeric-predicate unit tests.
% No scientific kernel is executed by this function. Historical sources are
% read-only. Full forward replay remains a separate, subsequently gated step.
repositoryRoot = char(repositoryRoot);
auditRoot = char(auditRoot);
assert(isfolder(repositoryRoot), 'Replay:MissingRepository', 'Repository missing.');
if ~isfolder(auditRoot), mkdir(auditRoot); end
report = struct('status','RUNNING_NATIVE_PREFLIGHT', ...
    'scientific_simulations_executed',false, ...
    'matlab_release',version('-release'), ...
    'code_files_checked',0,'native_analyzer_error_count',0, ...
    'negative_control_detected',false,'numeric_predicate_test_count',0);
try
    manifest = jsondecode(fileread(fullfile(repositoryRoot,'SOURCE_MANIFEST.json')));
    fileList = {};
    relativeNames = {};
    for k = 1:numel(manifest.files)
        rel = char(manifest.files(k).path);
        [~,~,ext] = fileparts(rel);
        if strcmpi(ext,'.m')
            fileList{end+1} = fullfile(repositoryRoot,rel); %#ok<AGROW>
            relativeNames{end+1} = rel; %#ok<AGROW>
        end
    end
    assert(~isempty(fileList), 'Replay:MissingSources', 'No MATLAB source files.');
    fprintf('native_matlab_preflight_started=True\n');
    issues = codeIssues(fileList, CodeAnalyzerConfiguration="factory", IncludeSubfolders=false);
    active = issues.Issues;
    suppressed = issues.SuppressedIssues;
    activeErrors = nnz(strcmpi(string(active.Severity),"error"));
    suppressedErrors = nnz(strcmpi(string(suppressed.Severity),"error"));
    report.native_analyzer_error_count = activeErrors + suppressedErrors;
    report.code_files_checked = numel(fileList);
    report.files = relativeNames;
    report.active_issue_count = height(active);
    report.suppressed_issue_count = height(suppressed);
    if ~isempty(active), writetable(active,fullfile(auditRoot,'MATLAB_CODE_ISSUES.csv')); end
    if ~isempty(suppressed), writetable(suppressed,fullfile(auditRoot,'MATLAB_SUPPRESSED_ISSUES.csv')); end
    if activeErrors > 0
        disp(active(strcmpi(string(active.Severity),"error"),:));
    end
    if suppressedErrors > 0
        disp(suppressed(strcmpi(string(suppressed.Severity),"error"),:));
    end
    assert(report.native_analyzer_error_count == 0, ...
        'Replay:NativeCodeErrors','Native Code Analyzer found an error; forward replay not started.');

    % Negative control: the exact original unmatched expression must be rejected.
    fixture = fullfile(auditRoot,'intentional_syntax_error_fixture.m');
    fid = fopen(fixture,'w','n','UTF-8');
    assert(fid >= 0,'Cannot create analyzer negative control.');
    fprintf(fid,'function ok = intentional_syntax_error_fixture(a,b,atol,rtol)\n');
    fprintf(fid,'bothNan=isnan(a)&isnan(b); finite=isfinite(a)&isfinite(b); d=abs(a-b);\n');
    fprintf(fid,'ok=all(bothNan(:) | (finite(:)&d(:)<=max(atol,rtol*max(abs(a(:)),abs(b(:)))));\n');
    fprintf(fid,'end\n');
    fclose(fid);
    negative = codeIssues(fixture, CodeAnalyzerConfiguration="factory", IncludeSubfolders=false);
    report.negative_control_detected = any(strcmpi(string(negative.Issues.Severity),"error"));
    assert(report.negative_control_detected,'Replay:AnalyzerControl', ...
        'Native analyzer did not reject the original malformed expression.');

    % Exercise the exact comparison block extracted from run_clean_replay.m,
    % rather than a reimplementation which could diverge from the real check.
    text = fileread(fullfile(repositoryRoot,'reproduction','run_clean_replay.m'));
    block = regexp(text, ...
        '(?s)bothNan=isnan\(a\)&isnan\(b\); finite=isfinite\(a\)&isfinite\(b\);.*?ok = all\(bothNan\(:\) \| withinTolerance\(:\)\);', ...
        'match');
    assert(numel(block)==1,'Replay:PredicateIdentity','Comparison block not unique.');
    block = block{1};
    cases = { ...
        [1 2 3], [1 2 3], 0, 0, true; ...
        0, 1e-9, 1e-9, 0, true; ...
        0, 2e-9, 1e-9, 0, false; ...
        100, 100.000001, 0, 1e-7, true; ...
        100, 100.001, 0, 1e-7, false; ...
        [NaN 1], [NaN 1], 0, 0, true; ...
        [NaN 1], [0 1], 1, 1, false; ...
        Inf, Inf, 0, 0, false; ...
        -Inf, -Inf, 0, 0, false; ...
        zeros(0,1), zeros(0,1), 0, 0, true; ...
        [1 2;3 4], [1 2;3 4], 0, 0, true; ...
        -100, -100.000001, 0, 1e-7, true};
    for k=1:size(cases,1)
        a=double(cases{k,1}); b=double(cases{k,2}); %#ok<NASGU>
        atol=cases{k,3}; rtol=cases{k,4}; %#ok<NASGU>
        eval(block); % Evaluates only the hash-verified comparison block, not a model.
        assert(isscalar(ok) && logical(ok)==cases{k,5}, ...
            'Replay:PredicateTest','Numeric predicate unit test %d failed.',k);
        report.numeric_predicate_test_count = k;
    end
    report.status = 'PASS_NATIVE_MATLAB_PREFLIGHT_AND_NUMERIC_PREDICATE';
    save_report(fullfile(auditRoot,'MATLAB_PREFLIGHT_RESULT.json'),report);
    fprintf('matlab_native_syntax_preflight=PASS\n');
    fprintf('matlab_code_files_checked=%d\n',report.code_files_checked);
    fprintf('original_syntax_error_negative_control=True\n');
    fprintf('native_numeric_predicate_tests=%d\n',report.numeric_predicate_test_count);
catch ME
    report.status = 'HOLD_NATIVE_MATLAB_PREFLIGHT';
    report.error_identifier = ME.identifier;
    save_report(fullfile(auditRoot,'MATLAB_PREFLIGHT_RESULT.json'),report);
    rethrow(ME);
end
end

function save_report(file,report)
fid=fopen(file,'w','n','UTF-8');
assert(fid>=0,'Cannot write native preflight result.');
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(report,'PrettyPrint',true));
end
