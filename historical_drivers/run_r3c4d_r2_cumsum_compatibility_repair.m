function run_r3c4d_r2_cumsum_compatibility_repair(artifactRoot, r3c3bRoot, analysisRasterPath, georefGateCsv, georefCenterCsv, officialRecordCsv, centerLat, centerLon)
% R3C4D-R2 cumulative-sum API compatibility repair plus corrected full-1000 population-only exposure propagation.
% Uses frozen R3C3B threshold radii and the source-bound native-grid analysis raster.
% Does not rerun entry, fragmentation, Monte Carlo source generation, or overpressure.
% Does not resample, reproject, redistribute, or modify the population raster.

fprintf('=== CHELYABINSK R3C4D-R2 CUMSUM API COMPATIBILITY REPAIR + FULL1000 POPULATION EXPOSURE ===\n');
fprintf('prior_R3C4D_R1_population_exposure_propagation_executed=false\n');
fprintf('root_cause=cumsum_does_not_accept_double_as_an_option\n');
fprintf('repair=explicit_double_cast_before_cumsum_with_parser_repair_preserved\n');
fprintf('scientific_entry_propagation_requested=false\n');
fprintf('fragmentation_propagation_requested=false\n');
fprintf('monte_carlo_source_propagation_requested=false\n');
fprintf('overpressure_propagation_requested=false\n');
fprintf('population_exposure_propagation_requested=true\n');
fprintf('population_raster_redistribution_requested=false\n');
fprintf('candidate_manuscript_modification_requested=false\n');
fprintf('canonical_modification_requested=false\n');
fprintf('exposure_inclusion_rule=CELL_CENTROID_WITHIN_GEODESIC_RADIUS\n');
fprintf('cell_area_rule=SPHERICAL_LATITUDE_DEPENDENT_NATIVE_GRID_CELL_AREA\n');
fprintf('artifact_root=%s\n', artifactRoot);

outRoot = fullfile(artifactRoot,'02_outputs');
figRoot = fullfile(outRoot,'figures');
tableRoot = fullfile(outRoot,'tables');
if ~exist(outRoot,'dir'); mkdir(outRoot); end
if ~exist(figRoot,'dir'); mkdir(figRoot); end
if ~exist(tableRoot,'dir'); mkdir(tableRoot); end

% -------------------------------------------------------------------------
% Bind passed geospatial eligibility.
% -------------------------------------------------------------------------
G = readGateCsvRobust(georefGateCsv);
requiredGeoref = [
    "rasterinfo_api_compatibility_pass"
    "native_grid_alignment_pass"
    "masked_extract_parity_pass"
    "metadata_binding_complete"
    "scenario_center_inside_analysis_bounds"
    "population_exposure_propagation_eligible"
    ];
fprintf('geospatial_gate_parser_row_count=%d\n',height(G));
fprintf('geospatial_gate_parser_unique_metric_count=%d\n',numel(unique(G.metric)));
if height(G) < 18
    error('Geospatial gate CSV contains fewer than 18 parsed metric rows.');
end
for k = 1:numel(requiredGeoref)
    v = gateValue(G,requiredGeoref(k));
    fprintf('bound_geospatial_gate_%s=%s\n',requiredGeoref(k),v);
    if ~asBool(v)
        error('Required geospatial gate is not true after robust parsing: %s=%s',requiredGeoref(k),v);
    end
end
geospatialGateParserRepairPass = true;

C = readtable(georefCenterCsv,'TextType','string');
if height(C) ~= 1
    error('Scenario-center binding must contain exactly one row.');
end
boundLat = numericScalar(C.scenario_center_lat_deg);
boundLon = numericScalar(C.scenario_center_lon_deg);
if abs(boundLat-centerLat) > 1e-10 || abs(boundLon-centerLon) > 1e-10
    error('Scenario-center input does not match frozen georeference binding.');
end
if ~asBool(string(C.inside_analysis_bounds(1)))
    error('Frozen scenario center is not inside the analysis bounds.');
end

% -------------------------------------------------------------------------
% Discover and normalize frozen R3C3B radius data.
% -------------------------------------------------------------------------
[Rlong, radiusSourcePath, radiusSchema, radiusMappingBasis] = discoverRadiusData(r3c3bRoot);
expectedThresholds = [3.45; 1.00; 0.50; 0.20];

if height(Rlong) ~= 4000
    error('Normalized threshold-radius table must contain 4000 rows; observed %d.',height(Rlong));
end
uniqueRuns = unique(Rlong.run_id);
if numel(uniqueRuns) ~= 1000
    error('Expected 1000 unique run IDs; observed %d.',numel(uniqueRuns));
end

thresholdSet = sort(unique(Rlong.threshold_kpa),'descend');
if numel(thresholdSet) ~= 4 || max(abs(thresholdSet-expectedThresholds)) > 1e-10
    error('Threshold set does not match [3.45, 1.00, 0.50, 0.20] kPa.');
end
if any(~isfinite(Rlong.radius_km)) || any(Rlong.radius_km <= 0)
    error('All threshold radii must be finite and positive.');
end

% Canonical row order: run_id ascending, pressure descending.
Rlong = sortrows(Rlong,{'run_id','threshold_kpa'},{'ascend','descend'});
runCountPass = true;
thresholdCountPass = true;
radiusOrderPass = true;
for i = 1:numel(uniqueRuns)
    idx = Rlong.run_id == uniqueRuns(i);
    if nnz(idx) ~= 4
        thresholdCountPass = false;
        radiusOrderPass = false;
        break
    end
    rr = Rlong.radius_km(idx);
    pp = Rlong.threshold_kpa(idx);
    [pp,ord] = sort(pp,'descend');
    rr = rr(ord);
    if max(abs(pp-expectedThresholds)) > 1e-10 || any(diff(rr) <= 0)
        radiusOrderPass = false;
        break
    end
end
if ~runCountPass || ~thresholdCountPass || ~radiusOrderPass
    error('R3C3B threshold-radius per-run contract failed.');
end

% -------------------------------------------------------------------------
% Read the georeferenced analysis raster without modifying it.
% -------------------------------------------------------------------------
[A,Ra] = readgeoraster(analysisRasterPath,'OutputType','double','StandardizeMissing',true);
if ndims(A)>2; A=A(:,:,1); end
if isempty(Ra) || ~contains(lower(string(class(Ra))),"geographic")
    error('Analysis raster must use a geographic raster-reference object.');
end

valid = isfinite(A) & A >= 0;
if nnz(valid) ~= 684872
    error('Analysis valid-cell count changed: expected 684872; observed %d.',nnz(valid));
end
if any(A(valid) < 0)
    error('Negative finite population values are not permitted.');
end

[validRows,validCols] = find(valid);
[latValid,lonValid] = intrinsicToGeographic(Ra,validCols,validRows);
populationValid = A(valid);

% Great-circle distance from the frozen scenario center to each valid-cell centroid.
earthRadiusKm = 6371.0088;
distanceKm = haversineKm(centerLat,centerLon,latValid,lonValid,earthRadiusKm);

% Native geographic-cell area varies with latitude.
dLat = double(Ra.CellExtentInLatitude);
dLon = double(Ra.CellExtentInLongitude);
latNorth = min(90,latValid + dLat/2);
latSouth = max(-90,latValid - dLat/2);
cellAreaKm2 = earthRadiusKm^2 * deg2rad(dLon) .* ...
    (sin(deg2rad(latNorth)) - sin(deg2rad(latSouth)));

if any(~isfinite(distanceKm)) || any(distanceKm < 0)
    error('Non-finite or negative cell-centroid distances detected.');
end
if any(~isfinite(cellAreaKm2)) || any(cellAreaKm2 <= 0)
    error('Non-finite or non-positive native cell areas detected.');
end

domainPopulation = sum(populationValid,'double');
domainAreaKm2 = sum(cellAreaKm2,'double');
if abs(domainPopulation - 87794916.9867744) > 1e-3
    error('Population-domain sum does not match the frozen numeric audit.');
end

[distanceSorted,sortOrder] = sort(distanceKm,'ascend');
populationSorted = populationValid(sortOrder);
areaSorted = cellAreaKm2(sortOrder);
populationSorted = double(populationSorted);
areaSorted = double(areaSorted);
cumulativePopulation = cumsum(populationSorted);
cumulativeArea = cumsum(areaSorted);
cumsumApiCompatibilityRepairPass = ...
    isa(cumulativePopulation,'double') && isa(cumulativeArea,'double') && ...
    numel(cumulativePopulation)==numel(populationSorted) && ...
    numel(cumulativeArea)==numel(areaSorted);
cumulativeCells = (1:numel(distanceSorted))';

% -------------------------------------------------------------------------
% Full-1000 exposure propagation using one sorted-distance sweep.
% -------------------------------------------------------------------------
queryRadii = Rlong.radius_km;
[queryPopulation,queryArea,queryCellCount] = cumulativeAtRadii( ...
    distanceSorted,cumulativePopulation,cumulativeArea,cumulativeCells,queryRadii);

exposureFraction = queryPopulation ./ domainPopulation;
E = table(Rlong.run_id,Rlong.threshold_kpa,Rlong.radius_km, ...
    queryArea,queryPopulation,exposureFraction,queryCellCount, ...
    'VariableNames',{'run_id','threshold_kpa','radius_km', ...
    'affected_area_km2','exposed_population','exposure_fraction','included_cell_count'});

fullResultsPath = fullfile(outRoot,'R3C4D_full1000_population_exposure_results.csv');
writetable(E,fullResultsPath);

% -------------------------------------------------------------------------
% Nominal G4 footprint. Prefer an upstream nominal-radius table; otherwise
% reproduce the frozen fixed-anchor contract.
% -------------------------------------------------------------------------
[nominalThreshold,nominalRadius,nominalSource] = discoverNominalRadii(r3c3bRoot);
[nomPop,nomArea,nomCells] = cumulativeAtRadii( ...
    distanceSorted,cumulativePopulation,cumulativeArea,cumulativeCells,nominalRadius);
nomFrac = nomPop ./ domainPopulation;
N = table(nominalThreshold,nominalRadius,nomArea,nomPop,nomFrac,nomCells, ...
    repmat(string(nominalSource),4,1), ...
    'VariableNames',{'threshold_kpa','radius_km','affected_area_km2', ...
    'exposed_population','exposure_fraction','included_cell_count','radius_source'});
nominalPath = fullfile(outRoot,'R3C4D_nominal_population_exposure.csv');
writetable(N,nominalPath);

% -------------------------------------------------------------------------
% Per-threshold ensemble summary.
% -------------------------------------------------------------------------
summaryRows = table();
for k = 1:numel(expectedThresholds)
    p = expectedThresholds(k);
    X = E(E.threshold_kpa==p,:);
    row = table();
    row.threshold_kpa = p;
    row.run_count = height(X);
    row.radius_min_km = min(X.radius_km);
    row.radius_p10_km = sampleQuantile(X.radius_km,0.10);
    row.radius_median_km = sampleQuantile(X.radius_km,0.50);
    row.radius_p90_km = sampleQuantile(X.radius_km,0.90);
    row.radius_max_km = max(X.radius_km);
    row.area_min_km2 = min(X.affected_area_km2);
    row.area_p10_km2 = sampleQuantile(X.affected_area_km2,0.10);
    row.area_median_km2 = sampleQuantile(X.affected_area_km2,0.50);
    row.area_p90_km2 = sampleQuantile(X.affected_area_km2,0.90);
    row.area_max_km2 = max(X.affected_area_km2);
    row.population_min = min(X.exposed_population);
    row.population_p10 = sampleQuantile(X.exposed_population,0.10);
    row.population_median = sampleQuantile(X.exposed_population,0.50);
    row.population_p90 = sampleQuantile(X.exposed_population,0.90);
    row.population_max = max(X.exposed_population);
    row.fraction_min = min(X.exposure_fraction);
    row.fraction_p10 = sampleQuantile(X.exposure_fraction,0.10);
    row.fraction_median = sampleQuantile(X.exposure_fraction,0.50);
    row.fraction_p90 = sampleQuantile(X.exposure_fraction,0.90);
    row.fraction_max = max(X.exposure_fraction);
    summaryRows = [summaryRows;row]; %#ok<AGROW>
end
summaryPath = fullfile(outRoot,'R3C4D_full1000_population_exposure_summary.csv');
writetable(summaryRows,summaryPath);

% -------------------------------------------------------------------------
% Domain and source bindings.
% -------------------------------------------------------------------------
D = table(string(analysisRasterPath),string(class(Ra)), ...
    double(Ra.RasterSize(1)),double(Ra.RasterSize(2)), ...
    nnz(valid),nnz(valid & A>0),nnz(valid & A==0), ...
    domainPopulation,domainAreaKm2,dLat,dLon,centerLat,centerLon, ...
    string("people_per_native_grid_cell"), ...
    string("cell centroid inside geodesic circular threshold"), ...
    string("spherical latitude-dependent cell area"), ...
    'VariableNames',{'analysis_raster_path','raster_reference_class', ...
    'raster_rows','raster_columns','valid_cell_count','positive_cell_count', ...
    'zero_cell_count','domain_population','domain_area_km2', ...
    'cell_extent_lat_deg','cell_extent_lon_deg','scenario_center_lat_deg', ...
    'scenario_center_lon_deg','population_semantics','inclusion_rule','cell_area_rule'});
writetable(D,fullfile(outRoot,'R3C4D_population_domain_summary.csv'));

B = table(string(radiusSourcePath),string(radiusSchema),string(radiusMappingBasis), ...
    height(Rlong),numel(uniqueRuns),numel(unique(Rlong.threshold_kpa)), ...
    min(Rlong.radius_km),max(Rlong.radius_km), ...
    string(officialRecordCsv), ...
    'VariableNames',{'radius_source_path','normalized_schema','threshold_mapping_basis', ...
    'normalized_row_count','unique_run_count','threshold_count', ...
    'minimum_radius_km','maximum_radius_km','official_population_record_path'});
writetable(B,fullfile(outRoot,'R3C4D_upstream_radius_source_binding.csv'));

% -------------------------------------------------------------------------
% Scientific gates.
% -------------------------------------------------------------------------
allExposureFinite = all(isfinite(E.affected_area_km2)) && ...
    all(isfinite(E.exposed_population)) && all(isfinite(E.exposure_fraction));
allExposureNonnegative = all(E.affected_area_km2>=0) && ...
    all(E.exposed_population>=0) && all(E.exposure_fraction>=0);
allExposureWithinDomain = all(E.affected_area_km2<=domainAreaKm2+1e-8) && ...
    all(E.exposed_population<=domainPopulation+1e-6) && all(E.exposure_fraction<=1+1e-12);

populationMonotonicPass = true;
areaMonotonicPass = true;
for i = 1:numel(uniqueRuns)
    X = E(E.run_id==uniqueRuns(i),:);
    X = sortrows(X,'threshold_kpa','descend');
    if any(diff(X.exposed_population)<-1e-8)
        populationMonotonicPass = false;
    end
    if any(diff(X.affected_area_km2)<-1e-8)
        areaMonotonicPass = false;
    end
end

nominalMonotonicPass = all(diff(N.exposed_population)>=-1e-8) && ...
    all(diff(N.affected_area_km2)>=-1e-8);
summaryOrderPass = all(diff(summaryRows.population_median)>=-1e-8) && ...
    all(diff(summaryRows.area_median_km2)>=-1e-8);

status = "FAIL_R3C4D_R2_CUMSUM_API_COMPATIBILITY_REPAIR_AND_FULL1000_POPULATION_EXPOSURE_PROPAGATION";
allGatesPass = geospatialGateParserRepairPass && cumsumApiCompatibilityRepairPass && height(E)==4000 && numel(uniqueRuns)==1000 && ...
    radiusOrderPass && allExposureFinite && allExposureNonnegative && ...
    allExposureWithinDomain && populationMonotonicPass && areaMonotonicPass && ...
    nominalMonotonicPass && summaryOrderPass;
if allGatesPass
    status = "PASS_R3C4D_R2_CUMSUM_API_COMPATIBILITY_REPAIR_AND_FULL1000_POPULATION_EXPOSURE_PROPAGATION";
end

gateMetric = [
    "status"
    "geospatial_gate_parser_repair_pass"
    "cumsum_api_compatibility_repair_pass"
    "population_exposure_propagation_executed"
    "scientific_entry_propagation_executed"
    "overpressure_propagation_executed"
    "population_raster_redistributed"
    "source_realization_count"
    "threshold_count"
    "exposure_result_row_count"
    "all_source_radii_finite_positive"
    "per_run_inverse_pressure_radius_order_pass"
    "all_exposure_values_finite"
    "all_exposure_values_nonnegative"
    "all_exposure_values_within_domain"
    "per_run_exposed_population_monotonic_pass"
    "per_run_affected_area_monotonic_pass"
    "nominal_exposure_monotonic_pass"
    "summary_median_monotonic_pass"
    "domain_population_frozen_audit_match"
    "population_statistics_publishable"
    "candidate_manuscript_modification_executed"
    "canonical_modification_executed"
    ];
gateValueOut = [
    status
    string(geospatialGateParserRepairPass)
    string(cumsumApiCompatibilityRepairPass)
    "true"
    "false"
    "false"
    "false"
    string(numel(uniqueRuns))
    string(numel(unique(Rlong.threshold_kpa)))
    string(height(E))
    string(all(isfinite(Rlong.radius_km)) && all(Rlong.radius_km>0))
    string(radiusOrderPass)
    string(allExposureFinite)
    string(allExposureNonnegative)
    string(allExposureWithinDomain)
    string(populationMonotonicPass)
    string(areaMonotonicPass)
    string(nominalMonotonicPass)
    string(summaryOrderPass)
    string(abs(domainPopulation-87794916.9867744)<=1e-3)
    string(allGatesPass)
    "false"
    "false"
    ];
GateOut = table(gateMetric,gateValueOut,'VariableNames',{'metric','value'});
writetable(GateOut,fullfile(outRoot,'R3C4D_population_exposure_gate_adjudication.csv'));
writelines(status,fullfile(outRoot,'R3C4D_status.txt'));

% -------------------------------------------------------------------------
% Publication-candidate tables. These are not inserted into the manuscript.
% -------------------------------------------------------------------------
writeNominalLatex(N,fullfile(tableRoot,'R3C4D_nominal_population_exposure_table.tex'));
writeSummaryLatex(summaryRows,fullfile(tableRoot,'R3C4D_full1000_population_exposure_summary_table.tex'));

% -------------------------------------------------------------------------
% Diagnostic/publication-candidate figures.
% -------------------------------------------------------------------------
makeNominalMap(A,Ra,N,centerLat,centerLon,figRoot);
makeExposureDistributionFigure(E,expectedThresholds,figRoot);
makeExposureSummaryFigure(summaryRows,figRoot);

fprintf('\n=== R3C4D-R2 RESULT ===\n');
fprintf('geospatial_gate_parser_repair_pass=%s\n',tf(geospatialGateParserRepairPass));
fprintf('cumsum_api_compatibility_repair_pass=%s\n',tf(cumsumApiCompatibilityRepairPass));
fprintf('status=%s\n',status);
fprintf('population_exposure_propagation_executed=true\n');
fprintf('source_realization_count=%d\n',numel(uniqueRuns));
fprintf('threshold_count=%d\n',numel(unique(Rlong.threshold_kpa)));
fprintf('exposure_result_row_count=%d\n',height(E));
fprintf('domain_population=%.10f\n',domainPopulation);
fprintf('domain_area_km2=%.10f\n',domainAreaKm2);
fprintf('all_source_radii_finite_positive=%s\n',tf(all(isfinite(Rlong.radius_km)) && all(Rlong.radius_km>0)));
fprintf('per_run_inverse_pressure_radius_order_pass=%s\n',tf(radiusOrderPass));
fprintf('all_exposure_values_finite=%s\n',tf(allExposureFinite));
fprintf('all_exposure_values_within_domain=%s\n',tf(allExposureWithinDomain));
fprintf('per_run_exposed_population_monotonic_pass=%s\n',tf(populationMonotonicPass));
fprintf('per_run_affected_area_monotonic_pass=%s\n',tf(areaMonotonicPass));
fprintf('nominal_exposure_monotonic_pass=%s\n',tf(nominalMonotonicPass));
fprintf('population_statistics_publishable=%s\n',tf(allGatesPass));
fprintf('population_raster_redistributed=false\n');
fprintf('candidate_manuscript_modification_executed=false\n');
fprintf('canonical_modification_executed=false\n');

for k=1:height(N)
    fprintf('nominal_threshold_%.2f_kpa_radius_km=%.10f\n',N.threshold_kpa(k),N.radius_km(k));
    fprintf('nominal_threshold_%.2f_kpa_exposed_population=%.10f\n',N.threshold_kpa(k),N.exposed_population(k));
    fprintf('nominal_threshold_%.2f_kpa_exposure_fraction=%.12f\n',N.threshold_kpa(k),N.exposure_fraction(k));
end

if ~allGatesPass
    error('R3C4D-R2 population-exposure scientific gates failed.');
end
end

% =========================================================================
% Radius discovery and normalization
% =========================================================================
function [Rlong,sourcePath,schema,mappingBasis] = discoverRadiusData(root)
files = dir(fullfile(root,'**','*.csv'));
bestScore = -Inf;
Rlong = table();
sourcePath = "";
schema = "";
mappingBasis = "";

for i=1:numel(files)
    path = fullfile(files(i).folder,files(i).name);
    nameLower = lower(string(files(i).name));
    if ~(contains(nameLower,"radius") || contains(nameLower,"overpressure") || contains(nameLower,"threshold"))
        continue
    end
    try
        T = readtable(path,'VariableNamingRule','preserve','TextType','string');
    catch
        continue
    end
    try
        [candidate,cSchema,cBasis,cScore] = normalizeRadiusTable(T,files(i).name);
    catch
        continue
    end
    if cScore > bestScore
        bestScore = cScore;
        Rlong = candidate;
        sourcePath = string(path);
        schema = cSchema;
        mappingBasis = cBasis;
    end
end

if isempty(Rlong) || height(Rlong)~=4000
    error('Unable to discover a valid 1000x4 R3C3B threshold-radius source table under %s.',root);
end
end

function [Rlong,schema,basis,score] = normalizeRadiusTable(T,fileName)
names = string(T.Properties.VariableNames);
canon = canonicalNames(names);
Rlong = table();
schema = "";
basis = "";
score = -Inf;

% Long schema.
radiusIdx = find(contains(canon,"radius"),1);
thresholdIdx = find(contains(canon,"threshold") | contains(canon,"pressurekpa") | canon=="kpa",1);
runIdx = find(contains(canon,"runid") | canon=="run" | canon=="realizationid",1);
if ~isempty(radiusIdx) && ~isempty(thresholdIdx) && height(T)>=4000
    radius = numericColumn(T{:,radiusIdx});
    threshold = numericColumn(T{:,thresholdIdx});
    if isempty(runIdx)
        run = repelem((1:(height(T)/4))',4);
        runBasis = "GENERATED_SEQUENTIAL";
    else
        run = numericColumn(T{:,runIdx});
        runBasis = "SOURCE_COLUMN";
    end
    good = isfinite(radius) & isfinite(threshold) & isfinite(run);
    radius=radius(good); threshold=threshold(good); run=run(good);
    if numel(radius)==4000 && numel(unique(run))==1000 && numel(unique(threshold))==4
        Rlong = table(run,threshold,radius,'VariableNames',{'run_id','threshold_kpa','radius_km'});
        schema = "LONG";
        basis = "EXPLICIT_THRESHOLD_COLUMN+"+runBasis;
        score = 100 + contains(lower(string(fileName)),"r3c3b")*20;
        return
    end
end

% Wide schema with 1000 rows and four identifiable radius columns.
if height(T)==1000
    radiusCandidates = find((contains(canon,"radius") | startsWith(canon,"r")) & ~contains(canon,"run"));
    mapped = nan(numel(radiusCandidates),1);
    for j=1:numel(radiusCandidates)
        mapped(j)=thresholdFromName(canon(radiusCandidates(j)));
    end
    explicitMask=isfinite(mapped);
    if nnz(explicitMask)>=4
        idx = radiusCandidates(explicitMask);
        th = mapped(explicitMask);
        [th,uniqueIdx] = unique(th,'stable');
        idx = idx(uniqueIdx);
        if numel(th)==4
            [th,ord]=sort(th,'descend'); idx=idx(ord);
            runIdx = find(contains(canon,"runid") | canon=="run" | canon=="realizationid",1);
            if isempty(runIdx)
                run=(1:1000)';
                runBasis="GENERATED_SEQUENTIAL";
            else
                run=numericColumn(T{:,runIdx});
                runBasis="SOURCE_COLUMN";
            end
            radiusMatrix=zeros(1000,4);
            for j=1:4
                radiusMatrix(:,j)=numericColumn(T{:,idx(j)});
            end
            [runLong,thLong,rLong]=wideToLong(run,th,radiusMatrix);
            Rlong=table(runLong,thLong,rLong,'VariableNames',{'run_id','threshold_kpa','radius_km'});
            schema="WIDE";
            basis="EXPLICIT_COLUMN_NAME+"+runBasis;
            score=95 + contains(lower(string(fileName)),"r3c3b")*20;
            return
        end
    end

    % Conservative fallback: exactly four radius-like numeric columns and
    % strict radius ordering in every row. Threshold mapping is by median order.
    numericRadiusIdx=[];
    for j=radiusCandidates(:)'
        x=numericColumn(T{:,j});
        if numel(x)==1000 && all(isfinite(x)) && all(x>0)
            numericRadiusIdx(end+1)=j; %#ok<AGROW>
        end
    end
    if numel(numericRadiusIdx)==4
        M=zeros(1000,4);
        for j=1:4; M(:,j)=numericColumn(T{:,numericRadiusIdx(j)}); end
        med=median(M,1);
        [~,ord]=sort(med,'ascend');
        M=M(:,ord);
        if all(all(diff(M,1,2)>0))
            runIdx = find(contains(canon,"runid") | canon=="run" | canon=="realizationid",1);
            if isempty(runIdx); run=(1:1000)'; runBasis="GENERATED_SEQUENTIAL";
            else; run=numericColumn(T{:,runIdx}); runBasis="SOURCE_COLUMN"; end
            th=[3.45;1.00;0.50;0.20];
            [runLong,thLong,rLong]=wideToLong(run,th,M);
            Rlong=table(runLong,thLong,rLong,'VariableNames',{'run_id','threshold_kpa','radius_km'});
            schema="WIDE";
            basis="STRICT_MEDIAN_ORDER_FALLBACK+"+runBasis;
            score=60 + contains(lower(string(fileName)),"r3c3b")*20;
            return
        end
    end
end
error('Not a valid threshold-radius source table.');
end

function [runLong,thLong,rLong]=wideToLong(run,thresholds,M)
n=numel(run);
runLong=repelem(run,4);
thLong=repmat(thresholds(:),n,1);
rLong=reshape(M.',[],1);
end

function t=thresholdFromName(name)
s=lower(regexprep(string(name),'[^a-z0-9]',''));
t=NaN;
if contains(s,"345"); t=3.45;
elseif contains(s,"100") || contains(s,"1kpa"); t=1.00;
elseif contains(s,"050") || contains(s,"05kpa") || contains(s,"50kpa"); t=0.50;
elseif contains(s,"020") || contains(s,"02kpa") || contains(s,"20kpa"); t=0.20;
end
end

function [threshold,radius,source] = discoverNominalRadii(root)
threshold=[3.45;1.00;0.50;0.20];
radius=[];
source="";
files=dir(fullfile(root,'**','*.csv'));
for i=1:numel(files)
    n=lower(string(files(i).name));
    if ~(contains(n,"nominal") && contains(n,"radius")); continue; end
    path=fullfile(files(i).folder,files(i).name);
    try
        T=readtable(path,'VariableNamingRule','preserve','TextType','string');
    catch
        continue
    end
    names=canonicalNames(string(T.Properties.VariableNames));
    rIdx=find(contains(names,"radius"),1);
    pIdx=find(contains(names,"threshold") | contains(names,"pressurekpa"),1);
    if ~isempty(rIdx) && ~isempty(pIdx) && height(T)>=4
        rr=numericColumn(T{:,rIdx}); pp=numericColumn(T{:,pIdx});
        good=isfinite(rr)&isfinite(pp);
        rr=rr(good); pp=pp(good);
        if numel(unique(pp))==4
            [pp,ord]=sort(pp,'descend'); rr=rr(ord);
            if max(abs(pp-threshold))<1e-10 && all(rr>0)
                radius=rr(:); source="UPSTREAM_NOMINAL_RADIUS_TABLE:"+string(path); return
            end
        end
    end
end

% Frozen fixed nominal anchor contract from R3C3A-R1/R3C3B.
C_op=1215.90123801035; % kPa km^3 kt^-1
E_c=331.6302;          % kt, corrected nominal G4 coupled burst energy
h_burst=23.75;         % km
slant=(C_op*E_c./threshold).^(1/3);
radius=sqrt(max(0,slant.^2-h_burst.^2));
source="FROZEN_GLOBAL_FIXED_NOMINAL_ANCHOR_CONTRACT";
end

% =========================================================================
% Exposure utilities
% =========================================================================
function [popAt,areaAt,countAt] = cumulativeAtRadii(dist,cumPop,cumArea,cumCount,radii)
[rSorted,ord]=sort(radii(:),'ascend');
popSorted=zeros(size(rSorted)); areaSorted=zeros(size(rSorted)); countSorted=zeros(size(rSorted));
j=0; n=numel(dist);
for q=1:numel(rSorted)
    while j<n && dist(j+1)<=rSorted(q)
        j=j+1;
    end
    if j>0
        popSorted(q)=cumPop(j);
        areaSorted(q)=cumArea(j);
        countSorted(q)=cumCount(j);
    end
end
popAt=zeros(size(rSorted)); areaAt=zeros(size(rSorted)); countAt=zeros(size(rSorted));
popAt(ord)=popSorted; areaAt(ord)=areaSorted; countAt(ord)=countSorted;
end

function d=haversineKm(lat0,lon0,lat,lon,R)
phi1=deg2rad(lat0); phi2=deg2rad(lat);
dphi=phi2-phi1; dlambda=deg2rad(lon-lon0);
a=sin(dphi/2).^2 + cos(phi1).*cos(phi2).*sin(dlambda/2).^2;
a=min(1,max(0,a));
d=2*R*asin(sqrt(a));
end

function q=sampleQuantile(x,p)
x=sort(x(isfinite(x)));
n=numel(x);
if n==0; q=NaN; return; end
if n==1; q=x; return; end
h=(n-1)*p+1;
lo=floor(h); hi=ceil(h);
if lo==hi; q=x(lo); else; q=x(lo)+(h-lo)*(x(hi)-x(lo)); end
end

% =========================================================================
% Figures and LaTeX
% =========================================================================
function makeNominalMap(A,R,N,lat0,lon0,figRoot)
fig=figure('Visible','off','Color','w','Position',[100 100 1000 700]);
popDisplay=A; popDisplay(~isfinite(popDisplay)|popDisplay<0)=NaN;
imagesc(R.LongitudeLimits,R.LatitudeLimits,log10(popDisplay+1));
set(gca,'YDir','normal','FontName','Times New Roman','FontSize',11);
axis tight; hold on;
xlabel('Longitude [deg]'); ylabel('Latitude [deg]');
title('Corrected nominal population-exposure screening zones');
cb=colorbar; ylabel(cb,'log_{10}(population per native grid cell + 1)');
for k=height(N):-1:1
    [latC,lonC]=smallCircle(lat0,lon0,N.radius_km(k),6371.0088);
    plot(lonC,latC,'LineWidth',1.4,'DisplayName',sprintf('%.2f kPa: %.2f km',N.threshold_kpa(k),N.radius_km(k)));
end
plot(lon0,lat0,'o','MarkerSize',6,'MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Scenario centre');
legend('Location','best');
grid on; box on;
exportgraphics(fig,fullfile(figRoot,'fig_R3C4D_nominal_population_exposure_map.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(figRoot,'fig_R3C4D_nominal_population_exposure_map.png'),'Resolution',240);
close(fig);
end

function makeExposureDistributionFigure(E,thresholds,figRoot)
fig=figure('Visible','off','Color','w','Position',[100 100 1100 760]);
tl=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
for k=1:4
    nexttile;
    x=E.exposed_population(E.threshold_kpa==thresholds(k))/1e6;
    histogram(x,28);
    xlabel('Exposed population [million]');
    ylabel('Realizations');
    title(sprintf('%.2f kPa',thresholds(k)));
    grid on; box on;
end
title(tl,'Corrected full-1000 exposed-population distributions');
exportgraphics(fig,fullfile(figRoot,'fig_R3C4D_full1000_exposed_population_distributions.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(figRoot,'fig_R3C4D_full1000_exposed_population_distributions.png'),'Resolution',240);
close(fig);
end

function makeExposureSummaryFigure(S,figRoot)
fig=figure('Visible','off','Color','w','Position',[100 100 900 620]);
x=1:height(S);
med=S.population_median/1e6;
lo=(S.population_median-S.population_p10)/1e6;
hi=(S.population_p90-S.population_median)/1e6;
errorbar(x,med,lo,hi,'o-','LineWidth',1.5,'MarkerSize',6);
set(gca,'XTick',x,'XTickLabel',compose('%.2f kPa',S.threshold_kpa), ...
    'FontName','Times New Roman','FontSize',11);
ylabel('Exposed population [million]');
xlabel('Overpressure threshold');
title('Corrected population exposure: median and P10--P90');
grid on; box on;
exportgraphics(fig,fullfile(figRoot,'fig_R3C4D_population_exposure_summary.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(figRoot,'fig_R3C4D_population_exposure_summary.png'),'Resolution',240);
close(fig);
end

function [lat,lon]=smallCircle(lat0,lon0,rKm,R)
bearing=linspace(0,2*pi,721);
delta=rKm/R;
phi1=deg2rad(lat0); lam1=deg2rad(lon0);
phi2=asin(sin(phi1).*cos(delta)+cos(phi1).*sin(delta).*cos(bearing));
lam2=lam1+atan2(sin(bearing).*sin(delta).*cos(phi1),cos(delta)-sin(phi1).*sin(phi2));
lat=rad2deg(phi2); lon=rad2deg(lam2);
end

function writeNominalLatex(N,path)
fid=fopen(path,'w');
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'\\begin{tabular}{rrrrr}\n');
fprintf(fid,'\\hline\n');
fprintf(fid,'Threshold [kPa] & Radius [km] & Area [km$^2$] & Exposed population & Fraction [\\%%] \\\\\n');
fprintf(fid,'\\hline\n');
for i=1:height(N)
    fprintf(fid,'%.2f & %.2f & %.2f & %.0f & %.4f \\\\\n',N.threshold_kpa(i),N.radius_km(i), ...
        N.affected_area_km2(i),N.exposed_population(i),100*N.exposure_fraction(i));
end
fprintf(fid,'\\hline\n\\end{tabular}\n');
end

function writeSummaryLatex(S,path)
fid=fopen(path,'w');
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'\\begin{tabular}{rrrrrrr}\n');
fprintf(fid,'\\hline\n');
fprintf(fid,'Threshold & Radius P10 & Radius med. & Radius P90 & Pop. P10 & Pop. med. & Pop. P90 \\\\\n');
fprintf(fid,'\\hline\n');
for i=1:height(S)
    fprintf(fid,'%.2f & %.2f & %.2f & %.2f & %.0f & %.0f & %.0f \\\\\n', ...
        S.threshold_kpa(i),S.radius_p10_km(i),S.radius_median_km(i),S.radius_p90_km(i), ...
        S.population_p10(i),S.population_median(i),S.population_p90(i));
end
fprintf(fid,'\\hline\n\\end{tabular}\n');
end

% =========================================================================
% Parsing helpers
% =========================================================================
function c=canonicalNames(names)
c=lower(regexprep(string(names),'[^a-zA-Z0-9]',''));
end

function x=numericColumn(v)
if isnumeric(v); x=double(v(:)); return; end
s=string(v(:));
x=str2double(regexprep(s,'[^0-9eE+\-.]',''));
end

function T=readGateCsvRobust(path)
raw = readlines(path);
raw = strip(erase(raw,char(65279)));
raw(raw=="") = [];
if numel(raw)<2
    error('Geospatial gate CSV is empty or has no data rows: %s',path);
end
header = lower(strip(split(raw(1),",")));
if numel(header)<2 || header(1)~="metric" || header(2)~="value"
    error('Unexpected geospatial gate CSV header: %s',raw(1));
end
metric = strings(numel(raw)-1,1);
value = strings(numel(raw)-1,1);
for i=2:numel(raw)
    parts = split(raw(i),",");
    if numel(parts)>=2
        metric(i-1) = strip(erase(erase(parts(1),'"'),''''));
        value(i-1) = strip(erase(erase(parts(2),'"'),''''));
    end
end
keep = strlength(metric)>0;
T = table(metric(keep),value(keep),'VariableNames',{'metric','value'});
end

function v=gateValue(T,metric)
keys = lower(strtrim(string(T.metric)));
target = lower(strtrim(string(metric)));
idx = find(keys==target,1);
if isempty(idx)
    v="";
else
    v=strtrim(string(T.value(idx)));
end
end

function b=asBool(v)
s=lower(strtrim(string(v)));
b=any(s==["true","1","yes"]);
end

function x=numericScalar(v)
if isnumeric(v); x=double(v(1)); else; x=str2double(string(v(1))); end
end

function s=tf(v)
if v; s='true'; else; s='false'; end
end
